# Preprocess to prepare for zUMI
# Mapped to both human and mouse
# William Ho

# Setup ------------------------------------------------------------------------

### WH: in here: mkdir output
### WH: move "helper_functions.R" to code dir
### WH: move samplesheet (Excel file) to data/samplesheets
### WH: update renv.lock (from similar project) and run renv::restore()

library(scPipe)
library(Rsubread)
library(here)
library(readxl)
library(dplyr)
library(janitor)
library(tidyr)

options("mc.cores" = 8L)

# Construct NN198 sample sheet -------------------------------------------------

file_nn198 <- here(
  "data",
  "sample_sheets",
  "JN_DZ_C051 NN198_SarahBest_Jin_SeqPrimer layout_Oct20.xlsx")

# NOTE: Header row is split across 2 lines, which I combine into 1 before
#       reading in the rest of the spreadsheet.
header_row <- read_excel(
  path = file_nn198,
  sheet = "Sample & Index",
  skip = 2,
  n_max = 1)

# NOTE: FACS data in columns >= "K"
facs_data_idx <- seq(which(LETTERS == "K"), ncol(header_row))
header_row <- c(
  paste0(colnames(header_row[, -facs_data_idx]), header_row[1, -facs_data_idx]),
  unlist(header_row[1, facs_data_idx], use.names = FALSE))
header_row <- gsub("^\\.\\.\\.[0-9]+", "", header_row)
sample_sheet_nn198 <- read_excel(
  path = file_nn198 ,
  sheet = "Sample & Index",
  skip = 4,
  col_names = header_row,
  # NOTE: Setting the max guess_max value avoids problems with incorrectly 
  #       guessed columns
  #       (https://github.com/tidyverse/readxl/issues/414#issuecomment-352437730)
  guess_max = 1048576)

# Restrict to specific plates
# sample_sheet_nn198 <- filter(
#   sample_sheet_nn198,
#   `Plate#` %in% c(
#     "LCE453"))

# Ensure FACS columns are stored as numeric (readxl sometimes fails, presumably 
# to weird pattern of empty cells).
sample_sheet_nn198 <- sample_sheet_nn198 %>%
  mutate_at(facs_data_idx, as.numeric)
# Tidy up names.
sample_sheet_nn198 <- bind_cols(
  clean_names(sample_sheet_nn198[, -facs_data_idx]),
  clean_names(sample_sheet_nn198[, facs_data_idx], case = "parsed"))

# Consolidate FACS markers. Originally, slightly different and/or incomplete 
# names were used for he FACS markers/fluorochromes and so Tracey used these
# unique names for each plate in the sample sheet. Subsequently, Jin provided 
# further information that allowed for the consolidation of these markers. This 
# information was manually extracted from the annotations in row 2 of the 
# 'Sample & Index' sheet.
markers <- c(
  "FSC_A", "FSC_H", "SSC_A", "V450_50_A_CD8_Vioblue", "V525_50_A_CD45_BV510",
  "B530_30_A_Epcam_FITC", "B576_26_A_CD31_140b_235a_PE", "B610_20_A_PI",
  "R780_60_A_CD3_APC_Cy7", "B710_50_A_NCAM_CD56_PERCP_vio700",
  "R660_20_A_MHC_I_APC", "V710_50_A_CD274_BV711")
markers <- setNames(markers, markers)
sample_sheet_nn198 <- bind_cols(
  clean_names(sample_sheet_nn198[, -facs_data_idx]),
  as.data.frame(
    lapply(markers, function(marker) {
      tmp <- as.matrix(
        cbind(sample_sheet_nn198[, grep(marker, colnames(sample_sheet_nn198))]))
      stopifnot(all(rowSums(!is.na(tmp)) <= 1))
      rowSums(tmp, na.rm = TRUE)
    })))

# Remove empty rows/columns.
sample_sheet_nn198 <- remove_empty(
  sample_sheet_nn198,
  which = c("rows", "cols"))

# Filter out those without a cell index sequence, with no cell, or that were 
# otherwise removed.
sample_sheet_nn198 <- sample_sheet_nn198 %>%
  filter(
    id != "removed",
    rd1_index_cell_index_index_sequence_as_in_c_rt1_primer != "removed",
    !is.na(rd1_index_cell_index_index_sequence_as_in_c_rt1_primer),
    !is.na(cell_type_sort_gate))

# Some final tidying.
sample_sheet_nn198 <- sample_sheet_nn198 %>%
  mutate(
    # NOTE: There are some wonky well_positions (e.g., 'I19=A1') that need to
    #       be fixed (these occur because it means well I19 with primer A1,
    #       in SCORE's terminology. I've asked for this to be avoided going
    #       forward.).
    well_position = gsub(" ", "", well_position),
    well_position = sapply(strsplit(well_position, "="), "[[", 1),
    well_position = factor(
      x = well_position,
      levels = unlist(
        lapply(LETTERS[1:16], function(x) paste0(x, 1:24)),
        use.names = TRUE)),
    # NOTE: `id` is `patient_id` in `donor_tbl` (below).
    sequencing_run = "NN198") %>%
  rename(patient_id = id) %>%
  arrange(plate_number, well_position)

# Construct data frame linking plates to donors --------------------------------

# NOTE: This table does not contain any of the patients from NN198. Instead, 
#       this information is included in the NN198 sample sheet.
# donor_tbl <- read_excel(
#   here("data", "sample_sheets", "Patient_Sort_Plate_overview.xlsx")) %>%
#   fill(`Patient ID`) %>%
#   mutate(plate_number = sub(" ", "", `Sort Plate ID`)) %>%
#   select(plate_number, patient_id = `Patient ID`)

# Construct final sample sheet -------------------------------------------------

### WH: skip this chunk as only deal with LCE453 of NN198 
# 
# sample_sheet <- full_join(
#   sample_sheet_nn147,
#   sample_sheet_nn157) %>%
#   left_join(donor_tbl) %>%
#   full_join(sample_sheet_nn180) %>%
#   full_join(sample_sheet_nn198) %>%
#   mutate(rowname = paste0(plate_number, "_", well_position)) %>%
#   tibble::column_to_rownames("rowname") %>%
#   DataFrame(., check.names = FALSE)
sample_sheet <- sample_sheet_nn198 %>%
  mutate(rowname = paste0(plate_number, "_", well_position)) %>%
  tibble::column_to_rownames("rowname") %>%
  DataFrame(., check.names = FALSE)

# NOTE: Check that there aren't any malformed well_positions (e.g., 'I19=A1').
stopifnot(!anyNA(sample_sheet$well_position))






# Key variables ----------------------------------------------------------------

plates <- unique(sample_sheet$plate_number)
names(plates) <- plates
sequencing_runs <- tapply(
  sample_sheet$sequencing_run,
  sample_sheet$plate_number,
  unique)

outdir <- here("data", "zUMIs")
dir.create(outdir, recursive = TRUE)

extdir <- here("extdata", sequencing_runs, "zUMIs", plates)
names(extdir) <- plates
sapply(extdir, dir.create, recursive = TRUE)

# # NOTE: Only using first 7 nt of barcode.
# read_structure <- get_read_str("CEL-Seq2")
# read_structure$bl2 <- 7
# # NOTE: Must be an element of biomaRt::listDatasets(), e.g.,
# #       biomaRt::listDatasets(biomaRt::useEnsembl("ensembl"))[["dataset"]]
# organism <- "mmusculus_gene_ensembl"
# # NOTE: Must be an element of biomaRt::listAttributes(), e.g.,
# #       biomaRt::listAttributes(biomaRt::useEnsembl("ensembl", organism))[["name"]]
# gene_id_type <- "ensembl_gene_id"

# Input files ------------------------------------------------------------------

# FASTQ files
# NOTE: Some plates have both single-cell samples as well as
#       'population/mini-bulk' samples. When both are present, there are two 
#       pairs of FASTQ files for that plate. We much combine these at the
#       plate-level before using with scPipe.
r1_fq <- c(
  # NN198
  grep(
    pattern = "Undetermined", 
    x = list.files(
      path = here("extdata", "NN198", "merged"),
      full.names = TRUE,
      pattern = glob2rx("*R1.fastq.gz")),
    invert = TRUE,
    value = TRUE))

r2_fq <- gsub("R1", "R2", r1_fq)
stopifnot(all(file.exists(r2_fq)))
tx_fq <- file.path(extdir, paste0(plates, ".R2.fastq.gz"))
names(tx_fq) <- plates
barcode_fq <- gsub("R2", "R1", tx_fq)







lapply(plates, function(plate) {
  message(plate)
  sequencing_run <- unique(
    sample_sheet$sequencing_run[sample_sheet$plate_number == plate])
  if (sequencing_run == "NN147") {
    cmd <- paste0(
      "cat ", 
      paste0(
        grep(
          plate, 
          grep(sequencing_run, r1_fq, value = TRUE), 
          value = TRUE), 
        collapse = " "),
      " > ", 
      barcode_fq[[plate]], 
      "\n",
      "cat ", 
      paste0(grep(
        plate, 
        grep(sequencing_run, r2_fq, value = TRUE), 
        value = TRUE), 
        collapse = " "),
      " > ", 
      tx_fq[[plate]])
  } else if (sequencing_run == "NN157") {
    rpi <- unique(
      sample_sheet$illumina_index_index_number_separate_index_read[sample_sheet$plate_number == plate])
    cmd <- paste0(
      "cat ", 
      paste0(
        sapply(rpi, function(x) {
          grep(
            paste0(sub(" ", "-", x), "_"),
            grep(sequencing_run, r1_fq, value = TRUE), 
            value = TRUE)
        }), 
        collapse = " "),
      " > ", 
      barcode_fq[[plate]], 
      "\n",
      "cat ", 
      paste0(
        sapply(rpi, function(x) {
          grep(
            paste0(sub(" ", "-", x), "_"),
            grep(sequencing_run, r2_fq, value = TRUE), 
            value = TRUE)
        }), 
        collapse = " "),
      " > ", 
      tx_fq[[plate]])
  } else if (sequencing_run == "NN198") {
    rpi <- unique(
      sample_sheet$illumina_index_index_number_separate_index_read[sample_sheet$plate_number == plate])
    cmd <- paste0(
      "cat ", 
      paste0(
        sapply(rpi, function(x) {
          grep(
            paste0(sub(" ", "-", x), "_"),
            grep(sequencing_run, r1_fq, value = TRUE), 
            value = TRUE)
        }), 
        collapse = " "),
      " > ", 
      barcode_fq[[plate]], 
      "\n",
      "cat ", 
      paste0(
        sapply(rpi, function(x) {
          grep(
            paste0(sub(" ", "-", x), "_"),
            grep(sequencing_run, r2_fq, value = TRUE), 
            value = TRUE)
        }), 
        collapse = " "),
      " > ", 
      tx_fq[[plate]])
  } else if (sequencing_run == "NN197") {
    rpi <- unique(
      sample_sheet$illumina_index_index_number_separate_index_read[sample_sheet$plate_number == plate])
    cmd <- paste0(
      "cat ", 
      paste0(
        sapply(rpi, function(x) {
          grep(
            paste0(sub(" ", "-", x), "_"),
            grep(sequencing_run, r1_fq, value = TRUE), 
            value = TRUE)
        }), 
        collapse = " "),
      " > ", 
      barcode_fq[[plate]], 
      "\n",
      "cat ", 
      paste0(
        sapply(rpi, function(x) {
          grep(
            paste0(sub(" ", "-", x), "_"),
            grep(sequencing_run, r2_fq, value = TRUE), 
            value = TRUE)
        }), 
        collapse = " "),
      " > ", 
      tx_fq[[plate]])
  }
  system(cmd)
})


# # Genome index
# #
# ### WH VIP: before running this line, the soft link pointing at the genome index needed to be properly set up !
# ### i.e. in directory of "extdata"
# ### ln -s /stornext/Projects/score/Indexes/
# #
# genome_index <- here("extdata", "Indexes", "GRCm38.p6", "GRCm38_with_ERCC")
# 
# # Genome annotation(s)
# annofn <- c(
#   here("extdata", "Indexes", "GRCm38.p6", "gencode.vM25.primary_assembly.annotation.gff3"),
#   system.file("extdata", "ERCC92_anno.gff3", package = "scPipe"))
# 
# Cell barcodes
bc_anno <- file.path(extdir, paste0(plates, ".barcode_annotation.csv"))
names(bc_anno) <- plates

for (plate in plates) {
  message(plate)
  tmp <- sample_sheet[sample_sheet$plate_number == plate, ]
  barcode_df <- data.frame(
    cell_id = row.names(tmp),
    # NOTE: For some reason the primer name and sequence columns have been
    #       reversed in this sample sheet.
    # NOTE: Only using first 7 nt of barcode.
    barcode = strtrim(tmp$c_rt1_primer_name, 7),
    stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(barcode_df$barcode))
  write.csv(
    x = barcode_df,
    file = bc_anno[[plate]],
    quote = FALSE,
    row.names = FALSE)
}
# 
# # Output files -----------------------------------------------------------------
# 
# combined_fq <- file.path(extdir, gsub("R[12]", "combined", basename(tx_fq)))
# names(combined_fq) <- names(tx_fq)
# subread_bam <- gsub("fastq.gz", "subread.bam", combined_fq, fixed = TRUE)
# exon_bam <- gsub("subread", "exon", subread_bam)
# 
# # FASTQ reformatting -----------------------------------------------------------
# 
# filter_settings <- list(rmlow = TRUE, rmN = FALSE, minq = 20, numbq = 2)
# # NOTE: Have to loop over files because sc_trim_barcode() is not vectorised.
# mclapply(seq_along(tx_fq), function(i) {
#   message(combined_fq[i])
#   sc_trim_barcode(
#     outfq = combined_fq[i],
#     r1 = tx_fq[i],
#     r2 = barcode_fq[i],
#     read_structure = read_structure,
#     filter_settings = filter_settings)
# })

# # Aligning reads to a reference genome -----------------------------------------
# 
# Rsubread::align(
#   index = genome_index,
#   readfile1 = combined_fq,
#   output_file = subread_bam,
#   nthreads = 8)
# 
# # Assigning reads to annotated exons -------------------------------------------
# 
# bam_tags <- list(am = "YE", ge = "GE", bc = "BC", mb = "OX")
# bc_len <- read_structure$bl1 + read_structure$bl2
# barcode_vector <- ""
# UMI_len <- read_structure$ul
# stnd <- TRUE
# fix_chr <- FALSE
# mclapply(seq_along(subread_bam), function(i) {
#   message(i)
#   sc_exon_mapping(
#     inbam = subread_bam[i],
#     outbam = exon_bam[i],
#     annofn = annofn,
#     bam_tags = bam_tags,
#     bc_len = bc_len,
#     barcode_vector = barcode_vector,
#     UMI_len = UMI_len,
#     stnd = stnd,
#     fix_chr = fix_chr)
# })
# 
# # De-multiplexing data ---------------------------------------------------------
# 
# max_mis <- 1
# has_UMI <- TRUE
# mito <- "chrM"
# mclapply(seq_along(exon_bam), function(i) {
#   message(i)
#   sc_demultiplex(
#     inbam = exon_bam[i],
#     outdir = extdir[i],
#     bc_anno = bc_anno[i],
#     max_mis = max_mis,
#     bam_tags = bam_tags,
#     mito = mito,
#     has_UMI = has_UMI)
# })
# 
# # Gene counting ----------------------------------------------------------------
# 
# UMI_cor <- 1
# gene_fl <- FALSE
# mclapply(seq_along(bc_anno), function(i) {
#   message(i)
#   sc_gene_counting(
#     outdir = extdir[i],
#     bc_anno = bc_anno[i],
#     UMI_cor = UMI_cor,
#     gene_fl = gene_fl)
# })
# 
# # Create and save SingleCellExperiment -----------------------------------------
# 
# list_of_sce <- lapply(plates, function(plate) {
#   message(plate)
#   create_sce_by_dir(
#     datadir = extdir[[plate]],
#     organism = organism,
#     gene_id_type = gene_id_type,
#     pheno_data = sample_sheet[sample_sheet$plate_number == plate, ],
#     # NOTE: Create the report separately for more fine-grained control.
#     report = FALSE)
# })
# 
# # source(here("code", "helper_functions.R"))
# # sce <- Reduce(function(x, y) .combine(x, y, rowData_by = NULL), list_of_sce)
# # assay(sce, withDimnames = FALSE) <- as(
# #   assay(sce, withDimnames = FALSE),
# #   "dgCMatrix")
# # sce <- splitAltExps(
# #   sce,
# #   ifelse(grepl("^ERCC-", rownames(sce)), "ERCC", "Endogenous"))
# # saveRDS(
# #   sce,
# #   file.path(outdir, "JinNg.xenograft.KnownHumanSample_ALLplates_Map2human.scPipe.SCE.rds"),
# #   compress = "xz")
# #
# # WH: somehow the chunk above returns an error of:
# #
# # the SCE object has to be saved individually as follows:
# #
# dss <- c("LCE437", "LCE438", "LCE439", "LCE440", "LCE441", "LCE442", "LCE443", 
#          "LCE444", "LCE446", "LCE447", "LCE448", "LCE449", "LCE453")
# 
# for(i in seq_along(dss)) {
#   ds <- dss[i]
#   single_sce <- list_of_sce[[i]]
#   
#   assay(single_sce, withDimnames = FALSE) <- as(
#     assay(single_sce, withDimnames = FALSE),
#     "dgCMatrix")
#   
#   single_sce <- splitAltExps(
#     single_sce,
#     ifelse(grepl("^ERCC-", rownames(single_sce)), "ERCC", "Endogenous"))
#   
#   saveRDS(
#     single_sce,
#     file.path(outdir, paste0("JinNg.xenograft.KnownHumanSample_Map2mouse_", ds, ".scPipe.SCE.rds")),
#     compress = "xz")
# }
# 
# # Create QC report -------------------------------------------------------------
# 
# library(readr)
# library(plotly)
# library(DT)
# library(scran)
# library(Rtsne)
# # NOTE: Needs a fix for https://github.com/LuyiTian/scPipe/issues/100.
# dir.create(here("output", "scPipe"))
# # NOTE: Tends to crap itself if using mclapply().
# lapply(plates, function(plate) {
#   try(create_report(
#     sample_name = plate,
#     outdir = extdir[[plate]],
#     r1 = tx_fq[[plate]],
#     r2 = barcode_fq[[plate]],
#     outfq = combined_fq[[plate]],
#     read_structure = read_structure,
#     filter_settings = filter_settings,
#     align_bam = subread_bam[[plate]],
#     genome_index = genome_index,
#     map_bam = exon_bam[[plate]],
#     exon_anno = annofn,
#     stnd = stnd,
#     fix_chr = fix_chr,
#     barcode_anno = bc_anno[[plate]],
#     max_mis = max_mis,
#     UMI_cor = UMI_cor,
#     gene_fl = gene_fl,
#     organism = organism,
#     gene_id_type = gene_id_type))
#   
#   # NOTE: Workaround bug in create_report() and stop output after 'Data summary'
#   #       section.
#   tmp <- readLines(file.path(extdir[[plate]], "report.Rmd"))
#   tmp <- c(tmp[1:161], "knitr::knit_exit()", tmp[162:length(tmp)])
#   writeLines(tmp, file.path(extdir[[plate]], "report.Rmd"))
#   knitr::wrap_rmd(
#     file = file.path(extdir[[plate]], "report.Rmd"), 
#     width = 120, 
#     backup = NULL)
#   rmarkdown::render(
#     input = file.path(extdir[[plate]], "report.Rmd"), 
#     output_file = file.path(extdir[[plate]], "report.html"), 
#     knit_root_dir = ".")
#   
#   # NOTE: Copy the QC report to the repository.
#   file.copy(
#     from = file.path(extdir[[plate]], "report.nb.html"),
#     to = here(
#       "output", 
#       "scPipe", 
#       paste0(plate, ".scPipe_QC_report.KnownHumanSample_Map2mouse.nb.html")),
#     overwrite = TRUE)
# })