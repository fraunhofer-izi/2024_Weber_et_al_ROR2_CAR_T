# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("devtools", "yaml", "dplyr")
.bioc_packages = c("recount", "GenomicRanges", "GenomicRanges")

## Install CRAN packages (if not already installed)
.inst = .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

# Install bioconductor packages (if not already installed)
.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  library(BiocManager)
  BiocManager::install(.bioc_packages[!.inst], ask = T)
}

list.of.packages = c(.cran_packages, .bioc_packages)

## Loading library
for (pack in list.of.packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# GDC Expression matrix
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# https://gdc.cancer.gov/about-data/publications/pancanatlas
gdc.tcga = readr::read_tsv(
  "data/RNASeqV2.geneExp.tsv"
) %>% data.frame()
gdc.tcga[is.na(gdc.tcga)] = 0
rownames(gdc.tcga) = gdc.tcga[, 1]
gdc.tcga[, 1] = NULL
gdc.tcga = data.frame(gdc.tcga)
colnames(gdc.tcga) = gsub("\\.", "-", colnames(gdc.tcga))

pdata.ids = openxlsx::read.xlsx(
  "data/TCGA-CDR-SupplementalTableS1.xlsx",
  sheet = 1, rowNames = F, startRow = 1, detectDates = T
)

mat = gdc.tcga
rc.clin.data = readRDS("data/tcga_recount_phenodata_gdc.Rds")
rc.clin.data = rc.clin.data[!duplicated(rc.clin.data$gdc_cases.samples.portions.analytes.aliquots.submitter_id), ]

shared.smpls = intersect(colnames(mat), rc.clin.data$gdc_cases.samples.portions.analytes.aliquots.submitter_id)
mat = mat[, colnames(mat) %in% shared.smpls]
rc.clin.data = rc.clin.data[rc.clin.data$gdc_cases.samples.portions.analytes.aliquots.submitter_id %in% shared.smpls, ]
mat = mat[, rc.clin.data$gdc_cases.samples.portions.analytes.aliquots.submitter_id]
stopifnot(identical(
  colnames(mat),
  rc.clin.data$gdc_cases.samples.portions.analytes.aliquots.submitter_id
))
mat = mat[!grepl("\\?", rownames(mat)), ]
mat = mat[!duplicated(gsub("\\|.+", "", rownames(mat))), ]
rownames(mat) = gsub("\\|.+", "", rownames(mat))
colnames(mat) = rownames(rc.clin.data)

pheno = rc.clin.data
stopifnot(identical(rownames(pheno), colnames(mat)))

pdata.ids.2 = pdata.ids[pdata.ids$bcr_patient_barcode %in% pheno$xml_bcr_patient_barcode, ]
pheno$DDS_TIME = pdata.ids.2$DSS.time[match(pheno$xml_bcr_patient_barcode, pdata.ids.2$bcr_patient_barcode)]
pheno$DDS_SURV = pdata.ids.2$DSS[match(pheno$xml_bcr_patient_barcode, pdata.ids.2$bcr_patient_barcode)]

se.work = SummarizedExperiment(assays = list("fpkm" = mat), colData = pheno)

saveRDS(se.work, file = "data/tcga_gene_gdc.Rds")
saveRDS(se.work["ROR2", ], file = "data/tcga_gene_ror2_gdc.Rds")
