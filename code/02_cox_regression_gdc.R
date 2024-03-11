# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries and helper functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "yaml", "ggplot2", "survival", "dplyr", "survminer", "tidyverse", "naturalsort"
)
.bioc_packages = c("SummarizedExperiment")

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

standardize = function(z) {
  rowmean = apply(z, 1, mean, na.rm=TRUE)
  rowsd = apply(z, 1, sd, na.rm=TRUE)
  rv = sweep(z, 1, rowmean,"-")
  rv = sweep(rv, 1, rowsd, "/")
  return(rv)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
tcga.ror = readRDS("data/tcga_gene_ror2_gdc.Rds")
# table(tcga.ror$gdc_cases.project.primary_site, tcga.ror$gdc_cases.samples.sample_type)
# table(tcga.ror$gdc_cases.project.project_id, tcga.ror$gdc_cases.samples.sample_type)

tcga.ror = tcga.ror[, tcga.ror$gdc_cases.samples.sample_type %in% c("Primary Tumor")]
colData(tcga.ror) = droplevels(colData(tcga.ror))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cox run
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# file = "TCGA-KIRC"
# file = "TCGA-PRAD"
cox.res = parallel::mclapply(unique(tcga.ror$gdc_cases.project.project_id), function(file) {

  TCGAType = gsub(".+-", "", file)
  print(TCGAType)

  se = tcga.ror[, tcga.ror$gdc_cases.project.project_id %in% file]
  pheno = data.frame(colData(se))

  organ = unique(se$gdc_cases.project.primary_site)

  # remove column contain only NAs
  pheno =  pheno[, which(apply(pheno, 2, function(x) !all(is.na(x))))]

  exprs.mat = data.frame(assays(se)$fpkm, check.names = F)
  exprs.mat = log(exprs.mat + 1)
  exprs.mat.tpm = exprs.mat
  exprs.mat = standardize(exprs.mat)

  # available sample types (normal, primary, mts, etc)
  print(table(pheno$gdc_cases.samples.sample_type))

  # sample with frozen and ffpe could be available, delete ffpe
  idx = which(pheno$gdc_cases.samples.is_ffpe == TRUE)
  if (length(idx) != 0){
    pheno = pheno[-idx,]
  }
  # looking for cases with duplicated id after ffpe deleted
  dup = unique(pheno$gdc_cases.submitter_id[duplicated(pheno$gdc_cases.submitter_id)])
  table(pheno$gdc_cases.submitter_id[pheno$gdc_cases.submitter_id %in% dup] )

  # Some cases have more than 1 primary tumor, choose the first occurrence
  # Alternative: tcga_replicateFilter.R
  for (duplicates in dup) {
    x = which(pheno$gdc_cases.submitter_id %in% duplicates)
    toDrop = x[2:length(x)]
    pheno = pheno[-toDrop, ]
  }
  table(pheno$gdc_cases.submitter_id[pheno$gdc_cases.submitter_id %in% dup] )

  sampleIDs = rownames(pheno)
  exprs.mat = exprs.mat[, colnames(exprs.mat) %in% sampleIDs]
  # Consistency check: Expression colnames are in the same order as pheno rownames
  stopifnot(identical(rownames(pheno), colnames(exprs.mat)))

  # Survival object
  survivalDOD = Surv(pheno$DDS_TIME,pheno$DDS_SURV)

  survModel = survfit(survivalDOD ~ 1, data = pheno)

  run_cox = function(obj) {

    exprs.mat.spl = split(as.matrix(obj), rownames(obj))

    coxModel.l = parallel::mclapply(exprs.mat.spl, function(gene) {
      res = summary(survival::coxph(survivalDOD ~ gene))
      log.rank =  tryCatch(
        survival::survdiff(survivalDOD ~ ifelse(gene > 0, 2, 1))$pvalue, error=function(e) "error"
      )
      if(log.rank == "error") {
        log.rank = 1
      }
      res = c(
        "HR" = res$coef[2],
        "logHR" = res$coef[1],
        "SE_logHR" = res$coef[3],
        "L95CI" = res$conf.int[,"lower .95"],
        "U95CI" = res$conf.int[,"upper .95"],
        "Rsquare" = res$rsq[[1]],
        "LR_Pval" = unname(res$logtest[3]),
        "LogRank_Pval" = log.rank,
        "Pval" = res$coef[5]
      )
    }, mc.cores = 2)

    coxModel = data.frame(do.call("rbind", coxModel.l))
    coxModel$Padj = stats::p.adjust(coxModel$Pval, method = "BH")
    coxModel
  }

  # Cox survival model fitting
  coxModel = run_cox(exprs.mat)
  coxModel = coxModel[!is.na(coxModel$Padj), ]
  coxModel$ORGAN = organ
  coxModel$TYPE = TCGAType
  map_signif_level = c("***" = 0.001, "**" = 0.01,  "*" = 0.05, "ns" = 1.1)

  signif.lvl = c()
  for (x in 1:length(coxModel$Pval)) {
    signif.lvl[x] = names(which.min(map_signif_level[which(map_signif_level > coxModel$Pval[x])]))
  }

  coxModel$Pval_lvl = signif.lvl

  exprs.mat.tpm = exprs.mat.tpm[rownames(coxModel), colnames(exprs.mat)]
  coxModel$AVE_EXPRS = rowMeans(exprs.mat.tpm)

  # Log Rank
  pheno$ROR2_EXPRS = as.numeric(exprs.mat[1, ])
  pheno$ROR2_GROUP = ifelse(pheno$ROR2_EXPRS > 0, 1, 0)
  coxModel$LogRank_Pval = survival::survdiff(survivalDOD ~ ROR2_GROUP, data = pheno)$pvalue

  l = list()
  l[[TCGAType]] = list(coxModel = coxModel, survModel = survModel)
  l
}, mc.cores = 20)

coxModelByTypes = list()
survModelByTypes = list()
for (i in 1:length(cox.res)) {
  coxModelByTypes[[names(cox.res[[i]])]] = cox.res[[i]][[1]]$coxModel
  survModelByTypes[[names(cox.res[[i]])]] = cox.res[[i]][[1]]$survModel
}

saveRDS(coxModelByTypes, file =paste0("data/tcga_cox_ror2.Rds"))
saveRDS(survModelByTypes, file =paste0("data/tcga_surv_ror2.Rds"))
