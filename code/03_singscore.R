# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "dplyr", "stringr", "yaml", "naturalsort",  "reshape2", "survival",
  "survminer", "data.table", "parallel","meta", "openxlsx",
  "foreach","doMC", "doParallel", "ggrepel", "data.table"
)
.bioc_packages = c(
  "Biobase", "singscore", "GSEABase", "genefilter", "SummarizedExperiment", "msigdbr"
)

# Install CRAN packages (if not already installed)
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

if (any(!"cancersea" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  devtools::install_github("https://github.com/camlab-bioml/cancersea")
}
library(cancersea)
data('available_pathways')

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
source("code/helper/singscore_helper.R")

standardize = function(z) {
  rowmean = apply(z, 1, mean, na.rm=TRUE)
  rowsd = apply(z, 1, sd, na.rm=TRUE)
  rv = sweep(z, 1, rowmean,"-")
  rv = sweep(rv, 1, rowsd, "/")
  return(rv)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Gene sets
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MSigDB Hallmakrs
gene.sets = msigdbr(species = "Homo sapiens", category = "H") %>% data.frame()
gene.sets.l = list()
for (i in unique(gene.sets$gs_name)) {
  # gene.set = GSEABase::GeneSet(gene.sets[gene.sets$gs_name == i, ]$human_ensembl_gene)
  gene.set = GSEABase::GeneSet(unique(gene.sets[gene.sets$gs_name == i, ]$human_gene_symbol))
  setName(gene.set) = gsub("HALLMARK_", "", i)
  gene.sets.l[[i]] = gene.set
}
msigdb.hm = split_gene_sets(gene.sets.l) # split in up/down gene sets (if available)

# CancerSEA
gene.sets.l = list()
for (i in available_pathways) {
  # gene.set = GSEABase::GeneSet( eval(as.name(i))$ensembl_gene_id )
  gene.set = GSEABase::GeneSet( unique(eval(as.name(i))$symbol) )
  setName(gene.set) = i
  gene.sets.l[[i]] = gene.set
}
cancersea.hm = split_gene_sets(gene.sets.l) # split in up/down gene sets (if available)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Prepare TCGA cohorts
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
TCGA.se = readRDS("data/tcga_gene_gdc.Rds")
TCGA.se = TCGA.se[, TCGA.se$gdc_cases.samples.sample_type %in% c("Primary Tumor")]

g.s.ftrs.1 = unique(unname(unlist(lapply(cancersea.hm$unpaired, function(i){i@geneIds}))))
g.s.ftrs.2 = unique(unname(unlist(lapply(msigdb.hm$unpaired, function(i){i@geneIds}))))
g.s.ftrs = unique(c(g.s.ftrs.1, g.s.ftrs.2))

tcga.prep = list()
for(i in 9:9){
# for(i in 1:length(unique(TCGA.se$gdc_cases.project.project_id))){

  TCGAType = gsub(".+-", "", unique(TCGA.se$gdc_cases.project.project_id)[i])
  print(TCGAType)

  se = TCGA.se[, TCGA.se$gdc_cases.project.project_id %in% unique(TCGA.se$gdc_cases.project.project_id)[i]]
  pheno = data.frame(colData(se))
  organ = unique(se$gdc_cases.project.primary_site)

  # remove column contain only NAs
  pheno =  pheno[, which(apply(pheno, 2, function(x) !all(is.na(x))))]

  exprs.mat = log(assays(se)$fpkm + 1)
  exprs.mat = genefilter::varFilter(as.matrix(exprs.mat), var.func = IQR,  var.cutoff = 0.25, filterByQuantile = TRUE)

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
  # table(pheno$gdc_cases.submitter_id[pheno$gdc_cases.submitter_id %in% dup] )

  sampleIDs = rownames(pheno)
  exprs.mat = exprs.mat[, colnames(exprs.mat) %in% sampleIDs]
  # Consistency check: Expression colnames are in the same order as pheno rownames
  stopifnot(identical(rownames(pheno), colnames(exprs.mat)))

  # Survival object
  survivalDOD = Surv(pheno$DDS_TIME, pheno$DDS_SURV)

  rankData = singscore::rankGenes(exprs.mat)

  pd = data.frame(pheno)
  pd = pd %>% dplyr::select(
    ORGAN = gdc_cases.project.primary_site, TYPEe = gdc_cases.samples.sample_type,
    COHORT = gdc_cases.project.project_id, CEP = DDS_SURV
  )
  pd$COHORT = gsub("TCGA-", "", pd$COHORT)
  pd = pd[rownames(pd) %in% colnames(rankData), ]
  pd = pd[colnames(rankData), ]
  stopifnot(identical(rownames(pd), colnames(rankData)))

  se = SummarizedExperiment(assays = list("log_fpkm" = exprs.mat), colData = pheno)
  se = se[rownames(se) %in% g.s.ftrs, ]

  tcga.prep[[TCGAType]] = list(
    rankData = rankData,
    p.data = pd,
    survivalDOD = survivalDOD,
    se.obj = se
  )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Function: SingScore (ssGSEA), LogRank & Cox with enrichment scores
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ssgsea_surv = function(
  .cohorts.prep,
  .gene.sets.l,
  .permut = F
){

  ss.surv.res.l <- list()
  ss.res.l = list()
  ss.res.unpaired.l = list()
  ss.res.paired.l = list()

  for(i in 1:length(.cohorts.prep)){

    study.id <- names(.cohorts.prep)[i]
    print(study.id)

    survivalDOD = .cohorts.prep[[study.id]]$survivalDOD
    rankData = .cohorts.prep[[study.id]]$rankData

    ss.res.unpaired = ssgsea(
      g.s.l = .gene.sets.l$unpaired,
      .rankData = rankData,
      paired = F,
      permut = .permut,
      alpha = 0.05,
      nbr.perm = 1000,
      ncores = 30,
      study.id = study.id
    )

    ss.res.paired = ssgsea(
      g.s.l = .gene.sets.l$paired,
      .rankData = rankData,
      paired = T,
      permut = .permut,
      alpha = 0.05,
      nbr.perm = 1000,
      ncores = 30,
      study.id = study.id
    )

    ss.sum = rbind(ss.res.paired$ss.sum, ss.res.unpaired$ss.sum)

    ss.res = rbind(ss.res.paired$ss.res, ss.res.unpaired$ss.res)
    var.non.zero = apply(ss.res, 1, function(x) { stats::var(x) > 0 })
    ss.res = ss.res[var.non.zero, ]
    ss.res.scaled = standardize(ss.res)

    ss.res.dicho = ss.res.scaled
    for (i in 1:nrow(ss.res.dicho)) {
      ss.res.dicho[i, ] = ifelse(ss.res.dicho[i, ] > median(ss.res.dicho[i, ]), 2, 1)
      # ss.res.dicho[i, ] = ifelse(ss.res.dicho[i, ] > 0, 2, 1)
    }

    # LogRank test
    logrank.res = runKaplanTab(dat = ss.res.dicho, surv = survivalDOD)[1, ]

    # Cox-Regr.
    surv.res = run_cox_gsea(ss.res.scaled, .surv = survivalDOD)
    colnames(surv.res) = paste0("COX_", colnames(surv.res))
    surv.res$GENE_SET = rownames(surv.res)
    surv.res$LogRank_Pval = logrank.res[match(surv.res$GENE_SET, names(logrank.res))]
    surv.res = cbind(surv.res, ss.sum[rownames(surv.res), ] %>% dplyr::select(-GENE_SET))

    ss.surv.res.l[[study.id]] = surv.res
    ss.res.l[[study.id]] = ss.res
    ss.res.unpaired.l[[study.id]] = ss.res.unpaired
    ss.res.paired.l[[study.id]] = ss.res.paired

    surv.res = NULL
  }
  list(
    ss.surv.res.l = ss.surv.res.l,
    ss.res.l = ss.res.l,
    ss.res.unpaired.l = ss.res.unpaired.l,
    ss.res.paired.l = ss.res.paired.l
  )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Run wit TCGA cohorts and save
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
do.permutation = F
ss.res.msigdb.hm = ssgsea_surv(.cohorts.prep = tcga.prep, .gene.sets.l = msigdb.hm, .permut = do.permutation)
ss.res.cancersea.hm = ssgsea_surv(.cohorts.prep = tcga.prep, .gene.sets.l = cancersea.hm, .permut = do.permutation)

df = rbind(
  do.call("rbind", ss.res.msigdb.hm$ss.surv.res.l) %>% mutate(PANEL = "MSigDB Hallmarks" ),
  do.call("rbind", ss.res.cancersea.hm$ss.surv.res.l) %>% mutate(PANEL = "CancerSEA" )
)
saveRDS(df, file = "data/singscore_tcga_surv.Rds")

l = list(
  "MSigDB Hallmarks" = ss.res.msigdb.hm$ss.res.l,
  "CancerSEA" = ss.res.cancersea.hm$ss.res.l
)
saveRDS(l, file = "data/singscore_tcga_scores.Rds")

l = list(
  "MSigDB Hallmarks" = ss.res.msigdb.hm,
  "CancerSEA" = ss.res.cancersea.hm
)
saveRDS(l, file = "data/singscore_tcga_full.Rds")

l = list(
  "msigdb.hm" = msigdb.hm,
  "cancersea.hm" = cancersea.hm
)
saveRDS(l, file = "data/singscore_gene_sets.Rds")

saveRDS(tcga.prep$KIRC$se.obj, file = "data/singscore_se_objects.Rds")

