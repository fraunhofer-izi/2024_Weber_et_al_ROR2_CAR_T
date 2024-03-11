# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries and helper functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "yaml", "ggplot2", "survival", "dplyr", "survminer", "tidyverse", "naturalsort",
  "cowplot", "scico", "forestmodel", "ggcorrplot", "msigdbr", "grid", "gridExtra", "scales"
)
.bioc_packages = c("SummarizedExperiment", "GSEABase")

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

source("code/helper/plotStyles.R")

if (any(!"cancersea" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  devtools::install_github("https://github.com/camlab-bioml/cancersea")
}
library(cancersea)
data('available_pathways')

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ROR2 (Uni Cox)
coxModelByTypes = readRDS("data/tcga_cox_ror2.Rds")
survModelByTypes = readRDS("data/tcga_surv_ror2.Rds")

tcga.anno = read.csv2("data/tcga_cancer_types_PMID_29625050.csv")

se.ror2 = readRDS("data/tcga_gene_ror2_gdc.Rds")

# Singscore
ss.scores = readRDS("data/singscore_tcga_scores.Rds")
ss.surv = readRDS("data/singscore_tcga_surv.Rds")
ss.gs = readRDS("data/singscore_gene_sets.Rds")
ss.se = readRDS("data/singscore_se_objects.Rds")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Filter samples
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
filter_obj = function(obj){
  obj = obj[, obj$gdc_cases.samples.sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]
  colData(obj) = droplevels(colData(obj))

  #  # select only normals with > 10 samples
  smpl.fiter = table(obj$gdc_cases.project.project_id, obj$gdc_cases.samples.sample_type)[, "Solid Tissue Normal"] < 10
  pd.tmp = colData(obj)[obj$gdc_cases.samples.sample_type == "Solid Tissue Normal", ]
  smpl.remove = rownames(pd.tmp[pd.tmp$gdc_cases.project.project_id %in% names(smpl.fiter[smpl.fiter]), ])
  obj = obj[, !colnames(obj) %in% smpl.remove]

  # sample with frozen and ffpe could be available, delete ffpe
  idx = which(obj$gdc_cases.samples.is_ffpe == TRUE)
  if (length(idx) != 0){
    obj = obj[, !obj$gdc_cases.samples.is_ffpe == TRUE]
  }

  # looking for cases with duplicated id after ffpe deleted
  dup = unique(obj$gdc_cases.submitter_id[duplicated(obj$gdc_cases.submitter_id)])
  table(obj$gdc_cases.submitter_id[obj$gdc_cases.submitter_id %in% dup] )

  obj
}
se.ror2 = filter_obj(obj = se.ror2)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
lvls = c(
  "CNS", "Eye", "Head/\nNeck", "Endocrine", "Thym.", "Thoracic", "Breast", "Core GI",
  "Dev. GI tract", "Genitourinary", "Gynecologic", "Skin", "Soft\nTissue",
  "Heme"
)
tcga.anno$DISEASE_TYPE[tcga.anno$DISEASE_TYPE=="Thymus"] <-"Thym."
tcga.anno$DISEASE_TYPE[tcga.anno$DISEASE_TYPE=="Head & Neck"] <-"Head/\nNeck"
tcga.anno$DISEASE_TYPE[tcga.anno$DISEASE_TYPE=="Soft Tissue"] <-"Soft\nTissue"
tcga.anno$DISEASE_TYPE[tcga.anno$DISEASE_TYPE=="Developmental GI tract"]<-"Dev. GI tract"
tcga.anno$DISEASE_TYPE = factor(tcga.anno$DISEASE_TYPE, levels = lvls)
my.cols = tcga.anno$DISEASE_TYPE_COLOUR[match(lvls, tcga.anno$DISEASE_TYPE)]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Figure 3A: TCGA: COX/Log-Rank test for ROR2 accross TCGA cohorts
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
obj = lapply(coxModelByTypes, function(x){
  x$TX = rownames(x)
  x
})
obj = do.call("rbind", obj)
obj$Pval
cox.alpha = 0.05
logrank.alpha = 0.05
pl.title = NULL

cox.leg.1 = paste0("HR >1 (p <", cox.alpha, ")")
cox.leg.2 = paste0("HR <=1 (p <", cox.alpha, ")")
cox.leg.3 = paste0("HR (p >=", cox.alpha, ")")
obj$HR_DIR = ifelse(obj$HR > 1, cox.leg.1, cox.leg.2)
obj$HR_DIR = ifelse(obj$Pval < cox.alpha, obj$HR_DIR, cox.leg.3)
obj$HR_DIR = factor(obj$HR_DIR, levels = c(cox.leg.1, cox.leg.2, cox.leg.3) )

obj$LogRank_Pval_sign = obj$LogRank_Pval < logrank.alpha
obj$X_LABEL = obj$TYPE
obj$ORGAN = tcga.anno$DISEASE_TYPE[match(obj$TYPE, tcga.anno$DISEASE_CODE)]

obj$HR_DIR = factor(obj$HR_DIR, levels = c(cox.leg.1, cox.leg.2, cox.leg.3) )

pl =
  ggplot(obj, aes(x=X_LABEL, y=TX), fill = TYPE) +
  scale_x_discrete() +
  scale_y_discrete() +
  geom_point(data = subset(obj, Padj <= 1), aes(size = HR, colour = HR_DIR )) +
  geom_point(aes(shape = ifelse(LogRank_Pval_sign == FALSE, NA, LogRank_Pval_sign), size = HR, stroke = 1.5), color = "#555555") +
  scale_size(range = c(1, 8)) +
  scale_shape_manual(values = c(1), name = paste0("Log-rank (p <", logrank.alpha, ")"), labels = c("yes", "")) +
  scale_colour_manual(values = c(
    "HR >1 (p <0.05)" = "#CC6677", "HR <=1 (p <0.05)" = "#6497bf", "HR (p >=0.05)" = "#DDDDDD"),
    drop = F
  ) +
  xlab(label = NULL) + ylab(label = NULL) +
  labs(size = "Hazard Ratio (HR)", col = "Cox-Regression") +
  mytheme_grid(base_size = 12) +
  theme(
    strip.text.x = element_text(size = rel(.8), colour = "white"),
    panel.spacing.y=unit(.75, "lines"),
    panel.spacing.x=unit(.25, "lines"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "bottom",
    legend.text = element_text(margin = margin(l = -15, unit = "pt"), size = rel(1)),
    legend.spacing.x = unit(20, 'pt')
    # plot.margin = margin(r=-20, l = 1, t = 1)
  ) +
  guides(
    size = guide_legend(order=1, title.position = "top", ncol = 2),
    color = guide_legend(order=2, title.position = "top", ncol = 1, override.aes = list(size = 3)),
    shape = guide_legend(order=3, title.position = "top", ncol = 1, override.aes = list(size = 3))
  ) +
  facet_grid( ~ ORGAN, space = "free", scales = "free")

g = ggplot_gtable(ggplot_build(pl))
strip_both = which(grepl('strip-', g$layout$name))
fills = tcga.anno$DISEASE_TYPE_COLOUR

k = 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- my.cols[k]
  k <- k+1
}

pan.cox.log = g

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig 3C & Fig S4B: Boxplot: Association between standardized ROR2 expression with
# tumor stage and histological grade.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.kidney = se.ror2[, se.ror2$gdc_cases.project.primary_site == "Kidney"]
se.kidney = se.kidney[, se.kidney$gdc_cases.samples.sample_type %in% c("Primary Tumor")]

colData(se.kidney) = droplevels(colData(se.kidney))
p.data = data.frame(colData(se.kidney))
p.data$ID = rownames(p.data)
p.data = p.data[, colSums(is.na(p.data)) != nrow(p.data)]
p.data = p.data[ , !unlist(lapply(1:ncol(p.data), function(x){is.list(p.data[, x])}))]

p.data = p.data %>% mutate(
  cgc_case_pathologic_t_groups = case_when(
    grepl("^T1|^T2", cgc_case_pathologic_t) ~ "Grp 1-2",
    grepl("^T3|^T4", cgc_case_pathologic_t) ~ "Grp 3-4",
    TRUE ~ cgc_case_pathologic_t
  )
)

expr = SummarizedExperiment::assays(se.kidney)$fpkm[1, ,drop = FALSE]
expr = log10(expr + 1)
expr = expr[, colnames(expr) %in% rownames(p.data), drop = F]

df.m = reshape2::melt(as.matrix(expr))
df.m$ID = df.m$Var2
df.m = merge(df.m, p.data, by = "ID")

cp_parameter = function(df, cp, x.title = NULL, max.y = NULL){

  df$target = df[[cp]]
  df = df[!is.na(df$target), ]
  df = subset(df, target != "not reported")
  df = subset(df, target != "GX")
  df = subset(df, target != "TX")
  df = subset(df, target != "NX")
  df$target = gsub("stage ", "", df$target)

  supp.labs <- c("ChRCC", "ccRCC", "pRCC")
  names(supp.labs) <- c("TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP")

  if(!is.null(max.y)){
    min.max = c(min(df$value) - 0.001, max.y)
  } else {
    min.max = c(min(df$value) - 0.001, max(df$value))
  }

pl =
  ggplot(df, aes(target, value)) +
    geom_boxplot(outlier.colour = NA, fill = "white") +
    geom_point(size = .25, alpha = 1, position = position_jitter(seed = 1, width = .2), color = "#555555") +
    xlab(x.title) +
    ylim(min.max) +
    ylab("log10(FPKM-UQ + 1)") +
    mytheme(base_size = 12) +
    theme(
      aspect.ratio = 1,
      panel.spacing = unit(1, "lines"),
      strip.text.x = element_text(size = rel(1), face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.y = element_text(vjust = + 3),
      axis.title.x = element_text(vjust = -0.75)
    ) +
    facet_wrap(
      ~ gdc_cases.project.project_id, scales = "free_x",
      labeller = labeller(gdc_cases.project.project_id = supp.labs)
    )
  if(length(unique(df$target)) == 2) {
    pl = pl + stat_compare_means(label.x = 1.2, vjust = 1, size = 3)
  } else {
    pl = pl + stat_compare_means(vjust = 1, label.x = 1.5, size = 4)
  }
  pl
}

# cp_parameter(df = df.m, cp = "gdc_cases.diagnoses.tumor_stage", x.title = "Tumor stage")
# cp_parameter(df.m, "cgc_case_pathologic_t_groups", x.title = "Pathological T stage")

bp.kirc.grade = cp_parameter(df.m, "xml_neoplasm_histologic_grade", x.title = "Histological grade", max.y = 4) +
  theme( strip.text.x = element_blank())

bp.kirc.stage = cp_parameter(df = df.m, cp = "gdc_cases.diagnoses.tumor_stage", x.title = "Tumor stage", max.y = 4)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig 3B Kaplan Meier curves
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.kidney = se.ror2[, se.ror2$gdc_cases.project.project_id == "TCGA-KIRC"]
se.kidney = se.kidney[, se.kidney$gdc_cases.samples.sample_type %in% c("Primary Tumor")]

pheno = data.frame(colData(se.kidney))
exprs.mat = data.frame(assays(se.kidney)$fpkm, check.names = F)
exprs.mat = log(exprs.mat + 1)
exprs.mat = standardize(exprs.mat)

# Define sample groups and make Kaplan meier curves
group.cutoff = quantile(as.numeric(exprs.mat[1, ]), c(.25, .5, .75))
pheno$ROR2_EXPRS = as.numeric(exprs.mat[1, ])
pheno = pheno %>% dplyr::mutate(
  ROR2_GROUP = dplyr::case_when(
    ROR2_EXPRS <= group.cutoff[1] ~ "low",
    (ROR2_EXPRS > group.cutoff[1]) & (ROR2_EXPRS < group.cutoff[3]) ~ "med",
    ROR2_EXPRS >= group.cutoff[3] ~ "high",
    TRUE ~ "something wrong"
  )
)

pheno$DDS_TIME = pheno$DDS_TIME / 365.25

x.breaks = 2
do.trunc = T
trunc.cut = 10

survivalDOD = Surv(pheno$DDS_TIME, pheno$DDS_SURV)
survModel = survfit(survivalDOD ~ ROR2_GROUP, data = pheno)

# ungrouped KM estimate to determine cutoff
survModel.ungrouped = survminer::surv_fit(survivalDOD ~ 1, data = pheno)

if (do.trunc == T) {
  cutoff  = min(survModel.ungrouped$time[survModel.ungrouped$n.risk <= trunc.cut])
}

cust.th = theme(
  plot.title = element_text(hjust = 0.5, face = "plain", size = rel(1)),
  panel.border = element_blank(),
  axis.line = element_line(colour = "black", size = .3),
  axis.title.y = element_text(vjust = + 3),
  axis.title.x = element_text(vjust = -0.75),
  plot.margin = margin(5.5, 5.5, .6, 8, unit = "pt")
)

surv.pl =
  ggsurvplot(
    survModel,
    conf.int = T,
    break.time.by = x.breaks,
    legend = "none",
    xlab = "Time (years)",
    ylab = "DSS",
    censor.size = 1.5,
    censor.shape = 124,
    size = 0.6,
    ggtheme = mytheme(base_size = 12) + cust.th,
    tables.theme = theme_cleantable(),
    risk.table.fontsize = 3.25,
    risk.table.title = "",
    risk.table.pos = "in",
    risk.table = "absolute",
    palette = c("#997700", "#009988", "#6699CC"),
    title = NULL
  )

if (do.trunc == T) {
  surv.pl$plot = surv.pl$plot +
    coord_cartesian(xlim=c(0, cutoff)) +
    scale_x_continuous(limits = c(0, cutoff), breaks = seq(0, cutoff, by = x.breaks))
}

surv.pl$table = surv.pl$table +
  theme(
    panel.border = element_blank(),
    legend.position = "none",
    axis.text.y = element_blank(),
    plot.margin =  margin(0, 0, 0, 0, unit = "pt")
  )

if (do.trunc == T) {
  surv.pl$table = surv.pl$table +
    coord_cartesian(xlim=c(0, cutoff)) +
    scale_x_continuous(limits = c(0, cutoff), breaks = seq(0, cutoff, by = x.breaks))
}

surv.pl = surv.pl$plot + patchwork::inset_element(surv.pl$table, left = 0, bottom = 0, right = 1, top =.28)

pheno.sub = subset(pheno, ROR2_GROUP != "med")
pheno.sub = subset(pheno, ROR2_GROUP != "high")
pheno.sub = subset(pheno, ROR2_GROUP != "low")
survival::survdiff(
  Surv(pheno.sub$DDS_TIME, pheno.sub$DDS_SURV) ~ ROR2_GROUP,
  data = pheno.sub
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig 3D Singscore
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ss.mat = ss.scores$CancerSEA$KIRC
# ss.mat = ss.scores$`MSigDB Hallmarks`$KIRC
ss.mat = ss.mat[rowMax(ss.mat) > .1, ]

se = se.ror2[, colnames(se.ror2) %in% colnames(ss.mat)]
ss.mat = ss.mat[, colnames(ss.mat) %in% colnames(se.ror2)]
stopifnot(identical(colnames(se), colnames(ss.mat)))
ss.mat = ss.mat[, colnames(se)]

res.cor = cor(
  t(log10(assays(se)$fpkm + 1)),
  t(ss.mat), method = "spearman"
)

res.cor = data.frame(t(res.cor))
res.cor = res.cor[abs(res.cor$ROR2) > 0.05, ,drop = F]
res.cor = res.cor[order(res.cor[, 1]), , drop = F]
colnames(res.cor) = "COR"
res.cor$PATHWAY = rownames(res.cor)
res.cor$PATHWAY = factor(res.cor$PATHWAY, levels = res.cor$PATHWAY)
res.cor$HR = ss.surv$COX_HR[match(res.cor$PATHWAY, ss.surv$GENE_SET)]
res.cor$COX_Pval = ss.surv$COX_Pval[match(res.cor$PATHWAY, ss.surv$GENE_SET)]
res.cor$COX_SIGN = ifelse(res.cor$COX_Pval < 0.05, "Y", NA)
max.cor = max(abs(res.cor$HR))

pl =
  ggplot(res.cor, aes(PATHWAY, COR, fill = HR, shape = COX_SIGN)) +
  coord_flip() +
  geom_bar(stat="identity") +
  geom_point(
    aes(y = ifelse(COR > 0, COR + .02, COR - .05)),
    position=position_dodge(0.9),
    size = 2, color = "#555555", show.legend=T
  ) +
  scale_shape_manual(values = c(8), name = paste0("Cox-Regression"), labels = c("p < 0.05", "")) +
  mytheme(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(vjust = -.5),
    plot.title = element_text(hjust = 0, face = "bold", colour = "black", size = rel(1)),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = .3)
  ) +
  guides(
    fill = guide_colorbar(
      barwidth = unit(.4,'lines'), order = 1,
      barheight = unit(5, 'lines'), ticks.linewidth = 1.5/.pt
    )
  ) +
  ylab("Spearman correlation between\nROR2 expression and gene set score") +
  xlab(NULL) + labs(fill = "Hazard Ratio") +
  ggtitle(NULL) +
  scale_fill_scico(palette = "vik", midpoint = 1, breaks = scales::breaks_pretty(4))
  # scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-max.cor, max.cor))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Figure 3
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ggsave2(
  cowplot::plot_grid(
    NULL,
    pan.cox.log,
    NULL,
    plot_grid(
      surv.pl,
      NULL,
      bp.kirc.grade +
        theme(
          aspect.ratio = NULL,
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", size = .3)
        ),
      NULL,
      pl + theme(plot.margin = margin(5.5, 5.5, -7, 5.5, unit = "pt")),
      ncol = 5, rel_widths = c(1.05, .15, .93, .15, 2.07),
      labels = c("B", "", "C", "", "D"), label_size = 22, label_y = 1.1
    ),
    NULL,
    nrow = 5, rel_heights = c(.1, .9, .15, 1, .05),
    labels = c("A", "", "", ""), label_size = 22, label_y = 1.1
  ),
  filename = "analysis/publication/Main_Fig_3.png",
  dpi = 300, width = 180, height = 85, units = "mm", bg = "white", scale = 2
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig S4A: Gene expression Boxplot TCGA (primary, normal)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
tcga_epxrs = function(obj, qry) {

  obj.ftr = obj[qry, ]

  expr = SummarizedExperiment::assays(obj.ftr)$fpkm[1, ]
  expr = as.numeric(expr)
  site = obj.ftr$gdc_cases.project.primary_site
  type = obj.ftr$gdc_cases.samples.sample_type
  tumor_perc = obj.ftr$cgc_slide_percent_tumor_cells
  id = gsub("TCGA-", "", obj.ftr$gdc_cases.project.project_id)

  df = data.frame(expr, site, type, id, tumor_perc)
  df$type_coarse = tcga.anno$DISEASE_TYPE[match(df$id, tcga.anno$DISEASE_CODE)]
  df$TX = qry
  df
}

tcga.plot = function(obj) {

  df = split(obj, obj$TX)[[1]]
  pr = df[df$type == "Primary Tumor", ]
  pr = data.frame(table(pr$id, pr$type))
  nrml = df[df$type == "Solid Tissue Normal", ]
  nrml = data.frame(table(nrml$id, nrml$type))
  pr$Freq2 = nrml$Freq[match(pr$Var1, nrml$Var1)]
  pr$label = paste0(pr$Var1, " (", pr$Freq, "/", pr$Freq2, ")")
  pr$label = gsub("/NA", "", pr$label)

  lbls = pr$label
  names(lbls) = pr$Var1

  obj$expr = log10(obj$expr + 1)

  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, Inf), symbols = c("****", "***", "**", "ns"))

  pl =
    ggplot(obj, aes(x=id, y = expr, fill = type)) +
    geom_boxplot(
      lwd = .3,
      outlier.size = .01,
      outlier.fill = NA,
      outlier.colour = "grey",
      position = position_dodge(preserve = "single")
    ) +
    mytheme(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x = element_text(size = rel(.8), colour = "white"),
      legend.position="bottom",
      axis.text.x = element_text(angle=45, hjust=1),
      axis.title.y = element_text(vjust = 2),
      plot.title = element_text(hjust = 0, colour = "black", face = "bold", size = rel(1)),
      legend.margin=margin(t=-5),
      panel.spacing.y=unit(.75, "lines"),
      panel.spacing.x=unit(.25, "lines"),
      legend.spacing.x = unit(1.5, 'pt'),
      legend.text = element_text(margin = margin(r = 15, unit = "pt")),
      legend.key.size = unit(20, "pt")
      # plot.margin = margin(r=-20, l = 5, t = 1)
    ) +
    xlab(NULL) +
    ylab("log10(FPKM-UQ + 1)") +
    guides(fill = guide_legend(title = NULL)) +
    scale_fill_manual(values = c("#997700", "#6699CC")) +
    scale_x_discrete(labels = lbls) +
    facet_grid( ~ type_coarse, space = "free_x", scales = "free_x") +
    stat_compare_means(
      aes(group = type), label = "p.signif", hide.ns = TRUE,
      vjust = 1, size = 2.5, symnum.args = symnum.args
    )

  g = ggplot_gtable(ggplot_build(pl))

  strip_both = which(grepl('strip-', g$layout$name))

  fills = tcga.anno$DISEASE_TYPE_COLOUR

  k = 1
  for (i in strip_both) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- my.cols[k]
    k <- k+1
  }

  return(cowplot::plot_grid(g))
}

ror2.tcga.exprs = tcga.plot(obj = tcga_epxrs(obj = se.ror2, qry = rownames(se.ror2)[1]))

ggsave2(
  ror2.tcga.exprs,
  filename = "analysis/publication/tcga_ror2_exprs.png",
  dpi = 300, width = 180, height = 45, units = "mm", bg = "white", scale = 2
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # Fig S4C: Singscore, highest correlated genes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se = se.ror2[, colnames(se.ror2) %in% colnames(ss.se)]
se = se[, colnames(ss.se)]

ftrs_cor = function(gs){
  gs.ftrs = geneIds(gs)
  gs.title = setName(gs)

  res.cor = cor(
    t(log10(assays(se)$fpkm + 1)),
    t(assays(ss.se[rownames(ss.se) %in% gs.ftrs, ])$log_fpkm),
    method = "spearman"
  )
  # eg = clusterProfiler::bitr(colnames(res.cor), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  # eg = eg[!duplicated(eg$ENSEMBL), ]

  df.cor = reshape2::melt(res.cor)
  df.cor$Var1 = NULL
  colnames(df.cor) = c("SYMBOL", "COR")
  # df.cor$SYMBOL = eg$SYMBOL[match(df.cor$ENSEMBL, eg$ENSEMBL)]
  df.cor = df.cor[order(df.cor$COR, decreasing = T), ]
  df.cor$SYMBOL = factor(df.cor$SYMBOL, levels = df.cor$SYMBOL)
  # df.cor = df.cor %>% dplyr::select(ENSEMBL, SYMBOL, COR)
  df.cor = df.cor %>% dplyr::select(SYMBOL, COR)
  df.cor$COR = round(df.cor$COR, 3)

  colnames(df.cor) = c("Gene", "r")
  t =
    ggpubr::ggtexttable(
    head(df.cor, 10),
    theme = ttheme(
      "light", base_size = 12,
      colnames.style = colnames_style(color = "black", face = "plain", fill = "white")
    ),
    rows = NULL
  )

  gridExtra::grid.arrange(
    top = grid::textGrob(gs.title,gp = gpar(fontsize=12, fontface = "bold")),
    t, heights=c(1, 0.12)
  )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supplement Figure 4
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ggsave2(
  cowplot::plot_grid(
    NULL,
    ror2.tcga.exprs,
    NULL,
    plot_grid(
      bp.kirc.stage,
      NULL,
      plot_grid(
        NULL,
        plot_grid(
          ftrs_cor(gs = ss.gs$cancersea.hm$unpaired$emt),
          NULL,
          ftrs_cor(ss.gs$cancersea.hm$unpaired$metastasis),
          NULL,
          ftrs_cor(ss.gs$cancersea.hm$unpaired$invasion),
          ncol = 5, rel_widths = c(1, .1, 1, .1, 1)
        ),
        nrow = 2, rel_heights = c(0, 1)
      ),
      ncol = 3, rel_widths = c(2, .15, 1.2),
      labels = c("B", "", "C"), label_size = 22, label_y = 1.075
    ),
    nrow = 4, rel_heights = c(.1, 1, .1, 1),
    labels = c("A", "", "", ""), label_size = 22, label_y = 1.1
  ),
  filename = "analysis/publication/Supp_Fig_4.png",
  dpi = 300, width = 180, height = 100, units = "mm", bg = "white", scale = 2
)


