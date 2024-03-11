# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Teile gene set in Hoch- und Runterregluierte sets auf
# Gene sets müssen durch Suffix _DN/_UP gekennzeichnet sein
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
split_gene_sets = function(gene.sets.l){
  nms = names(gene.sets.l)
  nms.dir = nms[grepl("_UP$|_DN$", nms)]
  up = nms.dir[grepl("_UP$", nms.dir)]
  names(up) = gsub("_UP$|_DN$", "", up)
  dn = nms.dir[grepl("_DN$", nms.dir)]
  names(dn) = gsub("_UP$|_DN$", "", dn)
  dn = dn[names(up)]
  stopifnot(identical(names(up), names(dn)))

  up.g.s = lapply(gene.sets.l[up], function(x){
    setName(x) = gsub("_UP$", "", setName(x))
    x
  })
  dn.g.s = lapply(gene.sets.l[dn], function(x){
    setName(x) = gsub("_DN$", "", setName(x))
    x
  })

  # List with up and down regulated gene sets
  gene.sets.l.paired = list(
    "UP" = up.g.s,
    "DN" = dn.g.s
  )
  gene.sets.l.unpaired = gene.sets.l[!names(gene.sets.l) %in% c(dn, up, names(dn), names(up))]
  list(paired = gene.sets.l.paired, unpaired = gene.sets.l.unpaired)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Log-Rank
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
runKaplanTab = function(
  dat,
  surv
) {
  out = t(apply(dat, 1, function(y, s) {
    res = survdiff(surv ~ y)
    p.val = pchisq(res$chisq, length(res$n)-1, lower.tail = FALSE)
    c("LOG_RANK_PVALUE" = p.val)
  }, s=surv))
  out
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cox
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
run_cox_gsea = function(
  obj,
  .case.weights = NULL,
  .surv =  NULL,
  .threads = 1
) {
  exprs.mat = obj

  exprs.mat.spl = split(as.matrix(exprs.mat), rownames(exprs.mat))

  coxModel.l = parallel::mclapply(exprs.mat.spl, function(gene) {
    if (is.null(.case.weights)) {
      res = summary(survival::coxph(.surv ~ gene))
    } else {
      res = summary(survival::coxph(.surv ~ gene, weights = .case.weights))
    }
    res = c(
      "HR" = signif(res$coef[2], digits=3),
      "logHR" = signif(res$coef[1], digits=3),
      "SE_logHR" = signif(res$coef[3], digits=3),
      # "L95CI" = signif(res$conf.int[,"lower .95"], 3),
      # "U95CI" = signif(res$conf.int[,"upper .95"], 3),
      "L95CI" = res$coef[1] - res$coef[3] * 1.96,
      "U95CI" = res$coef[1] + res$coef[3] * 1.96,
      "Pval" = res$coef[5]
    )
  }, mc.cores = .threads)

  coxModel = data.frame(do.call("rbind", coxModel.l))
  coxModel$Padj = stats::p.adjust(coxModel$Pval, method = "BH")
  coxModel
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Singscore & Permutation test
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ssgsea = function(
  g.s.l = NULL,
  .rankData = NULL,
  paired = F,
  permut = F,
  nbr.perm = 1000,
  ncores = 20,
  alpha = 0.05,
  study.id

){

  if(paired == T){
    if(length(g.s.l) != 2){
      stop("Error: Paired gene sets, but list has not length == 2")
    }
  }

  simpleScore.output = list()
  getPvals.output = list()
  permuteResult.output = list()
  pval.l = list()
  not.in.rankData = list()
  ss.res.l = list()

  if(paired == F){

    if(length(g.s.l) == 0) {
      return(list(ss.sum = data.frame(), ss.res = matrix(nrow = 0, ncol = ncol(.rankData))))
    }

    for (j in 1:length(g.s.l)) {

      g.s = g.s.l[[j]]
      g.s.name = GSEABase::setName(g.s)

      ftrs.avail.bool = GSEABase::geneIds(g.s) %in% rownames(.rankData)
      ftrs.not.avail = table(ftrs.avail.bool)["FALSE"]
      GSEABase::geneIds(g.s) = GSEABase::geneIds(g.s)[ftrs.avail.bool]

      scoredf = singscore::simpleScore(.rankData, upSet = g.s,  knownDirection = T)
      simpleScore.output[[g.s.name]] = scoredf

      if(permut == T){
        permuteResult = singscore::generateNull(
          upSet = g.s,
          rankData = .rankData,
          knownDirection = T,
          B = nbr.perm,
          seed = 1234,
          ncores = ncores
        )
        permuteResult.output[[g.s.name]] = permuteResult
        pvals = singscore::getPvals(permuteResult, scoredf)
        getPvals.output[[g.s.name]] = pvals
        pval.l[[g.s.name]] = t(data.frame(pvals))
      }

      not.in.rankData[[j]] = ftrs.not.avail
      scoredf = scoredf[, 1, drop = F]
      colnames(scoredf) = g.s.name
      ss.res.l[[g.s.name]] = scoredf
    }
  }
  if(paired == T){

    if(length(g.s.l[[1]]) == 0) {
      return(list(ss.sum = data.frame(), ss.res = matrix(nrow = 0, ncol = ncol(.rankData))))
    }

    for (j in 1:length(g.s.l[[1]])) {

      g.s.up = g.s.l[[1]][[j]]
      g.s.dn = g.s.l[[2]][[j]]
      g.s.name = GSEABase::setName(g.s.up)

      ftrs.avail.up.bool = GSEABase::geneIds(g.s.up) %in% rownames(.rankData)
      ftrs.avail.dn.bool = GSEABase::geneIds(g.s.dn) %in% rownames(.rankData)
      ftrs.not.avail = table(c(ftrs.avail.up.bool, ftrs.avail.dn.bool))["FALSE"]
      GSEABase::geneIds(g.s.up) = GSEABase::geneIds(g.s.up)[ftrs.avail.up.bool]
      GSEABase::geneIds(g.s.dn) = GSEABase::geneIds(g.s.dn)[ftrs.avail.dn.bool]

      scoredf = singscore::simpleScore(.rankData, upSet = g.s.up, downSet = g.s.dn,  knownDirection = T)
      simpleScore.output[[g.s.name]] = scoredf

      if(permut == T){
        permuteResult = singscore::generateNull(
          upSet = g.s.up, downSet = g.s.dn,
          rankData = .rankData,
          knownDirection = T,
          B = nbr.perm,
          seed = 1234,
          ncores = ncores
        )
        permuteResult.output[[g.s.name]] = permuteResult
        pvals = singscore::getPvals(permuteResult, scoredf)
        getPvals.output[[g.s.name]] = pvals
        # plotNull(permuteResult, scoredf, pvals, sampleNames = names(pvals)[3])
        pval.l[[g.s.name]] = t(data.frame(pvals))
      }

      not.in.rankData[[j]] = ftrs.not.avail
      scoredf = scoredf[, 1, drop = F]
      colnames(scoredf) = g.s.name
      ss.res.l[[g.s.name]] = scoredf
    }
  }

  if(permut == T){
    ss.sum = data.frame(t(sapply(pval.l, function (x) {
      sign = unname(table(x < alpha)["TRUE"])
      sign[is.na(sign)] = 0
      c(length(x) - sign, sign)
    })))
    ss.sum$X1 = NULL
    colnames(ss.sum) = "SIGN_ENRICHED_SAMPLES"
    ss.sum$SIGN_ENRICHED_SAMPLES_PERC = round((ss.sum$SIGN_ENRICHED_SAMPLES * 100) / ncol(.rankData), 3)
    ss.sum$SAMPLES = ncol(.rankData)
  } else {
    ss.sum = data.frame(
      SAMPLES = rep(ncol(.rankData), length(names(ss.res.l))),
      row.names = names(ss.res.l)
    )
  }

  if(paired == F){
    ss.sum$GENE_SET_FTRS = unname(unlist(lapply(g.s.l, function(x){length(GSEABase::geneIds(x))})))
  }else{
    ss.sum$GENE_SET_FTRS = unname(unlist(lapply(g.s.l[[1]], function(x){length(GSEABase::geneIds(x))}))) +
      unname(unlist(lapply(g.s.l[[2]], function(x){length(GSEABase::geneIds(x))})))
  }
  ss.sum$GENE_SET_FTRS_NOT_IN_FLTRD_DATA = unlist(not.in.rankData)
  ss.sum$GENE_SET_FTRS_NOT_IN_FLTRD_DATA[is.na(ss.sum$GENE_SET_FTRS_NOT_IN_FLTRD_DATA)] = 0
  ss.sum$GENE_SET = rownames(ss.sum)
  ss.sum$COHORT = study.id
  # print(ss.sum)
  list(
    ss.sum = ss.sum,
    ss.res = t(do.call("cbind", ss.res.l)),
    simpleScore.output = simpleScore.output,
    permuteResult.output = permuteResult.output,
    getPvals.output = getPvals.output
  )
}
