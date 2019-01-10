#' Parainflammation: A package for calculating the parainflammation score.
#'
#' @docType package
#' @name Parainflammation
NULL

#' CD45 curves learned from GTEx
#'
#' @format list:
#' \describe{
#'   \item{coeff}{a matrix with CD45 curves learned from GTEx}
#'   \item{tissues}{tissue names}
#' }
"fit_params"

#' The parainflammation signature
#'
#' @format oject:
"PI_signature"

#' Calculating the PI score
#'
#' \code{calculatePIscore} Returns the PI scores based on expression profiles.
#'
#' @param expr the gene expression data set. A matrix with row names as symbols and columns as samples.
#' @param fit_params the cruves of association with CD45 learned from GTEx.
#' @param tissue the tissue source of the expression profile. Must be from the list of tissues in fit_params$tissues
#'
#' @return the PI scores
calculatePIscore = function(expr,tissue) {
  require(gsva)
  require(GSEABase)
  id = which(fit_params$tissues==tissue)
  params = fit_params$coef[,(id*2)]
  expr.adjusted = adjustExpressionToCD45(expr,params)

  egc <- GeneSetCollection(PI_signature)
  pi = gsva(expr.adjusted,egc,method='ssgsea')
}

#' Adjusting expression to CD45
#'
#' \code{adjustExpressionToCD45} Returns the adjusted expression profile.
#'
#' @param expr the gene expression data set. A matrix with row names as symbols and columns as samples.
#' @param params a vector of the curves of association between genes and CD45 learned from GTEx.
#'
#' @return the adjusted expression profile
adjustExpressionToCD45 = function(expr,params) {
  A = intersect(rownames(expr),names(params))
  p = as.matrix(params[A])
  ex = as.matrix(expr[A,])
  cd45 = as.matrix(ex['PTPRC',])
  expr.adjusted = ex-p%*%t(cd45)
  expr.adjusted = expr.adjusted-apply(expr.adjusted,1,min)
}

#' #' Calculating the parameters
#' #'
#' #' \code{createParams} Returns a matrix of the parameters learned from GTEx.
#'
#' #' @return the parameters matrix
#' createParams = function() {
#'   load ('data/Toil_RSEM.Rdata')
#'   samples = read.table("data/TcgaTargetGTEX_phenotype.txt",sep="\t",header=TRUE,row.names=1, as.is=TRUE)
#'   rownames(samples) = gsub('-','.',rownames(samples))
#'   A = intersect(rownames(samples),colnames(toil))
#'   samples =samples[A,]
#'   toil = toil[,A]
#'   A = samples$X_study=="GTEX"
#'   gtex.samples =samples[A,]
#'   gtex = toil[,A]
#'   tissues = unique(gtex.samples$X_primary_site)
#'   tissues = sort(tissues[-which(tissues=="")])
#'   A = samples$X_study=="TCGA"
#'   tcga.samples =samples[A,]
#'   tcga = toil[,A]
#'
#'   cd45 = as.matrix(gtex['PTPRC',])
#'   cd45[cd45==0]=NA
#'   fit_params = list()
#'   fit_params$coef = matrix(NA,dim(gtex)[1],length(tissues)*2)
#'   fit_params$tissues = tissues
#'
#'   for (j in 1:length(tissues)) {
#'     print(paste(j,tissues[j]))
#'     A = gtex.samples$X_primary_site==tissues[j]
#'     x = as.matrix(cd45[A])
#'     for (i in 1:dim(gtex)[1]) {
#'       y = t(as.matrix(gtex[i,A]))
#'       z <- lm(y ~ x)
#'       fit_params$coef[i,(j*2)-1] = coef(z)[1]
#'       fit_params$coef[i,j*2] = coef(z)[2]
#'     }
#'   }
#'   rownames(fit_params$coef) = rownames(gtex)
#'   save(fit_params,file='data/fit_params.Rdata')
#' }

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/R-dev",
    devtools.install.args = "",
    devtools.name = "Dvir Aran",
    devtools.desc.author = '"Dvir Aran <dvir.aran@ucsf.edu> [aut, cre]"',
    devtools.desc.license = "GPL-3",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
  )
  toset <- !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])

  require(GSVA)
  require(GSEABase)

  invisible()
}
