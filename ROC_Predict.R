suppressMessages(suppressWarnings(install.packages("getopt", repos = "https://cloud.r-project.org/", 
 quiet = TRUE)))
suppressMessages(suppressWarnings(install.packages("optparse", repos = "https://cloud.r-project.org/", 
 quiet = TRUE)))
suppressMessages(suppressWarnings(library("getopt")))
suppressMessages(suppressWarnings(library("optparse")))

arguments = commandArgs(trailingOnly = TRUE)

option_list <- list(
make_option("--test.set", dest = "test.set"),
make_option("--roc.result", dest = "roc.result"),
make_option("--sig.threshold", dest = "sig.threshold", default = "0.05"),
make_option("--up.label", dest = "up.label", default = "1"),
make_option("--down.label", dest = "down.label", default = "2"),
make_option("--basename", dest = "basename", default = "ROC_prediction"),
make_option("--nperm", dest = "nperm", default = "1000"),
make_option("--random.seed", dest = "random.seed", default = as.integer(Sys.time()))
)

opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE, 
 args = arguments)$options

set.seed(seed = as.integer(opt$random.seed), kind = NULL)
nperm = as.integer(opt$nperm)

roc <- read.table(opt$roc.result, sep = "\t", header = TRUE, na = "", stringsAsFactors = FALSE)

n.orig <- nrow(roc)

if ("AUC.NOM.pValue" %in% colnames(roc)) {
 roc = roc[roc$AUC.NOM.pValue <= as.numeric(opt$sig.threshold), ]
 cat(paste0(nrow(roc), " gene sets remain after applying AUC.NOM.pValue<=0.05 threshold.\n"))
}

if ("MCC.NOM.pValue" %in% colnames(roc)) {
 roc = roc[roc$MCC.NOM.pValue <= as.numeric(opt$sig.threshold), ]
 cat(paste0(nrow(roc), " gene sets remain after applying MCC.NOM.pValue<=0.05 threshold.\n"))
}

if ("ssGSEA.Score.Wilcox.pValue" %in% colnames(roc)) {
 roc = roc[roc$ssGSEA.Score.Wilcox.pValue <= as.numeric(opt$sig.threshold), ]
 cat(paste0(nrow(roc), " gene sets remain after applying ssGSEA.Score.Wilcox.pValue<=0.05 threshold.\n"))
}

if ("MCC.FDR" %in% colnames(roc)) {
 roc = roc[roc$MCC.FDR <= as.numeric(opt$sig.threshold), ]
 cat(paste0(nrow(roc), " gene sets remain after applying MCC.FDR<=0.05 threshold.\n"))
}

n.filter <- nrow(roc)

if (n.orig == n.filter) {
 cat("Classifying using all gene sets.\n")
}
if (n.orig != n.filter) {
 cat(paste0("Classifying using ", n.filter, " gene sets that passed thresholds.\n"))

}

if (any(roc$AUC == round(0.5, 1))) {
 warning("Ambiguous AUC (0.5) Detected")
}

if (any(roc$Matthews.Correlation..MCC. == round(0))) {
 warning("Ambiguous MCC (0) Detected")
}

read.gct <- function(file) {
 if (is.character(file)) 
  if (file == "") 
   file <- stdin() else {
   file <- file(file, "r")
   on.exit(close(file))
  }
 if (!inherits(file, "connection")) 
  stop("argument `file' must be a character string or connection")

 # line 1 version
 version <- readLines(file, n = 1)

 # line 2 dimensions
 dimensions <- scan(file, what = list("integer", "integer"), nmax = 1, quiet = TRUE)
 rows <- dimensions[[1]]
 columns <- dimensions[[2]]
 # line 3 Name\tDescription\tSample names...
 column.names <- read.table(file, header = FALSE, quote = "", nrows = 1, sep = "\t", 
  fill = FALSE, comment.char = "")
 column.names <- column.names[3:length(column.names)]

 
 if (length(column.names) != columns) {
  stop(paste("Number of sample names", length(column.names), "not equal to the number of columns", 
   columns, "."))
 }

 colClasses <- c(rep(c("character"), 2), rep(c("double"), columns))

 x <- read.table(file, header = FALSE, quote = "", row.names = NULL, comment.char = "", 
  sep = "\t", colClasses = colClasses, fill = FALSE)
 row.descriptions <- as.character(x[, 2])
 data <- as.matrix(x[seq(from = 3, to = dim(x)[2], by = 1)])

 column.names <- column.names[!is.na(column.names)]

 colnames(data) <- column.names
 row.names(data) <- x[, 1]
 return(list(row.descriptions = row.descriptions, data = data))
}

if (all((roc$AUC >= round(0.5, 1)) == (roc$Matthews.Correlation..MCC. >= round(0)))) {

 samples <- read.gct(opt$test.set)

 # If MCC>0 Phenotype A > cutoff. If MCC<0 Phenotype A < cutoff
 sample.matrix <- as.data.frame(samples$data[roc$Name, , drop = FALSE])
 colnames(sample.matrix) <- colnames(samples$data)
 result.matrix <- matrix(rep("null"), nrow = ncol(samples$data), ncol = 3)
 rownames(result.matrix) = colnames(samples$data)
 colnames(result.matrix) = c("Phenotype_Call", "ES", "NOM pValue")

 pos.test <- as.integer(roc$Matthews.Correlation..MCC. >= round(0))
 neg.test <- (-1) * as.integer(roc$Matthews.Correlation..MCC. <= round(0))
 test.directions <- t(rbind(pos.test, neg.test))
 expected <- rowSums(test.directions)

 RES.expected <- cumsum(weight.expected <- roc$Matthews.Correlation..MCC. * expected)
  max.ES.expected <- max(RES.expected)

 phi <- matrix(nrow = length(colnames(samples$data)), ncol = nperm)
 p.vals <- matrix(0, nrow = length(colnames(samples$data)), ncol = 2)
 reshuffled.tests <- matrix(nrow = nperm, ncol = length(expected))

for (i in 1:nperm) {
reshuffled.tests[i,] <- sample(1:length(expected))
}

 for (i in 1:length(colnames(samples$data))) {

  sample.pos.test.up <- as.integer(sample.matrix[[i]] >= roc$cutoff.value..Youden.Index.)[as.logical(pos.test)]
  sample.pos.test.dn <- ((-1) * as.integer(sample.matrix[[i]] < roc$cutoff.value..Youden.Index.))[as.logical(pos.test)]
  sample.pos.test = sample.pos.test.up + sample.pos.test.dn

  sample.neg.test.up <- ((-1) * as.integer(sample.matrix[[i]] <= roc$cutoff.value..Youden.Index.))[as.logical(neg.test)]
  sample.neg.test.dn <- as.integer(sample.matrix[[i]] > roc$cutoff.value..Youden.Index.)[as.logical(neg.test)]
  sample.neg.test = sample.neg.test.up+sample.neg.test.dn

  observed <- c(sample.pos.test,sample.neg.test)

  RES <- cumsum(roc$Matthews.Correlation..MCC. * observed)
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > -min.ES) {
   # ES <- max.ES
   ES <- signif(max.ES, digits = 5)
  } else {
   # ES <- min.ES
   ES <- signif(min.ES, digits = 5)
  }

	for (j in 1:nperm) {
	 RES2 <- cumsum(roc$Matthews.Correlation..MCC. * observed[reshuffled.tests[j, 
	  ]])
	 max.ES2 <- max(RES2)
	 min.ES2 <- min(RES2)
	 if (max.ES2 > -min.ES2) {
	  # ES <- max.ES
	  ES2 <- signif(max.ES2, digits = 5)
	 } else {
	  # ES <- min.ES
	  ES2 <- signif(min.ES2, digits = 5)
	 }
	 phi[i, j] <- ES2
	}



 pos.phi <- NULL
 neg.phi <- NULL
 for (j in 1:nperm) {
	if (phi[i, j] >= 0) {
	 pos.phi <- c(pos.phi, phi[i, j])
	} else {
	 neg.phi <- c(neg.phi, phi[i, j])
	}
 }

 if (ES >= 0) {
	p.vals[i, 1] <- signif(sum(pos.phi >= ES)/length(pos.phi), digits = 5)
 } else {
	p.vals[i, 1] <- signif(sum(neg.phi <= ES)/length(neg.phi), digits = 5)
 }



  if (ES > 0) {
   result.matrix[i, 1] <- opt$up.label
#   result.matrix[i, 3] <- 1-(mcnemar.test(observed,expected)$p.value)
#   result.matrix[i, 3] <- chisq.test(table(observed,expected),simulate.p.value = TRUE)$p.value
   result.matrix[i, 3] <- signif(sum(pos.phi >= ES)/length(pos.phi), digits = 5)
  }
  if (ES < 0) {
   result.matrix[i, 1] <- opt$down.label
#   result.matrix[i, 3] <- 1-(mcnemar.test(observed,-expected)$p.value)
#   result.matrix[i, 3] <- chisq.test(table(observed,-expected),simulate.p.value = TRUE)$p.value
   result.matrix[i, 3] <- signif(sum(neg.phi <= ES)/length(neg.phi), digits = 5)
  }
  if (ES == 0) {
   result.matrix[i, 1] <- "NO_PREDICTION"
#   result.matrix[i, 3] <- "null"
   result.matrix[i, 3] <- "null"
  }

   result.matrix[i, 2] <- ES

 }
 result.matrix <- cbind(Sample_ID = rownames(result.matrix), result.matrix)
 write.table(result.matrix, paste0(opt$basename,"_result.txt"), sep = "\t", quote = FALSE, 
  col.names = TRUE, row.names = FALSE)

} else {
 message("Error: MCC and AUC classifiers conflict. Exiting.")
}
