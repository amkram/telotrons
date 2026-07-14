#!/usr/bin/env Rscript
# Run NuPoP predNuPoP(species=0, model=4) on each FASTA path listed in the file given as arg1.
# predNuPoP writes "<basename>_Prediction4.txt" into the current working directory, so we setwd to
# each input file's own directory before calling it. This (a) keeps the repo root clean and
# (b) prevents with/<lid>.fa and without/<lid>.fa (identical basenames) from overwriting each other.
suppressMessages(library(NuPoP))
files <- readLines(commandArgs(trailingOnly=TRUE)[1])
owd <- normalizePath(getwd())
for (f in files) {
  d <- dirname(f); b <- basename(f)
  outp <- file.path(d, paste0(b, "_Prediction4.txt"))
  if (file.exists(outp) && file.info(outp)$size > 0) next
  setwd(d)
  tryCatch(predNuPoP(b, species=0, model=4),
           error=function(e) cat("ERR", f, conditionMessage(e), "\n"))
  setwd(owd)
}
