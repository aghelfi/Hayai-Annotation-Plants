db <-
function() {
  dokohayai <- system.file(package="hayai")
  dokodata <- paste(dokohayai,"data", sep="/")
  dokoR <- paste(dokohayai,"R", sep="/")
  dokoftp <- "ftp://ftp.kazusa.or.jp/pub/hayai/hayai-annotation/plants"

  names <- c("phyta_aa.udb","phyta_ab.udb","phyta_ac.udb","phyta_ad.udb")
  l <- length(names)
  for (j in 1:l) {
  destfile <- paste(dokodata, names[j], sep="/")
  ftp <- paste(dokoftp, names[j], sep="/")
  download.file(ftp, destfile)
  }
}
