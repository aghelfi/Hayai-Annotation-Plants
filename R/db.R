db <-
function() {
  dokohayai <- system.file(package="hayai")
  dokohayai <- paste(dokohayai,"data", sep="/")
  dokoftp <- "ftp://ftp.kazusa.or.jp/pub/hayai/hayai-annotation/plants"
  names <- c("phyta_aa.udb","phyta_ab.udb","phyta_ac.udb","phyta_ad.udb")
  for (j in 1:4) {
  destfile <- paste(dokohayai, names[j], sep="/")
  ftp <- paste(dokoftp, names[j], sep="/")
  download.file(ftp, destfile)
  }
}
