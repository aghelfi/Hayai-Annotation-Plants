db <-
function() {
  dokohayai <- system.file(package="hayai")
  dokodata <- paste(dokohayai,"data", sep="/")
  dokoR <- paste(dokohayai,"R", sep="/")
  dokoftp <- "ftp://ftp.kazusa.or.jp/pub/hayai/hayai-annotation/plants"

  names2 <- c("sysdata.rda") 
  destfile2 <- paste(dokoR, names2, sep="/")
  ftp2 <- paste(dokoftp, names2, sep="/")
  download.file(ftp2, destfile2)
 
  names <- c("phyta_aa.udb","phyta_ab.udb","phyta_ac.udb","phyta_ad.udb","info_database.rda", "mf_table.rda","bp_table.rda","cc_table.rda")
  l <- length(names)
  for (j in 1:l) {
  destfile <- paste(dokodata, names[j], sep="/")
  ftp <- paste(dokoftp, names[j], sep="/")
  download.file(ftp, destfile)
  }
}
