\name{hayai}
\alias{hayai}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hayai-Annotation Plants}
\description{ayai-Annotation Plants is an R-package with a shiny-based interface for an ultra-fast and accurate gene annotation s
ystem for plant species using R. The pipeline combines the sequence-similarity searches, using USEARCH against UniProtKB (taxonomy
 Embryophyta), with a functional annotation step. Hayai-Annotation Plants provides five layers of annotation: 1) gene name; 2) gen
e ontology terms consisting of its three main domains (Biological Process, Molecular Function, and Cellular Component); 3) enzyme 
commission number; 4) protein existence level; 5) and evidence type.}
\usage{
hayai()
}
%- maybe also 'usage' for other objects documented here.
\details{After install.package ("hayai_1.0.tar.gz"), install database.
To install database, type:
db ()
This command will download Embryophyta/UniProt database from ftp.kazusa.or.jp into hayai/data directory permanentely.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{https://github.com/kdri-genomics/Hayai-Annotation-Plants/}
\author{Andrea Ghelfi
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function () 
{
    options(shiny.maxRequestSize = 1000 * 1024^2)
    loadData <- function() {
        if (exists("responses")) {
            userInputFile <- read.csv("df.csv")
            userInputFile <- as.character(userInputFile$datapath)
            path_hayai <- path.package("hayai")[1]
            tmp <- gsub("\\..*", "", userInputFile)
            pos <- gregexpr("/", tmp)
            pos <- summary(pos[[1]])[6][[1]]
            tmp <- substr(tmp, 1, pos)
            title_fasta <- "awk '{gsub(\"[|]\",\"_\"); print $1}' USERSFILE > TMPquery.fasta "
            query_title <- gsub("USERSFILE", userInputFile, title_fasta)
            query_title <- gsub("TMP", tmp, query_title)
            system(query_title)
            system("rm df.csv")
            if (align == "local") {
                temp_usearch1 <- "./usearch -usearch_useralignment TMPquery.fasta -db PHYTA -id userseqID -maxaccepts usermaxaccepts -blast6out TMPRESULTSUSEARCH -evalue 1e-USERSEVALUE -query_cov USERSQUERYCOV"
                usearch_run1 <- gsub("useralignment", "local", 
                  temp_usearch1)
            }
            if (align == "global") {
                temp_usearch1 <- "./usearch -usearch_useralignment TMPquery.fasta -db PHYTA -id userseqID -maxaccepts usermaxaccepts -blast6out TMPRESULTSUSEARCH -evalue 1e-USERSEVALUE -query_cov USERSQUERYCOV"
                usearch_run1 <- gsub("useralignment", "global", 
                  temp_usearch1)
            }
            usearch_run1 <- gsub("USERSEVALUE", e_value$e_value[[1]], 
                usearch_run1)
            seq_id <- as.numeric(levels(seq_id$seq_id))/100
            usearch_run1 <- gsub("userseqID", seq_id, usearch_run1)
            usearch_run1 <- gsub("usermaxaccepts", hits$hits[[1]], 
                usearch_run1)
            usearch_run1 <- gsub("TMP", tmp, usearch_run1)
            query_cov <- as.numeric(levels(query_cov$query_cov))/100
            usearch_run1 <- gsub("USERSQUERYCOV", query_cov, 
                usearch_run1)
            for (p in 1:4) {
                database <- paste(path_hayai, phyta[p], sep = "/")
                usearch_run2 <- gsub("PHYTA", database, usearch_run1)
                usearch_run2 <- gsub("RESULTSUSEARCH", results[p], 
                  usearch_run2)
                system(usearch_run2)
                print(usearch_run2)
            }
            write.csv(usearch_run1, "hayai_usearch.log", row.names = F)
            linecat <- "cat TMPresults_usearch_a*.b6 > TMPoutput_usearch.txt"
            linecat <- gsub("TMP", tmp, linecat)
            system(linecat)
            lineawk <- "awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' TMPoutput_usearch.txt | awk -F\"|\" '{print $3}' > TMPtable_usearch.txt"
            lineawk <- gsub("TMP", tmp, lineawk)
            system(lineawk)
            lineinput <- "TMPtable_usearch.txt"
            lineinput <- gsub("TMP", tmp, lineinput)
            usearch <- read.table(lineinput)
            colnames(usearch) <- c("uniprot", "query", "seqID", 
                "length", "mis", "gaps", "startquery", "endquery", 
                "starttarget", "endtarget", "evalue", "score")
            temp <- merge(usearch, info_database, by = "uniprot", 
                all.x = T)
            temp_bp <- merge(temp, bp_table, by = "uniprot", 
                all.x = T)
            temp_mf <- merge(temp, mf_table, by = "uniprot")
            temp_cc <- merge(temp, cc_table, by = "uniprot")
            temp_ec <- temp[!is.na(temp$EC), ]
            if (organism$Type_algorithm == "protein_existence") {
                if (align == "local") {
                  temp_bp <- temp_bp[order(temp_bp$Protein_Evidence, 
                    temp_bp$evidence_type, -temp_bp$score), ]
                  temp_mf <- temp_mf[order(temp_mf$Protein_Evidence, 
                    temp_mf$evidence_type, -temp_mf$score), ]
                  temp_cc <- temp_cc[order(temp_cc$Protein_Evidence, 
                    temp_cc$evidence_type, -temp_cc$score), ]
                  temp_ec <- temp_ec[order(temp_ec$Protein_Evidence, 
                    temp_ec$evidence_type, -temp_ec$score), ]
                }
                if (align == "global") {
                  temp_bp <- temp_bp[order(temp_bp$Protein_Evidence, 
                    temp_bp$evidence_type, -temp_bp$seqID), ]
                  temp_mf <- temp_mf[order(temp_mf$Protein_Evidence, 
                    temp_mf$evidence_type, -temp_mf$seqID), ]
                  temp_cc <- temp_cc[order(temp_cc$Protein_Evidence, 
                    temp_cc$evidence_type, -temp_cc$seqID), ]
                  temp_ec <- temp_ec[order(temp_ec$Protein_Evidence, 
                    temp_ec$evidence_type, -temp_ec$seqID), ]
                }
            }
            if (organism$Type_algorithm == "score") {
                temp_bp <- temp_bp[order(-temp_bp$seqID, temp_bp$Protein_Evidence, 
                  temp_bp$evidence_type), ]
                temp_mf <- temp_mf[order(-temp_mf$seqID, temp_mf$Protein_Evidence, 
                  temp_mf$evidence_type), ]
                temp_cc <- temp_cc[order(-temp_cc$seqID, temp_cc$Protein_Evidence, 
                  temp_cc$evidence_type), ]
                temp_ec <- temp_ec[order(-temp_ec$seqID, temp_ec$Protein_Evidence, 
                  temp_ec$evidence_type), ]
            }
            temp_bp <- temp_bp[!duplicated(temp_bp[, c("query")]), 
                ]
            temp_mf <- temp_mf[!duplicated(temp_mf[, c("query")]), 
                ]
            temp_mf <- temp_mf[, c("query", "uniprot", "GO_MF", 
                "GO_MF_name")]
            temp_mf2 <- temp_mf[, c("query", "GO_MF", "GO_MF_name")]
            temp_cc <- temp_cc[!duplicated(temp_cc[, c("query")]), 
                ]
            temp_cc <- temp_cc[, c("query", "uniprot", "GO_CC", 
                "GO_CC_name")]
            temp_cc2 <- temp_cc[, c("query", "GO_CC", "GO_CC_name")]
            anota <- merge(temp_bp, temp_mf2, by = "query", all.x = T, 
                all.y = T)
            anota <- merge(anota, temp_cc2, by = "query", all.x = T, 
                all.y = T)
            anota$evidence_type <- gsub("xIEA", "IEA", anota$evidence_type)
            write.table(anota, "hayai_annotation.csv", row.names = F, 
                col.names = T, sep = ",")
            temp_bp <- temp_bp[, c("query", "uniprot", "GO_BP")]
            temp_bp <- temp_bp[!is.na(temp_bp$GO_BP), ]
            temp_bp <- temp_bp[, c("query", "uniprot")]
            bp <- merge(temp_bp, bp_table, by = "uniprot")
            write.table(bp, "GO_BP_table.csv", row.names = F, 
                col.names = T, sep = ",")
            temp_mf <- temp_mf[!is.na(temp_mf$GO_MF), ]
            temp_mf <- temp_mf[, c("query", "uniprot")]
            mf <- merge(temp_mf, mf_table, by = "uniprot")
            write.table(mf, "GO_MF_table.csv", row.names = F, 
                col.names = T, sep = ",")
            temp_cc <- temp_cc[!is.na(temp_cc$GO_CC), ]
            temp_cc <- temp_cc[, c("query", "uniprot")]
            cc <- merge(temp_cc, cc_table, by = "uniprot")
            write.table(cc, "GO_CC_table.csv", row.names = F, 
                col.names = T, sep = ",")
            temp_ec <- temp_ec[!duplicated(temp_ec[, c("query")]), 
                c("uniprot", "query", "EC")]
            write.table(temp_ec, "EC_table.csv", row.names = F, 
                col.names = T, sep = ",")
            sum_bp <- sort(table(bp$GO_BP_name), decreasing = T)
            sum_bp <- as.data.frame(sum_bp)
            sum_bp <- sum_bp[sum_bp$Freq > 0, ]
            colnames(sum_bp) <- c("GO_BP_name", "Counts")
            sum_bp$GO_BP_name <- strtrim(sum_bp$GO_BP_name, 55)
            bar_bp <- sum_bp[50:1, ]
            y <- bar_bp[50, 2]
            z <- 1.25 * y
            pdf("GO_BP.pdf", width = 11, height = 11)
            par(mar = c(2, 2, 6, 2), oma = c(0.5, 16, 5, 0.5))
            barplot(bar_bp$Counts, names = bar_bp$GO_BP_name, 
                las = 1, horiz = T, cex.names = 0.75, space = 2, 
                xlim = c(0, z), axes = F)
            axis(3)
            mtext("Number of Genes", side = 3, line = 2, font = 1)
            title(main = "GO Biological Process - Top 50 - Gene Level", 
                font.main = 4, line = 4)
            dev.off()
            write.table(sum_bp, "GO_BP_counts.csv", row.names = F, 
                col.names = T, sep = ",")
            sum_mf <- sort(table(mf$GO_MF_name), decreasing = T)
            sum_mf <- as.data.frame(sum_mf)
            sum_mf <- sum_mf[sum_mf$Freq > 0, ]
            colnames(sum_mf) <- c("GO_MF_name", "Counts")
            sum_mf$GO_MF_name <- strtrim(sum_mf$GO_MF_name, 55)
            bar_mf <- sum_mf[50:1, ]
            y <- bar_mf[50, 2]
            z <- 1.25 * y
            pdf("GO_MF.pdf", width = 11, height = 11)
            par(mar = c(2, 2, 6, 2), oma = c(0.5, 16, 5, 0.5))
            barplot(bar_mf$Counts, names = bar_mf$GO_MF_name, 
                las = 1, horiz = T, cex.names = 0.75, space = 2, 
                xlim = c(0, z), axes = F)
            axis(3)
            mtext("Number of Genes", side = 3, line = 2, font = 1)
            title(main = "GO Molecular Function - Top 50 - Gene Level", 
                font.main = 4, line = 4)
            dev.off()
            write.table(sum_mf, "GO_MF_counts.csv", row.names = F, 
                col.names = T, sep = ",")
            sum_cc <- sort(table(cc$GO_CC_name), decreasing = T)
            sum_cc <- as.data.frame(sum_cc)
            sum_cc <- sum_cc[sum_cc$Freq > 0, ]
            colnames(sum_cc) <- c("GO_CC_name", "Counts")
            sum_cc$GO_CC_name <- strtrim(sum_cc$GO_CC_name, 55)
            bar_cc <- sum_cc[50:1, ]
            y <- bar_cc[50, 2]
            z <- 1.25 * y
            pdf("GO_CC.pdf", width = 11, height = 11)
            par(mar = c(2, 2, 6, 2), oma = c(0.5, 16, 5, 0.5))
            barplot(bar_cc$Counts, names = bar_cc$GO_CC_name, 
                las = 1, horiz = T, cex.names = 0.75, space = 2, 
                xlim = c(0, z), axes = F)
            axis(3)
            mtext("Number of Genes", side = 3, line = 2, font = 1)
            title(main = "GO Cellular Component - Top 50 - Gene Level", 
                font.main = 4, line = 4)
            dev.off()
            write.table(sum_cc, "GO_CC_counts.csv", row.names = F, 
                col.names = T, sep = ",")
            sum_ec <- sort(table(droplevels(temp_ec$EC)), decreasing = T)
            bar_ec <- sum_ec[50:1]
            y <- bar_ec[[50]]
            z <- 1.25 * y
            pdf("EC_codes.pdf", width = 11, height = 11)
            par(mar = c(2, 2, 6, 2), oma = c(0.5, 5, 5, 1))
            barplot(bar_ec, las = 1, horiz = T, cex.names = 0.75, 
                space = 2, xlim = c(0, z), axes = F)
            axis(3)
            mtext("Number of Genes", side = 3, line = 2, font = 1)
            title(main = "Enzyme Commission - Top 50", font = 2, 
                line = 4)
            dev.off()
            sum_ec <- as.data.frame(sum_ec)
            colnames(sum_ec) <- (c("EC", "Counts"))
            write.table(sum_ec, "EC_counts.csv", row.names = F, 
                col.names = T, sep = ",")
            write.table(sum_ec[, 1], "unique_EC.csv", row.names = F, 
                col.names = F, quote = F)
            anota
        }
    }
    write.table(anota, "hayai_annotation.csv", col.names = F, 
        row.names = F, sep = ",")
    saveData <- function(data) {
        data <- as.data.frame(t(data))
        if (exists("responses")) {
            align <<- data[1]
            hits <<- data[2]
            seq_id <<- data[3]
            organism <<- data[4]
            e_value <<- data[5]
            query_cov <<- data[6]
            responses <<- rbind(responses, data)
            write.csv(responses, "hayai_annotation.parameters", 
                row.names = F)
            loadData()
        }
        else {
            align <<- data[1]
            hits <<- data[2]
            seq_id <<- data[3]
            organism <<- data[4]
            e_value <<- data[5]
            query_cov <<- data[6]
            responses <<- data
            write.csv(responses, "hayai_annotation.parameters", 
                row.names = F)
            loadData()
        }
    }
    if (interactive()) {
        ui <- fluidPage(h1(id = "big-heading", "Hayai-Annotation Plants"), 
            tags$style(HTML("#big-heading{color: DarkGreen;}")), 
            br(), h5(id = "small-heading", "Kazusa DNA Research Institute"), 
            tags$style(HTML("#small-heading{color: green;}")), 
            br(), br(), br(), fluidRow(column(3, radioButtons("Align", 
                "Type of Alignment", choices = list(Local = "local", 
                  Global = "global"), selected = "local"), br(), 
                radioButtons("Type_algorithm", "Type of Algorithm", 
                  choices = list(`Protein Existence Level` = "protein_existence", 
                    `Alignment Score` = "score"), selected = "protein_existence")), 
                column(3, numericInput("hits", "Max hits per query", 
                  value = 1, min = 1, max = 20, step = 1), br(), 
                  numericInput("e_value", "Evalue 1e-", value = 6, 
                    min = 1, max = 100, step = 1)), column(3, 
                  numericInput("seq_id", "Minimum Sequence Identity (\%)", 
                    value = 80, min = 40, max = 100, step = 1), 
                  br(), numericInput("query_cov", "Minimum Query Coverage (\%)", 
                    value = 80, min = 20, max = 100, step = 1))), 
            br(), br(), fileInput("userInput", "Upload FASTA File", 
                multiple = F, accept = c(".fasta", ".fa", ".faa", 
                  ".fna")), actionButton("submit", "Submit"), 
            br(), br(), br(), downloadLink("downloadData", "Download"), 
            br(), br(), br(), fluidRow(DT::dataTableOutput("contents")))
        server <- function(input, output, session) {
            observeEvent(input$userInput, {
                inFile <- input$userInput
                if (is.null(inFile)) 
                  return(NULL)
                write.csv(inFile, "df.csv", row.names = F)
            })
            formData <- reactive({
                data <- sapply(fields, function(x) input[[x]])
                data
            })
            observeEvent(input$submit, {
                saveData(formData())
            })
            output$contents <- DT::renderDataTable(DT::datatable({
                input$submit
                anota <- read.table("hayai_annotation.csv", header = T, 
                  sep = ",")
                anota
            }))
            output$downloadData <- downloadHandler(filename = function() {
                paste("output_HayaiAnnotation", "zip", sep = ".")
            }, content = function(fname) {
                fs <- c("hayai_annotation.csv", "GO_BP_counts.csv", 
                  "GO_MF_counts.csv", "GO_CC_counts.csv", "GO_BP.pdf", 
                  "GO_MF.pdf", "GO_CC.pdf", "EC_counts.csv", 
                  "EC_codes.pdf", "GO_BP_table.csv", "GO_MF_table.csv", 
                  "GO_CC_table.csv", "EC_table.csv", "unique_EC.csv")
                zip(zipfile = fname, files = fs)
            }, contentType = "application/zip")
        }
        shinyApp(ui, server)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
