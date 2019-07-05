hayai <-
function() {
  options(shiny.maxRequestSize = 1000*1024^2)
  loadData <- function() {
    if (exists("responses")) {
      userInputFile <- read.csv("df.csv")
      userInputFile <- as.character(userInputFile$datapath)
#      path_hayai <- "/Users/andreaghelfi/projects/database" # test
      path_hayai <- path.package("hayai")[1] # r-package
      tmp <- gsub("\\..*","",userInputFile)
      pos <- gregexpr("/",tmp)
      pos <- summary(pos[[1]])[6][[1]]
      tmp <- substr(tmp, 1, pos)
      title_fasta <- 'awk \'{gsub("[|]","_"); print $1}\' USERSFILE > TMPquery.fasta '
      query_title <- gsub ("USERSFILE", userInputFile, title_fasta)
      query_title <- gsub ("TMP", tmp, query_title)
      system (query_title)
      system ('rm df.csv')
      if (align == "local"){
        temp_usearch1 <- 'usearch -usearch_useralignment TMPquery.fasta -db PHYTA -id userseqID -maxaccepts usermaxaccepts -blast6out TMPRESULTSUSEARCH -evalue 1e-USERSEVALUE -query_cov USERSQUERYCOV -threads 4'
        usearch_run1 <- gsub ("useralignment", "local", temp_usearch1)
        }
      if (align == "global"){
        temp_usearch1 <- 'usearch -usearch_useralignment TMPquery.fasta -db PHYTA -id userseqID -maxaccepts usermaxaccepts -blast6out TMPRESULTSUSEARCH -evalue 1e-USERSEVALUE -query_cov USERSQUERYCOV -threads 4'
        usearch_run1 <- gsub ("useralignment", "global", temp_usearch1)
      }
      usearch_run1 <- gsub ("USERSEVALUE", e_value$e_value[[1]], usearch_run1)
      seq_id <- as.numeric(levels(seq_id$seq_id))/100
      usearch_run1 <- gsub ("userseqID", seq_id, usearch_run1)
      usearch_run1 <- gsub ("usermaxaccepts", hits$hits[[1]], usearch_run1)
      usearch_run1 <- gsub ("TMP", tmp, usearch_run1) # changed filename
      query_cov <- as.numeric(levels(query_cov$query_cov))/100
      usearch_run1 <- gsub ("USERSQUERYCOV", query_cov, usearch_run1)
      for (p in 1:4) {
        database <- paste(path_hayai,phyta[p], sep="/")
        usearch_run2 <- gsub ("PHYTA", database, usearch_run1)
        usearch_run2 <- gsub ("RESULTSUSEARCH", results[p], usearch_run2)
        system(usearch_run2)
        print(usearch_run2)
      }
      write.csv(usearch_run1, "hayai_usearch.log", row.names=F)
      linecat <- 'cat TMPresults_usearch_a*.b6 > TMPoutput_usearch.txt'
      linecat <- gsub ("TMP", tmp, linecat) # changed filename
      system (linecat)
#      print (linecat)
      lineawk <- 'awk -F"\t" \'{gsub("[#]",""); print $1"|"$2"|"$3"|"$4"|"$5"|"$6"|"$7"|"$8"|"$9"|"$10"|"$11"|"$12 }\' TMPoutput_usearch.txt > TMPtable_usearch.txt'
      lineawk <- gsub("TMP", tmp, lineawk)
#      print(lineawk)
      system (lineawk)
      lineinput <- "TMPtable_usearch.txt"
      lineinput <- gsub("TMP", tmp, lineinput)
      usearch <- read.table(lineinput, sep="|")
#      print(usearch)
      colnames(usearch) <-  c("Query_ID", "UNIPROT_AC","UNIPROT_ID","Protein_name","Gene_name","EC","Protein_existence", "Sequence_identity", "Alignment_length", "Mismatches", "Gaps", "Start_query", "End_query", "Start_target", "End_target", "E_value", "Bitscore")
      temp_ec <- usearch[!is.na(usearch$EC), c("Query_ID", "EC", "Protein_existence", "Sequence_identity", "Bitscore")]
      temp <- usearch[,c("Query_ID", "UNIPROT_ID", "Protein_existence", "Sequence_identity", "Bitscore")]
      temp_bp <- merge(usearch, bp_table, by ="UNIPROT_ID", all.x=T)
      temp_mf <- merge(temp, mf_table, by ="UNIPROT_ID")
      temp_cc <- merge(temp, cc_table, by ="UNIPROT_ID")
      if (organism$Type_algorithm == "protein_existence"){
        if (align == "local"){
          temp_bp <- temp_bp[order(temp_bp$Protein_existence, temp_bp$Evidence, -temp_bp$Bitscore),]
          temp_mf <- temp_mf[order(temp_mf$Protein_existence, temp_mf$Evidence, -temp_mf$Bitscore),]
          temp_cc <- temp_cc[order(temp_cc$Protein_existence, temp_cc$Evidence, -temp_cc$Bitscore),]
          temp_ec <- temp_ec[order(temp_ec$Protein_existence, -temp_ec$Bitscore),]
        }
        if (align == "global"){
          temp_bp <- temp_bp[order(temp_bp$Protein_existence, temp_bp$Evidence, -temp_bp$Sequence_identity),]
          temp_mf <- temp_mf[order(temp_mf$Protein_existence, temp_mf$Evidence, -temp_mf$Sequence_identity),]
          temp_cc <- temp_cc[order(temp_cc$Protein_existence, temp_cc$Evidence, -temp_cc$Sequence_identity),]
          temp_ec <- temp_ec[order(temp_ec$Protein_existence, -temp_ec$Sequence_identity),]
        }
      }
      if (organism$Type_algorithm == "score"){
        if (align == "local"){
          temp_bp <- temp_bp[order(-temp_bp$Bitscore, temp_bp$Protein_existence, temp_bp$Evidence),]
          temp_mf <- temp_mf[order(-temp_mf$Bitscore, temp_mf$Protein_existence, temp_mf$Evidence),]
          temp_cc <- temp_cc[order(-temp_cc$Bitscore, temp_cc$Protein_existence, temp_cc$Evidence),]
          temp_ec <- temp_ec[order(-temp_ec$Bitscore, temp_ec$Protein_existence),]
        }
        if (align == "global"){
          temp_bp <- temp_bp[order(-temp_bp$Sequence_identity, temp_bp$Protein_existence, temp_bp$Evidence),]
          temp_mf <- temp_mf[order(-temp_bp$Sequence_identity, temp_bp$Protein_existence, temp_bp$Evidence),]
          temp_cc <- temp_cc[order(-temp_bp$Sequence_identity, temp_bp$Protein_existence, temp_bp$Evidence),]
          temp_ec <- temp_ec[order(-temp_ec$Sequence_identity, temp_ec$Protein_existence),]
        }
      }
      # bp
      temp_bp <- temp_bp[!duplicated(temp_bp[, c('Query_ID')]),]
      temp_bp <- merge(temp_bp, goid2term, by="GO_ID", all.x=T)
      temp_bp <- temp_bp[,c(2:19,1,20)]
      colnames(temp_bp)[19] <- "GO_BP_ID"
      colnames(temp_bp)[20] <- "GO_BP_Term"
      # mf
      temp_mf <- temp_mf[!duplicated(temp_mf[, c('Query_ID')]),]
      temp_mf <- temp_mf[,c("Query_ID","UNIPROT_ID", "GO_ID")]
      temp_mf <- merge(temp_mf, goid2term, by="GO_ID")
      temp_mf2 <- temp_mf[,c("Query_ID", "GO_ID", "GO_Term")]
      colnames(temp_mf2) <- c("Query_ID","GO_MF_ID", "GO_MF_Term")
      # cc
      temp_cc <- temp_cc[!duplicated(temp_cc[, c('Query_ID')]),]
      temp_cc <- temp_cc[,c("Query_ID","UNIPROT_ID", "GO_ID")]
      temp_cc <- merge(temp_cc, goid2term, by="GO_ID")
      temp_cc2 <- temp_cc[,c("Query_ID", "GO_ID", "GO_Term")]
      colnames(temp_cc2) <- c("Query_ID","GO_CC_ID", "GO_CC_Term")
      #join all layers
      anota <- merge (temp_bp, temp_mf2, by="Query_ID", all.x=T)
      anota <- merge (anota, temp_cc2, by="Query_ID", all.x=T)
      anota$Evidence <- gsub("xIEA", "IEA", anota$Evidence)
      write.table(anota, "hayai_annotation.csv", row.names=F, col.names=T, sep=",")

      # GO_BP
      temp_bp <- temp_bp [, c("Query_ID", "UNIPROT_ID", "GO_BP_ID")]
      temp_bp <- temp_bp [!is.na(temp_bp$GO_BP_ID), ]
      temp_bp <- temp_bp [, c("Query_ID", "UNIPROT_ID")]
      bp <- merge(temp_bp, bp_table, by ="UNIPROT_ID")
      bp <- bp[,1:3]
      bp <- merge(bp, goid2term, by="GO_ID", all.x=T)
      write.table(bp, "GO_BP_table.csv", row.names=F, col.names=T, sep=",")
      # GO_MF
      temp_mf <- temp_mf [!is.na(temp_mf$GO_ID), ]
      temp_mf <- temp_mf [, c("Query_ID", "UNIPROT_ID")]
      mf <- merge(temp_mf, mf_table, by ="UNIPROT_ID")
      mf <- mf[,1:3]
      mf <- merge(mf, goid2term, by="GO_ID", all.x=T)
      write.table(mf, "GO_MF_table.csv", row.names=F, col.names=T, sep=",")
      # GO_CC
      temp_cc <- temp_cc [!is.na(temp_cc$GO_ID), ]
      temp_cc <- temp_cc[,c("Query_ID", "UNIPROT_ID")]
      cc <- merge(temp_cc, cc_table, by ="UNIPROT_ID")
      cc <- cc[,1:3]
      cc <- merge(cc, goid2term, by="GO_ID", all.x=T)
      write.table(cc, "GO_CC_table.csv", row.names=F, col.names=T, sep=",")
      # EC
      temp_ec <- temp_ec[!duplicated(temp_ec[, c("Query_ID")]), c("Query_ID","EC", "Protein_existence")]
      write.table(temp_ec[,1:3], "EC_table.csv", row.names=F, col.names=T, sep=",")
      length_anota <- dim(anota)[1]

      if (length_anota > 500) { # if number of annotated genes are higher than 500, do graphics
        # Graphics GO_BP
        sum_bp <- sort(table(bp$GO_Term), decreasing=T)
        sum_bp <- as.data.frame(sum_bp)
        sum_bp <- sum_bp[sum_bp$Freq > 0,]
        colnames(sum_bp) <- c("GO_BP_name","Counts")
        sum_bp$GO_BP_name <- strtrim(sum_bp$GO_BP_name, 55)
        bar_bp <- sum_bp[50:1,]
        y <-bar_bp[50,2]
        z <- 1.25*y
        pdf('GO_BP.pdf', width = 11, height = 11 )
        par(mar=c(2,2,6,2), oma=c(0.5,16,5,0.5))
        barplot(bar_bp$Counts, names=bar_bp$GO_BP_name, las=1, horiz=T,  cex.names=0.75, space=2, xlim=c(0,z), axes=F)
        axis(3)
        mtext("Number of Genes", side=3, line=2, font=1)
        title(main = "GO Biological Process - Top 50 - Gene Level", font.main = 4, line=4)
        dev.off()
        write.table(sum_bp, "GO_BP_counts.csv", row.names=F, col.names=T, sep=",")
        # Graphics GO_MF
        sum_mf <- sort(table(mf$GO_Term), decreasing=T)
        sum_mf <- as.data.frame(sum_mf)
        sum_mf <- sum_mf[sum_mf$Freq > 0,]
        colnames(sum_mf) <- c("GO_MF_name","Counts")
        sum_mf$GO_MF_name <- strtrim(sum_mf$GO_MF_name, 55)
        bar_mf <- sum_mf[50:1,]
        y <-bar_mf[50,2]
        z <- 1.25*y
        pdf('GO_MF.pdf', width = 11, height = 11 )
        par(mar=c(2,2,6,2), oma=c(0.5,16,5,0.5))
        barplot(bar_mf$Counts, names=bar_mf$GO_MF_name, las=1, horiz=T,  cex.names=0.75, space=2, xlim=c(0,z), axes=F)
        axis(3)
        mtext("Number of Genes", side=3, line=2, font=1)
        title(main = "GO Molecular Function - Top 50 - Gene Level", font.main = 4, line=4)
        dev.off()
        write.table(sum_mf, "GO_MF_counts.csv", row.names=F, col.names=T, sep=",")
        # Graphics GO_CC
        sum_cc <- sort(table(cc$GO_Term), decreasing=T)
        sum_cc <- as.data.frame(sum_cc)
        sum_cc <- sum_cc[sum_cc$Freq > 0,]
        colnames(sum_cc) <- c("GO_CC_name","Counts")
        sum_cc$GO_CC_name <- strtrim(sum_cc$GO_CC_name, 55)
        bar_cc <- sum_cc[50:1,]
        y <-bar_cc[50,2]
        z <- 1.25*y
        pdf('GO_CC.pdf', width = 11, height = 11 )
        par(mar=c(2,2,6,2), oma=c(0.5,16,5,0.5))
        barplot(bar_cc$Counts, names=bar_cc$GO_CC_name, las=1, horiz=T,  cex.names=0.75, space=2, xlim=c(0,z), axes=F)
        axis(3)
        mtext("Number of Genes", side=3, line=2, font=1)
        title(main = "GO Cellular Component - Top 50 - Gene Level", font.main = 4, line=4)
        dev.off()
        write.table(sum_cc, "GO_CC_counts.csv", row.names=F, col.names=T, sep=",")
        # Graphics EC
        sum_ec <- sort(table(droplevels(temp_ec$EC)), decreasing=T)
        bar_ec <- sum_ec[50:1]
        y <-bar_ec[[50]]
        z <- 1.25*y
        pdf('EC_codes.pdf', width = 11, height = 11 )
        par(mar=c(2,2,6,2), oma=c(0.5,5,5,1))
        barplot(bar_ec, las=1, horiz=T, cex.names=0.75, space=2, xlim=c(0,z), axes=F)
        axis(3)
        mtext("Number of Genes", side=3, line=2, font=1)
        title(main = "Enzyme Commission - Top 50", font = 2, line=4)
        dev.off()
        sum_ec <- as.data.frame(sum_ec)
        colnames(sum_ec) <- (c("EC","Counts"))
        write.table(sum_ec, "EC_counts.csv", row.names=F, col.names=T, sep=",")
        write.table(sum_ec[,1], "unique_EC.csv", row.names=F, col.names=F, quote=F) # use on KEGG Mapper
        } # end of graphics
        anota
    }
  }
  write.table (anota, "hayai_annotation.csv", col.names=F, row.names=F, sep = ",")
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
      write.csv(responses, "hayai_annotation.parameters", row.names=F)
      loadData()

    } else {
      align <<- data[1]
      hits <<- data[2]
      seq_id <<- data[3]
      organism <<- data[4]
      e_value <<- data[5]
      query_cov <<- data[6]
      responses <<- data
      write.csv(responses, "hayai_annotation.parameters", row.names=F)
      loadData()
    }
  }
  if (interactive()) {
    ui <- fluidPage(
      shinyjs::useShinyjs(),
      h1(id="big-heading", "Hayai-Annotation Plants v1.0.2 - Functional Protein Annotation for Plant Species"),
      tags$style(HTML("#big-heading{color: DarkGreen;}")),
      br(),
      h5(id="small-heading", "Kazusa DNA Research Institute - Kisarazu - Japan"),
      tags$style(HTML("#small-heading{color: green;}")),
      br(),
      br(),
      br(),
      fluidRow (
          column (3,
          radioButtons("Align","Type of Alignment",
          choices= list("Local" = "local", "Global" = "global"),
          selected = "local"),
          br(),
          radioButtons("Type_algorithm","Type of Algorithm",
          choices= list("Protein Existence Level" = "protein_existence", "Alignment Score" = "score"),
          selected = "protein_existence")),
          column (3,
          numericInput("hits", "Max hits per query", value = 1, min = 1, max = 20, step = 1),
          br(),
          numericInput("e_value", "Evalue 1e-", value= 6, min = 1, max = 100, step = 1)
          ),
          column (3,
          numericInput("seq_id", "Minimum Sequence Identity (%)", value = 80, min = 40, max = 100, step = 1),
          br(),
          numericInput("query_cov", "Minimum Query Coverage (%)", value= 80, min = 20, max = 100, step = 1))
      ),
      br(),
      br(),
      fileInput("userInput", "Enter Query Sequence (FASTA format)", multiple = F, accept = c(".fasta", ".fa",".faa", ".fna")) ,
      actionButton("submit", "Submit"),
      # wait bottom information
      shinyjs::hidden(p(id = "text1", "Processing...")),
      # end of wait bottom information
      br(),
      br(),
      br(),
      downloadLink("downloadData", "Download"),
      br(),
      br(),
      br(),
      fluidRow(
        DT::dataTableOutput('contents')
      )
    )
    server <- function(input, output, session) {
    # this command controls submit button
    observe({
      shinyjs::toggleState("submit", !is.null(input$userInput) && input$userInput != "")
    })
    observeEvent(input$userInput , {
      inFile <- input$userInput
      if (is.null(inFile))
      return(NULL)
      write.csv(inFile, "df.csv", row.names=F) # changed here add row.namesF
    })
    formData <- reactive({
      data <- sapply(fields, function(x) input[[x]])
      data
    })
    # this starts block bottom action
    plotReady <- reactiveValues(ok = FALSE)
    observeEvent(input$submit, {
      # block submit bottom until annotation if finished
      shinyjs::disable("submit")
      shinyjs::show("text1")
      plotReady$ok <- FALSE
      saveData(formData())
      # unblock submit bottom
      Sys.sleep(1)
      plotReady$ok <- TRUE
    })
    output$contents <- DT::renderDataTable (DT::datatable({
      if (plotReady$ok) {
  #     shinyjs::enable("submit")
        shinyjs::hide("text1")
        input$submit
        anota <- read.table("hayai_annotation.csv", header=T, sep = ",")
        anota <- anota[, c("Query_ID", "Protein_name", "Gene_name", "EC", "GO_BP_ID", "GO_BP_Term", "GO_MF_ID","GO_MF_Term","GO_CC_ID","GO_CC_Term", "Protein_existence", "Sequence_identity", "E_value")]
        anota
      }
    }))
    output$downloadData <- downloadHandler(
        filename = function() {
          paste("output_HayaiAnnotation", "zip", sep=".")
        },
        content = function(fname) {
          fs <- c("hayai_annotation.csv","GO_BP_counts.csv","GO_MF_counts.csv", "GO_CC_counts.csv", "GO_BP.pdf", "GO_MF.pdf", "GO_CC.pdf", "EC_counts.csv", "EC_codes.pdf", "GO_BP_table.csv", "GO_MF_table.csv", "GO_CC_table.csv", "EC_table.csv", "unique_EC.csv")
          zip(zipfile=fname, files=fs)
          },
        contentType = "application/zip"
    )
    }
  shinyApp(ui,server)
  }
}
