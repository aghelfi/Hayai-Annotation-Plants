# Hayai-Annotation Plants v1.0.2

R-package for an Ultra-Fast and Comprehensive Gene Annotation in Plants

Description
-----------
Hayai-Annotation-Plants is an R-package with a shiny-browser interface for a highly accurate and fast gene annotation system for plant species (Embryophyta). 

Hayai-Annotation-Plants provides five levels of annotation: 

1) Protein Name; 
2) Gene Ontology (GO) consisting of three main domains (Biological Process, Molecular Function and Cellular Component); 
3) [Enzyme Commission Code - EC](https://enzyme.expasy.org/); 
4) [Protein Existence Level](https://www.uniprot.org/help/protein_existence); 
5) [Evidence Code for GO annotation](http://geneontology.org/docs/guide-go-evidence-codes/).


Documentation
-------------
For downloading, installing and running **Hayai-Annotation Plants** please see [Hayai-Annotation-Plants wiki](https://github.com/kdri-genomics/Hayai-Annotation-Plants/wiki) 

Reference
---------
Hayai-Annotation Plants: an ultra-fast and comprehensive functional gene annotation system in plants <br/>
Andrea Ghelfi, Kenta Shirasawa, Hideki Hirakawa, Sachiko Isobe <br/>
Bioinformatics, btz380, https://doi.org/10.1093/bioinformatics/btz380 

Updates
-------
2019/07/05. Update database on 2019/06/21.<br/>
Taxonomy Viridiplantae were implemented, now it includes algae. <br/>
Database were clustered to 99% sequence identity and a new algorithm were implemented to make Hayai-annotation faster <br/>

2019/06/20. [Hayai-annotation Plants on Docker Container beta version](https://hub.docker.com/r/kazusa005/hayai-annotation-plants). Recommended for up to 2,000 sequences. <br/>

2019/01/04. Low number of input sequences crash - Fixed. Now, graphics are generated if number of annotated queries are higher than 500. <br/>
Hayai-Annotation Plants run with even if USEARCH is installed in another directory (requires symbolic link).
