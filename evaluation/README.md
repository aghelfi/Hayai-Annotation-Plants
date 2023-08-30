## Benchmarking Hayai-Annotation Plants: A Re-evaluation Using Standard Evaluation Metrics 

Prepare gold standard dataset using goa_arabidopsis.gaf

#### Extract gene ids list from Araport11 
egrep ">" dataset/Araport11_genes.201606.pep.repr.fasta | awk '{sub(">",""); print $1}' | awk -F"." '{print $1}' > dataset/araport11_list.txt

#### Create a tsv file with gene_id go_id pairs for each gene from Araport11 gene list
python3 scripts/gene2go.py

#### Generate information accretion file using R-packages SemDist and org.At.tair.db
Rscript scripts/make_ia.R

#### Clone CAFA-evaluator "https://github.com/BioComputingUP/CAFA-evaluator"

cafa_dir="CAFA-evaluator/src"

python3 ${cafa_dir}/main.py dataset/go-basic_2019Jul.obo prediction_dir/ dataset/gold_standards.tsv -out_dir results/ -ia dataset/IA.tsv -prop fill -norm cafa -th_step 0.001 -max_terms 500 -no_orphans -threads 8

#### Generate Graphics 
jupyter notebook ${cafa_dir}/plot.ipynb
