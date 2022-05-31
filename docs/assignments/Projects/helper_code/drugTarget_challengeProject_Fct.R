#######################################################################
## Drug-target analysis of proteins encoded by genes in peak regions ##
#######################################################################

## The following should be executed under highest level of the ChIP-Seq
## workflow project. In addition, a subdirectory named 'drug_target_challenge_project'
## needs to be created. Most of the challenge project results will be written
## to this sub-directory.


## Load peak annotation table
library(systemPipeR)
sal <- SPRproject(resume=TRUE)
annofiles <- getColumn(sal, step = "annotation_ChIPseeker", "outfiles") # requires sal object; path in annofiles is: ./results/AP1_1_ChIPseeker_annotated.xls
peak_annot <- read.delim(annofiles) 

## Identify genes with peaks within their <=1kb promoter region
## Note, second line subsets to peaks that are fully contained in the 1kb promoter regions!
peak_annot_sub <- peak_annot[peak_annot$annotation=="Promoter (<=1kb)",]
peak_annot_sub <- peak_annot_sub[loc_within <- (peak_annot_sub$start >= peak_annot_sub$geneStart) & (peak_annot_sub$end <= peak_annot_sub$geneEnd), ]
geneIDs <- unique(peak_annot_sub$geneId)
length(geneIDs) # final result contains 289 unique genes 

## Obtain for geneIDs the corresponding protein sequences from Arabidopsis
library(GenomicFeatures); library(Biostrings); library(Rsamtools); library(BSgenome.Athaliana.TAIR.TAIR9)
txdb <- makeTxDbFromGFF(file = "data/tair10.gff", format = "gff", dataSource = "TAIR", organism = "Arabidopsis thaliana")
genome <- BSgenome.Athaliana.TAIR.TAIR9
cds <- cdsBy(txdb, by="tx", use.names=TRUE)
cds_seqs <- extractTranscriptSeqs(genome, cds)
protein_seqs <- translate(cds_seqs, if.fuzzy.codon="solve")
subset_index <- gsub("\\..*", "", names(protein_seqs)) %in% geneIDs
protein_seqs <- protein_seqs[subset_index]
writeXStringSet(protein_seqs, "drug_target_challenge_project/Arab_pep.fasta")

## Obtain human protein sequences
library(ensembldb); library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
prts <- proteins(edb, return.type = "AAStringSet")
prts <- prts[!duplicated(names(prts))]
writeXStringSet(prts, "drug_target_challenge_project/hs_pep.fasta")

## Identify human orthologs via BLASTP
mydir <- getwd()
setwd("drug_target_challenge_project")
moduleload("ncbi-blast/2.2.31+")
system("makeblastdb -in hs_pep.fasta -out hs_pep.fasta -dbtype prot -hash_index -parse_seqids")
system("blastp -query Arab_pep.fasta -db hs_pep.fasta -outfmt 6 -evalue 1e-6 -out blastp.tab")
setwd(mydir)

## Perform drug-target analysis
library(drugTargetInteractions); library(ensembldb); library(EnsDb.Hsapiens.v86)
blast_result <- read.delim("drug_target_challenge_project/blastp.tab", header=FALSE)
blast_result <- blast_result[!duplicated(as.character(blast_result[,1])),] # keeps for each query only best scoring (first) result
blast_result_filtered <- blast_result[blast_result[,11] <= 1e-10,]
id_mapping <- select(edb, keys = blast_result_filtered[,2], keytype = "PROTEINID", columns = c("GENEID", "GENENAME", "UNIPROTID"))
uniprot_query_ids <- id_mapping$UNIPROTID



