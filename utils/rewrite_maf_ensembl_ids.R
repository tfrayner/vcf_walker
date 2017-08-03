#!/usr/bin/env Rscript
##
## Usage: Rscript rewrite_maf_ensembl_ids.R <MAF filename> [Ensembl version (optional)]

library(biomaRt)

args <- commandArgs(trailingOnly=TRUE)

if ( length(args) > 1 ) {
  mart <- useEnsembl(version=as.numeric(args[2]), biomart='ensembl', dataset='mmusculus_gene_ensembl')
} else {
  mart <- useEnsembl(biomart='ensembl', dataset='mmusculus_gene_ensembl')
}

message("Reading in MAF file...")
maf <- read.table(args[1], sep="\t", header=TRUE, stringsAsFactors=FALSE)

message("Reannotating MAF file...")
## We can't really handle transcripts mapping to multiple Entrez IDs;
## here we just take the first Entrez ID in the file.
ann <- getBM(mart=mart, filter='ensembl_transcript_id', attr=c('entrezgene', 'external_gene_name', 'ensembl_transcript_id'), values=maf$Entrez_Gene_Id)
ann <- ann[ ! duplicated( ann$ensembl_transcript_id ), ]
rownames(ann) <- ann$ensembl_transcript_id

eg <- ann[ maf$Hugo_Symbol, 'entrezgene' ]

maf$Hugo_Symbol    <- ann[ maf$Hugo_Symbol, 'external_gene_name' ]
maf$Hugo_Symbol[is.na(eg)] <- 'Unknown'

maf$Entrez_Gene_Id <- eg
maf$Entrez_Gene_Id[is.na(eg)] <- 0

message("Writing out new MAF file...")
write.table(maf, file=paste(args[1], 'fixed', sep='.'), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

message("Done.")
