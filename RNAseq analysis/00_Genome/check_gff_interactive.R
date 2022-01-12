getwd()

setwd("C:/Users/RaoGrp_Anshu/Box Sync/202002_Methylococcus_RNAseq/Genome/")

library(rtracklayer)

rm(list = ls())

dir()

gff0 <- import("GCF_000008325.1_ASM832v1_genomic.gtf")

table(gff0$type)
# gene         CDS start_codon  stop_codon        exon 
# 3095        3043        3039        3039          56 

gff1 <- gff0[gff0$type == "gene"]

names(mcols(gff1))
# [1] "source"        "type"          "score"         "phase"         "gene_id"       "gbkey"         "gene"          "gene_biotype" 
# [9] "locus_tag"     "old_locus_tag" "transcript_id" "inference"     "product"       "protein_id"    "transl_table"  "pseudo"       
# [17] "note"          "partial"       "anticodon"     "exon_number"   "exception"     "gene_synonym"  "db_xref

mcols(gff1)
# DataFrame with 3095 rows and 23 columns
# source     type     score     phase     gene_id       gbkey        gene   gene_biotype   locus_tag old_locus_tag transcript_id
# <factor> <factor> <numeric> <integer> <character> <character> <character>    <character> <character>   <character>   <character>
#   1      RefSeq     gene        NA        NA MCA_RS00005        Gene        mnmG protein_coding MCA_RS00005       MCA0001            NA
# 2      RefSeq     gene        NA        NA MCA_RS00010        Gene        rsmG protein_coding MCA_RS00010       MCA0002            NA
# 3      RefSeq     gene        NA        NA MCA_RS00015        Gene          NA protein_coding MCA_RS00015       MCA0003            NA
# 4      RefSeq     gene        NA        NA MCA_RS00020        Gene          NA protein_coding MCA_RS00020       MCA0004            NA
# 5      RefSeq     gene        NA        NA MCA_RS00025        Gene          NA protein_coding MCA_RS00025       MCA0005            NA

length(unique(gff1$gene))
# 698

length(unique(gff1$gene_id))
# 3095

length(unique(gff1$transcript_id))
# 1

# If there doesn't appear to be an attribute name that contains good 
# gene-level IDs, then you need to do more testing using the S4_objects_29mar19.R
# code from here: https://go.illinois.edu/introR  inside the introBioC_29mar19.zip file

# If there is a good attribute name and the format was already .gft, no 
# modification is needed and you can simply exit R by typing q() at the prompt.

# If there is a good attribute name but the format was gff, you can quickly
# write out the original gff0 object to gtf format by:
export(gff0, "desired_new_name.gtf", format = "gtf")

# Don't forget to document how you changed the gene feature file!! To quit
# interactive R and go back to the biocluster command prompt, type in:
q()
  #when asked to save workspace, select no!







