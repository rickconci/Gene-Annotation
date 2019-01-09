####################################################################################
#                                      Functional domains
####################################################################################

setwd("~/Desktop/code/CompBio MPhil/GenomeInformatics/GIa3/Function and sequences/Function/")


human.sry.ppt <- ">SRY-201 peptide: ENSP00000372547 pep:protein_coding
MQSYASAMLSVFNSDDYSPAVQENIPALRRSSSFLCTESCNSKYQCETGENSKGNVQDRVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPKNCSLLPADPASVLCSEVQLDNRLYRDDCTKATHSRMEHQLGHLPPINAASSPQQRDRYSHWTKL"

human.sry.ppt.only <- strsplit(human.sry.ppt, split="\n")[[1]][2]

human.sry.ppt.NTD <- substr(human.sry.ppt.only, 0, 59)
human.sry.ppt.CTD <- substr(human.sry.ppt.only, 129, 204)

write.table(human.sry.ppt.NTD, "human.sry.ppt.NTD.txt", 
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(human.sry.ppt.CTD, "human.sry.ppt.CTD.txt", 
            row.names = FALSE, quote = FALSE, col.names = FALSE)



#### remove duplicated entries on full fasta sequences
setwd("~/Desktop/code/CompBio MPhil/GenomeInformatics/GIa3/
      Function and sequences/Phylogeny/")

library("seqinr")
homologues.fasta <- read.fasta("sryblastfullfasta")
hits <- names(homologues.fasta)
length(hits)
length(unique(hits))

#select top 200 fasta sequences
twohundred.homologues.fasta <- homologues.fasta[1:200]
length(unique(names(twohundred.homologues.fasta))) #check if duplicates
#save each fasta as 1 string
names(twohundred.homologues.fasta) <- paste0(">", names(twohundred.homologues.fasta))
for (name in names(twohundred.homologues.fasta)) {
  twohundred.homologues.fasta[[name]] <- toupper(paste(
    twohundred.homologues.fasta[[name]], collapse=""))
  twohundred.homologues.fasta[[name]] <- paste(
    name, twohundred.homologues.fasta[[name]], sep ="\n" )
}

twohundred.homologues.fasta <- unlist(twohundred.homologues.fasta)

write.table(twohundred.homologues.fasta, "twohundred.homologues.fasta.txt",sep="\n", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


####################################################################################
#                                      Variation
####################################################################################

setwd("~/Desktop/code/CompBio MPhil/GenomeInformatics/GIa3/Variation")

#from https://www.ncbi.nlm.nih.gov/gene/6736 
sry.start.GRCh37 <- 2654896
sry.start.GRCh38 <- 2786855
HMGbox.start <- 60*3 #starts at 60th AA * 3 nt per amino acid
HMGbox.end <- 128*3 
SRY.length <- 204*3

#downloaded from ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/chr_rpts/
#chr.Y <- read.delim("chr_Y.txt", header=TRUE, stringsAsFactors = FALSE)
#sry.snps <- chr.Y[grepl('SRY', chr.Y$local), ]
#sry.snps.select <- sry.snps[, c(12, 17, 24, 25)]
#sry.snps.select$chr.2 <- as.numeric(sry.snps.select$chr.2) - sry.start.GRCh38
#sry.snps.select <- sry.snps.select[sry.snps.select$vali. >0, ]

#downloaded from http://gnomad.broadinstitute.org/gene/ENSG00000184895
# gnomad data: GRCh37, not GRCh38... 
gnomad.sry.snps <- read.csv("gnomAD_v2.1_ENSG00000184895_2018_12_18_21_58_03.csv",
                            header=TRUE)
gnomad.sry.snps$Position <- gnomad.sry.snps$Position - sry.start.GRCh37

#polyphen compares to GRCh37
polyphen.input <- gnomad.sry.snps$rsID
write.table(polyphen.input,"polyphen.input.txt",sep="\n",
            row.names=FALSE, quote = FALSE)


#ensembl data: GRCh38 
#downloaded from https://www.ensembl.org/Homo_sapiens/Gene/Variation_Gene/...
#   Table?db=core;g=ENSG00000184895;r=Y:2786855-2787699;t=ENST00000383070
ensembl.sry.snp <- read.csv("ensembl-gene-variations.csv", header = TRUE)
ensembl.sry.snp <- ensembl.sry.snp[grepl("SNP", ensembl.sry.snp$Class), ]
ensembl.sry.snp <- ensembl.sry.snp[order(ensembl.sry.snp$sift_sort,
                                         decreasing = FALSE), ]
ensembl.sry.snp <- ensembl.sry.snp[grepl('missense', ensembl.sry.snp$Conseq..Type), ]
ensembl.sry.snp$Location <- gsub("Y:", "", ensembl.sry.snp$Location)

#VEP compares to GRCh38
vep.snps <- cbind("Y", ensembl.sry.snp$Location, ensembl.sry.snp$Location,
                  as.character(ensembl.sry.snp$Alleles), "+")
write.table(vep.snps,"VEP.input.snps.txt",sep="\t",row.names=FALSE,
            quote = FALSE, col.names = FALSE)

ensembl.sry.snp$Location <- as.numeric(ensembl.sry.snp$Location) - sry.start.GRCh38


# VEP results

vep.results <- read.delim("VEPresults.txt", header=TRUE)
vep.results.select <- vep.results[grep("protein_coding", vep.results$BIOTYPE), ]
vep.results.select <- vep.results.select[,c(1,4,17,18, 20, 28, 29, 47)]

write.csv(vep.results.select, "vep.results.select.csv", quote=FALSE, 
          col.names = FALSE, row.names = FALSE)


#convert correct PPT sequence to mutated one
locations <- vep.results.select$Protein_position
alleles <- vep.results.select$Amino_acids
for ( i in 1:length(locations)){
  substr(human.sry.ppt.only, locations[i], locations[i]) <- strsplit(as.character(alleles),
                                                                     split="/")[[i]][2]
}

substr(human.sry.ppt.only, 155, 157 )



####################################################################################
#                                      Expression
####################################################################################
setwd("~/Desktop/code/CompBio MPhil/GenomeInformatics/GIa3/Expression/")

#downloaded from Encode
rna.tissues <- read.delim("rna_tissue.tsv")
rna.tissues.sry <- rna.tissues[grep("SRY", rna.tissues$Gene.name), ]
rna.tissues.sry <- rna.tissues.sry[order(rna.tissues.sry$Value, decreasing = TRUE),]

library("ggplot2")
ggplot(data=rna.tissues.sry, aes(x = reorder(Sample, -Value), y = Value, fill = Gene.name)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
  
  


