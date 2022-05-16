######R script to generate a table summarizing EVE annotations (position in the genomes of interest, evidence for taxonomy assignment, to be checked by eye)

#before running the R script, please read and follow instructions on the taxonomizr webpage to install appropriate databases: https://github.com/sherrillmix/taxonomizr

#To obtain the "outputDiamondblastX_EVE_candidate.fas_Versus_nr" use the script available here https://github.com/CaroleBelliardo/EVE_annotation.git/EVE_annotation.sh 

setwd("location_directory")

library(data.table)
library(taxonomizr)

blast=fread("outputDiamondblastX_EVE_candidate.fas_Versus_nr", sep = "\t", na.strings=c(""," ","NA"))

ids<-accessionToTaxa(blast$V2,"/path/to/accessionTaxa_prot4thfeb2021.sql")
taxo=getTaxonomy(ids,"/path/to/accesion/accessionTaxa_prot4thfeb2021.sql")

blast=cbind(blast, taxo)
VirFamToDelete=c("Baculoviridae", "Metaviridae", "Nudiviridae", "Phycodnaviridae", "Poxviridae", "Mimiviridae", "Adintoviridae", "Iridoviridae", "Polydnaviridae", "Marseilleviridae", "Caudovirales", "Siphoviridae", "Myoviridae", "Podoviridae")

blast=blast[!(blast$family %in% VirFamToDelete)]
blast=blast[!(blast$order %in% VirFamToDelete)]

setorder(blast, V1, -V12)

virusHits=blast[blast$superkingdom=="Viruses"]
nonVirusHits=blast[blast$superkingdom!="Viruses"]
NAhits=blast[is.na(blast$superkingdom)]
nonVirusHits=rbind(nonVirusHits, NAhits)

tab1=table(nonVirusHits$V1)
tab2=table(virusHits$V1)
tab1=as.data.table(tab1)
tab2=as.data.table(tab2)
tab2$V=tab1$N[match(tab2$V1, tab1$V1)]
tab1$V=tab2$N[match(tab1$V1, tab2$V1)]
tab2$X=tab2$N/(tab2$V+tab2$N)*100
setorder(tab2, -X)

setorder(blast, V1, -V12)
bestHits=blast[!duplicated(blast$V1)]

tab2$bestHit_AccNum=bestHits$V2[match(tab2$V1, bestHits$V1)]
tab2$bestHit_SuperKingdom=bestHits$superkingdom[match(tab2$V1, bestHits$V1)]
tab2$bestHit_species=bestHits$species[match(tab2$V1, bestHits$V1)]
tab2$bestHit_Fam=bestHits$family[match(tab2$V1, bestHits$V1)]

setorder(tab2, -X)

setorder(nonVirusHits, V1, -V12)
head(nonVirusHits)
nonVirusHits_BestHits=nonVirusHits[!duplicated(nonVirusHits$V1)]

setorder(virusHits, V1, -V12)
virusHits_BestHits=virusHits[!duplicated(virusHits$V1)]

tab2$aLength=bestHits$V4[match(tab2$V1, bestHits$V1)]
tab2$scoreVirusHit=virusHits_BestHits$V12[match(tab2$V1, virusHits_BestHits$V1)]
tab2$scoreNonVirusHit=nonVirusHits_BestHits$V12[match(tab2$V1, nonVirusHits_BestHits$V1)]

tab2$ID_VirusHit=virusHits_BestHits$V3[match(tab2$V1, virusHits_BestHits$V1)]
tab2$ID_NonVirusHit=nonVirusHits_BestHits$V3[match(tab2$V1, nonVirusHits_BestHits$V1)]

tab2$evalueVirusHit=virusHits_BestHits$V11[match(tab2$V1, virusHits_BestHits$V1)]

tab2$VirusBestHit=virusHits_BestHits$V13[match(tab2$V1, virusHits_BestHits$V1)]
tab2$nonVirusBesstHit=nonVirusHits_BestHits$V13[match(tab2$V1, nonVirusHits_BestHits$V1)]

tab2$scoreDiff=tab2$scoreVirusHit-tab2$scoreNonVirusHit

tab2=tab2[tab2$evalueVirusHit<0.0001]

tab2[, c("chr", "coord") := tstrsplit(tab2$V1, ":", fixed=T)]

tab2[, c("start", "end") := tstrsplit(tab2$coord, "-", fixed=T)]

tab2$EVElength=as.numeric(tab2$end)-as.numeric(tab2$start)

tab2$percQuerAln=(100*(tab2$aLength*3))/tab2$EVElength


tab2$bestHitLineage=bestHits$V14[match(tab2$V1, bestHits$V1)]

tab2$V1=NULL

tab2$coord=NULL

tab2=tab2[,c("chr", "start", "end", "EVElength", "percQuerAln", "N", "V", "X", "bestHit_AccNum", "bestHit_SuperKingdom", "bestHit_species", "bestHit_Fam", "aLength", "scoreVirusHit", "scoreNonVirusHit", "ID_VirusHit", "ID_NonVirusHit", "evalueVirusHit", "scoreDiff", "bestHitLineage", "VirusBestHit", "nonVirusBesstHit")]

tab2=tab2[,c("chr", "start", "end", "EVElength", "percQuerAln", "N", "V", "X", "bestHit_AccNum", "bestHit_SuperKingdom", "bestHit_species", "bestHit_Fam", "aLength", "scoreVirusHit", "scoreNonVirusHit", "ID_VirusHit", "ID_NonVirusHit", "evalueVirusHit", "scoreDiff", "VirusBestHit", "nonVirusBesstHit")]

names(tab2)[names(tab2) == 'N'] <- 'nbVirHits'

names(tab2)[names(tab2) == 'V'] <- 'nbNonVirHits'

names(tab2)[names(tab2) == 'X'] <- 'percVirHits'

setorder(tab2, -bestHit_SuperKingdom, -percVirHits)

write.table(tab2, "EVEannotation_species_name.txt", quote = F, sep = "\t", row.names = F)


q()

