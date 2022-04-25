#!/usr/bin/Rscript
args <- commandArgs(TRUE)
if (length(args) == 0){
 print ('Rscript --vanilla ../absolutPosi.R inFile_path outDirectory_path extension_string') }
#if (length(args) == 1) {
# args[3] <- '.absP'}

df=read.table(args[1] , stringsAsFactors=F, sep='\t')
names(df)<-c("Contig", "start", "stop", "ligneage","stranded",'a','b','c')

#splits
infoC= strsplit(df$Contig,':')
df.infoC <- data.frame(matrix(unlist(infoC), nrow=length(infoC), byrow=T))
names(df.infoC)<- c('Contig','posi')
posiC= strsplit(as.character(df.infoC$posi),'-')
df.posiC <- data.frame(matrix(unlist(posiC), nrow=length(posiC), byrow=T))
names(df.posiC)<- c('start','stop')
# concatenate all info
info_Posi= data.frame(df.infoC$Contig,  strtoi(df.posiC$start),  strtoi(df.posiC$stop ))
names(info_Posi)<- c('Contig','start','stop')
starts=(info_Posi$start+df$start)
stops=(info_Posi$start+df$stop)
df_F= data.frame(info_Posi$Contig,starts,stops,df$ligneage,df$stranded,df$a,df$b,df$c)
print(head(df_F))
outt=paste(args[2],'/',args[1],sep='')
write.table(df_F, file=outt, row.names=FALSE, quote=FALSE, col.names=FALSE, sep='\t')
