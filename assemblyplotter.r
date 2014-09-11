#!/usr/bin/Rscript

#install.packages("seqinr")
library(seqinr)
###########
#file list#
###########
files=list.files(pattern=".fa")
#################
#load assemblies#
#################
writeLines("Loading assemblies")
assemblies=NULL
for (i in files){
  assemblies[[i]]<-read.fasta(i)
}
#######################
#contiglength function#
#######################
contiglength=function(nr,data){
  as.numeric(lapply(data[[nr]],length))
}
##########################
#assemblybarplot function#
##########################
assemblybarplot=function(a){
  x=list()
  for(i in 1:length(a)){
    x[[i]]=sort(contiglength(i,a),decreasing=T)
  }
  xlength=lapply(x,length)
  maxxlength=max(unlist(lapply(x,length)))
  y=matrix(nrow=length(a),ncol=maxxlength)
  for(i in 1:length(a)){
    x[[i]]=c(x[[i]],rep(0,maxxlength-as.numeric(xlength[i])))
    y[i,]=sort(as.numeric(x[[i]]))
  }
  #barplot(t(y),ylab="Assemblylength (bp)",xlab="Specimen",horiz=T)
  barplot(t(y),names.arg=1:length(assemblies),ylab="Assemblylength (bp)",xlab="Specimen")
}
###################################
#functions to remove short contigs#
###################################
min.length=1000
removeshortcontigs=function(nr,threshold,data){
  remove=which(contiglength(nr,data)<threshold)
  if(length(remove)==0){
    return(data[[nr]])
  }
  else{
    return(data[[nr]][-remove])
  }
}
corrected=NULL
for(i in 1:length(assemblies)){
  corrected[[i]]=removeshortcontigs(i,min.length,assemblies)
}
names(corrected)=names(assemblies)
for(i in 1:length(names(corrected))){
  names(corrected)[i]=paste(names(corrected)[i],"corrected",sep="_")
}
##########
#plotting#
##########
writeLines("Generating plot")
jpeg("plot",width=length(corrected)*60,height=500)
par(mfrow=c(1,2))
assemblybarplot(corrected)
assemblybarplot(assemblies)
dev.off()
#############################
#saving corrected assemblies#
#############################
writeLines("Saving corrected assemblies")
for(i in 1:length(corrected)){
  correctedname=paste0(strsplit(names(corrected)[i],".fa")[[1]],"_corrected.fa")
  write.fasta(corrected[[i]],names=names(corrected[[i]]),file.out=correctedname[1],open="w")
}

write.table(cbind(1:length(assemblies),names(assemblies)),file="plot_legend.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

