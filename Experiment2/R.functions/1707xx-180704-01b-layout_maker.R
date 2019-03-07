library(ggplot2)
library(reshape)
library(reshape2)
library(dplyr)
library(Hmisc)
library(plyr)



head.dir<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis'
head.dir.sys<-"/Users/anew/Downloads"
layout.dir<-'./layout.sample.data'
setwd(head.dir)


# make description of source plates from ordered clone name description file. 
# This file is a list of 144 items, composed of 3 sets of 48 allele / clone names of GAL3, 80 and 4 in that order.
clone.names<-read.table(paste0(layout.dir,'/1707xx-CloneNames.txt'),header=F)
GAL3<-paste(droplevels(clone.names[1:48,]))
GAL80<-paste(droplevels(clone.names[49:96,]))
GAL4<-paste(droplevels(clone.names[97:144,]))

first.set<-data.frame(GAL3.allele= paste(rep(GAL3,48)),
	GAL80.allele= paste(as.vector(sapply(GAL80,function(x){rep(x,48)}))),
	GAL4.allele= paste(rep('GAL4.WT',48*48)))
	
second.set<-data.frame(GAL3.allele= paste(rep('GAL3.WT',48*48)),
	GAL80.allele= paste(as.vector(sapply(GAL80,function(x){rep(x,48)}))),
	GAL4.allele= paste(rep(GAL4,48)))

third.set<-data.frame(GAL3.allele= paste(as.vector(sapply(GAL3,function(x){rep(x,48)}))),
	GAL80.allele= paste(rep('GAL80.WT',48*48)),
	GAL4.allele= paste(rep(GAL4,48)))


source.plate.genotypes<-as.data.frame(rbind(first.set,second.set,third.set))


write.table(source.plate.genotypes,paste0(layout.dir,'/1707xx-SourcePlateGenotypes.txt'),sep='\t',row.names=FALSE)



source.date<-170600 # date plates labeled in freezer
date<-'1707xx' # beginning date clones were measured

columns<-sapply(1:12,function(x){rep(x,8)})
rows<-t(sapply(1:8,function(x){rep(LETTERS[x],12)}))

rc<-data.frame(r=melt(t(rows))$value,c=melt(t(columns))$value)
rc$rc<-sapply((apply(rc,1,function(x){paste(x[1],x[2],sep='')})),function(x){gsub(" ","",x)})

# n is number of source plates
n <- 72
source.plates<-data.table()
for(i in 1:n){
	if(i<10){source.plate<-rep(paste("Plate_00",i,sep=''),96)}else(source.plate<-rep(paste("Plate_0",i,sep=''),96))
	platei<-data.frame(source.plate,rc)
	source.plates<-rbind(source.plates,platei)
	
}
source.plates[,source.plate.well:=apply(data.frame(source.plates$source.plate,source.plates$rc),1,paste,collapse='.')]

source.plate.merge1<-data.table(source.plate.genotypes, source.plate.well=source.plates$source.plate.well)
source.plate.merge1$background.strain<-'AN612'



# add background strain information
background.strain.info<-fread(paste0(layout.dir,'/1707xx-BackgroundStrainInfo.txt'),sep='\t',header=T)
source.plate.merge<-merge(source.plate.merge1, background.strain.info,by='background.strain')
table(source.plate.merge $background.strain)

# now begin merging the experimental plates with the source plates
# experimental plates are described in the PlateSetSummary.txt file
# this file has rows corresponding to each plate.
# columns describe the general features of the plate, including the source plate for the experimental plate, 
# media batch, date of measurement, as well as the folder name for where the plates with data are stored.

# all plates with data have been named according to the day they were measured. generally 24 plates were measured
# per day, therefore plates range from _001 - _024, so a typical name for the plate would be 170707-Plate_001

PlateSetSummary<-fread(paste0(layout.dir,'/1707xx-PlateSetSummary.txt'),header=T,sep='\t')


PlateSetSummary$samples.to.mask<-as.logical(!(as.logical(is.na(PlateSetSummary$mask.scoring) * is.na(PlateSetSummary$mask.transformation))))


head(PlateSetSummary)
#PlateDatRepeater = PDR . This repeats the information for each plate 96 times, corresponding to each row and column of the plate.
PDR <-function(x,nrow,ncol=96){ 
	# because the data frame gets translated the cols are the rows and the rows are the cols. sorry for the fucked up naming.
	as.data.frame(t(matrix(x,nrow=nrow,ncol=ncol)))}

numberrows<-ncol(PlateSetSummary)


ExpPlates1<-apply(PlateSetSummary,1,function(x){
	xa<-paste(sapply(x,as.character))
	
	# this stupid shit was the first way i could find for the apply fctn to interpret the TRUE or FALSES in the x variable as a logical
	stm<-unlist(as.vector(as.data.frame(as.character(x['samples.to.mask']))))
	ms<-unlist(as.vector(as.data.frame(as.character(x['mask.scoring']))))
	mt<-unlist(as.vector(as.data.frame(as.character(x['mask.transformation']))))

# TRYING TO FIGURE OUT HOW TO LABEL WELLS THAT NEED MASKING FROM ANALYSIS
# i think this works
	if(stm){
		ya<-as.vector(ms)
		za<-as.vector(mt)
		a<-!is.na(ya)
			if(a){yya<-strsplit(ya,"[ ]")[[1]]}else(yya<-NA)
		b<-!is.na(za)
			if(b){zza<-strsplit(za,"[ ]")[[1]]}else(zza<-NA)
		ab<-as.vector(unique(na.omit(c(zza,yya))))
		match(ab,rc$rc,)
		keep<-!(rc$rc %in% ab)		
	}else(keep<-rep(FALSE,96))
	df<-data.frame(PDR(xa,numberrows,ncol=96))
	df$keep<-keep
	df
})

ExpPlates2 <- data.table(ldply(ExpPlates1, data.frame))
colnames(ExpPlates2)<-c(colnames(PlateSetSummary),'keep')

ExperimentalPlates<-data.table(rc,ExpPlates2)
ExperimentalPlates

ExperimentalPlates[, source.plate.well:=apply(data.frame(ExperimentalPlates $source.plate, ExperimentalPlates $rc),1,paste,collapse='.')]
ExperimentalPlates$samp_id<-paste(ExperimentalPlates$measurement.date,"-",ExperimentalPlates$plate.name.in.source.fcs,".",ExperimentalPlates$rc,sep='')

# everything seems to work :)
merge(ExperimentalPlates, source.plate.merge, by='source.plate.well',all=TRUE)

# still need a way to mask samples that were problematic in the cloning / inoculation etc processes.
# this will come from another source file.

layout.expt1<-merge(ExperimentalPlates, source.plate.merge, by='source.plate.well',all.x=TRUE)
layout.expt<-layout.expt1[!is.na(layout.expt1$background.strain )]
sampling.data1<-fread('1707xx-SamplingData.txt',header=T,sep='\t')
a<-sampling.data1[grep('1709', sampling.data1$samp_id)]
b<-a[grep('Plate_0', a$samp_id)]
bglu<-b[1:(nrow(b)/2)]
bgal<-b[(nrow(b)/2+1):nrow(b)]
bgal$samp_id<-gsub('170927','170928',bgal$samp_id)
pre<-sampling.data1[!grepl('1709', sampling.data1$samp_id)==TRUE]
sampling.data<-rbind(pre,bglu,bgal)

# METADATA<-fread('1707xx-METADATA-ProcessedFromFCS.txt',header=T,sep='\t')
allele.description<-fread(paste0(layout.dir,"/1707xx-BigMixNMatch-AlleleDescriptions.txt"),sep='\t',header=T)
# test<-merge(sampling.data,METADATA,'samp_id',all=T)
# test[grep('1709',samp_id)]$sample_flow_rate
# test2<-merge(test,layout.expt,'samp_id',all=T)
# test2[grep('1709',samp_id)]$sample_flow_rate
# test2[is.na(sample_flow_rate)]$samp_id


together<-data.table(merge(layout.expt,sampling.data,'samp_id',all=T))
# together<-data.table(merge(together2,METADATA,by='samp_id',all=T))
together[grep('1709',samp_id)]$sample_flow_rate

colnames(allele.description)
keep.cols<-c('allele',"clone.named", "classic.allele","mutation.nickname","mutation.AA","mutant.class")

GAL3.descrip<-allele.description[!is.na(GAL3.allele)]
GAL3.descrip<-GAL3.descrip[,keep.cols,with=FALSE]
setnames(GAL3.descrip,colnames(GAL3.descrip),paste("GAL3",colnames(GAL3.descrip),sep='.'))

GAL80.descrip<-allele.description[!is.na(GAL80.allele)]
GAL80.descrip<-GAL80.descrip[,keep.cols,with=FALSE]
setnames(GAL80.descrip,colnames(GAL80.descrip),paste("GAL80",colnames(GAL80.descrip),sep='.'))

GAL4.descrip<-allele.description[!is.na(GAL4.allele)]
GAL4.descrip<-GAL4.descrip[,keep.cols,with=FALSE]
setnames(GAL4.descrip,colnames(GAL4.descrip),paste("GAL4",colnames(GAL4.descrip),sep='.'))

layout<-merge(merge(merge(together,GAL3.descrip,by='GAL3.allele',all=T),GAL80.descrip,by='GAL80.allele',all=T),GAL4.descrip,by='GAL4.allele',all=T)

write.table(layout,paste0(layout.dir,"/1707xx-180704-layout.txt"),row.names=FALSE,sep='\t')

