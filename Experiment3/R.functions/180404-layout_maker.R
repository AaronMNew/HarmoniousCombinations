# source plates were made by hand.

source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'
source(source.code.path.file)
load.all() # this function loads all the standard packages you use in this script

head.dir<-"/Users/anew/Google Drive/002_lehnerlabData/180404-VarGALK.With.G80xG4_2"
date<-'180404' # beginning date clones were measured
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/'
outputdata<-'./output.data/'
layoutdir<-'./layout.sample.data/'
setwd(head.dir)


###################################################
### make source plates
###################################################

sourceplate1<-'180324-SourcePlateLayoutPre.txt'

SP1<-fread(paste0(layoutdir,sourceplate1),header=T,sep='\t')
dim(SP1)

# contamsamps<-fread(paste0(layoutdir,'180302-ContaminatedWells.txt'))
# x<-contamsamps[1,c('contaminated.wells','low.tx.efficiency','contamination.second.score')]
# masksamp1<-data.table(ldply(apply(contamsamps,1,function(xd){
	# x<-xd[2:length(xd)]
	# xa<-as.character(strsplit(paste0(c(x[1],x[2],x[3]),collapse=' '),'\ ')[[1]])
	# DT<-data.table(source.plate=xd[1],well.to.mask=unique(xa[!grepl('NA',xa)]))
	# DT#paste0(DT,collapse='.')
# })))
# masksamp<-apply(masksamp1,1,function(x)paste(x,collapse='.'))
# keepsamp<-SP1$clone.id%w/o%masksamp
# sourceplate.layout<-SP1[clone.id%in%keepsamp]
sourceplate.layout<-SP1

write.table(sourceplate.layout,paste0(layoutdir,date,'-sourceplate.layout.txt'),row.names=F,sep='\t')


###################################################
### make sampling layouts
###################################################

sourceplate.layout<-fread(paste0(layoutdir,date,'-sourceplate.layout.txt'),header=T,sep='\t')
samplingdata<-fread(paste0(layoutdir,date,'-SamplingData.txt'),header=T,sep='\t')
# for some reason plate 24 is not getting exported. On 140404 sampling rate was 0.5 µl/s and on 140405 was 2.0 µl/s
Plate_24<-samplingdata[plate.num=='Plate_023'][1:96]
Plate_24$plate.num<-'Plate_024'
Plate_24$date<-180404
Plate_24[,samp_id:=paste0(date,'-',plate.num,'.',tube.name)]
Plate_24a<-Plate_24[,colnames(samplingdata),with=F]
Plate_24$date<-180405
Plate_24[,samp_id:=paste0(date,'-',plate.num,'.',tube.name)]
Plate_24$sample_flow_rate<-'2.0'
Plate_24b<-Plate_24[,colnames(samplingdata),with=F]
samplingdata<-rbind(samplingdata,Plate_24a,Plate_24b)

samplingdata$samp.plate<-sapply(samplingdata$samp_id,function(x){as.character(strsplit(x,"\\.")[[1]][1])})
samplingnotes<-fread(paste0(layoutdir,date,'-SamplePlateLayout.txt'),header=T,sep='\t')
# conditiondata<-fread(paste0(layoutdir,date,'-GrowthConditions.txt'),header=T,sep='\t')

rcs<-SP1[1:96,c('r','c','rc')]
rcDT<-data.table(ldply(lapply(1:nrow(samplingnotes),function(i){
	rcs$samp.plate<-samplingnotes$samp.plate[i]
	rcs
})))

samplingplatespre1<-merge(rcDT,samplingnotes,by='samp.plate');samplingplatespre1$samp.plate<-as.character(samplingplatespre1$samp.plate)
samplingplatespre1[,samp_id:=paste0(samp.plate,'.',rc)]
samplingplatespre<-merge(samplingplatespre1,samplingdata,by=c('samp_id','samp.plate'))
layout<-merge(samplingplatespre,sourceplate.layout,by=c('source.plate',colnames(rcs)))[order(samp.plate,r,c)]
table(layout$source.plate)/96
nrow(sourceplate.layout)
# still problems because plate_024 is missing from the samplingdata object....
write.table(layout, paste0(layoutdir,date,'-layout.txt'),row.names=F,sep='\t')
fread(paste0(layoutdir,date,'-layout.txt',sep=''))




