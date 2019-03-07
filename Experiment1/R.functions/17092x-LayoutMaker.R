# Make a layout file for experiment
source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'

head.dir<-"/Users/anew/Desktop/17092x-CopyNumberExperiment"
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/'
outputdata<-'./output.data/'
layoutdir<-'./layout.sample.data/'
date<-'170927' # beginning date clones were measured
setwd(head.dir)

br<-seq(-0.5,2.7,length=61)



###################################################
### custom functions for analysis
###################################################
source(source.code.path.file)
load.all() # this function loads all the standard packages you use in this experiment


###################################################
### custom functions for analysis
###################################################

a<-fread(paste0(layoutdir,'17092x-SamplingData.txt'))
a[,rc:=unlist(sapply(samp_id,function(x)strsplit(x,'\\.')[[1]][2]))]
a[,c('r','c','rc','measurement.date'):= data.table(t(a[,unlist(sapply(samp_id,function(x){
	rc<-strsplit(x,'\\.')[[1]][2]
	md<-strsplit(x,'-')[[1]][1]
	len<-nchar(rc)
	if(len==2){
		r<-substr(rc,1,1)
		c<-substr(rc,2,2)
	}
	if(len==3){
		r<-substr(rc,1,1)
		c<-substr(rc,2,3)
	}
	c(r,c,rc,md)
	}))]))
]

# run below to see the problems with sample id's in this set
# plateplot(a,'sample_flow_rate')+facet_grid(~plate.num)
# OK it's clear that the xml files had some problems. Plates 009 - 016 all have the correct data, and from my notes I know that the sampling rate was 0.5 µl per second.

b<-fread(paste0(layoutdir,'1707xx-17092x-LayoutByHand.txt'))[!is.na(genotype.strain)]
asub<-a[plate.num%in%c('Plate_009','Plate_010','Plate_011','Plate_012','Plate_013','Plate_014','Plate_015','Plate_016'),.(sample_flow_rate,samp_id)]
asub[,samp_id:=gsub('170927','170928',samp_id)]
DTm<-merge(asub,b,by='samp_id',all=T)
DTm[is.na(sample_flow_rate),sample_flow_rate:=0.5]


layout<-DTm[,.(
	samp_id,
	sample_flow_rate,
	gal,
	glu,
	plate,
	source.plate,
	r,
	c,
	rc,
	genotype.strain,
	genotype.plasmid,
	genotype.strain1,
	strain.arb,
	delta.GAL3,
	delta.GAL80,
	delta.GAL4,
	copies.GAL3.genome,
	copies.GAL80.genome,
	copies.GAL4.genome,
	tx.mix.arb,
	pRS415,
	copies.GAL3.plasmid,
	copies.GAL80.plasmid,
	copies.GAL4.plasmid,
	total.copies.GAL3,
	total.copies.GAL80,
	total.copies.GAL4
	)]


layout[,plasgeno:=applyPaste(
	data.frame(	copies.GAL3.plasmid,
	copies.GAL80.plasmid,
	copies.GAL4.plasmid,
	copies.GAL3.genome,
	copies.GAL80.genome,
	copies.GAL4.genome),' '
)]

layout[,plasmidcopies:=applyPaste(
	data.frame(	copies.GAL3.plasmid,
	copies.GAL80.plasmid,
	copies.GAL4.plasmid),' '
)]


layout[,genomecopies:=applyPaste(
	data.frame(
	copies.GAL3.genome,
	copies.GAL80.genome,
	copies.GAL4.genome),' '
)]

layout[,totalcopies:=applyPaste(data.frame(	total.copies.GAL3,
	total.copies.GAL80,
	total.copies.GAL4),' '
	)]

layout[,binarycount:=applyPaste(data.frame(	as.numeric(total.copies.GAL3>0),
	as.numeric(total.copies.GAL80>0),
	as.numeric(total.copies.GAL4>0)),' '
	)]

layout[,measurement.date:=unlist(sapply(samp_id,function(x)as.numeric(strsplit(x,'-')[[1]][1])))]

layout[measurement.date=='170927',c('gal','glu'):=list(0,1)]
layout[measurement.date=='170928',c('gal','glu'):=list(1,0)]


write.table(layout,paste0(layoutdir,'17092x-layout.txt'),row.names=F,sep='\t')
fread(paste0(layoutdir,'17092x-layout.txt'))
	
	
	
