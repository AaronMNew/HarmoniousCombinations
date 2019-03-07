
library(plyr)
library(XML)
library(data.table)


head.dir<-"/Users/anew/Google Drive/002_lehnerlabData/180404-VarGALK.With.G80xG4_2"
date<-'180404' # beginning date clones were measured
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/'
outputdata<-'./output.data/'
layoutdir<-'./layout.sample.data/'
setwd(head.dir)




# # # code below parses a plate
# data <- xmlParse('Plate_014.xml')
# xml_data <- xmlToList(data)

# # data for A1
# x<-xml_data$experiment$tray$specimen[3]
# test<-NULL

# names(x$tube$loader_user_info)
# # below returns a table with useful information about the experiment
# test<-c("sample_flow_rate","sample_volume","record_sample_volume","mixing_volume","mixing_speed","number_mixes", "wash_volume",'tube.name')
# for(i in 3:87){
	# x<-xml_data$experiment$tray$specimen[i]
	# a<-unlist(c(as.vector(x$tube$loader_user_info),as.vector(x$tube$.attrs)))
	# test<-rbind(test,a)
# }


# extend the individual plate parsing above into entire experiments
template.directory<-'/Templates'
xml.dir.index<-list.files(paste(head.dir,template.directory,sep=''))
xmls<-xml.dir.index[grep('.xml',xml.dir.index)]
output1<-lapply(xmls,function(x){
	xab<-paste(head.dir,template.directory,'/',x,sep='')
	data2 <- xmlParse(xab) 
	xml_data2 <- xmlToList(data2)
	#as.data.frame(names(xml_data2$experiment)[9:(length(xml_data2$experiment)-3)]) # check you have identified the "trays" which are the plates
	test<-c("sample_flow_rate","sample_volume","record_sample_volume","mixing_volume","mixing_speed","number_mixes", "wash_volume",'tube.name','plate.num')
	plates<-xml_data2$experiment[9:(length(xml_data2$experiment)-3)]
	for (ix in 1:length(plates)){
		# define plate
		plate<-plates[ix]
		# define the range of tubes included in this plate
		len<-(length(plate$tray$specimen)-1)
		# get each well of the experiment into a specimens object you cycle through extracting information
		specimens<-plate$tray$specimen[3:len]
		print(length(specimens))
		for(i in 1:length(specimens)){
			# create tube-specific object
			x<-specimens[i]
			# the numbers below for date.plate are tailored for the naming of the different plates 
			# (in this case, a 6-digit date, a dash, then the standard Plate_XXX)
			date.plate<-substr(xab,(nchar(xab)-22),(nchar(xab)-23+6))
			# pull out the information you want
			a<-unlist(c(as.vector(x$tube$loader_user_info),as.vector(x$tube$.attrs), as.vector(plate$tray$.attrs[1]),date.plate=date.plate))
			# append it onto the test dummy variable
			test<-rbind(test,a)
		}
	}
	rownames(test)<-NULL
	output<-as.data.frame(test[2:nrow(test),])

	colnames(output)<-test[1,]
	nchar(xab)
	
	
	output$samp_id<-paste(date.plate, '-',output$plate.num,'.',output$tube.name,sep='')
	output
})
lapply(output1,function(d)nrow(d)/96)
lapply(output1,function(d)table(d$plate.num))

sampling.data <- data.table(ldply(output1, data.table))
sampling.data$sample_flow_rate.1<-NULL
setwd(head.dir)
write.table(sampling.data,paste0(layoutdir,date,'-SamplingData.txt'),row.names=FALSE,sep='\t')
fread(paste0(layoutdir,date,'-SamplingData.txt'),header=T,sep='\t')






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






# Scripts for dealing with flow cytometry data.
source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'


head.dir<-"/Users/anew/Google Drive/002_lehnerlabData/180404-VarGALK.With.G80xG4_2"
date<-'180404' # beginning date clones were measured
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/'
outputdata<-'./output.data/'
layoutdir<-'./layout.sample.data/'
setwd(head.dir)

br<-seq(-0.5,2.7,length=61)

# this file generates these objects which can be used in downstream analysis
# merged.final<-fread(paste(date,"-DatasetWithGenerationsAndGeneExpression.txt",sep=''),sep='\t',header=T)


###################################################
### custom functions for analysis
###################################################
source(source.code.path.file)
load.all() # this function loads all the standard packages you use in this experiment

###################################################
### analysis: first filtering on FSC and SSC
###################################################

setwd(plates.dir)

dirs<-list.dirs()[grep('Plate',list.dirs())]

# this for loop dumps a bunch of filtered flowCore data files in .rData format

for (ix in 1:length(dirs)){
	wash=NULL
	filtsize<-4
	setwd(plates.dir)
	setwd(dirs[ix])
	dir<-paste(dirs[ix])
	setwd(plates.dir);setwd(dir)
	fcs.files<-list.files(pattern='../*.fcs')
	flowframes<-read.flowSet(fcs.files, name.keyword="TUBE NAME", phenoData=list(name="TUBE NAME", Filename="$FIL"))
#	head(pData(phenoData(flowframes)))
	
	filtered <- fsApply(flowframes,function(fcs){
		if(!is.null(wash)){strt<-0}
		if(is.null(wash)){strt<-0}	
		len<-length(exprs(fcs$'FSC-A'))
		if(len<(strt+1000)){filteredPlus<-NULL}else{
		ssc<-quantile(exprs(fcs[strt:len,3]),c(0.1,0.9))
		fsc<-quantile(exprs(fcs[strt:len,1]),c(0.1,0.9))
		rectGate <- rectangleGate("FSC-H"=fsc,"SSC-H"=ssc)
		filteredRec<-split(fcs,flowCore::filter(fcs,rectGate))$`defaultRectangleGate+`
		c2f <- curv2Filter("FSC-H","SSC-H",bwFac=filtsize)
		filtered<-split(filteredRec,flowCore::filter(filteredRec,c2f))
		filteredPlus<-filtered$`area 1`}

		filteredPlus
	})
	setwd(plates.dir)
	save(filtered,file=paste0(head.dir,'/rdats.output/', substr(dir,3,18),'.rData',sep=''))
	
}

###################################################
### analysis: crunch numbers on filtered data 
### and merge with layout
###################################################

setwd(paste(head.dir,'/rdats.output',sep=''))
rdats<-list.files(pattern='*.rData')

br<-seq(-0.5,2.7,length=61)
denslabels<-as.character(as.data.frame(table(cut(rnorm(1000),breaks=br)))$Var1)
summary.names<-c('yfp.min',
				'yfp.5th', 
				'yfp.25th',
				'yfp.median',
				'yfp.75th',
				'yfp.95th',
				'yfp.max',
				'yfp.mean',
				'yfp.variance',
				'fracon',
				'mean.norm.fsc',
				'var.norm.fsc',
				'N',
				'mean.on.fraction',
				'var.on.fraction',
				'cells.per.sec',
				'dip.stat',
				'dip.pval',
				'mean.fsc',
				'median.fsc',
				'mean.ssc',
				'median.ssc',
				'FITC.saturation.corrected',
				'pct.oversaturation',
				denslabels)

# this for() loop analyzes what you'd like within the fcs data
output<-NULL
nullsamps<-NULL
nullsamps.summary<-NULL
for(i in 1:length(rdats)){
	load(paste(rdats[i]))
	print(rdats[i])
	
	# in above analysis, if the number of cells did not match a certain value, then samples were scored as NULL. 
	# These have to be removed for fsApply to work, hence all this bullshit below.
	nullsamps<-c();for(ix in 1:length(filtered)){nullsamps<-c(nullsamps,(is.null(filtered[[ix]])))}
	if(sum(as.numeric(nullsamps))>0){
			nulls<-paste(substr(rdats[i],1,nchar(rdats[i])-6),names(filtered[nullsamps]),sep='.')
			if(is.null(nullsamps.summary)){nullsamps.summary<-nulls}
				else{nullsamps.summary<-c(nullsamps.summary, nulls)}
			filtered<-flowSet(filtered[!nullsamps])

		}
	nullsamps<-NULL
	
		analyzed<-as.data.frame(fsApply(filtered,function(x){
			# debugging:
			# x<-filtered[[ix]]
			# instead of log, a biexponential is used, which is pseudo-log except around 0. It avoids lost numbers due to log transformation.
			# these parameters below for the biexponential were taken from somewhere on the internet and are not optimized at all.
			biexp<-biexponentialTransform(transformationId="biexptrans", 
				a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0, 
				tol = .Machine$double.eps^0.25, 
				maxit = as.integer(5000))
			biexp.dat<-transform(x,transformList(c("FITC-A"),biexp))
			df<-data.frame(exprs(biexp.dat$'FITC-A'), exprs(x$'FITC-A'))
			# 6.65 is the empirically decided cutoff for autofluorescence
			df.on<-data.frame(df[df[,1]>6.65,]) 
			
			FITCA<-exprs(x$'FITC-A')
			PEA<-exprs(x$'PE-A')
			FSC<-exprs(x$'FSC-A')
			# correct FITC-A saturation by regression with PE-A, if necessary
			# if the FITC-A signal is saturated, you can't use the variance measures to determine noise
			
			max<-262000 # near the saturation in the machine, which I believe to be 262144 AU (2^18 AU)
			if(length(FITCA[FITCA>max])>5){FITC.saturation.corrected<-1}else{FITC.saturation.corrected<-0}
			if(FITC.saturation.corrected==1){

				df.fl<-data.frame(PEA,FITCA,FSC)
				df.subsat1<-as.data.frame(df.fl[df.fl$FITC.A < max,])
				df.subsat2<-as.data.frame(df.fl[df.fl$FITC.A >1*10^3,])

				df.sat<-as.data.frame(df.fl[df.fl$FITC.A> max,])
				mod1<-lm(FITC.A~PE.A,data=df.subsat2)
				corrected<-predict.lm(mod1,data.frame(PE.A=df.sat$PE.A))
				
				df.sat$FITC.A<-corrected
				df<-data.frame(rbind(df.sat,df.subsat1))
				pct.oversaturation<-nrow(df.sat)/nrow(df)
			}else{df<-data.frame(PEA,FITCA,FSC);pct.oversaturation<-0}
		

			mean.norm.fsc<-mean(df$FITC.A / df$FSC.A)
			yfp<-df$FITC.A
			yfp.mean<-mean(yfp)
			var.yfp<-var(yfp,na.rm=TRUE) 
			fracon<-nrow(df.on)/nrow(df) 
			mean.on.fraction<-mean(df.on[,2])
			var.on.fraction<-var(df.on[,2])

			fitcalog<-log(df$FITC.A+10^3,base=10)-3
			# hist(fitcalog,main=paste(x@description$'TUBE NAME'),xlim=c(-0.5,2.7),breaks=30)}
			densities<-table(cut(fitcalog,breaks=br))



# debug / look at what you have
#			norm<-exprs(x$'FITC-A') / exprs(x$'FSC-A')
#			plot(yfp,norm,log='xy')
#			smoothScatter(yfp, exprs(x$'FSC-A'),log='xy')

			var.norm.fsc<-var(exprs(x$'FITC-A') / exprs(x$'FSC-A'))
			N<-length(exprs(x$'FITC-A'))
			# for cell density per ml, cells.per.sec * 1/??l per second sampling rate * 1000 ??l / ml
			cells.per.sec<-lm(1:length(exprs(x$Time))~exprs(x$Time))$coefficients[2]*1000 
			neg.absent<-log(FITCA[FITCA>0])
			diptest<-dip.test(neg.absent)
			mean.fsca<-mean(exprs(x$'FSC-A'))
			median.fsca<-median(exprs(x$'FSC-A'))
			mean.ssca<-mean(exprs(x$'SSC-A'))
			median.ssca<-median(exprs(x$'SSC-A'))

						
			a<-c(quantile(yfp,c(0,0.05,0.25,0.5,0.75,0.95,1)),
					yfp.mean,
					var.yfp,
					fracon,
					mean.norm.fsc,
					var.norm.fsc,
					N,
					mean.on.fraction,
					var.on.fraction,
					cells.per.sec,
					diptest$statistic,
					diptest$p.value,
					mean.fsca,
					median.fsca,
					mean.ssca,
					median.ssca,
					FITC.saturation.corrected,
					pct.oversaturation,
					densities) # always keep densities at the end!)
			a
		}))

	colnames(analyzed)<-summary.names
	plate.name<-strsplit(rdats[i],'[.]')[[1]][1]
	analyzed$samp_id<-paste(plate.name,'.',rownames(analyzed),sep='')

	if(is.null(output)){output<-analyzed}else{output<-rbind(output,analyzed)}
	

}

unique(warnings())
# typical warning that i just ignore:
# > Warning message:
# > In FUN(if (use.exprs) exprs(y) else y, ...) : NaNs produced
setwd(head.dir)
# first time running script, write output.
write.table(output,paste(outputdata,date,'-analyzed.tab',sep=''),sep="\t",quote=FALSE, row.names=FALSE)
write.table(nullsamps.summary, paste(outputdata,date,'-null.samps.tab',sep=''),sep="\t",quote=FALSE, row.names=FALSE)




###################################################
### load data
###################################################

setwd(head.dir)
datloc<-paste(outputdata,date,'-analyzed.tab',sep='')
nullsamploc<-paste(outputdata,date,'-null.samps.tab',sep='')
# analyzed data have already been written from the Preprocessing.R file and here are fread output written the first time around.
output<-fread(paste(outputdata,date,'-analyzed.tab',sep=''),header=TRUE,sep='\t')[]
old.cols<-colnames(output)[c(grep("\\]",colnames(output)))]
dens.cols<-paste('dens',round(br[1:(length(br)-1)],2),sep='.')
setnames(output,old= old.cols,new= dens.cols)


nullsamps.summary<-fread(nullsamploc,header=TRUE,sep='\t')
# layout is a table with headers and at least a "samp_id" column either generated in the layout.maker script or by hand in excel.
layout<-fread(paste(layoutdir,date,'-layout.txt',sep=''),header=TRUE,sep='\t')

if(!is.null(nullsamps.summary)){
	#create a dummy data frame for null samples generated and add to output
	nullsamps.layout <-as.data.frame(t(apply(nullsamps.summary,1,function(x){
		c(rep(NA,(ncol(output)-1)),x)
	})))
	colnames(nullsamps.layout)<-colnames(output)
	output.dummy<-rbind(output,nullsamps.layout,make.row.names=FALSE)

	output.final<-data.table(t(apply(output.dummy[,1:(ncol(output.dummy)-1)],1,as.numeric)),output.dummy$samp_id)
	colnames(output.final)<-colnames(output)
	}else{output.final<-output}


merged<-merge(output.final,layout,by='samp_id')[!is.na(GAL3.clone.named)][order(samp.plate,r,c)]
merged[,measurement.date:=sapply(samp.plate,function(x)as.numeric(strsplit(x,'-')[[1]][1]))]
vars<-c('cond','GAL3.clone.named','GAL80.clone.named','GAL4.clone.named','GALK.clone.named')
# merged[,biol.rep:=1:length(yfp.mean),by=vars]
# data.frame(merged[,c('clone.id','samp_id.prefix','r.x','c.x','biol.rep')])


###################################################
### add a few things: compute number of generations + phenotypic index
###################################################

dens.norm<-data.table(t(apply(merged[, dens.cols,with=FALSE],1,row.normalization)))
colnames(dens.norm)<-dens.cols
merged[,dens.cols]<-dens.norm
merged2<-merged
# x<-c()
# x[1]<-merged2$clone.id[1]
# x[2]<-merged2$genotype[1]

# for cell density per ml, cells.per.sec * 1/µl per second sampling rate * 1000 µl / ml

merged2$cell.per.ml<-merged2$cells.per.sec*(1/as.numeric(merged2$sample_flow_rate))*1000*(1/merged2$fold.dilution)
DT.glu<-merged2[glu==0.1]
DT.gal<-merged2[gal==0.2]

# to couple direct measurements with each other (glucose > galactose) i need to create an id variable that will allow joining glu with gal measurements
DT.glu$next.day<-as.numeric(DT.glu$measurement.date)+1
DT.gal$next.day<-NA
DT.glu$id<-paste(DT.glu$source.plate,DT.glu$clone.id,DT.glu$next.day)
DT.gal$id<-paste(DT.gal$source.plate,DT.gal$clone.id,DT.gal$measurement.date)

DT.glu$id%in%DT.gal$id
table(DT.glu$id%in%DT.gal$id)
# quick check of how well they match
a<-DT.glu[match(DT.glu$id,DT.gal$id)]$id
b<-DT.gal$id
table(a==b)

DT.gal[!id %in% DT.glu$id]
DT.glu[!id %in% DT.gal$id]

# i re-merge them here and then calculate generations, growth rate and phenotypic index
# setnames(DT.glu, colnames(output),paste0(colnames(output),'.glu'))
# setnames(DT.gal, colnames(output),paste0(colnames(output),'.gal'))
setnames(DT.glu, colnames(DT.glu),paste0(colnames(DT.glu),'.glu'))
setnames(DT.gal, colnames(DT.gal),paste0(colnames(DT.gal),'.gal'))

setnames(DT.glu,'id.glu','id')
setnames(DT.gal,'id.gal','id')
DT.merge1<-merge(DT.glu,DT.gal,by='id')
duplicated.columns<-gsub('.glu','',gsub(".gal",'',colnames(DT.merge1[,duplicated(t(DT.merge1)),with=F])))
duplicated.layouts.glu<-unlist(lapply(duplicated.columns[duplicated.columns%in%colnames(layout)],function(x){paste0(x,'.glu',collapse='')}))
duplicated.layouts.gal<-unlist(lapply(duplicated.columns[duplicated.columns%in%colnames(layout)],function(x){paste0(x,'.gal',collapse='')}))
DT.merge<-merge(DT.glu[,colnames(DT.glu)%w/o%duplicated.layouts.glu,with=F],DT.gal,by='id')
gsub('.gal','',duplicated.layouts.gal)[!duplicated(gsub('.gal','',duplicated.layouts.gal))]
setnames(DT.merge,duplicated.layouts.gal[!duplicated(duplicated.layouts.gal)],gsub('.gal','',duplicated.layouts.gal[!duplicated(duplicated.layouts.gal)]))

dil.factor<-1;hours<-12
DT.merge$generations<-log((DT.merge$cell.per.ml.gal/((DT.merge$cell.per.ml.glu)*dil.factor)),base=2)
DT.merge$growth.rate<-log((DT.merge$cell.per.ml.gal/((DT.merge$cell.per.ml.glu)*dil.factor)))/hours
DT.merge$phenotypic.index<-DT.merge$fracon.glu+DT.merge$fracon.gal

write.table(DT.merge,paste0(outputdata,date,'-preprocessed_data.txt'),sep='\t',row.names=F)
fread(paste0(outputdata,date,'-preprocessed_data.txt'),sep='\t',header=T)

###################################################
### quick checks
###################################################

vars1<-c('GALK.clone.named',"GAL80.clone.named","GAL4.clone.named")
test<-merge(DT.merge[GAL3.clone.named=='GAL3.delta',summary.func(yfp.mean.gal),by=vars1],DT.merge[GAL3.clone.named=='GAL3.WT',summary.func(yfp.mean.gal),by=vars1],by=vars1)
test2<-merge(DT.merge[GAL3.clone.named=='GAL3.delta',summary.func(yfp.mean.glu),by=vars1],DT.merge[GAL3.clone.named=='GAL3.WT',summary.func(yfp.mean.glu),by=vars1],by=vars1)

plot(test$mean.y,test$mean.x);abline(a=0,b=1)

test$fold.induction.delta<-test$mean.x/test2$mean.x
test$fold.induction.WT<-test$mean.y/test2$mean.y
dev.new();plot(test$fold.induction.WT,test$fold.induction.delta,log='xy');abline(a=0,b=1)

franplot(test$fold.induction.delta,names.arg=paste(test$GAL80.clone.named,test$GAL4.clone.named))
dev.new()
franplot(test$fold.induction.WT,names.arg=paste(test$GAL80.clone.named,test$GAL4.clone.named))

######
test<-merge(DT.merge[GAL3.clone.named=='GAL3.delta',summary.func(fracon.gal,'GAL3.delta.gal'),by=vars1],DT.merge[GAL3.clone.named=='GAL3.WT',summary.func(fracon.gal,'GAL3.WT.gal'),by=vars1],by=vars1)
test2<-merge(DT.merge[GAL3.clone.named=='GAL3.delta',summary.func(fracon.glu,'GAL3.delta.glu'),by=vars1],DT.merge[GAL3.clone.named=='GAL3.WT',summary.func(fracon.glu,'GAL3.WT.glu'),by=vars1],by=vars1)

qplot(test$GAL3.WT.gal_mean,test$GAL3.delta.gal_mean,col=factor(test$GALK.clone.named))+geom_abline(intercept=0,slope=1)

test$fold.induction.delta<-test$GAL3.delta.gal_mean/test2$GAL3.delta.glu_mean
test$fold.induction.WT<-test$GAL3.WT.gal_mean/test2$GAL3.WT.glu_mean
dev.new();qplot(test$fold.induction.WT,test$fold.induction.delta,log='xy',col=test$GALK.clone.named)+geom_abline(intercept=0,slope=1)

test[fold.induction.WT>20]

######
test<-merge(DT.merge[GAL3.clone.named=='GAL3.delta',summary.func(yfp.mean.gal,'GAL3.delta.gal'),by=vars1],DT.merge[GAL3.clone.named=='GAL3.WT',summary.func(yfp.mean.gal,'GAL3.WT.gal'),by=vars1],by=vars1)
test2<-merge(DT.merge[GAL3.clone.named=='GAL3.delta',summary.func(yfp.mean.glu,'GAL3.delta.glu'),by=vars1],DT.merge[GAL3.clone.named=='GAL3.WT',summary.func(yfp.mean.glu,'GAL3.WT.glu'),by=vars1],by=vars1)

qplot(test$GAL3.WT.gal_mean,test$GAL3.delta.gal_mean,col=factor(test$GALK.clone.named))+geom_abline(intercept=0,slope=1)

test$fold.induction.delta<-test$GAL3.delta.gal_mean/test2$GAL3.delta.glu_mean
test$fold.induction.WT<-test$GAL3.WT.gal_mean/test2$GAL3.WT.glu_mean
dev.new();qplot(test$fold.induction.WT,test$fold.induction.delta,log='xy',col=test$GALK.clone.named)+geom_abline(intercept=0,slope=1)

test[fold.induction.WT>20]

