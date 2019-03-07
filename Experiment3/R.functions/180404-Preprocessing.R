

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

nullsamps.summary<-NULL
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

### BE VERY CAREFUL HERE TO BE SURE ID'S ACTUALLY MAKE SENSE WITH EXPERIMENT
DT.glu$id<-paste(DT.glu$plate.num,DT.glu$clone.id,DT.glu$next.day)
DT.gal$id<-paste(DT.gal$plate.num,DT.gal$clone.id,DT.gal$measurement.date)

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



######################################################################################################
### MORE PREPROCESSING THAN NORMAL HERE
######################################################################################################

###################################################
### read data for samples 
###################################################
DT.merge<-fread(paste0(outputdata,date,'-preprocessed_data.txt'),sep='\t',header=T)
DT00<-DT.merge[order(measurement.date.glu,samp.plate.glu,r)]
G3alleles<-na.exclude(unique(DT00 $GAL3.clone.named))
G80alleles<-na.exclude(unique(DT00 $GAL80.clone.named))
G4alleles<-na.exclude(unique(DT00 $GAL4.clone.named))
G3bin<-data.table(t(ldply(lapply(G3alleles,function(x){
	grepl(x, DT00 $GAL3.clone.named)
}))))
table(apply(G3bin,1,sum))
G80bin<-data.table(t(ldply(lapply(G80alleles,function(x){
	grepl(x, DT00 $GAL80.clone.named)
}))))
table(apply(G80bin,1,sum))

G4bin<-data.table(t(ldply(lapply(G4alleles,function(x){
	grepl(x, DT00 $GAL4.clone.named)
}))))
table(apply(G4bin,1,sum))

binaryalleles<-data.table(G3bin,G80bin,G4bin)
colnames(binaryalleles)<-c(G3alleles,G80alleles,G4alleles)
DT0<-data.table(DT00,binaryalleles)

###################################################
### set variables
###################################################
phenocols<-c('yfp.mean','mean.on.fraction','var.on.fraction','fracon','dip.stat')
galcols<-c(c('growth.rate','phenotypic.index'),paste(phenocols,'.gal',sep=''))
glucols<-paste(phenocols,'.glu',sep='')
coordcols<-colnames(DT0)[ grep('Dim',colnames(DT0))]
dependents<-c(galcols,glucols,coordcols)
genocols<-c('genotype','clone.named','GAL4.allele', 'GAL80.allele', 'GAL3.allele')
clone.named.set<-colnames(DT0)[grep('clone.name',colnames(DT0))]
groupings<-c('GAL3.clone.named','GAL80.clone.named','GAL4.clone.named','clone.named')
GALgenes<-c('GAL3','GAL80','GAL4')
groupings2<-c('GAL3.clone.named','GAL80.clone.named','GAL4.clone.named')


batchcols1<-c('replicate.glu', # 
	'number.mutants.glu', # this scores the number of mutations present compared to WT (0,1, or 2)
	'id',# id joins the glucose sample with the same one the next day that was measured in galactose (date term is the date gal measurement was made)
	'samp_id.glu',
	'measurement.date.glu',
	'measurement.date.gal',
	'media.batch.date.glu',
	'source.plate.glu',
	'plate.name.in.source.fcs.glu',
	'plate.glu',
	'plate.gal',
	'biol.rep.glu',
	'onemut.glu',
	'twomut.glu',
	'samp_id.gal',
	'biol.rep.gal',
	'replicate.gal')
	
batchcols2<-c('measurement.date.glu',
	'GAL3.clone.named',
	'GAL80.clone.named',
	'GAL4.clone.named',
	'GAL3.allele',
	'GAL80.allele',
	'GAL4.allele',
	'genotype',
	'measurement.date.gal',
	'source.plate',
	'clone.named',
	'tx.date'
	)

###################################################
### get orders right for clone.named factors
###################################################

G3ctrl<-c('GAL3.delta','GAL3.WT')
G80ctrl<-c('GAL80.delta','GAL80.WT')
G4ctrl<-c('GAL4.delta','GAL4.WT')

GAL3s<-factor(unique(c(DT0$GAL3.clone.named,G3ctrl)),levels=c(G3ctrl,unique(DT0$GAL3.clone.named%w/o%G3ctrl)))
GAL80s<-factor(unique(c(DT0$GAL80.clone.named,G80ctrl)),levels=c(G80ctrl,unique(DT0$GAL80.clone.named%w/o%G80ctrl)))
GAL4s<-factor(unique(c(DT0$GAL4.clone.named,G4ctrl)),levels=c(G4ctrl,unique(DT0$GAL4.clone.named%w/o%G4ctrl)))

DT0$GAL3.clone.named <-factor(DT0$GAL3.clone.named,levels=c(G3ctrl,unique(DT0$GAL3.clone.named%w/o%G3ctrl)))
DT0$GAL80.clone.named <-factor(DT0$GAL80.clone.named,levels=c(G80ctrl,unique(DT0$GAL80.clone.named%w/o%G80ctrl)))
DT0$GAL4.clone.named <-factor(DT0$GAL4.clone.named,levels=c(G4ctrl,unique(DT0$GAL4.clone.named%w/o%G4ctrl)))

###################################################
### MUST RUN to correct mean signal in glucose problems with first 4 plates of the day ()
###################################################
# save the old data in _old columns
DT0[,yfp.mean.glu_old:=yfp.mean.glu]
DT0[,mean.on.fraction.glu_old:=mean.on.fraction.glu]
DT0[,biol.rep:=duplicated(clone.id)]
table(DT0$biol.rep)
reps1<-DT0[biol.rep==FALSE]
reps2<-DT0[biol.rep==TRUE]
reps<-merge(reps1,reps2,by='clone.id')
byfac1<-c(clone.named.set, 'samp.plate.glu')
byfac<-c(paste0(byfac1,'.x'))
rep1_mean<-reps[,lapply(.SD,mean),by= byfac,.SDcols=paste0(dependents,'.x')]
rep2_mean<-reps[,lapply(.SD,mean),by= byfac,.SDcols= paste0(dependents,'.y')]
rep1_upper<-reps[,lapply(.SD,function(x)summary.func(x)$upper),by= byfac,.SDcols= paste0(dependents,'.x')]
rep2_upper<-reps[,lapply(.SD,function(x)summary.func(x)$upper),by= byfac,.SDcols= paste0(dependents,'.y')]
rep1_lower<-reps[,lapply(.SD,function(x)summary.func(x)$lower),by= byfac,.SDcols= paste0(dependents,'.x')]
rep2_lower<-reps[,lapply(.SD,function(x)summary.func(x)$lower),by= byfac,.SDcols= paste0(dependents,'.y')]
colnames(rep1_mean)<-colnames(rep2_mean)<-colnames(rep1_upper)<-colnames(rep2_upper)<-colnames(rep1_lower)<-colnames(rep2_lower)<-c(byfac1,dependents)


problemplates<-c('180404-Plate_001','180404-Plate_002','180404-Plate_003','180404-Plate_004')

# this is an overview of what you're doing
x<-problemplates[1]
# growth rate looks fine
plot(rep1_mean[samp.plate.glu==x]$growth.rate,rep2_mean[samp.plate.glu==x]$growth.rate);abline(a=0,b=1)
# systematic bias for rep1 vs rep2 in log space
plot(rep1_mean[samp.plate.glu==x]$yfp.mean.glu,rep2_mean[samp.plate.glu==x]$yfp.mean.glu,log='xy');abline(a=0,b=1)
# define 'bad' samples to be corrected
bad<-log(rep1_mean[samp.plate.glu==x]$yfp.mean.glu)
# define 'good' samples to be corrected
good<-log(rep2_mean[samp.plate.glu==x]$yfp.mean.glu)
# re-plot to just make sure you're not crazy
plot(bad,good);abline(a=0,b=1)
# determine the mean effect in log space
meaneffect<-mean(good-bad)
# calculate the corrected values of the bad samples by adding the mean effect
corrected<-bad+meaneffect
# plot to be sure you're not crazy
plot(corrected,good);abline(a=0,b=1)
# move measurements out of log space and plot to be sure you're not crazy
plot(exp(corrected), rep2_mean[samp.plate.glu==x]$yfp.mean.glu);abline(a=0,b=1)
# correct the values in the yfp.mean col
DT0[samp.plate.glu==x, exp(log(yfp.mean.glu)+meaneffect)]

# Do above for the yfp.mean.glu
for(i in 1:length(problemplates)){
	x<-problemplates[i]
	bad<-log(rep1_mean[samp.plate.glu==x]$yfp.mean.glu)
	good<-log(rep2_mean[samp.plate.glu==x]$yfp.mean.glu)
	# plot(good,bad);abline(a=0,b=1)
	meaneffect<-mean(good-bad,na.rm=T)
	corrected<-bad+meaneffect
	test<-t.test(good,bad,paired=T)
	print(test);print(test$p.value)
	# same as above but correct the data directly in the data.table
	if(test$p.value<0.001)DT0[samp.plate.glu==x, yfp.mean.glu:=exp(log(yfp.mean.glu)+meaneffect)]
	}
# Do above for the yfp.mean.gal
for(i in 1:length(problemplates)){
	x<-problemplates[i]
	bad<-log(rep1_mean[samp.plate.glu==x]$yfp.mean.gal)
	good<-log(rep2_mean[samp.plate.glu==x]$yfp.mean.gal)
	meaneffect<-mean(good-bad,na.rm=T)
	corrected<-bad+meaneffect
	test<-t.test(good,bad,paired=T)
	print(test);print(test$p.value)
	# same as above but correct the data directly in the data.table
	if(test$p.value<0.001)DT0[samp.plate.glu==x, yfp.mean.gal:=exp(log(yfp.mean.gal)+meaneffect)]
	}


# repeat plotting as above to be sure you've removed the plate-specific effects

DT0[,biol.rep:=duplicated(clone.id)]
table(DT0$biol.rep)
reps1<-DT0[biol.rep==FALSE]
reps2<-DT0[biol.rep==TRUE]
reps<-merge(reps1,reps2,by='clone.id')
byfac1<-c(clone.named.set, 'samp.plate.glu')
byfac<-c(paste0(byfac1,'.x'))
rep1_mean<-reps[,lapply(.SD,mean),by= byfac,.SDcols=paste0(dependents,'.x')]
rep2_mean<-reps[,lapply(.SD,mean),by= byfac,.SDcols= paste0(dependents,'.y')]
rep1_upper<-reps[,lapply(.SD,function(x)summary.func(x)$upper),by= byfac,.SDcols= paste0(dependents,'.x')]
rep2_upper<-reps[,lapply(.SD,function(x)summary.func(x)$upper),by= byfac,.SDcols= paste0(dependents,'.y')]
rep1_lower<-reps[,lapply(.SD,function(x)summary.func(x)$lower),by= byfac,.SDcols= paste0(dependents,'.x')]
rep2_lower<-reps[,lapply(.SD,function(x)summary.func(x)$lower),by= byfac,.SDcols= paste0(dependents,'.y')]
colnames(rep1_mean)<-colnames(rep2_mean)<-colnames(rep1_upper)<-colnames(rep2_upper)<-colnames(rep1_lower)<-colnames(rep2_lower)<-c(byfac1,dependents)


i<-1
plot.list<-list()
plot.list<-lapply(dependents,function(x){
	dat<-data.table(a=rep1_mean[,x,with=F],b=rep2_mean[,x,with=F])
	setnames(dat,colnames(dat),paste0(x,c('_rep1','_rep2')))
	uppers<-data.table(a=rep1_upper[,x,with=F],b=rep2_upper[,x,with=F])
	lowers<-data.table(a=rep1_lower[,x,with=F],b=rep2_lower[,x,with=F])
	dat<-data.frame(dat)
	ggplot(dat,aes(x=dat[,1],y=dat[,2]))+
	geom_point(alpha=0.2,aes(colour=factor(rep1_mean$samp.plate.glu=='180404-Plate_001'| rep1_mean$samp.plate.glu=='180404-Plate_002'| rep1_mean$samp.plate.glu=='180404-Plate_003'| rep1_mean$samp.plate.glu=='180404-Plate_004')))+
	geom_errorbarh(xmax=uppers$a,xmin=lowers$a,alpha=0.1)+geom_errorbar(ymax=uppers$b,ymin=lowers$b,alpha=0.1)+
		ylab(colnames(dat)[2])+xlab(colnames(dat)[1])+geom_abline(intercept=0,slope=1,linetype='dotted',colour='indianred')+theme(legend.position='none')
	
})
do.call(grid.arrange,plot.list)

###################################################
### check replicate plates
###################################################
tx1_mean<-DT0[tx.rep==1,lapply(.SD,mean),by=c(clone.named.set),.SDcols=dependents]
tx2_mean<-DT0[tx.rep==2,lapply(.SD,mean),by=c(clone.named.set),.SDcols=dependents]
tx1_upper<-DT0[tx.rep==1,lapply(.SD,function(x)summary.func(x)$upper),by=c(clone.named.set),.SDcols=dependents]
tx2_upper<-DT0[tx.rep==2,lapply(.SD,function(x)summary.func(x)$upper),by=c(clone.named.set),.SDcols=dependents]
tx1_lower<-DT0[tx.rep==1,lapply(.SD,function(x)summary.func(x)$lower),by=c(clone.named.set),.SDcols=dependents]
tx2_lower<-DT0[tx.rep==2,lapply(.SD,function(x)summary.func(x)$lower),by=c(clone.named.set),.SDcols=dependents]

i<-1
plot.list<-list()
plot.list<-lapply(dependents,function(x){
	dat<-data.table(a=tx1_mean[,x,with=F],b=tx2_mean[,x,with=F])
	setnames(dat,colnames(dat),paste0(x,c('_tx1','_tx2')))
	uppers<-data.table(a=tx1_upper[,x,with=F],b=tx2_upper[,x,with=F])
	lowers<-data.table(a=tx1_lower[,x,with=F],b=tx2_lower[,x,with=F])
	dat<-data.frame(dat)
	ggplot(dat,aes(x=dat[,1],y=dat[,2]))+geom_point(alpha=0.2,colour='cornflowerblue')+geom_errorbarh(xmax=uppers$a,xmin=lowers$a,alpha=0.1)+geom_errorbar(ymax=uppers$b,ymin=lowers$b,alpha=0.1)+
		ylab(colnames(dat)[2])+xlab(colnames(dat)[1])+geom_abline(intercept=0,slope=1,linetype='dotted',colour='indianred')
	
})
w<-12;h<-w*1
ggsave(paste0(figout,'180404-ReproducibilityAcrossTXreps.pdf'),do.call(grid.arrange,plot.list),width=w,height=h)

# do.call(grid.arrange,plot.list)

DT0[,biol.rep:=duplicated(clone.id)]
table(DT0$biol.rep)
reps1<-DT0[biol.rep==FALSE]
reps2<-DT0[biol.rep==TRUE]
reps<-merge(reps1,reps2,by='clone.id')
byfac1<-c(clone.named.set, 'samp.plate.glu')
byfac<-c(paste0(byfac1,'.x'))
rep1_mean<-reps[,lapply(.SD,mean),by= byfac,.SDcols=paste0(dependents,'.x')]
rep2_mean<-reps[,lapply(.SD,mean),by= byfac,.SDcols= paste0(dependents,'.y')]
rep1_upper<-reps[,lapply(.SD,function(x)summary.func(x)$upper),by= byfac,.SDcols= paste0(dependents,'.x')]
rep2_upper<-reps[,lapply(.SD,function(x)summary.func(x)$upper),by= byfac,.SDcols= paste0(dependents,'.y')]
rep1_lower<-reps[,lapply(.SD,function(x)summary.func(x)$lower),by= byfac,.SDcols= paste0(dependents,'.x')]
rep2_lower<-reps[,lapply(.SD,function(x)summary.func(x)$lower),by= byfac,.SDcols= paste0(dependents,'.y')]
colnames(rep1_mean)<-colnames(rep2_mean)<-colnames(rep1_upper)<-colnames(rep2_upper)<-colnames(rep1_lower)<-colnames(rep2_lower)<-c(byfac1,dependents)


i<-1
plot.list<-list()
plot.list<-lapply(dependents,function(x){
	dat<-data.table(a=rep1_mean[,x,with=F],b=rep2_mean[,x,with=F])
	setnames(dat,colnames(dat),paste0(x,c('_rep1','_rep2')))
	uppers<-data.table(a=rep1_upper[,x,with=F],b=rep2_upper[,x,with=F])
	lowers<-data.table(a=rep1_lower[,x,with=F],b=rep2_lower[,x,with=F])
	dat<-data.frame(dat)
	# ggplot(dat,aes(x=dat[,1],y=dat[,2]))+geom_point(alpha=0.2,colour='cornflowerblue')+geom_errorbarh(xmax=uppers$a,xmin=lowers$a,alpha=0.1)+geom_errorbar(ymax=uppers$b,ymin=lowers$b,alpha=0.1)+
		# ylab(colnames(dat)[2])+xlab(colnames(dat)[1])+geom_abline(intercept=0,slope=1,linetype='dotted',colour='indianred')
	ggplot(dat,aes(x=dat[,1],y=dat[,2]))+geom_point(alpha=0.2,aes(colour=factor(rep1_mean$samp.plate.glu=='180404-Plate_001'| rep1_mean$samp.plate.glu=='180404-Plate_002'| rep1_mean$samp.plate.glu=='180404-Plate_003'| rep1_mean$samp.plate.glu=='180404-Plate_004')))+geom_errorbarh(xmax=uppers$a,xmin=lowers$a,alpha=0.1)+geom_errorbar(ymax=uppers$b,ymin=lowers$b,alpha=0.1)+
		ylab(colnames(dat)[2])+xlab(colnames(dat)[1])+geom_abline(intercept=0,slope=1,linetype='dotted',colour='indianred')+theme(legend.position='none')
	
})
w<-12;h<-w*1
plot.list$top='blue colors are the first plates measured at the beginning of the day, which had a batch effect due to immature fluorophores'
ggsave(paste0(figout,'180404-ReproducibilityAcrossBiolReps.pdf'),do.call(grid.arrange,plot.list),width=w,height=h)
# do.call(grid.arrange,plot.list)


###################################################
### calculate total protein production
###################################################
background.gal<-median(DT0[GAL4.clone.named=='GAL4.delta']$yfp.mean.gal)
background.glu<-median(DT0[GAL4.clone.named=='GAL4.delta']$yfp.mean.glu)
DT0[,yfp.gal:=sapply(yfp.mean.gal-background.gal,function(x)max(x,1))]
DT0[,yfp.glu:=sapply(yfp.mean.glu-background.glu,function(x)max(x,1))]
DT0[,log.protein.production:=log2((yfp.gal)/(yfp.glu*1/(2^generations)))]
DT0[,fac:=paste0(GALK.clone.named,'+',GAL3.clone.named)]
ggplot(DT0,aes(x=phenotypic.index,y= log.protein.production))+geom_point(aes(colour=fac))

DTm<-merge(DT0[,summary.func(phenotypic.index),by=c(clone.named.set,colnames(binaryalleles))],DT0[,summary.func(log.protein.production),by=c(clone.named.set,colnames(binaryalleles))],by=c(clone.named.set,colnames(binaryalleles)))
DTm[,GAL4.delta:=GAL4.clone.named=='GAL4.delta']
DTm[,fac:=paste0(GALK.clone.named,'+',GAL3.clone.named)]
# DTm$fac<-DTm $GAL4.WT==T & DTm $GAL3.WT==T & DTm $GAL80.WT==T


xlabel<-'phenotypic.index'
ylabel<-'log 2 protein produced'
ggplot(DTm,aes(mean.x,mean.y))+geom_point(aes(col=fac,shape=GALK.clone.named))+
	geom_errorbarh(xmax=DTm$upper.x,xmin=DTm$lower.x,aes(col=fac),alpha=0.2)+
	geom_errorbar(ymax=DTm$upper.y,ymin=DTm$lower.y,aes(col=fac),alpha=0.2)+
	xlab(xlabel)+ylab(ylabel)
	

###################################################
### quick look at GAL4-L868C vs GAL4.WT: they behave almost exactly the same except in fracon.gal.
###################################################

DTm[mean.y>11]
qplot(DTm[GAL4.clone.named=='GAL4-L868C']$mean.y,DTm[GAL4.clone.named=='GAL4.WT']$mean.y)
qplot(DTm[GAL4.clone.named=='GAL4-L868C']$mean.x,DTm[GAL4.clone.named=='GAL4.WT']$mean.x)
L868C<-DT0[GAL4.clone.named=='GAL4-L868C',dependents,with=F]
WT<-DT0[GAL4.clone.named=='GAL4.WT',dependents,with=F]
plot.list<-lapply(1:ncol(WT),function(i){
	DT<-data.table(WT[,i,with=F], L868C[,i,with=F])
	colnames(DT)<-paste0(c('GAL4.WT-','GAL4-L8686C-'),colnames(DT))
	df<-as.data.frame(DT)
	ggplot(df,aes(x= df[,colnames(df)[1]],y= df[,colnames(df)[2]]))+geom_point()+ylab(colnames(df)[2])+xlab(colnames(df)[1])+
		geom_abline(intercept=0,slope=1,col='red')
})
do.call(grid.arrange,plot.list)

###################################################
### quick look at GALK.Can_alb: clone pAMN51.1 - 1 - A7 is problematic
###################################################
unique(DT0[GALK.clone.named=='GALK.Can_abl']$GALK.allele)
DT0[GALK.clone.named=='GALK.Can_abl',summary.func(growth.rate),by='GALK.allele']

DT1<-DT0[GALK.clone.named=='GALK.Can_abl',summary.func(growth.rate),by=c('GALK.allele',groupings2)]
DTSPLIT<-split(DT1,by='GALK.allele')
DT2<-merge(DTSPLIT[[1]],DTSPLIT[[2]],by=groupings2)
qplot(DT2$mean.x,DT2$mean.y)+geom_abline(intercept=0,slope=1)+ylab(names(DTSPLIT)[2])+xlab(names(DTSPLIT)[1])

DT1<-DT0[GALK.clone.named=='GALK.Can_abl',summary.func(yfp.mean.gal),by=c('GALK.allele',groupings2)]
DTSPLIT<-split(DT1,by='GALK.allele')
DT2<-merge(DTSPLIT[[1]],DTSPLIT[[2]],by=groupings2)
qplot(DT2$mean.x,DT2$mean.y)+geom_abline(intercept=0,slope=1)+ylab(names(DTSPLIT)[2])+xlab(names(DTSPLIT)[1])

DT1<-DT0[GALK.clone.named=='GALK.Can_abl',summary.func(yfp.mean.glu),by=c('GALK.allele',groupings2)]
DTSPLIT<-split(DT1,by='GALK.allele')
DT2<-merge(DTSPLIT[[1]],DTSPLIT[[2]],by=groupings2)
qplot(DT2$mean.x,DT2$mean.y)+geom_abline(intercept=0,slope=1)+ylab(names(DTSPLIT)[2])+xlab(names(DTSPLIT)[1])

DT1<-DT0[GALK.clone.named=='GALK.Can_abl',summary.func(fracon.gal),by=c('GALK.allele',groupings2)]
DTSPLIT<-split(DT1,by='GALK.allele')
DT2<-merge(DTSPLIT[[1]],DTSPLIT[[2]],by=groupings2)
qplot(DT2$mean.x,DT2$mean.y)+geom_abline(intercept=0,slope=1)+ylab(names(DTSPLIT)[2])+xlab(names(DTSPLIT)[1])

DT1<-DT0[GALK.clone.named=='GALK.Can_abl',summary.func(fracon.glu),by=c('GALK.allele',groupings2)]
DTSPLIT<-split(DT1,by='GALK.allele')
DT2<-merge(DTSPLIT[[1]],DTSPLIT[[2]],by=groupings2)
qplot(DT2$mean.x,DT2$mean.y)+geom_abline(intercept=0,slope=1)+ylab(names(DTSPLIT)[2])+xlab(names(DTSPLIT)[1])

# based on all this I suspect that Clone pAMN51.1 - 1 - A7 is the true WT. There is something wrong with the Clone pAMN51.2 - 1 - A8

###################################################
### quick look at GALK.Can_alb: clone pAMN51.1 - 1 - A7 is problematic
###################################################

DTFINAL<-DT0[GALK.allele!='pAMN51.2 - 1 - A8']

write.table(DTFINAL,paste0(outputdata,date,'-BatchEffectsFixedAndClonesCorrected_data.txt'),sep='\t',row.names=F)
fread(paste0(outputdata,date,'-BatchEffectsFixedAndClonesCorrected_data.txt'),sep='\t',header=T)








