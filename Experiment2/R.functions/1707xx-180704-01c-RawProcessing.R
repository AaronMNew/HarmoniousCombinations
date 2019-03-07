

# Scripts for dealing with flow cytometry data.
source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'

head.dir<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis'
date<-'1707xx' # beginning date clones were measured
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/'
outputdata<-'./output.data/'
layoutdir<-'./layout.sample.data/'
setwd(head.dir)

br<-seq(-0.5,2.7,length=61)

# this file generates these objects which can be used in downstream analysis


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


