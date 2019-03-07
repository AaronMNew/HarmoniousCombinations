

# Scripts for dealing with flow cytometry data.
source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'
source.code.path.file2<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis/R.functions/1707xx-180704-ReanalysisCustomScripts.R'

head.dir<-"/Users/anew/Google Drive/002_lehnerlabData/17092x-CopyNumberExperiment"
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-paste0(head.dir,'/figure.output/publication/')
outputdata<-'./output.data/'
layoutdir<-'./layout.sample.data/'
date<-'170927' # beginning date clones were measured
setwd(head.dir)
source(source.code.path.file)
source(source.code.path.file2)
aesthetics.dir<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/181119-Aesthetics_Harmonious_Combinations.R'
source(aesthetics.dir)
writeplot=F
br<-seq(-0.5,2.7,length=61)

# this file generates these objects which can be used in downstream analysis
# merged.final<-fread(paste(date,"-DatasetWithGenerationsAndGeneExpression.txt",sep=''),sep='\t',header=T)


###################################################
### analysis: illustrative plots of initial filtering done on smaples based on size and shape characteristics determined by FCS and SSC channels.
###################################################

setwd(plates.dir)

dirs<-list.dirs()[grep('Plate',list.dirs())][1]
ix<-1 # this is normally a for loop which dumps a bunch of filtered flowCore data files in .rData format
wash=NULL # use if a blank sample was taken before every measurement
filtsize<-4 # filter size set for the curv2Filter
setwd(plates.dir) # directory where the plates are
setwd(dirs[ix]) # move to the directory where the plates are
dir<-paste(dirs[ix]) # make a list of all the directories with plates
setwd(plates.dir);setwd(dir) # move to the first directory in the list
fcs.files<-list.files(pattern='../*.fcs') # get all teh fcs files in the folder
# fsapply is normally used to go through each sample
flowframes<-read.flowSet(fcs.files, name.keyword="TUBE NAME", phenoData=list(name="TUBE NAME", Filename="$FIL")) # make a flowframe object of the fcs files
#	head(pData(phenoData(flowframes)))
fcs<-flowframes[[1]] # within fsapply, go through each item in the flow frame
if(!is.null(wash)){strt<-0} # establish a starting cell population 
if(is.null(wash)){strt<-500} # set a variable of number of cells to discard if the previous reading was not washed
len<-length(exprs(fcs$'FSC-A')) # determine the length of the list acquired
if(len<(strt+1000)){filteredPlus<-NULL # if you have fewer than 1000 cells plus the amount you are discarding then count the sample as a NULL
	}else{
	ssc<-quantile(exprs(fcs[strt:len,3]),c(0.1,0.9)) # determine the quantiles to filter the ssc
	fsc<-quantile(exprs(fcs[strt:len,1]),c(0.1,0.9)) # determine the quantiles to filter the fsc
	rectGate <- rectangleGate("FSC-H"=fsc,"SSC-H"=ssc) # create a gate for the samples based on these quantiles
	filteredRec<-split(fcs,flowCore::filter(fcs,rectGate))$`defaultRectangleGate+` # actually filter the samples based on the rectangular gate
	c2f <- curv2Filter("FSC-H","SSC-H",bwFac=filtsize) # take rectangle filtered samples and apply the curve2Filter gate
	filtered<-split(filteredRec,flowCore::filter(filteredRec,c2f)) # actually filter the samples based on the curve2Filter
	filteredPlus<-filtered$`area 1`} # extract the cells that wall within the central area of the gate

	filteredPlus # return the final filtered sample

# after filtering, samples are re-analyzed for expression distribution characteristics
# make a plot of the raw data through filtered final
DT0<-data.table(FSC=exprs(fcs$'FSC-A'),SSC=exprs(fcs$'SSC-A'),fac='original')
DTrect<-data.table(FSC=exprs(filteredRec $'FSC-A'),SSC=exprs(filteredRec $'SSC-A'),fac='rectangle filtered')
DTcurve<-data.table(FSC=exprs(filteredPlus $'FSC-A'),SSC=exprs(filteredPlus $'SSC-A'),fac='curve2Filt filtered')
filtlevs<-c('original','rectangle filtered','curve2Filt filtered')
DT<-rbind(DT0,DTrect,DTcurve)
DT[,`filtered data`:=factor(fac,levels=filtlevs)]
setnames(DT,c('FSC.FSC-A', 'SSC.SSC-A'),c('FSC','SSC'))
DT[,c('log10FSC','log10SSC'):=list(log10(FSC+1000)-3,log10(SSC+1000)-3)]

# plot filtered samples
p<-ggplot(DT,aes(log10FSC,log10SSC,col=`filtered data`))
p+geom_point(size=0.1)

library(sparr)
X<-log10(as.numeric(as.vector(unlist(DT0[,1])))+1000)-3
Y<-log10(as.numeric(as.vector(unlist(DT0[,2])))+1000)-3

dtp<-ppp(X,Y,window=owin(range(X),range(Y)))
test<-bivariate.density(dtp,h0=0.05)
zmat<-as.matrix(test$z)
x<-test$q$xcol
y<-test$q$yrow
DTsm<-data.table(Y=y,zmat)
DTm<-melt(DTsm,id='Y')
dum<-data.table(variable=paste0('V',seq(1,128)),X=x)
DTf<-merge(dum,DTm,by='variable')
sample(letters,4)

pIUDX<-p+geom_point(size=0.3,shape=21,stroke=0.1)+geom_contour(data=DTf,aes(x=X,y=Y,z=value),inherit.aes=F,col='black',alpha=0.2)+
	theme_NewPub+
	xlab('log10 FSC signal')+ylab('log10 SSC signal')

if(writeplot==T){
	w<-3.9;h<-2.85
	ggsave(paste0(figout,'17092x-190110-pIUDX-FilteringSSC,FSC,illustration.pdf'), pIUDX,width=w,height=h)
}


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




