
# Scripts for dealing with flow cytometry data.
source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'
source.code.path.file2<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis/R.functions/1707xx-180704-ReanalysisCustomScripts.R'

head.dir<-"/Users/anew/Desktop/17092x-CopyNumberExperiment"
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/'
outputdata<-'./output.data/'
layoutdir<-'./layout.sample.data/'
date<-'170927' # beginning date clones were measured
setwd(head.dir)
source(source.code.path.file)
source(source.code.path.file2)

br<-seq(-0.5,2.7,length=61)

# this file generates these objects which can be used in downstream analysis
# merged.final<-fread(paste(date,"-DatasetWithGenerationsAndGeneExpression.txt",sep=''),sep='\t',header=T)


###################################################
### load data
###################################################

setwd(head.dir)
datloc<-paste(outputdata,date,'-analyzed.tab',sep='')
nullsamploc<-paste(outputdata,date,'-null.samps.tab',sep='')
# analyzed data have already been written from the Preprocessing.R file and here are fread output written the first time around.
output<-fread(paste(outputdata,date,'-analyzed.tab',sep=''),header=TRUE,sep='\t')
old.cols<-colnames(output)[c(grep("\\]",colnames(output)))]
dens.cols<-paste('dens',round(br[1:(length(br)-1)],2),sep='.')
setnames(output,old= old.cols,new= dens.cols)

nullsamps.summary<-NULL
nullsamps.summary<-fread(nullsamploc,header=TRUE,sep='\t')
# layout is a table with headers and at least a "samp_id" column either generated in the layout.maker script or by hand in excel.
layout<-fread(paste(layoutdir,'17092x-Layout.txt',sep=''),header=TRUE,sep='\t')

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


merged<-merge(output.final,layout,by='samp_id')


###################################################
### add a few things: compute number of generations + phenotypic index
###################################################

dens.norm<-data.table(t(apply(merged[, dens.cols,with=FALSE],1,function(x)round(row.normalization(x),3))))
merged[,dens.cols]<-dens.norm

merged[,c('measurement.date','samp.plate.rc'):=colsplit(samp_id,'-',c('measurement.date','samp.plate.rc'))]
merged[order(samp_id),arb:=1:.N,by=measurement.date]

# for cell density per ml, cells.per.sec * 1/µl per second sampling rate * 1000 µl / ml
dim(merged)

merged[,cell.per.ml:=cells.per.sec*(1/as.numeric(sample_flow_rate))*1000]
DT.glu<-merged[glu==1]
DT.gal<-merged[gal==1]

# to couple direct measurements with each other (glucose > galactose) i need to create an id variable that will allow joining glu with gal measurements
DT.glu[,next.day:=as.numeric(measurement.date)+1]
DT.gal[,next.day:=NA]
DT.glu[,id:=arb]
DT.gal[,id:=arb]
DT.gal[,gal:=NULL]
DT.gal[,glu:=NULL]
DT.glu[,gal:=NULL]
DT.glu[,glu:=NULL]

#variables definition
dil.factor<-9/150
hours <-12
galphenos<-colnames(output)
# i re-merge them here and then calculate generations, growth rate and phenotypic index
phenos1<-c(colnames(output.final),'cell.per.ml','measurement.date','plate')
setnames(DT.gal,phenos1,paste0(phenos1,'.gal'))
setnames(DT.glu,phenos1,paste0(phenos1,'.glu'))
DT.merge1<-merge(DT.glu[,!colnames(layout),with=F],DT.gal,by=c('id'),suffixes=c('','___y'),all=TRUE)
DT.merge<-DT.merge1[,!grepincols(DT.merge1,'___'),with=F]
DT.merge[,generations:=log((cell.per.ml.gal/((cell.per.ml.glu)*dil.factor)),base=2)]
DT.merge[,growth.rate:=log((cell.per.ml.gal/((cell.per.ml.glu)*dil.factor)))/hours]
DT.merge[,phenotypic.index:=fracon.glu+fracon.gal]
mDT<-melt(DT.merge[,c('generations','growth.rate','phenotypic.index')])
ggplot(mDT,aes(value))+geom_histogram()+facet_grid(~variable,scales='free')

write.table(DT.merge,paste0(outputdata,'17092x-Dataset.txt'),row.names=F,sep='\t')
