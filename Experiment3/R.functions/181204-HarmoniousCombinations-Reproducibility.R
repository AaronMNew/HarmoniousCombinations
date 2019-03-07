# Scripts for dealing with flow cytometry data.
source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'
source.code.path.file2<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis/R.functions/1707xx-180704-ReanalysisCustomScripts.R'
aesthetics.dir<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/181119-Aesthetics_Harmonious_Combinations.R'
source(source.code.path.file)
source(source.code.path.file2)
source(aesthetics.dir)
library(tidyr)
writeplot<-F # when true runs all plotting

###################################################
### Experiment 1
###################################################


	head.dir<-"/Users/anew/Google Drive/002_lehnerlabData/17092x-CopyNumberExperiment"
	plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
	figout<-'./figure.output/publication/'
	outputdata<-'./output.data/'
	layoutdir<-'./layout.sample.data/'
	date<-'170927' # beginning date clones were measured
	setwd(head.dir)
	br<-seq(-0.5,2.7,length=61)
	
	# this file generates these objects which can be used in downstream analysis
	# merged.final<-fread(paste(date,"-DatasetWithGenerationsAndGeneExpression.txt",sep=''),sep='\t',header=T)
	
	###################################################
	### custom functions
	###################################################
	
	reverse_lev<-function(x,cust_order=NULL){
		if(!is.factor(x)){stop('\nx not a factor')}
		if(is.null(cust_order))cust_order<-length(levels(x)):1
		if(!is.null(cust_order))if(length(unique(cust_order))!=length(levels(x)))stop('\nyou screwed up the custom order')
		return(factor(x,levels=levels(x)[cust_order]))
	}
	
	
	###################################################
	### load data
	###################################################
	DT00<-fread(paste0(outputdata,'17092x-Dataset.txt'),sep='\t',header=T)
	DT0<-DT00[genotype.strain!='KV447'&source.plate%in%c("Plate_073",'Plate_074')]
	
	
	phenocols<-c('growth.rate','phenotypic.index','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal','cell.per.ml.glu','cell.per.ml.gal')
	ensembles<-c('plasgeno','binarycount','totalcopies')
	dens.cols<-grepincols(DT0,'dens')
	factorDefault<-c('orange1','cornflowerblue','darkgreen','red','chartreuse3')
	DT0[,c('a','b','c'):=list(1-total.copies.GAL4,1-total.copies.GAL3,1-total.copies.GAL80)]
	DT0[,WT:=plasgeno=='0 0 0 1 1 1']
	genos<-c(ensembles,letters[1:3],'WT')
	dum1<-list(data.table(a=c(1,0),aa=c('GAL4.delta','GAL4.WT')),data.table(b=c(1,0), bb=c('GAL3.delta','GAL3.WT')),data.table(c=c(1,0), cc=c('GAL80.delta','GAL80.WT')))
	phenDT<-merge_recurse2(DT0[plasgeno=='0 0 0 1 1 1'|genomecopies=='0 0 0',c(genos,phenocols),with=F],dum1,by.list=c('a','b','c'))
	phenDT[,clone.named:=applyPaste(data.frame(bb,cc,aa),' ')]
	phenDT[plasgeno=='0 0 0 1 1 1',clone.named:='WT_control']
	
	phenocols<-c('growth.rate','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')
	Exp1Phen<-copy(phenDT[,c('clone.named',phenocols),with=F])
	Exp1Phen[,expt:='Experiment 1']


###################################################
### Experiment 2
###################################################

	head.dir<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis'
	layout.dir<-'./layout.sample.data'
	date<-'1707xx' # beginning date clones were measured
	plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
	figout<-'./figure.output/publication/'
	outputdata<-'./output.data/'
	setwd(head.dir)
	br<-seq(-0.5,2.7,length=61)
	
	###################################################
	### load data
	###################################################
	DT00<-fread(paste0(outputdata,'1707xx-180704-OutliersRemoved.txt'),sep='\t',header=T)
	
	###################################################
	### cust functions
	###################################################
	
	denscolconv<-function(DT, vnull =T){
		DT[,variable:=gsub('dens__','',variable)]
		DT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
		DT[,sugar:=paste0('g',sugar)]
		if(vnull ==T)DT[,variable:=NULL]
	}
	
	
	
	###################################################
	### set variables and add labels
	###################################################
	
	# pull in columns to keep
	keepcolDT1<-fread(paste0(layout.dir,'/1707xx-180710-ColumnsToKeep.txt'),sep='\t',header=T)
	keepcolDT<-keepcolDT1[keepcols==T]
	
	phenocols<-c('growth.rate','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')
	# phenosh<-c('GR','PI','fONglu','fONgal','ymglu','ymgal')
	factorDefault<-c('orange1','cornflowerblue','darkgreen','red','chartreuse3')
	
	dens.cols<-grepincols(DT00,'dens')
	phenDT<-CloneNamedSummary(DT00[,c('clone.named','allele.named',phenocols),with=F])
	densDT<-CloneNamedSummary(DT00[,c('clone.named','allele.named',dens.cols),with=F])
	genocols<-colnames(CloneNamedSummary(DT00[,c('clone.named','allele.named'),with=F]))
	
	Exp2Phen<-copy(phenDT[!is.na(clone.named),c('clone.named',phenocols),with=F])
	Exp2Phen[,expt:='Experiment 2']

###################################################
### Experiment 3
###################################################

	head.dir<-"/Users/anew/Google Drive/002_lehnerlabData/180404-VarGALK.With.G80xG4_2"
	date<-'180404' # beginning date clones were measured
	plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
	figout<-'./figure.output/publication/'
	outputdata<-'./output.data/'
	layoutdir<-'./layout.sample.data/'
	setwd(head.dir)
	br<-seq(-0.5,2.7,length=61)
	
	###################################################
	### custom functions
	###################################################
	
	reverse_lev<-function(x,cust_order=NULL){
		if(!is.factor(x)){stop('\nx not a factor')}
		if(is.null(cust_order))cust_order<-length(levels(x)):1
		if(!is.null(cust_order))if(length(unique(cust_order))!=length(levels(x)))stop('\nyou screwed up the custom order')
		return(factor(x,levels=levels(x)[cust_order]))
	}
	
	###################################################
	### read data for samples 
	###################################################
	
	DT.merge<-fread(paste0(outputdata,date,'-BatchEffectsFixedAndClonesCorrected_data.txt'),sep='\t',header=T)
	DT.merge[GALK.clone.named=='GALK.Can_abl',GALK.clone.named:='GALK.Can_alb']
	DT00<-DT.merge[order(measurement.date.glu,samp.plate.glu,r)]
	genos<-c('aa','bb','cc','dd')
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
	### get orders right for clone.named factors
	###################################################
	
	G3ctrl<-c('GAL3.delta','GAL3.WT')
	G80ctrl<-c('GAL80.delta','GAL80.WT','GAL80.07','GAL80S-0','GAL80S-2','GAL80S-1')
	G4ctrl<-c('GAL4.delta','GAL4.WT','GAL4-L868C','GAL4-L868G','GAL4-L868K','GAL4-L868P')
	GKctrl<-c('HIS5.Sch_pom','GALK.Sac_cer')
	
	
	GAL3s<-factor(unique(c(DT0$GAL3.clone.named,G3ctrl)),levels=c(G3ctrl,unique(DT0$GAL3.clone.named%w/o%G3ctrl)))
	GAL80s<-factor(unique(c(DT0$GAL80.clone.named,G80ctrl)),levels=c(G80ctrl,unique(DT0$GAL80.clone.named%w/o%G80ctrl)))
	GAL4s<-factor(unique(c(DT0$GAL4.clone.named,G4ctrl)),levels=c(G4ctrl,unique(DT0$GAL4.clone.named%w/o%G4ctrl)))
	
	DT0$GAL3.clone.named <-factor(DT0$GAL3.clone.named,levels=c(G3ctrl,unique(DT0$GAL3.clone.named%w/o%G3ctrl)))
	DT0$GAL80.clone.named <-factor(DT0$GAL80.clone.named,levels=c(G80ctrl,unique(DT0$GAL80.clone.named%w/o%G80ctrl)))
	DT0$GAL4.clone.named <-factor(DT0$GAL4.clone.named,levels=c(G4ctrl,unique(DT0$GAL4.clone.named%w/o%G4ctrl)))
	DT0$GALK.clone.named <-factor(DT0$GALK.clone.named,levels=c(GKctrl,unique(DT0$GALK.clone.named%w/o%GKctrl)))
	
	DT0[,aa:=GAL4.clone.named]
	DT0[,bb:=GAL3.clone.named]
	DT0[,cc:=GAL80.clone.named]
	DT0[,dd:=GALK.clone.named]
	
	phenocols<-c('growth.rate','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')
	dens.cols<-grepincols(DT00,'dens')
	phenDT<-DT0[,c(genos,phenocols),with=F]
	phenDT[,clone.named:=applyPaste(data.frame(bb,cc,aa),' ')]
	
	Exp3Phen<-copy(phenDT[dd =='GALK.Sac_cer',c('clone.named',phenocols),with=F])
	Exp3Phen[,expt:='Experiment 3']



###################################################
### reproducibility within experiments
###################################################

DT0<-rbind(Exp1Phen,Exp2Phen,Exp3Phen)
bycols<-c('clone.named','expt') 
mDT<-melt(DT0,id=bycols)
mDT[,scaled:=scale(value),by=c('variable','expt')]
summ<-mDT[,summary.func(value),by=c(bycols,'variable')]
DTm<-merge(mDT,summ,by=c('clone.named','variable','expt'))
phenolevels<-data.table(
	pheno=c('growth.rate','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal'),
	variable=c('growth.rate','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal'),
	p1=factor(c('growth rate','fraction ON glucose','fraction ON galactose','mean YFP expression\nglucose, AU','mean YFP expression\ngalactose, AU'),
		levels=c('growth rate','fraction ON glucose','fraction ON galactose','mean YFP expression\nglucose, AU','mean YFP expression\ngalactose, AU')))


dum1<-merge(phenolevels,DTm[,round(varexplained(value,mean),3),by=c('expt','variable')])[order(p1)]
dum1[,phenotype:=applyPaste(data.frame(p1,V1),collapse='\nVarExp= ')]
dum1[,phenotype:=factor(phenotype,levels=dum1$phenotype)]
DTm2<-merge(DTm,dum1,by=c('variable','expt'))
DTm2[,alpha:=1/.N+0.1,by=c('expt','variable')]

sample(LETTERS,4)
pFTIQ<-ggplot(DTm2,aes(mean, value))+
	geom_point(shape=21,size=0.2,aes(alpha=alpha))+
	facet_wrap(expt~phenotype,scales='free',ncol=5)+
	theme_NewPub+
	geom_abline(col='grey70')+theme(legend.position='none')+
	xlab('within-genotype mean')+ylab('individual observations')

if(writeplot==T){
	w<-8.5*0.8;h<-6.4*0.8
	ggsave(paste0(figout,date,'-181203-pFTIQ-IntraExperimentalReproducibility.pdf'), pFTIQ,width=w,height=h)	
}

###################################################
### reproducibility between experiments
###################################################

phenocols<-c('growth.rate','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')

summ1<-summary.func.all(Exp1Phen[clone.named!='WT_control'],phenocols,bylist=c('clone.named'))
summ2<-summary.func.all(Exp2Phen,phenocols,bylist=c('clone.named'))
summ3<-summary.func.all(Exp3Phen,phenocols,bylist=c('clone.named'))

mDT1<-melt(summ1,id='clone.named')
mDT2<-melt(summ2,id='clone.named')
mDT3<-melt(summ3,id='clone.named')

E1E2<-merge(mDT1,mDT2,by=c('clone.named','variable'))
E1E3<-merge(mDT1,mDT3,by=c('clone.named','variable'))
E2E3<-merge(mDT2,mDT3,by=c('clone.named','variable'))

E1E2[,expt:='Experiment 1 vs Experiment 2']
E1E3[,expt:='Experiment 1 vs Experiment 3']
E2E3[,expt:='Experiment 2 vs Experiment 3']

DTL<-list(E1E2,E1E3,E2E3)

DT1<-data.table(ldply(lapply(DTL,function(DTx){
	#DTx<-DTL[[1]]
	DTx[,c('pheno','stat'):=colsplit(variable,'_',c('pheno','stat'))]
	setnames(DTx,c('value.x','value.y'),c('x','y'))
	form<-as.formula(castform(c('expt','clone.named','pheno'),c('stat')))
	cDT<-dcast(data=DTx,formula=form,value.var=c('x','y'))
	dum1<-merge(cDT[,round(varexplained(x_mean,y_mean),3),by=c('pheno')],phenolevels,by='pheno')[order(p1)]
	dum1[,phenotype1:=applyPaste(data.frame(p1,V1),collapse='\nVarExp= ')]
	dum1[,phenotype:=factor(phenotype1,levels=dum1$phenotype1)]
	merge(cDT,dum1[,.(pheno,phenotype)],by='pheno')
})))

sample(LETTERS,4)

pMFKB<-ggplot(DT1,aes(x_mean,y_mean,xmin=x_lower,xmax=x_upper,ymin=y_lower,ymax=y_upper))+
	geom_errorbar(col='grey70')+geom_errorbarh(col='grey70')+geom_point(shape=21,col='black')+
	geom_abline(col='grey70')+theme_NewPub+xlab('experiment x')+ylab('experiment y')+
	facet_wrap(expt~phenotype,ncol=5,scales='free')


if(writeplot==T){
	w<-8.5*0.8;h<-6.4*0.8
	ggsave(paste0(figout,date,'-181203-pMFKB-InterExperimentalReproducibility.pdf'), pMFKB,width=w,height=h)	
}




























