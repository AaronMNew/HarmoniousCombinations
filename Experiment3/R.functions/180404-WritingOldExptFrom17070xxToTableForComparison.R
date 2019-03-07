#########################################################################################################################################################
### Old dataset
#########################################################################################################################################################	

source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'
source.code.path.file2<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis/R.functions/1707xx-180704-ReanalysisCustomScripts.R'
head.dir<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis'
layout.dir<-'./layout.sample.data'
date<-'1707xx' # beginning date clones were measured
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/'
outputdata<-'./output.data/'
setwd(head.dir)
source(source.code.path.file)
source(source.code.path.file2)
rm(letters)

br<-seq(-0.5,2.7,length=61)
library(tidyr)


###################################################
### load data
###################################################
DT00<-fread(paste0(outputdata,'1707xx-180704-OutliersRemoved.txt'),sep='\t',header=T)
# load "clusterframe" a data table that was made of the hdbscan expression density clusters
expression_clusters<-paste0(outputdata,'1707xx-180704-HDB_Clusters_Based_On_GAL_Distributions.rData')
load(expression_clusters) # loads "clusterframe" data.table			


###################################################
### write overlapping data to new table
###################################################

# pull in columns to keep
keepcolDT1<-fread(paste0(layout.dir,'/1707xx-180710-ColumnsToKeep.txt'),sep='\t',header=T)
keepcolDT<-keepcolDT1[keepcols==T]

phenocols<-c('growth.rate','phenotypic.index','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')
phenosh<-c('GR','PI','fONglu','fONgal','ymglu','ymgal')

dens.cols<-grepincols(DT00,'dens')
phenDT<-CloneNamedSummary(DT00[,c('clone.named','allele.named',phenocols),with=F])
densDT<-CloneNamedSummary(DT00[,c('clone.named','allele.named',dens.cols),with=F])
genocols<-colnames(CloneNamedSummary(DT00[,c('clone.named','allele.named'),with=F]))

DTX<-merge(phenDT,clusterframe,by=c('aa','bb','cc','clone.named'))
aas<-c('GAL4.WT','GAL4.delta','GAL4-L868P','GAL4-L868G','GAL4-L868K','GAL4-L868C')
bbs<-c('GAL80.WT','GAL80.delta','GAL80.07','GAL80S-0','GAL80S-1','GAL80S-2')
ccs<-c('GAL3.WT','GAL3.delta')
bb_clusts<-'inducible'
dSub<-DTX[aa%in%aas|cc%in%ccs&(bb=='GAL3.WT'|bb=='GAL3.delta')]
write.table(dSub,file=paste0(outputdata,'180906-SubsetOfDataForComparisonWith180404Data.txt'),sep='\t',row.names=F)
'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis/output.data/180906-SubsetOfDataForComparisonWith180404Data.txt'



#########################################################################################################################################################
### New experiment
#########################################################################################################################################################

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
source.code.path.file2<-'/Users/anew/Google Drive/002_lehnerlabData/180404-VarGALK.With.G80xG4_2/R.functions/180404-CustomScripts.R'
source(source.code.path.file2)
writeplot<-F
br<-seq(-0.5,2.7,length=61)
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

# BELOW I classified genotypes based on background switching rate in the absence of GAL sensing.
# "inducible" is TRUE when the grwoth rate is faster than expectation (derived from no-sensor backgrounds) of given genotype's expression in glucose
fname<-paste0(outputdata,'180404-180904-InducibleGenotypesBasedOnGluExpression.rdata')
load(fname)
# InducibleFrame # <- name of data.table

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

phenocols<-c('growth.rate','phenotypic.index','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')
dens.cols<-grepincols(DT00,'dens')
phenDT<-DT0[,c(genos,phenocols),with=F]
densDT<-DT0[,c(genos,dens.cols),with=F]

###################################################
### compare results with first experiment
###################################################

phenSplit<-function(DT,varcol='variable',splitfac='_'){
	DT<-data.table(DT)
	setnames(DT,varcol,'variable')
	DT[,c('pheno','stat'):=colsplit(variable,splitfac,c('pheno','stat'))]
	DT[,variable:=NULL]
	data.table(data.frame(DT))
}

genosfirst<-c('aa','bb','cc')
# GAL80.07 behaved differently here compared to previous expt
FirstExperimentDataLoc<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis/output.data/180906-SubsetOfDataForComparisonWith180404Data.txt'
firstexpt<-fread(FirstExperimentDataLoc)
# load "clusterframe" a data table that was made of the hdbscan expression density clusters
expression_clusters<-paste0('/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis/output.data/1707xx-180704-HDB_Clusters_Based_On_GAL_Distributions.rData')
load(expression_clusters) # loads "clusterframe" data.table			

summ1<-summary.func.all(firstexpt,phenocols,genosfirst)
mDT1<-phenSplit(melt(summ1,id=genosfirst))
mDT1[,expt:='a']
summ2<-summary.func.all(phenDT[dd=='GALK.Sac_cer'],phenocols,genos)[,!'dd']
mDT2<-phenSplit(melt(summ2,id=genosfirst))
mDT2[,expt:='b']

form1<-castform(c(genosfirst,'pheno'),c('expt','stat'))
cDT1<-dcast(mDT1,form1,value.var='value')
cDT2<-dcast(mDT2,form1,value.var='value')
DTm<-merge(cDT1,cDT2,by=c(genosfirst,'pheno'))
DTm[,G80:=cc=='GAL80.07']
DTm[,varexp:=varexplained(a_mean,b_mean),by='pheno']
DTm[,pct:=percent(round(varexp,2))]
DTm[,fac:=applyPaste(data.frame(gsub('\\.',' ',pheno),pct),'\nvar. explained: ')]
DTm[,varexp2:=summary(lm(a_mean~b_mean))$r.squared,by='pheno']
DTm[,pct2:=percent(round(varexp2,2))]
DTm[,fac2:=applyPaste(data.frame(fac,pct2),'\nR2: ')]

sample(LETTERS,4)
pOKJG<-ggplot(DTm,aes(a_mean,b_mean,xmax=a_upper,xmin=a_lower,ymax=b_upper,ymin=b_lower))+geom_abline(col='grey70')+theme_minimal()+
	geom_errorbar(col='grey70')+geom_errorbarh(col='grey70')+geom_point(shape=21)+scale_colour_manual(values=c('black'))+
	ylab('phenotype value for second experiment')+xlab('phenotype value for first experiment')+
	facet_wrap(~fac2,scales='free')
w<-7;h<-4
ggsave(paste0(figout,date,'-180906-pOKJG-ComparisonOfFirstExperimentWithSecond-AllPhenotypes.png'), pOKJG,width=w,height=h)



###################################################
### compare relationship between growth rate and glucose expression predicted growth rate and yfp.mean.glu from 
###################################################

# annotate a smaller version of DT0
dSub1<-DT0[,.(aa,cc,bb,dd,yfp.mean.glu,growth.rate)]
dum<-data.table(dd=c('GALK.Esc_col','GALK.Can_alb','GALK.Sac_cer','HIS5.Sch_pom'),GALK.cat=c('GALK (E. coli or C. albicans)\nsensor(-) kinase(+)','GALK (E. coli or C. albicans)\nsensor(-) kinase(+)','GAL1\nsensor(+) kinase(+)','HIS5\nsensor(-) kinase(-)'))
dum[,GALK.cat:=factor(GALK.cat,levels=unique(dum$GALK.cat))]
dSub<-merge(dSub1,dum,by='dd')
dSub[,x:= yfp.mean.glu]
dSub[,y:=growth.rate]
dSub[GALK.cat=='GALK (E. coli or C. albicans)\nsensor(-) kinase(+)'&bb=='GAL3.delta',fac:='no sensors']


# fit model and calculate deviation from expectation
dat<-dSub[fac=='no sensors',.(x,y)]
mod<-nls(y~SSlogis(x,Asym,xmid,scal),data=dat)
xes<-data.table(x=seq(1,50000,length=100))
lineframe<-data.table(xes,y=predict(mod,xes))
datcoli<-dSub[dd=='GALK.Esc_col'&bb=='GAL3.delta',.(x,y)]
modcoli<-nls(y~SSlogis(x,Asym,xmid,scal),data=datcoli)
datalb<-dSub[dd=='GALK.Can_alb'&bb=='GAL3.delta',.(x,y)]
modalb<-nls(y~SSlogis(x,Asym,xmid,scal), data =datalb)
dSub[dd=='GALK.Sac_cer'|dd=='HIS5.Sch_pom',pred:=predict(mod)]
dSub[dd=='GALK.Esc_col',pred:=predict(modcoli)]
dSub[dd=='GALK.Can_alb',pred:=predict(modalb)]


genos<-c('aa','bb','cc','clone.named')
DTX<-merge(firstexpt,clusterframe,by=c(genos,grepincols(firstexpt,c('HDB','clust'))))
summ<-CloneNamedSummary(summary.func.all(DTX,phenocols,c('clone.named',grepincols(DTX,c('HDB','clust')))),allele=F)
summ[,`single mutant\nlocus`:=mutant_gene]
summ[mutant_gene=='GAL3 GAL80 GAL4',`single mutant\nlocus`:='WT']
loessmod<-loess(b_mean~a_mean,DTm[pheno=='yfp.mean.glu'])
summ[,a_mean:=yfp.mean.glu_mean]
summ[,loessmod_mean:=predict(loessmod,a_mean)]
sample(LETTERS,4)
summ[,`HDB phenotype\nclassification`:=HDB2]
pBZRK<-ggplot(summ,aes(loessmod_mean,growth.rate_mean))+geom_point(aes(col=`HDB phenotype\nclassification`),shape=21,size=0.5)+
	geom_line(data=lineframe,aes(x,y),inherit.aes=F)+
	xlab('loess-corrected mean YFP signal in glucose')+ylab(expression(paste('growth rate ',mu,' ',hr^-1)))+
	ggtitle('Growth-rate-in-absence-of-sensor-prediction from second experiment\ncompared to growth rate and mean YFP sinal from first experiment')+
	theme(plot.title = element_text(size=10,face='plain'))
w<-5.5;h<-3
ggsave(paste0(figout,date,'-180906-pOKJG-ComparisonOfFirstExperimentWithSecond-GrowthRatePredictionBasedOnMeanYFPGluExpression.png'), pBZRK,width=w,height=h)









