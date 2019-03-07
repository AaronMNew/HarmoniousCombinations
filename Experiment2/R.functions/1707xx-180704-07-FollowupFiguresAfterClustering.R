# Scripts for analysis.
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

phenocols<-c('growth.rate','phenotypic.index','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')
phenosh<-c('GR','PI','fONglu','fONgal','ymglu','ymgal')

dens.cols<-grepincols(DT00,'dens')
phenDT<-CloneNamedSummary(DT00[,c('clone.named','allele.named',phenocols),with=F])
densDT<-CloneNamedSummary(DT00[,c('clone.named','allele.named',dens.cols),with=F])
genocols<-colnames(CloneNamedSummary(DT00[,c('clone.named','allele.named'),with=F]))

###################################################
### remake fracon gal fracon.glu plots
###################################################
genos<-c('aa','bb','cc','clone.named')
DTX<-merge(phenDT,clusterframe,by=genos)
summ<-CloneNamedSummary(summary.func.all(DTX,phenocols,c('clone.named',grepincols(DTX,c('HDB','clust')))),allele=F)
summ[,`single mutant\nlocus`:=mutant_gene]
summ[mutant_gene=='GAL3 GAL80 GAL4',`single mutant\nlocus`:='WT']

ggplot(summ,aes(yfp.mean.glu_mean,growth.rate_mean))+geom_point(aes(col=HDB2),size=0.1)

sample(letters,4)
DT<-summ
pGAVJ<-ggplot(DT,aes(fracon.glu_mean,fracon.gal_mean,ymax=fracon.gal_upper,ymin=fracon.gal_lower,xmax=fracon.glu_upper,xmin=fracon.glu_lower))+
	geom_abline(col='red')+
	geom_errorbar(col='grey70',alpha=0.3)+geom_errorbar(col='grey70',alpha=0.3)+geom_point(col='grey70',alpha=0.3)+
	geom_errorbar(data= DT[single_lv==T],aes(ymax=fracon.gal_upper,ymin=fracon.gal_lower,col=`single mutant\nlocus`))+
	geom_errorbarh(data= DT[single_lv==T],aes(xmax=fracon.glu_upper,xmin=fracon.glu_lower,col=`single mutant\nlocus`))+
	geom_point(data= DT[single_lv==T],aes(fracon.glu_mean,fracon.gal_mean,col=`single mutant\nlocus`),shape=21,size=3,stroke=1)+
	scale_colour_manual(values=c('cornflowerblue','orange1','darkgreen','black'))+
	ylab('fraction ON in galactose')+xlab('fraction ON in glucose')+
	xlim(c(0,1))+ylim(c(0,1))

DTX[,c('x','y')]
mod1<-lm(growth.rate~fracon.gal+fracon.glu,DTX)
summ[,c('growth.rate','fracon.gal','fracon.glu'):=list(growth.rate_mean,fracon.gal_mean,fracon.glu_mean)]
summ[,pred_mean:=predict(mod1,data.frame(growth.rate,fracon.gal,fracon.glu))]
summ[,pred_se:=predict(mod1,data.frame(growth.rate,fracon.gal,fracon.glu),se.fit=T)$se]

ggplot(DT,aes(pred_mean,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower,xmax= pred_mean+ pred_se,xmin= pred_mean-pred_se))+
	geom_abline(col='red')+
	geom_errorbar(col='grey70',alpha=0.3)+geom_errorbar(col='grey70',alpha=0.3)+geom_point(col='grey70',alpha=0.3)+
	geom_errorbar(data= DT[single_lv==T],aes(ymax=growth.rate_upper,ymin=growth.rate_lower,col=`single mutant\nlocus`))+
	geom_errorbarh(data= DT[single_lv==T],aes(xmax= pred_mean+ pred_se,xmin= pred_mean-pred_se,col=`single mutant\nlocus`))+
	geom_point(data= DT[single_lv==T],aes(pred_mean,growth.rate_mean,col=`single mutant\nlocus`),shape=21,size=3,stroke=1)+
	scale_colour_manual(values=c('cornflowerblue','orange1','darkgreen','black'))+
	ylab(expression(paste('observed growth rate ',mu,' ',hr^-1)))+xlab(expression(paste('predicted growth rate ',mu,' ',hr^-1)))+
	xlim(c(0,0.3))+ylim(c(0,0.3))

DT<-summ[cc_clust!='const, low expr'&aa_clust!='const, low expr']
DT<-summ[cc_clust!='const, low expr'&aa_clust!='const, low expr'&aa_clust!='weakly inducible']

pQUFD<-ggplot(DT,aes(fracon.glu_mean,fracon.gal_mean,ymax=fracon.gal_upper,ymin=fracon.gal_lower,xmax=fracon.glu_upper,xmin=fracon.glu_lower))+
	geom_abline(col='red')+
	geom_errorbar(col='grey70',alpha=0.3)+geom_errorbar(col='grey70',alpha=0.3)+geom_point(col='grey70',alpha=0.3)+
	geom_errorbar(data= DT[single_lv==T],aes(ymax=fracon.gal_upper,ymin=fracon.gal_lower,col=`single mutant\nlocus`))+
	geom_errorbarh(data= DT[single_lv==T],aes(xmax=fracon.glu_upper,xmin=fracon.glu_lower,col=`single mutant\nlocus`))+
	geom_point(data= DT[single_lv==T],aes(fracon.glu_mean,fracon.gal_mean,col=`single mutant\nlocus`),shape=21,size=3,stroke=1)+
	scale_colour_manual(values=c('cornflowerblue','orange1','darkgreen','black'))+
	ylab('fraction ON in galactose')+xlab('fraction ON in glucose')+
	xlim(c(0,1))+ylim(c(0,1))

DT<-summ[cc_clust!='const, low expr'&aa_clust!='const, low expr'&aa_clust!='weakly inducible']
DT[,varexplained(growth.rate_mean,pred_mean)]
summ[,varexplained(growth.rate_mean,pred_mean)]
ggplot(DT,aes(pred_mean,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower,xmax= pred_mean+ pred_se,xmin= pred_mean-pred_se))+
	geom_abline(col='red')+
	geom_errorbar(col='grey70',alpha=0.3)+geom_errorbar(col='grey70',alpha=0.3)+geom_point(col='grey70',alpha=0.3)+
	geom_errorbar(data= DT[single_lv==T],aes(ymax=growth.rate_upper,ymin=growth.rate_lower,col=`single mutant\nlocus`))+
	geom_errorbarh(data= DT[single_lv==T],aes(xmax= pred_mean+ pred_se,xmin= pred_mean-pred_se,col=`single mutant\nlocus`))+
	geom_point(data= DT[single_lv==T],aes(pred_mean,growth.rate_mean,col=`single mutant\nlocus`),shape=21,size=3,stroke=1)+
	scale_colour_manual(values=c('cornflowerblue','orange1','darkgreen','black'))+
	ylab(expression(paste('observed growth rate ',mu,' ',hr^-1)))+xlab(expression(paste('predicted growth rate ',mu,' ',hr^-1)))+
	xlim(c(0,0.3))+ylim(c(0,0.3))

w<-6.5;h<-5
ggsave(paste0(figout,'1707xx-180907-pGAVJ-FraconGluVsFraconGal.png'), pGAVJ,width=w,height=h)













