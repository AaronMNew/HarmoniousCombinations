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
### geometric model: growth rate
###################################################

DT0<-copy(phenDT)
backg<-mean(DT0[grepl('GAL4.delta',clone.named)]$growth.rate,na.rm=T)
backsd<-sd(DT0[grepl('GAL4.delta',clone.named)]$growth.rate,na.rm=T)
credmin<-backg-backsd
pheno<-'growth.rate'
genos<-c('aa','bb','cc')


pheno<-'growth.rate'
summname<-'pheno'
DT0[,x:=lapply(.SD,function(x)x),.SDcols=pheno]
bygenos<-c('clone.named','aa','bb','cc')
summ<-DT0[,summary.func(x,summname),by=bygenos]
WT<-DT0[WT_lv==T,summary.func(x)]
refmu<-WT$mean-backg
backg<-mean(DT0[grepl('GAL4.delta',clone.named)]$x,na.rm=T)
backsd<-sd(DT0[grepl('GAL4.delta',clone.named)]$x,na.rm=T)

G4<-DT0[G4.single==T,summary.func(x),by='clone.named']
G3<-DT0[G3.single==T,summary.func(x),by='clone.named']
G80<-DT0[G80.single==T,summary.func(x),by='clone.named']
bycols<-colnames(DT0)
dt.list<-list(G4,G3,G80)
boo<-c('aa','bb','cc')
dt.list1<-lapply(1:length(dt.list),function(i){
	xdt<-dt.list[[i]]
	oldcols<-notin(colnames(xdt),'clone.named')
	newcols<-paste0(boo[i],'__',oldcols)
	setnames(xdt,oldcols,newcols)
	xdt2<-(CloneNamedSummary(xdt,allele=F))
	return(xdt2[,c(boo[i],newcols),with=F])
})
DTm<-merge_recurse3(summ,dt.list1,by.list=c('aa','bb','cc'))
DTm[,mult_pred_mean:=apply(data.frame(aa__mean,bb__mean,cc__mean),1,function(x){
	prod(x-backg)/refmu^2+backg
})]
DTm[,mult_pred_se:=apply(data.frame(mult_pred_mean,aa__SE_cov,bb__SE_cov,cc__SE_cov),1,function(x){
	abs(x[1])*sqrt(sum(x[2:4]^2))
})]
DTm[,mult_pred_N:=apply(data.frame(aa__N,bb__N,cc__N),1,function(x){
	sum(x)
})]
DTm[,c('mult_pred_upper','mult_pred_lower'):=list(mult_pred_mean+ mult_pred_se, mult_pred_mean-mult_pred_se)]

geom_model<-copy(DTm)

DTm1<-merge(DTm,DTm[,t.test2(pheno_mean, mult_pred_mean,pheno_se, mult_pred_se,pheno_N, mult_pred_N,list=T),by='clone.named'],by='clone.named')
DTm1[,p.adj:=p.adjust(p.value,'fdr')]
sigcutoff<-0.05
diffcutoff<-0.025
DTm1[,sig:=p.adj<sigcutoff&abs(diff.of.means)> diffcutoff& pheno_mean >credmin]

DTm2<-CloneNamedSummary(DTm1[,!c('aa','bb','cc')],allele=F)
datcols<-grepincols(DTm2,c('mult_pred','pheno'))
DTm2[,pairing:=gsub('\\ ',' v ',mutant_gene)]
pairlevs<-c('GAL3 v GAL80','GAL3 v GAL4','GAL80 v GAL4')
restlevs<-unique(notin(DTm2$pairing,pairlevs))
DTm2[,pairing:=factor(pairing,levels=c(pairlevs,restlevs))]
plotDT<-DTm2[!(single_lv==T|clone.named=='GAL3.WT GAL80.WT GAL4.WT')&!is.na(sig),c('clone.named','sig','pairing',datcols),with=F]
setnames(plotDT,datcols,gsub('pheno','y',gsub('mult_pred','x',datcols)))

plotDT[,`FDR < 0.05`:=sig]
sample(LETTERS,4)

pGQPU<-ggplot(plotDT[order(sig)],aes(x_mean,y_mean,xmax=x_upper,xmin=x_lower,ymax=y_upper,ymin=y_lower))+
	geom_errorbar(col='grey50',alpha=0.2)+geom_errorbarh(col='grey50',alpha=0.2)+geom_point(aes(col=`FDR < 0.05`),alpha=0.3)+
	geom_abline(col='indianred')+scale_colour_manual(values=c('cornflowerblue','orange1'))+
	xlab(expression(paste('predicted growth rate ',mu,' ',hr^-1)))+ylab(expression(paste('observed growth rate ',mu,' ',hr^-1)))+
	facet_grid(~ pairing)+theme_light()
# pGQPU
if(writeplot==T){
	w<-7;h<-2.5
	ggsave(paste0(figout,'1707xx-180907-pGQPU-MultiplicativePredictionExpectedVsObserved-FacetByPairingOnly.png'), pGQPU,width=w,height=h)	
}

p1<-ggplot(plotDT[order(sig)],aes(x_mean,y_mean,xmax=x_upper,xmin=x_lower,ymax=y_upper,ymin=y_lower))+geom_errorbar(col='grey50',alpha=0.2)+geom_errorbarh(col='grey50',alpha=0.2)+geom_point(aes(col=sig),alpha=0.3)+
	geom_abline(col='indianred')+scale_colour_manual(values=c('cornflowerblue','orange1'))+
	facet_grid(sig~ pairing)+theme_light()
# p1
#do.call(grid.arrange,list(p1,p0,ncol=1))

plotDTa<-copy(plotDT)
plotDTa[,sig1:=as.character(sig)][,sig1:='ALL TOGETHER']
plotDT[,sig1:=as.character(sig)]
plotDT2<-rbind(plotDT,plotDTa)
sample(letters,4)
pVLYF<-ggplot(plotDT2[order(sig)],aes(x_mean,y_mean,xmax=x_upper,xmin=x_lower,ymax=y_upper,ymin=y_lower))+
	geom_errorbar(col='grey50',alpha=0.2)+geom_errorbarh(col='grey50',alpha=0.2)+geom_point(aes(col=sig),alpha=0.3)+
	geom_abline(col='indianred')+scale_colour_manual(values=c('cornflowerblue','orange1'))+
	facet_grid(sig1~ pairing)+theme_light()+
	ylab(expression(paste('observed growth rate ',hr^-1)))+xlab(expression(paste('geometric model predicted growth rate ',hr^-1)))
# w<-6.5;h<-5
# ggsave(paste0(figout,'1707xx-180823-pVLYF-MultiplicativePredictionExpectedVsObserved.png'), pVLYF,width=w,height=h)

hdbtest<-hdbscan(plotDT[,.(x_mean,y_mean)],minPts=40)
plotDT[,cluster:=hdbtest$cluster]
clust<-c('unclustered','pred high & obs high','pred low & obs low','pred low & obs high')
dumm<-data.table(cluster=0:3,cluster_new=factor(clust,levels=clust[c(2:4,1)]))
plotDT3<-merge(plotDT,dumm,by='cluster')
sample(letters,4)
pKCGL<-ggplot(plotDT3[order(sig)],aes(x_mean,y_mean,xmax=x_upper,xmin=x_lower,ymax=y_upper,ymin=y_lower))+
	geom_errorbar(col='grey50',alpha=0.2)+geom_errorbarh(col='grey50',alpha=0.2)+geom_point(aes(col=sig),alpha=0.3)+
	geom_abline(col='indianred')+scale_colour_manual(values=c('cornflowerblue','orange1'))+
	facet_grid(cluster_new~ pairing)+theme_light()+
	ylab(expression(paste('observed growth rate ',hr^-1)))+xlab(expression(paste('geometric model predicted growth rate ',hr^-1)))
# w<-7.15;h<-5.9
# ggsave(paste0(figout,'1707xx-180823-pKCGL-HDBClustered-MultiplicativePredictionExpectedVsObserved.png'), pKCGL,width=w,height=h)


hdbtest<-hdbscan(plotDT[sig==T,.(x_mean,y_mean)],minPts=20)
plotDT[sig==T,cluster2:=as.factor(hdbtest$cluster)]
plotDT[cluster2==0|cluster2==1,cluster2:=NA]
plotDT[,`sig epist cluster`:=as.factor(cluster2)]
sample(letters,4)
pUZDI<-ggplot(plotDT[order(sig)],aes(x_mean,y_mean,xmax=x_upper,xmin=x_lower,ymax=y_upper,ymin=y_lower))+
	geom_errorbar(col='grey50',alpha=0.2)+geom_errorbarh(col='grey50',alpha=0.2)+geom_point(aes(col=`sig epist cluster`),alpha=0.3)+
	geom_abline(col='indianred')+#scale_colour_manual(values=c('cornflowerblue','orange1'))+
	facet_grid(~ pairing)+theme_light()+
	ylab(expression(paste('observed growth rate ',hr^-1)))+xlab(expression(paste('geometric model predicted growth rate ',hr^-1)))
# w<-8.1;h<-3.8
# ggsave(paste0(figout,'1707xx-180823-pUZDI-HDBClusteredInSignificantEpistasis-MultiplicativePredictionExpectedVsObserved.png'), pUZDI,width=w,height=h)



############
#### show enrichment of different clones for being involved with epistatic interactions
############
# not really very informative
DT11<-CloneNamedSummary(plotDT,allele=F)
mDT<-melt(DT11,id=notin(colnames(DT11),c('aa','bb','cc')))
mDT[,pairing:=as.character(pairing)]
dSub<-mDT[grepl('v',pairing)&!grepl('WT',value)]
summ1<-dSub[,.N,by=c('sig','value')][order(N)]
summ<-dSub[,.N,by=c('pairing','sig','value')]
allelelev<-summ1[sig==T]$value
summ[,allele:=factor(value,levels= allelelev)]
summ[,epistasis:=sig]
# ggplot(summ,aes(epistasis, allele))+geom_tile(aes(fill=N))+scale_fill_gradientn(colours=c('white','indianred'))+
	# theme(axis.text.y=element_text(size=6))+
	# facet_wrap(~pairing,scales='free',drop=T)


############
#### show underlying phenotypes behind epistasis clusters: raw phenotypes
############

test<-CloneNamedSummary(merge(plotDT3[,.(clone.named,cluster_new,pairing)],
	na.exclude(plotDT[,grepincols(plotDT,c('clone.named','sig epis','pairing')),with=F]),
	by=c('clone.named','pairing')),allele=F)
test2<-test[,grepincols(test,c('pairing','clust','aa','bb','cc')),with=F]
test3<-melt(test2,id=grepincols(test2,c('clust','pair')))
test4<-test3[!grepl('WT',value)]

singles<-phenDT[single_lv==T]
#singles[,c("growth.rate","phenotypic.index","fracon.glu","fracon.gal","yfp.mean.glu","yfp.mean.gal"):=lapply(.SD,zscore),.SDcols=phenocols]
singles2<-melt(singles[,notin(colnames(singles),c('aa1','bb1','cc1')),with=F],id= notin(colnames(singles),c('aa','bb','cc')))
singles3<-singles2[!grepl('WT',value)]

DTm<-merge(test4,singles3,by=c('variable','value'),allow.cartesian=T)
mDT<-melt(DTm,id=notin(colnames(DTm),phenocols),variable.name='phenotype',value.name='obs')
mDT[,locus:=gsub('aa','GAL4',gsub('bb','GAL3',gsub('cc','GAL80',variable)))]
mDT[,fac:=applyPaste(data.frame(pairing,locus),': ')]
mDT[,dup:=applyPaste(mDT)]
mDT1<-mDT[!duplicated(dup)];mDT1[,dup:=NULL]
sample(letters,4)

# all epistasis clusters
pKHWE<-ggplot(mDT1[cluster_new!='unclustered'],aes(fac,obs))+geom_jitter(alpha=0.2,width=0.2,shape=21,size=0.1,col='grey70')+
	geom_violin(aes(col=pairing),alpha=0.5)+coord_flip()+geom_hline(yintercept=0,col='grey50',alpha=0.5)+
	facet_grid(cluster_new~phenotype,scales='free')+theme(axis.text.x=element_text(angle=40,hjust=1))+xlab('observed phenotype')+ylab('')
# w<-12.4;h<-6.45
# ggsave(paste0(figout,'1707xx-180823-pKHWE-PhenotypesOfEpistasisClusters.png'), pKHWE,width=w,height=h)

# all clones within pred low obs high clusters
pIUXT<-ggplot(mDT1[cluster_new=='pred low & obs high'],aes(fac,obs))+geom_jitter(alpha=0.2,width=0.2,shape=21,size=0.1,col='grey70')+
	geom_violin(aes(col=pairing),alpha=0.5)+coord_flip()+geom_hline(yintercept=0,col='grey50',alpha=0.5)+
	facet_grid(cluster_new~phenotype,scales='free')+theme(axis.text.x=element_text(angle=40,hjust=1))+xlab('observed phenotype')+ylab('')
# w<-12.4;h<-3
# ggsave(paste0(figout,'1707xx-180823-pIUXT-PhenotypesOfEpistasisClusters-PredLowObsHigh.png'), pIUXT,width=w,height=h)

# significant epistasis clusters
pHANF<-ggplot(mDT1,aes(fac,obs))+geom_jitter(alpha=0.2,width=0.2,shape=21,size=0.1,col='grey70')+
	geom_violin(aes(col=pairing),alpha=0.5)+coord_flip()+geom_hline(yintercept=0,col='grey50',alpha=0.5)+
	facet_grid(`sig epist cluster`~phenotype,scales='free')+theme(axis.text.x=element_text(angle=40,hjust=1))

# w<-12.4;h<-6.7
# ggsave(paste0(figout,'1707xx-180823-pHANF-PhenotypesOfSignificantEpistasisClusters.png'), pHANF,width=w,height=h)

############
#### show underlying phenotypes behind epistasis clusters: scaled phenotypes
############

test<-CloneNamedSummary(merge(plotDT3[,.(clone.named,cluster_new,pairing)],
	na.exclude(plotDT[,grepincols(plotDT,c('clone.named','sig epis','pairing')),with=F]),
	by=c('clone.named','pairing')),allele=F)
test2<-test[,grepincols(test,c('pairing','clust','aa','bb','cc')),with=F]
test3<-melt(test2,id=grepincols(test2,c('clust','pair')))
test4<-test3[!grepl('WT',value)]

singles<-phenDT[single_lv==T]
singles[,c("growth.rate","phenotypic.index","fracon.glu","fracon.gal","yfp.mean.glu","yfp.mean.gal"):=lapply(.SD,zscore),.SDcols=phenocols]
singles2<-melt(singles[,notin(colnames(singles),c('aa1','bb1','cc1')),with=F],id= notin(colnames(singles),c('aa','bb','cc')))
singles3<-singles2[!grepl('WT',value)]

DTm<-merge(test4,singles3,by=c('variable','value'),allow.cartesian=T)
mDT<-melt(DTm,id=notin(colnames(DTm),phenocols),variable.name='phenotype',value.name='obs')
mDT[,locus:=gsub('aa','GAL4',gsub('bb','GAL3',gsub('cc','GAL80',variable)))]
mDT[,fac:=applyPaste(data.frame(pairing,locus),': ')]
mDT[,dup:=applyPaste(mDT)]
mDT1<-mDT[!duplicated(dup)];mDT1[,dup:=NULL]
sample(letters,4)
pIJLN<-ggplot(mDT1,aes(fac,obs))+geom_jitter(alpha=0.2,width=0.2,shape=21,size=0.1,col='grey70')+
	geom_violin(aes(col=pairing),alpha=0.5)+coord_flip()+geom_hline(yintercept=0,col='grey50',alpha=0.5)+
	facet_grid(`sig epist cluster`~phenotype)

# w<-12.4;h<-6.7
# ggsave(paste0(figout,'1707xx-180823-pIJLN-PhenotypesOfEpistasisClusters-ScaledPhenotype.png'), pIJLN,width=w,height=h)


############
### show cluster of constitutive GAL4 with GAL80: phenotypic index and growth rate
############

DTm<-merge(phenDT,clusterframe,by=genos)
dSub<-DTm[((cc_clust=='uninducible'|cc=='GAL80.WT'|cc=='GAL80.delta')&(aa_clust=='constitutive'|aa=='GAL4.WT'|aa=='GAL4.delta')&bb=='GAL3.WT')|WT==T]
credmin<-dSub[aa=='GAL4.delta',mean(growth.rate)]
dSub[,GR:=growth.rate]
dSub[growth.rate<credmin,GR:=credmin]
G80s<-unique(dSub$cc)
G80lev<-G80s[c(2,3,4,5,6,1)]
G4s<-unique(dSub$aa)
G4lev<-G4s[c(7,8,1,6,5,2,3,4)]
dSub[,cc:=factor(cc,levels=G80lev[6:1])]
dSub[,aa:=factor(aa,levels=G4lev[8:1])]
dSub[,PI:=phenotypic.index]
sample(letters,4)
collim<-c(round(credmin,2),0.28)
pZOMT<-ggplot(dSub,aes(cc,aa))+geom_tile(aes(fill=GR))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),limits=collim)+
	theme(axis.text.x=element_text(angle=30,hjust=1))+ylab('GAL4 alleles')+xlab('GAL80 alleles')

collim<-c(0,2)
pMDBJ<-ggplot(dSub,aes(cc,aa))+geom_tile(aes(fill=phenotypic.index))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),limits=collim)+
	theme(axis.text.x=element_text(angle=30,hjust=1))+ylab('GAL4 alleles')+xlab('GAL80 alleles')

collim<-c(0,1)
mDT<-melt(dSub[,.(fracon.gal,fracon.glu,cc,aa)],id=c('cc','aa'))
pZUOM<-ggplot(mDT,aes(cc,aa))+geom_tile(aes(fill=value))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),limits=collim)+
	theme(axis.text.x=element_text(angle=30,hjust=1))+ylab('GAL4 alleles')+xlab('GAL80 alleles')+facet_grid(~variable)

collim<-c(0,1)
dSub[,GRnorm:=(GR-credmin)/max((GR-credmin),na.rm=T)]
mDT<-melt(dSub[,.(fracon.gal,fracon.glu,GRnorm,cc,aa)],id=c('cc','aa'))
pJXDV<-ggplot(mDT,aes(cc,aa))+geom_tile(aes(fill=value))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),limits=collim)+
	theme(axis.text.x=element_text(angle=30,hjust=1))+ylab('GAL4 alleles')+xlab('GAL80 alleles')+facet_grid(~variable)

# w<-7.7;h<-3.3
# ggsave(paste0(figout,'1707xx-180823-pJXDV-GoF-FraconGluGalAndGRnorm.png'), pJXDV,width=w,height=h)

dSub[,ctrl:=grepl('WT',aa)|grepl('delta',aa)|grepl('WT',cc)|grepl('delta',cc)]
# ggplot(dSub[order(ctrl,decreasing=T)],aes(fracon.glu,GR,col=ctrl))+geom_point(shape=21)

bycols<-c('clone.named',grepincols(dSub,'clust'),'aa','bb','cc','ctrl')
summ<-summary.func.all(dSub,c('growth.rate','GR','fracon.glu','yfp.mean.glu','yfp.mean.gal'),bycols)
DTm<-merge(dSub,summ,by=bycols)
summ[ctrl==F,genotype:='Combinations of GOF mutants']
summ[ctrl==T,genotype:='GOF mutants with control alleles']
summ[clone.named=='GAL3.WT GAL80.WT GAL4.WT',genotype:='WT']
labs<-c('GOF mutants with control alleles','Combinations of GOF mutants','WT')
summ[,genotype:=factor(genotype,levels=labs)]

sub<-dSub[ctrl==F]
mod1<-lm(fracon.glu~growth.rate,sub)
rsq1<-paste0('R2: ',round(summary(mod1)$r.squared,2))
sub<-dSub[ctrl==T&clone.named!='GAL3.WT GAL80.WT GAL4.WT']
mod2<-lm(fracon.glu~growth.rate,sub)
rsq2<-paste0('R2: ',round(summary(mod2)$r.squared,2))

sub<-dSub[ctrl==F]
mod1<-lm(yfp.mean.glu ~growth.rate,sub)
rsq1a<-paste0('R2: ',round(summary(mod1)$r.squared,2))
sub<-dSub[ctrl==T&clone.named!='GAL3.WT GAL80.WT GAL4.WT']
mod2<-lm(yfp.mean.glu ~growth.rate,sub)
rsq2a<-paste0('R2: ',round(summary(mod2)$r.squared,2))


#labs<-c('GOF mutants with control alleles','Combinations of GOF mutants')
txt<-data.table(genotype=factor(labs,levels=labs),label=c(rsq2,rsq1,NA),x=c(0.3,0.8,NA),y=c(0.22,0.05,NA))

ylimit<-c(0,0.31)
DTm<-merge(txt,summ,by='genotype')
pUEHX<-ggplot(data= DTm,aes(y=growth.rate_mean,x=fracon.glu_mean, ymax=growth.rate_upper,ymin=growth.rate_lower, xmax=fracon.glu_upper,xmin=fracon.glu_lower,col=genotype),size=3)+
#	geom_point(data=dSub[order(ctrl,decreasing=T)],aes(fracon.glu,GR,col=ctrl),shape=21)+
	geom_point()+geom_errorbar()+geom_errorbarh()+geom_smooth(method='lm')+
	scale_colour_manual(values=c('cornflowerblue','orange1','red'))+
	ylab(expression(paste(mu,' ',hr^-1)))+xlab('fraction ON in glucose')+ylim(ylimit)+
	geom_text(aes(x=x,y=y,label=label))
w<-6.7;h<-2.8
ggsave(paste0(figout,'1707xx-180827-pUEHX-FraconVsGrowthRate-GoFSubset.png'), pUEHX,width=w,height=h)

#labs<-c('GOF mutants with control alleles','Combinations of GOF mutants')
txt<-data.table(genotype=factor(labs,levels=labs),label=c(rsq2a,rsq1a,NA),x=c(8000,10000,NA),y=c(0.28,0.05,NA))

DTm<-merge(txt,summ,by='genotype')
pMFKB<-ggplot(data= DTm,aes(y=growth.rate_mean,x=yfp.mean.glu_mean, ymax=growth.rate_upper,ymin=growth.rate_lower, xmax=yfp.mean.glu_upper,xmin=yfp.mean.glu_lower,col=genotype),size=3)+
#	geom_point(data=dSub[order(ctrl,decreasing=T)],aes(fracon.glu,GR,col=ctrl),shape=21)+
	geom_point()+geom_errorbar()+geom_errorbarh()+geom_smooth(method='lm')+
	scale_colour_manual(values=c('cornflowerblue','orange1','red'))+
	ylab(expression(paste(mu,' ',hr^-1)))+xlab('YFP signal in glucose')+ylim(ylimit)+
	geom_text(aes(x=x,y=y,label=label))
w<-6.7;h<-2.8
ggsave(paste0(figout,'1707xx-180827-pMFKB-YFPMeanGluVsGrowthRate-GoFSubset.png'), pMFKB,width=w,height=h)

txt<-data.table(genotype=factor(labs,levels=labs),label=c(rsq2a,rsq1a,NA),x=c(8000,10000,NA),y=c(0.28,60000,NA))
DTm<-merge(txt,summ,by='genotype')


pDMJT<-ggplot(data= DTm,aes(y=yfp.mean.gal_mean,x=yfp.mean.glu_mean, ymax=yfp.mean.gal_upper,ymin=yfp.mean.gal_lower, xmax=yfp.mean.glu_upper,xmin=yfp.mean.glu_lower,col=genotype),size=3)+
#	geom_point(data=dSub[order(ctrl,decreasing=T)],aes(fracon.glu,GR,col=ctrl),shape=21)+
	geom_point()+geom_errorbar()+geom_errorbarh()+geom_smooth(method='lm')+
	scale_colour_manual(values=c('cornflowerblue','orange1','red'))+
	ylab('YFP signal in galactose')+xlab('YFP signal in glucose')+
	geom_text(aes(x=x,y=y,label=label))
w<-6.7;h<-2.8
ggsave(paste0(figout,'1707xx-180827-pDMJT-YFPMeanGluVsYFPmeanGal-GoFSubset.png'), pDMJT,width=w,height=h)


###################################################
### genotypic ensemble model
###################################################


############
### use archetypal expression clusters to predict growth rate 
############

boo<-c('aa','bb','cc')
clustcols<-unlist(lapply(1:3,function(i){
	clustcol<-paste0(boo[i],'_','clust')
}))
pheno<-'growth.rate'
genos<-c('clone.named','aa','bb','cc')

DTm<-merge(phenDT,clusterframe,by=genos)
#DTm<-merge_recurse3(qDT=phenDT,dt.list=dts,by.list=boo)
qDT<-DTm[,c(pheno,clustcols,genos,'mutant_gene'),with=F]
setnames(qDT,pheno,'x')
if(pheno=='growth.rate'){
	credMin<-qDT[aa=='GAL4.delta',mean(x,na.rm=T)-sd(x,na.rm=T)]
	qDT[x<credMin,x:=credMin]
}
pred1<-qDT[,summary.func(x,'pred'),by=c(clustcols,'mutant_gene')][order(pred_N,decreasing=T)]

pred<-qDT[,summary.func(x,'pred'),by=c(clustcols)][order(pred_N,decreasing=T)]
genosumm<-qDT[,summary.func(x,'obs'),by=genos]
DTm2<-merge(merge(qDT,pred,by=clustcols),genosumm,by=genos)
varexpA<-DTm2[,varexplained(obs_mean,pred_mean)]
varexpB<-DTm2[,varexplained(x,pred_mean)]
mod_genotypic_ensemble<-lm(obs_mean~pred_mean,DTm2) # this comes into play below
p0<-ggplot(DTm2[obs_CoV<0.5],aes(pred_mean,obs_mean,xmin=pred_lower,xmax=pred_upper,ymin=obs_lower,ymax=obs_upper))+geom_errorbar(col='grey70',alpha=0.2)+geom_errorbarh(col='grey70',alpha=0.2)+
	geom_point(alpha=0.1)+geom_abline()
p1<-ggplot(DTm2[obs_CoV<0.5],aes(pred_mean,obs_mean,xmin=pred_lower,xmax=pred_upper,ymin=obs_lower,ymax=obs_upper))+
	geom_point(alpha=0.1)+geom_abline()+annotate('text',label=paste0('R2: ',round(varexpA,2)),x=0.07,y=0.25)

summ<-DTm2[,lapply(.SD,unique),by='clone.named',.SDcols=c(grepincols(DTm2,c('obs','pred')),clustcols,notin(genos,'clone.named'))]
summ[,diff_mean:=obs_mean-pred_mean]
summ[,diff_se:=sqrt(obs_se^2+pred_se^2)]
summ[,diff_df:=obs_N+pred_N-2]
summ2<-data.table(summ,summ[,t.test3(obs_mean,pred_mean,se=diff_se,df=diff_df)])
summ2[,p.adj:=p.adjust(p.value,'fdr')]
summ2[,sig:=p.adj<0.01&abs(diff_mean)>0.05]
summ2[,`FDR < 0.01`:=p.adj<0.01]
# ggplot(summ2[obs_CoV<0.5],aes(pred_mean,obs_mean,xmin=pred_lower,xmax=pred_upper,ymin=obs_lower,ymax=obs_upper))+geom_errorbar(col='grey70',alpha=0.2)+geom_errorbarh(col='grey70',alpha=0.2)+
	# geom_point(alpha=0.1,aes(col=sig))+geom_abline()+scale_colour_manual(values=c('orange1','steelblue'))
p2<-ggplot(summ2[obs_CoV<0.5][order(sig,decreasing=T)],aes(pred_mean,obs_mean))+
	geom_point(alpha=0.1,aes(col=`FDR < 0.01`))+geom_abline()+scale_colour_manual(values=c('steelblue','orange1'))

jitfac<-0.0002
summ2[sig==T,x:=pred_mean+jitfac]
summ2[sig==F,x:=pred_mean-jitfac]
p3<-ggplot(summ2[obs_CoV<0.2][order(sig,decreasing=T)],aes(x,obs_mean))+
	geom_point(alpha=0.1,aes(col=`FDR < 0.01`))+geom_abline()+scale_colour_manual(values=c('steelblue','orange1'))

DT1<-CloneNamedSummary(summ2[,!c('aa','bb','cc'),with=F],allele=F)
p4<-ggplot(DT1[obs_CoV<0.5],aes(x,obs_mean))+
	geom_point(alpha=0.1,aes(col=`FDR < 0.01`))+geom_abline()+scale_colour_manual(values=c('steelblue','orange1'))+
	facet_grid(~mutant_gene)


############
### plot predictions of multiplicative model vs genotypic ensemble method
############

dTemp1<-merge(geom_model,clusterframe,by=c(boo,'clone.named'),all.x=T)
dTemp1[,pred_mean:=mean(pheno_mean,na.rm=T),by=c(grepincols(dTemp1,'clust'))]
dTemp<-CloneNamedSummary(dTemp1[,.(clone.named,pheno_mean,mult_pred_mean,pred_mean)],allele=F)
dTemp[,pairing:=gsub('\\ ',' v ',mutant_gene)]
mDT<-melt(na.exclude(dTemp[,.(pheno_mean,pairing, mult_pred_mean, pred_mean)]),id=c('pheno_mean','pairing'))
mDT[,model:=gsub('pred_mean','genotypic ensemble prediction',gsub('mult_pred_mean','multiplicative model',variable))]
sample(letters,4)
mDT[,varexp:=paste0('R2: ',round(varexplained(pheno_mean,value),2)),by=model]
mDT[,x:=0.2]
mDT[,y:=-0.05]
p1<-ggplot(mDT,aes(value,pheno_mean))+geom_point(shape=21,size=0.3,alpha=0.5)+geom_abline(col='red')+
	facet_grid(~model)+xlab(expression(paste('observed ',mu,' ',hr^-1)))+ylab(expression(paste('predicted ',mu,' ',hr^-1)))+
	geom_text(aes(label=varexp,x=x,y=y))

mDT[,varexp:=paste0('R2: ',round(varexplained(pheno_mean,value),2)),by=c('model','pairing')]
mDT[,x:=0.2]
mDT[,y:=0.05]
mDT[,keep:=grepl('v',pairing)]
mDT[nchar(pairing)>15,keep:=F]
pHMFA<-ggplot(mDT[keep==T],aes(value,pheno_mean))+geom_point(shape=21,size=0.3,alpha=0.5)+geom_abline(col='red')+
	facet_grid(model~pairing)+xlab(expression(paste('observed ',mu,' ',hr^-1)))+ylab(expression(paste('predicted ',mu,' ',hr^-1)))+
	geom_text(aes(label=varexp,x=x,y=y))

w<-8.8;h<-5.8
ggsave(paste0(figout,'1707xx-180827-pHMFA-CompareGenotEnsembleVsMultModel.png'), pHMFA,width=w,height=h)

loglikbyhand<-function(x,y){
	-length(y)/2*log(2*pi) - length(y)/2*log(var(y)) - 1/(2*var(y)) - 1/(2*var(y))*sum((y-x)^2)
}
mDT[, loglikbyhand(pheno_mean,value),by=model]
# doesn't really work

mod1<-lm(pheno_mean~value,mDT[model=='multiplicative model'])
mod2<-lm(pheno_mean~value,mDT[model=='genotypic ensemble prediction'])

l1<-logLik(mod1)
l2<-logLik(mod2)

par1<-phenDT[,length(unique(c(aa,bb,cc)))]
par2<-nrow(na.exclude(dTemp1[,mean(pheno_mean),by=c('aa_clust','bb_clust','cc_clust')]))

aic1<-as.numeric(par1*2-2*log(l1))
aic2<-as.numeric(par2*2-2*log(l2))
exp((aic2-aic1)/2)
exp((AIC(mod2,k=par2)-AIC(mod1,k=par1))/2)

################
### cluster epistasis across clone.named using hdbscan
################

clusterfun2<-function(DTr,left,right,singleGeno,PI,interp=T,min.pts=2){
	valuevar<-'value'
	if(interp==T){
		DTr[,temp:=DTr[,singleGeno,with=F]]
		mod1<-lm(value~temp+query,DTr)
		predmod<-predict(mod1,DTr)
		DTr[,dirty_pred:=predmod]
		DTr[,interp_val:=value]
		DTr[is.na(value),interp_val:=dirty_pred]
		valuevar<-'interp_val'
		varexp1<-DTr[,varexplained(value,dirty_pred)]
		DTr[,temp:=NULL]
		print(paste0('interpolated values var explained were ',round(varexp1,2)))
	}
	form<-as.formula(castform(left,right))
	cDT<-dcast(DTr,form,value.var=valuevar)
	sDT1<-cDT[,lapply(.SD,scale),.SDcols=c(notin(colnames(cDT),singleGeno))]
	sDT1[,x:=cDT[,singleGeno,with=F]]
	setnames(sDT1,'x',singleGeno)
	lv<-as.logical(paste0(unlist(sDT1[,(lapply(.SD,function(x)as.logical(prod(unique(!is.na(x))))))])))
	nacols<-colnames(sDT1)[!lv]
	notnacols<-colnames(sDT1)[lv]
	sDT<-sDT1[,notnacols,with=F]
	setcolorder(sDT, singleGeno)
	clusters<-hdbscan(sDT[,!singleGeno,with=F],minPts=min.pts)
	sDT[,arb:=1:nrow(sDT)]
	clust<-clusters$cluster
	setnames(PI,singleGeno,'x')
	ordframe <-data.table(arb=1:nrow(sDT),X1=1:nrow(sDT),X2=1:nrow(sDT),x=unlist(sDT[,singleGeno,with=F]),A=unlist(sDT[,singleGeno,with=F]),B=unlist(sDT[,singleGeno,with=F]),cluster1=clust)
	boo<-merge(ordframe,PI,by='x')[order(cluster1,V1)]
	boo[,ord:=1:nrow(boo)]
	distmat<-data.table(melt(as.matrix(cophenetic(clusters$hc))))
	dMT<-merge(merge(distmat,boo[,c('A','X1'),with=F],by='X1'), boo[,c('B','X2'),with=F],by='X2')
	dMT[,A:=factor(A,levels= boo$A)]
	dMT[,B:=factor(B,levels= boo$B)]
	p<-ggplot(dMT,aes(A,B))+geom_tile(aes(fill=value))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'))+theme(axis.text.x=element_text(angle=30,hjust=1))
	setnames(boo,c('x','V1'),c(singleGeno,'PI_mean'))
	boo_out<-boo[,c(singleGeno,'PI_mean','cluster1','ord'),with=F]
	return(list(plot=p,
		cluster_order=boo_out,
		clust_obj=clusters
		))
}

DTx<-copy(phenDT)
summx<-summary.func.all(DTx,c('fracon.gal','fracon.glu','yfp.mean.gal','yfp.mean.glu','growth.rate'),by='clone.named')
DT0x<-copy(DT1)#copy(phenDT)#
DT0<-merge(summx,DT0x,by='clone.named')
phenocols1<-'mag'#dens.cols#c('growth.rate')#
phenocols<-notin(grepincols(DT0,c('p.adj','_mean')),c('obs_mean','pred_mean'))#c('diff_mean','obs_mean')#dens.cols#
phenoexc<-'p.adj'#'diff_mean'
genocols<-c(colnames(CloneNamedSummary(DT00[,c('clone.named','allele.named'),with=F])))
subcols<-c(phenocols,'clone.named','aa','bb','cc')
a<-DT0[,c(subcols),with=F]
bycols2<-c('clone.named','aa','bb','cc')
genocols<-c(colnames(CloneNamedSummary(DT00[,c('clone.named','allele.named'),with=F])))
subcols<-c(phenocols,'clone.named','aa','bb','cc','single_lv','double_lv','sig')
a<-DT0[,c(subcols),with=F]
bycols2<-c('aa','bb','cc','clone.named','single_lv','double_lv','sig')
ord_within_cluster<-'obs_mean'
# ord_within_cluster<-'phenotypic.index'
summ1<-a
summ2<-summ1[,lapply(.SD,function(x){
	scale(x)
}),.SDcols=c(phenocols)]
lv<-!is.nan(as.numeric(paste0(summ2[1])))
keep<-colnames(summ2)[lv]
summ<-summ2[,keep,with=F]

# ord_within_cluster<-'phenotypic.index'


mDT1<-melt(a,id=c(bycols2))
varcols<-notin(phenocols,phenoexc)
#varcols<-phenocols
dens.arb<-paste0('X',1:length(varcols))
DT.dens.arb<-data.table(variable=varcols,dens.arb=factor(dens.arb,levels=dens.arb))
mDT<-merge(mDT1,DT.dens.arb,by='variable')

# get singles data across all genetic backgrounds together for GAL4
G4<-mDT[!(aa=='GAL4.WT')]
aabb<-G4[cc=='GAL80.WT',!'cc']
setnames(aabb,'bb','query')
aacc<-G4[bb=='GAL3.WT',!'bb']
setnames(aacc,'cc','query')
DTr<-rbind(aabb,aacc)
left<-c('aa')
right<-c('query','dens.arb')
singleGeno<-c('aa')
# PI<-DT0[,sqrt(sum(diff_mean^2,na.rm=T)),by=c(singleGeno)][order(V1)]
PI<-DT0[G4.single==T,c(singleGeno,'obs_mean'),with=F][order(obs_mean)]
PI[,V1:=obs_mean]
test<-clusterfun2(DTr,left,right,singleGeno,PI,interp=T,min.pts=2)
# test$plot;dev.new();plot(test$clust_obj);test$cluster_order

# get singles data across all genetic backgrounds together for GAL3
G3<-mDT[!(bb=='GAL3.WT')]
bbcc<-G3[aa=='GAL4.WT',!'aa']
setnames(bbcc,'cc','query')
bbaa<-G3[cc=='GAL80.WT',!'cc']
setnames(bbaa,'aa','query')
DTr<-rbind(bbcc, bbaa)
left<-c('bb')
right<-c('query','dens.arb')
singleGeno <-c('bb')
singleGeno2<-c('bb')
# PI<-DT0[,sqrt(sum(diff_mean^2,na.rm=T)),by=c(singleGeno)][order(V1)]
PI<-DT0[G3.single==T,c(singleGeno,'obs_mean'),with=F][order(obs_mean)]
PI[,V1:=obs_mean]
test2<-clusterfun2(DTr,left,right,singleGeno,PI,interp=T,min.pts=2)
# test2$plot;dev.new();plot(test2$clust_obj);test2$cluster_order

# get singles data across all genetic backgrounds together for GAL80
G80<-mDT[!(cc=='GAL80.WT')]
ccbb<-G80[aa=='GAL4.WT',!'aa']
setnames(ccbb,'bb','query')
ccaa<-G80[bb=='GAL3.WT',!'bb']
setnames(ccaa,'aa','query')
DTr<-rbind(ccbb, ccaa)
left<-c('cc')
right<-c('query','dens.arb')
singleGeno<-c('cc')
singleGeno2<-c('cc')
# PI<-DT0[,sqrt(sum(diff_mean^2,na.rm=T)),by=c(singleGeno)][order(V1)]
PI<-DT0[G80.single==T,c(singleGeno,'obs_mean'),with=F][order(obs_mean)]
PI[,V1:=obs_mean]
test3<-clusterfun2(DTr,left,right,singleGeno,PI,interp=T,min.pts=2)
# test3$plot;dev.new();plot(test3$clust_obj);test3$cluster_order

rankClust<-function(DT){
	DT[,cluster.orig:=cluster1]
	DT[,V1:=median(PI_mean),by='cluster1']
	ttt<-DT[,median(PI_mean),by='cluster1'][order(V1)]
	ttt[,cluster1:=1:nrow(ttt)]
	out<-merge(ttt,DT[,!'cluster1'],by='V1')
	out[,!'V1'][order(cluster1,PI_mean)]
}

G4DT<-rankClust(copy(test$cluster_order))
G3DT<-rankClust(copy(test2$cluster_order))
G80DT<-rankClust(copy(test3$cluster_order))

summ<-a

G80G3<-summ[aa=='GAL4.WT']
G80G3[,c('A','B'):=list(cc,bb)]
G3G4<-summ[cc=='GAL80.WT']
G3G4[,c('A','B'):=list(bb,aa)]
G80G4<-summ[bb=='GAL3.WT']
G80G4[,c('A','B'):=list(cc,aa)]

G80G3x<-data.table(A=c(c('GAL80.WT',G80DT$cc), c('GAL3.WT',G3DT$bb)))
G4G3x<-data.table(B=c(c('GAL4.WT',G4DT$aa), c('GAL3.WT',G3DT$bb)))

DTr<-rbind(G80G3,G3G4,G80G4)
DTr[,sig:=p.adj<0.05&abs(diff_mean)>0.035]
DTr[,A:=factor(A,levels= G80G3x$A[85:1])]
DTr[,B:=factor(B,levels= G4G3x$B[89:1])]
DTr[sig==T,`obs - pred`:=diff_mean]
DTr[sig==F,`obs - pred`:=0]
#DTr[abs(diff_mean)>0.13,`obs - pred`:=NA]
absmax<-max(abs(range(DTr$`obs - pred`,na.rm=T)))
collim<-c(-absmax,absmax)
center0colors<-c('blue','cornflowerblue','white','indianred','red')
pRTBU<-ggplot(DTr,aes(B, A))+geom_tile(aes(fill=`obs - pred`))+scale_fill_gradientn(colours=center0colors,na.value='grey70',limits=collim)+
	theme(axis.text.x=element_text(angle=50,hjust=1,size=6))+
	theme(axis.text.y=element_text(size=6))+
	scale_y_discrete(position='right')+xlab('')+ylab('')+
	theme(plot.margin = unit(c(1,1,1,1), "cm"))

tricolor<-c('cornflowerblue','yellow','indianred')
pRZWG<-ggplot(DTr,aes(B, A))+geom_tile(aes(fill=growth.rate_mean))+scale_fill_gradientn(colours= tricolor,na.value='grey70')+
	theme(axis.text.x=element_text(angle=50,hjust=1,size=6))+
	theme(axis.text.y=element_text(size=6))+
	scale_y_discrete(position='right')+xlab('')+ylab('')+
	theme(plot.margin = unit(c(1,1,1,1), "cm"))
bicolor<-c('white','indianred')
pRZWGb<-ggplot(DTr,aes(B, A))+geom_tile(aes(fill=growth.rate_mean))+scale_fill_gradientn(colours= bicolor,na.value='grey70')+
	theme(axis.text.x=element_text(angle=50,hjust=1,size=6))+
	theme(axis.text.y=element_text(size=6))+
	scale_y_discrete(position='right')+xlab('')+ylab('')+
	theme(plot.margin = unit(c(1,1,1,1), "cm"))

w<-13;h<-12.5
ggsave(paste0(figout,'1707xx-180710-pRTBU-ClusteredBy5PhenotypesAcrossGeneticBackgrounds-GrowthRateDeviationFromGenotypicEnsemble.png'), pRTBU,width=w,height=h)
w<-13;h<-12.5
ggsave(paste0(figout,'1707xx-180710-pRZWG-ClusteredBy5PhenotypesAcrossGeneticBackgrounds-GrowthRate.png'), pRZWG,width=w,height=h)
w<-13;h<-12.5
ggsave(paste0(figout,'1707xx-180710-pRZWGb-ClusteredBy5PhenotypesAcrossGeneticBackgrounds-WhiteToRedColors-GrowthRate.png'), pRZWGb,width=w,height=h)


####################################
#### linear regression using log-transformed values
####################################

phenoplot2<-function(summ,phenoy='GR',phenox='PI',hlight=NULL,plot=T,bycols=genocols,abline=F,facet=T,singles=T){
	mDT1<-melt(summ,notcols(summ,' v '),variable.name='pairing',value.name='boo')
	mDT1[pairing=='GAL3 v GAL80',test:=apply(data.frame(G3_mut_lv,G80_mut_lv,G4_mut_lv,WT_lv),1,function(x)as.logical(((x[1]|x[2])&!x[3])|x[4]))]
	mDT1[pairing=='GAL4 v GAL80',test:=apply(data.frame(G3_mut_lv,G80_mut_lv,G4_mut_lv, WT_lv),1,function(x)as.logical((!x[1]&(x[2]|x[3])|x[4])))]
	mDT1[pairing=='GAL3 v GAL4',test:=apply(data.frame(G3_mut_lv,G80_mut_lv,G4_mut_lv, WT_lv),1,function(x)as.logical(((x[1]|x[3])&!x[2])|x[4]))]
	mDT1[,highlight:=NA]
	if(!is.null(hlight))mDT1[,highlight:=unlist(grepl(hlight,clone.named))]
	mDT<-mDT1[test==T]
	oldcols<-grepincols(mDT,c(phenox,phenoy))
	newcols<-gsub('_mean','',gsub(phenox,'x',gsub(phenoy,'y', grepincols(mDT, oldcols))))
	setnames(mDT,oldcols,newcols)
	p0<-ggplot(mDT,aes(x,y))+geom_point(col='grey70',alpha=0.3)+
		xlab(phenox)+ylab(phenoy)
	if(singles==T){p<-p+
		geom_point(data=mDT[single_lv==T],aes(x,y,col= single_locus),shape=21,stroke=1,size=2)+
		geom_errorbar(data= mDT[single_lv==T],aes(ymin=y_lower,ymax=y_upper, col= single_locus),width=0)+
		geom_errorbarh(data= mDT[single_lv==T],aes(xmin=x_lower,xmax=x_upper, col= single_locus),height=0)+
		geom_point(data= mDT[WT_lv==T],aes(x,y),col='red',size=3,shape=21,stroke=1)+
		geom_errorbar(data= mDT[WT_lv==T],aes(ymin=y_lower,ymax=y_upper),col='red',width=0)+
		geom_errorbarh(data= mDT[WT_lv==T],aes(xmin=x_lower,xmax=x_upper),col='red',height=0)}

	if(facet==F)p<-p0
	if(facet==T)p<-p0+facet_grid(~ pairing)+theme_light()
	if(abline==T)p<-p+geom_abline(col='grey70')
	if(!is.null(hlight)){
		p<-p+geom_point(data= mDT[highlight==T],aes(x,y),col='orange1',size=1,alpha=0.6)+
			geom_point(data= mDT[highlight==T&single_lv==T],aes(x,y),col='black',size=4)+
			geom_point(data= mDT[highlight==T&single_lv==T],aes(x,y),col='orange1',size=3)
		}
	
	if(plot==F){return(mDT)}
	if(plot==T){return(p)}
	#geom_point(data=mDT[single_lv==T&test==T],aes(xmin=PI_lower,xmax=PI_upper,ymax=GR_upper,ymin=GR_lower))+geom_errorbar(col='red',alpha=0.3)+geom_errorbarh(col='red',alpha=0.3)+geom_point()
}

DT0<-copy(phenDT)
backg<-mean(DT0[grepl('GAL4.delta',clone.named)]$growth.rate,na.rm=T)
backsd<-sd(DT0[grepl('GAL4.delta',clone.named)]$growth.rate,na.rm=T)
pheno<-'growth.rate'
genos<-c('aa','bb','cc')
dsub<-DT0[,c(pheno,genos),with=F]
wt<-'GAL4.WT'
dsub[,aa:=factor(aa,levels=c(wt,notin(unique(aa),wt)))]
wt<-'GAL3.WT'
dsub[,bb:=factor(bb,levels=c(wt,notin(unique(bb),wt)))]
wt<-'GAL80.WT'
dsub[,cc:=factor(cc,levels=c(wt,notin(unique(cc),wt)))]


setnames(dsub,pheno,'x')
lDT<-logconv2(dsub,backg=backg,backgsd=backsd,genos=genos)
mod1<-lm(lx~aa+bb+cc, lDT)
predDT<-expMod(dsub,mod1,backg=backg,genos=genos)
DTm<-merge(dsub,predDT,by=genos)
varexplained(DTm$x,DTm$fit_mean)
AIC(mod_genotypic_ensemble)
AIC(mod1)
exp((AIC(mod_genotypic_ensemble)-AIC(mod1))/2)
exp((BIC(mod_genotypic_ensemble)-BIC(mod1))/2)


DTm[,clone.named:=applyPaste(data.frame(bb,cc,aa),' ')]
DTm2<-merge(CloneNamedSummary(DTm[,!genos,with=F],allele=F), DTm[,summary.func(x,'GR'),by='clone.named'],by='clone.named')

mDT1<-melt(DTm2,notcols(DTm2,'v G'),variable.name='pairing',value.name='boo')
mDT1[pairing=='G3vG80_lv',test:=apply(data.frame(G3_mut_lv,G80_mut_lv,G4_mut_lv,WT_lv),1,function(x)as.logical(((x[1]|x[2])&!x[3])|x[4]))]
mDT1[pairing=='G4vG80_lv',test:=apply(data.frame(G3_mut_lv,G80_mut_lv,G4_mut_lv, WT_lv),1,function(x)as.logical((!x[1]&(x[2]|x[3])|x[4])))]
mDT1[pairing=='G3vG4_lv',test:=apply(data.frame(G3_mut_lv,G80_mut_lv,G4_mut_lv, WT_lv),1,function(x)as.logical(((x[1]|x[3])&!x[2])|x[4]))]
mDT1[,highlight:=NA]
hlight<-NULL
if(!is.null(hlight))mDT1[,highlight:=unlist(grepl(hlight,clone.named))]

p<-ggplot(mDT1,aes(fit_mean,GR_mean))+geom_point()+facet_grid(~pairing)+geom_abline(col='red')


DT0<-copy(phenDT)
phenocols<-c('growth.rate','phenotypic.index','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')
phenosh<-c('GR','PI','fONglu','fONgal','ymglu','ymgal')
genocols<-colnames(CloneNamedSummary(DT00[,c('clone.named','allele.named'),with=F]))
setnames(DT0,phenocols,phenosh)

summ1<-summary.func.all(DT0,phenosh,c('clone.named','allele.named'))
summ2<-CloneNamedSummary(summ1)
summ<-merge(summ2,predDT,by=c('aa','bb','cc'))

phenoplot2(summ[GR_CoV<0.5],phenoy='GR',phenox='fit',bycols=genocols,abline=T,facet=T,singles=F)













