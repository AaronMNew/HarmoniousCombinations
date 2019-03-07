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
writeplot<-F
br<-seq(-0.5,2.7,length=61)
library(tidyr)
factorDefault<-c('orange1','cornflowerblue','darkgreen','red','chartreuse3')



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

phenocols<-c('growth.rate','phenotypic.index','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')
phenosh<-c('GR','PI','fONglu','fONgal','ymglu','ymgal')

dens.cols<-grepincols(DT00,'dens')
phenDT<-CloneNamedSummary(DT00[,c('clone.named','allele.named',phenocols),with=F])
densDT<-CloneNamedSummary(DT00[,c('clone.named','allele.named',dens.cols),with=F])
genocols<-colnames(CloneNamedSummary(DT00[,c('clone.named','allele.named'),with=F]))

###################################################
### consider fold induction
###################################################
phenDT[,fold.induction:=yfp.mean.gal/yfp.mean.glu]
phenDT[,summary.func(fold.induction),by=clone.named]
phenDT[clone.named==c('GAL3.WT GAL80.delta GAL4.WT'),]
summ<-summary.func.all(phenDT,c('fracon.gal','fracon.glu','fold.induction'),by='clone.named')
ggplot(summ[order(fold.induction_mean)],aes(fracon.glu_mean,fracon.gal_mean,col=log2(fold.induction_mean)))+geom_point()+scale_color_gradientn(colours=c('cornflowerblue','yellow','indianred'))
library(plotly)
p <- plot_ly(summ, x = ~ fracon.glu_mean, y = ~ fracon.gal_mean, z = ~ log2(fold.induction_mean), color = ~ log2(fold.induction_mean)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'fracon glu'),
                     yaxis = list(title = 'fracon gal'),
                     zaxis = list(title = 'fold induction')))

###################################################
### cluster across glu,gal expression density with clone.named rows using hdbscan
###################################################
DT0<-copy(densDT)#copy(phenDT)#
phenocols<-dens.cols#c('growth.rate')#
genocols<-c(colnames(CloneNamedSummary(DT00[,c('clone.named','allele.named'),with=F])))
subcols<-c(phenocols,'clone.named','aa','bb','cc')
a<-DT0[,c(subcols),with=F]
bycols2<-c('clone.named','aa','bb','cc')
# ord_within_cluster<-'growth.rate'
ord_within_cluster<-'phenotypic.index'
summ1<-na.exclude(a[,lapply(.SD,function(x)mean(log10(x+0.001),na.rm=T)),by=bycols2])
summ2<-summ1[,lapply(.SD,function(x){
	scale(x)
}),.SDcols=c(phenocols)]
lv<-!is.nan(as.numeric(paste0(summ2[1])))
keep<-colnames(summ2)[lv]
summ<-summ2[,keep,with=F]
clusters<-hclust(dist(summ))


PI4<-phenDT[G4.single==T,mean(phenotypic.index,na.rm=T),by='aa']
PI3<-phenDT[G3.single==T,mean(phenotypic.index,na.rm=T),by='bb']
PI80<-phenDT[G80.single==T,mean(phenotypic.index,na.rm=T),by='cc']

dt.list=list(PI4,PI3,PI80)
boo1<-c('aa','bb','cc')
lapply(1:length(dt.list),function(i){
	x<-dt.list[[i]]
	setnames(x,c('V1'),c(paste0(boo1[i],'_','PI_mean')))
})

DTm1<-merge_recurse3(qDT=summ1,dt.list=dt.list,by.list=boo1)
DTm1[,cluster:=cutree(clusters,5)]
PI<-phenDT[,mean(phenotypic.index,na.rm=T),by='clone.named']
DTm11<-merge(DTm1,PI,by='clone.named')
hdbtest<-hdbscan((grepcols(DTm1,'dens')),minPts=20)
DTm1[,cluster:=hdbtest$cluster]
DTm1[,cluster_outlier_score:=hdbtest$outlier_score]

PI<-phenDT[,mean(phenotypic.index,na.rm=T),by='clone.named']
DTm11<-merge(DTm1,PI,by='clone.named')

rankClust<-function(DT){
	DT[,cluster.orig:=cluster]
	DT[,V2:=median(V1),by='cluster']
	ttt<-DT[,median(V1),by='cluster'][order(V1)]
	ttt[,cluster:=1:nrow(ttt)]
	setnames(ttt,"V1",'V2')
	out<-merge(ttt,DT[,!'cluster'],by='V2')
	out[,!c('V2','cluster.orig'),with=F][order(cluster,V1)]
}

# DTm<-rankClust(DTm11)
# mDT1<-melt(DTm[,!c('cluster.orig','V2'),with=F],id=c(bycols2,grepincols(DTm,'PI'),'cluster','V1'))[order(cluster,V1)]


DTm<-DTm11
lev<-DTm[,mean(V1),by=c('clone.named','cluster')][order(cluster,V1)]
mDT1<-melt(DTm,id=c(bycols2,grepincols(DTm,'PI'),'cluster','V1','cluster_outlier_score'))[order(cluster,V1)]
denscolconv(mDT1,vnull=F)
mDT<-CloneNamedSummary(mDT1[,!c('aa','bb','cc'),with=F],allele=F)
mDT[,Y:=factor(clone.named,levels=lev$clone.named)]
mDT[,X:=factor(round(AU.FL,3))]
mDT[sugar=='gal',X:=factor(round(AU.FL,3)+4)]
mDT[,dumm:=125]
mDT[,single:=1:.N,by=c('AU.FL','sugar')]
mDT[single_lv==F,single:=NA]
mDT[,`density, log10`:=value]
mDT[,`single mutant\nlocus`:=single_locus]
uniqfun<-function(x,fac=1,shift=1)(unique(x)[order(unique(x))][seq(shift,length(unique(x)),by=fac)])
collim<-range(mDT$value)
singM<-mDT[single_lv==T,.N,by=c('Y','dumm','single','single mutant\nlocus')]
xlab<-c(uniqfun(mDT$AU.FL,fac=4,shift=3),uniqfun(round(mDT$AU.FL,3),shift=3,fac=4))
xbreak<-uniqfun(mDT$X,fac=4,shift=3)
sample(letters,4)
pSZJE<-ggplot(mDT,aes(X,Y))+geom_tile(aes(fill=`density, log10`))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),na.value='white',limits=collim)+
	theme(axis.text.x=element_text(angle=90,vjust=0.4,hjust=1,size=8))+
	theme(axis.text.y=element_text(size=6))+
	geom_vline(xintercept=61,col='black',size=2)+
	geom_vline(xintercept=61,col='red',size=1.5)+
	geom_point(data=data.frame(x=126,y=1),aes(x,y),col='white')+
	geom_point(data=singM,aes(dumm-runif(length(dumm))*2,Y,colour=`single mutant\nlocus`),size=1,shape=21)+
	scale_colour_manual(values=c('orange1','indianred','chartreuse3','blue'),na.value='transparent')+
	scale_x_discrete(breaks=xbreak,labels=xlab)+
	theme(axis.title.y=element_blank(),axis.text.y=element_blank(),	axis.ticks.y=element_blank())+
	xlab('Expression in glucose               Exression in galactose   \n pseudo log10 A.U. Fluorescence units')+
	theme(axis.title=element_text(size=10,face="plain"))


if(writeplot==T){
	w<-5;h<-8
	ggsave(paste0(figout,'1707xx-180710-pSZJE-ExpressionDensityTilePlot-AllCloneNamed.png'), pSZJE,width=w,height=h)	
}




############
### plot archetypal expression densities: use all clusters
############

mDT1<-copy(mDT)
mDT1[,density:=exp(value)-0.001]
setnames(mDT1,'V1','PhenotypeForRanking')
bycols<-c('cluster','sugar','AU.FL')
summ1<-mDT1[,mean(density),by=bycols]
summ2<-mDT[,.N,by=bycols]
summ2[,prop:=N/length(unique(mDT$clone.named))]
summ<-merge(summ1,summ2,by=bycols)
summ[,facet:=percent(round(prop,3))]
summ[cluster%in%c(0),cluster_new:='outlying']
summ[cluster%in%c(1),cluster_new:='weakly inducible']
summ[cluster%in%c(2),cluster_new:='const, low expr 1']
summ[cluster%in%c(3),cluster_new:='const, low expr 2']
summ[cluster%in%c(4),cluster_new:='const, low expr 3']
summ[cluster%in%c(5),cluster_new:='uninducible']
summ[cluster%in%c(6),cluster_new:='inducible']
summ[cluster%in%c(7),cluster_new:='constitutive']


dumm<-summ[,unique(cluster_new),by='cluster']
summ[,cluster_new:=factor(cluster_new,levels=dumm$V1[c(5,6,7,1,2,3,4,0)+1])]
mDTm<-merge(summ,mDT1,by=bycols)
mDTm[,V2:=mean(density),by=c(bycols,'mutant_gene')]
mDTm[,dummyline:=paste0('z',sugar)]
mDTm[,mutant_facet:=mutant_gene]
mDTm[mutant_gene=='GAL3 GAL80 GAL4',mutant_facet:='Control']
mDTm[single_lv==T&aa=='GAL4.delta',mutant_facet:='Control']
mDTm[single_lv==T&bb=='GAL3.delta',mutant_facet:='Control']
mDTm[single_lv==T&cc=='GAL80.delta',mutant_facet:='Control']
dumm2a<-mDTm[,unique(mutant_facet),by='mutcount'][order(mutcount,V1)]
dumm2<-dumm2a[!(mutcount==1&V1=='Control')]
mDTm[,mutant_facet:=factor(mutant_facet,levels=dumm2$V1)]
labelframe<-mDTm[,.N/length(dens.cols),by=c('cluster_new','facet', 'mutant_facet')]
labelframe[,x:=2]
labelframe[,y:=0.75]
sample(LETTERS,4)
pBYGF<-ggplot(mDTm,aes(AU.FL,V2))+geom_point(data=mDTm,aes(AU.FL,density,col=sugar),size=0.1,alpha=0.3,shape=21)+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(cluster_new+facet~ mutant_facet)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density of observations') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))+
	geom_text(data=labelframe,aes(x=x,y=y,label=V1),size=3.5,col='grey50')
if(writeplot==T){
	w<-14;h<-10
	ggsave(paste0(figout,'1707xx-180720-pBYGF-ExpressionDensityClusters-LinePlot-FacetByMutatedGene.png'), pBYGF,width=w,height=h)
}
C1<-mDTm[,list(HDB1=unique(cluster_new)),by=c('aa','bb','cc','clone.named')]

############
### plot archetypal expression densities: combine all const low expr 
############

mDT1<-copy(mDT)
mDT1[,density:=exp(value)-0.001]
setnames(mDT1,'V1','PhenotypeForRanking')
mDT1[cluster%in%c(0),cluster_new:='outlying']
mDT1[cluster%in%c(1),cluster_new:='weakly inducible']
mDT1[cluster%in%c(2,3,4),cluster_new:='const, low expr']
mDT1[cluster%in%c(5),cluster_new:='uninducible']
mDT1[cluster%in%c(6),cluster_new:='inducible']
mDT1[cluster%in%c(7),cluster_new:='constitutive']

bycols<-c('cluster_new','sugar','AU.FL')
summ1<-mDT1[,mean(density),by=bycols]
summ2<-mDT1[,.N,by=bycols]
summ2[,prop:=N/length(unique(mDT$clone.named))]
summ<-merge(summ1,summ2,by=bycols)
summ[,facet:=percent(round(prop,3))]
dumm<-data.frame(arb=1:6,V1=summ[,unique(cluster_new)])
summ[,cluster_new:=factor(cluster_new,levels=dumm$V1[c(5,3,2,6,1,4)])]
mDTm<-merge(summ,mDT1,by=bycols)
mDTm[,V2:=mean(density),by=c(bycols,'mutant_gene')]
mDTm[,dummyline:=paste0('z',sugar)]
mDTm[,mutant_facet:=mutant_gene]
mDTm[mutant_gene=='GAL3 GAL80 GAL4',mutant_facet:='Control']
mDTm[single_lv==T&aa=='GAL4.delta',mutant_facet:='Control']
mDTm[single_lv==T&bb=='GAL3.delta',mutant_facet:='Control']
mDTm[single_lv==T&cc=='GAL80.delta',mutant_facet:='Control']
dumm2a<-mDTm[,unique(mutant_facet),by='mutcount'][order(mutcount,V1)]
dumm2<-dumm2a[!(mutcount==1&V1=='Control')]
mDTm[,mutant_facet:=factor(mutant_facet,levels=dumm2$V1)]
newlevels<-dumm$V1[c(5,3,2,6,1,4)]
mDTm[,cluster_new:=factor(cluster_new,levels=newlevels)]
labelframe<-mDTm[,.N/length(dens.cols),by=c('cluster_new','facet', 'mutant_facet')]
labelframe[,x:=2]
labelframe[,y:=0.75]
sample(LETTERS,4)

pUNWE<-ggplot(mDTm,aes(AU.FL,V2))+geom_point(data=mDTm,aes(AU.FL,density,col=sugar),size=0.1,alpha=0.3,shape=21)+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(cluster_new+facet~ mutant_facet)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density of observations') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))+
	geom_text(data=labelframe,aes(x=x,y=y,label=V1),size=3.5,col='grey50')

if(writeplot==T){
	w<-8;h<-6.5
	ggsave(paste0(figout,'1707xx-180720-pUNWE-ExpressionDensityClusters-AggregatedConstitutiveOutliers-LinePlot-FacetByMutatedGene.png'), pUNWE,width=w,height=h)
}

C2<-mDTm[,list(HDB2=unique(cluster_new)),by=c('aa','bb','cc','clone.named')]

############
### write clusters to new file for use downstream
############
clusterframe<-merge(C1,C2,by=c('aa','bb','cc','clone.named'))

DT0<-merge(phenDT,clusterframe,by=c('aa','bb','cc','clone.named'))
singles<-DT0[single_lv==T,list(cluster_new=unique(HDB2)),by=c('aa','bb','cc')]
G4<-singles[bb=='GAL3.WT'&cc=='GAL80.WT']
G3<-singles[aa=='GAL4.WT'&cc=='GAL80.WT']
G80<-singles[bb=='GAL3.WT'&aa=='GAL4.WT']
dt.list<-list(G4,G3,G80)
boo<-c('aa','bb','cc')
dts<-lapply(1:3,function(i){
	clustcol<-paste0(boo[i],'_','clust')
	setnames(dt.list[[i]],'cluster_new',clustcol)
	dt.list[[i]][,c(boo[i],clustcol),with=F]
})

clusterframe<-merge_recurse3(clusterframe,dts,by.list=boo)
clusterframe[clone.named=='GAL3.WT GAL80.07 GAL4.WT',c('HDB1','HDB2','cc_clust','cluster_new'):=list('inducible','inducible','inducible','inducible')]

if(writeplot==T){
	save(clusterframe,file=paste0(outputdata,'1707xx-180704-HDB_Clusters_Based_On_GAL_Distributions.rData'))
}


###################################################
### analysis of clusters based on pathway distributions
###################################################
rm(clusterframe)
datloc<-paste0(outputdata,'1707xx-180704-HDB_Clusters_Based_On_GAL_Distributions.rData')
load(datloc)

################
### show basic clusters
################

denscolconv<-function(DT, vnull =T){
	DT[,variable:=gsub('dens__','',variable)]
	DT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
	DT[,sugar:=paste0('g',sugar)]
	if(vnull ==T)DT[,variable:=NULL]
}


DT1<-merge(densDT,clusterframe,by=c('aa','bb','cc','clone.named'))
setnames(DT1,'HDB2','cluster_new')
mDT<-melt(DT1[,c('cluster_new',dens.cols),with=F],id='cluster_new')
denscolconv(mDT)
mDT[,dummyline:=paste0('z',sugar)]
mDT[,pheno:='fluorescence']

ylimits<-c(0,0.3)
ybreaks<-c(0.0,0.1,0.2,0.3)
ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
mDT1<-mDT[,mean(value,na.rm=T),by=c('cluster_new','dummyline','AU.FL','sugar')];setnames(mDT1,'V1','density_mean')
mDT1[density_mean>=max(ylimits),density_mean:=max(ylimits)][order(sugar,AU.FL,decreasing=T)]

pPXTE <-ggplot(mDT1,aes(AU.FL,density_mean))+theme_minimal()+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(~ cluster_new)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
	theme(legend.title=element_blank())+
	theme(strip.text.x = element_text(size = 10))

if(writeplot==T){
	w=7.8;h=1.6
	ggsave(paste0(figout,'1707xx-180724-pPXTE-CorePhenotypes-ExpressionDistributionLinePlots.png'), pPXTE,width=w,height=h)
}

################
### look at outlying clusters
################
genos<-c('clone.named','aa','bb','cc')
DT1<-merge(densDT,clusterframe,by=genos)
setnames(DT1,'HDB2','cluster_new')
mDT<-melt(DT1[,c('cluster_new',genos,dens.cols),with=F],id=c(genos,'cluster_new'))
denscolconv(mDT)
mDT[,dummyline:=paste0('z',sugar)]
mDT[,pheno:='fluorescence']

mDT0<-mDT[cluster_new=='outlying']
ylimits<-c(0,0.3)
ybreaks<-c(0.0,0.1,0.2,0.3)
ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
mDT1<-mDT0[,mean(value,na.rm=T),by=c('cluster_new',genos,'dummyline','AU.FL','sugar')];setnames(mDT1,'V1','density_mean')
mDT1[density_mean>=max(ylimits),density_mean:=max(ylimits)][order(sugar,AU.FL,decreasing=T)]


pOUTLYINGCLUSTERS<-ggplot(mDT1,aes(AU.FL,density_mean))+theme_minimal()+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_wrap(~clone.named,ncol=10)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
	theme(legend.title=element_blank())+
	theme(strip.text.x = element_text(size = 6))
if(writeplot==T){
	w<-15;h<-40
	ggsave(paste0(figout,'1707xx-180724-pOUTLYINGCLUSTERS-allOutlying.png'), pOUTLYINGCLUSTERS,width=w,height=h)

}

scored<-fread(paste0(layout.dir,'/1707xx-181025-Scoring_Outlying_By_Hand.txt'))
clones<-unique(mDT1$clone.named); clones<-clones[order(clones)]
clones[grepl('GAL3.13',clones)]

mScore<-na.exclude(melt(scored,id='asdf',value.name='cluster_by_hand')[order(asdf,variable)])
mScore[,clone.named:=clones]
mScore[,cluster_new:=cluster_by_hand]

mDT2<-merge(mScore[,!'cluster_new'],mDT1,by='clone.named')
mSplit<-split(mDT2,by='cluster_by_hand')

# check it's OK
fac_orig<-sqrt(prod(dim(mDT2)))/40
if(writeplot==T){
	lapply(mSplit,function(DTx){
		p<-ggplot(DTx,aes(AU.FL,density_mean))+theme_minimal()+
			geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
			facet_wrap(~clone.named,ncol=10)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
			ylab('density') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
			scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
			theme(legend.title=element_blank())+
			theme(strip.text.x = element_text(size = 6))
	
		w<-15;h<-h<-sqrt(prod(dim(DTx)))/fac_orig
		ggsave(paste0(figout,'1707xx-180724-', unique(DTx$cluster_by_hand),'-allOutlying.png'), p,width=w,height=h)
	})
}


mScore[cluster_by_hand=='partially inducible',cluster_new:='inducible']
clusterframe[,cluster_new:=HDB2]
tmp<-clusterframe[HDB2!='outlying']
tmp[,cluster_by_hand:=NA]
clusterframe2<-data.table(rbind(tmp,merge(clusterframe[HDB2=='outlying',!'cluster_new'],mScore[,c('clone.named','cluster_by_hand','cluster_new'),with=F],by='clone.named')))
clusterframe2[cluster_new =='partially inducible and partially constitutive'&aa_clust!='weakly inducible',cluster_new:='inducible']
DTm<-merge(clusterframe2,summary.func.all(phenDT,c('fracon.gal','fracon.glu','yfp.mean.gal','yfp.mean.glu'),'clone.named'),by='clone.named')
DTm[,fold.induction:=yfp.mean.gal_mean/yfp.mean.glu_mean]
DTm[cluster_new=='inducible'&(fracon.gal_mean<0.25|fracon.glu_mean>0.25)]
DTm[cluster_new=='inducible'&(fracon.gal_mean<0.25|fracon.glu_mean>0.25&fold.induction<9),cluster_new:='uninducible']
DTm[cluster_new=='trimodal in gal',cluster_new:='inducible']
consstclones<-unique(c("GAL3.24 GAL80.WT GAL4-L868E","GAL3.WT GAL80S-0 GAL4-L868E","GAL3.24 GAL80.WT GAL4-L868K","GAL3.WT GAL80.34 GAL4-L868K","GAL3.WT GAL80S-0 GAL4-L868K","GAL3.WT GAL80.34 GAL4-L868P","GAL3.WT GAL80S-0 GAL4-L868P","GAL3.WT GAL80S-2 GAL4-L868P","GAL3.WT GAL80.03 GAL4-L868F","GAL3.WT GAL80.43 GAL4-L868E","GAL3.WT GAL80.43 GAL4-L868K","GAL3.WT GAL80.43 GAL4-L868P","GAL3.24 GAL80.WT GAL4-L868E","GAL3.24 GAL80.WT GAL4-L868K","GAL3.24 GAL80.WT GAL4-L868P","GAL3.WT GAL80.34 GAL4-L868K","GAL3.WT GAL80.34 GAL4-L868P","GAL3.WT GAL80S-0 GAL4-L868E","GAL3.WT GAL80S-0 GAL4-L868K","GAL3.WT GAL80S-0 GAL4-L868P","GAL3.WT GAL80S-2 GAL4-L868K","GAL3.WT GAL80S-2 GAL4-L868P","GAL3.WT GAL80.WT GAL4-L868E","GAL3.WT GAL80.WT GAL4-L868K","GAL3.WT GAL80.WT GAL4-L868P","GAL3.WT GAL80.01 GAL4.WT","GAL3.WT GAL80.02 GAL4.WT","GAL3.WT GAL80.03 GAL4.WT","GAL3.WT GAL80.04 GAL4.WT","GAL3.WT GAL80.07 GAL4.WT","GAL3.WT GAL80.23 GAL4.WT","GAL3.WT GAL80.25 GAL4.WT","GAL3.WT GAL80.27 GAL4.WT","GAL3.WT GAL80.28 GAL4.WT","GAL3.WT GAL80.40 GAL4.WT","GAL3.WT GAL80.delta GAL4.WT","GAL3.WT GAL80S-2 GAL4-L868K",'GAL3.24 GAL80.WT GAL4-L868P','GAL3.33 GAL80.25 GAL4.WT'))

DTm[clone.named%in%clones,cluster_new:='constitutive']
uninclones<-c('GAL3.12 GAL80.34 GAL4.WT',"GAL3.13 GAL80.WT GAL4.29","GAL3.06 GAL80.34 GAL4.WT","GAL3.08 GAL80.15 GAL4.WT","GAL3.09 GAL80.10 GAL4.WT","GAL3.09 GAL80.15 GAL4.WT","GAL3.09 GAL80.22 GAL4.WT","GAL3.20 GAL80.33 GAL4.WT")
DTm[clone.named%in% uninclones,cluster_new:='uninducible']

clusterframe3<-DTm[,colnames(DTm)%in%colnames(clusterframe2),with=F]


################
### show how clusters lead to different fracon glu / gal bits
################
genos<-c('clone.named','aa','bb','cc')
phenocols<-c('growth.rate','phenotypic.index','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')

Aphens<-c('constitutive','inducible','uninducible')
Bphens<-c('constitutive','inducible','uninducible')
ABphens<-c('constitutive','inducible','uninducible')

clustset<-c('inducible','constitutive','uninducible')
abset<-c('inducible.constitutive','constitutive.uninducible','inducible.uninducible','inducible.inducible','uninducible.uninducible','constitutive.constitutive')
# DTemp<-merge(clusterframe3,phenDT[,summary.func(growth.rate),by=c(genos,'single_lv')],by=genos)
DT<-CloneNamedSummary(clusterframe3[,!c('aa','bb','cc')],allele=F)
#DT<-DTa[apply(data.frame(DTa[,grepincols(DTa,'clust'),with=F]),1,function(x)as.logical(prod(x%in%clustset))),]
ABgrid<-data.table(expand.grid(clustset,abset))
ABgrid[,c('Var2','Var3'):=colsplit(Var2,'\\.',c('Var2','Var3'))]
ABgrid1<-data.table(t(apply(ABgrid,1,as.character)))
landscape<-lapply(1: nrow(ABgrid),function(i){	
	Aphen<-'inducible'
	Bphen<-'constitutive'
	ABphen<-'uninducible'
	Aphen<-as.character(ABgrid[i,]$Var2)
	Bphen<-as.character(ABgrid[i,]$Var3)
	ABphen<-as.character(ABgrid[i,]$Var1)


	ABframe<-data.table(DT,t(DT[,apply(data.frame(aa_clust,bb_clust,cc_clust,aa,bb,cc),1,function(x){
		#x<-DT[,data.frame(aa_clust,bb_clust,cc_clust,aa,bb,cc)][4682,]
		df1<-data.table(pheno=as.character(unlist(x[1:3])),clone=as.character(unlist(x[4:6])))
		df<-df1[!grepl('.WT',clone)]
		A<-df[pheno%in%c(Aphen)]$clone[1]
		B<-df[clone!=A&pheno%in%(Bphen)]$clone
		if(length(A)==0)A<-NA
		if(length(B)==0)B<-NA
		list(A,B)
	})]))
	
	
	ABframe[,double_mutant:=unlist(sapply(mutant_gene,function(x)length(strsplit(x,'\\ ')[[1]])==2))]
#	print(ABframe[!is.na(V2)])
	dSub<-ABframe[double_mutant==T&cluster_new== ABphen&!is.na(V1)&!is.na(V2)]
	
	x1<-dSub[mutant_gene == 'GAL3 GAL4']
	x2<-dSub[mutant_gene == 'GAL3 GAL80']
	x3<-dSub[mutant_gene == 'GAL80 GAL4']
	aas<-unique(c(x1$aa,x3$aa))
	bbs<-unique(c(x1$bb,x2$bb))
	ccs<-unique(c(x2$cc,x3$cc))
	
	DTm<-merge(phenDT,clusterframe2,by=c(genos))
	singleframe<-DTm[single_lv==T|WT==T]
	WT<-summary.func.all(DTm[WT==T],phenocols,bylist=c('clone.named','mutant_gene','single_lv','cluster_new'))


	
	dSub2<-rbind(singleframe[(aa%in%aas|bb%in%bbs|cc%in%ccs)],DTm[clone.named%in%dSub$clone.named])
	summ<-summary.func.all(dSub2,c('fracon.gal','fracon.glu'),by=c('cluster_new','mutant_gene','single_lv','clone.named'))
	
	pre<-summ[,lapply(.SD,mean),.SDcols=c('fracon.glu_mean','fracon.gal_mean'),by='cluster_new']
	xWT<-WT$fracon.glu_mean
	yWT<-WT$fracon.gal_mean
	xendx<-pre[cluster_new==ABphen][1,2]
	yendx<-pre[cluster_new==ABphen][1,3]
	x0A<-pre[cluster_new==Aphen][1,2]
	y0A<-pre[cluster_new==Aphen][1,3]
	x0B<-pre[cluster_new==Bphen][1,2]
	y0B<-pre[cluster_new==Bphen][1,3]
	pre1<-data.table(x0=c(xWT,xWT,x0A,x0B),y0=c(yWT,yWT,y0A,y0B),xend=c(x0A,x0B,xendx,xendx),yend=c(y0A,y0B,yendx,yendx))
	arrowframe<-data.frame(apply(pre1,2,function(x)as.numeric(unlist(x))))
	arrowframe$col<-factor(c('single','single','double','double'),c('WT','single','double'))
	genotypelevels<-c('WT','GAL4','GAL3','GAL80','GAL3 GAL80','GAL80 GAL4','GAL3 GAL4')
	
	print(paste('made it here ',i))
	summ<-summary.func.all(dSub2,c('fracon.gal','fracon.glu'),by=c('cluster_new','mutant_gene','single_lv','clone.named'))
	summ[,fac:=factor(gsub('TRUE','single',gsub('FALSE','double',single_lv)),levels=c('WT','single','double'))]
	summ[,pairing:=factor(mutant_gene,levels=genotypelevels)]
	WT[,fac:=factor('WT',levels=c('WT','single','double'))]
	WT[,pairing:=factor('WT',levels=genotypelevels)]
	summ1<-rbind(summ,WT[,colnames(WT)%in%colnames(summ),with=F])
	summ1[,faceting:=paste0(Aphen,' + ',Bphen,' = ',ABphen)]
	arrowframe$faceting<-paste0(Aphen,' + ',Bphen,' = ',ABphen)
	return(list(summ1,arrowframe))
})

summs<-data.table(ldply(lapply(landscape,function(x)x[[1]])))
arrows<-(ldply(lapply(landscape,function(x)x[[2]])))
facy<-1
limits<-c(0,1.02)
ptsize<-4
strokesize<-1
custshapes<-c(1,3,4,0,2,5,6)
facx<-2


alphatab<-summs[,2/(1+exp(0.01*(.N-1))),by=c('fac','faceting')]
#alphatab<-data.table(pairing=factor(levels(summs$pairing),levels(summs$pairing)),alphas=c(1,0.8,0.8,0.8,0.3,0.3,0.3))
setnames(alphatab,'V1','alphas')
summs1<-merge(summs,alphatab,by=c('fac','faceting'))[order(pairing,decreasing=T)]
sample(letters,4)
pOKJG<-ggplot(summs1,aes(fracon.glu_mean,fracon.gal_mean,xmin=fracon.glu_lower,xmax=fracon.glu_upper,ymin=fracon.gal_lower,ymax=fracon.gal_upper))+
	theme_minimal()+
	geom_errorbar(col='grey70')+geom_errorbarh(col='grey70')+
	geom_point(aes(col= fac,shape= pairing,alpha=alphas),size=ptsize,stroke=strokesize)+geom_abline(col='grey70')+
	scale_shape_manual(values=custshapes)+
	scale_color_manual(values=factorDefault)+
	geom_segment(data= arrows,aes(y=y0,yend=yend,x=x0,xend=xend),size=0.1,col='black',arrow=arrow(type='closed',length=unit(0.3,'cm')),inherit.aes=F)+
	geom_segment(data= arrows,aes(y=y0,yend=yend,x=x0,xend=xend,col=col),size=0.5,inherit.aes=F)+
	ylim(limits)+xlim(limits)+ylab('fraction ON cells in galactose')+xlab('fraction ON cells in glucose')+
	facet_wrap(~faceting,ncol=3)

# w<-7;h<-4.5
# ggsave(paste0(figout,'1707xx-180724-pDYHQ-GeneExpressionClusterHierarchy-Rectangles.png'), pDYHQ,width=w,height=h)


# look at weirdly categorized clones
combo<-'inducible + inducible = constitutive'
clones1<-unique(summs[faceting%in%combo]$clone.named)
clones0<-unique(summs[faceting%in%combo&single_lv==F]$clone.named)
clones<-clones1#clones1%w/o%c('GAL3.WT GAL80.WT GAL4.WT','GAL3.WT GAL80.29 GAL4.WT')
# clones<-c(dSub$clone.named)
# DTxx<-merge(phenDT[clone.named%in%clones],clusterframe3[clone.named%in%clones],by=c('aa','bb','cc','clone.named'))
# clones<-unique(DTxx[fracon.gal>0.5]$clone.named)
DT1<-merge(densDT[clone.named%in%clones],clusterframe3[clone.named%in%clones],by=c('aa','bb','cc','clone.named'))
mDT<-melt(DT1[,c('cluster_new','clone.named',dens.cols,'aa','bb','cc','single_lv'),with=F],id=c('cluster_new','clone.named','aa','bb','cc','single_lv'))
denscolconv(mDT)
mDT[,dummyline:=paste0('z',sugar)]
mDT[,pheno:='fluorescence']

ylimits<-c(0,0.3)
ybreaks<-c(0.0,0.1,0.2,0.3)
ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
mDT1<-mDT[,mean(value,na.rm=T),by=c('clone.named','dummyline','AU.FL','sugar','aa','bb','cc','single_lv')];setnames(mDT1,'V1','density_mean')
mDT1[density_mean>=max(ylimits),density_mean:=max(ylimits)][order(sugar,AU.FL,decreasing=T)]

if(writeplot==T){
	ggplot(mDT1,aes(AU.FL,density_mean))+theme_minimal()+
		geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
		facet_grid(bb+aa~cc)+scale_colour_manual(values=c('red','blue','orange1','steelblue','grey90','grey90'),na.value='transparent')+
		ylab('density') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
		scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
		theme(legend.title=element_blank())+
		geom_point(data=data.table(mDT1[,unique(single_lv),by=c('aa','bb','cc')][,c('x','y'):=list(1,0.2)]),aes(x=x,y=y,col=V1),size=3,inherit.aes=F)+
		theme(strip.text.x = element_text(size = 10))
	
	ggplot(summs1[faceting==combo],aes(fracon.glu_mean,fracon.gal_mean,xmin=fracon.glu_lower,xmax=fracon.glu_upper,ymin=fracon.gal_lower,ymax=fracon.gal_upper))+
		theme_minimal()+
		geom_errorbar(col='grey70')+geom_errorbarh(col='grey70')+
		geom_point(aes(col= fac,shape= pairing,alpha=alphas),size=ptsize,stroke=strokesize)+geom_abline(col='grey70')+
		scale_shape_manual(values=custshapes)+
		scale_color_manual(values=factorDefault)+
		geom_segment(data= arrows[arrows$faceting%in%combo,],aes(y=y0,yend=yend,x=x0,xend=xend),size=0.1,col='black',arrow=arrow(type='closed',length=unit(0.3,'cm')),inherit.aes=F)+
		geom_segment(data= arrows[arrows$faceting%in%combo,],aes(y=y0,yend=yend,x=x0,xend=xend,col=col),size=0.5,inherit.aes=F)+
		ylim(limits)+xlim(limits)+ylab('fraction ON cells in galactose')+xlab('fraction ON cells in glucose')+
		facet_wrap(~faceting,ncol=3)
}




# look at GAL3 var
g3var<-'GAL80.01'
ab<-phenDT[cc%in%c(g3var,'GAL80.WT'),summary.func(yfp.mean.gal),by=genos]
a<-ab[cc=='GAL80.WT']
b<-ab[cc==g3var]
abm<-merge(a,b,by=c('aa','bb'))
ggplot(abm,aes(mean.x,mean.y))+geom_point()+geom_abline()
tests<-(data.table(apply(t(abm[,apply(data.frame(m1=mean.x,m2=mean.y,s1=sd.x,s2=sd.y,n1=N.x,n2=N.y),1,function(x){
	
	t.test2(m1=x[1],m2=x[2],s1=x[3],s2=x[4],n1=x[5],n2=x[6],list=T)})
	]),2,as.numeric)))
	
colnames(tests)<-names(t.test2(1,1,1,1,1,1))
tests[,p.adj:=p.adjust(p.value,'fdr')]
tests[,sig:=p.adj<0.05]
abmt<-data.table(abm,tests)
ggplot(abmt,aes(mean.x,mean.y,xmin=lower.x,xmax=upper.x,ymin=lower.y,ymax=upper.y,col=sig))+geom_point()+geom_errorbar()+geom_errorbarh()+geom_abline()
abmt[order(p.adj)]
aas<-unique(c(abmt[sig==T]$aa,'GAL4.WT'))
ccs<-unique(c(abmt[sig==T]$cc,'GAL80.WT','GAL80.delta'))

DTxx<-merge(densDT, clusterframe3,by=genos)
DT1<-DTxx[bb%in%c(g3var,'GAL3.WT')&(aa%in%aas)&cc=='GAL80.WT']
mDT<-melt(DT1[,c('cluster_new','clone.named',dens.cols,'aa','bb','cc'),with=F],id=c('cluster_new','clone.named','aa','bb','cc'))
denscolconv(mDT)
mDT[,dummyline:=paste0('z',sugar)]
mDT[,pheno:='fluorescence']

ylimits<-c(0,0.3)
ybreaks<-c(0.0,0.1,0.2,0.3)
ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
mDT1<-mDT[,mean(value,na.rm=T),by=c('clone.named','dummyline','AU.FL','sugar','aa','bb','cc')];setnames(mDT1,'V1','density_mean')
mDT1[density_mean>=max(ylimits),density_mean:=max(ylimits)][order(sugar,AU.FL,decreasing=T)]
if(writeplot==T){
	ggplot(mDT1,aes(AU.FL,density_mean))+theme_minimal()+
		geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
		facet_grid(bb+cc~aa)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
		ylab('density') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
		scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
		theme(legend.title=element_blank())+
		theme(strip.text.x = element_text(size = 10))
}


p2<-ggplot(summ1,aes(yfp.mean.glu_mean,fracon.gal_mean,xmin=yfp.mean.glu_lower,xmax=yfp.mean.glu_upper,ymin=fracon.gal_lower,ymax=fracon.gal_upper))+
	theme_minimal()+
	geom_errorbar()+geom_errorbarh()+geom_vline(xintercept=1000)+
	geom_point(aes(col= single_lv,shape=mutant_gene))
#	ylim(limits)+xlim(limits)
do.call(grid.arrange,list(p1,p2))

################
### show hierarchy of clusters
################

PhenoHierarchyRect<-function(notframe,clustCols,clusterlevels,plot=T,varname=NULL){
	cCols<-c('A_clust','B_clust')
	test1<-notframe[,.N,by=clustCols]
	test2<-notframe[,.N,by=c(clustCols,'cluster_new')]
	test3a <-merge(test1,test2,by=clustCols)
	setnames(test3a,clustCols,cCols)
	test3a[, A_clust:=factor(A_clust,levels=clusterlevels)]
	test3a[, B_clust:=factor(B_clust,levels=clusterlevels)]
	test3a[, cluster_new:=factor(cluster_new,levels=clusterlevels)]
	test3<-test3a[order(A_clust, B_clust,cluster_new,decreasing=T)]
	test3[,A:=.N,by='A_clust']
	test3[,B:=.N,by=cCols]
	test3[,C:=.N,by=c(cCols,'cluster_new')]
	test4<-test3[order(A_clust,B_clust,cluster_new,decreasing=T)]
	test4[,seq:=1:nrow(test4)]
	test4[,dupA:=!duplicated(A_clust)]
	test4[,dupB:=!duplicated(B_clust),by=A_clust]
	test4[,dupC:=!duplicated(cluster_new),by=cCols]
	test4[,A2:=seq+A]
	test4[,B2:=seq+B]
	test4[,C2:=seq+C]
	A<-test4[dupA==T]
	A[,c('xmin','xmax','ymin','ymax','single_clust'):=list(0,1,seq,A2,A_clust)]
	B<-test4[dupB==T]
	B[,c('xmin','xmax','ymin','ymax','single_clust'):=list(1,2,seq,B2,B_clust)]
	C<-test4[dupC==T]
	C[,Ctest2:=(C2-seq)*N.y/N.x*.N,by=cCols]
	C[,C2:=cumsum(Ctest2)+1]
	newseq<-c(1,C$C2[1:(nrow(C)-1)])
	C[,seq:=newseq]
	C[,Ctest2:=NULL]
	C[,c('xmin','xmax','ymin','ymax','single_clust'):=list(2,3,seq,C2,cluster_new)]
	DTr<-rbind(A,B,C)
	miny<-min(DTr[,'ymin',with=F])
	DTr[,c('ymin','ymax'):=lapply(.SD,function(x)(x-miny)),.SDcols=c('ymin','ymax')]
	maxy<-max(DTr[,'ymax',with=F])
	DTr[,c('ymin','ymax'):=lapply(.SD,function(x)(x/maxy)),.SDcols=c('ymin','ymax')]
	p<-ggplot(DTr,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))+geom_rect(aes(fill= single_clust),col='grey40')+scale_fill_manual(values=c(rainbow(5),'grey70'),na.value='white')+
		theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
		theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),panel.border=element_blank(),axis.line=element_blank())+NULL
	if(!is.null(varname))DTr[,var:=varname]
	if(plot==T)	return(p)
	if(plot==F)return(DTr)
}


DT2<-merge(phenDT,clusterframe,by=c('aa','bb','cc','clone.named'))
setnames(DT2,'HDB2','cluster_new')

clustCols<-c('bb_clust','cc_clust')
notaa<-(DT2[aa_clust=='inducible',c(clustCols,'cluster_new'),with=F])
clusterlevels<-c('uninducible','inducible','constitutive','weakly inducible','const, low expr','outlying')
notframe<-notaa
DTr1<-PhenoHierarchyRect(notframe=notframe,clustCols=clustCols,clusterlevels=clusterlevels,plot=F,varname='GAL3 v GAL80')

clustCols<-c('aa_clust','cc_clust')
notbb<-(DT2[bb_clust=='inducible',c(clustCols,'cluster_new'),with=F])
clusterlevels<-c('uninducible','inducible','constitutive','weakly inducible','const, low expr','outlying')
notframe<-notbb
DTr2<-PhenoHierarchyRect(notframe=notframe,clustCols=clustCols,clusterlevels=clusterlevels,plot=F,varname='GAL4 v GAL80')

clustCols<-c('aa_clust','bb_clust')
notcc<-(DT2[cc_clust=='inducible',c(clustCols,'cluster_new'),with=F])
clusterlevels<-c('uninducible','inducible','constitutive','weakly inducible','const, low expr','outlying')
notframe<-notcc
DTr3<-PhenoHierarchyRect(notframe=notframe,clustCols=clustCols,clusterlevels=clusterlevels,plot=F,varname='GAL4 v GAL3')

DTr<-rbind(DTr1,DTr2,DTr3)
DTr[,`gene expression cluster`:=single_clust]
sample(letters,4)

# # failed attempt to make axis labels
# xlabs<-c('GAL3','GAL80','GAL80 + GAL3','GAL4','GAL3','GAL4 + GAL3','GAL4','GAL80','GAL4 + GAL80')
# dum<-DTr[,unique(var),by=c('xmin','xmax')];setnames(dum,'V1','var')
# dum[,xlab:=xlabs]

pDYHQ<-ggplot(DTr,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))+geom_rect(aes(fill= `gene expression cluster`),col='grey40')+scale_fill_manual(values=c(rainbow(5),'grey70'),na.value='white')+
	theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+
	theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),panel.border=element_blank(),axis.line=element_blank())+
	facet_grid(~var)
# w<-7;h<-4.5
# ggsave(paste0(figout,'1707xx-180724-pDYHQ-GeneExpressionClusterHierarchy-Rectangles.png'), pDYHQ,width=w,height=h)



################
### constitutive GAL4 + unindubible GAL80 yield a lot of uncharacterized profiles. look closer
################
denscol_fun<-function(DT, vnull=T, meltframe=F,bycolsi=NULL){
	if(meltframe ==F){mDT<-DT}else{
		if(is.null(bycolsi))stop('need to provide columns to melt on')
		mDT<-melt(DT,id=bycolsi)
	}
	mDT[,variable:=gsub('dens__','',variable)]
	mDT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
	mDT[,sugar:=paste0('g',sugar)]
	setnames(mDT,'value','density')
	if(vnull ==T)mDT[,variable:=NULL]
	return(data.table(mDT))
}

# DT<-densDT[,lapply(.SD,mean,na.rm=T),.SDcols=dens.cols,by='clone.named']
# mDT<-denscol_fun(DT,meltframe=T,bycolsi='clone.named')
# mDT1<-CloneNamedSummary(mDT,allele=F)
DT2<-merge(phenDT,clusterframe,by=c('aa','bb','cc','clone.named'))
setnames(DT2,'HDB2','cluster_new')

clones1<-DT2[aa_clust=='constitutive'&cc_clust=='uninducible',c('aa','cc')]
ctrls<-data.table(aa=c('GAL4.WT','GAL4.delta'),cc=c('GAL80.WT','GAL80.delta'))
clones2<-rbind(clones1,ctrls)
clones2[,dumm:=1]
clones3<-data.table(expand.grid(clones2[,unique(dumm),by=c('aa','cc')]))[,V1:=NULL]
clones3[,bb:='GAL3.WT']
clones3[,clone.named:=applyPaste(data.frame(bb,cc,aa),' ')]
clones<-unique(clones3$clone.named)

mDT2<-(melt(densDT[clone.named%in%clones,c(dens.cols,'clone.named'),with=F],id='clone.named'))
mDT1<-merge(mDT2,mDT2[,summary.func(value),by=c('clone.named','variable')],by=c('clone.named','variable'))
denscolconv(mDT1)
mDT<-CloneNamedSummary(mDT1,allele=F)
mDT[,dummyline:=paste0('z',sugar)]
mDT[,aa:=factor(aa,levels=c('GAL4.WT','GAL4.delta',unique(clones1$aa)))]
mDT[,cc:=factor(cc,levels=c('GAL80.WT','GAL80.delta',unique(clones1$cc)))]
mDT[mean>0.3,mean:=0.3]
sample(letters,4)
pIWSZ<-ggplot(mDT,aes(AU.FL,mean))+
	geom_point(data= mDT,aes(AU.FL,value,col=sugar),size=0.1,shape=21)+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(aa~cc)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density of observations') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))+
	ylim(c(0,0.3))
# w<-7;h<-4.5
# ggsave(paste0(figout,'1707xx-180724-pIWSZ-GoF_Unclustered-GeneExpressionLinePlots.png'), pIWSZ,width=w,height=h)

################
### weakly inducbile GAL4 + constitutive GAL80 yield a lot of uncharacterized profiles. look closer
################
DT2<-merge(phenDT,clusterframe,by=c('aa','bb','cc','clone.named'))
setnames(DT2,'HDB2','cluster_new')

clones1<-DT2[aa_clust=='weakly inducible'&cc_clust=='constitutive',c('aa','cc')]
ctrls<-data.table(aa=c('GAL4.WT','GAL4.delta'),cc=c('GAL80.WT','GAL80.delta'))
clones2<-rbind(clones1,ctrls)
clones3<-data.table(expand.grid(clones2))
clones3[,bb:='GAL3.WT']
clones3[,clone.named:=applyPaste(data.frame(bb,cc,aa),' ')]
clones<-unique(clones3$clone.named)
G80Rank<-notin(phenDT[clone.named%in%clones&aa=='GAL4.WT',mean(fracon.glu),by=c('cc','clone.named')][order(V1)]$cc,c('GAL80.WT','GAL80.delta'))
G80lev<-c('GAL80.WT',G80Rank,'GAL80.delta')

mDT2<-(melt(densDT[clone.named%in%clones,c(dens.cols,'clone.named'),with=F],id='clone.named'))
mDT1<-merge(mDT2,mDT2[,summary.func(value),by=c('clone.named','variable')],by=c('clone.named','variable'))
denscolconv(mDT1)
mDT<-CloneNamedSummary(mDT1,allele=F)
mDT[,dummyline:=paste0('z',sugar)]
mDT[,aa:=factor(aa,levels=c('GAL4.WT','GAL4.delta',unique(clones1$aa)))]
mDT[,cc:=factor(cc,levels=G80lev)]
mDT[mean>0.3,mean:=0.3]
sample(letters,4)
pMKJV<-ggplot(mDT,aes(AU.FL,mean))+
	geom_point(data= mDT,aes(AU.FL,value,col=sugar),size=0.1,shape=21)+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(aa~cc)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density of observations') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))+
	ylim(c(0,0.3))
	

# w<-10;h<-4.5
# ggsave(paste0(figout,'1707xx-180724-pMKJV-GAL4weakInduceAndConstGAL80-GeneExpressionLinePlots.png'), pMKJV,width=w,height=h)


###################################################
### analysis of GAL80 repression scores variants in GAL3 backgrounds: allelic series
###################################################

rm(clusterframe)
datloc<-paste0(outputdata,'1707xx-180704-HDB_Clusters_Based_On_GAL_Distributions.rData')
load(datloc)
phenDT0<-merge(clusterframe,phenDT,by=c('aa','bb','cc','clone.named'))

################
### get GAL80 and GAL3 functional scores and plot in 2d
################

######
### GAL80's in different GAL3 backgrounds
######
dTemp<-phenDT0[aa_clust=='inducible'&(cc_clust!='uninducible'&cc_clust!='const, low expr')]
Y<-dTemp[,summary.func(logit(fracon.gal+0.001),'Y'),by=c('bb','cc','bb_clust','cc_clust')]
X<-dTemp[,summary.func(logit(fracon.glu+0.001),'X'),by=c('bb','cc','bb_clust','cc_clust')]
XY<-merge(X,Y,by=c('bb','cc','cc_clust','bb_clust'))
XY[,`GAL3 GAL80 combination`:=applyPaste(data.frame('GAL3',bb_clust,'GAL80',cc_clust),' ')]
XY1<-copy(XY)
sample(letters,4)
pVRHT<-ggplot(XY1,aes(X_mean,Y_mean))+geom_point(aes(col=`GAL3 GAL80 combination`),shape=21)+
	ylab('logit fraction ON') + xlab('')+
	theme(legend.position=c(0.5,0.2),legend.title=element_text(size=6),legend.text=element_text(size=6))

dTemp<-phenDT0[aa_clust=='inducible'&(cc_clust!='uninducible'&cc_clust!='const, low expr')]
Y<-dTemp[,summary.func(growth.rate,'Y'),by=c('bb','cc','bb_clust','cc_clust')]
X<-dTemp[,summary.func(logit(fracon.glu+0.001),'X'),by=c('bb','cc','bb_clust','cc_clust')]
XY<-merge(X,Y,by=c('bb','cc','cc_clust','bb_clust'))
XY[,`GAL3 GAL80 combination`:=applyPaste(data.frame('GAL3',bb_clust,'GAL80',cc_clust),' ')]
XY2<-copy(XY)
sample(letters,4)
pYIOV<-ggplot(XY2,aes(X_mean,Y_mean))+geom_point(aes(col=`GAL3 GAL80 combination`),shape=21)+
	ylab(expression(paste('growth rate',h^-1))) + xlab('logit fraction ON in glucose')+ theme(legend.position='none')

XY1[,fac:='logit fraction ON']
XY2[,fac:='growth rate']

XY3<-rbind(XY1,XY2)
pUOBL <-ggplot(XY3,aes(X_mean,Y_mean))+geom_point(aes(col=`GAL3 GAL80 combination`),shape=21)+
	ylab('phenotypes in galactose') + xlab('logit fraction ON in glucose')+facet_wrap(~fac,ncol=1,scales='free',strip.position = "left")+
	theme(legend.position=c(0.5,0.1),legend.title=element_text(size=6),legend.text=element_text(size=6))
w<-4.6;h<-6
ggsave(paste0(figout,'1707xx-180809-pUOBL-AllelicSeriesGAL80vsGAL3-GalPhenotypesVsLogitFraconGlu.png'), pUOBL,width=w,height=h)

######
### GAL80 repression with GAL3 vs no GAL3
######
dTemp<-phenDT0[aa_clust=='inducible'&(cc_clust!='uninducible'&cc_clust!='const, low expr')]
A<-dTemp[,summary.func(fracon.glu),by=c('cc','bb','bb_clust','cc_clust')]

X<-A[bb_clust=='uninducible',summary.func(mean,'X'),by='cc']
Y<-A[bb_clust=='inducible',summary.func(mean,'Y'),by='cc']
XY<-merge(X,Y,by='cc')
ggplot(XY,aes(X_mean,Y_mean,ymax=Y_upper,ymin=Y_lower,xmax=X_upper,xmin=X_lower))+geom_errorbarh()+geom_errorbar()+geom_point()+geom_abline()

attach(XY)
test<-cbind(XY,XY[,t(apply(data.frame(X_mean,Y_mean,X_sd,Y_sd,X_N,Y_N),1,function(x){
	t.test2(x[1],x[2],x[3],x[4],x[5],x[6])}))])
test[,sig:=p.value<0.05]
ggplot(test,aes(X_mean,Y_mean,ymax=Y_upper,ymin=Y_lower,xmax=X_upper,xmin=X_lower,col=sig))+geom_errorbarh()+geom_errorbar()+geom_point()+geom_abline()


################
### sensitized background plotting
################


######
### get data together
######
# GAL80
clonesub1<-phenDT0[(aa_clust=='inducible'|aa%in%c('GAL4.01','GAL4.03','GAL4.06'))&bb=='GAL3.WT'&cc_clust!='uninducible'&cc_clust!='const, low expr']
clonesub1[,outlier:=clone.named%in%c('GAL3.WT GAL80.01 GAL4.01','GAL3.WT GAL80.02 GAL4.01')]
clonesub<-clonesub1[outlier==F]
clonesub[,lyfpglu:=log10(yfp.mean.glu)]
WT<-clonesub[aa_clust=='inducible',c('cc',grepincols(clonesub,c('fracon.glu','yfp.mean.glu','lyfpglu'))),with=F]
MUT<-clonesub[aa%in%c('GAL4.01','GAL4.03','GAL4.06'),c('cc',grepincols(clonesub,c('fracon.glu','yfp.mean.glu','lyfpglu'))),with=F]
G80_mut<-copy(MUT)
G80_wt<-copy(WT)
# DTm<-merge(WT,MUT,by='cc',allow.cartesian=T)
# PCAmod<-PCA(grepcols(DTm,'glu'))
# summary(PCAmod)
# dims<-data.table(cc=DTm$cc,PCAmod$ind$coord)
# dumsumm<-dims[,summary.func(Dim.1,'cc_Dim.1'),by='cc']


# GAL3
clonesub1<-phenDT0[(aa_clust=='inducible'|aa%in%c('GAL4.01','GAL4.03','GAL4.06'))&cc=='GAL80.WT']
clonesub1[,outlier:=clone.named%in%c('GAL3.WT GAL80.01 GAL4.01','GAL3.WT GAL80.02 GAL4.01')]
clonesub<-clonesub1[outlier==F]
clonesub[,lyfpgal:=log10(yfp.mean.gal)]
WT<-clonesub[aa_clust=='inducible',c('bb',grepincols(clonesub,c('fracon.gal','yfp.mean.gal','lyfpgal'))),with=F]
MUT<-clonesub[aa%in%c('GAL4.01','GAL4.03','GAL4.06'),c('bb',grepincols(clonesub,c('fracon.gal','yfp.mean.gal','lyfpgal'))),with=F]
G3_mut<-copy(MUT)
G3_wt<-copy(WT)



######
### overview of GAL80 in 2 different inducible-classified GAL4 backgrounds
######

clonesub1<-phenDT0[(aa_clust=='inducible'|aa%in%c('GAL4.01','GAL4.03','GAL4.06'))&bb=='GAL3.WT'&cc_clust!='uninducible'&cc_clust!='const, low expr']
clonesub1[,outlier:=clone.named%in%c('GAL3.WT GAL80.01 GAL4.01','GAL3.WT GAL80.02 GAL4.01')]
clonesub<-clonesub1[outlier==F]
# non-linear relationship between the two indubible classes with useful information for each
WT<-clonesub[aa_clust=='inducible',summary.func(fracon.glu,'WT'),by=c('cc')]
MUT<-clonesub[aa%in%c('GAL4.01','GAL4.03','GAL4.06'),summary.func(fracon.glu,'MUT'),by=c('cc')]
tog1<-merge(WT,MUT,by=c('cc'))
p1<-ggplot(tog1,aes(WT_mean,MUT_mean,xmin=WT_lower,xmax=WT_upper,ymin=MUT_lower,ymax=MUT_upper))+geom_errorbarh()+geom_errorbar()+geom_point()+
	ylab('fracon.glu')+xlab('fracon.glu')+geom_abline()+xlim(c(0,1))+ylim(c(0,1))

#asdfasdf
WT<-clonesub[aa_clust=='inducible',summary.func(logit(fracon.glu+ 0.001),'WT'),by=c('cc')]
MUT<-clonesub[aa%in%c('GAL4.01','GAL4.03','GAL4.06'),summary.func(logit(fracon.glu+0.001),'MUT'),by=c('cc')]
tog<-merge(WT,MUT,by=c('cc'))
logTog<-copy(tog)
p2<-ggplot(tog,aes(WT_mean,MUT_mean,xmin=WT_lower,xmax=WT_upper,ymin=MUT_lower,ymax=MUT_upper))+geom_errorbarh()+geom_errorbar()+geom_point()+geom_smooth(method='lm')+
	ylab('logit(fracon.glu)')+xlab('logit(fracon.glu)')
left<-'GAL80 variants in weakly inducible GAL4 background'
bottom<-'GAL80 variants in WT GAL4 background'
pMKRL<-list(p1,p2,left=left,bottom=bottom)
do.call(grid.arrange,pMKRL)
# sample(letters,4)
# w<-4;h<-5.3
# ggsave(paste0(figout,'1707xx-180710-pMKRL-AllelicSeriesGAL80vsGAL3.png'), do.call(grid.arrange,pMKRL),width=w,height=h)

WT<-clonesub[aa_clust=='inducible',summary.func((yfp.mean.glu),'WT'),by=c('cc')]
MUT<-clonesub[aa%in%c('GAL4.01','GAL4.03','GAL4.06'),summary.func((yfp.mean.glu),'MUT'),by=c('cc')]
tog<-merge(WT,MUT,by=c('cc'))
yfpTog<-copy(tog)
# ggplot(tog,aes(WT_mean,MUT_mean,xmin=WT_lower,xmax=WT_upper,ymin=MUT_lower,ymax=MUT_upper))+geom_errorbarh()+geom_errorbar()+geom_point()+
	# ylab('yfp sig in glu\nGAL80 variants in weakly inducible GAL4 background')+xlab('yfp sig in glu\nGAL80 variants in WT GAL4 background')


WT<-clonesub[aa_clust=='inducible',summary.func(log10(yfp.mean.glu),'WT'),by=c('cc')]
MUT<-clonesub[aa%in%c('GAL4.01','GAL4.03','GAL4.06'),summary.func(log10(yfp.mean.glu),'MUT'),by=c('cc')]
tog<-merge(WT,MUT,by=c('cc'))
yfplogTog<-copy(tog)
# ggplot(tog,aes(WT_mean,MUT_mean,xmin=WT_lower,xmax=WT_upper,ymin=MUT_lower,ymax=MUT_upper))+geom_errorbarh()+geom_errorbar()+geom_point()+
	# ylab('log yfp glu\nGAL80 variants in weakly inducible GAL4 background')+xlab('log yfp glu\nGAL80 variants in WT GAL4 background')



######
### OK way of showing GAL80 across GAL3 different backgrounds
######

# PCA on GAL3
DTm<-merge(G3_wt, G3_mut,by='bb',allow.cartesian=T)
PCAmod<-PCA(DTm[,!'bb',with=F],graph=F)
summary(PCAmod)
dims<-data.table(bb=DTm$bb,PCAmod$ind$coord)
dimsumm<-dims[,summary.func(Dim.1,'bb_Dim.1'),by='bb']


setnames(G80_mut,c('fracon.glu','yfp.mean.glu','lyfpglu'),c('G80.fON','G80.mean.yfp','G80.lyfp'))
G80_mut[,G80.logit:=logit(G80.fON+0.001)]
G80_mut_summ<-summary.func.all(G80_mut,notin(colnames(G80_mut),'cc'),bylist='cc')
setnames(G3_mut,c('fracon.gal','yfp.mean.gal','lyfpgal'),c('G3.fON','G3.mean.yfp','G3.lyfp'))
G3_mut[,G3.logit:=logit(G3.fON+0.001)]
DTm<-merge(merge(phenDT0,dimsumm,by='bb'), G80_mut_summ,by='cc')
dSumm<-DTm[aa_clust=='inducible',summary.func(growth.rate,'GR'),by=c('clone.named','bb','cc',grepincols(DTm,c('Dim','80')))]
dSumm[,x:=rnorm(.N,mean=(-1*G80.logit_mean),sd= G80.logit_se),by='clone.named']
dSumm[,y:=rnorm(.N,mean=bb_Dim.1_mean,sd=bb_Dim.1_se*1.96),by='clone.named']
dSumm[GR_mean<0.06,GR_mean:=0.04]
# ggplot(dSumm,aes(x,y))+geom_point(aes(col= GR_mean),shape=21)+scale_colour_gradientn(colours=c('cornflowerblue','yellow','indianred'))

# classify GAL3
dSumm[bb_Dim.1_mean <0,GAL3_cat:='GAL3 deletion']
dSumm[round(bb_Dim.1_mean,2)>0.53,GAL3_cat:='GAL3 wild-type']
dSumm[bb_Dim.1_mean<0.54&bb_Dim.1_mean>-1,GAL3_cat:='GAL3 mildly detrimental']
# ggplot(dSumm[order(GAL3_cat)],aes(x=x,y=GR_mean,col= GAL3_cat))+geom_point()+geom_smooth(se=F)
# ggplot(dSumm[order(GAL3_cat)],aes(x=x,y=GR_mean,col=bb))+geom_point()+geom_smooth(se=F)
ggplot(dSumm,aes(x=x,y=GR_mean))+geom_point()+geom_smooth(se=F)+ 
	facet_wrap(bb~.,ncol=6)+
	xlab('GAL80 repression strength')

# make plot
dSumm2<-dSumm[,summary.func(GR_mean),by=c('G80.logit_mean','G80.logit_upper','G80.logit_lower','GAL3_cat')]
ggplot(dSumm2,aes(-1* G80.logit_mean,mean,xmax=(-1* G80.logit_upper),xmin=-1* G80.logit_lower,ymax=upper,ymin=lower))+geom_point()+geom_errorbarh()+geom_errorbar()+
	geom_smooth()+ylab(expression(paste('mean growth rate ',h^-1)))+xlab('GAL80 repression strength\n (-) log odds (fraction ON in glucose)')+
	facet_grid(~GAL3_cat)



######
### Do above a different way
######

G3_mut2<-data.table(bb=G3_mut$bb,G3_mut[,lapply(.SD,scale),.SDcols=notin(colnames(G3_mut),'bb')])
G80_mut2<-data.table(cc=G80_mut$cc,G80_mut[,lapply(.SD,scale),.SDcols=notin(colnames(G80_mut),'cc')])
G3<-summary.func.all(G3_mut2,notin(colnames(G3_mut2),'bb'),'bb')
G80<-summary.func.all(G80_mut2,notin(colnames(G80_mut2),'cc'),'cc')
mG80<-melt(G80,id='cc')
mG80[,c('pheno','stat'):=colsplit(variable,'_',c('pheno','stat'))]; mG80[,variable:=NULL]
left<-c('cc','pheno')
right<-'stat'
cG80<-dcast(mG80,castform(left,right),value.var='value')
clonesub1<-phenDT0[(aa_clust=='inducible'|aa%in%c('GAL4.01','GAL4.03','GAL4.06'))&cc_clust!='uninducible'&cc_clust!='const, low expr']
clonesub1[,outlier:=clone.named%in%c('GAL3.WT GAL80.01 GAL4.01','GAL3.WT GAL80.02 GAL4.01')]
clonesub<-clonesub1[outlier==F]
grDT<-(clonesub[,summary.func(growth.rate,'GR'),by=c('cc','bb')])
grDT<-clonesub[,c('growth.rate','cc','bb'),with=F]
DTm1<-merge(grDT,cG80,by='cc',allow.cartesian=T)
qdt<-dSumm[,unique(GAL3_cat),by='bb'];setnames(qdt,'V1','GAL3_cat')
DTm<-merge(DTm1,qdt,by='bb')
sDTm<-DTm[,summary.func(growth.rate,'GR'),by=c('cc','pheno','GAL3_cat',names(summary.func(1:2)))]
mutsDTm<-copy(sDTm)
# ggplot(sDTm,aes(mean,GR_mean))+geom_point()+facet_grid(GAL3_cat~pheno)



setnames(G80_wt,c('fracon.glu','yfp.mean.glu','lyfpglu'),c('G80.fON','G80.mean.yfp','G80.lyfp'))
G80_wt[,G80.logit:=logit(G80.fON+0.001)]
setnames(G3_wt,c('fracon.gal','yfp.mean.gal','lyfpgal'),c('G3.fON','G3.mean.yfp','G3.lyfp'))
G3_wt[,G3.logit:=logit(G3.fON+0.001)]
G3_wt2<-data.table(bb=G3_wt$bb,G3_wt[,lapply(.SD,scale),.SDcols=notin(colnames(G3_wt),'bb')])
G80_wt2<-data.table(cc=G80_wt$cc,G80_wt[,lapply(.SD,scale),.SDcols=notin(colnames(G80_wt),'cc')])
G3<-summary.func.all(G3_wt2,notin(colnames(G3_wt2),'bb'),'bb')
G80<-summary.func.all(G80_wt2,notin(colnames(G80_wt2),'cc'),'cc')
mG80<-melt(G80,id='cc')
mG80[,c('pheno','stat'):=colsplit(variable,'_',c('pheno','stat'))]; mG80[,variable:=NULL]
left<-c('cc','pheno')
right<-'stat'
cG80<-dcast(mG80,castform(left,right),value.var='value')
clonesub1<-phenDT0[(aa_clust=='inducible'|aa%in%c('GAL4.01','GAL4.03','GAL4.06'))&cc_clust!='uninducible'&cc_clust!='const, low expr']
clonesub1[,outlier:=clone.named%in%c('GAL3.WT GAL80.01 GAL4.01','GAL3.WT GAL80.02 GAL4.01')]
clonesub<-clonesub1[outlier==F]
grDT<-(clonesub[,summary.func(growth.rate,'GR'),by=c('cc','bb')])
grDT<-clonesub[,c('growth.rate','cc','bb'),with=F]
DTm1<-merge(grDT,cG80,by='cc',allow.cartesian=T)
qdt<-dSumm[,unique(GAL3_cat),by='bb'];setnames(qdt,'V1','GAL3_cat')
DTm<-merge(DTm1,qdt,by='bb')
sDTm<-DTm[,summary.func(growth.rate,'GR'),by=c('cc','pheno','GAL3_cat',names(summary.func(1:2)))]
WTsDTm<-copy(sDTm)
mutsDTm[,`GAL80 phenotype background`:='GAL4 weakly inducible']
WTsDTm[,`GAL80 phenotype background`:='GAL4 WT-like-inducible']
combDT1<-rbind(mutsDTm,WTsDTm)
qdt<-data.table(pheno=c('G80.fON','G80.logit','G80.mean.yfp','G80.lyfp'),pheno2=c('phenotypes in glucose:\nfraction ON','\nlog odds fraction ON','\nmean YFP signal, AU','\nlog(mean YFP signal, AU)'))
combDT<-merge(combDT1,qdt,by='pheno')
combDT[,pheno2:=factor(pheno2,levels=qdt$pheno2)]
sample(letters,4)
pDJFT<-ggplot(combDT,aes(mean,GR_mean,xmin=lower,xmax=upper,ymin=GR_upper,ymax=GR_lower))+geom_errorbar(col='grey50')+geom_errorbarh(col='grey50')+
	geom_point(aes(col=`GAL80 phenotype background`),shape=21)+
	facet_grid(GAL3_cat~ pheno2)+scale_colour_manual(values=c('orange1','cornflowerblue'))+
	theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))+
	ylab(expression(paste('mean growth rate ',hr^-1)))+xlab('scaled repression phenotype for GAL80 across given inducible GAL4 backgrounds')+
	NULL
w<-10;h<-5
ggsave(paste0(figout,'1707xx-180710-pDJFT-AllelicSeriesGAL80vsGAL3.png'), pDJFT,width=w,height=h)


