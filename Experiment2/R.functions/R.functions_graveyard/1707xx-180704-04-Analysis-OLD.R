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
### technical overview
###################################################
DT0<-copy(phenDT)
phenosh<-c('GR','PI','frac ON glu','frac ON gal','mean YFP glu','mean YFP gal')
setnames(DT0,phenocols,phenosh)
# use this for factors: scale_colour_manual(values=c('orange1','cornflowerblue','chartreuse3','blue'))
sample(letters,4)
pBRDL<-pTechrep2(DT0,'clone.named',phenosh=phenosh,plot=T)+theme(axis.text.x=element_text(angle=30,hjust=1,size=8),axis.text.y=element_text(size=8))
w<-8.4;h<-5.6
ggsave(paste0(figout,'1707xx-180817-pBRDL-pTechRep-PairwiseCombGenetics.png'), pBRDL,width=w,height=h)

###################################################
### simple regression models with growth rate as a fucntion of expression
###################################################
DT0<-copy(phenDT)
summ1<-CloneNamedSummary(summary.func.all(DT0, phenocols,'clone.named'),allele=F)
mod1<-lm(growth.rate~fracon.gal+fracon.glu,data=DT0)
summary(mod1)
mod2<-lm(growth.rate_mean~fracon.gal_mean+fracon.glu_mean,data=summ1)
summary(mod2)
pred<-predict(mod2,summ1,se.fit=T)
summ1[,c('GR.predicted_mean','GR.predicted_sd'):=list(pred$fit,pred$se.fit)]
summ1[,c('pred_mean','pred_sd'):=list(pred$fit,pred$se.fit)]
summ1[,c('xmax','xmin'):=list(pred_mean+pred_sd,pred_mean-pred_sd)]
singles<-summ1[!is.na(single_locus)]
singles[,`single mutated\nlocus`:=single_locus]
sample(letters,4)
pNCDK<-ggplot(summ1,aes(pred_mean,growth.rate_mean))+
	geom_point(shape=21,col='grey70')+
	geom_point(data=singles,aes(x=GR.predicted_mean,y=growth.rate_mean,col=`single mutated\nlocus`),shape=21,size=3,stroke=1)+
	geom_errorbar(data=singles,aes(ymax=growth.rate_upper,ymin=growth.rate_lower,col=`single mutated\nlocus`))+
	geom_errorbarh(data=singles,aes(xmax=xmax,xmin=xmin,col=`single mutated\nlocus`))+
	geom_abline(col='red')+xlab(expression(atop('predicted growth rate '*h^-1,'lm(growth.rate ~ fracon.gal+fracon.glu) ')))+
	ylab(expression(paste('observed growth rate ',h^-1)))
varexplained(summ1$growth.rate_mean,summ1$pred_mean)
w<-6;h<-5.25
ggsave(paste0(figout,'1707xx-180808-pNCDK-lm(gr~fracon.gal+fracon.glu).png'), pNCDK,width=w,height=h)


###################################################
### basic phenotype X-Y plotting
###################################################
DT0<-copy(phenDT)
phenocols<-c('growth.rate','phenotypic.index','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')
phenosh<-c('GR','PI','fONglu','fONgal','ymglu','ymgal')
genocols<-colnames(CloneNamedSummary(DT00[,c('clone.named'),with=F]))
setnames(DT0,phenocols,phenosh)

summ<-summary.func.all(DT0,phenosh,genocols)

phenoplot1<-function(summ,phenoy='GR',phenox='PI',hlight=NULL,plot=T,bycols=genocols,abline=F,facet=T){
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
		geom_point(data=mDT[single_lv==T],aes(x,y,col= single_locus),shape=21,stroke=1,size=2)+
		geom_errorbar(data= mDT[single_lv==T],aes(ymin=y_lower,ymax=y_upper, col= single_locus),width=0)+
		geom_errorbarh(data= mDT[single_lv==T],aes(xmin=x_lower,xmax=x_upper, col= single_locus),height=0)+
		geom_point(data= mDT[WT_lv==T],aes(x,y),col='red',size=3,shape=21,stroke=1)+
		geom_errorbar(data= mDT[WT_lv==T],aes(ymin=y_lower,ymax=y_upper),col='red',width=0)+
		geom_errorbarh(data= mDT[WT_lv==T],aes(xmin=x_lower,xmax=x_upper),col='red',height=0)+
		xlab(phenox)+ylab(phenoy)
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

sample(letters,4)
pYDOP<-phenoplot1(summ)
pFZMH<-phenoplot1(summ,phenoy='fONgal',phenox='fONglu',abline=T)
pYDOPa<-phenoplot1(summ,facet=F)
pFZMH<-phenoplot1(summ,phenoy='fONgal',phenox='fONglu',abline=T)
pFZMHa<-phenoplot1(summ,phenoy='fONgal',phenox='fONglu',abline=T,facet=F)

w<-10;h<-4
ggsave(paste0(figout,'1707xx-180710-pYDOP-PIvsGR-LibraryWithSinglesHighlighted.png'), pYDOP,width=w,height=h)
w<-5;h<-4
ggsave(paste0(figout,'1707xx-180710-pYDOPa-PIvsGR-NotFacetted-LibraryWithSinglesHighlighted.png'), pYDOPa,width=w,height=h)
w<-10;h<-4
ggsave(paste0(figout,'1707xx-180710-pFZMH-FONgluvsFONgal-LibraryWithSinglesHighlighted.png'), pFZMH,width=w,height=h)
w<-5;h<-4
ggsave(paste0(figout,'1707xx-180710-pFZMHa-FONgluvsFONgal-NotFacetted-LibraryWithSinglesHighlighted.png'), pFZMHa,width=w,height=h)


###################################################
### cluster growth rate across genetic backgrounds
###################################################

clusterfun<-function(DTr,left,right,singleGeno,singleGeno2,nclust,PI){
	form<-as.formula(castform(left,right))
	cDT<-dcast(DTr,form,value.var='lx')
	sDT1<-cDT[,lapply(.SD,scale),.SDcols=c(notin(colnames(cDT),singleGeno))]
	sDT1[,x:=cDT[,singleGeno,with=F]]
	setnames(sDT1,'x',singleGeno)
	lv<-as.logical(paste0(unlist(sDT1[,(lapply(.SD,function(x)as.logical(prod(unique(!is.nan(x))))))])))
	notnacols<-colnames(sDT1)[lv]
	sDT<-sDT1[,notnacols,with=F]
#	tidbit(sDT)
	setcolorder(sDT, singleGeno)
	clusters<-hclust(dist(sDT[,! singleGeno,with=F]))
#	plot(clusters)
	sDT[,arb:=1:nrow(sDT)]
#	str(clusters)
	clust<-cutree(clusters,nclust)
	setnames(PI,singleGeno,'x')
	ordframe <-data.table(arb=1:nrow(sDT),X1=1:nrow(sDT),X2=1:nrow(sDT),x=unlist(sDT[,singleGeno,with=F]),A=unlist(sDT[,singleGeno,with=F]),B=unlist(sDT[,singleGeno,with=F]),ord=clusters$order,cluster1=clust)
	boo<-merge(ordframe,PI,by='x')[order(cluster1,V1)]
	boo[,ord:=1:nrow(boo)]
	distmat<-data.table(melt(as.matrix(cophenetic(clusters))))
	dMT<-merge(merge(distmat,boo[,c('A','X1'),with=F],by='X1'), boo[,c('B','X2'),with=F],by='X2')
	dMT[,A:=factor(A,levels= boo$A)]
	dMT[,B:=factor(B,levels= boo$B)]
	p<-ggplot(dMT,aes(A,B))+geom_tile(aes(fill=value))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'))+theme(axis.text.x=element_text(angle=30,hjust=1))
	setnames(boo,c('x','V1'),c(singleGeno,'PI_mean'))
	return(list(plot=p,
		cluster_order=boo[,c(singleGeno,singleGeno2,'PI_mean','cluster1','ord'),with=F],
		clust_obj=clusters
		))
}


DT0<-copy(phenDT)#copy(densDT)

phenocols<-c('growth.rate')#dens.cols#
genocols<-c(colnames(CloneNamedSummary(DT00[,c('clone.named','allele.named'),with=F])))
subcols<-c(phenocols,'clone.named','aa','bb','cc','aa1','bb1','cc1','allele.named','single_lv','double_lv')
a<-DT0[,c(subcols),with=F]
bycols2<-c('aa1','bb1','cc1','aa','bb','cc','clone.named','allele.named','single_lv','double_lv')
ord_within_cluster<-'growth.rate'
# ord_within_cluster<-'phenotypic.index'


mDT1<-melt(a,id=c(bycols2))
dens.arb<-paste0('X',1:length(phenocols))
DT.dens.arb<-data.table(variable=phenocols,dens.arb=factor(dens.arb,levels=dens.arb))
mDT2<-merge(mDT1,DT.dens.arb,by='variable')

# get expression data together
subcols3<-c('aa1','bb1','cc1','dens.arb')
mDT1<-mDT2[,mean(value,na.rm=T),by=c(subcols3)]
mDT1[V1<0,V1:=0]
mDT1[,lx:=log10(V1+0.001)]
mDT<-mDT1[!is.na(lx)]
nclust<-6

# get singles data across all genetic backgrounds together for GAL4
G4<-mDT[!(aa1=='GAL4.WT')]
aabb<-G4[cc1=='GAL80.WT',!'cc1']
setnames(aabb,'bb1','query')
aacc<-G4[bb1=='GAL3.WT',!'bb1']
setnames(aacc,'cc1','query')
DTr<-rbind(aabb,aacc)
left<-c('aa1')
right<-c('query','dens.arb')
singleGeno<-c('aa1')
singleGeno2<-c('aa')
if(ord_within_cluster=='phenotypic.index')PI<-phenDT[G4.single==T,mean(phenotypic.index,na.rm=T),by=c(singleGeno2, singleGeno)]
if(ord_within_cluster=='growth.rate')PI<-phenDT[G4.single==T,mean(growth.rate,na.rm=T),by=c(singleGeno2, singleGeno)]
test<-clusterfun(DTr,left,right,singleGeno,singleGeno2,nclust,PI)
# test$plot;dev.new();plot(test$clust_obj);test$cluster_order

# get singles data across all genetic backgrounds together for GAL3
G3<-mDT[!(bb1=='GAL3.WT')]
bbcc<-G3[aa1=='GAL4.WT',!'aa1']
setnames(bbcc,'cc1','query')
bbaa<-G3[cc1=='GAL80.WT',!'cc1']
setnames(bbaa,'aa1','query')
DTr<-rbind(bbcc, bbaa)
left<-c('bb1')
right<-c('query','dens.arb')
singleGeno<-c('bb1')
singleGeno2<-c('bb')
if(ord_within_cluster=='phenotypic.index')PI<-phenDT[G3.single==T,mean(phenotypic.index,na.rm=T),by=c(singleGeno2, singleGeno)]
if(ord_within_cluster=='growth.rate')PI<-phenDT[G3.single==T,mean(growth.rate,na.rm=T),by=c(singleGeno2, singleGeno)]
test2<-clusterfun(DTr,left,right,singleGeno,singleGeno2,nclust, PI)
# test2$plot;dev.new();plot(test2$clust_obj);test2$cluster_order

# get singles data across all genetic backgrounds together for GAL80
G80<-mDT[!(cc1=='GAL80.WT')]
ccbb<-G80[aa1=='GAL4.WT',!'aa1']
setnames(ccbb,'bb1','query')
ccaa<-G80[bb1=='GAL3.WT',!'bb1']
setnames(ccaa,'aa1','query')
DTr<-rbind(ccbb, ccaa)
left<-c('cc1')
right<-c('query','dens.arb')
singleGeno<-c('cc1')
singleGeno2<-c('cc')
if(ord_within_cluster=='phenotypic.index')PI<-phenDT[G80.single==T,mean(phenotypic.index,na.rm=T),by=c(singleGeno2, singleGeno)]
if(ord_within_cluster=='growth.rate')PI<-phenDT[G80.single==T,mean(growth.rate,na.rm=T),by=c(singleGeno2, singleGeno)]
test3<-clusterfun(DTr,left,right,singleGeno,singleGeno2,nclust, PI)
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


phenomap<-'growth.rate'
dSub<-phenDT[,c('aa1','bb1','cc1','aa','bb','cc',phenomap),with=F]
dSub[,aa1:=factor(aa1,c('GAL4.WT',G4DT$aa1))]
dSub[,bb1:=factor(bb1,c('GAL3.WT',G3DT$bb1))]
dSub[,cc1:=factor(cc1,c('GAL80.WT',G80DT$cc1))]

summ<-dSub[,summary.func(growth.rate),by=c('aa1','bb1','cc1','aa','bb','cc')][order(aa1,bb1,cc1)]
backg_summ<-dSub[aa1=='GAL4.46'|aa1=='GAL4.45',summary.func(growth.rate)]
backg<-backg_summ$mean
backgsd<-backg_summ$sd
summ[mean<backg-backgsd,mean:=backg-sd]
G80G3<-summ[aa1=='GAL4.WT']
G80G3[,c('A','B'):=list(cc1,bb1)]
G3G4<-summ[cc1=='GAL80.WT']
G3G4[,c('A','B'):=list(bb1,aa1)]
G80G4<-summ[bb1=='GAL3.WT']
G80G4[,c('A','B'):=list(cc1,aa1)]

G80G3[,c('Axx','Bxx'):=list(cc,bb)]
G3G4[,c('Axx','Bxx'):=list(bb,aa)]
G80G4[,c('Axx','Bxx'):=list(cc,aa)]

DTr<-rbind(G80G3,G3G4,G80G4)
G80G3x<-data.table(A=c(c('GAL80.WT',G80DT$cc1), c('GAL3.WT',G3DT$bb1)),Ax=factor(1:98,levels=1:98))
G4G3x<-data.table(B=c(c('GAL4.WT',G4DT$aa1), c('GAL3.WT',G3DT$bb1)),Bx= factor(1:98,levels=1:98))
G80G3xx<-data.table(Axx=c(c('GAL80.WT',G80DT$cc), c('GAL3.WT',G3DT$bb)),Ax=factor(1:98,levels=1:98))
G4G3xx<-data.table(Bxx=c(c('GAL4.WT',G4DT$aa), c('GAL3.WT',G3DT$bb)),Bx= factor(1:98,levels=1:98))

DTr[,A:=factor(A,levels= G80G3x$A[98:1])]
DTr[,B:=factor(B,levels= G4G3x$B[98:1])]
DTr[,Axx:=factor(Axx,levels=G80G3xx $Axx[85:1])]
DTr[,Bxx:=factor(Bxx,levels= G4G3xx $Bxx[89:1])]

DTr[,`growth rate, 1/h`:=mean]

collim<-c(backg-backgsd-0.01,0.3) # collim<-c(0,2)# 
sample(letters,4)
pIOZG<-ggplot(DTr,aes(B, A))+geom_tile(aes(fill=`growth rate, 1/h`))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),na.value='white',limits=collim)+
	theme(axis.text.x=element_text(angle=50,hjust=1,size=6))+
	theme(axis.text.y=element_text(size=6))+
	# scale_x_discrete(labels= G4G3xx $Bxx[98:1])+
	# scale_y_discrete(labels= G80G3xx $Axx[98:1],position='right')+xlab('')+ylab('')+
	scale_y_discrete(position='right')+xlab('')+ylab('')+
	theme(plot.margin = unit(c(1,1,1,1), "cm"))
w<-13;h<-12.5
ggsave(paste0(figout,'1707xx-180710-pIOZG-ClusteredAcrossGeneticBackgroudnds-GrowthRate.png'), pIOZG,width=w,height=h)



###################################################
### cluster expression density across genetic backgrounds
###################################################

DT0<-copy(densDT)#copy(phenDT)#

phenocols<-dens.cols#c('growth.rate')#
genocols<-c(colnames(CloneNamedSummary(DT00[,c('clone.named','allele.named'),with=F])))
subcols<-c(phenocols,'clone.named','aa','bb','cc','aa1','bb1','cc1','allele.named','single_lv','double_lv')
a<-DT0[,c(subcols),with=F]
bycols2<-c('aa1','bb1','cc1','aa','bb','cc','clone.named','allele.named','single_lv','double_lv')
# ord_within_cluster<-'growth.rate'
ord_within_cluster<-'phenotypic.index'


mDT1<-melt(a,id=c(bycols2))
dens.arb<-paste0('X',1:length(phenocols))
DT.dens.arb<-data.table(variable=phenocols,dens.arb=factor(dens.arb,levels=dens.arb))
mDT2<-merge(mDT1,DT.dens.arb,by='variable')

# get expression data together
subcols3<-c('aa1','bb1','cc1','dens.arb')
mDT1<-mDT2[,mean(value,na.rm=T),by=c(subcols3)]
mDT1[,lx:=log10(V1+0.001)]
mDT<-mDT1[!is.na(lx)]
nclust<-3

# get singles data across all genetic backgrounds together for GAL4
G4<-mDT[!(aa1=='GAL4.WT')]
aabb<-G4[cc1=='GAL80.WT',!'cc1']
setnames(aabb,'bb1','query')
aacc<-G4[bb1=='GAL3.WT',!'bb1']
setnames(aacc,'cc1','query')
DTr<-rbind(aabb,aacc)
left<-c('aa1')
right<-c('query','dens.arb')
singleGeno<-c('aa1')
singleGeno2<-c('aa')
if(ord_within_cluster=='phenotypic.index')PI<-phenDT[G4.single==T,mean(phenotypic.index,na.rm=T),by=c(singleGeno2, singleGeno)]
if(ord_within_cluster=='growth.rate')PI<-phenDT[G4.single==T,mean(growth.rate,na.rm=T),by=c(singleGeno2, singleGeno)]
test<-clusterfun(DTr,left,right,singleGeno,singleGeno2,nclust,PI)
# test$plot;dev.new();plot(test$clust_obj);test$cluster_order

# get singles data across all genetic backgrounds together for GAL3
G3<-mDT[!(bb1=='GAL3.WT')]
bbcc<-G3[aa1=='GAL4.WT',!'aa1']
setnames(bbcc,'cc1','query')
bbaa<-G3[cc1=='GAL80.WT',!'cc1']
setnames(bbaa,'aa1','query')
DTr<-rbind(bbcc, bbaa)
left<-c('bb1')
right<-c('query','dens.arb')
singleGeno<-c('bb1')
singleGeno2<-c('bb')
if(ord_within_cluster=='phenotypic.index')PI<-phenDT[G3.single==T,mean(phenotypic.index,na.rm=T),by=c(singleGeno2, singleGeno)]
if(ord_within_cluster=='growth.rate')PI<-phenDT[G3.single==T,mean(growth.rate,na.rm=T),by=c(singleGeno2, singleGeno)]
test2<-clusterfun(DTr,left,right,singleGeno,singleGeno2,nclust, PI)
# test2$plot;dev.new();plot(test2$clust_obj);test2$cluster_order

# get singles data across all genetic backgrounds together for GAL80
G80<-mDT[!(cc1=='GAL80.WT')]
ccbb<-G80[aa1=='GAL4.WT',!'aa1']
setnames(ccbb,'bb1','query')
ccaa<-G80[bb1=='GAL3.WT',!'bb1']
setnames(ccaa,'aa1','query')
DTr<-rbind(ccbb, ccaa)
left<-c('cc1')
right<-c('query','dens.arb')
singleGeno<-c('cc1')
singleGeno2<-c('cc')
if(ord_within_cluster=='phenotypic.index')PI<-phenDT[G80.single==T,mean(phenotypic.index,na.rm=T),by=c(singleGeno2, singleGeno)]
if(ord_within_cluster=='growth.rate')PI<-phenDT[G80.single==T,mean(growth.rate,na.rm=T),by=c(singleGeno2, singleGeno)]
test3<-clusterfun(DTr,left,right,singleGeno,singleGeno2,nclust, PI)
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


phenomap<-'phenotypic.index'
dSub<-phenDT[,c('aa1','bb1','cc1',phenomap),with=F]
dSub[,aa1:=factor(aa1,c('GAL4.WT',G4DT$aa1))]
dSub[,bb1:=factor(bb1,c('GAL3.WT',G3DT$bb1))]
dSub[,cc1:=factor(cc1,c('GAL80.WT',G80DT$cc1))]

summ<-dSub[,summary.func(phenotypic.index),by=c('aa1','bb1','cc1')][order(aa1,bb1,cc1)]
# backg_summ<-dSub[aa1=='GAL4.46'|aa1=='GAL4.45',summary.func(phenomap)]
# backg<-backg_summ$mean
# backgsd<-backg_summ$sd
# summ[mean<backg-backgsd,mean:=backg-sd]
G80G3<-summ[aa1=='GAL4.WT']
G80G3[,c('A','B'):=list(cc1,bb1)]
G3G4<-summ[cc1=='GAL80.WT']
G3G4[,c('A','B'):=list(bb1,aa1)]
G80G4<-summ[bb1=='GAL3.WT']
G80G4[,c('A','B'):=list(cc1,aa1)]
DTr<-rbind(G80G3,G3G4,G80G4)
G80G3x<-data.table(A=c(c('GAL80.WT',G80DT$cc1), c('GAL3.WT',G3DT$bb1)),Ax=factor(1:98,levels=1:98))
G4G3x<-data.table(B=c(c('GAL4.WT',G4DT$aa1), c('GAL3.WT',G3DT$bb1)),Bx= factor(1:98,levels=1:98))
DTr[,A:=factor(A,levels= G80G3x$A[98:1])]
DTr[,B:=factor(B,levels= G4G3x$B[98:1])]
dev.new()
collim<-collim<-c(0,2)# c(backg-backgsd-0.01,0.3) # 
# need to reorder clusters to make smoother transitions
ggplot(DTr,aes(B, A))+geom_tile(aes(fill=mean))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),na.value='white',limits=collim)+
	theme(axis.text.x=element_text(angle=40,hjust=1,size=6))+
	theme(axis.text.y=element_text(size=6))+
	scale_y_discrete(position='right')


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
denscolconv<-function(DT, vnull =T){
	DT[,variable:=gsub('dens__','',variable)]
	DT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
	DT[,sugar:=paste0('g',sugar)]
	if(vnull ==T)DT[,variable:=NULL]
}

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
mDT[,`expression density, log10`:=value]
uniqfun<-function(x,fac=1,shift=1)(unique(x)[order(unique(x))][seq(shift,length(unique(x)),by=fac)])
collim<-range(mDT$value)
singM<-mDT[single_lv==T,.N,by=c('Y','dumm','single','single_locus')]
xlab<-c(uniqfun(mDT$AU.FL,fac=4,shift=3),uniqfun(round(mDT$AU.FL,3),shift=3,fac=4))
xbreak<-uniqfun(mDT$X,fac=4,shift=3)
sample(letters,4)
pSZJE<-ggplot(mDT,aes(X,Y))+geom_tile(aes(fill=`expression density, log10`))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),na.value='white',limits=collim)+
	theme(axis.text.x=element_text(angle=40,hjust=1,size=8))+
	theme(axis.text.y=element_text(size=6))+
	geom_vline(xintercept=61,col='black',size=2)+
	geom_vline(xintercept=61,col='red',size=1.5)+
	geom_point(data=data.frame(x=126,y=1),aes(x,y),col='white')+
	geom_point(data=singM,aes(dumm-runif(length(dumm))*2,Y,colour=single_locus),size=1,shape=21)+scale_colour_manual(values=c('orange1','indianred','chartreuse3','blue'),na.value='transparent')+
	scale_x_discrete(breaks=xbreak,labels=xlab)+
	theme(axis.title.y=element_blank(),axis.text.y=element_blank(),	axis.ticks.y=element_blank())+
	xlab('Expression in glucose                            Exression in galactose       \n pseudo log10 A.U. Fluorescence units')+
	theme(axis.title=element_text(size=10,face="plain"))


# # w<-7;h<-9
# ggsave(paste0(figout,'1707xx-180710-pSZJE-ExpressionDensityTilePlot-AllCloneNamed.png'), pSZJE,width=w,height=h)




############
### plot archetypal expression densities: use all clusters
############

mDT1<-mDT[cluster_outlier_score<1]
mDT1[,density:=exp(value)-0.001]
setnames(mDT1,'V1','PhenotypeForRanking')
bycols<-c('cluster','sugar','AU.FL')
summ1<-mDT1[,mean(density),by=bycols]
summ2<-mDT[,.N,by=bycols]
summ2[,prop:=N/length(unique(mDT$clone.named))]
summ<-merge(summ1,summ2,by=bycols)
summ[,facet:=percent(round(prop,2))]
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
sample(LETTERS,4)
pBYGF<-ggplot(mDTm,aes(AU.FL,V2))+geom_point(data=mDTm,aes(AU.FL,density,col=sugar),size=0.1,alpha=0.3,shape=21)+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(cluster_new+facet~ mutant_facet)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density of observations') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))
# w<-14;h<-10
# ggsave(paste0(figout,'1707xx-180720-pBYGF-ExpressionDensityClusters-LinePlot-FacetByMutatedGene.png'), pBYGF,width=w,height=h)


############
### plot archetypal expression densities: combine all const low expr 
############

mDT1<-mDT[cluster_outlier_score<1]
mDT1[,density:=exp(value)-0.001]
setnames(mDT1,'V1','PhenotypeForRanking')
mDT1[cluster%in%c(0),cluster_new:='outlying']
mDT1[cluster%in%c(1),cluster_new:='weakly inducible']
mDT1[cluster%in%c(2),cluster_new:='const, low expr']
mDT1[cluster%in%c(3),cluster_new:='const, low expr']
mDT1[cluster%in%c(4),cluster_new:='const, low expr']
mDT1[cluster%in%c(5),cluster_new:='uninducible']
mDT1[cluster%in%c(6),cluster_new:='inducible']
mDT1[cluster%in%c(7),cluster_new:='constitutive']
bycols1<-c('cluster_new','sugar','AU.FL')
# bycols2<-c('cluster_new','cluster','sugar','AU.FL')
summ<-mDT1[,.N,by= bycols1]
summ[,prop:=N/length(unique(mDT$clone.named))]
summ[,facet:=percent(round(prop,2))]
mDTm<-merge(summ,mDT1,by=bycols1)
dumm<-mDTm[,unique(cluster_new)]
data.table(dumm,1:length(dumm))
cluster_new_vec<-factor(mDTm$cluster_new,levels=dumm[c(5,3,2,6,1,4)])
mDTm[,cluster_new:=cluster_new_vec]
mDTm[,V2:=mean(density),by=bycols1]
mDTm[,dummyline:=paste0('z',sugar)]
mDTm[,mutant_facet:=mutant_gene]
mDTm[mutant_gene=='GAL3 GAL80 GAL4',mutant_facet:='Control']
mDTm[single_lv==T&aa=='GAL4.delta',mutant_facet:='Control']
mDTm[single_lv==T&bb=='GAL3.delta',mutant_facet:='Control']
mDTm[single_lv==T&cc=='GAL80.delta',mutant_facet:='Control']
dumm2a<-mDTm[,unique(mutant_facet),by='mutcount'][order(mutcount,V1)]
dumm2<-dumm2a[!(mutcount==1&V1=='Control')]
mDTm[,mutant_facet:=factor(mutant_facet,levels=dumm2$V1)]
sample(LETTERS,4)
pUNWE<-ggplot(mDTm,aes(AU.FL,V2))+geom_point(data=mDTm,aes(AU.FL,density,col=sugar),size=0.1,alpha=0.3,shape=21)+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(cluster_new+facet~ mutant_facet)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density of observations') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	theme(strip.text.x = element_text(size = 6),strip.text.y = element_text(size = 6))
# w<-6.8;h<-5.7
# ggsave(paste0(figout,'1707xx-180720-pUNWE-ExpressionDensityClusters-AggregatedConstitutiveOutliers-LinePlot-FacetByMutatedGene2.png'), pUNWE,width=w,height=h)

# ## quick look at the weakly inducible class
# table(mDTm[cluster_new=='weakly inducible'&single_lv==T]$clone.named)
# weak_inducible_clone.named<-unique(mDTm[cluster_new=='weakly inducible'&single_lv==T]$clone.named)
# weak_inducible<-paste0(sapply(weak_inducible_clone.named,function(x)strsplit(x,'\\ ')[[1]][3]))
# phenDT[clone.named%in%weak_inducible_clone.named,mean(growth.rate,na.rm=T),by=clone.named]
# qsumm<-phenDT[aa%in%weak_inducible,mean(growth.rate,na.rm=T),by=c('clone.named','bb','cc','aa')]
# qsumm[,val:=V1]
# qsumm[V1<0.04,val:=0.04]
# ggplot(qsumm,aes(bb,cc))+geom_tile(aes(fill=val))+
	# scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),limits=c(0.04,0.14))+facet_wrap(~aa)

############
### use archetypal expression clusters to predict growth rate 
############
singles<-mDTm[single_lv==T,unique(cluster_new),by=c('aa','bb','cc')]
G4<-singles[bb=='GAL3.WT'&cc=='GAL80.WT']
G3<-singles[aa=='GAL4.WT'&cc=='GAL80.WT']
G80<-singles[bb=='GAL3.WT'&aa=='GAL4.WT']
dt.list<-list(G4,G3,G80)
boo<-c('aa','bb','cc')
dts<-lapply(1:3,function(i){
	clustcol<-paste0(boo[i],'_','clust')
	setnames(dt.list[[i]],'V1',clustcol)
	dt.list[[i]][,c(boo[i],clustcol),with=F]
})
clustcols<-unlist(lapply(1:3,function(i){
	clustcol<-paste0(boo[i],'_','clust')
}))

pheno<-'growth.rate'
genos<-c('clone.named','aa','bb','cc')
DTm<-merge_recurse3(qDT=phenDT,dt.list=dts,by.list=boo)
qDT<-DTm[,c(pheno,clustcols,genos),with=F]
setnames(qDT,pheno,'x')
if(pheno=='growth.rate'){
	credMin<-qDT[aa=='GAL4.delta',mean(x,na.rm=T)-sd(x,na.rm=T)]
	qDT[x<credMin,x:=credMin]
}
pred<-qDT[,summary.func(x,'pred'),by=clustcols][order(pred_N,decreasing=T)]
genosumm<-qDT[,summary.func(x,'obs'),by=genos]
DTm2<-merge(merge(qDT,pred,by=clustcols),genosumm,by=genos)
varexpA<-DTm2[,varexplained(obs_mean,pred_mean)]
varexpB<-DTm2[,varexplained(x,pred_mean)]
mod_genotypic_ensemble<-lm(obs_mean~pred_mean,DTm2) # this comes into play below
p0<-ggplot(DTm2[obs_CoV<0.5],aes(pred_mean,obs_mean,xmin=pred_lower,xmax=pred_upper,ymin=obs_lower,ymax=obs_upper))+geom_errorbar(col='grey70',alpha=0.2)+geom_errorbarh(col='grey70',alpha=0.2)+
	geom_point(alpha=0.1)+geom_abline()
p1<-ggplot(DTm2[obs_CoV<0.5],aes(pred_mean,obs_mean,xmin=pred_lower,xmax=pred_upper,ymin=obs_lower,ymax=obs_upper))+
	geom_point(alpha=0.1)+geom_abline()

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

pred
DT2<-merge(DT1, mDTm[,unique(cluster_new),by='clone.named'],by='clone.named')
setnames(DT2,'V1','cluster_new')

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
pDYHQ<-ggplot(DTr,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))+geom_rect(aes(fill= `gene expression cluster`),col='grey40')+scale_fill_manual(values=c(rainbow(5),'grey70'),na.value='white')+
	theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
	theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),panel.border=element_blank(),axis.line=element_blank())+
	facet_grid(~var)
# w<-7;h<-4.5
# ggsave(paste0(figout,'1707xx-180724-pDYHQ-GeneExpressionClusterHierarchy-Rectangles.png'), pDYHQ,width=w,height=h)

################
### constitutive GAL4 + unindubible GAL80 yield a lot of uncharacterized profiles. look closer
################
DT2<-merge(DT1, mDTm[,unique(cluster_new),by='clone.named'],by='clone.named')
setnames(DT2,'V1','cluster_new')

clones1<-DT2[aa_clust=='constitutive'&cc_clust=='uninducible',c('aa','cc')]
ctrls<-data.table(aa=c('GAL4.WT','GAL4.delta'),cc=c('GAL80.WT','GAL80.delta'))
clones2<-rbind(clones1,ctrls)
clones3<-data.table(expand.grid(clones2))
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
DT2<-merge(DT1, mDTm[,unique(cluster_new),by='clone.named'],by='clone.named')
phenDT0<-merge(DT2[,c('clone.named',grepincols(DT2,'clust')),with=F],phenDT,by='clone.named')
setnames(DT2,'V1','cluster_new')

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


################
### get GAL80 and GAL3 functional scores and plot in 2d
################
DT2<-merge(DT1, mDTm[,unique(cluster_new),by='clone.named'],by='clone.named')
phenDT0<-merge(DT2[,c('clone.named',grepincols(DT2,'clust')),with=F],phenDT,by='clone.named')

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


######
### sensitized background plotting
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
ggplot(tog,aes(WT_mean,MUT_mean,xmin=WT_lower,xmax=WT_upper,ymin=MUT_lower,ymax=MUT_upper))+geom_errorbarh()+geom_errorbar()+geom_point()+
	ylab('yfp sig in glu\nGAL80 variants in weakly inducible GAL4 background')+xlab('yfp sig in glu\nGAL80 variants in WT GAL4 background')


WT<-clonesub[aa_clust=='inducible',summary.func(log10(yfp.mean.glu),'WT'),by=c('cc')]
MUT<-clonesub[aa%in%c('GAL4.01','GAL4.03','GAL4.06'),summary.func(log10(yfp.mean.glu),'MUT'),by=c('cc')]
tog<-merge(WT,MUT,by=c('cc'))
yfplogTog<-copy(tog)
ggplot(tog,aes(WT_mean,MUT_mean,xmin=WT_lower,xmax=WT_upper,ymin=MUT_lower,ymax=MUT_upper))+geom_errorbarh()+geom_errorbar()+geom_point()+
	ylab('log yfp glu\nGAL80 variants in weakly inducible GAL4 background')+xlab('log yfp glu\nGAL80 variants in WT GAL4 background')


######
### characterize GAL80 in glucose and GAL3 in galactose in GAL4 inducible backgrounds
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


# PCA on GAL3
DTm<-merge(G3_wt, G3_mut,by='bb',allow.cartesian=T)
PCAmod<-PCA(DTm[,!'bb',with=F])
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



# do above a different way
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
ggplot(sDTm,aes(mean,GR_mean))+geom_point()+facet_grid(GAL3_cat~pheno)



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

################
### plot growth rate as function of GAL80 facetted across GAL3 and GAL4 functional categories
################

# # This code is greyed out but it works fine. I prefer the plots above

# ########
# ### get GAL80 function single values
# ########
# clonesub1<-phenDT0[(aa_clust=='inducible'|aa%in%c('GAL4.01','GAL4.03','GAL4.06'))&bb_clust=='inducible'&cc_clust!='uninducible'&cc_clust!='const, low expr']
# clonesub1[,outlier:=clone.named%in%c('GAL3.WT GAL80.01 GAL4.01','GAL3.WT GAL80.02 GAL4.01')]
# clonesub<-clonesub1[outlier==F]
# clonesub[,lyfpglu:=log10(yfp.mean.glu)]
# G80<-copy(clonesub)
# setnames(G80,c('fracon.glu','yfp.mean.glu','lyfpglu'),c('G80.fON','G80.mean.yfp','G80.lyfp'))
# G80[,G80.logit:=logit(G80.fON+0.001)]
# subG80<-G80[,c('aa','bb','cc','clone.named','growth.rate',grepincols(G80,'G80.')),with=F][,!c('G80_mut_lv','G80.single'),with=F]
# sumG80<-summary.func.all(subG80,grepincols(subG80,c('80')),bylist=notin(colnames(subG80),c('growth.rate','aa','bb','clone.named','80')))


# # sumG80<-summary.func.all(subG80,grepincols(subG80,c('growth.rate','80')),bylist=notin(colnames(subG80),c('growth.rate','80')))
# # ggplot(sumG80,aes(G80.logit_mean))+geom_histogram()
# # ggplot(sumG80,aes(G80.logit_mean,growth.rate_mean))+geom_point()


# ########
# ### categorize GAL4 and GAL3 "singles" by their behavior in WT backgrounds
# ########
# phenDT[,logit.fon.gal:=logit(fracon.gal+0.001)]
# hist(phenDT[single_lv==T& mutant_gene=='GAL4',fracon.gal])
# hist(phenDT[single_lv==T& mutant_gene=='GAL4',logit.fon.gal])
# G3tab<-table(cut(phenDT[single_lv==T& mutant_gene=='GAL3',fracon.gal],4))
# # G4tab<-table(cut(phenDT[single_lv==T& mutant_gene=='GAL4', logit.fon.gal],4))
# # G3tab<-table(cut(phenDT[single_lv==T& mutant_gene=='GAL3', logit.fon.gal],4))

# G3cut<-strsplit(gsub(']','',gsub('\\(','',names(G3tab))),',')

# sub<-phenDT0[cc=='GAL80.WT'&bb_clust=='inducible',c('aa','bb','cc','fracon.gal'),with=F]
# summ<-sub[,summary.func(fracon.gal),by=c('aa')]
# G4tab<-table(cut(summ[,mean],4))
# G4cut<-strsplit(gsub(']','',gsub('\\(','',names(G4tab))),',')
# te<-as.list(as.numeric(lapply(unlist(G4cut),function(x)x)))
# summ[,c('A1','A2','B1','B2','C1','C2','D1','D2'):=te]
# x<-summ[1,c('A1','A2','B1','B2','C1','C2','D1','D2','mean')]
# attach(summ);boo<-data.table(t(apply(data.frame(A1,A2,B1,B2,C1,C2,D1,D2,mean),1,function(x){
	# x<-sapply(x,function(x)x)
	# fon<-round(x[9],3)
	# A<-0<=fon&&fon<x[2]
	# B<-x[3]<=fon&&fon<x[4]
	# C<-x[5]<=fon&&fon<x[6]
	# D<-x[7]<=fon&&fon<=1
	# c(A,B,C,D)
# })));detach(summ)
# summ[,c('low','mid-low','mid-high','high'):=boo]
# out<-summ[,c('aa','low','mid-low','mid-high','high'),with=F]
# mout<-melt(out,c('aa'))
# G4_cat<-mout[value==T,!'value',with=F]


# sub<-phenDT0[cc=='GAL80.WT'&aa_clust=='inducible',c('aa','bb','cc','fracon.gal'),with=F]
# summ<-sub[,summary.func(logit(fracon.gal+0.001)),by=c('bb')]
# G3tab<-table(cut(summ1[,mean],5))
# G3cut<-strsplit(gsub(']','',gsub('\\(','',names(G3tab[G3tab>0]))),',')
# te<-as.list(as.numeric(lapply(unlist(G3cut),function(x)x)))
# summ[,c('A1','A2','B1','B2','C1','C2','D1','D2'):=te]
# x<-summ[1,c('A1','A2','B1','B2','C1','C2','D1','D2','mean')]
# attach(summ);boo<-data.table(t(apply(data.frame(A1,A2,B1,B2,C1,C2,D1,D2,mean),1,function(x){
	# x<-sapply(x,function(x)x)
	# fon<-round(x[9],3)
	# A<--5<=fon&&fon<x[2]
	# B<-x[3]<=fon&&fon<x[4]
	# C<-x[5]<=fon&&fon<x[6]
	# D<-x[7]<=fon&&fon<=5
	# c(A,B,C,D)
# })));detach(summ)

# summ[,c('low','mid-low','mid-high','high'):=boo]
# out<-summ[,c('bb','low','mid-low','mid-high','high'),with=F]
# mout<-melt(out,c('bb'))
# G3_cat<-mout[value==T,!'value',with=F]

# setnames(G4_cat,'variable','aa_clust2')
# setnames(G3_cat,'variable','bb_clust2')
# grSum<-phenDT0[aa_clust!='constitutive',summary.func(growth.rate,'GR'),by=c('aa','bb','cc')]
# DTm<-merge(merge(merge(grSum, G4_cat,by='aa'),G3_cat,by='bb'), sumG80,by='cc')
# ggplot(DTm,aes(G80.logit_mean,GR_mean))+geom_point()+
	# facet_grid(aa_clust2~bb_clust2)+theme_minimal()+geom_hline(yintercept=0.04)


# DTm2<-merge(merge(merge(phenDT0[aa_clust!='constitutive',c('aa','bb','cc','growth.rate'),with=F], G4_cat,by='aa'),G3_cat,by='bb'), sumG80,by='cc')
# sDT<-DTm2[,summary.func(growth.rate,'GR'),by=c(grepincols(DTm,c('G80.logit','clust')),'aa','cc')]
# ggplot(sDT,aes(G80.logit_mean,GR_mean))+geom_point()+
	# facet_grid(aa_clust2~bb_clust2)+theme_minimal()+geom_hline(yintercept=0.04)

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

tricolor<-c('cornflowerblue','cornflowerblue','yellow','indianred')
pRZWG<-ggplot(DTr,aes(B, A))+geom_tile(aes(fill=growth.rate_mean))+scale_fill_gradientn(colours= tricolor,na.value='grey70')+
	theme(axis.text.x=element_text(angle=50,hjust=1,size=6))+
	theme(axis.text.y=element_text(size=6))+
	scale_y_discrete(position='right')+xlab('')+ylab('')+
	theme(plot.margin = unit(c(1,1,1,1), "cm"))
bicolor<-c('cornflowerblue','indianred')
pRZWGb<-ggplot(DTr,aes(B, A))+geom_tile(aes(fill=growth.rate_mean))+scale_fill_gradientn(colours= bicolor,na.value='grey70')+
	theme(axis.text.x=element_text(angle=50,hjust=1,size=6))+
	theme(axis.text.y=element_text(size=6))+
	scale_y_discrete(position='right')+xlab('')+ylab('')+
	theme(plot.margin = unit(c(1,1,1,1), "cm"))

w<-13;h<-12.5
ggsave(paste0(figout,'1707xx-180710-pRTBU-ClusteredBy5PhenotypesAcrossGeneticBackgrounds-GrowthRateDeviationFromGenotypicEnsemble.png'), pRTBU,width=w,height=h)
w<-13;h<-12.5
ggsave(paste0(figout,'1707xx-180710-pRZWG-ClusteredBy5PhenotypesAcrossGeneticBackgrounds-GrowthRate.png'), pRZWG,width=w,height=h)






####################################
#### first try at linear regression
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


p<-ggplot(mDT1[test==T],aes(fit_mean,x))+geom_errorbar(col='grey70',alpha=0.3)+geom_errorbarh(col='grey70',alpha=0.3)+geom_point(col='grey70',alpha=0.3)+
	geom_point(data= mDT1[single_lv==T],aes(PI_mean,GR_mean),col='red')+
	geom_point(data= mDT1[WT_lv==T],aes(PI_mean,GR_mean),col='blue')+
	geom_point(data= mDT1[WT==T],aes(PI_mean,GR_mean),col='red',size=3,shape=21,stroke=1)+facet_grid(~ pairing)+theme_light()

DT0<-copy(phenDT)
phenocols<-c('growth.rate','phenotypic.index','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')
phenosh<-c('GR','PI','fONglu','fONgal','ymglu','ymgal')
genocols<-colnames(CloneNamedSummary(DT00[,c('clone.named','allele.named'),with=F]))
setnames(DT0,phenocols,phenosh)

summ1<-summary.func.all(DT0,phenosh,c('clone.named','allele.named'))
summ2<-CloneNamedSummary(summ1)
summ<-merge(summ2,predDT,by=c('aa','bb','cc'))

phenoplot2(summ[GR_CoV<0.5],phenoy='GR',phenox='fit',bycols=genocols,abline=T,facet=T,singles=F)



####################################
#### try multiplicative model
####################################

DT0<-copy(phenDT)


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

DTm1<-merge(DTm,DTm[,t.test2(pheno_mean, mult_pred_mean,pheno_se, mult_pred_se,pheno_N, mult_pred_N,list=T),by='clone.named'],by='clone.named')
DTm1[,p.adj:=p.adjust(p.value,'fdr')]
sigcutoff<-0.05
diffcutoff<-0.025
DTm1[,sig:=p.adj<sigcutoff&abs(diff.of.means)> diffcutoff]

DTm2<-CloneNamedSummary(DTm1[,!c('aa','bb','cc')],allele=F)

p1<-phenoplot2(DTm2[pheno_CoV<0.5& single_lv!=T],phenoy='pheno',phenox='mult_pred',bycols=genocols,abline=T,facet=T,singles=F)
plotDT<-phenoplot2(DTm2[pheno_CoV<0.5& single_lv!=T],phenoy='pheno',phenox='mult_pred',bycols=genocols,abline=T,facet=T,singles=F,plot=F)

	
p0<-ggplot(plotDT[order(sig)],aes(x,y,xmax=x_upper,xmin=x_lower,ymax=y_upper,ymin=y_lower))+geom_errorbar(col='grey50',alpha=0.2)+geom_errorbarh(col='grey50',alpha=0.2)+geom_point(aes(col=sig),alpha=0.3)+
	geom_abline(col='indianred')+scale_colour_manual(values=c('orange1','cornflowerblue'))+
	facet_grid(~ pairing)+theme_light()
p0

p1<-ggplot(plotDT[order(sig)],aes(x,y))+geom_point(aes(col=sig),alpha=0.3)+
	geom_abline(col='indianred')+scale_colour_manual(values=c('orange1','cornflowerblue'))+
	facet_grid(~ pairing)+theme_light()
p1

do.call(grid.arrange,list(p1,p0,ncol=1))

varexplained(DTm2$PI_mean,DTm2$mult_pred_mean)



















function(summ,phenoy='GR',phenox='PI',hlight=NULL,plot=T,bycols=genocols,abline=F,facet=T){
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
		geom_point(data=mDT[single_lv==T],aes(x,y,col= single_locus),shape=21,stroke=1,size=2)+
		geom_errorbar(data= mDT[single_lv==T],aes(ymin=y_lower,ymax=y_upper, col= single_locus),width=0)+
		geom_errorbarh(data= mDT[single_lv==T],aes(xmin=x_lower,xmax=x_upper, col= single_locus),height=0)+
		geom_point(data= mDT[WT_lv==T],aes(x,y),col='red',size=3,shape=21,stroke=1)+
		geom_errorbar(data= mDT[WT_lv==T],aes(ymin=y_lower,ymax=y_upper),col='red',width=0)+
		geom_errorbarh(data= mDT[WT_lv==T],aes(xmin=x_lower,xmax=x_upper),col='red',height=0)+
		xlab(phenox)+ylab(phenoy)
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











sample(letters,4)
pYDOP<-phenoplot1(summ)
pFZMH<-phenoplot1(summ,phenoy='fONgal',phenox='fONglu',abline=T)
pYDOPa<-phenoplot1(summ,facet=F)
pFZMH<-phenoplot1(summ,phenoy='fONgal',phenox='fONglu',abline=T)
pFZMHa<-phenoplot1(summ,phenoy='fONgal',phenox='fONglu',abline=T,facet=F)

w<-10;h<-4
ggsave(paste0(figout,'1707xx-180710-pYDOP-PIvsGR-LibraryWithSinglesHighlighted.png'), pYDOP,width=w,height=h)
w<-5;h<-4
ggsave(paste0(figout,'1707xx-180710-pYDOPa-PIvsGR-NotFacetted-LibraryWithSinglesHighlighted.png'), pYDOPa,width=w,height=h)
w<-10;h<-4
ggsave(paste0(figout,'1707xx-180710-pFZMH-FONgluvsFONgal-LibraryWithSinglesHighlighted.png'), pFZMH,width=w,height=h)
w<-5;h<-4
ggsave(paste0(figout,'1707xx-180710-pFZMHa-FONgluvsFONgal-NotFacetted-LibraryWithSinglesHighlighted.png'), pFZMHa,width=w,height=h)






























########
### GAL80S-V352E continues to be interesting: low rates of growth for ∆GAL3 and 2 other alleles
########

# there are three GAL3 variants taht activate transcription less than the rest in the GAL80S-V352E background
mDT<-clone.plot3(DT0,hlight='GAL80S-V352E',F)
mDT[PI_mean>1&GR_mean<0.1&single_lv==T]
hlight<-'GAL80S-V352E'
mDT<-clone.plot3(DT0,hlight=hlight,F)
weird<-DT0[clone.named%in%mDT[double_lv==T&highlight==T&PI_mean<1& G3_mut_lv==T]$clone.named]$allele.named
weirdplates<-unique(DT0[allele.named%in%weird]$source.plate)
DTemp<-copy(DT0[source.plate%in%weirdplates])
DTemp[allele.named%in%weird,highlight:=T]
plateplot(DTemp[source.plate==weirdplates[1]],'phenotypic.index')+geom_point(aes(col=highlight))
plateplot(DTemp[source.plate==weirdplates[2]],'phenotypic.index')+geom_point(aes(col=highlight))
DTemp[highlight==T]
weirdG3<-unique(DT0[clone.named%in%mDT[double_lv==T&highlight==T&PI_mean<1& G3_mut_lv==T]$clone.named]$bb)
mDT[single_lv==T&bb%in%weirdG3]

########
### basic allele-phenotype plotting 
########

DTemp<-copy(DT0)
DTemp[GAL80.allele=='GAL80.36',GAL80.clone.named:='GAL80S-2']
DTemp[,G80fac:=applyPaste(data.frame(GAL80.clone.named,GAL80.allele),' ')]
DTemp[,G4fac:=applyPaste(data.frame(GAL4.clone.named,GAL4.allele),' ')]
ggplot(DTemp[GAL3.clone.named=='GAL3.WT'],aes(GAL4.clone.named,G80fac))+geom_tile(aes(fill=phenotypic.index))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'))+theme(axis.text.x=element_text(angle=30,hjust=1))
ggplot(DTemp[GAL4.clone.named=='GAL4.WT'],aes(GAL3.clone.named,G80fac))+geom_tile(aes(fill=phenotypic.index))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'))+theme(axis.text.x=element_text(angle=30,hjust=1))
ggplot(DTemp[GAL80.clone.named=='GAL80.WT'],aes(GAL3.clone.named,G4fac))+geom_tile(aes(fill=phenotypic.index))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'))+theme(axis.text.x=element_text(angle=30,hjust=1))


PI<-DT0[G80.single==T|G4.single==T|G3.single==T,summary.func(phenotypic.index),by=c('clone.named','aa','bb','cc','G80.single','G4.single','G3.single')]
mDT1<-melt(PI,id=notcols(PI,'sing'),variable.name='single')
mDT2<-melt(mDT1,id=notcols(mDT1,c('aa','bb','cc'),grepl=F),value.name='allele')
mDT2[value==T]
mDT<-mDT1[value==T]
lDT2<-subsetList(lDT= split(mDT,by='single'), commonCols='mean',return.list=list('cc','aa','bb'),newnames=T)

wc2<-merge_recurse3(weirdclones,lDT2, list('cc','aa','bb'),all_lv=F)
wc2[order(aa)]

# Plate_055 D1 and D2 are failed transformations
plateplot(DT0[source.plate=='Plate_055'],'phenotypic.index')
DT0<-DT0[!(source.plate=='Plate_055'&(rc=='D1'|rc=='D2'))]



lapply(lDT,function(x){
	x[,]
})
merge(merge(merge(weirdclones, mDT[single=='G4.single',c('aa','mean'),with=F],by='aa'), mDT[single=='G3.single',c('bb','mean'),with=F],by='bb'), mDT[single=='G80.single',c('cc','mean'),with=F],by='cc')


DTm<-merge(weirdclones,PIsing)

DT<-mDT1
x<-c('aa','bb','cc')



####################################
#### first try at linear regression
####################################

backg<-mean(DT0[grepl('GAL4.delta',clone.named)]$growth.rate,na.rm=T)
backsd<-sd(DT0[grepl('GAL4.delta',clone.named)]$growth.rate,na.rm=T)
pheno<-'growth.rate'
genos<-c('aa','bb','cc')
dsub<-DT0[,c(pheno,genos),with=F]
setnames(dsub,pheno,'x')
lDT<-logconv2(dsub,backg=backg,backgsd=backsd,genos=genos)
mod1<-lm(lx~aa+bb+cc,lDT)
predDT<-expMod(dsub,mod1,backg=backg,genos=genos)
DTm<-merge(dsub,predDT,by=genos)
DTm[,clone.named:=applyPaste(data.frame(bb,cc,aa),' ')]
DTm2<-CloneNamedSummary(DTm[,!genos,with=F])
mDT1<-melt(DTm2,notcols(DTm2,'vG'),variable.name='pairing',value.name='boo')
mDT1[pairing=='G3vG80_lv',test:=apply(data.frame(G3_mut_lv,G80_mut_lv,G4_mut_lv,WT_lv),1,function(x)as.logical(((x[1]|x[2])&!x[3])|x[4]))]
mDT1[pairing=='G4vG80_lv',test:=apply(data.frame(G3_mut_lv,G80_mut_lv,G4_mut_lv, WT_lv),1,function(x)as.logical((!x[1]&(x[2]|x[3])|x[4])))]
mDT1[pairing=='G3vG4_lv',test:=apply(data.frame(G3_mut_lv,G80_mut_lv,G4_mut_lv, WT_lv),1,function(x)as.logical(((x[1]|x[3])&!x[2])|x[4]))]
mDT1[,highlight:=NA]
if(!is.null(highlight))mDT1[,highlight:=unlist(grepl(hlight,clone.named))]


p<-ggplot(mDT[test==T],aes(fit_mean,x))+geom_errorbar(col='grey70',alpha=0.3)+geom_errorbarh(col='grey70',alpha=0.3)+geom_point(col='grey70',alpha=0.3)+
	geom_point(data=mDT[single_lv==T],aes(PI_mean,GR_mean),col='red')+
	geom_point(data=mDT[WT_lv==T],aes(PI_mean,GR_mean),col='blue')+
	geom_point(data= mDT[WT==T],aes(PI_mean,GR_mean),col='red',size=3,shape=21,stroke=1)+facet_grid(~ pairing)+theme_light()






DT0<-copy(phenDT)
phenocols<-c('growth.rate','phenotypic.index','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')
phenosh<-c('GR','PI','fONglu','fONgal','ymglu','ymgal')
genocols<-colnames(CloneNamedSummary(DT00[,c('clone.named'),with=F]))
setnames(DT0,phenocols,phenosh)

summ<-summary.func.all(DT0,phenosh,genocols)

phenoplot1<-function(summ,phenoy='GR',phenox='PI',hlight=NULL,plot=T,bycols=genocols,abline=F,facet=T){
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
		geom_point(data=mDT[single_lv==T],aes(x,y,col= single_locus),shape=21,stroke=1,size=2)+
		geom_errorbar(data= mDT[single_lv==T],aes(ymin=y_lower,ymax=y_upper, col= single_locus),width=0)+
		geom_errorbarh(data= mDT[single_lv==T],aes(xmin=x_lower,xmax=x_upper, col= single_locus),height=0)+
		geom_point(data= mDT[WT_lv==T],aes(x,y),col='red',size=3,shape=21,stroke=1)+
		geom_errorbar(data= mDT[WT_lv==T],aes(ymin=y_lower,ymax=y_upper),col='red',width=0)+
		geom_errorbarh(data= mDT[WT_lv==T],aes(xmin=x_lower,xmax=x_upper),col='red',height=0)+
		xlab(phenox)+ylab(phenoy)
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

sample(letters,4)
pYDOP<-phenoplot1(summ)
pFZMH<-phenoplot1(summ,phenoy='fONgal',phenox='fONglu',abline=T)
pYDOPa<-phenoplot1(summ,facet=F)
pFZMH<-phenoplot1(summ,phenoy='fONgal',phenox='fONglu',abline=T)
pFZMHa<-phenoplot1(summ,phenoy='fONgal',phenox='fONglu',abline=T,facet=F)

w<-10;h<-4
ggsave(paste0(figout,'1707xx-180710-pYDOP-PIvsGR-LibraryWithSinglesHighlighted.png'), pYDOP,width=w,height=h)
w<-5;h<-4
ggsave(paste0(figout,'1707xx-180710-pYDOPa-PIvsGR-NotFacetted-LibraryWithSinglesHighlighted.png'), pYDOPa,width=w,height=h)
w<-10;h<-4
ggsave(paste0(figout,'1707xx-180710-pFZMH-FONgluvsFONgal-LibraryWithSinglesHighlighted.png'), pFZMH,width=w,height=h)
w<-5;h<-4
ggsave(paste0(figout,'1707xx-180710-pFZMHa-FONgluvsFONgal-NotFacetted-LibraryWithSinglesHighlighted.png'), pFZMHa,width=w,height=h)



