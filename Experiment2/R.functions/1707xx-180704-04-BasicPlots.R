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
# phenosh<-c('GR','PI','fONglu','fONgal','ymglu','ymgal')
factorDefault<-c('orange1','cornflowerblue','darkgreen','red','chartreuse3')

dens.cols<-grepincols(DT00,'dens')
phenDT<-CloneNamedSummary(DT00[,c('clone.named','allele.named',phenocols),with=F])
densDT<-CloneNamedSummary(DT00[,c('clone.named','allele.named',dens.cols),with=F])
genocols<-colnames(CloneNamedSummary(DT00[,c('clone.named','allele.named'),with=F]))

###################################################
### technical overview
###################################################
phenDT<-CloneNamedSummary(DT00[,c('clone.named','allele.named',phenocols),with=F])
DT0<-copy(phenDT)
phenosh<-phenocols#c('GR','PI','frac ON glu','frac ON gal','mean YFP glu','mean YFP gal')
setnames(DT0,phenocols,phenosh)
# use this for factors: scale_colour_manual(values=c('orange1','cornflowerblue','chartreuse3','blue'))
sample(letters,4)
# pTechrep2(DT0,'clone.named',phenosh=phenosh,plot=T)
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
### basic phenotype X-Y plotting: fracon glu vs fracon gal. Newer code than that below that matches aesthetics from first experiment.
###################################################
DT0<-copy(phenDT)
phenocols<-c('growth.rate','phenotypic.index','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')
phenosh<-c('GR','PI','fONglu','fONgal','ymglu','ymgal')
genocols<-colnames(CloneNamedSummary(DT00[,c('clone.named'),with=F]))


summ<-summary.func.all(DT0,phenocols,notin(genocols,c('allele.named','aa1','bb1','cc1')))
summ[,fac:=single_locus]
summ[is.na(single_locus),fac:='double mutant']
summ[,`Single mutant\nlocus`:=factor(fac,levels=c('GAL4','GAL3','GAL80','double mutant'))]
cols<-c(factorDefault[1:3],'grey70')
limits<-c(0,1)
pOCES<-ggplot(summ[order(`Single mutant\nlocus`,decreasing=T)],aes(x=fracon.glu_mean,y=fracon.gal_mean,xmax=fracon.glu_upper,xmin=fracon.glu_lower,ymax=fracon.gal_upper,ymin=fracon.gal_upper))+theme_minimal()+
	geom_abline(col='grey70')+
	geom_errorbarh(aes(col=`Single mutant\nlocus`))+geom_errorbar(aes(col=`Single mutant\nlocus`))+geom_point(aes(col=`Single mutant\nlocus`),shape=21,size=3,stroke=1)+
	scale_colour_manual(values=alpha(cols,c(1,1,1,0.3)))+
	theme(strip.text.x=element_blank(),strip.text.y=element_blank())+
	ylim(limits)+xlim(limits)+
	xlab('fraction ON in glucose')+ylab('fraction ON in galactose')
if(writeplot==T){
	w<-3.8;h<-2.2
	ggsave(paste0(figout,date,'-181010-pOCES-FraconGluVsFraconGal.png'), pOCES,width=w,height=h)	
}



###################################################
### basic phenotype X-Y plotting. Old, not useful for now.
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

# A<-data.table(allele=G4DT$aa1,locus='GAL4')
# B<-data.table(allele=G3DT$bb1,locus='GAL3')
# C<-data.table(allele=G80DT$cc1,locus='GAL80')

# GR<-copy(rbind(A,B,C))
# GR[,phenotype:='growth rate rank']
A<-data.table(allele=G4DT$aa1,locus='GAL4',GR_rank=1:length(G4DT$aa1))
B<-data.table(allele=G3DT$bb1,locus='GAL3',GR_rank=1:length(G3DT$bb1))
C<-data.table(allele=G80DT$cc1,locus='GAL80',GR_rank=1:length(G80DT$cc1))

GR <-copy(rbind(A,B,C))


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

sample(letters,4)
pCOUK<-ggplot(DTr,aes(B, A))+geom_tile(aes(fill=mean))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),na.value='white',limits=collim)+
	theme(axis.text.x=element_text(angle=40,hjust=1,size=6))+
	theme(axis.text.y=element_text(size=6))+
	scale_y_discrete(position='right')

w<-13;h<-12.5
ggsave(paste0(figout,'1707xx-180710-pCOUK-ClusteredAcrossGeneticBackgroudnds-PI.png'), pCOUK,width=w,height=h)

A<-data.table(allele=G4DT$aa1,locus='GAL4',PI_rank=1:length(G4DT$aa1))
B<-data.table(allele=G3DT$bb1,locus='GAL3',PI_rank=1:length(G3DT$bb1))
C<-data.table(allele=G80DT$cc1,locus='GAL80',PI_rank=1:length(G80DT$cc1))

PI_rank<-copy(rbind(A,B,C))



###################################################
### compare ranks
###################################################

DTm<-merge(GR,PI_rank,by=c('allele','locus'))
txt<-DTm[,round(varexplained(GR_rank,PI_rank),2),by=locus]
txt[,c('x','y'):=list(12,40)]
txt[,label:=paste0('Spear R2: ',V1)]
sample(letters,4)
pIMTR<-ggplot(DTm,aes(PI_rank,GR_rank))+geom_point(aes(col=locus))+facet_grid(~locus)+
	geom_abline(col='grey70')+geom_text(data=txt,aes(x=x,y=y,label=label))+
	xlab('phenotypic index rank for allele')+ylab('growth rate rank for allele')
w<-8.5;h<-2.7
ggsave(paste0(figout,'1707xx-180710-pIMTR-ClusterRankCorrelations.png'), pIMTR,width=w,height=h)


