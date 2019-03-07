

# Scripts for analysis.
source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'
source.code.path.file2<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis/R.functions/1707xx-180704-ReanalysisCustomScripts.R'
head.dir<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis'
layout.dir<-'./layout.sample.data'
date<-'1707xx' # beginning date clones were measured
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/publication/'
outputdata<-'./output.data/'
setwd(head.dir)
source(source.code.path.file)
source(source.code.path.file2)
aesthetics.dir<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/181119-Aesthetics_Harmonious_Combinations.R'
source(aesthetics.dir)
writeplot<-F
rm(letters)
br<-seq(-0.5,2.7,length=61)
library(tidyr)

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
# phenosh<-c('GR','PI','fONglu','fONgal','ymglu','ymgal')
factorDefault<-c('orange1','cornflowerblue','darkgreen','red','chartreuse3')

dens.cols<-grepincols(DT00,'dens')
phenDT<-CloneNamedSummary(DT00[,c('clone.named','allele.named',phenocols),with=F])
densDT<-CloneNamedSummary(DT00[,c('clone.named','allele.named',dens.cols),with=F])
genocols<-colnames(CloneNamedSummary(DT00[,c('clone.named','allele.named'),with=F]))



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
summ[,`Single mutant locus`:=factor(fac,levels=c('GAL4','GAL3','GAL80','double mutant'))]
cols<-c(brewer.pal(5, "Set1")[1:3],'grey70')

limits<-c(0,1)
sample(LETTERS,4)
pOKJG<-ggplot(summ[order(`Single mutant\nlocus`,decreasing=T)],aes(x=fracon.glu_mean,y=fracon.gal_mean,xmax=fracon.glu_upper,xmin=fracon.glu_lower,ymax=fracon.gal_upper,ymin=fracon.gal_upper))+
	theme_NewPub+
	geom_abline(col='grey70')+
	geom_errorbarh(aes(col=`Single mutant\nlocus`))+geom_errorbar(aes(col=`Single mutant\nlocus`))+geom_point(aes(col=`Single mutant\nlocus`),shape=21,size=2,stroke=1)+
	scale_colour_manual(values=alpha(cols,c(1,1,1,0.3)))+
	theme(strip.text.x=element_blank(),strip.text.y=element_blank())+
	ylim(limits)+xlim(limits)+
	xlab('fraction ON in glucose')+ylab('fraction ON in galactose')+
	theme(legend.position='bottom')
if(writeplot==T){
	w<-2.8;h<-2.6
	ggsave(paste0(figout,date,'-181119-pOKJG-FraconGluVsFraconGal.pdf'), pOKJG,width=w,height=h)	
}


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

test<-a[,lapply(.SD,function(x)mean(log10(x+0.001),na.rm=T)),by=bycols2]
cnameds<-a[is.na(dens__2.195.gal)]$clone.named
table(DT00[clone.named%in%cnameds&is.na(growth.rate)]$plate.gal)

summ2<-summ1[,lapply(.SD,function(x){
	scale(x)
}),.SDcols=c(phenocols)]
lv<-!is.nan(as.numeric(paste0(summ2[1])))
keep<-colnames(summ2)[lv]
summ<-summ2[,keep,with=F]

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
PI<-phenDT[,mean(phenotypic.index,na.rm=T),by='clone.named']
DTm11<-merge(DTm1,PI,by='clone.named')
hdbtest<-hdbscan((grepcols(DTm1,'dens')),minPts=20)
DTm1[,cluster:=hdbtest$cluster]
DTm1[,cluster_outlier_score:=hdbtest$outlier_score]


############
### classify outlying
############

denscolconv<-function(DT, vnull =T){
	DT[,variable:=gsub('dens__','',variable)]
	DT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
	DT[,sugar:=paste0('g',sugar)]
	if(vnull ==T)DT[,variable:=NULL]
}

densMelt<-function(DTX,meltby=c('clone.named','cluster_test'),summby=c('cluster_test','sugar','AU.FL'),ylimits=c(0,0.35)){
	mDT<-melt(DTX,id=meltby)#;mDT<-mDT[order(sugar,AU.FL,decreasing=T)]
	denscolconv(mDT)
	mDT[,dummyline:=paste0('z',sugar)]
	mDT[,density_mean:=mean(value,na.rm=T),by=summby]
	mDT[,density_sd:=sd(value,na.rm=T),by=summby]
	# mDT[density_mean>=max(ylimits),density_mean:=max(ylimits)]
	# mDT[,bb:=reverse_lev(bb)]
	return(mDT)
}
ExpDensLinePlot<-function(mDT1,ylimits=c(0,0.35)){
	xlimits<-c(-0.3,3)
	xbreaks<-c(0,1,2)
	xlabs<-c('0','1','2')
	ybreaks<-c(0.0,0.1,0.2,0.3)
	ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
	ggplot(mDT1[order(sugar)],aes(AU.FL,density_mean))+#geom_point(data=mDT1,aes(AU.FL,value,col=sugar),size=0.1,alpha=0.3,shape=21)+
		theme_NewPub+
		geom_ribbon(aes(x=AU.FL,ymax=density_mean+density_sd,ymin=density_mean-density_sd,fill=sugar),alpha=0.3)+
		geom_line(aes(col=sugar))+
		scale_colour_manual(values= factorTwoSugar,na.value='transparent')+
		ylab('density') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
#		scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
#		scale_x_continuous(limits=xlimits,breaks=xbreaks,labels=xlabs)+
		theme(legend.position='none')
}

dSub1<-DTm1[cluster=='0'|cluster=='1',c(grepincols(DTm1,'dens__'),'cluster','clone.named'),with=F]
hdbtest0<-hdbscan(dSub1[, grepincols(dSub1,'dens__'),with=F],minPts=8)
dSub1$cluster_test<-hdbtest0$cluster
DTX<-merge(densDT[clone.named%in%dSub1$clone.named,c(grepincols(densDT,'dens__'),'clone.named'),with=F],dSub1[,.(clone.named,cluster_test)],by='clone.named')
mDT1<-densMelt(DTX)
p1<-ExpDensLinePlot(mDT1)+facet_grid(~cluster_test)
dSub1[,cluster_new:=as.character(cluster_test)]
dSub1[cluster==0&cluster_test=='1',cluster_new:='inducible']
dSub1[cluster==0&cluster_test=='2',cluster_new:='constitutive']
dSub1[cluster==0&cluster_test=='3',cluster_new:='const, low expr']
dSub1[cluster==0&cluster_test=='4',cluster_new:='weak expr, induc']
dSub1[cluster==0&cluster_test=='5',cluster_new:='weak expr, const']
DTm2<-merge(DTm1,dSub1[,.(clone.named,cluster_new)],by='clone.named',all=T)
DTm2[is.na(cluster_new),cluster_new:=as.character(cluster)]
DTm2[cluster==0]


dSub1<-DTm2[cluster_new=='0'|cluster=='7',c(grepincols(DTm1,'dens__'),'cluster_new','clone.named'),with=F]
hdbtest0<-hdbscan(dSub1[, grepincols(dSub1,'dens__'),with=F],minPts=8)
dSub1$cluster_test<-hdbtest0$cluster
DTX<-merge(densDT[clone.named%in%dSub1$clone.named,c(grepincols(densDT,'dens__'),'clone.named'),with=F],dSub1[,.(clone.named,cluster_test)],by='clone.named')
mDT1<-densMelt(DTX)
p1<-ExpDensLinePlot(mDT1)+facet_grid(~cluster_test)
dSub1[,cluster_new2:=as.character(cluster_test)]
dSubx<-dSub1[,.(clone.named,cluster_new,cluster_test,cluster_new2)]
dSubx[cluster_new =='0'&cluster_test=='1',cluster_new2:='partially inducible']
dSubx[cluster_new =='0'&cluster_test=='2',cluster_new2:='constitutive']
DTm3<-merge(DTm2,dSubx,by='clone.named',all=T)
DTm3[is.na(cluster_new2),cluster_new2:=as.character(cluster_new)]
DTm3[cluster==0]


dSub1<-DTm3[cluster_new2=='0'|cluster=='6',c(grepincols(DTm1,'dens__'),'cluster_new2','clone.named'),with=F]
hdbtest0<-hdbscan(dSub1[, grepincols(dSub1,'dens__'),with=F],minPts=4)
dSub1$cluster_test<-hdbtest0$cluster
DTX<-merge(densDT[clone.named%in%dSub1$clone.named,c(grepincols(densDT,'dens__'),'clone.named'),with=F],dSub1[,.(clone.named,cluster_test)],by='clone.named')
mDT1<-densMelt(DTX)
p1<-ExpDensLinePlot(mDT1)+facet_grid(~cluster_test)
dSub1[,cluster_new3:=as.character(cluster_test)]
dSub1[cluster_test=='1',cluster_new3:='inducible']
dSub1[cluster_test=='2',cluster_new3:='constitutive']
dSub1[cluster_test=='3',cluster_new3:='weak expr, induc']
dSub1[cluster_test=='4',cluster_new3:='weak expr, const']
dSub1[cluster_test=='5',cluster_new3:='const, low expr']
DTm4<-merge(DTm3,dSub1[,.(clone.named,cluster_new3)],by='clone.named',all=T)
DTm4[is.na(cluster_new3),cluster_new3:=as.character(cluster_new2)]
table(DTm4$cluster_new3)


DTm2[,cluster_A:=cluster_new]
# categorize all clusters
DTm2[cluster_new==c('0'),cluster_A:='outlying']
DTm2[cluster_new==c('1'),cluster_A:='weak expr, induc']
DTm2[cluster_new%in%c('2'),cluster_A:='const, low expr']
DTm2[cluster_new%in%c('3'),cluster_A:='const, low expr']
DTm2[cluster_new%in%c('4'),cluster_A:='const, low expr']
DTm2[cluster_new%in%c('5'),cluster_A:='uninducible']
DTm2[cluster_new%in%c('6'),cluster_A:='inducible']
DTm2[cluster_new%in%c('7'),cluster_A:='constitutive']
table(DTm2$cluster_new)
table(DTm2 $cluster_A)


dSub1<-DTm2[cluster_A=='outlying',c(grepincols(DTm1,'dens__'),'cluster','clone.named'),with=F]
hdbtest0<-hdbscan(dSub1[, grepincols(dSub1,'dens__'),with=F],minPts=8)
dSub1$cluster_test<-hdbtest0$cluster
DTX<-merge(densDT[clone.named%in%dSub1$clone.named,c(grepincols(densDT,'dens__'),'clone.named'),with=F],dSub1[,.(clone.named,cluster_test)],by='clone.named')
mDT1<-densMelt(DTX)
p1<-ExpDensLinePlot(mDT1)+facet_grid(~cluster_test)

mDT2<-densMelt(DTX,summby=c('cluster_test','clone.named','sugar','AU.FL'))
ExpDensLinePlot(mDT2)+facet_wrap(~clone.named+cluster_test,ncol=9)


dSub2<-DTm2[cluster_new=='0'|cluster_new=='6',c(grepincols(DTm2,'dens__'),'cluster_new','clone.named'),with=F]
hdbtest0<-hdbscan(dSub1[, grepincols(dSub1,'dens__'),with=F],minPts=8)
dSub2$cluster_test<-hdbtest0$cluster
DTX<-merge(densDT[clone.named%in%dSub$clone.named,c(grepincols(densDT,'dens__'),'clone.named'),with=F], dSub2[,.(clone.named,cluster_test)],by='clone.named')
mDT1<-densMelt(DTX)
p1<-ExpDensLinePlot(mDT1)+facet_grid(~cluster_test)



dSub2<-DTm1[cluster=='0'|cluster=='7',c(grepincols(DTm1,'dens__'),'cluster','clone.named'),with=F]
hdbtest0<-hdbscan(dSub[, grepincols(dSub,'dens__'),with=F],minPts=8)
dSub$cluster_test<-hdbtest0$cluster
DTX<-merge(densDT[clone.named%in%dSub$clone.named,c(grepincols(densDT,'dens__'),'clone.named'),with=F],dSub[,.(clone.named,cluster_test)],by='clone.named')
mDT2<-densMelt(DTX)
p2<-ExpDensLinePlot(mDT2)+facet_grid(~cluster_test)

do.call(grid.arrange,list(p1,p2))

dSub[,cluster_new:=cluster]
dSub[cluster_test==]

dSub<-DTm1[cluster=='0'|cluster=='7',c(grepincols(DTm1,'dens__'),'cluster'),with=F]
hdbtest0<-hdbscan(dSub[, grepincols(dSub,'dens__'),with=F],minPts=8)

dSub<-DTm1[cluster=='0'|cluster=='7'|cluster=='6',c(grepincols(DTm1,'dens__'),'cluster'),with=F]
hdbtest0<-hdbscan(dSub[, grepincols(dSub,'dens__'),with=F],minPts=8)


############
### classify clusters
############

clusterframe0<-DTm1[,.(clone.named,aa,bb,cc,cluster,cluster_outlier_score)]

# categorize all clusters
clusterframe0[cluster%in%c(0),cluster_A:='outlying']
clusterframe0[cluster%in%c(1),cluster_A:='weakly inducible']
clusterframe0[cluster%in%c(2),cluster_A:='const, low expr 1']
clusterframe0[cluster%in%c(3),cluster_A:='const, low expr 2']
clusterframe0[cluster%in%c(4),cluster_A:='const, low expr 3']
clusterframe0[cluster%in%c(5),cluster_A:='uninducible']
clusterframe0[cluster%in%c(6),cluster_A:='inducible']
clusterframe0[cluster%in%c(7),cluster_A:='constitutive']


# combine const low expression clusters
clusterframe0[cluster%in%c(0),cluster_B:='outlying']
clusterframe0[cluster%in%c(1),cluster_B:='weakly inducible']
clusterframe0[cluster%in%c(2),cluster_B:='const, low expr']
clusterframe0[cluster%in%c(3),cluster_B:='const, low expr']
clusterframe0[cluster%in%c(4),cluster_B:='const, low expr']
clusterframe0[cluster%in%c(5),cluster_B:='uninducible']
clusterframe0[cluster%in%c(6),cluster_B:='inducible']
clusterframe0[cluster%in%c(7),cluster_B:='constitutive']

clusterframe0[,cluster_C:=cluster_B]

genos<-c('clone.named','aa','bb','cc')
clusterframe0[cluster_A=='constitutive']
DT0<-copy(phenDT)
DT0[,foldinduction:=yfp.mean.gal/yfp.mean.glu]
DT1<-merge(DT0,clusterframe0,by=genos)
summ<-DT1[,summary.func(foldinduction),by=c(genos,'cluster_C')]
ggplot(summ,aes(y=log10(mean),x=cluster_C))+geom_violin()+geom_hline(yintercept=0.92)+geom_hline(yintercept=1.4)



############
### sub-classify constitutive
############
summ[cluster_A=='constitutive' &mean<10^0.92,cluster_test:='constitutive']
summ[cluster_A=='constitutive' &mean>=10^0.92&mean<10^1.4,cluster_test:='partially constitutive']
summ[cluster_A=='constitutive' &mean>=10^1.4,cluster_test:='leaky glucose']


ylimits<-c(0,0.35)

denscolconv<-function(DT, vnull =T){
	DT[,variable:=gsub('dens__','',variable)]
	DT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
	DT[,sugar:=paste0('g',sugar)]
	if(vnull ==T)DT[,variable:=NULL]
}
DTA<-copy(densDT[clone.named%in%summ[!is.na(cluster_test)]$clone.named,c(genos,grepincols(densDT,'dens__')),with=F])
DTX<-merge(DTA,summ[,.(clone.named,cluster_test)],by='clone.named')
meltby<-c(genos,'cluster_test')

mDT<-melt(DTX,id=meltby)#;mDT<-mDT[order(sugar,AU.FL,decreasing=T)]
denscolconv(mDT)
mDT[,dummyline:=paste0('z',sugar)]
mDT[,density_mean:=mean(value,na.rm=T),by=c('cluster_test','sugar','AU.FL')]
mDT[,density_sd:=sd(value,na.rm=T),by=c('cluster_test','sugar','AU.FL')]
mDT[density_mean>=max(ylimits),density_mean:=max(ylimits)]
mDT[,bb:=reverse_lev(bb)]

mDT1<-mDT

ExpDensLinePlot<-function(mDT1,ylimits=c(0,0.35)){
	xlimits<-c(-0.3,3)
	xbreaks<-c(0,1,2)
	xlabs<-c('0','1','2')
	ybreaks<-c(0.0,0.1,0.2,0.3)
	ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
	ggplot(mDT1[order(sugar)],aes(AU.FL,density_mean))+#geom_point(data=mDT1,aes(AU.FL,value,col=sugar),size=0.1,alpha=0.3,shape=21)+
		theme_NewPub+
		geom_ribbon(aes(x=AU.FL,ymax=density_mean+density_sd,ymin=density_mean-density_sd,fill=sugar),alpha=0.3)+
		geom_line(aes(col=sugar))+
		scale_colour_manual(values= factorTwoSugar,na.value='transparent')+
		ylab('density') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
		scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
		scale_x_continuous(limits=xlimits,breaks=xbreaks,labels=xlabs)+
		theme(legend.position='none')
}

ExpDensLinePlot(mDT1)+facet_grid(~cluster_test)








############
### plot out clusters
############

PI<-phenDT[,mean(phenotypic.index,na.rm=T),by='clone.named']
DTm <-merge(DTm1,PI,by='clone.named')

rankClust<-function(DT){
	DT[,cluster.orig:=cluster]
	DT[,V2:=median(V1),by='cluster']
	ttt<-DT[,median(V1),by='cluster'][order(V1)]
	ttt[,cluster:=1:nrow(ttt)]
	setnames(ttt,"V1",'V2')
	out<-merge(ttt,DT[,!'cluster'],by='V2')
	out[,!c('V2','cluster.orig'),with=F][order(cluster,V1)]
}

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
pLBXR<-ggplot(mDT,aes(X,Y))+geom_tile(aes(fill=`density, log10`))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),na.value='white',limits=collim)+
	theme(axis.text.x=element_text(angle=90,vjust=0.4,hjust=1,size=8))+
	theme(axis.text.y=element_text(size=6))+
	geom_vline(xintercept=61,col='black',size=2)+
	geom_vline(xintercept=61,col='white',size=1.5)+
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
pYIAD<-ggplot(mDTm,aes(AU.FL,V2))+#geom_point(data=mDTm,aes(AU.FL,density,col=sugar),size=0.1,alpha=0.3,shape=21)+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(cluster_new+facet~ mutant_facet)+scale_colour_manual(values=c(factorTwoSugar,'grey90','grey90'),na.value='transparent')+
	ylab('density of observations') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	theme_NewPub+
	geom_text(data=labelframe,aes(x=x,y=y,label=V1),size=2,col='grey50')
if(writeplot==T){
	w<-5;h<-4.2
	ggsave(paste0(figout,'1707xx-181119-pYIAD-ExpressionDensityClusters-LinePlot-FacetByMutatedGene.pdf'), pYIAD,width=w,height=h)
}

############
### plot archetypal expression densities: combine all const low expr 
############

ylimits<-c(0,0.35)

denscolconv<-function(DT, vnull =T){
	DT[,variable:=gsub('dens.','',variable)]
	DT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
	DT[,sugar:=paste0('g',sugar)]
	if(vnull ==T)DT[,variable:=NULL]
}
DTX<-copy(densDT)
mDT<-melt(DTX,id=genos)#;mDT<-mDT[order(sugar,AU.FL,decreasing=T)]
denscolconv(mDT)
mDT[,dummyline:=paste0('z',sugar)]
mDT[,density_mean:=mean(value,na.rm=T),by=c(genos,'sugar','AU.FL')]
mDT[,density_sd:=sd(value,na.rm=T),by=c(genos,'sugar','AU.FL')]
mDT[density_mean>=max(ylimits),density_mean:=max(ylimits)]
mDT[,bb:=reverse_lev(bb)]

mDT1<-mDT

ExpDensLinePlot<-function(mDT1,ylimits=c(0,0.35)){
	xlimits<-c(-0.3,3)
	xbreaks<-c(0,1,2)
	xlabs<-c('0','1','2')
	ybreaks<-c(0.0,0.1,0.2,0.3)
	ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
	ggplot(mDT1[order(sugar)],aes(AU.FL,density_mean))+#geom_point(data=mDT1,aes(AU.FL,value,col=sugar),size=0.1,alpha=0.3,shape=21)+
		theme_NewPub+
		geom_ribbon(aes(x=AU.FL,ymax=density_mean+density_sd,ymin=density_mean-density_sd,fill=sugar),alpha=0.3)+
		geom_line(aes(col=sugar))+
		scale_colour_manual(values= factorTwoSugar,na.value='transparent')+
		ylab('density') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
		scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
		scale_x_continuous(limits=xlimits,breaks=xbreaks,labels=xlabs)+
		theme(legend.position='none')
}

G4<-c('GAL4.WT','GAL4-L868P','GAL4-L868K','GAL4-L868G') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80S-1','GAL80.07')
GK<-c('GALK.Sac_cer','GALK.Esc_col')
G3<-c('GAL3.WT','GAL3.delta')
mDT1<-mDT[aa%in%G4&cc%in%G80&dd%in%GK&bb%in%G3]
sample(LETTERS,4)




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
mDTm[,V2:=mean(density,na.rm=T),by=c(bycols,'mutant_gene')]
mDTm[,V3:=sd(density,na.rm=T),by=c(bycols,'mutant_gene')]
mDTm[,upper:=V2+V3]
mDTm[,lower:=V2-V3]
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
mDTm[,alphas:=0.1/.N,by=c('cluster_new','facet','mutant_facet')]
pUNWE<-ggplot(mDTm,aes(AU.FL,V2))+geom_point(data=mDTm,aes(AU.FL,density,col=sugar,alpha=alphas),size=0.1,shape=21)+theme_minimal()+
	geom_ribbon(aes(x=AU.FL,ymax=V2+0.0025,ymin=V2-0.0025,col=dummyline),size=2,alpha=0.3,inherit.aes=F)+geom_line(aes(col=sugar))+
	facet_grid(cluster_new+facet~ mutant_facet)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density of observations') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))+
	geom_text(data=labelframe,aes(x=x,y=y,label=V1),size=3.5,col='grey50')

if(writeplot==T){
	w<-8;h<-6.5
	ggsave(paste0(figout,'1707xx-180720-pUNWE-ExpressionDensityClusters-AggregatedConstitutiveOutliers-LinePlot-FacetByMutatedGene.png'), pUNWE,width=w,height=h)
}

