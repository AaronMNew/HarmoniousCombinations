# code for publication

source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'
source(source.code.path.file)
# load.all() # this function loads all the standard packages you use in this script

head.dir<-"/Users/anew/Google Drive/002_lehnerlabData/180404-VarGALK.With.G80xG4_2"
date<-'180404' # beginning date clones were measured
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/publication/'
outputdata<-'./output.data/'
layoutdir<-'./layout.sample.data/'
setwd(head.dir)
source.code.path.file2<-'/Users/anew/Google Drive/002_lehnerlabData/180404-VarGALK.With.G80xG4_2/R.functions/180404-CustomScripts.R'
source(source.code.path.file2)
aesthetics.dir<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/181119-Aesthetics_Harmonious_Combinations.R'
source(aesthetics.dir)

writeplot<-F
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
### supplement: Comparison of GALKs
###################################################
DTX<-data.table(DTX[,lapply(.SD,scale),.SDcols= notin(phenocols,'phenotypic.index')],phenDT[,genos,with=F])#copy(phenDT)
summ<-summary.func.all(DTX,notin(phenocols,'phenotypic.index'),genos)
ggplot(phenDT[dd=='HIS5.Sch_pom'],aes(yfp.mean.gal, growth.rate))+geom_point()+geom_smooth()
summ[,clone.named:=applyPaste(data.frame(aa,bb,cc,dd),' ')]
bylist<-c(genos,'clone.named')

meltby<-c('aa','bb','cc','dd','clone.named')
leftside<-meltby
meltcast.summary<-function(summ,meltby,leftside,rightside='stat'){
	# summ is a summary data.table made by summary.func.all() custom function
	mDT<-melt(summ,id=meltby)
	mDT[,c('pheno','stat'):=colsplit(variable,'_',c('pheno','stat'))]
	form<-castform(c(leftside,'pheno'),rightside)
	return(dcast(mDT,form,value=stat))
}


cDT<-meltcast.summary(summ,meltby=meltby,leftside=meltby)
cDT[,NotGK:=applyPaste(data.frame(aa,bb,cc))]
cDT[,GK:=dd]

DTl<-split(cDT,by='GK')
DTl2<-lapply(DTl,function(DTx){
	xlabel<-unique(DTx$dd)
	print(xlabel)
	DTm<-merge(DTx,cDT[dd!=xlabel],by=c('NotGK','pheno'),allow.cartesian=T)
})

DT1<-data.table(ldply(DTl2))

ggplot(DT1[pheno=='growth.rate'],aes(mean.x,mean.y,xmax=upper.x,xmin=lower.x,ymin=lower.y,ymax=upper.y))+
	geom_errorbarh(col='grey70')+geom_errorbar(col='grey70')+
	geom_point(shape=21)+facet_grid(dd.y~dd.x)+geom_abline(col='grey70')

sample(LETTERS,4)

# plot all the phenotypes paired between GALK orthologues or HIS5
pBAGQ<-ggplot(DT1,aes(mean.x,mean.y,xmax=upper.x,xmin=lower.x,ymin=lower.y,ymax=upper.y))+
	theme_NewPub+
	geom_errorbarh(aes(col=pheno))+geom_errorbar(aes(col=pheno))+
	geom_point(shape=21,size=1,stroke=0.1,aes(col=pheno))+facet_grid(dd.y~dd.x)+geom_abline(col='grey70')+
	xlab('scaled within-genotype phenotypic mean x')+ylab('scaled within-genotype phenotyipic mean y')

# C. albicans grows more slowly at higher growth rates. do a kolmogorov-smirnov test and anderson-darling test to put some "significance" number on this and to show with a plot of ECDF
library(kSamples)
statset<-sapply(phenocols,function(x){
	xx<-unlist(phenDT[dd=='GALK.Can_alb',eval(x),with=F])
	xy<-unlist(phenDT[dd=='GALK.Esc_col',x,with=F])
	ao<-ks.test(xx,xy)
	bo<-ad.test(xx,xy)
	round(c(ks.stat=ao$statistic,ks.pval=ao$p.value,ad.stat=data.frame(bo$ad)[,1][1],ad.pval=data.frame(bo$ad)[,3][1]),3)
})
print('table of non-parametric stat comparisons for comparison of phenotypic distributions between coli and albicans GALK values');print(statset)


ColVsAlb<-phenDT[dd%in%c('GALK.Esc_col','GALK.Can_alb'),mean(growth.rate,na.rm=T),by=genos]
pUZVC<-ggplot(ColVsAlb,aes(col=dd,x=V1))+stat_ecdf()+
	theme_NewPub+
	xlab(expression(paste('growth rate ',mu,' ',hr^-1)))+
	ylab('cumulative disribution')+
	theme(legend.position=c(0.7,0.25))+
	guides(colour=guide_legend(title='GALK ortholog'))+
	scale_colour_manual(values=factorDefault)


if(writeplot==T){
	w=5;h=4
	ggsave(paste0(figout,'180404-181214-pBAGQ-GALK_ortholog_comparisons.pdf'), pBAGQ,width=w,height=h)
	w=2.15;h=1.9
	ggsave(paste0(figout,'180404-181214-pUZVC-ECDF_growth.rate_ColiVSAlb.pdf'), pUZVC,width=w,height=h)

}


###################################################
### Figure 3: Scatterplot illustrating clones that can grow without GAL3
###################################################

DTx<-copy(phenDT)
DTx[dd%in%c('GALK.Sac_cer'),GK:='GAL1']
DTx[dd%in%c('GALK.Esc_col','GALK.Can_alb'),GK:='GALK']
DTx1<-DTx[dd!= 'HIS5.Sch_pom']
genosx<-c('aa','bb','cc','GK')

summ<-DTx1[,summary.func(growth.rate),by=genosx]
form<-castform(c('aa','cc','bb'),c('GK'))
summnames<-names(summary.func(1:4))
cDT1<-dcast(summ,form,value.var=summnames)
ttests<-data.table(t(apply(cDT1,1,function(xa){
	xb<-as.list(xa)
	x<-lapply(xb,as.numeric)
	names(x)<-names(xa)
	t.test2(x$mean_GAL1,x$mean_GALK,x$sd_GAL1,x$sd_GALK,x$N_GAL1,x$N_GALK)
	})))
ttests[,padj:=p.adjust(p.value,'fdr')]
ttests[,`FDR < 0.05`:=as.character(padj<=0.05& abs(diff.of.means)>0.03)]
cDT<-data.table(cDT1,ttests)
# cDT[aa=='GAL4.WT'&cc=='GAL80.WT',`FDR < 0.05`:='WT GAL80 and GAL4']
sample(letters,4)
cDT[,bb:=reverse_lev(bb)]
pJTVD<-ggplot(cDT[order(`FDR < 0.05`,decreasing=T)],aes(mean_GAL1,mean_GALK,xmin=lower_GAL1,xmax=upper_GAL1,ymin=lower_GALK,ymax=upper_GALK))+
	theme_NewPub+
	geom_errorbar(col='grey50')+geom_errorbarh(col='grey50')+geom_point(aes(col=`FDR < 0.05`),size=2,shape=21,stroke=0.2)+
	scale_colour_manual(values=c('black','red','cornflowerblue'))+
	geom_abline(col='grey70')+	
	geom_point(data=cDT[aa=='GAL4.WT'&cc=='GAL80.WT'], aes(mean_GAL1,mean_GALK),col='cornflowerblue',size=3,shape=21,stroke=1.5)+
	theme(legend.position='bottom')+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+
	ylab(expression(paste('GALK - sensor (-) ',mu,' ',hr^-1)))+
	xlab(expression(paste('GAL1 - sensor (+) ',mu,' ',hr^-1)))+
	facet_grid(~bb)

if(writeplot==T){
	w=5;h=3
	ggsave(paste0(figout,'180404-181113-pJTVD-GAL1VsGALK-FacetByGAL3.pdf'), pJTVD,width=w,height=h)

}

# GAL4 alleles involved with inducible phenotypes: 
aas<-unique(cDT[`FDR < 0.05`==T]$aa)
# GAL80 alleles invovled with inducible phenotypes:
ccs<-unique(cDT[`FDR < 0.05`==T]$cc)


###################################################
### Figure 3: quantify effect of leaky expression on GAL pathway
###################################################

# summarize data
summ<-summary.func.all(DT0,c('growth.rate','fracon.glu','yfp.mean.glu','yfp.mean.gal','mean.on.fraction.glu'),bylist=genos)
summ[aa=='GAL4.WT'&bb=='GAL3.WT']
summ[,clone.named:=applyPaste(data.frame(aa,cc),' ')]
summ[,WT:=grepl(c('GAL4.WT GAL80.WT'),clone.named)]


############
### Figure 3: fit growth rate in galactose as a logistic function of GAL pathway expression in glucose
############

# annotate a smaller version of DT0
dSub1<-DT0[,.(aa,cc,bb,dd,yfp.mean.glu,growth.rate)]
dum<-data.table(dd=c('GALK.Esc_col','GALK.Can_alb','GALK.Sac_cer','HIS5.Sch_pom'),GALK.cat=c('GALK\nsensor(-)','GALK\nsensor(-)','GAL1\nsensor(+)','HIS5\nsensor(-)'))
dum[,GALK.cat:=factor(GALK.cat,levels=unique(dum$GALK.cat)[c(2,1,3)])]
dSub<-merge(dSub1,dum,by='dd')
dSub[,x:= yfp.mean.glu]
dSub[,y:=growth.rate]
dSub[GALK.cat=='GALK\nsensor(-)'&bb=='GAL3.delta',fac:='no sensors']

# fit logistic model and calculate deviation from expectation
dat<-dSub[fac=='no sensors',.(x,y)]
mod<-nls(y~SSlogis(x,Asym,xmid,scal),data=dat)
xes<-data.table(x=seq(1,50000,length=100))
lineframe<-data.table(xes,y=predict(mod,xes))
datcoli<-dSub[dd=='GALK.Esc_col'&bb=='GAL3.delta',.(x,y)]
modcoli<-nls(y~SSlogis(x,Asym,xmid,scal),data=datcoli)
datalb<-dSub[dd=='GALK.Can_alb'&bb=='GAL3.delta',.(x,y)]
modalb<-nls(y~SSlogis(x,Asym,xmid,scal), data =datalb)
dSub[dd=='GALK.Sac_cer'|dd=='HIS5.Sch_pom',pred:=predict(mod,newdata=data.frame(x,y))]
dSub[dd=='GALK.Esc_col',pred:=predict(modcoli, newdata=data.frame(x,y))]
dSub[dd=='GALK.Can_alb',pred:=predict(modalb, newdata=data.frame(x,y))]

# merge summary and raw data
DTm<-merge(summ,dSub,by=genos)
DTm[,clone.named:=applyPaste(data.frame(aa,bb,cc,dd),' ')]
length(unique(DTm$clone.named))
DTm[,pred:=predict(mod,DTm[,.(x)])]
DTm[,estimate:=y-pred]
# ggplot(DTm,aes(pred,y,col=GALK.cat))+geom_point()+facet_wrap(~bb)+geom_abline()
tframe<-DTm[,t.test(estimate,alternative='greater'),by=genos]
tframe[,p.adj:=p.adjust(p.value)]
tframe[,sig:=p.adj<0.05]
tframe[,`FDR < 0.05`:=sig=='TRUE'&abs(estimate)>0.03]
tframe[,dup:=applyPaste(tframe[,.(aa,bb,cc,dd)],'')]
DTm2<-merge(DTm,tframe[!duplicated(dup),c(genos,'sig','FDR < 0.05'),with=F],by=genos)
DTm2[,bb:=reverse_lev(bb)]
sample(LETTERS,4)
# smaller plot of expression vs growth rate without HIS5 for main fig - strip text normal
pALUN<-ggplot(DTm2[order(`FDR < 0.05`)][dd!='HIS5.Sch_pom'],aes(yfp.mean.glu_mean,growth.rate_mean,xmax= yfp.mean.glu_upper,xmin=yfp.mean.glu_lower,ymax= growth.rate_upper,ymin=growth.rate_lower))+
	geom_errorbar(col='grey50')+geom_errorbarh(col='grey50')+geom_point(aes(col=`FDR < 0.05`),size=2,shape=21,stroke=0.2)+
	geom_line(data=lineframe,aes(x,y),inherit.aes=F)+
	ylab(expression(paste('mean ',mu,' ',hr^-1,' in galactose')))+
	xlab('mean GAL pathway expression in glucose, A.U.')+
	facet_grid(GALK.cat~bb)+
	theme_NewPub+ 
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+
	scale_colour_manual(values=c('black','red'))+
	theme(legend.position='bottom')

# GAL4 alleles involved with inducible phenotypes: 
aas<-unique(DTm2[`FDR < 0.05`=='TRUE']$aa)
# GAL80 alleles involved with inducible phenotypes:
ccs<-unique(DTm2[`FDR < 0.05`=='TRUE']$cc)


############
### Figure 3: combine two plots into one
############
library(ggplot2)   
library(gtable)    
library(grid)
library(gridExtra) 
gA<-ggplotGrob(pJTVD)
gB<-ggplotGrob(pALUN)
# Set the widths
gA$widths <- gB$widths

# Arrange the charts.
grid.newpage()
grid.arrange(gA, gB,ncol=2)
if(writeplot==T){
	w<-5.25;h<-2.5
	ggsave(paste0(figout,date,'-180828-pJTVD&pALUN-SensorEffects-Figure3.pdf'), grid.arrange(gA, gB,ncol=2),width=w,height=h)
	}

# plot of background GAL pathway expression in glucose vs growth rate for different sensor / GALK backgrounds 
p1<-ggplot(DTm2[order(`FDR < 0.05`)],aes(yfp.mean.glu_mean,growth.rate_mean,xmax= yfp.mean.glu_upper,xmin=yfp.mean.glu_lower,ymax= growth.rate_upper,ymin=growth.rate_lower))+
#	geom_point(data=DTm,aes(x,y),shape=21,size=1,col='grey50',inherit.aes=F)+
	geom_errorbar(col='grey50')+geom_errorbarh(col='grey50')+geom_point(aes(col=`FDR < 0.05`),size=2,shape=21,stroke=1)+
	geom_line(data=lineframe,aes(x,y),inherit.aes=F)+ylab(expression(paste('mean ',mu,' ',hr^-1,' in galactose')))+xlab('mean expression of GAL pathway in glucose, A.U.')+
	facet_grid(GALK.cat~bb)+theme_minimal()+ theme(axis.text.x = element_text(angle = 30, hjust = 1))+scale_colour_manual(values=c('black','red'))+
	theme(strip.text.y = element_text(angle =0),legend.position='bottom')

# tile plot showing differences and significance between deviation of observed growth rate from expectation based on no-sensor background growth rate vs YFP in glucose
p2<-ggplot(DTm2,aes(aa,cc))+geom_tile(aes(fill=`observed -\nexpected`))+scale_fill_gradientn(colours=c('cornflowerblue','white','indianred'),limits=round(collim,2))+
	facet_grid(GALK.cat~bb)+geom_point(aes(aa,cc,col=`FDR < 0.05`))+scale_colour_manual(values=c('transparent','red'))+theme_minimal()+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+xlab('GAL4 alleles')+ylab('GAL80 alleles')+
	theme(strip.text.y = element_text(angle =0),legend.position='bottom')


###################################################
### Figure 3: expression density line plots of genotypes
###################################################

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
pZVOQ<-ExpDensLinePlot(mDT1)+facet_grid(bb+dd~cc+aa)
pZVOQLegend<-ExpDensLinePlot(mDT1)+facet_grid(bb+dd~cc+aa)+theme(legend.position='bottom')
if(writeplot==T){
	w<-8;h<-4
	ggsave(paste0(figout,date,'-181113-pZVOQ-ExpressionDistributions-WithLegend-SubsetForPub.pdf'), pZVOQLegend,width=w,height=h)	
}

pLETD<-ExpDensLinePlot(mDT)+facet_grid(cc+aa~bb+dd)
if(writeplot==T){
	w<-8;h<-35
	ggsave(paste0(figout,date,'-181113-pLETD-ExpressionDistributions-AllGenotypes.pdf'), pLETD,width=w,height=h, limitsize = FALSE)
}


###################################################
### Figure 3: pairwise fitness landscapes for GAL3, GAL1 across different GAL4, GAL80 backgrounds
###################################################

############
### pairwise fitness landscapes: growth rate
############

# summarize data and get together
dSub<-DT0[!dd%in%c('HIS5.Sch_pom','GALK.Can_alb')]
pheno<-'growth.rate'
backg<-summary.func(unlist(DT0[aa=='GAL4.delta'|dd=='HIS5.Sch_pom',pheno,with=F]))
wt<-summary.func(unlist(DT0[aa=='GAL4.WT'&bb=='GAL3.WT'&cc=='GAL80.WT'&dd=='GALK.Sac_cer',pheno,with=F]))
nogrowth<-data.table(x=seq(0,4,length=100),upper=backg$upper,lower=backg$lower)
wtgrowth<-data.table(x=seq(0,4,length=100),upper=wt$upper,lower=wt$lower)
midval1<-2.05
midval2<-1.95 
summ<-summary.func.all(dSub,c('growth.rate','phenotypic.index','yfp.mean.gal','yfp.mean.glu'),c('aa','bb','cc','dd'))
summ[dd=='GALK.Esc_col'&bb=='GAL3.delta',c('x','shape','shape2'):=list(3,factor(11),'d.GALS')]
summ[dd=='GALK.Sac_cer'&bb=='GAL3.delta',c('x','shape','shape2'):=list(midval2,factor(6),'d.GAL3')]
summ[dd=='GALK.Esc_col'&bb=='GAL3.WT',c('x','shape','shape2'):=list(midval1,factor(2),'d.GAL1')]
summ[dd=='GALK.Sac_cer'&bb=='GAL3.WT',c('x','shape','shape2'):=list(1,factor(1),'WT GALS')]
senslev<-c("WT GALS","d.GAL3","d.GAL1","d.GALS")
summ[,shape2:=factor(shape2,levels=senslev)]

arrowframe<-data.table(x0=c(1,1,midval1,midval2),xend=c(midval1,midval2,3,3),shape2=senslev,lev=c(0,1,1,2))
summ2<-merge(summ,arrowframe,by='shape2')
setnames(summ2,paste0(pheno,'_mean'),'phenoi')
sub<-summ2[,.(aa,cc,bb,dd,phenoi,shape2,x0,xend,lev)]
sub0<-sub[lev==0&dd=='GALK.Sac_cer',.(aa,cc,x0=1,y0= phenoi)]
sub01a<-sub[lev==1&dd=='GALK.Esc_col',.(aa,cc,xend=midval1,shape2,y0= phenoi)]
sub01b<-sub[lev==1&dd=='GALK.Sac_cer',.(aa,cc,xend=midval2,shape2,y0= phenoi)]
sub01<-rbind(sub01a,sub01b)
setnames(sub01,'y0','yend')
sub01c<-merge(sub0,sub01,by=c('aa','cc'))
sub02<-sub[lev==2,.(aa,cc,xend=3,shape2,yend= phenoi)]
sub02a<-copy(sub01)
setnames(sub02a,c('xend','yend'),c('x0','y0'))
sub02b<-merge(sub02,sub02a[,!'shape2'],by=c('aa','cc'))
dum<-sub02b[1]
dum[,notin(colnames(dum),c('aa','cc'))]#<-0
dum[,shape2:='d.GALS']
arrowframe3<-rbind(sub01c,sub02b,dum)
arrowframe3[,shape3:=factor(shape2,levels=senslev)]


facy<-1#0.96
facx<-0
xlabels<-c('background','single mutants','double mutants')
ynudge<-0.02
DT<-copy(summ)
arrowframe2<-copy(arrowframe3)
arrowframe2[,shape3:=shape2]
df<-merge(InducibleFrame,summ,by=genos)[,.(aa,bb,cc,dd,inducible,shape2)][!is.na(shape2)]
pointframe<-merge(df,arrowframe2,by=c('aa','cc','shape2'))
arrowframe2[,grp:=shape3]
DT[,grp:=shape2]

myColors <- c('cornflowerblue','orange1','grey50','darkgreen','red','black')
names(myColors) <- c('d.GAL3','d.GAL1','WT GALS','d.GALS','inducible','uninducible')
legname<-'GAL sensor genotypes'
colScale <- scale_colour_manual(name = "GAL sensor genotypes",values = myColors)
sample(LETTERS,4)
DTx<-DT
textopt<-geom_text(data=DTx,nudge_y= ynudge,aes(x=x,y=growth.rate_mean,label=shape2,col=(shape2)),size=2)
FitnessLandscapePlot<-function(DTx,arrowframe2,nogrowth,wtgrowth){

	ggplot(DTx,aes(x,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+
		theme_NewPub+
		geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
		geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
		geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='open',length=unit(0.3,'cm'),angle=20),size=0.5,alpha=0.7,inherit.aes=F)+
		geom_errorbar(col='grey50',width=0.4)+#geom_point(col='grey50')+
		colScale+
		xlab('')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
		scale_y_continuous(breaks=c(round(backg$mean,2),round((backg$mean+wt$mean)/2,2),round(wt$mean,2)),limits=c(0.05,0.27))+
		scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+
		facet_grid(cc~aa)+
		theme(legend.position='bottom',
			axis.line.x=element_line(size=0.3,colour='grey50'),
			axis.text.x = element_text(angle = 30, hjust = 1),
			panel.spacing = unit(0.7, "lines"))

}
G4<-c('GAL4.WT','GAL4-L868P','GAL4-L868K','GAL4-L868G') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80S-1','GAL80.07')
GK<-c('GALK.Sac_cer','GALK.Esc_col')
G3<-c('GAL3.WT','GAL3.delta')
ynudge<-0.02
DTx<-DT[aa%in%G4&cc%in%G80]
arrowframe2<-arrowframe3[aa%in%G4&cc%in%G80]
sample(letters,4)
pIKLE<-FitnessLandscapePlot(DTx,arrowframe2,nogrowth,wtgrowth)

DTx<-copy(DT)
arrowframe2<-copy(arrowframe3)
sample(letters,4)
pBNAC<-FitnessLandscapePlot(DTx,arrowframe2,nogrowth,wtgrowth)


if(writeplot==T){
	w<-3.7;h<-4.2
#	ggsave(paste0(figout,date,'-181113-pIKLE-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-ForPubFinalFigure-growth.rate.png'), pIKLE,width=w,height=h)
	ggsave(paste0(figout,date,'-181113-pIKLE-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-ForPubFinalFigure-growth.rate.pdf'), pIKLE,width=w,height=h)
	w<-5;h<-4.5
	ggsave(paste0(figout,date,'-181113-pBNAC-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-Supplemenet-All-growth.rate.pdf'), pBNAC,width=w,height=h)
}




###################################################
### statistical analyses for publication
###################################################

phenDT[aa=='GAL4.WT'&cc=='GAL80.07',mean(yfp.mean.gal/yfp.mean.glu),by=c('bb','dd')]
phenDT[aa=='GAL4.WT'&cc=='GAL80.07',mean(yfp.mean.gal),by=c('bb','dd')]
phenDT[aa=='GAL4.WT'&cc=='GAL80.07',mean(yfp.mean.glu),by=c('bb','dd')]

dSub<-copy(phenDT)
dSub[,cor:=fracon.glu]
dSub[fracon.glu==0,cor:=rnorm(1,mean=0.0001,sd=0.0001)]
dSub[fracon.glu==1,cor:= rnorm(1,mean=0.9999,sd=0.0001)]
dSub[,logfon:=logit(cor)]

padj='fdr'
bylist=c('aa','bb','cc','dd')
pheno=c('logfon')
sig.cutoff = 0.05
ttestsDT<-all.ttest(dSub,pheno,c('aa','bb','cc','dd'),padj=padj,sig.cutoff=sig.cutoff)


# analysis @statHSEI test of effect of GAL3 and GAL1 in GAL80.07 background

aas<-unique(phenDT$aa)#c('GAL4.WT')
bbs<-c('GAL3.WT','GAL3.delta')
ccs<-c('GAL80.07')
dds<-c('GALK.Sac_cer','GALK.Esc_col','GALK.Can_alb')

dSub1<-ttestsDT[aa.x%in%aas&bb.x%in%bbs&cc.x%in%ccs&dd.x%in%dds&aa.y%in%aas&bb.y%in%bbs&cc.y%in%ccs&dd.y%in%dds]
ttests.plot(DT=dSub1,padj=padj,sig.cutoff=sig.cutoff)



dt<-dSub[,.(fracon.glu,aa,bb,cc,dd)]
dt[,cor:=fracon.glu]
dt[fracon.glu<0.001,cor:=0.001]
dt[fracon.glu>0.999,cor:=0.999]
dt[,val:=logit(cor)]
dt1<-dt
tt1<-TukeyHSD(aov(lm(val~(bb+dd)^2, dt1)))
dt2<-(tt1[[length(tt1)]])
rownames(dt2)<-gsub('4-','4.',gsub('S-','S.',gsub('\\:','\\+',rownames(dt2))))
TukeyInteractionFlip2(dt2)
TukeyInteractionFlip2(dt2,plot=F)
summ1<-dt[,summary.func(cor),by=c('bb','dd')]
summ1[,diff:=(summ1[bb=='GAL3.WT'&dd=='GALK.Sac_cer']$mean-mean)][order(bb,dd)]
p1<-ggplot(summ1,aes(bb,dd))+geom_tile(aes(fill=mean))+scale_fill_gradientn(colours=c('black','indianred'))+ggtitle('fraction ON cells')
statHSEI <-glm(cor~bb+dd,dt,family=quasibinomial(link='logit'))
summary(statHSEI)
summary(aov(statHSEI))
data.frame(exp(statHSEI $coefficients)/(1+exp(statHSEI $coefficients)))

dt<-dSub[,.(yfp.mean.glu,aa,bb,cc,dd)]
dt[,val:=(yfp.mean.glu)]
dt1<-dt
tt1<-TukeyHSD(aov(lm(val~(bb+dd)^2, dt1)))
dt2<-(tt1[[length(tt1)]])
rownames(dt2)<-gsub('4-','4.',gsub('S-','S.',gsub('\\:','\\+',rownames(dt2))))
TukeyInteractionFlip2(dt2)
TukeyInteractionFlip2(dt2,plot=F)
summ2<-dt[,summary.func(val),by=c('bb','dd')]
p2<-ggplot(summ2,aes(bb,dd))+geom_tile(aes(fill=mean))+scale_fill_gradientn(colours=c('black','indianred'))+ggtitle('mean YFP expression, A.U.')
statHSEI <-lm(val~bb*dd,dt)
summary(statHSEI)

pHSEI<-list(p1+theme(axis.title.x=element_blank(),axis.title.y=element_blank()),p2+theme(axis.title.x=element_blank(),axis.title.y=element_blank()),ncol=2)
# do.call(grid.arrange,pHSEI)


# @statGKMW test of effect of GAL3 and GAL1 in GAL4C backgrounds on leaky glucose expression
library(betareg)
aas<-c('GAL4-L868G','GAL4-L868K','GAL4-L868P')
bbs<-c('GAL3.WT','GAL3.delta')
ccs<-c('GAL80.WT')
dds<-c('GALK.Sac_cer','GALK.Esc_col','GALK.Can_alb')
dSub<-phenDT[aa%in%aas&bb%in%bbs&cc%in%ccs&dd%in%dds]


dt<-dSub[,.(fracon.glu,aa,bb,cc,dd)]
dt[,cor:=fracon.glu]
dt[fracon.glu<0.001,cor:=0.001]
dt[fracon.glu>0.999,cor:=0.999]
dt[,val:=logit(cor)]
statGKMW <-glm(cor~aa*bb*dd,dt,family=quasibinomial(link='logit')) # betareg(cor~aa*bb*dd,dt)
data.frame(exp(statGKMW$coefficients$mean)/(1+exp(statGKMW$coefficients$mean)))
summary(aov(statGKMW))
tt1<-TukeyHSD(aov(statGKMW))
dt2<-(tt1[[length(tt1)]])
rownames(dt2)<-gsub('4-','4.',gsub('S-','S.',gsub('\\:','\\+',rownames(dt2))))
TukeyInteractionFlip2(dt2)
TukeyInteractionFlip2(dt2,plot=F)
summ1<-dt[,summary.func(cor),by=c('bb','dd','aa')]
p1<-ggplot(summ1,aes(bb,dd))+geom_tile(aes(fill=mean))+scale_fill_gradientn(colours=c('black','indianred'))+ggtitle('fraction ON cells in glucose')+
	facet_grid(~aa)

dt<-dSub[,.(yfp.mean.glu,aa,bb,cc,dd)]
dt[,cor:= yfp.mean.glu]
dt[,val:=(cor)]
statGKMW <-lm(cor~aa*bb*dd,dt) 
summary(aov(statGKMW))
tt1<-TukeyHSD(aov(statGKMW))
dt2<-(tt1[[length(tt1)]])
rownames(dt2)<-gsub('4-','4.',gsub('S-','S.',gsub('\\:','\\+',rownames(dt2))))
TukeyInteractionFlip2(dt2)
TukeyInteractionFlip2(dt2,plot=F)
summ2<-dt[,summary.func(cor),by=c('bb','dd','aa')]
p2<-ggplot(summ2,aes(bb,dd))+geom_tile(aes(fill=mean))+scale_fill_gradientn(colours=c('black','indianred'))+ggtitle('mean YFP expression, A.U.')+
	facet_grid(~aa)

pGKMW<-list(p1+theme(axis.title.x=element_blank(),axis.title.y=element_blank()),p2+theme(axis.title.x=element_blank(),axis.title.y=element_blank()),ncol=2)

# do.call(grid.arrange,pGKMW)

# 
aas<-c('GAL4.WT','GAL4-L868G','GAL4-L868K')
bbs<-c('GAL3.WT','GAL3.delta')
ccs<-c('GAL80.WT','GAL80.07','GAL80S-1')
dds<-c('GALK.Sac_cer','GALK.Esc_col')
dSub<-phenDT[aa%in%aas&bb%in%bbs&cc%in%ccs&dd%in%dds] # phenDT[!aa=='GAL4.WT'&!dd=='HIS5.Sch_pom']
#dSub<-phenDT[!aa=='GAL4.WT'&!dd=='HIS5.Sch_pom']
# fracon.glu
dt<-dSub[,.(fracon.glu,aa,bb,cc,dd)]
dt[,cor:=fracon.glu]
dt[fracon.glu<0.001,cor:=0.001]
dt[fracon.glu>0.999,cor:=0.999]
dt[,val:=logit(cor)]
dt1<-dt
dt2<-TukeyHSD(aov(lm(val~(aa+bb+cc+dd)^4, dt1)))$`aa:bb:cc:dd`
rownames(dt2)<-gsub('4-','4.',gsub('S-','S.',gsub('\\:','\\+',rownames(dt2))))
pFraconGlu<-TukeyInteractionFlip2(dt2,plot=T)
DTFraconGlu<-TukeyInteractionFlip2(dt2,plot=F,degfreedom=14,tukeyPoint=F,fdr= 0.05)

colsA<-paste0('A_',c('aa','bb','cc','dd'))
colsB<-paste0('B_',c('aa','bb','cc','dd'))
DTFraconGlu[,paste0('A_',c('aa','bb','cc','dd')):=colsplit(A,'\\+',paste0('A_',c('aa','bb','cc','dd')))]
DTFraconGlu[,paste0('B_',c('aa','bb','cc','dd')):=colsplit(B,'\\+',paste0('B_',c('aa','bb','cc','dd')))]
DTFraconGlu[A_aa=='GAL4.WT'&B_aa=='GAL4.WT'&A_cc=='GAL80.07'&B_cc=='GAL80.07']

boo<-table(data.table(DTFraconGlu)$sig2)
boo[1]/(sum(boo))

# fracon.gal
dt<-dSub[,.(fracon.gal,aa,bb,cc,dd)]
dt[,cor:=fracon.gal]
dt[fracon.gal <0.001,cor:=0.001]
dt[fracon.gal>0.999,cor:=0.999]
dt[,val:=logit(cor)]
dt1<-dt
dt2<-TukeyHSD(aov(lm(val~(aa+bb+cc+dd)^4, dt1)))$`aa:bb:cc:dd`
rownames(dt2)<-gsub('4-','4.',gsub('S-','S.',gsub('\\:','\\+',rownames(dt2))))
pFraconGal<-TukeyInteractionFlip2(dt2)

# yfp.mean.glu
dt<-dSub[,.(val=yfp.mean.glu,aa,bb,cc,dd)]
dt1<-dt
dt2<-TukeyHSD(aov(lm(val~(aa+bb+cc+dd)^4, dt1)))$`aa:bb:cc:dd`
rownames(dt2)<-gsub('4-','4.',gsub('S-','S.',gsub('\\:','\\+',rownames(dt2))))
pYFPmeanGlu<-TukeyInteractionFlip2(dt2)

# yfp.mean.gal
dt<-dSub[,.(val=yfp.mean.gal,aa,bb,cc,dd)]
dt1<-dt
dt2<-TukeyHSD(aov(lm(val~(aa+bb+cc+dd)^4, dt1)))$`aa:bb:cc:dd`
rownames(dt2)<-gsub('4-','4.',gsub('S-','S.',gsub('\\:','\\+',rownames(dt2))))
pYFPmeanGal<-TukeyInteractionFlip2(dt2)

# growth.rate
dt<-dSub[,.(val= growth.rate,aa,bb,cc,dd)]
dt1<-dt
dt2<-TukeyHSD(aov(lm(val~(aa+bb+cc+dd)^4, dt1)))$`aa:bb:cc:dd`
rownames(dt2)<-gsub('4-','4.',gsub('S-','S.',gsub('\\:','\\+',rownames(dt2))))
pGrowthRate<-TukeyInteractionFlip2(dt2)

# fold.induction
dt<-dSub[,.(val=yfp.mean.gal/yfp.mean.glu,aa,bb,cc,dd)]
dt1<-dt
dt2<-TukeyHSD(aov(lm(val~(aa+bb+cc+dd)^4, dt1)))$`aa:bb:cc:dd`
rownames(dt2)<-gsub('4-','4.',gsub('S-','S.',gsub('\\:','\\+',rownames(dt2))))
pFoldInduction<-TukeyInteractionFlip2(dt2)


###################################################
### sandbox
###################################################
boo<-wt$mean-backg$mean
mut<-mean(phenDT[aa=='GAL4-L868P'&bb=='GAL3.delta'&cc=='GAL80S-1'&dd=='GALK.Sac_cer']$growth.rate)-backg$mean
mut/boo
mut<-mean(phenDT[aa=='GAL4-L868K'&bb=='GAL3.delta'&cc=='GAL80S-1'&dd=='GALK.Sac_cer']$growth.rate)-backg$mean
mut/boo
mut<-mean(phenDT[aa=='GAL4-L868G'&bb=='GAL3.delta'&cc=='GAL80S-1'&dd=='GALK.Sac_cer']$growth.rate)-backg$mean
mut/boo


mut<-mean(phenDT[aa=='GAL4-L868P'&bb=='GAL3.WT'&cc=='GAL80S-1'&dd=='GALK.Esc_col']$growth.rate)-backg$mean
mut/boo
mut<-mean(phenDT[aa=='GAL4-L868K'&bb=='GAL3.WT'&cc=='GAL80S-1'&dd=='GALK.Esc_col']$growth.rate)-backg$mean
mut/boo
mut<-mean(phenDT[aa=='GAL4-L868G'&bb=='GAL3.WT'&cc=='GAL80S-1'&dd=='GALK.Esc_col']$growth.rate)-backg$mean
mut/boo

