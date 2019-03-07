

# Scripts for dealing with flow cytometry data.
source.code.path.file<-'/Users/anew/Downloads/HarmoniousCombinations-master/CustomFunctions.R'
aesthetics.dir<-'/Users/anew/Downloads/HarmoniousCombinations-master/Aesthetics_Harmonious_Combinations.R'
head.dir<-"/Users/anew/Downloads/HarmoniousCombinations-master/Experiment1"
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/publication/'
outputdata<-'./output.data/'
layoutdir<-'./layout.sample.data/'
date<-'170927' # beginning date clones were measured
setwd(head.dir)
source(source.code.path.file)
source(aesthetics.dir)
br<-seq(-0.5,2.7,length=61)
sessionInfo()
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
DT0<-fread(paste0(outputdata,'17092x-Dataset.txt'),sep='\t',header=T)

phenocols<-c('growth.rate','phenotypic.index','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal','cell.per.ml.glu','cell.per.ml.gal')
ensembles<-c('plasgeno','binarycount','totalcopies')
dens.cols<-grepincols(DT0,'dens')
factorDefault<-c('orange1','cornflowerblue','darkgreen','red','chartreuse3')
DT0[,c('a','b','c'):=list(1-total.copies.GAL4,1-total.copies.GAL3,1-total.copies.GAL80)]
DT0[,WT:=plasgeno=='0 0 0 1 1 1']
genos<-c(ensembles,letters[1:3],'WT')
phenDT<-DT0[plasgeno=='0 0 0 1 1 1'|genomecopies=='0 0 0',c(genos,phenocols),with=F]
densDT<-DT0[plasgeno=='0 0 0 1 1 1'|genomecopies=='0 0 0',c(genos,dens.cols),with=F]


###################################################
### basic plots: generations with time, expression in glucose and galactose
###################################################


################
### Figure 1 show basic clusters
################

denscolconv<-function(DT, vnull =T){
	DT[,variable:=gsub('dens.','',variable)]
	DT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
	DT[,sugar:=paste0('g',sugar)]
	if(vnull ==T)DT[,variable:=NULL]
}

dum<-phenDT[,mean(growth.rate,n_distincta.rm=T),by=c(letters[1:3])]
clusts<-c('uninducible','inducible','constitutive')
dum[,cluster_new:=factor(c('uninducible','inducible','uninducible','uninducible','constitutive','uninducible','constitutive','uninducible'),levels=clusts)];dum[,V1:=NULL]
DT1<-merge(densDT, dum,by=c('a','b','c'))
mDT<-melt(DT1[,c('cluster_new',dens.cols),with=F],id='cluster_new')
denscolconv(mDT)
mDT[,dummyline:=paste0('z',sugar)]
mDT[,pheno:='fluorescence']
ylimits<-c(0,0.3)
ybreaks<-c(0.0,0.1,0.2,0.3)
ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
bycols<-c('cluster_new','dummyline','AU.FL','sugar')
mDT1a<-mDT[,mean(value,na.rm=T),by=bycols];setnames(mDT1a,'V1','density_mean')
mDT1b<-mDT[,sd(value,na.rm=T),by= bycols];setnames(mDT1b,'V1','density_sd')
mDT1<-merge(mDT1a,mDT1b,by=bycols)
mDT1[density_mean>=max(ylimits),density_mean:=max(ylimits)][order(sugar,AU.FL,decreasing=T)]
sample(LETTERS,4)


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

pUEHX<-ExpDensLinePlot(mDT1,ylimit=c(0,0.35))+ facet_grid(~ cluster_new)

if(writeplot==T){
	w=2.85;h=1.5
	ggsave(paste0(figout,'17092x-181116-pUEHX-CorePhenotypes-ExpressionDistributionLinePlots.pdf'), pUEHX,width=w,height=h)
}



###################################################
### correlations genotypes
###################################################
genos<-c(ensembles,letters[1:3],'WT')

summ<-summary.func.all(phenDT,c('growth.rate','fracon.glu','fracon.gal'),genos)
summ[,ord:=apply(data.frame(a,b,c),1,sum)]
summ[,`Mutation order`:=factor(ord)]
summ[WT==T,`Mutation order`:='WT control']
genos2<-c(genos,'Mutation order','ord')
mDT<-melt(summ,id=genos2)
mDT[,c('pheno','stat'):=colsplit(variable,'_',c('pheno','stat'))][,variable:=NULL]
form1<-castform(c(genos2,'pheno'),c('stat'))
cDT<-dcast(mDT, plasgeno+binarycount+totalcopies+a+b+c+WT+`Mutation order`+ord+pheno~stat)

cx<-merge(cDT[pheno=='fracon.glu'], cDT[pheno=='fracon.gal'],by=genos2)
cx[,facx:='3a']
cx[,facy:='a']

sample(letters,4)

pFonGluVsFonGal<-ggplot(cx[order(ord)],aes(x=mean.x,y=mean.y,xmax=upper.x,xmin=lower.x,ymax=upper.y,ymin=lower.y))+
	theme_NewPub+
	geom_abline(col='grey70')+
	geom_errorbarh(col='grey70')+geom_errorbar(col='grey70')+geom_point(aes(col=`Mutation order`), size=2,shape=21,stroke=0.5)+
	factorDefault2+
	theme(strip.text.x=element_blank(),strip.text.y=element_blank())+xlim(c(0,1))+ylim(c(0,1))+
	xlab('fraction ON in glucose')+ylab('fraction ON in galactose')+
	theme(legend.position='none')

dx<-merge(cDT[pheno=='fracon.glu'], cDT[pheno=='growth.rate'],by=genos2)
dx[,facx:='2a']
dx[,facy:='a']
sample(letters,4)

pFonGluVsGrowthRate <-ggplot(dx[order(ord)],aes(x=mean.x,y=mean.y,xmax=upper.x,xmin=lower.x,ymax=upper.y,ymin=lower.y))+
	theme_NewPub+
	geom_errorbarh(col='grey70')+geom_errorbar(col='grey70')+geom_point(aes(col=`Mutation order`), size=2,shape=21,stroke=0.5)+
	factorDefault2+geom_smooth(method='lm',size=0.5,fill='cornflowerblue',alpha=0.1)+
	theme(strip.text.x=element_blank(),strip.text.y=element_blank())+xlim(c(0,1))+ylim(c(0,0.3))+
	xlab('fraction ON in glucose')+ ylab(expression(paste('growth rate ',mu,' ',hr^-1)))+
	theme(legend.position='none')

ex<-merge(cDT[pheno=='fracon.gal'], cDT[pheno=='growth.rate'],by=genos2)
ex[,facx:='1a']
ex[,facy:='a']
sample(letters,4)

pFonGalVsGrowthRate <-ggplot(ex[order(ord)],aes(x=mean.x,y=mean.y,xmax=upper.x,xmin=lower.x,ymax=upper.y,ymin=lower.y))+
	theme_NewPub+
	geom_errorbarh(col='grey70')+geom_errorbar(col='grey70')+geom_point(aes(col=`Mutation order`), size=2,shape=21,stroke=0.5)+
	factorDefault2+geom_smooth(method='lm',size=0.5,fill='cornflowerblue',alpha=0.1)+
	theme(strip.text.x=element_blank(),strip.text.y=element_blank())+xlim(c(0,1))+ylim(c(0,0.3))+
	xlab('fraction ON in galactose')+ylab(expression(paste('growth rate ',mu,' ',hr^-1)))+
	theme(legend.position='none')


pFonGluVsFonGal
pFonGluVsGrowthRate
pFonGalVsGrowthRate


library(ggplot2)   
library(gtable)    
library(grid)
library(gridExtra) 
gA<-ggplotGrob(pFonGluVsFonGal)
gB<-ggplotGrob(pFonGluVsGrowthRate)
gC<-ggplotGrob(pFonGalVsGrowthRate)
# Set the widths
gA$widths <- gB$widths<-gC$widths

# Arrange the charts.
grid.newpage()
pNWHA<-list(gA, gB,gC,ncol=3)

summary(lm(mean.x~mean.y,cx))
summary(lm(mean.x~mean.y,dx))
summary(lm(mean.x~mean.y,ex))

pNWHA

if(writeplot==T){
	w=4.8;h=1.5
	ggsave(paste0(figout,'17092x-181116-pNWHA-FraconGluVsFraconGal.pdf'), do.call(grid.arrange,pNWHA),width=w,height=h)
}


