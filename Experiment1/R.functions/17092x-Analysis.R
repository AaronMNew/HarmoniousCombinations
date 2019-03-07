# Scripts for dealing with flow cytometry data.
source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'
source.code.path.file2<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis/R.functions/1707xx-180704-ReanalysisCustomScripts.R'
head.dir<-"/Users/anew/Google Drive/002_lehnerlabData/17092x-CopyNumberExperiment"
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/'
outputdata<-'./output.data/'
layoutdir<-'./layout.sample.data/'
date<-'170927' # beginning date clones were measured
setwd(head.dir)
source(source.code.path.file)
source(source.code.path.file2)
br<-seq(-0.5,2.7,length=61)

# this file generates these objects which can be used in downstream analysis
# merged.final<-fread(paste(date,"-DatasetWithGenerationsAndGeneExpression.txt",sep=''),sep='\t',header=T)


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
phenDT<-DT0[plasgeno=='0 0 0 1 1 1'|genomecopies=='0 0 0',c(genos,phenocols),with=F]
densDT<-DT0[plasgeno=='0 0 0 1 1 1'|genomecopies=='0 0 0',c(genos,dens.cols),with=F]


###################################################
### basic plots: generations with time, expression in glucose and galactose
###################################################


################
### show basic clusters
################

denscolconv<-function(DT, vnull =T){
	DT[,variable:=gsub('dens.','',variable)]
	DT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
	DT[,sugar:=paste0('g',sugar)]
	if(vnull ==T)DT[,variable:=NULL]
}

dum<-phenDT[,mean(growth.rate,na.rm=T),by=c(letters[1:3])]
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
mDT1<-mDT[,mean(value,na.rm=T),by=c('cluster_new','dummyline','AU.FL','sugar')];setnames(mDT1,'V1','density_mean')
mDT1[density_mean>=max(ylimits),density_mean:=max(ylimits)][order(sugar,AU.FL,decreasing=T)]
sample(LETTERS,4)

pXIMG<-ggplot(mDT1,aes(AU.FL,density_mean))+theme_minimal()+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(~ cluster_new)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
	theme(legend.title=element_blank())+
	theme(strip.text.x = element_text(size = 10))

if(writeplot==T){
	w=6;h=1.75
	ggsave(paste0(figout,'17092x-pXIMG-CorePhenotypes-ExpressionDistributionLinePlots.png'), pXIMG,width=w,height=h)
}



############
### expression density line plots of all genotypes
############

denscolconv<-function(DT, vnull =T){
	DT[,variable:=gsub('dens.','',variable)]
	DT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
	DT[,sugar:=paste0('g',sugar)]
	if(vnull ==T)DT[,variable:=NULL]
}


DTX<-copy(densDT)
mDT<-melt(DTX,id=genos)
denscolconv(mDT)
mDT[,dummyline:=paste0('z',sugar)]
mDT[,pheno:='fluorescence']

ylimits<-c(0,0.3)
ybreaks<-c(0.0,0.1,0.2,0.3)
ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
mDT[,density_mean:=mean(value),by=c(genos,'AU.FL','sugar')]
mDT1<-mDT[density_mean>=max(ylimits),density_mean:=max(ylimits)][order(sugar,AU.FL,decreasing=T)]
sample(LETTERS,4)
mDT1[,ord:=apply(data.frame(a,b,c),1,sum)]
mDT1[,c('A','B','C'):=lapply(.SD,function(x)gsub('0','(+)',gsub('1','(-)',x))),.SDcols=c('a','b','c')]
mDT1[,genotype:=applyPaste(data.frame(A,B,C,ord),'\n')]
mDT1[WT==T,genotype:='GAL4:\nGAL3:\nGAL80:\nMutation order']
mDT1[WT==T,ord:=0]
dum<-mDT1[,unique(ord),by='genotype'][order(V1)]
mDT1[,genotype:=factor(genotype,levels=dum$genotype)]


pKJFC<-ggplot(mDT1,aes(AU.FL,density_mean))+geom_point(data=mDT1,aes(AU.FL,value,col=sugar),size=0.1,alpha=0.3,shape=21)+theme_minimal()+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(~ genotype)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
	theme(strip.text.x = element_text(size = 10))


DTX2<-copy(phenDT)
mDT2<-melt(DTX2[,c(genos,grepincols(DTX2,c('cell.per.ml'))),with=F],id=genos)
mDT2[,c('boo','sugar'):=colsplit(variable,'ml\\.',c('boo','sugar'))][,c('variable','boo'):=list(NULL,NULL)]
mDT2[sugar=='glu',value:=value*9/150]
mDT2[,c('dummyline','density_mean','AU.FL'):=list(NA,NA,NA)]
#mDT2[,value:=mean]
mDT2[,pheno:='cell density']

mDT2[,ord:=apply(data.frame(a,b,c),1,sum)]
mDT2[,c('A','B','C'):=lapply(.SD,function(x)gsub('0','(+)',gsub('1','(-)',x))),.SDcols=c('a','b','c')]
mDT2[,genotype:=applyPaste(data.frame(A,B,C,ord),'\n')]
mDT2[WT==T,genotype:='GAL4:\nGAL3:\nGAL80:\nMutation order']
mDT2[WT==T,ord:=0]
dum<-mDT2[,unique(ord),by='genotype'][order(V1)]
mDT2[,genotype:=factor(genotype,levels=dum$genotype)]
mDT2[,generations:=log2(value),by=c(genos,'genotype','sugar')]
summ1<-mDT2[,summary.func(generations),by=c(genos,'genotype','sugar')]
summ2<-merge(summ1,summ1[sugar=='glu',.(init=mean,genotype)],by='genotype')
summ2[,mean:=mean-init]
summ2[,upper:=upper-init]
summ2[,lower:=lower-init]
summ2[sugar=='glu',time:=0]
summ2[sugar=='gal',time:=12]
sample(LETTERS,4)
pTWCM<-ggplot(summ2,aes(x=time,y=mean,ymax=upper,ymin=lower))+theme_minimal()+
	geom_line(col='grey70')+geom_errorbar(col='black',size=1)+geom_point(aes(col=sugar),shape=21,size=3,stroke=1)+
	facet_grid(~ genotype)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('generations') + xlab('time, hours')+
	scale_x_continuous(limits=c(-0.3,12.3),breaks=c(0,12),labels=c(0,12))+	scale_y_continuous(limits=c(-0.3,4),breaks=c(0,1,2,3,4),labels=c('0.0','1.0','2.0','3.0','4.0'))+
	theme(strip.text.x = element_blank(),panel.grid.minor=element_blank())	

ROWS<-1000;COLS<-ROWS
lmatblank<-data.frame(matrix(nrow=ROWS,ncol=COLS))
a<-0.6
b<-0.015
X<-1
Y<-1
lmatblank[1:(a*ROWS),1:(COLS*X)]<-1
lmatblank[(a*ROWS):ROWS,1:(COLS*X)]<-2
#lmatblank[(a*ROWS):ROWS,1:(COLS*b)]<-3
pBLANK<-ggplot()
lmat<-as.matrix(lmatblank)            

pFJCZ<-list(pKJFC,pTWCM,layout_matrix=lmat)
do.call(grid.arrange,pFJCZ)

if(writeplot==T){
	w<-10.4;h<-3.2
	ggsave(paste0(figout,date,'-180828-pFJCZ-ExpressionDistributionsAndGenerationsGrowth-SimpleSingleDoubleAndTriples.png'), do.call(grid.arrange,pFJCZ),width=w,height=h)
}



###################################################
### correlations genotypes
###################################################

DTX<-phenDT
mod1<-lm(growth.rate~fracon.gal+fracon.glu,phenDT)
# # bootstrap
# varexps<-(sapply(1:1000,function(i){
	# x<-sample_n(phenDT,replace=T,nrow(phenDT))
	# varexplained(x$growth.rate,predict(lm(growth.rate~fracon.gal+fracon.glu,x)))
# }))
summary(mod1)
varexpbyexpression<-mean(varexps)
sd(varexps)
varexpbymean<-DTm[,varexplained(growth.rate,growth.rate_mean)]
varexpbyexpression/varexpbymean

summary(mod1)$r.squared
summ<-summary.func.all(phenDT,c('growth.rate','fracon.glu','fracon.gal'),genos)
DTm<-merge(DTX,summ,by=genos)
varexpbymean<-DTm[,varexplained(growth.rate,growth.rate_mean)]
summ[,c('fracon.glu','fracon.gal'):=list(fracon.glu_mean,fracon.gal_mean)]
summ[,pred.expression_mean:=predict(mod1,newdata=data.frame(fracon.glu,fracon.gal))]
summ[,pred.expression_se:=predict(mod1,newdata=data.frame(fracon.glu,fracon.gal),se.fit=T)$se]
summ[,pred.expression_upper:=pred.expression_mean-pred.expression_se]
summ[,pred.expression_upper:=pred.expression_mean-pred.expression_se]
factorDefault<-c('orange1','cornflowerblue','darkgreen','red','chartreuse3')

summ[,ord:=apply(data.frame(a,b,c),1,sum)]
summ[,`Mutation order`:=factor(ord)]
summ[WT==T,`Mutation order`:='WT control']
# ggplot(summ[order(ord)],aes(fracon.glu_mean,growth.rate_mean,xmin=fracon.glu_lower,xmax=fracon.glu_upper,ymin=growth.rate_lower,ymax=growth.rate_lower))+
	# geom_errorbarh(col='grey70')+geom_errorbar(col='grey70')+geom_point(aes(col=`Mutation order`),shape=21,size=3,stroke=1)+scale_colour_manual(values=factorDefault)
	
# ggplot(summ[order(ord)],aes(fracon.gal_mean,growth.rate_mean,xmin=fracon.gal_lower,xmax=fracon.gal_upper,ymin=growth.rate_lower,ymax=growth.rate_lower))+
	# geom_errorbarh(col='grey70')+geom_errorbar(col='grey70')+geom_point(aes(col=`Mutation order`),shape=21,size=3,stroke=1)+scale_colour_manual(values=factorDefault)

# ggplot(summ,aes(fracon.glu_mean,fracon.gal_mean,xmin=fracon.glu_lower,xmax=fracon.glu_upper,ymin=fracon.gal_lower,ymax=fracon.gal_lower))+
	# geom_errorbarh(col='grey70')+geom_errorbar(col='grey70')+geom_point(aes(col=`Mutation order`),shape=21,size=3,stroke=1)+scale_colour_manual(values=factorDefault)+
	# geom_abline(col='grey70')

genos2<-c(genos,'Mutation order','ord')
mDT<-melt(summ,id=genos2)
mDT[,c('pheno','stat'):=colsplit(variable,'_',c('pheno','stat'))][,variable:=NULL]
form1<-castform(c(genos2,'pheno'),c('stat'))
cDT<-dcast(mDT, plasgeno+binarycount+totalcopies+a+b+c+WT+`Mutation order`+ord+pheno~stat)
a<-merge(cDT[pheno=='fracon.glu'],cDT[pheno=='growth.rate'],by=genos2)
a[,facx:='a']
a[,facy:='a']

d<-merge(cDT[pheno=='pred.expression'],cDT[pheno=='growth.rate'],by=genos2)
d[,facx:='a']
d[,facy:='c']


b<-merge( cDT[pheno=='fracon.gal'],cDT[pheno=='growth.rate'],by=genos2)
b[,facx:='a']
b[,facy:='b']

cx<-merge(cDT[pheno=='fracon.glu'], cDT[pheno=='fracon.gal'],by=genos2)
cx[,facx:='b']
cx[,facy:='a']

DTr<-rbind(a,b,d)

sample(letters,4)
DTr[,rsq:=summary(lm(mean.y~mean.x))$r.squared,by=c('facx','facy')]
labelframe<-DTr[,paste0('R2 = ',unique(round(rsq,3))),by=c('facx','facy')]
labelframe[,c('x','y'):=list(c(0.2,0.2,0.075),c(0.3,0.3,0.3))]
pKQSY<-ggplot(DTr[order(ord)],aes(x=mean.x,y=mean.y,xmax=upper.x,xmin=lower.x,ymax=upper.y,ymin=lower.y))+theme_minimal()+
	geom_errorbarh(col='grey70')+geom_errorbar(col='grey70')+geom_point(aes(col=`Mutation order`),shape=21,size=3,stroke=1)+
	scale_colour_manual(values=factorDefault)+
	theme(strip.text.x=element_blank(),strip.text.y=element_blank())+
	xlab('')+ylab('')+geom_smooth(method='lm',fill='grey90',colour='grey70')+
	geom_text(data=labelframe,aes(x,y,label=V1),inherit.aes=F,size=3)+theme(legend.position='none')+
	facet_grid(facx~facy,scales='free')

pLOGK<-ggplot(cx[order(ord)],aes(x=mean.x,y=mean.y,xmax=upper.x,xmin=lower.x,ymax=upper.y,ymin=lower.y))+theme_minimal()+
	geom_abline(col='grey70')+
	geom_errorbarh(col='grey70')+geom_errorbar(col='grey70')+geom_point(aes(col=`Mutation order`),shape=21,size=3,stroke=1)+
	scale_colour_manual(values=factorDefault)+
	theme(strip.text.x=element_blank(),strip.text.y=element_blank())+xlim(c(0,1))+ylim(c(0,1))+
	xlab('fraction ON in glucose')+ylab('fraction ON in galactose')

if(writeplot==T){
	w<-6.4;h<-2.2
	ggsave(paste0(figout,date,'-180828-pKQSY-GrowthRateVsFracONs.png'), pKQSY,width=w,height=h)
	w<-3.8;h<-2.2
	ggsave(paste0(figout,date,'-180828-pLOGK-FraconGluVsFraconGal.png'), pLOGK,width=w,height=h)
	
}



###################################################
### prediction of phenotypes
###################################################

pheno<-'growth.rate'
genos<-c(letters[1:3],'binarycount')
dSub<-phenDT[,c(letters[1:3],'binarycount',pheno),with=F]
setnames(dSub,pheno,'phenox')
backg<-dSub[a=='1',summary.func(phenox)]
backg0<-0#backg$mean
varx<-10^-16
credmin<-backg0+backg$sd*1.96
dSub[,x:=phenox]
dSub[phenox<credmin,x:=backg0]
dSub[,lx:=log(x-backg0)]
summ<-dSub[,summary.func(phenox,pheno),by=genos]
f1<-as.formula(paste0('lx~a+b+c'))
f2<-as.formula(paste0('lx~(a+b+c)^2'))
f3<-as.formula(paste0('lx~(a+b+c)^3'))
f4<-as.formula(paste0('lx~a+b*c'))
forms<-list(f1,f2,f3,f4)
mods<-lapply(forms,function(formx){
	lm(formx,data=dSub)
})
preds<-lapply(mods,function(modx){
	predict.lm(modx,newdata= summ,se.fit=T)
})
pred_means<-data.table(t(ldply(lapply(preds,function(predx)exp(predx$fit)+backg0))))
colnames(pred_means)<-paste0(colnames(pred_means),'_mean')
logSEConv<-function(predx){
	# have to have a lm() model's prediction with se.fit=T
	exp(predx$fit)*predx$se.fit
	}
pred_ses<-data.table(t(ldply(lapply(preds,function(predx)logSEConv(predx)))))
colnames(pred_ses)<-paste0(colnames(pred_ses),'_se')
out<-data.table(summ[,c(genos,grepincols(summ,c('_mean','_se'))),with=F],pred_means,pred_ses)
out[,ord:=as.factor(apply(data.frame(a,b,c),1,sum))]
# mDT<-melt(out[,c(grepincols(out,'V'),'pheno',genos,'ord'),with=F],id=c('pheno',genos,'ord'))
mDT<-melt(out,id=c(genos,'ord'))
mDT[,c('pred','stat'):=colsplit(variable,'_',c('pred','stat'))][,variable:=NULL]
form<-castform(leftside=c(genos,'ord','pred'),rightside=c('stat'))
cDT1<-dcast(mDT,form,value.var='value')
cDT<-merge(cDT1[pred!=pheno],cDT1[pred==pheno],by=c(genos,'ord'))
levs<-c('first-order\nsingles','second-order\nsingles + doubles','third-order\nsingles + doubles + triple','first-order\nsingles + GAL80:GAL3 double')

dum<-data.table(pred.x=paste0('V',1:length(forms)),fac=factor(levs,levels=levs))
DTm<-merge(cDT,dum,by='pred.x')
fac2<-applyPaste(DTm[,percent(round(varexplained(mean.y,mean.x),2)),by=fac],'\n')
DTm[,fac2:=factor(paste0(fac,'\n',percent(round(varexplained(mean.y,mean.x),2))),levels=fac2),by='fac']

limits<-c(0,0.35)
ptotal<-ggplot(DTm,aes(mean.x,mean.y,xmax=mean.x+se.x,xmin=mean.x-se.x,ymax=mean.y+se.y,ymin=mean.y-se.y))+theme_minimal()+
	geom_abline(col='grey70')+
	geom_errorbarh(col='grey70')+geom_errorbar(col='grey70')+geom_point(aes(col=ord),shape=21,size=2)+
	theme(strip.text.x=element_text(size=8))+
#	theme(legend.position='none')+
	scale_colour_manual(values=factorDefault)+
	ylim(limits)+xlim(limits)+
	ylab(expression(paste('observed growth rate ',mu,' ',hr^-1)))+
	xlab(expression(paste('predicted growth rate ',mu,' ',hr^-1)))+
	facet_grid(~fac2)

pQXHT<-ggplot(DTm[fac!='third-order\nsingles + doubles + triple'],aes(mean.x,mean.y,xmax=mean.x+se.x,xmin=mean.x-se.x,ymax=mean.y+se.y,ymin=mean.y-se.y))+theme_minimal()+
	geom_abline(col='grey70')+
	geom_errorbarh(col='grey70')+geom_errorbar(col='grey70')+geom_point(aes(col=ord),shape=21,size=2)+
	theme(strip.text.x=element_text(size=8))+
	theme(legend.position='none')+
	scale_colour_manual(values=factorDefault)+
	ylim(limits)+xlim(limits)+
	ylab(expression(paste('observed growth rate ',mu,' ',hr^-1)))+
	xlab(expression(paste('predicted growth rate ',mu,' ',hr^-1)))+
	facet_grid(~fac2)

if(writeplot==T){
	w<-7.4;h<-2.6
	ggsave(paste0(figout,date,'-181009-pQXHT-MultiplicativeModelPrediction.png'), pQXHT,width=w,height=h)
}

ggplot(DTm,aes(mean.x,mean.y,xmax=mean.x+se.x,xmin=mean.x-se.x,ymax=mean.y+se.y,ymin=mean.y-se.y))+theme_minimal()+
	geom_abline(col='grey70')+
	geom_errorbarh()+geom_point(aes(col=ord),shape=21,size=2)+
	theme(strip.text.x=element_text(size=8))+
	scale_colour_manual(values=factorDefault)+
	facet_grid(~fac,scales='free')





