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
### custom ggplot themes
###################################################

factorDefault<-c('orange1','steelblue','darkgreen','red','chartreuse3')
factorTwoSugar<-c('orangered','deepskyblue')
heatmapColRange<-c('cornflowerblue','yellow','indianred')
heatmapColBounded<-c('white','cornflowerblue','yellow','indianred')
heatmapColCentered<-c('cornflowerblue','white','indianred')

cust_font_facet<-theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))
cust_font_axis<-theme(axis.text = element_text(size = 8))
cust_font_axislab<-theme(axis.title = element_text(size = 12))
cust_font_legend<-theme(legend.text=element_text(size=8),legend.title=element_text(size=10))

cust_font<-cust_font_facet+cust_font_axis+cust_font_axislab+ cust_font_legend

cust_panel<-theme(panel.background=element_blank(),panel.border=element_blank())
cust_strip<-theme(strip.background=element_rect(fill='grey90',colour='transparent'))
cust_axis<-theme(axis.line=element_line(arrow=arrow(angle=20,length=unit(0.15,'cm')),
	size=0.3,colour='grey50'),
	axis.ticks=element_blank())

theme_NewPub<-theme_classic()+cust_strip+ cust_axis + cust_font+ cust_panel

# booasdfawfgb<-data.table(x=1:10,y=1:10,fac=factor(rep(c('A','B'),5)))

# ggplot(booasdfawfgb,aes(x,y))+geom_point()+geom_line(aes(col=fac))+
	# theme_NewPub +
	# scale_colour_manual(values=factorDefault)+
	# facet_grid(~fac)

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
### quick overview of phenotypes
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

############
### tile plots of growth rate across all genotypes
############

genos<-c('aa','bb','cc','dd')
summ<-DT0[,summary.func(growth.rate),by=genos]
summ[,`growth rate 1/hr`:=mean]
sample(letters,4)

pDGMY<-ggplot(summ,aes(aa,cc))+theme_NewPub+
	geom_tile(aes(fill=`growth rate 1/hr`))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'))+
	facet_grid(dd~bb)+ theme(axis.text.x = element_text(angle = 30,margin=ggplot2::margin(0.4,10,-0.5,-1,'cm')))+xlab('GAL4 alleles')+ylab('GAL80 alleles')




pVEGX<-ggplot(summ,aes(aa,cc))+geom_tile(aes(fill=`growth rate 1/hr`))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'))+
	facet_grid(dd~bb)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))+xlab('GAL4 alleles')+ylab('GAL80 alleles')
if(writeplot==T){
	w<-7;h<-6.2
	ggsave(paste0(figout,date,'-180828-pVEGX-TilePlotsAllGenotypes-GrowthRate.png'), pVEGX,width=w,height=h)
	ggsave(paste0(figout,date,'-180828-pVEGX-TilePlotsAllGenotypes-GrowthRate.pdf'), pVEGX,width=w,height=h)
	
}

############
### expression density line plots of all genotypes
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
# Old -- switched to new R script
# # ExpDensLinePlot<-function(mDT1,ylimits=c(0,0.35)){
	# xlimits<-c(-0.3,3)
	# xbreaks<-c(0,1,2)
	# xlabs<-c('0','1','2')
	# ybreaks<-c(0.0,0.1,0.2,0.3)
	# ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
	# ggplot(mDT1[order(sugar)],aes(AU.FL,density_mean))+#geom_point(data=mDT1,aes(AU.FL,value,col=sugar),size=0.1,alpha=0.3,shape=21)+
		# theme_NewPub+
		# geom_ribbon(aes(x=AU.FL,ymax=density_mean+density_sd,ymin=density_mean-density_sd,fill=sugar),alpha=0.3)+
		# geom_line(aes(col=sugar))+
		# scale_colour_manual(values= factorTwoSugar,na.value='transparent')+
		# ylab('density') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
		# scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
		# scale_x_continuous(limits=xlimits,breaks=xbreaks,labels=xlabs)+
		# theme(legend.position='none')
# }

# G4<-c('GAL4.WT','GAL4-L868P','GAL4-L868K','GAL4-L868G') # G4<-'GAL4-L868P'
# G80<-c('GAL80.WT','GAL80S-1','GAL80.07')
# GK<-c('GALK.Sac_cer','GALK.Esc_col')
# G3<-c('GAL3.WT','GAL3.delta')
# mDT1<-mDT[aa%in%G4&cc%in%G80&dd%in%GK&bb%in%G3]
# sample(LETTERS,4)
# pAFMX<-ExpDensLinePlot(mDT1)+facet_grid(bb+dd~cc+aa)
# pAFMXLegend<-ExpDensLinePlot(mDT1)+facet_grid(bb+dd~cc+aa)+theme(legend.position='bottom')
# if(writeplot==T){
	# w<-11;h<-4.15
	# ggsave(paste0(figout,date,'-181109-pAFMX-ExpressionDistributions-SubsetForPub.png'), pAFMX,width=w,height=h)
	# ggsave(paste0(figout,date,'-181109-pAFMX-ExpressionDistributions-SubsetForPub.pdf'), pAFMX,width=w,height=h)
	# ggsave(paste0(figout,date,'-181109-pAFMX-ExpressionDistributions-WithLegend-SubsetForPub.pdf'), pAFMXLegend,width=w,height=h)
	# w<-8;h<-4
	# ggsave(paste0(figout,date,'-181109-pAFMX-ExpressionDistributions-WithLegend-SubsetForPub.pdf'), pAFMXLegend,width=w,height=h)	
# }

# pDIRG<-ExpDensLinePlot(mDT)+facet_grid(cc+aa~bb+dd)
# if(writeplot==T){
	# w<-8;h<-35
	# ggsave(paste0(figout,date,'-181109-pDIRG-ExpressionDistributions-AllGenotypes.pdf'), pDIRG,width=w,height=h, limitsize = FALSE)
# }


G4<-c('GAL4.WT','GAL4-L868P','GAL4-L868K','GAL4-L868G') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80S-1')
GK<-c('GALK.Sac_cer','GALK.Esc_col')
G3<-c('GAL3.WT','GAL3.delta')
mDT1<-mDT[aa%in%G4&cc%in%G80&dd%in%GK&bb%in%G3]
sample(LETTERS,4)
pWIPK<-ExpDensLinePlot(mDT1)+facet_grid(cc+dd~aa+bb)
pWIPK
if(writeplot==T){
	w<-8.9;h<-4.3
	ggsave(paste0(figout,date,'-180828-pWIPK-ExpressionDistributions-SubsetIncludingNewInducible.png'), pWIPK,width=w,height=h)
}
mDT1<-mDT[density_mean>=max(ylimits),density_mean:=max(ylimits)][order(sugar,AU.FL,decreasing=T)]
mDT1[,bb:=reverse_lev(bb)]
pYDTE<-ggplot(mDT1,aes(AU.FL,density_mean))+geom_point(data=mDT1,aes(AU.FL,value,col=sugar),size=0.1,alpha=0.3,shape=21)+theme_minimal()+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(cc+dd~aa+bb)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density of observations') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
	theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))
if(writeplot==T){
	w<-14;h<-24
	ggsave(paste0(figout,date,'-180828-pYDTE-ExpressionDistributions-FullDataset.png'), pYDTE,width=w,height=h)
}

G4<-c('GAL4.WT','GAL4-L868G') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80S-1')
GK<-c('GALK.Sac_cer','GALK.Esc_col')
G3<-c('GAL3.WT','GAL3.delta')
ybreaks<-c(0.0,0.1,0.2,0.3)
ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
mDT1<-mDT[density_mean>=max(ylimits),density_mean:=max(ylimits)][order(sugar,AU.FL,decreasing=T)][aa%in%G4&cc%in%G80&dd%in%GK&bb%in%G3]
sample(LETTERS,4)
mDT1[,bb:=reverse_lev(bb)]
pYRXI<-ggplot(mDT1,aes(AU.FL,density_mean))+geom_point(data=mDT1,aes(AU.FL,value,col=sugar),size=0.1,alpha=0.3,shape=21)+theme_minimal()+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(cc+dd~aa+bb)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density of observations') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
	cust_font+
	theme(legend.position='none')
if(writeplot==T){
	w<-4.4;h<-4.2
	ggsave(paste0(figout,date,'-180828-pYRXI-ExpressionDistributions-SmallSubsetIncludingNewInducible.png'), pYRXI,width=w,height=h)
}

G4<-c('GAL4.WT') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80.07')
GK<-c('GALK.Sac_cer','GALK.Esc_col')
G3<-c('GAL3.WT','GAL3.delta')
ybreaks<-c(0.0,0.1,0.2,0.3)
ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
mDT1<-mDT[density_mean>=max(ylimits),density_mean:=max(ylimits)][order(sugar,AU.FL,decreasing=T)][aa%in%G4&cc%in%G80&dd%in%GK&bb%in%G3]
sample(LETTERS,4)
mDT1[,bb:=reverse_lev(bb)]
pSMGP<-ExpDensLinePlot(mDT1)
if(writeplot==T){
	w<-2.9;h<-4.14
	ggsave(paste0(figout,date,'-180828-pSMGP-ExpressionDistributions-SmallSubsetIncludingGAL80.07.png'), pSMGP,width=w,height=h)
	ggsave(paste0(figout,date,'-180828-pSMGP-ExpressionDistributions-SmallSubsetIncludingGAL80.07.pdf'), pSMGP,width=w,height=h)

}



###################################################
### Scatterplot illustrating clones that can grow without GAL3
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
ttests[,`FDR < 0.05`:=as.character(padj<=0.05)]
cDT<-data.table(cDT1,ttests)
# cDT[aa=='GAL4.WT'&cc=='GAL80.WT',`FDR < 0.05`:='WT GAL80 and GAL4']
sample(letters,4)
cDT[,bb:=reverse_lev(bb)]
pZDSN<-ggplot(cDT[order(`FDR < 0.05`,decreasing=T)],aes(mean_GAL1,mean_GALK,xmin=lower_GAL1,xmax=upper_GAL1,ymin=lower_GALK,ymax=upper_GALK))+theme_minimal()+
	geom_errorbar(col='grey50')+geom_errorbarh(col='grey50')+geom_point(aes(col=`FDR < 0.05`),size=2,shape=21,stroke=1)+
	scale_colour_manual(values=c('black','red','cornflowerblue'))+
	geom_abline(col='grey70')+	
	geom_point(data=cDT[aa=='GAL4.WT'&cc=='GAL80.WT'], aes(mean_GAL1,mean_GALK),col='cornflowerblue',size=3,shape=21,stroke=1.5)+
	theme(legend.position='bottom')+
	ylab(expression(paste('GALK sensor (-) ',mu,' ',hr^-1)))+	xlab(expression(paste('GAL1 sensor (+) ',mu,' ',hr^-1)))+
	facet_grid(~bb)

if(writeplot==T){
	w=5;h=3
	ggsave(paste0(figout,'180404-181011-pZDSN-GAL1VsGALK-FacetByGAL3.png'), pZDSN,width=w,height=h)

}

aas<-unique(cDT[`FDR < 0.05`==T&abs(diff.of.means)>0.05]$aa)
ccs<-unique(cDT[`FDR < 0.05`==T& abs(diff.of.means)>0.05]$cc)


###################################################
### Scatterplot illustrating clones that can grow without GAL3
###################################################
denscolconv<-function(DT, vnull =T){
	DT[,variable:=gsub('dens.','',variable)]
	DT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
	DT[,sugar:=paste0('g',sugar)]
	if(vnull ==T)DT[,variable:=NULL]
}
aas<-c('GAL4.WT','GAL4-L868G')
ccs<-c('GAL80.WT','GAL80S-1')
mDT1<-melt(densDT,id=c(genos))[aa%in%aas&cc%in%ccs&bb=='GAL3.WT'&dd=='GALK.Sac_cer']
denscolconv(mDT1,vnull=T)
mDT<-mDT1[,density_mean:=mean(value,na.rm=T),by=c(notin(colnames(mDT1),'value'))]
mDT[,dummyline:=paste0('z',sugar)]
ylimits<-c(0,0.3)
mDT[density_mean>max(ylimits),density_mean:=max(ylimits)]
ybreaks<-c(0.0,0.1,0.2,0.3)
ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
sample(LETTERS,4)

pBCNO<-ggplot(mDT,aes(AU.FL,density_mean))+geom_point(data=mDT1,aes(AU.FL,value,col=sugar),size=0.1,alpha=0.3,shape=21)+theme_minimal()+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(cc~aa)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density of observations') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
	theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))


if(writeplot==T){
	w=4;h=3
	ggsave(paste0(figout,'180404-181011-pBCNO-NovelInduciblePhenotypes-LinePlotsExpressionDensity.png'), pBCNO,width=w,height=h)

}

###################################################
### cluster across glu,gal expression density with clone.named rows using hdbscan
###################################################

denscolconv<-function(DT, vnull =T){
	DT[,variable:=gsub('dens.','',variable)]
	DT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
	DT[,sugar:=paste0('g',sugar)]
	if(vnull ==T)DT[,variable:=NULL]
}

DTX<-copy(densDT[dd!='HIS5.Sch_pom'])#copy(phenDT)#
phenocols<-dens.cols#c('growth.rate')#
genocols<-genos
subcols<-c(phenocols,genos)
a<-DTX[,c(subcols),with=F]
bycols2<-c(genos)
# ord_within_cluster<-'growth.rate'
ord_within_cluster<-'phenotypic.index'
summ1<-na.exclude(a[,lapply(.SD,function(x)mean(log10(x+0.001),na.rm=T)),by=bycols2])

booDT<-phenDT[dd!='HIS5.Sch_pom']
PI4<-booDT[,median(phenotypic.index,na.rm=T),by='aa']
PI3<-booDT[,median(phenotypic.index,na.rm=T),by='bb']
PI80<-booDT[,median(phenotypic.index,na.rm=T),by='cc']
PIGK<-booDT[,median(phenotypic.index,na.rm=T),by='dd']

dt.list=list(PI4,PI3,PI80,PIGK)
boo1<-c('aa','bb','cc','dd')
lapply(1:length(dt.list),function(i){
	x<-dt.list[[i]]
	setnames(x,c('V1'),c(paste0(boo1[i],'_','PI_mean')))
})

DTm1<-merge_recurse3(qDT=summ1,dt.list=dt.list,by.list=boo1)
# pts<-4:20
# lapply(pts,function(x){
	# hdbtest<-hdbscan((grepcols(DTm1,'dens')),minPts=x)
	# DTm1[,cluster:=hdbtest$cluster]
	# c(paste0('minPts = ',x,' '),table(DTm1$cluster)/sum(table(DTm1$cluster)))
# })
# running above shows that 10 or 12 pts are best. 10 gives 4 clusters and 12 gives three, with ~10% uncategorized
hdbtest<-hdbscan((grepcols(DTm1,'dens')),minPts=9)
DTm1[,cluster:=hdbtest$cluster]
DTm1[,cluster_outlier_score:=hdbtest$outlier_score]

PI<-phenDT[,mean(phenotypic.index,na.rm=T),by=genos]
DTm <-merge(DTm1,PI,by=genos)
DTm[,clone.named:=applyPaste(data.frame(aa,bb,cc,dd),' ')]
PCl<-DTm[,summary.func(V1),by=cluster][order(mean)]
pointframe2<-data.table(DTm[cluster==0,'clone.named'],WT=apply(10^DTm[cluster==0,(dens.cols[15:60]),with=F],1,sum)<0.1)
DTm[cluster==2, cluster_new:='constitutive']
DTm[cluster==1, cluster_new:='inducible - new class']
DTm[cluster==0, cluster_new:='constitutive - new class']
DTm[cluster==3, cluster_new:='uninducible']
DTm[clone.named%in%pointframe2[WT==T]$clone.named,cluster_new:='inducible - WT']
clustlev<-c('constitutive','constitutive - new class','inducible - WT','inducible - new class','uninducible')
DTm[, cluster_new:=factor(cluster_new,levels=clustlev)]

############
### plot all clustered samples as logged density tile plots
############

lev<-DTm[,mean(V1),by=c(genos,'clone.named','cluster_new','cluster_outlier_score')][order(cluster_new, -V1,cluster_outlier_score)]
mDT1<-melt(DTm,id=c(bycols2,grepincols(DTm,'PI'),'clone.named','cluster_new','V1','cluster_outlier_score','cluster'))[order(cluster_new,V1)]
denscolconv(mDT1,vnull=T)
mDT<-mDT1
mDT[,Y:=factor(clone.named,levels=lev$clone.named)]
mDT[,X:=factor(round(AU.FL,3))]
mDT[sugar=='gal',X:=factor(round(AU.FL,3)+4)]
mDT[,`density, log10`:=value]
uniqfun<-function(x,fac=1,shift=1)(unique(x)[order(unique(x))][seq(shift,length(unique(x)),by=fac)])
collim<-range(mDT$value)
xlab<-c(uniqfun(mDT$AU.FL,fac=4,shift=3),uniqfun(round(mDT$AU.FL,3),shift=3,fac=4))
xbreak<-uniqfun(mDT$X,fac=4,shift=3)
sample(letters,4)
pointframe<-lev[,.(Y=clone.named,X=125, cluster_new=factor(cluster_new))]
sample(LETTERS, 4)
pJOPS<-ggplot(mDT,aes(X,Y))+geom_tile(aes(fill=`density, log10`))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),na.value='white',limits=collim)+
	theme(axis.text.x=element_text(angle=90,vjust=0.4,hjust=1,size=8))+
	theme(axis.text.y=element_text(size=6))+
	geom_vline(xintercept=61,col='black',size=2)+
	geom_vline(xintercept=61,col='white',size=1.5)+
	geom_point(data=pointframe,aes(X,Y,col= cluster_new))+
#	geom_point(data= pointframe2,aes(126,clone.named,col=WT))+
	scale_colour_manual(values=c('orange1','cornflowerblue','darkgreen','red','chartreuse3'),na.value='transparent')+
	scale_x_discrete(breaks=xbreak,labels=xlab)+
	theme(axis.title.y=element_blank(),axis.text.y=element_blank(),	axis.ticks.y=element_blank())+
	xlab('Expression in glucose               Exression in galactose   \n pseudo log10 A.U. Fluorescence units')+
	theme(axis.title=element_text(size=10,face="plain"))
if(writeplot==T){
	w<-5;h<-4
	ggsave(paste0(figout,'180404-180905-pJOPS-ExpressionDensityTilePlot-AllClones.png'), pJOPS,width=w,height=h)
}
############
### plot archetypal expression densities
############

mDT1<-copy(mDT)
mDT1[,density:=exp(value)-0.001]
setnames(mDT1,'V1','PhenotypeForRanking')
bycols<-c('cluster_new','sugar','AU.FL')
summ1<-mDT1[,mean(density),by=bycols]
summ2<-mDT[,.N,by=bycols]
summ2[,prop:=N/length(unique(mDT$clone.named))]
summ<-merge(summ1,summ2,by=bycols)
summ[,facet:=percent(round(prop,3))]
#summ[,cluster_new:=cluster]

# dumm<-summ[,unique(cluster_new),by='cluster']
mDTm<-merge(summ,mDT1,by=bycols)
mDTm[,V2:=mean(density),by=c(bycols)]
mDTm[,dummyline:=paste0('z',sugar)]
sample(LETTERS,4)
pUNQT<-ggplot(mDTm,aes(AU.FL,V2))+geom_point(data=mDTm,aes(AU.FL,density,col=sugar),size=0.1,alpha=0.3,shape=21)+
	geom_line(aes(col=dummyline),size=2)+geom_line(aes(col=sugar))+
	facet_grid(~cluster_new+facet)+scale_colour_manual(values=c('orange1','steelblue','grey90','grey90'),na.value='transparent')+
	ylab('density of observations') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
	theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))

if(writeplot==T){
	w<-8;h<-2.6
	ggsave(paste0(figout,'180404-180905-pUNQT-ExpressionDensityClusters-LinePlot-FacetByClass.png'), pUNQT,width=w,height=h)
}



###################################################
### measure effect of leaky expression on GAL pathway
###################################################

# summarize data
summ<-summary.func.all(DT0,c('growth.rate','fracon.glu','yfp.mean.glu','yfp.mean.gal','mean.on.fraction.glu'),bylist=genos)
summ[aa=='GAL4.WT'&bb=='GAL3.WT']
summ[,clone.named:=applyPaste(data.frame(aa,cc),' ')]
summ[,WT:=grepl(c('GAL4.WT GAL80.WT'),clone.named)]


############
### fit growth rate in galactose as a logistic function of GAL pathway expression in glucose
############

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


# merge summary and raw data
DTm<-merge(summ,dSub,by=genos)
DTm[,clone.named:=applyPaste(data.frame(aa,bb,cc,dd),' ')]
length(unique(DTm$clone.named))
DTm[,pred:=predict(mod,DTm[,.(x)])]
DTm[,estimate:=y-pred]
# ggplot(DTm,aes(pred,y,col=GALK.cat))+geom_point()+facet_wrap(~bb)+geom_abline()
tframe<-DTm[,t.test(estimate,alternative='greater'),by=genos]
tframe[,p.adj:=p.adjust(p.value)]
tframe[,`FDR < 0.05`:=p.adj<0.05]
tframe[,sig:=`FDR < 0.05`]
tframe[,dup:=applyPaste(tframe[,.(aa,bb,cc,dd)],'')]
DTm2<-merge(DTm,tframe[!duplicated(dup),c(genos,'sig','FDR < 0.05'),with=F],by=genos)

# plot of background GAL pathway expression in glucose vs growth rate for different sensor / GALK backgrounds 
pOYIC<-ggplot(DTm2[order(`FDR < 0.05`)],aes(yfp.mean.glu_mean,growth.rate_mean,xmax= yfp.mean.glu_upper,xmin=yfp.mean.glu_lower,ymax= growth.rate_upper,ymin=growth.rate_lower))+
#	geom_point(data=DTm,aes(x,y),shape=21,size=1,col='grey50',inherit.aes=F)+
	geom_errorbar(col='grey50')+geom_errorbarh(col='grey50')+geom_point(aes(col=`FDR < 0.05`),size=2,shape=21,stroke=1)+
	geom_line(data=lineframe,aes(x,y),inherit.aes=F)+ylab(expression(paste('mean ',mu,' ',hr^-1,' in galactose')))+xlab('mean expression of GAL pathway in glucose, A.U.')+
	facet_grid(GALK.cat~bb)+theme_minimal()+ theme(axis.text.x = element_text(angle = 30, hjust = 1))+scale_colour_manual(values=c('black','red'))+
	theme(strip.text.y = element_text(angle =0),legend.position='bottom')

# smaller plot of expression vs growth rate without HIS5 for main fig - strip text normal
pHZMVa<-ggplot(DTm2[order(`FDR < 0.05`)][dd!='HIS5.Sch_pom'],aes(yfp.mean.glu_mean,growth.rate_mean,xmax= yfp.mean.glu_upper,xmin=yfp.mean.glu_lower,ymax= growth.rate_upper,ymin=growth.rate_lower))+
#	geom_point(data=DTm,aes(x,y),shape=21,size=1,col='grey50',inherit.aes=F)+
	geom_errorbar(col='grey50')+geom_errorbarh(col='grey50')+geom_point(aes(col=`FDR < 0.05`),size=2,shape=21,stroke=1)+
	geom_line(data=lineframe,aes(x,y),inherit.aes=F)+ylab(expression(paste('mean ',mu,' ',hr^-1,' in galactose')))+xlab('mean expression of GAL pathway in glucose, A.U.')+
	facet_grid(GALK.cat~bb)+theme_minimal()+ theme(axis.text.x = element_text(angle = 30, hjust = 1))+scale_colour_manual(values=c('black','red'))+
	theme(legend.position='bottom')


DTemp<-copy(DTm2)
DTemp[,bb:=reverse_lev(bb)]
DTemp[,GALK.cat:=reverse_lev(GALK.cat,c(2,1,3))]
# smaller plot of expression vs growth rate without HIS5 for main fig - strip text normal
pHZMVb<-ggplot(DTemp[order(`FDR < 0.05`)][dd!='HIS5.Sch_pom'],aes(yfp.mean.glu_mean,growth.rate_mean,xmax= yfp.mean.glu_upper,xmin=yfp.mean.glu_lower,ymax= growth.rate_upper,ymin=growth.rate_lower))+
#	geom_point(data=DTm,aes(x,y),shape=21,size=1,col='grey50',inherit.aes=F)+
	geom_errorbar(col='grey50')+geom_errorbarh(col='grey50')+geom_point(aes(col=`FDR < 0.05`),size=2,shape=21,stroke=1)+
	geom_line(data=lineframe,aes(x,y),inherit.aes=F)+ylab(expression(paste('mean ',mu,' ',hr^-1,' in galactose')))+xlab('mean expression of GAL pathway in glucose, A.U.')+
	facet_grid(GALK.cat~bb)+theme_minimal()+ theme(axis.text.x = element_text(angle = 30, hjust = 1))+scale_colour_manual(values=c('black','red'))+
	theme(legend.position='bottom')


# smaller plot of expression vs growth rate without HIS5 for main fig - strip text rotated out
pHZMV<-ggplot(DTm2[order(`FDR < 0.05`)&dd!='HIS5.Sch_pom'],aes(yfp.mean.glu_mean,growth.rate_mean,xmax= yfp.mean.glu_upper,xmin=yfp.mean.glu_lower,ymax= growth.rate_upper,ymin=growth.rate_lower))+
#	geom_point(data=DTm,aes(x,y),shape=21,size=1,col='grey50',inherit.aes=F)+
	geom_errorbar(col='grey50')+geom_errorbarh(col='grey50')+geom_point(aes(col=`FDR < 0.05`),size=2,shape=21,stroke=1)+
	geom_line(data=lineframe,aes(x,y),inherit.aes=F)+ylab(expression(paste('mean ',mu,' ',hr^-1,' in galactose')))+xlab('mean expression of GAL pathway in glucose, A.U.')+
	facet_grid(GALK.cat~bb)+theme_minimal()+ theme(axis.text.x = element_text(angle = 30, hjust = 1))+scale_colour_manual(values=c('black','red'))+
	theme(strip.text.y = element_text(angle =0),legend.position='bottom')

collim<-c(-max(abs(range(DTm$estimate))),max(abs(range(DTm$estimate))))
DTm2[,`observed -\nexpected`:=estimate]
sample(letters,4)
# tile plot showing differences and significance between deviation of observed growth rate from expectation based on no-sensor background growth rate vs YFP in glucose
pUNPY<-ggplot(DTm2,aes(aa,cc))+geom_tile(aes(fill=`observed -\nexpected`))+scale_fill_gradientn(colours=c('cornflowerblue','white','indianred'),limits=round(collim,2))+
	facet_grid(GALK.cat~bb)+geom_point(aes(aa,cc,col=`FDR < 0.05`))+scale_colour_manual(values=c('transparent','red'))+theme_minimal()+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+xlab('GAL4 alleles')+ylab('GAL80 alleles')+
	theme(strip.text.y = element_text(angle =0),legend.position='bottom')


if(writeplot==T){
	w<-6;h<-3.5
	ggsave(paste0(figout,date,'-180828-pHZMV-GLUExpressionVsGrowthRate-DoseResponse-Small.png'), pHZMV,width=w,height=h)
	w<-6*0.8;h<-5*0.8
	ggsave(paste0(figout,date,'-180828-pHZMVa-GLUExpressionVsGrowthRate-DoseResponse-Small.png'), pHZMVa,width=w,height=h)
	w<-6*0.8;h<-5*0.8
	ggsave(paste0(figout,date,'-180828-pHZMVb-GLUExpressionVsGrowthRate-DoseResponse-Small.png'), pHZMVb,width=w,height=h)
	w<-5.3;h<-6
	ggsave(paste0(figout,date,'-180828-pOYIC-GLUExpressionVsGrowthRate-DoseResponse.png'), pOYIC,width=w,height=h)
	w<-5.3;h<-6
	ggsave(paste0(figout,date,'-180828-pUNPY-DeviationGrowthRateFromExpectationByYFPGLU-DoseResponse.png'), pUNPY,width=w,height=h)
}

# save your inducible classifications as a data.table
InducibleFrame<-tframe[!duplicated(dup),c(genos,'FDR < 0.05'),with=F];setnames(InducibleFrame,'FDR < 0.05','inducible')
save(InducibleFrame,file=paste0(outputdata,'180404-180904-InducibleGenotypesBasedOnGluExpression.rdata'))
fname<-paste0(outputdata,'180404-180904-InducibleGenotypesBasedOnGluExpression.rdata')
rm(InducibleFrame)
load(fname)

############
### fit GAL pathway expression in galactose as a logistic function of GAL pathway expression in glucose
############

# annotate a smaller version of DT0
dSub1<-DT0[,.(aa,cc,bb,dd,yfp.mean.glu, yfp.mean.gal)]
dum<-data.table(dd=c('GALK.Esc_col','GALK.Can_alb','GALK.Sac_cer','HIS5.Sch_pom'),GALK.cat=c('GALK (E. coli or C. albicans)\nsensor(-) kinase(+)','GALK (E. coli or C. albicans)\nsensor(-) kinase(+)','GAL1\nsensor(+) kinase(+)','HIS5\nsensor(-) kinase(-)'))
dum[,GALK.cat:=factor(GALK.cat,levels=unique(dum$GALK.cat))]
dSub<-merge(dSub1,dum,by='dd')
dSub[,x:= yfp.mean.glu]
dSub[,y:= yfp.mean.gal]
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


# merge summary and raw data
DTm<-merge(summ,dSub,by=genos)
DTm[,clone.named:=applyPaste(data.frame(aa,bb,cc,dd),' ')]
length(unique(DTm$clone.named))
DTm[,pred:=predict(mod,DTm[,.(x)])]
DTm[,estimate:=y-pred]
# ggplot(DTm,aes(pred,y,col=GALK.cat))+geom_point()+facet_wrap(~bb)+geom_abline()
tframe<-DTm[,t.test(estimate,alternative='greater'),by=genos]
tframe[,p.adj:=p.adjust(p.value)]
tframe[,`FDR < 0.05`:=p.adj<0.05]
tframe[,sig:=`FDR < 0.05`]
tframe[,dup:=applyPaste(tframe[,.(aa,bb,cc,dd)],'')]
DTm2<-merge(DTm,tframe[!duplicated(dup),c(genos,'sig','FDR < 0.05'),with=F],by=genos)

# similar plot but of GAL pathway expression in glucose vs GAL pathway expression in galactose
pKREZ<-ggplot(DTm2,aes(yfp.mean.glu_mean,yfp.mean.gal_mean,xmax= yfp.mean.glu_upper,xmin=yfp.mean.glu_lower,ymax= yfp.mean.gal_upper,ymin=yfp.mean.gal_lower))+
	geom_line(data=lineframe,aes(x,y),inherit.aes=F)+
	geom_errorbar(col='grey50')+geom_errorbarh(col='grey50')+geom_point(aes(col=`FDR < 0.05`),size=2,shape=21,stroke=1)+
	scale_colour_manual(values=c('black','red'))+theme_minimal()+ylab('mean expression of GAL pathway\nin galactose, A.U.')+xlab('mean expression of GAL pathway in glucose, A.U.')+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+
	geom_abline(col='grey50')+facet_grid(GALK.cat~bb)+
	theme(strip.text.y = element_text(angle =0),legend.position='bottom')

if(writeplot==T){
	w<-5.3;h<-6
	ggsave(paste0(figout,date,'-180828-pKREZ-GLUExpressionVsYFPGal-DoseResponse.png'), pKREZ,width=w,height=h)	
}


############
### fit growth rate as a weibull function of fraction ON in glucose
############

# annotate a smaller version of DT0
dSub1<-DT0[,.(aa,cc,bb,dd,fracon.glu, growth.rate)]
dum<-data.table(dd=c('GALK.Esc_col','GALK.Can_alb','GALK.Sac_cer','HIS5.Sch_pom'),GALK.cat=c('GALK (E. coli or C. albicans)\nsensor(-) kinase(+)','GALK (E. coli or C. albicans)\nsensor(-) kinase(+)','GAL1\nsensor(+) kinase(+)','HIS5\nsensor(-) kinase(-)'))
dum[,GALK.cat:=factor(GALK.cat,levels=unique(dum$GALK.cat))]
dSub<-merge(dSub1,dum,by='dd')
dSub[,x:= fracon.glu]
dSub[,y:= growth.rate]
dSub[GALK.cat=='GALK (E. coli or C. albicans)\nsensor(-) kinase(+)'&bb=='GAL3.delta',fac:='no sensors']

ctrl<-nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/(2^20),
            printEval = FALSE, warnOnly = T)

# fit model and calculate deviation from expectation
cutoff<-0.2
dat<-dSub[x> cutoff&fac=='no sensors',.(x,y)]
mod<-nls(y~SSweibull(x, Asym, Drop, lrc, pwr),data=dat)
xes<-data.table(x=seq(0.0001,1,length=100))
lineframe<-data.table(xes,y=predict(mod,xes))
datcoli<-dSub[x> cutoff&dd=='GALK.Esc_col'&bb=='GAL3.delta',.(x,y)]
modcoli<-nls(y~SSlogis(x,Asym,xmid,scal),data=datcoli)
datalb<-dSub[x> cutoff&dd=='GALK.Can_alb'&bb=='GAL3.delta',.(x,y)]
modalb<-nls(y~SSlogis(x,Asym,xmid,scal), data =datalb)
dSub[x<= cutoff,x:= cutoff]
dSub[dd=='GALK.Sac_cer'|dd=='HIS5.Sch_pom',pred:=predict(mod)]
dSub[dd=='GALK.Esc_col',pred:=predict(modcoli)]
dSub[dd=='GALK.Can_alb',pred:=predict(modalb)]


# merge summary and raw data
DTm<-merge(summ,dSub,by=c(genos))
DTm[,clone.named:=applyPaste(data.frame(aa,bb,cc,dd),' ')]
length(unique(DTm$clone.named))
DTm[,pred:=predict(mod,DTm[,.(x)])]
DTm[,estimate:=y-pred]
# ggplot(DTm,aes(pred,y,col=GALK.cat))+geom_point()+facet_wrap(~bb)+geom_abline()
tframe<-DTm[,t.test(estimate,alternative='greater'),by=genos]
tframe[,p.adj:=p.adjust(p.value)]
tframe[,`FDR < 0.05`:=p.adj<0.05]
tframe[,sig:=`FDR < 0.05`]
tframe[,dup:=applyPaste(tframe[,.(aa,bb,cc,dd)],'')]
DTm2<-merge(DTm,tframe[!duplicated(dup),c(genos,'sig','FDR < 0.05'),with=F],by=genos)

# fracon in gucose vs growth rate with model fit above
pADGB<-ggplot(DTm2,aes(fracon.glu_mean,growth.rate_mean,xmax= fracon.glu_upper,xmin=fracon.glu_lower,ymax= growth.rate_upper,ymin=growth.rate_lower))+
	geom_line(data=lineframe,aes(x,y),inherit.aes=F)+
	geom_errorbar(col='grey50')+geom_errorbarh(col='grey50')+geom_point(aes(col=`FDR < 0.05`),size=2,shape=21,stroke=1)+
	scale_colour_manual(values=c('black','red'))+theme_minimal()+ylab(expression(paste('mean ',mu,' ',hr^-1,' in galactose')))+xlab('fraction of cells ON in glucose')+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+
	facet_grid(GALK.cat~bb)+
	theme(strip.text.y = element_text(angle =0),legend.position='bottom')


# show that these cells with low fracon values have low mean expression
backg_gr<-summary.func(unlist(DT0[aa=='GAL4.delta'|dd=='HIS5.Sch_pom',pheno,with=F]))
summ01<-data.table(apply(t(summ[,apply(data.frame(growth.rate_mean,growth.rate_sd,growth.rate_N),1,function(x){
	t.test2(m1=x[1],s1=x[2],n1=x[3],m2= backg_gr$mean,s2= backg_gr$sd, n2=backg_gr$N,list=T)	
})]),2,as.numeric))
colnames(summ01)<-c("diff.of.means", "std.err", "t", "p.value")
summ01[,padj:=p.adjust(p.value)]
summ01[,fac2:=padj<0.05]
summ02<-cbind(summ[,.(aa,bb,cc,dd)],summ01)
DTm22<-merge(DTm2,summ02[,.(aa,bb,cc,dd,fac2)])
DTm22[,`able to grow above background`:=fac2]
sample(LETTERS,4)

pKCTY<-ggplot(DTm22[order(fac2)],aes(fracon.glu_mean,yfp.mean.glu_mean,xmax= fracon.glu_upper,xmin=fracon.glu_lower,ymax= yfp.mean.glu_upper,ymin=yfp.mean.glu_lower))+
#	geom_line(data=lineframe,aes(x,y),inherit.aes=F)+
	geom_errorbar(col='grey50')+geom_errorbarh(col='grey50')+geom_point(aes(col=`able to grow above background`),size=2,shape=21,stroke=1,alpha=0.5)+
	scale_colour_manual(values=c('black','red'))+theme_minimal()+ylab('mean expression of GAL pathway in glucose, A.U.')+xlab('fraction of cells ON in glucose')+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+
	facet_grid(GALK.cat~bb)+
	theme(strip.text.y = element_text(angle =0),legend.position='bottom')

if(writeplot==T){
	w<-5.3;h<-6
	ggsave(paste0(figout,date,'-180828-pADGB-FraconGLUVsGrowthRate-DoseResponse.png'), pADGB,width=w,height=h)	
	ggsave(paste0(figout,date,'-180828-pKCTY-FraconGLUVsYFPMeanGlu.png'), pKCTY,width=w,height=h)	
}


###################################################
### pairwise fitness landscapes for GAL3, GAL1 across different GAL4, GAL80 backgrounds
###################################################
writeplot<-F
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
pDSEX<-FitnessLandscapePlot(DTx,arrowframe2,nogrowth,wtgrowth)

if(writeplot==T){
	w<-3.7;h<-4.2
	ggsave(paste0(figout,date,'-181109-pDSEX-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-ForPubFinalFigure-growth.rate.png'), pDSEX,width=w,height=h)
	ggsave(paste0(figout,date,'-181109-pDSEX-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-ForPubFinalFigure-growth.rate.pdf'), pDSEX,width=w,height=h)
}





pCTPO<-ggplot(DTx,aes(x,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+theme_minimal()+
	geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
	geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
	geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	geom_errorbar(col='grey50')+
	colScale+
	xlab('')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+theme(panel.grid.minor.x = element_blank())+
	scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+
	facet_grid(cc~aa)

if(writeplot==T){
	w<-12;h<-9.5
	ggsave(paste0(figout,date,'-181105-pCTPO-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-all-growth.rate.png'), pCTPO,width=w,height=h)	
	w<-12;h<-9.5
	ggsave(paste0(figout,date,'-181105-pCTPOopt-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-all-growth.rate.png'), pCTPO+textopt,width=w,height=h)	

}


G4<-c('GAL4.WT','GAL4-L868P','GAL4-L868K','GAL4-L868G') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80S-1')
GK<-c('GALK.Sac_cer','GALK.Esc_col')
G3<-c('GAL3.WT','GAL3.delta')
ynudge<-0.02
DTx<-DT[aa%in%G4&cc%in%G80]
arrowframe2<-arrowframe3[aa%in%G4&cc%in%G80]
sample(letters,4)
pBWRG<-ggplot(DTx,aes(x,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+theme_minimal()+
	geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
	geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
	geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	geom_errorbar(col='grey50')+
	colScale+
	xlab('')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+theme(panel.grid.minor.x = element_blank())+
	scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+
	facet_grid(cc~aa)	

if(writeplot==T){
	w<-6.9;h<-4.3
	ggsave(paste0(figout,date,'-181105-pBWRG-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-Subset-growth.rate.png'), pBWRG,width=w,height=h)
}

# even smallser subset
G4<-c('GAL4.WT','GAL4-L868G') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80S-1')
GK<-c('GALK.Sac_cer','GALK.Esc_col')
G3<-c('GAL3.WT','GAL3.delta')
ynudge<-0.02
DTx<-DT[aa%in%G4&cc%in%G80]
arrowframe2<-arrowframe3[aa%in%G4&cc%in%G80]
sample(letters,4)
pTYLE<-ggplot(DTx,aes(x,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+theme_minimal()+
	geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
	geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
	geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	geom_errorbar(col='grey50')+#geom_point(col='grey50')+
	colScale+
	xlab('')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+theme(panel.grid.minor.x = element_blank())+
	scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+
	facet_grid(cc~aa)+theme(legend.position='none')
if(writeplot==T){
	w<-3.9;h<-4.5
	ggsave(paste0(figout,date,'-181105-pTYLE-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-RepresentativeSubset-growth.rate.png'), pTYLE,width=w,height=h)
}



G4<-c('GAL4.WT') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80.07')
GK<-c('GALK.Sac_cer','GALK.Esc_col')
G3<-c('GAL3.WT','GAL3.delta')
ynudge<-0.02
DTx<-DT[aa%in%G4&cc%in%G80]
arrowframe2<-arrowframe3[aa%in%G4&cc%in%G80]
sample(letters,4)
textopt<-geom_text(data=DTx,nudge_y= ynudge,aes(x=x,y=growth.rate_mean,label=shape2,col=(shape2)),size=2)

pADHM<-ggplot(DTx,aes(x,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+theme_minimal()+
	geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
	geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
	geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	geom_errorbar(col='grey50')+#geom_point(col='grey50')+
	colScale+
	xlab('')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+theme(panel.grid.minor.x = element_blank())+
	scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+
	facet_grid(cc~aa)+theme(legend.position='none')

if(writeplot==T){
	w<-2.6;h<-4.5
	ggsave(paste0(figout,date,'-181105-pADHM-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-SubsetGAL80.07-growth.rate.png'), pADHM,width=w,height=h)
}

library(ggplot2)   
library(gtable)    
library(grid)
library(gridExtra) 
gA<-ggplotGrob(pSMGP)
gB<-ggplotGrob(pYRXI)
gC<-ggplotGrob(pADHM)
gD<-ggplotGrob(pTYLE)
# Set the widths
gA$widths <- gB$widths
gC$widths <- gD$widths

# Arrange the charts.
grid.newpage()
grid.arrange(gB, gA, gD, gC,ncol=2)



###################################################
### OLD BUT SOME STILL GOOD
###################################################




###################################################
### OLD: reverse ordered GALS backgrounds pairwise fitness landscapes for GAL3, GAL1 across different GAL4, GAL80 backgrounds
###################################################

############
### growth rate
############

# summarize data and get together
dSub<-DT0[dd!='HIS5.Sch_pom']
pheno<-'growth.rate'
backg<-summary.func(unlist(DT0[aa=='GAL4.delta'|dd=='HIS5.Sch_pom',pheno,with=F]))
wt<-summary.func(unlist(DT0[aa=='GAL4.WT'&bb=='GAL3.WT'&cc=='GAL80.WT'&dd=='GALK.Sac_cer',pheno,with=F]))
nogrowth<-data.table(x=seq(0,4,length=100),upper=backg$upper,lower=backg$lower)
wtgrowth<-data.table(x=seq(0,4,length=100),upper=wt$upper,lower=wt$lower)
midval1<-1.85
midval2<-2.15
summ<-summary.func.all(dSub,c('growth.rate','phenotypic.index','yfp.mean.gal','yfp.mean.glu'),c('aa','bb','cc','dd'))
summ[dd=='GALK.Esc_col'&bb=='GAL3.delta',c('x','shape','shape2'):=list(1,factor(1),'no sensor')]
summ[dd=='GALK.Sac_cer'&bb=='GAL3.delta',c('x','shape','shape2'):=list(midval1,factor(6),'GAL1')]
summ[dd=='GALK.Esc_col'&bb=='GAL3.WT',c('x','shape','shape2'):=list(midval2,factor(2),'GAL3')]
summ[dd=='GALK.Sac_cer'&bb=='GAL3.WT',c('x','shape','shape2'):=list(3,factor(11),'GAL1+GAL3')]
senslev<-c('no sensor','GAL1','GAL3','GAL1+GAL3')
summ[,shape2:=factor(shape2,levels=senslev)]

arrowframe<-data.table(x0=c(1,1,midval1,midval2),xend=c(midval1,midval2,3,3),shape2=c('no sensor','GAL1','GAL3','GAL1+GAL3'),lev=c(0,1,1,2))
summ2<-merge(summ,arrowframe,by='shape2')
setnames(summ2,paste0(pheno,'_mean'),'phenoi')
sub<-summ2[,.(aa,cc,bb,dd,phenoi,shape2,x0,xend,lev)]
sub0<-sub[lev==0&dd=='GALK.Esc_col',.(aa,cc,x0=1,y0= phenoi)]
sub01a<-sub[lev==1&dd=='GALK.Esc_col',.(aa,cc,xend=midval2,shape2,y0= phenoi)]
sub01b<-sub[lev==1&dd=='GALK.Sac_cer',.(aa,cc,xend=midval1,shape2,y0= phenoi)]
sub01<-rbind(sub01a,sub01b)
setnames(sub01,'y0','yend')
sub01c<-merge(sub0,sub01,by=c('aa','cc'))
sub02<-sub[lev==2,.(aa,cc,xend=3,shape2,yend= phenoi)]
sub02a<-copy(sub01)
setnames(sub02a,c('xend','yend'),c('x0','y0'))
sub02b<-merge(sub02,sub02a[,!'shape2'],by=c('aa','cc'))
dum<-sub02b[1]
dum[,notin(colnames(dum),c('aa','cc'))]<-0
dum[,shape2:='no sensor']
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
pointframe[inducible==T,grp:='inducible']
pointframe[inducible==F,grp:='uninducible']
pointframe[shape2=='no sensor',xend:=1]

myColors <- c('cornflowerblue','orange1','darkgreen','grey50','red','black')
names(myColors) <- c('GAL1','GAL3','GAL1+GAL3','no sensor','inducible','uninducible')
legname<-'GAL sensor genotypes'
colScale <- scale_colour_manual(name = "GAL sensor genotypes",values = myColors)

pAHJE<-ggplot(DT,aes(x,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+theme_minimal()+
	geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
	geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
	geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	geom_errorbar(col='grey50')+#geom_point(col='grey50')+
	# geom_point(data=na.exclude(pointframe),aes(x=xend,y=0.05,col=grp),inherit.aes=F)+
	geom_text(data=DT,nudge_y= ynudge,aes(x=x,y=growth.rate_mean,label=shape2,col=(shape2)),size=2)+
	colScale+
	xlab('')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+theme(panel.grid.minor.x = element_blank())+
	scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+
	facet_grid(cc~aa)

if(writeplot==T){
	w<-12;h<-9.5
	ggsave(paste0(figout,date,'-180828-pAHJE-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-all-growth.rate.png'), pAHJE,width=w,height=h)	
}


xlabels<-c('background','single mutants','double mutants')
G4<-c('GAL4.WT','GAL4-L868P','GAL4-L868K','GAL4-L868G') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80S-1')
GK<-c('GALK.Sac_cer','GALK.Esc_col')
G3<-c('GAL3.WT','GAL3.delta')
ynudge<-0.02
DT<-summ[aa%in%G4&cc%in%G80&dd%in%GK&bb%in%G3]
arrowframe2<-arrowframe3[aa%in%G4&cc%in%G80]
sample(letters,4)
pKRGM<-ggplot(DT,aes(x,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+theme_minimal()+
	geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
	geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
	geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	geom_errorbar(col='grey50')+#geom_point(col='grey50')+
	geom_text(data=DT,nudge_y= ynudge,aes(x=x,y=growth.rate_mean,label=shape2,col=(shape2)),size=1.5)+
#	scale_colour_manual(values=c('grey50','orange1','cornflowerblue','darkgreen'),name='GAL sensor genotypes')+
	scale_colour_manual(values=c('cornflowerblue','darkgreen','orange1','grey50'),name='GAL sensor genotypes')+
	xlab('')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+theme(panel.grid.minor.x = element_blank())+
	scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+
	facet_grid(cc~aa)	

if(writeplot==T){
	w<-6.9;h<-4.3
	ggsave(paste0(figout,date,'-180828-pKRGM-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-Subset-growth.rate.png'), pKRGM,width=w,height=h)
}

# even smallser subset
xlabels<-c('background','single mutants','double mutants')
G4<-c('GAL4.WT','GAL4-L868G') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80S-1')
GK<-c('GALK.Sac_cer','GALK.Esc_col')
G3<-c('GAL3.WT','GAL3.delta')
ynudge<-0.02
DT<-summ[aa%in%G4&cc%in%G80&dd%in%GK&bb%in%G3]
arrowframe2<-arrowframe3[aa%in%G4&cc%in%G80]
sample(letters,4)
pRMOF<-ggplot(DT,aes(x,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+theme_minimal()+
	geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
	geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
	geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	geom_errorbar(col='grey50')+#geom_point(col='grey50')+
	geom_text(data=DT,nudge_y= ynudge,aes(x=x,y=growth.rate_mean,label=shape2,col=(shape2)),size=1.5)+
#	scale_colour_manual(values=c('grey50','orange1','cornflowerblue','darkgreen'),name='GAL sensor genotypes')+
	scale_colour_manual(values=c('cornflowerblue','darkgreen','orange1','grey50'),name='GAL sensor genotypes')+
	xlab('')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+theme(panel.grid.minor.x = element_blank())+
	scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+
	facet_grid(cc~aa)+theme(legend.position='none')

if(writeplot==T){
	w<-3.9;h<-4.5
	ggsave(paste0(figout,date,'-180828-pRMOF-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-RepresentativeSubset-growth.rate.png'), pRMOF,width=w,height=h)
}




xlabels<-c('background','single mutants','double mutants')
G4<-c('GAL4.WT') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80.07')
GK<-c('GALK.Sac_cer','GALK.Esc_col')
G3<-c('GAL3.WT','GAL3.delta')
ynudge<-0.02
DT<-summ[aa%in%G4&cc%in%G80&dd%in%GK&bb%in%G3]
arrowframe2<-arrowframe3[aa%in%G4&cc%in%G80]
sample(letters,4)
pVMNU<-ggplot(DT,aes(x,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+theme_minimal()+
	geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
	geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
	geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	geom_errorbar(col='grey50')+#geom_point(col='grey50')+
	geom_text(data=DT,nudge_y= ynudge,aes(x=x,y=growth.rate_mean,label=shape2,col=(shape2)),size=1.5)+
#	scale_colour_manual(values=c('grey50','orange1','cornflowerblue','darkgreen'),name='GAL sensor genotypes')+
	scale_colour_manual(values=c('cornflowerblue','darkgreen','orange1','grey50'),name='GAL sensor genotypes')+
	xlab('')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+theme(panel.grid.minor.x = element_blank())+
	scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+
	facet_grid(cc~aa)+theme(legend.position='none')

if(writeplot==T){
	w<-2.6;h<-4.5
	ggsave(paste0(figout,date,'-180828-pVMNU-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-SubsetGAL80.07-growth.rate.png'), pVMNU,width=w,height=h)
}



############
### YFP in glucose
############

dSub<-DT0[dd!='HIS5.Sch_pom']
pheno<-'yfp.mean.glu'
backg<-summary.func(unlist(DT0[aa=='GAL4.delta',pheno,with=F]))
wt<-summary.func(unlist(DT0[aa=='GAL4.WT'&bb=='GAL3.WT'&cc=='GAL80.WT'&dd=='GALK.Sac_cer',pheno,with=F]))
nogrowth<-data.table(x=seq(0,4,length=100),upper=backg$upper,lower=backg$lower)
wtgrowth<-data.table(x=seq(0,4,length=100),upper=wt$upper,lower=wt$lower)
midval1<-1.85
midval2<-2.15
summ<-summary.func.all(dSub,c('growth.rate','phenotypic.index','yfp.mean.gal','yfp.mean.glu'),c('aa','bb','cc','dd'))
summ[dd=='GALK.Esc_col'&bb=='GAL3.delta',c('x','shape','shape2'):=list(1,factor(1),'no sensor')]
summ[dd=='GALK.Sac_cer'&bb=='GAL3.delta',c('x','shape','shape2'):=list(midval1,factor(6),'GAL1')]
summ[dd=='GALK.Esc_col'&bb=='GAL3.WT',c('x','shape','shape2'):=list(midval2,factor(2),'GAL3')]
summ[dd=='GALK.Sac_cer'&bb=='GAL3.WT',c('x','shape','shape2'):=list(3,factor(11),'GAL1+GAL3')]
senslev<-c('no sensor','GAL1','GAL3','GAL1+GAL3')
summ[,shape2:=factor(shape2,levels=senslev)]
phenocols<-grepincols(summ,pheno)
setnames(summ,phenocols,gsub(pheno,'phenoi',phenocols))

arrowframe<-data.table(x0=c(1,1,midval1,midval2),xend=c(midval1,midval2,3,3),shape2=c('no sensor','GAL1','GAL3','GAL1+GAL3'),lev=c(0,1,1,2))
summ2<-merge(summ,arrowframe,by='shape2')
sub<-summ2[,.(aa,cc,bb,dd,phenoi_mean,shape2,x0,xend,lev)]
sub0<-sub[lev==0&dd=='GALK.Esc_col',.(aa,cc,x0=1,y0= phenoi_mean)]
sub01a<-sub[lev==1&dd=='GALK.Esc_col',.(aa,cc,xend=midval2,shape2,y0= phenoi_mean)]
sub01b<-sub[lev==1&dd=='GALK.Sac_cer',.(aa,cc,xend=midval1,shape2,y0= phenoi_mean)]
sub01<-rbind(sub01a,sub01b)
setnames(sub01,'y0','yend')
sub01c<-merge(sub0,sub01,by=c('aa','cc'))
sub02<-sub[lev==2,.(aa,cc,xend=3,shape2,yend= phenoi_mean)]
sub02a<-copy(sub01)
setnames(sub02a,c('xend','yend'),c('x0','y0'))
sub02b<-merge(sub02,sub02a[,!'shape2'],by=c('aa','cc'))
dum<-sub02b[1]
dum[,notin(colnames(dum),c('aa','cc'))]<-0
dum[,shape2:='no sensor']
arrowframe3<-rbind(sub01c,sub02b,dum)
arrowframe3[,shape3:=factor(shape2,levels=senslev)]

facy<-1#0.96
facx<-0
xlabels<-c('background','single mutants','double mutants')
ynudge<-5000
DT<-copy(summ)
arrowframe2<-copy(arrowframe3)
sample(letters,4)
xlabel<-''
ylabel<-'Mean YFP signal in glucose, AU'#expression(paste('growth rate ',mu,' ',hr^-1,sep=' '))
pCFKD<-ggplot(DT,aes(x,phenoi_mean,ymax=phenoi_upper,ymin=phenoi_lower))+theme_minimal()+
	geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='red',col='transparent',alpha=1,inherit.aes=F)+
	geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='red',col='transparent',alpha=0.5,inherit.aes=F)+
	geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	geom_errorbar(col='grey50')+#geom_point(col='grey50')+
	geom_text(data=DT,nudge_y= ynudge,aes(x=x,y= phenoi_mean,label=shape2,col=(shape2)),size=2)+
#	scale_colour_manual(values=c('grey50','orange1','cornflowerblue','darkgreen'),name='GAL sensor genotypes')+
	scale_colour_manual(values=c('cornflowerblue','darkgreen','orange1','grey50'),name='GAL sensor genotypes')+
	xlab(xlabel)+ylab(ylabel)+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+theme(panel.grid.minor.x = element_blank())+
	scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+#	scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+
	facet_grid(cc~aa)	
if(writeplot==T){
	w<-12;h<-9.5
	ggsave(paste0(figout,date,'-180828-pCFKD-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-all-yfp.mean.glu.png'), pCFKD,width=w,height=h)
}

xlabels<-c('background','single mutants','double mutants')
G4<-c('GAL4.WT','GAL4-L868P','GAL4-L868K','GAL4-L868G') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80S-1')
GK<-c('GALK.Sac_cer','GALK.Esc_col')
G3<-c('GAL3.WT','GAL3.delta')
ynudge<-5000
DT<-summ[aa%in%G4&cc%in%G80&dd%in%GK&bb%in%G3]
arrowframe2<-arrowframe3[aa%in%G4&cc%in%G80]
xlabel<-''
ylabel<-'Mean YFP signal in glucose, AU'#expression(paste('growth rate ',mu,' ',hr^-1,sep=' '))

sample(letters,4)
pJTYP<-ggplot(DT,aes(x,phenoi_mean,ymax=phenoi_upper,ymin=phenoi_lower))+theme_minimal()+
	geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='red',col='transparent',alpha=1,inherit.aes=F)+
	geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='red',col='transparent',alpha=0.5,inherit.aes=F)+
	geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	geom_errorbar(col='grey50')+#geom_point(col='grey50')+
	geom_text(data=DT,nudge_y= ynudge,aes(x=x,y= phenoi_mean,label=shape2,col=(shape2)),size=1.5)+
#	scale_colour_manual(values=c('grey50','orange1','cornflowerblue','darkgreen'),name='GAL sensor genotypes')+
	scale_colour_manual(values=c('cornflowerblue','darkgreen','orange1','grey50'),name='GAL sensor genotypes')+
	xlab(xlabel)+ylab(ylabel)+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+theme(panel.grid.minor.x = element_blank())+
	scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+#	scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+
	facet_grid(cc~aa)	

if(writeplot==T){
	w<-6.9;h<-4.3
	ggsave(paste0(figout,date,'-180828-pJTYP-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-Subset-yfp.mean.glu.png'), pJTYP,width=w,height=h)
}
############
### pairwise fitness landscapes for GAL3, GAL1 across different GAL4, GAL80 backgrounds: YFP in galactose
############

dSub<-DT0[dd!='HIS5.Sch_pom']
pheno<-'yfp.mean.gal'
backg<-summary.func(unlist(DT0[aa=='GAL4.delta',pheno,with=F]))
wt<-summary.func(unlist(DT0[aa=='GAL4.WT'&bb=='GAL3.WT'&cc=='GAL80.WT'&dd=='GALK.Sac_cer',pheno,with=F]))
nogrowth<-data.table(x=seq(0,4,length=100),upper=backg$upper,lower=backg$lower)
wtgrowth<-data.table(x=seq(0,4,length=100),upper=wt$upper,lower=wt$lower)
midval1<-1.85
midval2<-2.15
summ<-summary.func.all(dSub,c('growth.rate','phenotypic.index','yfp.mean.gal','yfp.mean.glu'),c('aa','bb','cc','dd'))
summ[dd=='GALK.Esc_col'&bb=='GAL3.delta',c('x','shape','shape2'):=list(1,factor(1),'no sensor')]
summ[dd=='GALK.Sac_cer'&bb=='GAL3.delta',c('x','shape','shape2'):=list(midval1,factor(6),'GAL1')]
summ[dd=='GALK.Esc_col'&bb=='GAL3.WT',c('x','shape','shape2'):=list(midval2,factor(2),'GAL3')]
summ[dd=='GALK.Sac_cer'&bb=='GAL3.WT',c('x','shape','shape2'):=list(3,factor(11),'GAL1+GAL3')]
senslev<-c('no sensor','GAL1','GAL3','GAL1+GAL3')
summ[,shape2:=factor(shape2,levels=senslev)]
phenocols<-grepincols(summ,pheno)
setnames(summ,phenocols,gsub(pheno,'phenoi',phenocols))



arrowframe<-data.table(x0=c(1,1,midval1,midval2),xend=c(midval1,midval2,3,3),shape2=c('no sensor','GAL1','GAL3','GAL1+GAL3'),lev=c(0,1,1,2))
summ2<-merge(summ,arrowframe,by='shape2')
sub<-summ2[,.(aa,cc,bb,dd,phenoi_mean,shape2,x0,xend,lev)]
sub0<-sub[lev==0&dd=='GALK.Esc_col',.(aa,cc,x0=1,y0= phenoi_mean)]
sub01a<-sub[lev==1&dd=='GALK.Esc_col',.(aa,cc,xend=midval2,shape2,y0= phenoi_mean)]
sub01b<-sub[lev==1&dd=='GALK.Sac_cer',.(aa,cc,xend=midval1,shape2,y0= phenoi_mean)]
sub01<-rbind(sub01a,sub01b)
setnames(sub01,'y0','yend')
sub01c<-merge(sub0,sub01,by=c('aa','cc'))
sub02<-sub[lev==2,.(aa,cc,xend=3,shape2,yend= phenoi_mean)]
sub02a<-copy(sub01)
setnames(sub02a,c('xend','yend'),c('x0','y0'))
sub02b<-merge(sub02,sub02a[,!'shape2'],by=c('aa','cc'))
dum<-sub02b[1]
dum[,notin(colnames(dum),c('aa','cc'))]<-0
dum[,shape2:='no sensor']
arrowframe3<-rbind(sub01c,sub02b,dum)
arrowframe3[,shape3:=factor(shape2,levels=senslev)]

# arrowframe4<-rbind(arrowframe3,arrowframe3[nrow(arrowframe3)+1])
# arrowframe4[nrow(arrowframe4),shape2:='no sensor']
# arrowframe4[,shape2:=factor(shape2,levels=senslev)]
# arrowframe2<-na.exclude(arrowframe4)

facy<-1#0.96
facx<-0
xlabels<-c('background','single mutants','double mutants')
ynudge<-5000
DT<-copy(summ)
arrowframe2<-copy(arrowframe3)
sample(letters,4)
xlabel<-''
ylabel<-'Mean YFP signal in galactose, AU'#expression(paste('growth rate ',mu,' ',hr^-1,sep=' '))
pILNQ<-ggplot(DT,aes(x,phenoi_mean,ymax=phenoi_upper,ymin=phenoi_lower))+theme_minimal()+
#	geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='red',col='transparent',alpha=1,inherit.aes=F)+
#	geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='red',col='transparent',alpha=0.5,inherit.aes=F)+
	geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	geom_errorbar(col='grey50')+#geom_point(col='grey50')+
	geom_text(data=DT,nudge_y= ynudge,aes(x=x,y= phenoi_mean,label=shape2,col=(shape2)),size=2)+
#	scale_colour_manual(values=c('grey50','orange1','cornflowerblue','darkgreen'),name='GAL sensor genotypes')+
	scale_colour_manual(values=c('cornflowerblue','darkgreen','orange1','grey50'),name='GAL sensor genotypes')+
	xlab(xlabel)+ylab(ylabel)+
	theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+theme(panel.grid.minor.x = element_blank())+
	scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+#	scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+
	facet_grid(cc~aa)	

if(writeplot==T){
	w<-12;h<-9.5
	ggsave(paste0(figout,date,'-180828-pILNQ-PairwiseFitnessGAL3GAL1AcrossGAL4GAL80-Subset-yfp.mean.gal.png'), pILNQ,width=w,height=h)
}













###################################################
### make expression plot figures of different backgrounds based on a GALK background
###################################################
GALKVars<-c('HIS5.Sch_pom','GALK.Can_alb','GALK.Esc_col','GALK.Sac_cer' )
for(i in 1:length(GALKVars)){	
	GALKvar<-GALKVars[i]
	DTsub<-DT0[GALK.clone.named== GALKvar]
	
	GAL3wt<-DTsub[GAL3.clone.named=='GAL3.WT']
	GAL3delta<-DTsub[GAL3.clone.named=='GAL3.delta']
	

	dens.cols<-colnames(DT0)[grep('dens',colnames(DT0))]
	WTdens<-GAL3wt[order(GAL4.clone.named,GAL80.clone.named),lapply(.SD,mean,na.rm=T),by=c('GAL80.clone.named','GAL4.clone.named'),.SDcols=dens.cols]
	DELTAdens<-GAL3delta[order(GAL4.clone.named,GAL80.clone.named),lapply(.SD,mean,na.rm=T),by=c('GAL80.clone.named','GAL4.clone.named'),.SDcols=dens.cols]
	
	labs<-melt(WTdens[,c('GAL4.clone.named','GAL80.clone.named'),with=F])
	table(labs)
	labs$variable<-paste0('V',rownames(WTdens))
	
	
	# GAL3.WT DENSITY PLOTS
	
	df<-as.data.frame(t(WTdens[,3:ncol(WTdens)]))
	df$sugar<-sapply(rownames(df),function(x){
		strsplit(x,'\\.')[[1]][4]
	})
	mDT<-mdf<-data.frame(melt(df,id='sugar'))
	mdf<-merge(mDT,labs,by='variable')
	mdf$variable<-as.numeric(gsub('V','',mdf$variable))
	mdf$pseudoLog10.Fluorescence.AU<-br[-1]
	
	
	mdf$variable2<-paste0(mdf$GAL4.clone.named,'x',mdf$GAL80.clone.named)
	mdf$variable2<-factor(mdf$variable2,levels=unique(data.table(mdf)[order(variable)]$variable2))
	p5<-ggplot(mdf,aes(x=pseudoLog10.Fluorescence.AU,y=value))+geom_line(aes(colour=sugar))+facet_wrap(~variable2,ncol= length(unique(mdf$GAL80.clone.named)))+
		scale_colour_manual(values=c('cornflowerblue','indianred'))+coord_cartesian(ylim=c(0,0.2))+
		theme( strip.background = element_blank(), strip.text.x = element_text(size=5),
			 panel.spacing=unit(0,'lines'),axis.text.y=element_blank(),
			 axis.ticks=element_blank(), axis.title.x=element_blank())
	
	w=10;h=w*0.9
	ggsave(paste0(figout,date,'-',GALKvar,'-GAL3.WT.background.pdf'),grid.arrange(p5,top=paste0('Background: GAL3.WT in ', GALKvar)),width=w,height=h)
	
	
	# GAL3.delta DENSITY PLOTS
	
	df<-as.data.frame(t(DELTAdens[,3:ncol(DELTAdens)]))
	df$sugar<-sapply(rownames(df),function(x){
		strsplit(x,'\\.')[[1]][4]
	})
	mDT<-mdf<-data.frame(melt(df,id='sugar'))
	mdf<-merge(mDT,labs,by='variable')
	mdf$variable<-as.numeric(gsub('V','',mdf$variable))
	mdf$pseudoLog10.Fluorescence.AU<-br[-1]
	
	
	mdf$variable2<-paste0(mdf$GAL4.clone.named,'x',mdf$GAL80.clone.named)
	mdf$variable2<-factor(mdf$variable2,levels=unique(data.table(mdf)[order(variable)]$variable2))
	p6<-ggplot(mdf,aes(x=pseudoLog10.Fluorescence.AU,y=value))+geom_line(aes(colour=sugar))+facet_wrap(~variable2,ncol= length(unique(mdf$GAL80.clone.named)))+
		scale_colour_manual(values=c('cornflowerblue','indianred'))+coord_cartesian(ylim=c(0,0.2))+
		theme( strip.background = element_blank(), strip.text.x = element_text(size=5),
			 panel.spacing=unit(0,'lines'),axis.text.y=element_blank(),
			 axis.ticks=element_blank(), axis.title.x=element_blank())
	#dev.new();p6
	w=10;h=w*0.9		 
	ggsave(paste0(figout,date,'-',GALKvar,'-GAL3.delta.background.pdf'),grid.arrange(p6,top=paste0('Background: GAL3.delta in', GALKvar)),width=w,height=h)
	
	
	# put GAL3.WT and GAL3.delta in same graphs
	
	
	df<-as.data.frame(t(DELTAdens[,3:ncol(DELTAdens)]))
	df$sugar<-sapply(rownames(df),function(x){
		strsplit(x,'\\.')[[1]][4]
	})
	mDT<-mdf<-data.frame(melt(df,id='sugar'))
	mdf<-merge(mDT,labs,by='variable')
	mdf$variable<-as.numeric(gsub('V','',mdf$variable))
	mdf$pseudoLog10.Fluorescence.AU<-br[-1]
	mdf1<-mdf
	df<-as.data.frame(t(WTdens[,3:ncol(WTdens)]))
	df$sugar<-sapply(rownames(df),function(x){
		strsplit(x,'\\.')[[1]][4]
	})
	mDT<-mdf<-data.frame(melt(df,id='sugar'))
	mdf<-merge(mDT,labs,by='variable')
	mdf$variable<-as.numeric(gsub('V','',mdf$variable))
	mdf$pseudoLog10.Fluorescence.AU<-br[-1]
	mdf2<-mdf
	mdf1$GAL3.background<-"GAL3.delta"
	mdf2$GAL3.background<-"GAL3.WT"
	
	mdf<-rbind(mdf1,mdf2)
	
	mdf$variable2<-paste0(mdf$GAL4.clone.named,'x',mdf$GAL80.clone.named)
	mdf$variable2<-factor(mdf$variable2,levels=unique(data.table(mdf)[order(variable)]$variable2))
	p7<-ggplot(mdf,aes(x=pseudoLog10.Fluorescence.AU,y=value))+geom_line(aes(colour=sugar,linetype= GAL3.background))+facet_wrap(~variable2,ncol= length(unique(mdf$GAL80.clone.named)))+
		scale_colour_manual(values=c('cornflowerblue','indianred'))+coord_cartesian(ylim=c(0,0.2))+
		theme( strip.background = element_blank(), strip.text.x = element_text(size=5),
			 panel.spacing=unit(0,'lines'),axis.text.y=element_blank(),
			 axis.ticks=element_blank(), axis.title.x=element_blank())
	
	w=10;h=w*0.9		 
	ggsave(paste0(figout,date,'-',GALKvar,'-GAL3.WTandDelta.pdf'),grid.arrange(p7,top=paste0('GALK background: ',GALKvar)),width=w,height=h)
}

###################################################
### quick look at GAL4-L868P
###################################################

DT0[GAL4.clone.named=='GAL4-L868P',summary.func(growth.rate),by=c('GALK.clone.named','GAL3.clone.named','GAL80.clone.named')][order(GALK.clone.named,GAL3.clone.named,GAL80.clone.named)]
DT0[GAL4.clone.named=='GAL4-L868P',summary.func(yfp.mean.gal/yfp.mean.glu),by=c('GALK.clone.named','GAL3.clone.named','GAL80.clone.named')][order(GALK.clone.named,GAL3.clone.named,GAL80.clone.named)]

DT0[GAL4.clone.named=='GAL4-L868P',summary.func((yfp.mean.gal*2^(generations))/yfp.mean.glu),by=c('GALK.clone.named','GAL3.clone.named','GAL80.clone.named')][order(GALK.clone.named,GAL3.clone.named,GAL80.clone.named)]


###################################################
### quick look at GAL80.07 vs GAL80.WT
###################################################


vars1<-c('GALK.clone.named','GAL3.clone.named','GAL4.clone.named','GAL80.clone.named')
x<-'growth rate 1/hr'
GAL8007<-DT0[GAL80.clone.named=='GAL80.07',summary.func(growth.rate),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
GAL80WT<-DT0[GAL80.clone.named=='GAL80.WT',summary.func(growth.rate),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
G80m<-merge(GAL80WT,GAL8007,by=vars1[1:3])
ggplot(G80m,aes(x=mean.x,y=mean.y))+geom_errorbarh(xmax= G80m $upper.x,xmin= G80m $lower.x)+geom_errorbar(ymax= G80m $upper.y,ymin= G80m $lower.y)+
	geom_point(aes(shape=GAL3.clone.named),size=3)+
	geom_abline(intercept=0,slope=1,col='red') + xlab(paste0('GAL80.WT ',x)) + ylab(paste0('GAL80.07 ',x))
G80m[mean.y>0.11&mean.x<0.1]

x<-'mean expresssion in glucose'
GAL8007<-DT0[GAL80.clone.named=='GAL80.07',summary.func(yfp.mean.glu),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
GAL80WT<-DT0[GAL80.clone.named=='GAL80.WT',summary.func(yfp.mean.glu),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
G80m<-merge(GAL80WT,GAL8007,by=vars1[1:3])
ggplot(G80m,aes(x=mean.x,y=mean.y))+geom_errorbarh(xmax= G80m $upper.x,xmin= G80m $lower.x)+geom_errorbar(ymax= G80m $upper.y,ymin= G80m $lower.y)+
	geom_point(aes(shape=GAL3.clone.named),size=3)+
	geom_abline(intercept=0,slope=1,col='red') + xlab(paste0('GAL80.WT ',x)) + ylab(paste0('GAL80.07 ',x))

x<-'fraction ON in glucose'
GAL8007<-DT0[GAL80.clone.named=='GAL80.07',summary.func(fracon.glu),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
GAL80WT<-DT0[GAL80.clone.named=='GAL80.WT',summary.func(fracon.glu),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
G80m<-merge(GAL80WT,GAL8007,by=vars1[1:3])
ggplot(G80m,aes(x=mean.x,y=mean.y))+geom_errorbarh(xmax= G80m $upper.x,xmin= G80m $lower.x)+geom_errorbar(ymax= G80m $upper.y,ymin= G80m $lower.y)+
	geom_point(aes(shape=GAL3.clone.named),size=3)+
	geom_abline(intercept=0,slope=1,col='red') + xlab(paste0('GAL80.WT ',x)) + ylab(paste0('GAL80.07 ',x))

x<-'growth rate 1/hr'
GAL8007<-DT0[GAL80.clone.named=='GAL80.07',summary.func(growth.rate),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
GAL80WT<-DT0[GAL80.clone.named=='GAL80.WT',summary.func(growth.rate),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
GAL8007a<-DT0[GAL80.clone.named=='GAL80.07',summary.func(fracon.glu,'a'),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
#GAL80WTa<-DT0[GAL80.clone.named=='GAL80.WT',summary.func(fracon.glu),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
G80m<-merge(merge(GAL80WT,GAL8007,by=vars1[1:3]), GAL8007a,by=vars1[1:3])

collim<-signif(as.vector(c(quantile(G80m $a_mean,0.005,na.rm=T),quantile(G80m $a_mean,0.995,na.rm=T))),3)
G80m[a_mean<=collim[1]]$a_mean<-collim[1]
G80m[a_mean>=collim[2]]$a_mean<-collim[2]
brks<-signif(seq(collim[1],collim[2],length=5),3)
labs<-round(brks,1)
	
ggplot(G80m,aes(x=mean.x,y=mean.y))+geom_errorbarh(xmax= G80m $upper.x,xmin= G80m $lower.x)+geom_errorbar(ymax= G80m $upper.y,ymin= G80m $lower.y)+
	geom_point(aes(shape=GAL3.clone.named,col= a_mean),size=3)+
	scale_colour_gradientn(colours=c('cornflowerblue','yellow','indianred'),breaks=brks, labels=labs,limits = collim)+
	geom_abline(intercept=0,slope=1,col='red') + xlab(paste0('GAL80.WT ',x)) + ylab(paste0('GAL80.07 ',x))

x<-'growth rate 1/hr'
GAL8007<-DT0[GAL80.clone.named=='GAL80.07',summary.func(growth.rate),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
GAL80WT<-DT0[GAL80.clone.named=='GAL80.WT',summary.func(growth.rate),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
GAL8007a<-DT0[GAL80.clone.named=='GAL80.07',summary.func(yfp.glu,'a'),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
#GAL80WTa<-DT0[GAL80.clone.named=='GAL80.WT',summary.func(fracon.glu),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
G80m<-merge(merge(GAL80WT,GAL8007,by=vars1[1:3]), GAL8007a,by=vars1[1:3])

collim<-signif(as.vector(c(quantile(G80m $a_mean,0.005,na.rm=T),quantile(G80m $a_mean,0.995,na.rm=T))),3)
G80m[a_mean<=collim[1]]$a_mean<-collim[1]
G80m[a_mean>=collim[2]]$a_mean<-collim[2]
brks<-signif(seq(collim[1],collim[2],length=5),3)
labs<-round(brks,1)
	
ggplot(G80m,aes(x=mean.x,y=mean.y, ymax=upper.y,ymin=lower.y, xmax= upper.x,xmin=lower.x))+geom_errorbarh()+geom_errorbar()+
	geom_point(col='black',aes(shape=GALK.clone.named),size=3.5)+
	geom_point(aes(shape=GALK.clone.named,col= a_mean),size=2)+
	scale_colour_gradientn(colours=c('cornflowerblue','yellow','indianred'),breaks=brks, labels=labs,limits = collim)+
	geom_abline(intercept=0,slope=1,col='red') + xlab(paste0('GAL80.WT ',x)) + ylab(paste0('GAL80.07 ',x))+facet_wrap(~GAL3.clone.named)

x<-'growth rate 1/hr'

###################################################
### quick look at GAL80s vs GAL80.WT
###################################################

collim<-signif(as.vector(c(quantile(DT0 $yfp.glu,0.005,na.rm=T),quantile(DT0 $yfp.glu,0.995,na.rm=T))),3)
brks<-signif(seq(collim[1],collim[2],length=3),3)
labs<-round(brks,1)
GAL80s<-levels(DT0$GAL80.clone.named)
plot.list<-lapply(1:length(GAL80s),function(i){
	
	GAL8007<-DT0[GAL80.clone.named==GAL80s[i],summary.func(growth.rate),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
	GAL80WT<-DT0[GAL80.clone.named=='GAL80.WT',summary.func(growth.rate),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
	GAL8007a<-DT0[GAL80.clone.named==GAL80s[i],summary.func(yfp.glu,'a'),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
	#GAL80WTa<-DT0[GAL80.clone.named=='GAL80.WT',summary.func(fracon.glu),by=vars1][order(GALK.clone.named,GAL3.clone.named,GAL4.clone.named)]
	G80m<-merge(merge(GAL80WT,GAL8007,by=vars1[1:3]), GAL8007a,by=vars1[1:3])
	G80m[a_mean<=collim[1]]$a_mean<-collim[1]
	G80m[a_mean>=collim[2]]$a_mean<-collim[2]
	setnames(G80m,'a_mean','intial expression level in glucose')
		
	ggplot(G80m,aes(x=mean.x,y=mean.y, ymax=upper.y,ymin=lower.y, xmax= upper.x,xmin=lower.x))+geom_errorbarh()+geom_errorbar()+
		geom_point(col='black',aes(shape=GALK.clone.named),size=3.5)+
		geom_point(aes(shape=GALK.clone.named,col= `intial expression level in glucose`),size=2)+
		scale_colour_gradientn(colours=c('cornflowerblue','yellow','indianred'),breaks=brks, labels=labs,limits = collim)+
		geom_abline(intercept=0,slope=1,col='red') + xlab(paste0('GAL80.WT ',x)) + ylab(paste0(GAL80s[i],' ',x))+facet_wrap(~GAL3.clone.named)
	})

h<-16;w<-h*0.9
ggsave(paste0(figout,date,'-GrowthRateInDifferentGAL80,GAL3backgrounds.pdf'),
	grid_arrange_shared_legend(plot.list[[1]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],ncol=2,nrow=3, for.new.dev=F),width=w,height=h)
grid_arrange_shared_legend(plot.list[[1]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],ncol=2,nrow=3)
###################################################
### linear modeling
###################################################
DT1<-DT0[!is.na(growth.rate)]
mod1<-lm(growth.rate~GAL80.clone.named+GAL4.clone.named+GALK.clone.named+GAL3.clone.named,data=DT1)
summary(mod1)
DT1$predicted_linear.model_growth.rate<-mod1$fitted.values
ggplot(DT1,aes(x= predicted_linear.model_growth.rate,y=growth.rate))+geom_point()
aovmod1<-aov(mod1)
TukeyHSD(aovmod1)

mod2<-lm(growth.rate~GAL80.clone.named*GAL4.clone.named*GALK.clone.named*GAL3.clone.named,data=DT1)
summary(mod2)
DT1$predicted_linear.model.interaction_growth.rate<-mod2$fitted.values
ggplot(DT1,aes(x= predicted_linear.model.interaction_growth.rate,y=growth.rate))+geom_point()
aovmod2<-aov(mod2)
summary(aovmod2)
TukeyHSD(aovmod2)

mod3<-lm(yfp.gal~GAL80.clone.named*GAL4.clone.named*GALK.clone.named*GAL3.clone.named,data=DT1)
summary(mod3)
DT1$predicted_linear.model.interaction_yfp.gal<-mod3$fitted.values
ggplot(DT1,aes(x= predicted_linear.model.interaction_yfp.gal,y=yfp.gal))+geom_point()
aovmod3<-aov(mod3)
summary(aovmod3)



mod4<-lm(yfp.gal~GAL80.clone.named*GAL4.clone.named*GAL3.clone.named,data=DT1)
summary(mod4)
DT1$predicted_linear.model.interaction_noGALK_yfp.gal<-mod4 $fitted.values
ggplot(DT1,aes(x= predicted_linear.model.interaction_noGALK_yfp.gal,y=yfp.gal))+geom_point()
aovmod4<-aov(mod4)
summary(aovmod3)

DT1[,surprise_GALK:=yfp.gal-predicted_linear.model.interaction_noGALK_yfp.gal]
DT1$surprise_GALK

DT<-DT1

collim<-signif(as.vector(c(quantile(DT$surprise_GALK,0.005,na.rm=T),quantile(DT$surprise_GALK,0.995,na.rm=T))),3)
brks<-signif(seq(collim[1],collim[2],length=5),3)
labs<-round(brks,1)
DT[surprise_GALK <=collim[1],'surprise_GALK'] <-collim[1]
DT[surprise_GALK>=collim[2],'surprise_GALK'] <-collim[2]

x<-'epistasis'
ggplot(DT,aes(x=GAL80.clone.named,y=GAL4.clone.named))+ geom_tile(aes(fill=surprise_GALK))+
	scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),guide = guide_legend(title = NULL),breaks=brks, labels=labs,limits = collim) +
	facet_wrap(~GALK.clone.named)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))+ggtitle(paste0('Phenotype: ',x,' in GAL3 WT backround'))


mod5<-lm(growth.rate~(GAL80.clone.named*GAL4.clone.named)+GAL3.clone.named+GALK.clone.named,data=DT1[GALK.clone.named!='HIS5.Sch_pom'])
summary(mod5)
DT1$predicted_linear.model.interaction_noGALK_growth.rate<-mod5 $fitted.values
ggplot(DT1,aes(x= predicted_linear.model.interaction_noGALK_growth.rate,y=yfp.gal))+geom_point()
# aovmod5<-aov(mod5)
# summary(aovmod5)

DT1[,surprise.GALK_growth.rate:= growth.rate -predicted_linear.model.interaction_noGALK_growth.rate]

DT<-DT1[GALK.clone.named!='HIS5.Sch_pom']
DT[,pval:=t.test(surprise.GALK_growth.rate)$p.value,by=vars1]
DT[,pval.bonf:=p.adjust(pval),by=vars1]
collim<-signif(as.vector(c(quantile(DT$surprise.GALK_growth.rate,0.005,na.rm=T),quantile(DT$surprise.GALK_growth.rate,0.995,na.rm=T))),3)
brks<-signif(seq(collim[1],collim[2],length=4),3)
labs<-round(brks,2)

# custom breaks and labs
brks<-labs<-c(collim[1],-0.06,-0.03,0,0.03,0.06,collim[2])

DT[surprise.GALK_growth.rate <=collim[1],'surprise.GALK_growth.rate'] <-collim[1]
DT[surprise.GALK_growth.rate>=collim[2],'surprise.GALK_growth.rate'] <-collim[2]
DT[pval.bonf<0.1&growth.rate>0.07,epistasis:=mean(surprise.GALK_growth.rate),by=vars1]
x<-'epistasis'
ggplot(DT[,unique(epistasis),by=vars1],aes(x=GAL80.clone.named,y=GAL4.clone.named))+ geom_tile(aes(fill= V1))+geom_text(aes(label = round(V1, 2)),col='grey50')+
	scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),na.value='white',guide = guide_legend(title = NULL),breaks=brks, labels=labs,limits = collim) +
	facet_wrap(~GAL3.clone.named+GALK.clone.named)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))+ggtitle(paste0('Epistasis when varoius things are ignored in predictive interaction model'))

DT[GALK.clone.named=='GALK.Esc_col',mean(growth.rate),by=vars1][order(GAL3.clone.named,GAL80.clone.named,GAL4.clone.named)]
DT0[GAL80.clone.named=='GAL80.WT'&GAL4.clone.named=='GAL4-L868G',summary.func(growth.rate),by=vars1]


###################################################
### calculate "fold induction"
###################################################
vars1<-c('GALK.clone.named','GAL3.clone.named','GAL4.clone.named','GAL80.clone.named')
DT0[,fold.induction:=yfp.mean.gal/yfp.mean.glu]
qplot(DT0$yfp.mean.glu,DT0$yfp.mean.gal,col=DT0$GALK.clone.named)



###################################################
### make tile plots
###################################################
i<-1
for(i in 1:length(dependents2)){
	x<-dependents2[i]
	DTx<-DT0[,c(vars1,x),with=F]
	colnames(DTx)<-c(vars1,'V1')
	summ<-DTx[,summary.func(V1),by=c(vars1[1:4])]
	collim<-signif(as.vector(c(quantile(summ $mean,0.005,na.rm=T),quantile(summ $mean,0.995,na.rm=T))),3)
	summ[mean<=collim[1]]$mean<-collim[1]
	summ[mean>=collim[2]]$mean<-collim[2]
	brks<-signif(seq(collim[1],collim[2],length=5),3)
	labs<-round(brks,3)
	
	p1<-ggplot(summ[GAL3.clone.named=='GAL3.WT'],aes(x=GAL4.clone.named,y=GAL80.clone.named))+geom_tile(aes(fill=mean))+
		scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),guide = guide_legend(title = NULL),breaks=brks, labels=labs,limits = collim) +
		facet_wrap(~GALK.clone.named)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))+ggtitle(paste0('Phenotype: ',x,' in GAL3 WT backround'))
	p2<-ggplot(summ[GAL3.clone.named=='GAL3.delta'],aes(x=GAL4.clone.named,y=GAL80.clone.named))+geom_tile(aes(fill=mean))+
		scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),guide = guide_legend(title = NULL),breaks=brks, labels=labs,limits = collim) +
		facet_wrap(~GALK.clone.named)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))+ggtitle(paste0('Phenotype: ',x,' in GAL3 Deletion backround'))
	
	w<-12*1.5*1.333;h=w/3.33
	
	summerge<-merge(summ[GAL3.clone.named=='GAL3.WT'],summ[GAL3.clone.named=='GAL3.delta'],by=vars1[c(1,3:4)])
	p3<-ggplot(summerge,aes(mean.x, mean.y))+geom_point()+
		geom_errorbarh(xmax= summerge $upper.x,xmin= summerge $lower.x)+geom_errorbar(ymax= summerge $upper.y,ymin= summerge $lower.y)+
		geom_abline(intercept=0,slope=1,col='red') + xlab(paste0('GAL3 WT ',x)) + ylab(paste0('GAL3 delta ',x))+
		facet_wrap(~GALK.clone.named)
	
	ggsave(paste0(figout,date,'-TilePlot-',x,'.pdf'),do.call(grid.arrange,list(p1,p2,p3,ncol=3)),width=w,height=h)
}
	
ggsave(paste0(figout,date,'-FoldInductionScatterplot-GAL3WTvsGAL3delta.pdf'),p1,width=5,height=5)

###################################################
### take means Scer vs Eco
###################################################
GALKVars<-c('HIS5.Sch_pom','GALK.Can_alb','GALK.Esc_col')#,'GALK.Sac_cer' )
i<-2
GALKvar<-GALKVars[i]
vars1<-c('GALK.clone.named','GAL3.clone.named','GAL4.clone.named','GAL80.clone.named')
Scer<-DT0[GALK.clone.named=='GALK.Sac_cer']
Eco<-DT0[GALK.clone.named==GALKvar]

summ<-merge(Scer[,summary.func(yfp.mean.gal,'Scer'),by=vars1], Eco[,summary.func(yfp.mean.gal,'Eco'),by=vars1],by=vars1[2:4])
p1<-ggplot(summ,aes(Scer_mean,Eco_mean))+geom_point()+geom_errorbarh(xmax= summ$Scer_upper,xmin=summ$Scer_lower)+geom_errorbar(ymax=summ$Eco_upper,ymin=summ$Eco_lower)+geom_abline(intercept=0,slope=1,col='red')

summ<-merge(Scer[,summary.func(fracon.gal,'Scer'),by=vars1], Eco[,summary.func(fracon.gal, GALKvar),by=vars1],by=vars1[2:4])
p2<-ggplot(summ,aes(Scer_mean,Eco_mean))+geom_point()+geom_errorbarh(xmax= summ$Scer_upper,xmin=summ$Scer_lower)+geom_errorbar(ymax=summ$Eco_upper,ymin=summ$Eco_lower)+geom_abline(intercept=0,slope=1,col='red')

summ<-merge(Scer[,summary.func(fold.induction,'Scer'),by=vars1], Eco[,summary.func(fold.induction,'Eco'),by=vars1],by=vars1[2:4])
p3<-ggplot(summ,aes(Scer_mean,Eco_mean))+geom_point()+geom_errorbarh(xmax= summ$Scer_upper,xmin=summ$Scer_lower)+geom_errorbar(ymax=summ$Eco_upper,ymin=summ$Eco_lower)+geom_abline(intercept=0,slope=1,col='red')

summ<-merge(Scer[,summary.func(growth.rate,'Scer'),by=vars1], Eco[,summary.func(growth.rate,'Eco'),by=vars1],by=vars1[2:4])
p4<-ggplot(summ,aes(Scer_mean,Eco_mean))+geom_point()+geom_errorbarh(xmax= summ$Scer_upper,xmin=summ$Scer_lower)+geom_errorbar(ymax=summ$Eco_upper,ymin=summ$Eco_lower)+geom_abline(intercept=0,slope=1,col='red')

w<-12;h=w/2.8
ggsave(paste0(figout,date,'-KeyStatComparisonsGALK-ScervsEcol.pdf'),do.call(grid.arrange,list(p1,p2,p3,p4,top='mean expression galactose, fraction ON in galactose, the fold induction, and growth rate',ncol=4)),height=h,width=w)

###################################################
### take means Scer vs HIS3
###################################################

Scer<-DT0[GALK.clone.named=='Scer.GALK']
HIS<-DT0[GALK.clone.named=='Spomb.HIS5']

summ<-merge(Scer[,summary.func(yfp.mean.gal,'Scer'),by=vars1], HIS[,summary.func(yfp.mean.gal,'HIS'),by=vars1],by=vars1[2:4])
p1<-ggplot(summ,aes(Scer_mean,HIS_mean))+geom_point()+geom_errorbarh(xmax= summ$Scer_upper,xmin=summ$Scer_lower)+geom_errorbar(ymax=summ$HIS_upper,ymin=summ$HIS_lower)+geom_abline(intercept=0,slope=1,col='red')

summ<-merge(Scer[,summary.func(fracon.gal,'Scer'),by=vars1], HIS[,summary.func(fracon.gal,'HIS'),by=vars1],by=vars1[2:4])
p2<-ggplot(summ,aes(Scer_mean,HIS_mean))+geom_point()+geom_errorbarh(xmax= summ$Scer_upper,xmin=summ$Scer_lower)+geom_errorbar(ymax=summ$HIS_upper,ymin=summ$HIS_lower)+geom_abline(intercept=0,slope=1,col='red')

summ<-merge(Scer[,summary.func(fold.induction,'Scer'),by=vars1], HIS[,summary.func(fold.induction,'HIS'),by=vars1],by=vars1[2:4])
p3<-ggplot(summ,aes(Scer_mean,HIS_mean))+geom_point()+geom_errorbarh(xmax= summ$Scer_upper,xmin=summ$Scer_lower)+geom_errorbar(ymax=summ$HIS_upper,ymin=summ$HIS_lower)+geom_abline(intercept=0,slope=1,col='red')

summ<-merge(Scer[,summary.func(growth.rate,'Scer'),by=vars1], HIS[,summary.func(growth.rate,'HIS'),by=vars1],by=vars1[2:4])
p4<-ggplot(summ,aes(Scer_mean,HIS_mean))+geom_point()+geom_errorbarh(xmax= summ$Scer_upper,xmin=summ$Scer_lower)+geom_errorbar(ymax=summ$HIS_upper,ymin=summ$HIS_lower)+geom_abline(intercept=0,slope=1,col='red')

ggsave(paste0(figout,date,'-KeyStatComparisonsGALK-ScervsHIS.pdf'),do.call(grid.arrange,list(p1,p2,p3,p4,top='mean expression galactose, fraction ON in galactose, the fold induction, and growth rate',ncol=4)), height=h,width=w)




###################################################
### epistasis from linear model
###################################################

lv<-as.logical(paste0(unlist(DT0[,lapply(.SD,is.logical)])))
ids1<-DT0[,colnames(DT0)[lv][-1],with=F]
ids2<-ids1[,!duplicated(colnames(ids1)),with=F]
ids<-ids2[,!'biol.rep',with=F]

ids<-DT0[,c('GAL3.clone.named','GAL80.clone.named','GAL4.clone.named','GALK.clone.named')]

DAT1<-data.table(x=DT0$growth.rate,ids)
summx<-DAT1[,summary.func(x),by=c(colnames(ids))]
boo<-(merge(DAT1,summx,by=c(colnames(ids))))
totvar<-varexplained(boo$x,boo$mean)
qplot(boo$mean,boo$x)+geom_abline(intercept=0,slope=1)+geom_vline(xintercept=credible.growth.rate.background,col='grey70')+
	geom_hline(yintercept=credible.growth.rate.background,col='grey70')+xlab('within-genotype means')+ylab('observed across N=28 observations')+
	xlim(c(0,0.25))+ylim(c(0,0.25))+ggtitle('overview of technical variation in higher-order epistasis experiment')+
	annotate('text',x=0.07,y=0.20,label=paste0('total variation explained =\n',round(totvar,3)))

mod1<-lm(x~.^4, DAT1)
summary(mod1)$coefficients
accuracy(mod1)
qplot(mod1$fitted.values, DAT1$x)+geom_abline()
varexplained(DAT1$x,mod1$fitted.values)

mod1df<-as.data.frame(summary(mod1)$coefficients)
mod1df$genotype<-rownames(mod1df)
DTmod1<-data.table(mod1df[c(2:nrow(mod1df),1),])
dfaov<-(as.data.frame(anova(mod1)))
dfaov$genotype<-rownames(dfaov)
DTaov<-data.table(dfaov)
ggplot(DTaov[order(var_exp)&!genotype=='(Intercept)'],aes(log10(var_exp),genotype))+geom_point(aes(colour=significant))+xlab('log10 variance explained')

DTmm<-na.exclude(merge(DTmod1,DTaov,by=c('genotype'),all=T))
DTmm[,p.adj.coef:=p.adjust(`Pr(>|t|)`,'fdr')]
DTmm[,p.adj.aov:=p.adjust(`Pr(>F)`,'fdr')]

DTmm[,var_exp:=sapply(`Sum Sq`,function(x)x/sum(`Sum Sq`,na.rm=T))]
DTmm[,genotype:=factor(genotype,levels=DTmm[order(var_exp)]$genotype)]
DTmm[,significant:=`p.adj.aov`<0.05]
DTmm[,genotype:=factor(genotype,levels=DTmm[order(Estimate)]$genotype)]
DTmm[,significant:=`p.adj.coef`<0.05]


ggplot(DTmm[order(Estimate)&!genotype=='Residuals'],aes(Estimate,genotype))+geom_point(aes(colour=significant))+xlab('mean effect')





###################################################
### calculate multiplicative epistasis: relative to "WT"
###################################################

groupings3a<-c('GALK.clone.named','GAL3.clone.named','GAL80.clone.named','GAL4.clone.named')

qDT<-DT0[,c(groupings3a,'growth.rate'),with=F]
qDT[,clone.named:=paste0(GALK.clone.named,GAL3.clone.named,GAL80.clone.named,GAL4.clone.named)]
groupings3<-c('clone.named','GALK.clone.named','GAL3.clone.named','GAL80.clone.named','GAL4.clone.named')
backg <-qDT[GAL4.clone.named=='GAL4.delta',summary.func(growth.rate)] # background growthrate
qDT$gr.actual<-as.numeric(as.character(DT0$growth.rate))
qDT[,growth.rate:=gr.actual-backg$mean]
plot(qDT$growth.rate,qDT$gr.actual);abline(a=0,b=1)


qDT3<-qDT[,c(groupings3,'growth.rate','gr.actual'),with=F]
DT3<-na.exclude(qDT3[,summary.func(growth.rate),by=groupings3][,!c('upper','lower')])

refvar<-c('GALK.Sac_cer','GAL3.WT','GAL80.WT','GAL4.WT')
WT<-DT3[GALK.clone.named==refvar[1]&GAL3.clone.named==refvar[2]&GAL80.clone.named==refvar[3]&GAL4.clone.named==refvar[4]]
asdf<-WT$clone.named
#setnames(WT,varis,paste0('WTexpect_',varis))


varis<-c('mean','sd','CoV','N')

dGKa<-DT3[GAL3.clone.named==refvar[2]&GAL80.clone.named==refvar[3]&GAL4.clone.named==refvar[4],!c(groupings3)[c(groupings3)!='GALK.clone.named'],with=F]
setnames(dGKa,varis,paste0('GALKexpect_',varis))
dG3a<-DT3[GALK.clone.named==refvar[1]&GAL80.clone.named==refvar[3]&GAL4.clone.named==refvar[4],!c(groupings3)[c(groupings3)!='GAL3.clone.named'],with=F]
setnames(dG3a,varis,paste0('GAL3expect_',varis))
dG80a<-DT3[GALK.clone.named==refvar[1]&GAL3.clone.named==refvar[2]&GAL4.clone.named==refvar[4],!c(groupings3)[c(groupings3)!='GAL80.clone.named'],with=F]
setnames(dG80a,varis,paste0('GAL80expect_',varis))
dG4a<-DT3[GALK.clone.named==refvar[1]&GAL3.clone.named==refvar[2]&GAL80.clone.named==refvar[3],!c(groupings3)[c(groupings3)!='GAL4.clone.named'],with=F]
setnames(dG4a,varis,paste0('GAL4expect_',varis))

merge_recurse2<-function(qDT,dt.list,by.list){
	if(length(dt.list)!=(length(by.list)))stop('list lengths look funny')
	x<-qDT
	i<-1
	while(i <= length(dt.list)){
		x<-merge(x,dt.list[[i]],by=by.list[[i]],all.x=T)
		i<-i+1
	}
	x
}
dt.list<-(list(dGKa,dG3a,dG80a,dG4a))
by.list<-(groupings3[-1])
merge(qDT,dt.list[[1]],by=by.list[[1]],all.x=T)
expectset<-merge_recurse2(qDT,dt.list,by.list)[!duplicated(clone.named),!c('growth.rate','gr.actual')]
	
xdt<-expectset[1,]
multiplicativemodel<-data.table(ldply(lapply(1:nrow(expectset),function(i){
	xdt<-expectset[i,]
	meannorm<-grepcols(xdt,'_mean')/WT$mean
	setnames(meannorm,colnames(meannorm),gsub('expect_mean','norm_mean',colnames(meannorm)))
	meannorm$expected.norm.multiplicative_mean<-prod(meannorm)
	meannorm$expected.gr.multiplicative_mean<-meannorm$expected.norm.multiplicative_mean*WT$mean
	meannorm$expected.gr.actual.multiplicative_mean<-meannorm$expected.gr.multiplicative_mean +backg$mean

	stdev<-sqrt((grepcols(xdt,'CoV')^2+WT$CoV^2))*abs(grepcols(meannorm,'norm_mean'))
	setnames(stdev,colnames(stdev),gsub('norm_mean','norm_sd',colnames(stdev)))
	stdev$expected.norm.multiplicative_sd<-sqrt(sum((grepcols(stdev,'norm_sd')/grepcols(meannorm,'norm_mean'))^2))*meannorm$expected.norm.multiplicative_mean
	stdev$expected.gr.multiplicative_sd<-sqrt(sum((stdev$expected.norm.multiplicative_sd/meannorm$expected.norm.multiplicative_mean)^2+(WT$sd/WT$mean)^2))*abs(meannorm$expected.gr.multiplicative_mean)
	stdev$expected.gr.actual.multiplicative_sd<-sqrt(stdev$expected.gr.multiplicative_sd^2 + backg$sd^2)

#	stderr.pre<-(sqrt(((grepcols(xdt,'sd')/sqrt(grepcols(xdt,'N')))/grepcols(xdt,'mean'))^2))
#	stderr<-data.table(stderr.pre*meannorm)
#	setnames(stderr,colnames(stderr),gsub('norm_mean','norm_stderr',colnames(stderr)))
	cbind(xdt,meannorm,stdev)
})))

summ<-qDT3[,mean(gr.actual),by=groupings3]
quick<-merge(summ,multiplicativemodel,by=groupings3)[,]
ggplot(quick,aes(x= expected.gr.actual.multiplicative_mean,y=V1))+geom_point()+geom_abline(intercept=0,slope=1,col='red')
varexplained(quick$V1 ,quick$expected.gr.actual.multiplicative_mean)
###################################################
### calculate multiplicative epistasis: mean across all backgrounds
###################################################

DT0[,clone.named:=paste0(GALK.clone.named,GAL3.clone.named,GAL80.clone.named,GAL4.clone.named)]
groupings3<-c('clone.named','GALK.clone.named','GAL3.clone.named','GAL80.clone.named','GAL4.clone.named')
qDT<-DT0[,c(groupings3,'growth.rate'),with=F]


backg <-qDT[GAL4.clone.named=='GAL4.delta',summary.func(growth.rate)] # background growthrate
qDT$gr.actual<-as.numeric(as.character(DT0$growth.rate))
qDT[,growth.rate:=gr.actual-backg$mean]
plot(qDT$growth.rate,qDT$gr.actual);abline(a=0,b=1)


qDT3<-qDT[,c(groupings3,'growth.rate','gr.actual'),with=F]
DT3<-na.exclude(qDT3[,summary.func(growth.rate),by=groupings3][,!c('upper','lower')])

refvar<-c('GALK.Sac_cer','GAL3.WT','GAL80.WT','GAL4.WT')
WT<-DT3[GALK.clone.named==refvar[1]&GAL3.clone.named==refvar[2]&GAL80.clone.named==refvar[3]&GAL4.clone.named==refvar[4]]
asdf<-WT$clone.named
#setnames(WT,varis,paste0('WTexpect_',varis))


varis<-c('mean','sd','CoV','N')

dGKa<-qDT[,summary.func(growth.rate),by='GALK.clone.named'][,!c('upper','lower'),with=F]
setnames(dGKa,varis,paste0('GALKexpect_',varis))
dG3a<-qDT[,summary.func(growth.rate),by='GAL3.clone.named'][,!c('upper','lower'),with=F]
setnames(dG3a,varis,paste0('GAL3expect_',varis))
dG80a<-qDT[,summary.func(growth.rate),by='GAL80.clone.named'][,!c('upper','lower'),with=F]
setnames(dG80a,varis,paste0('GAL80expect_',varis))
dG4a<-qDT[,summary.func(growth.rate),by='GAL4.clone.named'][,!c('upper','lower'),with=F]
setnames(dG4a,varis,paste0('GAL4expect_',varis))

merge_recurse2<-function(qDT,dt.list,by.list){
	if(length(dt.list)!=(length(by.list)))stop('list lengths look funny')
	x<-qDT
	i<-1
	while(i <= length(dt.list)){
		x<-merge(x,dt.list[[i]],by=by.list[[i]],all.x=T)
		i<-i+1
	}
	x
}
dt.list<-(list(dGKa,dG3a,dG80a,dG4a))
by.list<-(groupings3[-1])
summ<-qDT[,summary.func(growth.rate),by=groupings3]
expectset<-merge_recurse2(summ,dt.list,by.list)
expected<-apply(grepcols(expectset,'expect_mean')/WT$mean,1,prod)*WT$mean
qplot(expected,expectset$mean)+geom_abline(intercept=0,slope=1)+ylab('observed growth rate')+xlab('predicted growth rate based on background-averaged effects for each genotype')

colnames(grepcols(expectset,'expect_mean'))
mod1<-lm(mean~GALKexpect_mean+GAL3expect_mean+GAL80expect_mean+GAL4expect_mean,expectset)
mod1<-lm(mean~.,grepcols(expectset,'mean'))
summary(mod1)
qplot(mod1$fitted.values,expectset$mean)

testset<-data.table(t(apply(grepcols(expectset,'mean')/WT$mean,1,function(x)log(abs(x)))))
mod2<-lm(mean~.,testset)
summary(mod2)
qplot(mod2$fitted.values,testset$mean)

###################################################
### calculate multiplicative considering all possible single backgrounds "WT"
###################################################


groupings3<-c('clone.named','GALK.clone.named','GAL3.clone.named','GAL80.clone.named','GAL4.clone.named')
qDT<-DT0[,c(groupings3[-1],'growth.rate'),with=F]

qDT[,clone.named:=paste0(GALK.clone.named,GAL3.clone.named,GAL80.clone.named,GAL4.clone.named)]
backg <-qDT[GAL4.clone.named=='GAL4.delta',summary.func(growth.rate)] # background growthrate
qDT$gr.actual<-as.numeric(as.character(DT0$growth.rate))
qDT[,growth.rate:=gr.actual-backg$mean]
plot(qDT$growth.rate,qDT$gr.actual);abline(a=0,b=1)


qDT3<-qDT[,c(groupings3,'growth.rate','gr.actual'),with=F]
DT3<-na.exclude(qDT3[,summary.func(growth.rate),by=groupings3][,!c('upper','lower')])
stateset<-data.frame(qDT3[!duplicated(clone.named),lapply(.SD,function(x)paste(as.character(x))),by='clone.named',.SDcols=groupings3[groupings3!='clone.named']][,!'clone.named'])


merge_recurse2<-function(qDT,dt.list,by.list){
	if(length(dt.list)!=(length(by.list)))stop('list lengths look funny')
	x<-qDT
	i<-1
	while(i <= length(dt.list)){
		x<-merge(x,dt.list[[i]],by=by.list[[i]],all.x=T)
		i<-i+1
	}
	x}
		

MultiplicativeModeling<-lapply(1:nrow(stateset),function(i){	
	
	refvar<-unlist(stateset[i,])
	WT<-DT3[GALK.clone.named==refvar[1]&GAL3.clone.named==refvar[2]&GAL80.clone.named==refvar[3]&GAL4.clone.named==refvar[4]]
	asdf<-WT$clone.named
	varis<-c('mean','sd','CoV','N')
	dGKa<-DT3[GAL3.clone.named==refvar[2]&GAL80.clone.named==refvar[3]&GAL4.clone.named==refvar[4],!c(groupings3)[c(groupings3)!='GALK.clone.named'],with=F]
	setnames(dGKa,varis,paste0('GALKexpect_',varis))
	dG3a<-DT3[GALK.clone.named==refvar[1]&GAL80.clone.named==refvar[3]&GAL4.clone.named==refvar[4],!c(groupings3)[c(groupings3)!='GAL3.clone.named'],with=F]
	setnames(dG3a,varis,paste0('GAL3expect_',varis))
	dG80a<-DT3[GALK.clone.named==refvar[1]&GAL3.clone.named==refvar[2]&GAL4.clone.named==refvar[4],!c(groupings3)[c(groupings3)!='GAL80.clone.named'],with=F]
	setnames(dG80a,varis,paste0('GAL80expect_',varis))
	dG4a<-DT3[GALK.clone.named==refvar[1]&GAL3.clone.named==refvar[2]&GAL80.clone.named==refvar[3],!c(groupings3)[c(groupings3)!='GAL4.clone.named'],with=F]
	setnames(dG4a,varis,paste0('GAL4expect_',varis))
	
	dt.list<-(list(dGKa,dG3a,dG80a,dG4a))
	by.list<-(groupings3[-1])
	merge(qDT,dt.list[[1]],by=by.list[[1]],all.x=T)
	expectset<-merge_recurse2(qDT,dt.list,by.list)[!duplicated(clone.named),!c('growth.rate','gr.actual')]
		
	multiplicativemodel<-data.table(ldply(lapply(1:nrow(expectset),function(i){
		xdt<-expectset[i,]
		meannorm<-grepcols(xdt,'_mean')/WT$mean
		setnames(meannorm,colnames(meannorm),gsub('expect_mean','norm_mean',colnames(meannorm)))
		meannorm$expected.norm.multiplicative_mean<-prod(meannorm)
		meannorm$expected.gr.multiplicative_mean<-meannorm$expected.norm.multiplicative_mean*WT$mean
		meannorm$expected.gr.actual.multiplicative_mean<-meannorm$expected.gr.multiplicative_mean +backg$mean
	
		stdev<-sqrt((grepcols(xdt,'CoV')^2+WT$CoV^2))*abs(grepcols(meannorm,'norm_mean'))
		setnames(stdev,colnames(stdev),gsub('norm_mean','norm_sd',colnames(stdev)))
		stdev$expected.norm.multiplicative_sd<-sqrt(sum((grepcols(stdev,'norm_sd')/grepcols(meannorm,'norm_mean'))^2))*meannorm$expected.norm.multiplicative_mean
		stdev$expected.gr.multiplicative_sd<-sqrt(sum((stdev$expected.norm.multiplicative_sd/meannorm$expected.norm.multiplicative_mean)^2+(WT$sd/WT$mean)^2))*abs(meannorm$expected.gr.multiplicative_mean)
		stdev$expected.gr.actual.multiplicative_sd<-sqrt(stdev$expected.gr.multiplicative_sd^2 + backg$sd^2)
		cbind(xdt,meannorm,stdev)
	})))
	
	
	quick<-merge(qDT3,multiplicativemodel,by=groupings3)[,]
	mod<-lm(gr.actual~ expected.gr.actual.multiplicative_mean,quick)
	list(refvar,summary(mod)$r.squared,mod,quick)
})
# save(MultiplicativeModeling,file=paste0(outputdata,'180403-180413-ModelingDifferentBackgroundsMultiplicativeModel.rData'))
load(paste0(outputdata,'180403-180413-ModelingDifferentBackgroundsMultiplicativeModel.rData'))


test<-data.table(ldply(lapply(MultiplicativeModeling,function(x)c(x[[1]],rsqd=as.numeric(x[[2]])))))
test$rsqd<-as.numeric(test$rsqd)
# graph of best-performing background datasets
ggplot(test,aes(x=GAL80.clone.named,y=GAL4.clone.named))+geom_tile(aes(fill=rsqd))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'))+
	facet_wrap(~GAL3.clone.named+ GALK.clone.named,ncol=4)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))
refvar<-c('GALK.Sac_cer','GAL3.WT','GAL80.WT','GAL4.WT')


summ<-qDT[,summary.func(gr.actual),by=eval(groupings3[-1])]
#setnames(summ,'V1','mean.gr')
tmerge<-merge(summ,test,by=groupings3[-1])
ggplot(tmerge,aes(mean,CoV))+geom_point()
tmerge[,replicate.noise:=CoV/mean]
mtm<-melt(tmerge)
mtm$rsqd<-mtm[variable=='rsqd']$value
mtm$replicate.noise<-mtm[variable=='replicate.noise']$value
mtm2<-mtm[variable!='rsqd'&variable!='N']
ggplot(mtm2,aes(value,rsqd,col=replicate.noise))+geom_point()+facet_wrap(~variable,scales='free')+scale_colour_gradientn(colours=c('cornflowerblue','yellow','indianred'))
ggplot(mtm2,aes(value,rsqd,col=GAL80.clone.named))+geom_point()+facet_wrap(~variable,scales='free')#+scale_colour_gradientn(colours=c('cornflowerblue','yellow','indianred'))

mod1<-lm(rsqd~.,tmerge[,c('rsqd','mean','CoV','replicate.noise'),with=F])
summary(mod1)
mod2<-lm(rsqd~.,tmerge[,c('rsqd','mean'),with=F])
summary(mod2)
mod3<-lm(rsqd~.,tmerge[,c('rsqd','mean','GAL80.clone.named'),with=F])
summary(mod3)



best<-paste(unlist(test[rsqd==max(rsqd),groupings3[2:5],with=F]))
best<-refvar

grepl(best, lapply(MultiplicativeModeling,function(x)x[[1]]))

ldply(lapply(MultiplicativeModeling,function(x)as.logical(prod(best%in%x[[1]]))))
bestI<-MultiplicativeModeling[[145]] # 145 is WT, 138 is the best
DT<-bestI[[4]]

test0<-merge(DT,DT0[,summary.func(yfp.glu),by=c('clone.named')],by='clone.named')
mod0<-lm(gr.actual~ expected.gr.actual.multiplicative_mean *mean,test0)
summary(mod0)
test0$mod0<-mod0$fitted.values
test0[,deltaGAL80GAL3:=apply(data.frame(GAL3.clone.named,GAL80.clone.named,GAL4.clone.named),1,function(x){x[1]=='GAL3.delta'&x[2]=='GAL80.delta'&x[3]=='GAL4.WT'})]
qplot(mod0$fitted.values,test0$gr.actual,col=test0$deltaGAL80GAL3)
ggplot(test0,aes(mod0,gr.actual))+geom_point(aes(col=deltaGAL80GAL3))+geom_abline(intercept=0,slope=1,col='red')

ggplot(test0,aes(mod0,gr.actual))+geom_point(aes(col=deltaGAL80GAL3))+
	geom_point(data=test0,aes(expected.gr.actual.multiplicative_mean,gr.actual),col='blue',shape=21)+
	geom_abline(intercept=0,slope=1,col='red')

ggplot(test0,aes(expected.gr.actual.multiplicative_mean,gr.actual))+geom_point(aes(col=deltaGAL80GAL3))+geom_abline(intercept=0,slope=1,col='red')


bestI<-MultiplicativeModeling[[138]]
DT<-bestI[[4]]

test2<-data.table(ldply(lapply(MultiplicativeModeling,function(x)c(x[[1]],rsqd=as.numeric(x[[2]])))))
test$rsqd<-as.numeric(test$rsqd)

test2<-data.table(ldply(lapply(MultiplicativeModeling,function(x){
	DT<-x[[4]]
	c(x[[1]],rmsd=sqrt(mean(DT$gr.actual-DT$expected.gr.multiplicative_mean)^2))
	})))
test3<-merge(test,test2,by=groupings3[2:5])
test3$rmsd<-as.numeric(test3$rmsd)
ggplot(test3,aes(x=GAL80.clone.named,y=GAL4.clone.named))+geom_tile(aes(fill=log10(rmsd)))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'))+
	facet_wrap(~GAL3.clone.named+ GALK.clone.named,ncol=4)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggplot(test3[rmsd<1],aes(rsqd,log10(rmsd)))+geom_point()

###################################################
### look at epistasis correlations across dataset
###################################################

load(paste0(outputdata,'180403-180413-ModelingDifferentBackgroundsMultiplicativeModel.rData'))
MultiplicativeModeling[[1]][[4]][,groupings3[-1],with=F]
cortab<-(((lapply(MultiplicativeModeling,function(x){
	
	xDT<-x[[4]]
	xDT[,epistasis:=gr.actual-expected.gr.actual.multiplicative_mean]
	asdf<-xDT[,mean(epistasis),by=groupings3][,c('clone.named','V1'),with=F]
	colnames(asdf)<-c('clone.named',paste0(x[[1]],collapse=''))
	asdf
}))))
dt<-cortab[[1]][,'clone.named',with=F]
cortab0<-merge_recurse2(dt,cortab,lapply(1:length(cortab),function(x){return('clone.named')}))
cortab1<-cortab0
nameset<-cortab0$clone.named
colnames(cortab0)<-c('clone.named',nameset)
df1<-data.frame(cortab1);
rownames(df1)<-nameset
colnames(df1)<-c('clone.named',nameset)
library(plotly)
library(reshape2)

# Do hierarchical clustering

clust <- df1 %>%  dist() %>%  hclust()

# Get order
ord <- clust$order

# Re-arrange based on order
df <- df1[ord,]
# melt
mdf<-melt(df) 
# plot
ggplot(mdf,aes(clone.named,variable))+geom_tile(aes(fill=log10(abs(value))))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'))
mdf[abs(value)<100]


###################################################
### correct for error that is unreasonably low in observed growth rate
###################################################
DAT<-na.exclude(qDT3[,summary.func(gr.actual),by=c('plasgeno','binary.GAL3','binary.GAL80','binary.GAL4')])
DAT1<-merge(qDT3[,c('gr.actual','plasgeno','binary.GAL3','binary.GAL80','binary.GAL4')],DAT,by=c('plasgeno','binary.GAL3','binary.GAL80','binary.GAL4'))
plot(DAT1[,c('gr.actual','CoV')],log='xy')
DTsub<-DAT1[,c('gr.actual','mean','CoV')]
mod1<-lm(log(CoV)~log(mean), DAT1[,c('mean','CoV')])
summary(mod1)
# mod2<-nls(DTsub$CoV ~ (DTsub$growth.rate ^ b), start = c(b = -0.75), trace = T) # doesn't work
int<-(mod1$coefficients[1])
slope<-(mod1$coefficients[2])

plot(DTsub[,lapply(.SD,log)]);abline(a=int,b=slope)
DTnew<-data.table(DAT1[,c('plasgeno','binary.GAL3','binary.GAL80','binary.GAL4')],DTsub[,lapply(.SD,log)],fit=mod1$fitted.values,resid=mod1$residuals)
hist(DTnew$resid)
cutoff<-(-sd(DTnew$resid)*1.96)
DTnew[,lowerr:=CoV<fit]
DTnew[,lowresid:=resid<cutoff]
plot(DTnew[,c('mean','CoV')],col=as.factor(DTnew$lowresid))
DTnew[lowresid==T,length(gr.actual),by='plasgeno']
DTnew[,sdFit:=exp(fit)*exp(gr.actual)]


sdCorrected<-as.character(DTnew$lowresid)
sdCorrected[sdCorrected=='TRUE']<-DTnew[lowresid ==T,sdFit]
sdCorrected[sdCorrected=='FALSE']<-DTnew[lowresid ==F,exp(mean)*exp(CoV)]
DTnew$sdCorrected<-as.numeric(sdCorrected)
DATout<-merge(DAT,DTnew[,c('plasgeno','binary.GAL3','binary.GAL80','binary.GAL4','sdCorrected')],by=c('plasgeno','binary.GAL3','binary.GAL80','binary.GAL4'),all.y=F)[!duplicated(plasgeno)]
plot(DATout$sd,DATout$sdCorrected)
DATout$sd.old<-DATout$sd
DATout$sd<-DATout$sdCorrected
observed<-DATout
setnames(observed,colnames(observed)[-(1:4)],paste0('observed_',colnames(observed)[-(1:4)]))


###################################################
### analyze multiplicative epistasis
###################################################

DAT2<-merge(observed,multiplicativemodel,by=c('binary.GAL3','binary.GAL80','binary.GAL4'))
# minimum number of observations for each expectation for single mutants
expected_N<-min(unlist(lapply(list(WT,dG3,dG80,dG4),function(x){
	min(grepcols(x,'N'))})))

vars<-c('observed_mean','expected.gr.actual.multiplicative_mean','observed_sd','expected.gr.actual.multiplicative_sd','observed_N')

dt3<-as.data.frame(grepcols(DAT2,vars))
x<-dt3[1,]
ttests<-data.table(t(apply(DAT2[,c(vars),with=F],1,function(x){
	xa<-unlist(c(x,expected_N))
	names(xa)<-c('m1','m2','s1','s2','n1','n2')
	with(as.list(xa),{t.test2(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)})
})))
ttests$epistasis_p.adj<-p.adjust(ttests$p.value,'fdr')
oldcols<-c('diff.of.means','std.err','t','p.value')
newcols<-c('epistasis_mean','epistasis_stderr','epistasis_t','epistasis_p.value')
setnames(ttests,oldcols,newcols)

summarymerge<-merge(DAT2,data.table(DAT2[,c('binary.GAL3','binary.GAL80','binary.GAL4','plasgeno')],ttests),by=c('binary.GAL3','binary.GAL80','binary.GAL4','plasgeno'))
summarymerge$epistasis_significant<-summarymerge[,epistasis_p.adj<0.1&observed_mean>backg]
summarymerge$expected_stderr<-summarymerge$expected.gr.actual.multiplicative_sd/sqrt(32)
summarymerge$expected_mean_jitter<-jitter(summarymerge$expected.gr.actual.multiplicative_mean,120)
summarymerge[,expected_upper:= expected_mean_jitter + expected_stderr]
summarymerge[,expected_lower:= expected_mean_jitter-expected_stderr]
clusterset<-hdbscan(summarymerge[epistasis_significant==T,c('observed_mean','expected.gr.actual.multiplicative_mean'),with=F],minPts=3)
tmp<-summarymerge[epistasis_significant==T]
tmp$cluster<-as.factor(clusterset$cluster+1)
summarymerge2<-merge(summarymerge, tmp[,c('cluster','plasgeno'),with=F], by='plasgeno',all.x=T)

DT<-merge(qDT3, summarymerge2,by=c('binary.GAL3','binary.GAL80','binary.GAL4','plasgeno'))
DT$epistasis.multiplicative<-DT$gr.actual-DT$expected.gr.actual.multiplicative_mean
DT$residuals.lm.obs.vs.expected<-lm(gr.actual~expected.gr.actual.multiplicative_mean,data=DT)$residuals
DT$residuals.lm.allele.identities<-lm(gr.actual~.,data=DT[,c('gr.actual','binary.GAL4','binary.GAL80','binary.GAL3'),with=F])$residuals
DT$expected.lm.allele.identities<-lm(gr.actual~.,data=DT[,c('gr.actual','binary.GAL4','binary.GAL80','binary.GAL3'),with=F])$fitted.values
library(forecast)
a<-accuracy(lm(gr.actual~.,DT[,c('gr.actual','binary.GAL4','binary.GAL80','binary.GAL3'),with=F]))
b<-accuracy(x=DT$gr.actual,f=DT$expected.gr.actual.multiplicative_mean)
dm.test(DT$epistasis.multiplicative, DT$residuals.lm.allele.identities,alternative='greater')
plot(ecdf(abs(DT$residuals.lm.allele.identities)));lines(ecdf(abs(DT$epistasis.multiplicative)),col='red')
dt<-data.table(residuals=c(DT$epistasis.multiplicative, DT$residuals.lm.allele.identities),model=gl(2,nrow(DT),labels=c('multiplicative','linear')))

DT1<-merge(merge(qDT[,c('plasgeno','growth.rate','phenotypic.index','yfp.mean.glu','yfp.mean.gal','fracon.glu','fracon.gal')], tmp[,c('cluster','plasgeno'),with=F], by='plasgeno',all.x=T), summarymerge2[,c('plasgeno','epistasis_mean'),with=F],by='plasgeno')
DT1$yfp.mean.glu<-DT1$yfp.mean.glu/10000
DT1$yfp.mean.gal<-DT1$yfp.mean.gal/10000
DT1$cluster<-as.factor(DT1$cluster)
mDT1<-melt(DT1[!is.na(cluster)])







####################################
####OLD
####################################
# # 


# facy<-1#0.96
# facx<-0
# xlabels<-c('background','single mutants','double mutants')
# G4<-c('GAL4.WT','GAL4-L868P','GAL4-L868K','GAL4-L868G') # G4<-'GAL4-L868P'
# G80<-c('GAL80.WT','GAL80S-1')
# GK<-c('GALK.Sac_cer','GALK.Esc_col')
# G3<-c('GAL3.WT','GAL3.delta')
# ynudge<-0.02
# DT<-summ[aa%in%G4&cc%in%G80&dd%in%GK&bb%in%G3]
# arrowframe2<-arrowframe3[aa%in%G4&cc%in%G80]
# p1a<-ggplot(DT,aes(x,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+theme_minimal()+
	# geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
	# geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
	# geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	# geom_errorbar(col='grey50')+#geom_point(col='grey50')+
	# geom_text(data=DT,nudge_y= ynudge,aes(x=x,y=growth.rate_mean,label=shape2,col=(shape2)),size=2)+
# #	scale_colour_manual(values=c('grey50','orange1','cornflowerblue','darkgreen'),name='GAL sensor genotypes')+
	# scale_colour_manual(values=c('cornflowerblue','darkgreen','orange1','grey50'),name='GAL sensor genotypes')+
	# xlab('')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
	# theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+theme(panel.grid.minor.x = element_blank())+
	# scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+
	# facet_grid(cc~aa)	

# facy<-1#0.96
# facx<-0
# xlabels<-c('background','single mutants','double mutants')
# G4<-c('GAL4.WT') # G4<-'GAL4-L868P'
# G80<-c('GAL80.WT','GAL80.07')
# GK<-c('GALK.Sac_cer','GALK.Esc_col')
# G3<-c('GAL3.WT','GAL3.delta')
# ynudge<-0.02
# DT<-summ[aa%in%G4&cc%in%G80&dd%in%GK&bb%in%G3]
# arrowframe2<-arrowframe3[aa%in%G4&cc%in%G80]
# p1b<-ggplot(DT,aes(x,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+theme_minimal()+
	# geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
	# geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
	# geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	# geom_errorbar(col='grey50')+#geom_point(col='grey50')+
	# geom_text(data=DT,nudge_y= ynudge,aes(x=x,y=growth.rate_mean,label=shape2,col=(shape2)),size=2)+
# #	scale_colour_manual(values=c('grey50','orange1','cornflowerblue','darkgreen'),name='GAL sensor genotypes')+
	# scale_colour_manual(values=c('cornflowerblue','darkgreen','orange1','grey50'),name='GAL sensor genotypes')+
	# xlab('')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
	# theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+theme(panel.grid.minor.x = element_blank())+
	# scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+
	# facet_grid(cc~aa)	


# p1b<-ggplot(summ,aes(x,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+theme_minimal()+
	# geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
	# geom_ribbon(data=wtgrowth,aes(x,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
	# geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), 
		# arrow=arrow(type='closed',length=unit(0.3,'cm')),size=0.25,alpha=0.5,inherit.aes=F)+
	# geom_errorbar(col='grey50')+#geom_point(col='grey50')+
# #	geom_text(data=summ,nudge_y=0.04,aes(x=x,y=growth.rate_mean,label=shape2,col=(shape2)),size=2)+
	# scale_colour_manual(values=c('grey50','orange1','cornflowerblue','darkgreen'),name='GAL sensor genotypes')+
# #	scale_colour_manual(values=c('cornflowerblue','darkgreen','orange1','grey50'),name='GAL sensor genotypes')+
	# xlab('')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
	# theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.grid.minor.x = element_blank())+
	# scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3),labels=xlabels)+
	# facet_grid(cc~aa)	

	
# p2<-ggplot(summ,aes(yfp.mean.glu_mean,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower,xmax=yfp.mean.glu_upper,xmin=yfp.mean.glu_lower))+theme_minimal()+
	# geom_ribbon(data=nogrowth,aes(x*50000,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
	# geom_ribbon(data=wtgrowth,aes(x*50000,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
# #	geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	# geom_errorbar(col='black')+geom_errorbarh(col='black')+geom_point(aes(col=shape2))+
	# geom_text(data=summ,nudge_y=0.02,aes(x=yfp.mean.glu_mean,y=growth.rate_mean,label=shape2,col=(shape2)),size=2)+
	# scale_colour_manual(values=c('grey50','orange1','cornflowerblue','darkgreen'),name='GAL sensor genotypes')+
	# xlab('mean YFP expression in glucose')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
	# theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+
	# scale_y_continuous(breaks=seq(0.05,0.25,by=0.05),limits=c(0.05,0.27))+scale_x_continuous(limits=c(-3000,40000))+
	# facet_grid(cc~aa)	


# p3<-ggplot(summ,aes(yfp.mean.glu_mean,yfp.mean.gal_mean,ymax=yfp.mean.gal_upper,ymin=yfp.mean.gal_lower,xmax=yfp.mean.glu_upper,xmin=yfp.mean.glu_lower))+theme_minimal()+
	# # geom_ribbon(data=nogrowth,aes(x*50000,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=1,inherit.aes=F)+
	# # geom_ribbon(data=wtgrowth,aes(x*50000,ymax=upper,ymin=lower),fill='grey70',col='transparent',alpha=0.5,inherit.aes=F)+
	# # geom_segment(data= arrowframe2,aes(y=y0,yend=yend*facy,x=x0,xend=xend,col= shape3), arrow=arrow(type='closed',length=unit(0.3,'cm')),size=1,alpha=0.5,inherit.aes=F)+
	# geom_errorbar(col='black')+geom_errorbarh(col='black')+geom_point(aes(col=shape2))+
	# geom_text(data=summ,nudge_y=0.02,aes(x=yfp.mean.glu_mean,y=yfp.mean.gal_mean,label=shape2,col=(shape2)),size=2)+
	# scale_colour_manual(values=c('grey50','orange1','cornflowerblue','darkgreen'),name='GAL sensor genotypes')+
	# xlab('mean YFP expression in glucose')+ylab(expression(paste('growth rate ',mu,' ',hr^-1,sep=' ')))+
	# theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(panel.spacing = unit(1, "lines"))+
	# scale_y_continuous(limits=c(-3000,60000))+scale_x_continuous(limits=c(-3000,40000))+geom_abline(col='grey50')+
	# facet_grid(cc~aa)	

# # do.call(grid.arrange,list(p1,p2))
# sample(letters,4)
# pUEHX<-list(p1,p2)
# w<-10;h<-8.1
# ggsave(paste0(figout,date,'-180831-EffectOfGAL3AndGAL1FlippedInNewInducibleBackgrounds.png'),do.call(grid.arrange,pUEHX),width=w,height=h)

# ggplot(summ,aes(x,phenotypic.index_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+
	# #geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),col='grey70',alpha=0.3,inherit.aes=F)+
	# ylim(c(0,2))+xlim(c(0.5,3.5))+geom_hline(yintercept=0,col='cornflowerblue')+geom_hline(yintercept=1,col='yellow')+geom_hline(yintercept=2,col='indianred')+
	# geom_text(aes(label=shape2,col=(shape)))+facet_grid(cc~aa)

# ggplot(summ,aes(yfp.mean.glu_mean,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+
	# ylim(c(0.05,0.25))+xlim(c(-5000,40000))+ geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),col='grey70',alpha=0.3,inherit.aes=F) +
	# geom_text(aes(label=shape2,col=(shape)),size=3)+facet_grid(cc~aa)+theme(axis.text.x = element_text(angle = 30, hjust = 1))

# G4<-c('GAL4.WT','GAL4-L868P') # G4<-'GAL4-L868G'
# G80<-c('GAL80.WT','GAL80S-1')
# GK<-c('GALK.Sac_cer','GALK.Esc_col')
# G3<-c('GAL3.WT','GAL3.delta')
# dSub<-DT0[aa%in%G4&cc%in%G80&dd%in%GK&bb%in%G3]
# backg<-DT0[aa=='GAL4.delta'|dd=='HIS5.Sch_pom',summary.func(growth.rate)]
# nogrowth<-data.table(x=seq(0,4,length=100),upper=backg$upper,lower=backg$lower)
# summ<-summary.func.all(dSub,c('growth.rate','yfp.mean.gal'),c('aa','bb','cc','dd'))
# summ[dd=='GALK.Esc_col'&bb=='GAL3.delta',c('x','shape','shape2'):=list(1,factor(1),'no sensor')]
# summ[dd=='GALK.Sac_cer'&bb=='GAL3.delta',c('x','shape','shape2'):=list(1.9,factor(6),'GAL1')]
# summ[dd=='GALK.Esc_col'&bb=='GAL3.WT',c('x','shape','shape2'):=list(2.1,factor(2),'GAL3')]
# summ[dd=='GALK.Sac_cer'&bb=='GAL3.WT',c('x','shape','shape2'):=list(3,factor(11),'GAL1+GAL3')]

# ggplot(summ,aes(x,growth.rate_mean,ymax=growth.rate_upper,ymin=growth.rate_lower))+
	# geom_ribbon(data=nogrowth,aes(x,ymax=upper,ymin=lower),col='grey70',alpha=0.3,inherit.aes=F)+ylim(c(0.05,0.25))+xlim(c(0.5,3.5))+
	# geom_text(data=summ,aes(x=x,y=growth.rate_mean,label=shape2,col=(shape)),size=3)+facet_grid(cc~aa)

