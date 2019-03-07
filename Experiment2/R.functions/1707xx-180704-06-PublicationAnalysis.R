

# this code was used to generate figures related to experiment 1 for the manuscript
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
### analysis: illustrative plots of initial filtering done on smaples based on size and shape characteristics determined by FCS and SSC channels.
###################################################
# The majority of this code is for analyzing already-processed FCS data. see code in "preprocessing" scripts for that. 
# For publication, I made an illustrative example of how a typical FCS file was analyzed

setwd(plates.dir)

dirs<-list.dirs()[grep('Plate',list.dirs())][1]
ix<-1 # this is normally a for loop which dumps a bunch of filtered flowCore data files in .rData format
wash=NULL # use if a blank sample was taken before every measurement
filtsize<-4 # filter size set for the curv2Filter
setwd(plates.dir) # directory where the plates are
setwd(dirs[ix]) # move to the directory where the plates are
dir<-paste(dirs[ix]) # make a list of all the directories with plates
setwd(plates.dir);setwd(dir) # move to the first directory in the list
fcs.files<-list.files(pattern='../*.fcs') # get all teh fcs files in the folder
# fsapply is normally used to go through each sample
flowframes<-read.flowSet(fcs.files, name.keyword="TUBE NAME", phenoData=list(name="TUBE NAME", Filename="$FIL")) # make a flowframe object of the fcs files
#	head(pData(phenoData(flowframes)))
fcs<-flowframes[[1]] # within fsapply, go through each item in the flow frame
if(!is.null(wash)){strt<-0} # establish a starting cell population 
if(is.null(wash)){strt<-500} # set a variable of number of cells to discard if the previous reading was not washed
len<-length(exprs(fcs$'FSC-A')) # determine the length of the list acquired
if(len<(strt+1000)){filteredPlus<-NULL # if you have fewer than 1000 cells plus the amount you are discarding then count the sample as a NULL
	}else{
	ssc<-quantile(exprs(fcs[strt:len,3]),c(0.1,0.9)) # determine the quantiles to filter the ssc
	fsc<-quantile(exprs(fcs[strt:len,1]),c(0.1,0.9)) # determine the quantiles to filter the fsc
	rectGate <- rectangleGate("FSC-H"=fsc,"SSC-H"=ssc) # create a gate for the samples based on these quantiles
	filteredRec<-split(fcs,flowCore::filter(fcs,rectGate))$`defaultRectangleGate+` # actually filter the samples based on the rectangular gate
	c2f <- curv2Filter("FSC-H","SSC-H",bwFac=filtsize) # take rectangle filtered samples and apply the curve2Filter gate
	filtered<-split(filteredRec,flowCore::filter(filteredRec,c2f)) # actually filter the samples based on the curve2Filter
	filteredPlus<-filtered$`area 1`} # extract the cells that wall within the central area of the gate

	filteredPlus # return the final filtered sample

# after filtering, samples are re-analyzed for expression distribution characteristics
# make a contour plot, which I guess is a publication standard for 2-D representations of flow cytometry data.

# put data together
DT0<-data.table(FSC=exprs(fcs$'FSC-A'),SSC=exprs(fcs$'SSC-A'),fac='original')
DTrect<-data.table(FSC=exprs(filteredRec $'FSC-A'),SSC=exprs(filteredRec $'SSC-A'),fac='rectangle filtered')
DTcurve<-data.table(FSC=exprs(filteredPlus $'FSC-A'),SSC=exprs(filteredPlus $'SSC-A'),fac='curve2Filt filtered')
filtlevs<-c('original','rectangle filtered','curve2Filt filtered')
DT<-rbind(DT0,DTrect,DTcurve)
DT[,`filtered data`:=factor(fac,levels=filtlevs)]
setnames(DT,c('FSC.FSC-A', 'SSC.SSC-A'),c('FSC','SSC'))
DT[,c('log10FSC','log10SSC'):=list(log10(FSC+1000)-3,log10(SSC+1000)-3)]

# make a 2d kernel density estimate of the two dimensions
library(sparr)
X<-log10(as.numeric(as.vector(unlist(DT0[,1])))+1000)-3
Y<-log10(as.numeric(as.vector(unlist(DT0[,2])))+1000)-3
dtp<-ppp(X,Y,window=owin(range(X),range(Y)))
test<-bivariate.density(dtp,h0=0.05)
zmat<-as.matrix(test$z)
x<-test$q$xcol
y<-test$q$yrow
DTsm<-data.table(Y=y,zmat)
DTm<-melt(DTsm,id='Y')
dum<-data.table(variable=paste0('V',seq(1,128)),X=x)
DTf<-merge(dum,DTm,by='variable')
sample(letters,4)

# plot filtered samples
p<-ggplot(DT,aes(log10FSC,log10SSC,col=`filtered data`))
pIUDX<-p+geom_point(size=0.3,shape=21,stroke=0.1)+geom_contour(data=DTf,aes(x=X,y=Y,z=value),inherit.aes=F,col='black',alpha=0.2)+
	theme_NewPub+
	xlab('log10 FSC signal')+ylab('log10 SSC signal')

if(writeplot==T){
	w<-3.9;h<-2.85
	ggsave(paste0(figout,'17092x-190110-pIUDX-FilteringSSC,FSC,illustration.pdf'), pIUDX,width=w,height=h)
}


###################################################
### load data
###################################################
DT00<-fread(paste0(outputdata,'1707xx-180704-OutliersRemoved.txt'),sep='\t',header=T)

# code to generate this file with clusters is below.
load(paste0(layout.dir,'/1707xx-181120-Clustered.rdata'))

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
### basic phenotype X-Y plotting: fracon glu vs fracon gal. 
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
	theme_NewPub(arrow=T)+
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
### all expression distribution phenotypes
###################################################

ylimits<-c(0,0.35)
genos<-c('aa','bb','cc','pairing')
denscolconv<-function(DT, vnull =T){
	DT[,variable:=gsub('dens.','',variable)]
	DT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
	DT[,sugar:=paste0('g',sugar)]
	if(vnull ==T)DT[,variable:=NULL]
}

DTx<-copy(densDT[!is.na(clone.named)])
DTX<-merge(DTx,clusterframe0,by='clone.named')

DTX[,pairing:=mutant_gene]
WTFrame<-DTX[WT==T]
a<-copy(WTFrame);b<-copy(WTFrame);cx<-copy(WTFrame);aa<-copy(WTFrame);bb<-copy(WTFrame);cxc<-copy(WTFrame)
a[,pairing:='GAL80 GAL4'];b[,pairing:='GAL3 GAL80'];cx[,pairing:='GAL80 GAL4'];aa[,pairing:='GAL3 GAL4'];bb[,pairing:='GAL3 GAL4'];cxc[,pairing:='GAL3 GAL80']
DTX1<-rbind(DTX[WT!=T],a,b,cx,aa,bb,cxc)

singleframe<-DTX[single_lv==T]
a<-copy(singleframe[WT!=T]);b<-copy(singleframe[WT!=T]);cx<-copy(singleframe[WT!=T]);aa<-copy(singleframe[WT!=T]);bb<-copy(singleframe[WT!=T]);cxc<-copy(singleframe[WT!=T])
a[,pairing:='GAL80 GAL4'];b[,pairing:='GAL3 GAL80'];cx[,pairing:='GAL80 GAL4'];aa[,pairing:='GAL3 GAL4'];bb[,pairing:='GAL3 GAL4'];cxc[,pairing:='GAL3 GAL80']
DTX2<-rbind(DTX1[single_lv!=T],a,b,cx,aa,bb,cxc)

singleframe2<-merge(phenDT[single_lv==T],clusterframe0,by='clone.named')
singleframe2[grepl('GAL4',mutant_gene),gene:='GAL4']
singleframe2[grepl('GAL3',mutant_gene),gene:='GAL3']
singleframe2[grepl('GAL80',mutant_gene),gene:='GAL80']
summ<-singleframe2[,summary.func(phenotypic.index),by=c('mutants','cluster_Simple','gene')]
summ[,var1:=mean(mean),by=cluster_Simple]
G3<-c(c('GAL3.WT','GAL3.delta'),notin(gsub('GAL3.WT GAL80.WT GAL4.WT','GAL3.WT',summ[gene=='GAL3'|mutants=='GAL3.WT GAL80.WT GAL4.WT'][order(var1,mean)]$mutants),c('GAL3.WT','GAL3.delta')))
G80<-c(c('GAL80.WT','GAL80.delta'),notin(gsub('GAL3.WT GAL80.WT GAL4.WT','GAL80.WT',summ[gene=='GAL80'|mutants=='GAL3.WT GAL80.WT GAL4.WT'][order(var1,mean)]$mutants),c('GAL80.WT','GAL80.delta')))
G4<-c(c('GAL4.WT','GAL4.delta'),notin(gsub('GAL3.WT GAL80.WT GAL4.WT','GAL4.WT',summ[gene=='GAL4'|mutants=='GAL3.WT GAL80.WT GAL4.WT'][order(var1,mean)]$mutants),c('GAL4.WT','GAL4.delta')))


DTX2[,aa:=factor(aa,levels=G4)]
DTX2[,bb:=factor(bb,levels=G3)]
DTX2[,cc:=factor(cc,levels=G80)]


mDT0<-melt(DTX2[,c(dens.cols,genos,'cluster_Simple'),with=F],id=c(genos,'cluster_Simple'))#;mDT<-mDT[order(sugar,AU.FL,decreasing=T)]
mDT<-mDT0[, list(density_mean=mean(value,na.rm=T),density_sd=sd(value,na.rm=T)),by=c(genos,'variable')]
denscolconv(mDT)
mDT[,AU.FL:=as.numeric(gsub('_','',AU.FL))]
mDT[density_mean>=max(ylimits),density_mean:=max(ylimits)]

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
	#	scale_x_continuous(limits=xlimits,breaks=xbreaks,labels=xlabs)+
		theme(legend.position='none')
}

mDT1<-mDT#[aa%in%G4&cc%in%G80&bb%in%G3]
sample(LETTERS,4)
pXFWV_GAL80vGAL4<-ExpDensLinePlot(mDT1[pairing=='GAL80 GAL4'])+facet_grid(cc~aa)+theme(axis.title=element_text(size=30))
pXFWV_GAL3vGAL4<-ExpDensLinePlot(mDT1[pairing=='GAL3 GAL4'])+facet_grid(bb~aa)+theme(axis.title=element_text(size=30))
pXFWV_GAL3vGAL80<-ExpDensLinePlot(mDT1[pairing=='GAL3 GAL80'])+facet_grid(bb~cc)+theme(axis.title=element_text(size=30))

if(writeplot==T){
	w<-27;h<-27
#	ggsave(paste0(figout,'1707xx-181203-pXFWV_GAL80vGAL4-ExpressionLineDensityPlots-AllClones-Experiment2.pdf'), pXFWV_GAL80vGAL4,width=w,height=h,limitsize=F)
	ggsave(paste0(figout,'1707xx-181203-pXFWV_GAL80vGAL4-ExpressionLineDensityPlots-AllClones-Experiment2.png'), pXFWV_GAL80vGAL4,width=w,height=h,limitsize=F)
	ggsave(paste0(figout,'1707xx-181203-pXFWV_GAL3vGAL4-ExpressionLineDensityPlots-AllClones-Experiment2.png'), pXFWV_GAL3vGAL4,width=w,height=h,limitsize=F)
	ggsave(paste0(figout,'1707xx-181203-pXFWV_GAL3vGAL80-ExpressionLineDensityPlots-AllClones-Experiment2.png'), pXFWV_GAL3vGAL80,width=w,height=h,limitsize=F)
	
	
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
	mDT[,upper:=density_mean+density_sd]
	mDT[,lower:=density_mean-density_sd]
	mDT[density_mean>=max(ylimits),density_mean:=max(ylimits)-0.01]
	mDT[upper>=max(ylimits), upper:=max(ylimits)-0.01]
	mDT[lower>=max(ylimits), lower:=max(ylimits)-0.01]
	mDT[lower<=min(ylimits), lower:=min(ylimits)+0.01]
	mDT[upper<=min(ylimits), upper:=min(ylimits)+0.01]

	return(mDT)
}
ExpDensLinePlot<-function(mDT1,ylimits=c(0,0.35),linethickness=0.3){
	xlimits<-c(-0.3,3)
	xbreaks<-c(0,1,2)
	xlabs<-c('0','1','2')
	ybreaks<-c(0.0,0.1,0.2,0.3)
	ylabs<-c('0.0',paste0(ybreaks[2:3]),paste0('>',max(ybreaks)))
	ggplot(mDT1[order(sugar)],aes(AU.FL,density_mean))+#geom_point(data=mDT1,aes(AU.FL,value,col=sugar),size=0.1,alpha=0.3,shape=21)+
		theme_NewPub+
		geom_ribbon(aes(x=AU.FL,ymax=upper,ymin=lower,fill=sugar),alpha=0.3)+
		geom_line(aes(col=sugar),size=linethickness)+
		scale_colour_manual(values= factorTwoSugar,na.value='transparent')+
		ylab('density') + xlab('pseudo-log10 fluorescence intensity, A.U.')+
		scale_y_continuous(limits=ylimits,breaks=ybreaks,labels=ylabs)+
		scale_x_continuous(limits=xlimits,breaks=xbreaks,labels=xlabs)+
		theme(legend.position='none')
}

# first go at reducing outliers and properly categorizing expression distributions
# cluster 1 is the low max expression GAL4 backgrounds. many were mis-characterized, so go in here and give them names in cluster_new
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

# many unclustered clones are constitutitve-like. by clustering with the constitutive cluster 7 we are able to put them in the right broad category
# creating cluster_new2 carries these over
dSub1<-DTm2[cluster_new=='0'|cluster=='7',c(grepincols(DTm1,'dens__'),'cluster_new','clone.named'),with=F]
hdbtest0<-hdbscan(dSub1[, grepincols(dSub1,'dens__'),with=F],minPts=8)
dSub1$cluster_test<-hdbtest0$cluster
DTX<-merge(densDT[clone.named%in%dSub1$clone.named,c(grepincols(densDT,'dens__'),'clone.named'),with=F],dSub1[,.(clone.named,cluster_test)],by='clone.named')
mDT1<-densMelt(DTX)
p1<-ExpDensLinePlot(mDT1)+facet_grid(~cluster_test)
dSub1[,cluster_new2:=as.character(cluster_test)]
dSubx<-dSub1[cluster_new=='0',.(clone.named,cluster_new,cluster_test,cluster_new2)]
dSubx[cluster_test=='1',cluster_new2:='partially inducible']
dSubx[cluster_test=='2',cluster_new2:='constitutive']
DTm3<-merge(DTm2,dSubx[,.(cluster_new2,clone.named)],by='clone.named',all=T)
DTm3[is.na(cluster_new2),cluster_new2:=as.character(cluster_new)]
DTm3[cluster_new==0]

# many unclustered clones are inducible-like. by clustering with the inducible cluster 6 we are able to put them in the right broad category
# creating cluster_new3 carries these over
dSub1<-DTm3[cluster_new2=='0'|cluster=='6',c(grepincols(DTm1,'dens__'),'cluster_new2','clone.named'),with=F]
hdbtest0<-hdbscan(dSub1[, grepincols(dSub1,'dens__'),with=F],minPts=4)
dSub1$cluster_test<-hdbtest0$cluster
DTX<-merge(densDT[clone.named%in%dSub1$clone.named,c(grepincols(densDT,'dens__'),'clone.named'),with=F],dSub1[,.(clone.named,cluster_test)],by='clone.named')
mDT1<-densMelt(DTX)
p1<-ExpDensLinePlot(mDT1)+facet_grid(~cluster_test)
DTX2<-merge(densDT[clone.named%in%dSub1[cluster_new2==0]$clone.named,c(grepincols(densDT,'dens__'),'clone.named'),with=F],dSub1[,.(clone.named,cluster_test)],by='clone.named')
mDT2<-densMelt(DTX2)
p2<-ExpDensLinePlot(mDT2)+facet_wrap(~clone.named+cluster_test,ncol=5)
dSubx<-dSub1[cluster_new2=='0',.(clone.named,cluster_new2,cluster_test)]
dSubx[,cluster_new3:=as.character(cluster_test)]
dSubx[cluster_test=='0',cluster_new3:='const, low expr']
dSubx[cluster_test=='1',cluster_new3:='leaky glucose']
dSubx[cluster_test=='2',cluster_new3:='inducible']
DTm4<-merge(DTm3,dSubx[,cluster_new3,clone.named],by='clone.named',all=T)
DTm4[is.na(cluster_new3),cluster_new3:=as.character(cluster_new2)]
table(DTm4$cluster_new3)
table(DTm4$cluster)

# all done - ready to characterize the ones that originally were scored well
DTm4[,cluster_A:= cluster_new3]
# categorize all clusters
DTm4[cluster_new3=='0',cluster_A:='outlying']
DTm4[cluster_new3=='1',cluster_A:='weak expr, induc']
DTm4[cluster_new3%in%c('2','3','4'),cluster_A:='const, low expr']
DTm4[cluster_new3=='5',cluster_A:='uninducible']
DTm4[cluster_new3=='6',cluster_A:='inducible']
DTm4[cluster_new3=='7',cluster_A:='constitutive']
table(DTm4$cluster_new3)
table(DTm4$cluster_A)


############
### many remain mis-categorized by mean effects.
############
# by taking the log of expression and amplifying this signal with the log of the density we were able to pull out fine differences between expression distributions
# this means we're missing some gross differences evident at the the mean level.
# a simple statistic to take is the fold induction. as you can see, the log of the fold induction is bi- or tri-modally distributed within certain clusters.
# below i go through these and characterize them by hand to the right clusters.

# Example: GAL4-L868E + GAL80S-1 is misclassified as partially inducible when it's basically inducible

############
### sub-classify based on fold indcution
############


clusterframe0<-DTm4[,.(clone.named,aa,bb,cc,cluster_A)]
genos<-c('clone.named','aa','bb','cc')
clusterframe0[cluster_A=='constitutive']
DT0<-copy(phenDT)
DT0[,foldinduction:=yfp.mean.gal/yfp.mean.glu]
DT1<-merge(DT0,clusterframe0,by=genos)
summ<-DT1[,summary.func(foldinduction),by=c(genos,'cluster_A')]
# ggplot(summ,aes(y=log10(mean),x= cluster_A))+geom_violin()+geom_hline(yintercept=0.92)+geom_hline(yintercept=1.4)



############
### sub-classify constitutive
############
cutoffA<-0.92
cutoffB<-1.4
# ggplot(summ,aes(y=log10(mean),x= cluster_A))+geom_violin()+geom_hline(yintercept= cutoffA)+geom_hline(yintercept= cutoffB)

summ[cluster_A=='constitutive' &mean<10^cutoffA,cluster_test:='constitutive']
summ[cluster_A=='constitutive' &mean>=10^cutoffA&mean<10^cutoffB,cluster_test:='partially constitutive']
summ[cluster_A=='constitutive' &mean>=10^cutoffB,cluster_test:='leaky glucose']

DTA<-copy(densDT[clone.named%in%summ[!is.na(cluster_test)]$clone.named,c('clone.named',grepincols(densDT,'dens__')),with=F])
DTX<-merge(DTA,summ[,.(clone.named,cluster_test)],by='clone.named')
meltby<-c(genos,'cluster_test')
mDT1<-densMelt(DTX)
# ExpDensLinePlot(mDT1)+facet_grid(~cluster_test)

############
### sub-classify leaky glucose
############
cutoffA<-1.1

# ggplot(summ,aes(y=log10(mean),x= cluster_A))+geom_violin()+geom_hline(yintercept= cutoffA)

summ[cluster_A=='leaky glucose' &mean<10^1.1,cluster_test:='const, low expr']
summ[cluster_A=='leaky glucose' &mean>=10^1.1,cluster_test:='weak expr, induc']

DTA<-copy(densDT[clone.named%in%summ[!is.na(cluster_test)]$clone.named,c('clone.named',grepincols(densDT,'dens__')),with=F])
DTX<-merge(DTA,summ[,.(clone.named,cluster_test)],by='clone.named')
meltby<-c(genos,'cluster_test')
mDT1<-densMelt(DTX)
# ExpDensLinePlot(mDT1)+facet_grid(~cluster_test)

############
### sub-classify uninducible
############
cutoffA<-1.0
# ggplot(summ,aes(y=log10(mean),x= cluster_A))+geom_violin()+geom_hline(yintercept=cutoffA)

summ[cluster_A=='uninducible' &mean<10^cutoffA,cluster_test:='uninducible']
summ[cluster_A=='uninducible' &mean>=10^cutoffA,cluster_test:='partially inducible']

DTA<-copy(densDT[clone.named%in%summ[!is.na(cluster_test)]$clone.named,c('clone.named',grepincols(densDT,'dens__')),with=F])
DTX<-merge(DTA,summ[,.(clone.named,cluster_test)],by='clone.named')
meltby<-c(genos,'cluster_test')
mDT1<-densMelt(DTX)
# ExpDensLinePlot(mDT1)+facet_grid(~cluster_test)



############
### sub-classify const, low expr
############
cutoffA<-1.2
# ggplot(summ,aes(y=log10(mean),x= cluster_A))+geom_violin()+geom_hline(yintercept=cutoffA)

summ[cluster_A=='const, low expr' &mean<10^cutoffA,cluster_test:='const, low expr']
summ[cluster_A=='const, low expr' &mean>=10^cutoffA,cluster_test:='weak expr, induc']#

DTA<-copy(densDT[clone.named%in%summ[!is.na(cluster_test)]$clone.named,c('clone.named',grepincols(densDT,'dens__')),with=F])
DTX<-merge(DTA,summ[,.(clone.named,cluster_test)],by='clone.named')
meltby<-c('clone.named','cluster_test')
mDT1<-densMelt(DTX,meltby=meltby, summby = c("cluster_test", "sugar", "AU.FL"))
# ExpDensLinePlot(mDT1)+facet_grid(~cluster_test)



############
### check partially inducibles
############
# rationale here is that some partially inducible clones are obviuosly inducible, perhaps with longer lag phase, others are uninducible
clones<-summ[cluster_A =='partially inducible']$clone.named
summ1<-merge(phenDT[clone.named%in%clones,summary.func(growth.rate,'gr'),by=genos],phenDT[clone.named%in%clones,summary.func(fracon.gal),by=genos],by=genos)
clones2<-summ1[mean>=0.4]$clone.named
clones3<-summ1[mean<0.4&mean>=0.2]$clone.named
clones4<-summ1[mean<0.2]$clone.named
DTA<-copy(densDT[clone.named%in%clones2,c('clone.named',grepincols(densDT,'dens__')),with=F])
DTX<-merge(DTA,summ[,.(clone.named, cluster_test)],by='clone.named')
meltby<-c('clone.named','cluster_test')
mDT1<-densMelt(DTX,meltby=meltby, summby = c("cluster_test",'clone.named', "sugar", "AU.FL"))
# ExpDensLinePlot(mDT1)+facet_grid(~cluster_test)+facet_wrap(~clone.named)
summ[clone.named%in%clones2,cluster_test:='inducible']
summ[clone.named%in%clones3,cluster_test:='partially inducible']
summ[clone.named%in%clones4,cluster_test:='uninducible']


############
### check how things went
############

summ[,cluster_C:=cluster_test]
summ[is.na(cluster_test),cluster_C:=cluster_A]
bylist<-c('clone.named','cluster_C')
summ0<-summary.func.all(merge(summ[,.(clone.named,cluster_C)],phenDT[,log10.fold.induction:=log10(yfp.mean.gal/yfp.mean.glu)],by='clone.named'),
	c('growth.rate','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal','log10.fold.induction'),bylist=bylist)

mDT<-melt(summ0[!is.na(clone.named)],id=bylist)
mDT[,c('pheno','stat'):=colsplit(variable,'_',c('pheno','stat'))]
form<-castform(c(bylist,'pheno'),c('stat'))
cDT<-dcast(mDT,form)

sample(LETTERS,4)
ggplot(cDT,aes(y=mean,x= cluster_C))+geom_jitter(shape=21,width=0.1,size=0.5)+geom_violin(alpha=0.7)+
	theme_NewPub_noarrow+theme(axis.text.x=element_text(angle=30,hjust=1))+facet_wrap(~pheno,scales='free')+	
	xlab('cluster')+ylab('within-genotype mean values')

############
### THIS IS THE FINAL CLASSIFICATION
############

dum<-data.table(cluster_C=unique(summ$cluster_C),ord=c(3,5,6,4,7,1,2,8,9))[order(ord)]
dum[,cluster:=factor(cluster_C,levels=dum[order(ord)]$cluster_C)]

# i renamed the clusters by hand based on broad categories (uninducible, inducible, constitutive), and then more narrow sub-categories
cluster_annotated<-fread(paste0(layout.dir,'/1707xx-181119-ClusterOrdersByHand.txt'))
cluster_annotated[,cluster_Broad:=factor(cluster_Broad,levels=unique(cluster_annotated$cluster_Broad))]
cluster_annotated[,cluster_Narrow:=factor(cluster_Narrow,levels=unique(cluster_annotated$cluster_Narrow))]

clusterframe0<-merge(summ[,.(cluster=cluster_C,clone.named)],cluster_annotated,by='cluster')
clusterframe0[,cluster_Simple:=factor(cluster_Simple,levels=c('inducible','uninducible','constitutive','leaky','weak expr'))]


if(writeplot==T){
	save(clusterframe0,file=paste0(layout.dir,'/1707xx-181120-Clustered.rdata'))
}


###################################################
### plot based on clusters
###################################################

load(paste0(layout.dir,'/1707xx-181120-Clustered.rdata'))
# clusterframe0 # new object loaded

############
### plot phentoypic mean values across clusters
############

bylist<-c('clone.named','cluster_Narrow','cluster_Simple','cluster_Broad')
summ0<-summary.func.all(merge(clusterframe0[,bylist,with=F],phenDT[,log10.fold.induction:=log10(abs(yfp.mean.gal/yfp.mean.glu))],by='clone.named'),
	c('growth.rate','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal','log10.fold.induction'),bylist=bylist)

mDT<-melt(summ0[!is.na(clone.named)],id=bylist)
mDT[,c('pheno','stat'):=colsplit(variable,'_',c('pheno','stat'))]
form<-castform(c(bylist,'pheno'),c('stat'))
cDT<-dcast(mDT,form)


mod1<-lm(mean~cluster_Simple*pheno,data=cDT)
summary(mod1)
mod2<-lm(mean~cluster_Broad*pheno,data=cDT)
summary(mod2)
mod3<-lm(mean~cluster*pheno,data=cDT)
summary(mod3)
BIC(mod1,mod2,mod3)
anova(mod1,mod2)
anova(mod1,mod3)
anova(mod3,mod2)


cDT[,lab:=applyPaste(data.frame(substr(cluster_Broad,1,4),': ',substr(cluster_Simple,1,4),': ',cluster_Narrow))]

sample(LETTERS,4)
pCNEO<-ggplot(cDT,aes(y=mean,x= lab))+geom_jitter(shape=21,width=0.1,size=0.5)+geom_violin(alpha=0.7)+
	theme_NewPub_noarrow+theme(axis.text.x=element_text(angle=30,hjust=1))+facet_wrap(~pheno,scales='free')+	
	xlab('cluster')+ylab('within-genotype mean phenotypic values')

pWHRL<-ggplot(cDT,aes(y=mean,x= cluster_Simple))+geom_jitter(shape=21,width=0.1,size=0.5)+geom_violin(alpha=0.7)+
	theme_NewPub_noarrow+theme(axis.text.x=element_text(angle=30,hjust=1))+facet_wrap(~pheno,scales='free')+	
	xlab('cluster')+ylab('within-genotype mean phenotypic values')

if(writeplot==T){
	w<-5.3;h<-4.5
	ggsave(paste0(figout,'1707xx-181213-pCNEO-PhenotypesAcrossClusters-MoreClusters.pdf'), pCNEO,width=w,height=h)	
	ggsave(paste0(figout,'1707xx-181213-pWHRL-PhenotypesAcrossClusters.pdf'), pWHRL,width=w,height=h)	

}

############
### plot archetypal expression densities as line plots: use all clusters
############

genos<-c('clone.named','aa','bb','cc','single_lv')
DTA<-copy(densDT[,c(genos,'mutant_gene',grepincols(densDT,'dens__')),with=F])
DTX<-merge(DTA,clusterframe0,by='clone.named')
DTX[,mutant_facet:=mutant_gene]
DTX[mutant_gene=='GAL3 GAL80 GAL4',mutant_facet:='Control']
DTX[single_lv==T&aa=='GAL4.delta',mutant_facet:='Control']
DTX[single_lv==T&bb=='GAL3.delta',mutant_facet:='Control']
DTX[single_lv==T&cc=='GAL80.delta',mutant_facet:='Control']
mut_facet_levels<-c('Control','GAL3','GAL80','GAL4','GAL3 GAL4','GAL3 GAL80','GAL80 GAL4')
DTX[,mutant_facet:=factor(mutant_facet,levels=mut_facet_levels)]

meltby<-c(genos,'cluster','mutant_facet','cluster_Broad','cluster_Narrow','cluster_Simple')
mDT1<-densMelt(DTX[,!'mutant_gene'],meltby=meltby, summby = c("cluster_Broad",'cluster_Narrow','mutant_facet', "sugar", "AU.FL"))

sample(letters,4)
len<-length(unique(mDT1$clone.named))
labelframe<-mDT1[,length(unique(clone.named)),by=c('cluster_Broad','cluster_Narrow', 'mutant_facet')]
labelframe[,frac:=(round(V1/len,4))]
labelframe[,x:=2]
labelframe[,y:=0.25]
labelframe[,label:=applyPaste(data.frame(V1,'\n',frac*100,'%'))]

pYTRK<-ExpDensLinePlot(mDT1,linethickness=0.2)+
	geom_text(data=labelframe,aes(x=x,y=y,label=label),size=2,col='grey50')+
	facet_grid(cluster_Broad+cluster_Narrow~ mutant_facet)
	
if(writeplot==T){
	w<-5;h<-4.2
	ggsave(paste0(figout,'1707xx-181119-pYTRK-ExpressionDensityClusters-LinePlot-FacetByMutatedGene.pdf'), pYTRK,width=w,height=h)	
}


############
### plot archetypal expression densities as line plots: use simple clustering
############

genos<-c('clone.named','aa','bb','cc','single_lv')
DTA<-copy(densDT[,c(genos,'mutant_gene',grepincols(densDT,'dens__')),with=F])
DTX<-merge(DTA,clusterframe0,by='clone.named')
DTX[,mutant_facet:=mutant_gene]
DTX[mutant_gene=='GAL3 GAL80 GAL4',mutant_facet:='Control']
DTX[single_lv==T&aa=='GAL4.delta',mutant_facet:='Control']
DTX[single_lv==T&bb=='GAL3.delta',mutant_facet:='Control']
DTX[single_lv==T&cc=='GAL80.delta',mutant_facet:='Control']
mut_facet_levels<-c('Control','GAL3','GAL80','GAL4','GAL3 GAL4','GAL3 GAL80','GAL80 GAL4')
DTX[,mutant_facet:=factor(mutant_facet,levels=mut_facet_levels)]

meltby<-c(genos,'cluster','mutant_facet','cluster_Broad','cluster_Narrow','cluster_Simple')
mDT1<-densMelt(DTX[,!'mutant_gene'],meltby=meltby, summby = c('cluster_Simple','mutant_facet', "sugar", "AU.FL"))

sample(letters,4)
len<-length(unique(mDT1$clone.named))
labelframe<-mDT1[,length(unique(clone.named)),by=c('cluster_Simple', 'mutant_facet')]
labelframe[,frac:=(round(V1/len,4))]
labelframe[,x:=2]
labelframe[,y:=0.25]
labelframe[,label:=applyPaste(data.frame(V1,'\n',frac*100,'%'))]

pKWUB<-ExpDensLinePlot(mDT1,linethickness=0.2)+
	geom_text(data=labelframe,aes(x=x,y=y,label=label),size=2,col='grey50')+
	facet_grid(cluster_Simple~ mutant_facet)

if(writeplot==T){
	w<-4.5;h<-3.2
	ggsave(paste0(figout,'1707xx-181119-pKWUB-ExpressionDensityClusters-LinePlot-FacetByMutatedGene.pdf'), pKWUB,width=w,height=h)	
}




############
### plot archetypal expression densities as heatmap tile plots: use simple clustering
############
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gg_color_hue(100)

genos<-c('clone.named','aa','bb','cc','single_lv')
DTA<-copy(densDT[,c(genos,'mutant_gene',grepincols(densDT,'dens__')),with=F])
DTX<-merge(DTA,clusterframe0,by='clone.named')
DTX[,mutant_facet:=mutant_gene]
DTX[mutant_gene=='GAL3 GAL80 GAL4',mutant_facet:='Control']
DTX[single_lv==T&aa=='GAL4.delta',mutant_facet:='Control']
DTX[single_lv==T&bb=='GAL3.delta',mutant_facet:='Control']
DTX[single_lv==T&cc=='GAL80.delta',mutant_facet:='Control']
mut_facet_levels<-c('Control','GAL3','GAL80','GAL4','GAL3 GAL4','GAL3 GAL80','GAL80 GAL4')
DTX[,mutant_facet:=factor(mutant_facet,levels=mut_facet_levels)]

clustlevsimp<-c('uninducible','weak expr','inducible','leaky','constitutive')
clustlevbroad<-c('uninducible','inducible','constitutive')
PI<-merge(clusterframe0,na.exclude(phenDT[,list(PI=mean(phenotypic.index),FI=mean(log(abs(yfp.mean.gal/yfp.mean.glu)))),by='clone.named']),by='clone.named')
PI[,clustfac1:=factor(cluster_Simple,levels=clustlevsimp)]
PI[,clustfac2:=factor(cluster_Broad,levels=clustlevbroad)]
PI[,c('PIM','FIM'):=list(mean(PI),mean(FI)),by='cluster']
cnames<-PI[order(clustfac1,clustfac2,FIM,PIM,FI,PI)]$clone.named

meltby<-c(genos,'cluster','mutant_facet','cluster_Broad','cluster_Narrow','cluster_Simple')
mDT1<-densMelt(DTX[DTX$clone.named%in%paste(cnames),!'mutant_gene'],meltby=meltby, summby = c('cluster_Simple','mutant_facet', "sugar", "AU.FL"))

clustlevs1<-c('')
mDT1[sugar=='glu',X:=factor(round(AU.FL,2))]
mDT1[sugar=='gal',X:=factor(round(AU.FL,2)+max(AU.FL)+abs(min(AU.FL)))]
mDT1[,Y:=factor(clone.named,levels=cnames)]
mDT1[,`density+0.001, log10`:= log10(round(density_mean,3)+0.001)]
collim<-round(range(mDT1$`density+0.001, log10`),2)+c(0,0.01)
uniqfun<-function(x,fac=1,shift=1)(unique(x)[order(unique(x))][seq(shift,length(unique(x)),by=fac)])
xlab<-round(c(uniqfun(mDT1$AU.FL,fac=4,shift=3),uniqfun(mDT1$AU.FL,shift=3,fac=4)),2)
xbreak<-uniqfun(mDT1$X,fac=4,shift=3)

mDT1[,dumm:=125]
singM<-mDT1[single_lv==T&mutant_facet!='Control',.N,by=c('Y','mutant_facet','dumm')]
cfx<-clusterframe0[clone.named%in%paste(cnames)]
clustPaste<-cfx[,applyPaste(data.frame(cluster_Broad,'\n',cluster_Narrow))]
clustM<-melt(cfx[,.(clone.named,cluster_Simple, cluster_Broad,clustPaste)],id='clone.named')
clustM[,dumm:=4+c(rep(128,nrow(cfx)),rep(131,nrow(cfx)),rep(134,nrow(cfx)))]
#clustM[,dumm:=4+c(rep(128,nrow(cfx)),rep(131,nrow(cfx)),rep(NA,nrow(cfx)))]
clustM[,Y:=factor(clone.named,levels=cnames)]
custCols<-c(factorDefault[1:3],gg_color_hue(length(unique(clustM$value))))

PI[,clone.named:=factor(clone.named,levels=cnames)]
dupdum<-!duplicated(PI[order(clustfac1,clustfac2,FIM,PIM,FI,PI)]$cluster_Simple)#for identifying breakpoints for plotting
dupdumm<-PI[order(clustfac1,clustfac2,FIM,PIM,FI,PI)][dupdum]
dupdumm[,y0:=clone.named]
dupdumm[,yend:=clone.named]
dupdumm[,x0:=0]
dupdumm[,xend:=130]

sample(letters,4)



pSTVX_PDF<-ggplot(mDT1,aes(X,Y))+geom_tile(aes(fill=`density+0.001, log10`))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),na.value='white',limits=collim)+
	# same as below with no points to avoid crashing illustrator
	theme_NewPub_noarrow+
	theme(axis.text.x=element_text(angle=90,vjust=0.4,hjust=1,size=8))+
	theme(axis.text.y=element_text(size=6))+
	geom_vline(xintercept=61,col='black',size=2)+
	geom_vline(xintercept=61,col='white',size=1.5)+
	scale_x_discrete(breaks=xbreak,labels=xlab)+
	theme(axis.title.y=element_blank(),axis.text.y=element_blank(),	axis.ticks.y=element_blank())+
	xlab('Expression in glucose               Exression in galactose   \n pseudo log10 A.U. GAL1p-YFP fluorescence units')+
	theme(axis.title=element_text(size=10,face="plain"))
pSTVX_PNG <- pSTVX_PDF+
	# same as above with points indicating class and single mutants which will crash illustrator
	geom_point(data=data.frame(x=126,y=1),aes(x,y),col='white')+
	geom_point(data=singM,aes(dumm-runif(length(dumm))*3,Y,colour=`mutant_facet`),size=1,alpha=0.5,shape=21,inherit.aes=F)+
	geom_point(data=clustM,aes(dumm-runif(length(dumm))*1,Y,colour=`value`),size=1,shape=15,alpha=0.5,inherit.aes=F)+
	geom_segment(data= dupdumm,aes(x=x0,xend=xend,y=y0,yend=yend),col='black',size=0.2)+
	scale_colour_manual(values=custCols)


g_legend <- function(a.gplot){ 
	# Created on 2018-05-31 by the reprex package (v0.2.0). https://reprex.tidyverse.org/
	tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
	legend <- tmp$grobs[[leg]] 
	return(legend)
} 

legend <- g_legend(pSTVX) 
grid.arrange(legend)

if(writeplot==T){
	w<-4.6;h<-3
	ggsave(paste0(figout,'1707xx-190109-pSTVX-ExpressionDensityTilePlot-AllCloneNamed.pdf'), pSTVX_PDF+theme(legend.position='none'),width=w,height=h)	
	w<-4.6;h<-5
	ggsave(paste0(figout,'1707xx-190109-pSTVX-ExpressionDensityTilePlot-AllCloneNamed.png'), pSTVX_PNG+theme(legend.position='none',panel.border=element_blank()),width=w,height=h)	
	w<-5.5;h<-7.5
	ggsave(paste0(figout,'1707xx-190109-pSTVX-ExpressionDensityTilePlot-AllCloneNamed-ONLY_LEGEND.pdf'), grid.arrange(g_legend(pSTVX_PNG)),width=w,height=h)	
}



############
### plot how single-mutant clustering is predictive of growth rate
############
clustermerge<-function(clusterset){
	DTA<-copy(phenDT)
	DTX<-merge(DTA,clusterframe0,by='clone.named')
	setnames(DTX,pheno,'pheno')
	summ1<-DTX[,summary.func(pheno),by=c(genos,clusterset)]
	
	phen4<-DTX[G4.single==T,mean(pheno,na.rm=T),by=c('aa',clusterset)]
	phen3<-DTX[G3.single==T,mean(pheno,na.rm=T),by=c('bb',clusterset)]
	phen80<-DTX[G80.single==T,mean(pheno,na.rm=T),by=c('cc',clusterset)]
	
	dt.list=list(phen4, phen3, phen80)
	boo1<-c('aa','bb','cc')
	lapply(1:length(dt.list),function(i){
		x<-dt.list[[i]]
		setnames(x,c('V1'),c(paste0(boo1[i],'_','pheno_mean')))
		setnames(x,clusterset,paste0(boo1[i],'_',clusterset,'_Cluster'))
	})
	DTm1<-merge_recurse3(qDT=DTX,dt.list=dt.list,by.list=boo1)
	clusterset2<-grepincols(DTm1,'Cluster')
	summ0<-DTm1[,summary.func(pheno,'mu'),by=clusterset2]
	
	
	# Brief interlude (don't run):
	if(writeplot==T&clusterset=='cluster_Simple'){
		# for the intermediate clustering, cluster_Simple
		# you can see here that the missing factors from the 32 possible combinations of 4 * 4 * 2 cluster combinations are that triple mutants
		# are not included. Specifically, "uninducible" GAL3 was never paired with the three classes of GAL80 and GAL4 that do not include "inducible". 
		# These three classes are 3 X 3 X 1 GAL4, GAL80 and GAL3 alleles = 9. 32 - 9 = 23 observed combinations
		tabPossible<-data.table(expand.grid(unique(summ0$aa_cluster_Simple_Cluster),unique(summ0$bb_cluster_Simple_Cluster),unique(summ0$cc_cluster_Simple_Cluster)))
		cols<-c('aa_cluster_Simple_Cluster','bb_cluster_Simple_Cluster','cc_cluster_Simple_Cluster')
		colnames(tabPossible)<-cols
		tabObserved<-merge(tabPossible,summ0,by=cols,all=T)
		# these tables include overlapping unique combinations, therefore appearing to have 32 unique combinations.
		setnames(tabObserved,cols,c('GAL4','GAL3','GAL80'))
		tabObserved[,arb:=1:.N]
		tabObserved[GAL4=='inducible']
		tabObserved[GAL3=='inducible']
		tabObserved[GAL80=='inducible']
	
	}
	
	DTm2<-merge(summ0,DTm1,by=clusterset2)
	DTm3<-merge(DTm2[,summary.func(pheno,'pheno'),by='clone.named'],DTm2,by='clone.named')
	
	DTm3[aa%in%c('GAL4-L868G','GAL4-L868K','GAL4-L868P','GAL4-L868E','GAL4-L868S')&cc%in%c('GAL80S-1'),harm_comb:=T]
	# ggplot(DTm3[order(harm_comb)],aes(jitter(mu_mean,factor=10),pheno))+geom_point(aes(col=harm_comb),shape=21,size=2)+geom_abline(col='grey70')
	DTm3[,varexplained(pheno,mu_mean),by=clusterset2]
	
	varexpTotal<-DTm3[,paste0('var. explained: ',round(varexplained(pheno,mu_mean),3))]
	varexp<-DTm3[single_lv==F,paste0('var. explained: ',round(varexplained(pheno,mu_mean),2)),by='mutant_gene']
	setnames(varexp,'V1','fac')
	DTm4<-merge(DTm3[!duplicated(clone.named),!'pheno'],varexp,by='mutant_gene')
	DTm4[,fac2:=applyPaste(data.table(mutant_gene,fac),'\n')]
	return(list(data.table(DTm4),summ0, varexpTotal, varexp))
}


pheno<-'growth.rate'
genos<-c('clone.named','aa','bb','cc','single_lv')
clusterset<-c('cluster_Simple')# var explained: 0.906
# clusterset<-c('cluster_Broad')# var explained: 0.692
# clusterset<-c('cluster_Narrow')# var explained: 0.214
# clusterset<-c('cluster_Narrow','cluster_Broad')# var explained: 0.91
simple<-clustermerge('cluster_Simple')
broad<-clustermerge('cluster_Broad')
complex<-clustermerge(c('cluster_Narrow','cluster_Broad'))

DTm4<-simple[[1]]
pKCDH<-ggplot(DTm4[single_lv==F&order(cluster_Simple)],aes(jitter(mu_mean,factor=10),pheno_mean))+
#	geom_point(data= DTm4[harm_comb==T],aes(mu_mean,pheno_mean), col='black',shape=21,size=2,inherit.aes=F)+
	geom_point(aes(col=cluster_Simple),shape=21,size=1,stroke=0.2)+geom_abline(col='grey70')+
	facet_grid(~fac2)+
	theme_NewPub+
	factorDefault2+
	theme(legend.position='bottom')+
	xlim(c(0,0.3))+ylim(c(0,0.3))+
	xlab(expression(paste('predicted growth rate ',mu,' ',hr^-1)))+ylab(expression(paste('observed growth rate ',mu,' ',hr^-1)))

sample(letters,4)
DTm4<-broad[[1]]
DTm4[,`Broad Cluster`:=cluster_Broad]
pXGYL<-ggplot(DTm4[single_lv==F&order(cluster_Simple)],aes(jitter(mu_mean,factor=10),pheno_mean))+
#	geom_point(data= DTm4[harm_comb==T],aes(mu_mean,pheno_mean), col='black',shape=21,size=2,inherit.aes=F)+
	geom_point(aes(col=`Broad Cluster`),shape=21,size=1,stroke=0.2)+geom_abline(col='grey70')+
	facet_grid(~fac2)+
	theme_NewPub+
	factorDefault2+
	theme(legend.position='bottom')+
	xlim(c(0,0.3))+ylim(c(0,0.3))+
	xlab(expression(paste('predicted growth rate ',mu,' ',hr^-1)))+ylab(expression(paste('observed growth rate ',mu,' ',hr^-1)))


DTm4<-complex[[1]]
DTm4[,`Narrow Cluster`:=applyPaste(data.frame(cluster_Broad,'\n',cluster_Narrow))]
pQVKS<-ggplot(DTm4[single_lv==F&order(cluster_Simple)],aes(jitter(mu_mean,factor=10),pheno_mean))+
#	geom_point(data= DTm4[harm_comb==T],aes(mu_mean,pheno_mean), col='black',shape=21,size=2,inherit.aes=F)+
	geom_point(aes(col= `Narrow Cluster`),shape=21,size=1,stroke=0.2)+geom_abline(col='grey70')+
	facet_grid(~fac2)+
	theme_NewPub+
	factorDefault2+
	theme(legend.position='bottom')+
	xlim(c(0,0.3))+ylim(c(0,0.3))+
	xlab(expression(paste('predicted growth rate ',mu,' ',hr^-1)))+ylab(expression(paste('observed growth rate ',mu,' ',hr^-1)))


if(writeplot==T){
	w<-4;h<-2.3
	ggsave(paste0(figout,'1707xx-181119-pKCDH-PredictabilityBySimpleClusterMeans.pdf'), pKCDH,width=w,height=h)	
	ggsave(paste0(figout,'1707xx-181119-pKCDH-PredictabilityBySimpleClusterMeans-Supp.pdf'), pKCDH+ggtitle(paste0('Intermediate classification\nnumber of parameters: ',nrow(simple[[2]]),'\nTotal ',simple[[3]]))+
		theme(plot.title=element_text(size=8)),width=w,height=h+0.5)
	ggsave(paste0(figout,'1707xx-181228-pXGYL-PredictabilityByBroadClusterMeans.pdf'), pXGYL+ggtitle(paste0('Broad classification\nnumber of parameters: ',nrow(broad[[2]]),'\nTotal ',broad[[3]]))+
		theme(plot.title=element_text(size=8)),width=w,height=h+0.5)	
	ggsave(paste0(figout,'1707xx-181228-pQVKS-PredictabilityByComplexClusterMeans.pdf'), pQVKS+ggtitle(paste0('Narrow classification\nnumber of parameters: ',nrow(complex[[2]]),'\nTotal ',complex[[3]]))+
		theme(plot.title=element_text(size=8)),width=w,height=h+0.7)	
}

############
### plot proportions of single mutant combinations that lead to given double-mutant phenotypes
############
# you have to run the code above for this to work
DTm4<-simple[[1]]
clusterset2<-grepincols(DTm4,'Cluster')
DTN<-merge(DTm4[single_lv!=T,.N,by=c(clusterset2,'cluster_Simple')],DTm4[,.N,by='cluster_Simple'],by='cluster_Simple')
DTN[,frac:=N.x/N.y]
DTN[,`AB combination`:=factor(cluster_Simple,levels=c('inducible','uninducible','constitutive','leaky','weak expr'))]
DTN[order(`AB combination`,frac),ord:=seq_len(nrow(DTN))]
DTN[order(`AB combination`,frac)]
DTN[,identity:=applyPaste(data.frame("G4:",aa_cluster_Simple_Cluster,'\nG3:',bb_cluster_Simple_Cluster,'\nG80:',cc_cluster_Simple_Cluster),' ')]
DTN$dum<-runif(nrow(DTN))*0.001
DTN[,fac:=rank(frac+dum),by='AB combination']
DTN[,fac2:=abs(fac-max(fac))+1,by='AB combination']
DTN[,`top ranking underlying\nsingle-mutant phenotypes`:=as.factor(fac2)]
DTN[fac2>4, `top ranking underlying\nsingle-mutant phenotypes`:=factor(paste0('>',4))]
DTN2<-DTN[order(`AB combination`,`top ranking underlying\nsingle-mutant phenotypes`,decreasing=T)]
DTN2[,y:=cumsum(frac),by=`AB combination`]
DTN2[frac>0.1,label:= identity]

sample(LETTERS,4)
pNBVH<-ggplot(DTN2,aes(x= `AB combination`, y= frac,fill=`top ranking underlying\nsingle-mutant phenotypes`))+
	geom_bar(stat='identity',alpha=0.3,width=0.4)+
#	scale_alpha_manual(values=(rep(0.1,5)))+
	theme_NewPub_noarrow+
	xlab('double mutant combination')+ylab('proportion')+
	theme(legend.position='bottom')+
	geom_label(aes(x=`AB combination`,y=y-0.07,label=label),size=1.7,lwd=NA,fill=NA)


DTN2[,identity2:=applyPaste(data.frame(ord,'\n',"G4:",aa_cluster_Simple_Cluster,'\nG3:',bb_cluster_Simple_Cluster,'\nG80:',cc_cluster_Simple_Cluster),' ')]
DTN2[order(ord),identity2:=factor(identity2,levels=DTN2[order(ord),identity2])]
DTN2[frac>0.1,label2:= identity2]


pIKDL<-ggplot(DTN2[order(ord)],aes(x= `AB combination`, y= frac,fill= `top ranking underlying\nsingle-mutant phenotypes`))+
	geom_bar(stat='identity',width=0.4,col='grey60')+
#	scale_alpha_manual(values=(rep(0.1,5)))+
	scale_fill_manual(values=c('grey80','grey80','grey80','grey80','black'))+
	theme_NewPub_noarrow+
	xlab('double mutant combination')+ylab('proportion')+
	theme(legend.position='bottom')+
	geom_label(aes(x=`AB combination`,y=y-0.07,label= label),size=1.7,lwd=NA,fill=NA)
	

if(writeplot==T){
	w<-3.8;h<-3.9
	ggsave(paste0(figout,'1707xx-181203-pNBVH-FractionInABFromA,B.pdf'), pNBVH,width=w,height=h)
	ggsave(paste0(figout,'1707xx-190208-pIKDL-FractionInABFromA,B-GreyScale.pdf'), pIKDL,width=w,height=h)
}



###################################################
### prediction of growth rate based on multiplicative model
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
predVar<-merge(plotDT[,.(x_mean,pairing,clone.named)],phenDT[,.(growth.rate,clone.named)],by='clone.named')

# overall variance explained by multiplicative model
predVar[,varexplained(growth.rate,x_mean)]
varExpx<-round(predVar[,varexplained(growth.rate,x_mean)],3)
# table of within-pairing variance explained by the model
predVar2<-predVar[,round(varexplained(growth.rate,x_mean),3),by='pairing']
setnames(predVar2,'V1','varexp')

# merge it together for plotting
plotDT2<-merge(plotDT,predVar2,by='pairing')
plotDT2[,lab:=applyPaste(data.table(pairing,'\n',varexp))]
sample(LETTERS,4)
lims<-c(-0.03,0.27)
plotDT2[x_upper>=max(lims),x_upper:=max(lims)]
plotDT2[x_lower<=min(lims),x_lower:=min(lims)]
plotDT2[y_upper>=max(lims),y_upper:=max(lims)]
plotDT2[y_lower<=min(lims),y_lower:=min(lims)]
plotDT3<-merge(plotDT2,clusterframe0,by='clone.named')

pLPRC <-ggplot(plotDT3[order(sig)],aes(x_mean,y_mean,xmax=x_upper,xmin=x_lower,ymax=y_upper,ymin=y_lower))+theme_NewPub+
	geom_errorbar(col='grey50',alpha=0.2)+geom_errorbarh(col='grey50',alpha=0.2)+geom_point(col='black',alpha=0.3)+
	geom_abline(col='grey70')+xlim(lims)+ylim(lims)+
	xlab(expression(paste('predicted growth rate ',mu,' ',hr^-1)))+ylab(expression(paste('observed growth rate ',mu,' ',hr^-1)))+
	facet_wrap(lab~.,ncol=1)

plotDT3[,`Double mutant 'Intermediate' cluster`:=cluster_Simple]
pPGOT <-ggplot(plotDT3[order(sig)],aes(x_mean,y_mean,xmax=x_upper,xmin=x_lower,ymax=y_upper,ymin=y_lower))+theme_NewPub+
	geom_errorbar(col='grey50',alpha=0.2)+geom_errorbarh(col='grey50',alpha=0.2)+geom_point(aes(col=`Double mutant 'Intermediate' cluster`),shape=21,size=1,stroke=0.2)+
	geom_abline(col='grey70')+
	factorDefault2+
	theme(legend.position='none')+
	xlab(expression(paste('predicted growth rate ',mu,' ',hr^-1)))+ylab(expression(paste('observed growth rate ',mu,' ',hr^-1)))+
	xlim(c(0,0.3))+ylim(c(0,0.3))+
	facet_grid(~lab)+ 
	ggtitle(paste0('Geometric model predictions based on single-mutant phenotypes\nTotal var. explained: ',varExpx))+
	theme(plot.title=element_text(size=8),legend.title=element_blank())


if(writeplot==T){
	w<-2.2;h<-4.5
	ggsave(paste0(figout,'1707xx-181213-pLPRC-multiplicative_model_predictedVsObserved.png'), pLPRC,width=w,height=h)	
	w<-4;h<-2.3
	ggsave(paste0(figout,'1707xx-181228-pPGOT-multiplicative_model_predictedVsObserved-ColByCluster.png'), pPGOT,width=w,height=h)	

}




###################################################
### plotting for presentations
###################################################

ylimits<-c(0,0.35)
genos<-c('aa','bb','cc','pairing')
denscolconv<-function(DT, vnull =T){
	DT[,variable:=gsub('dens.','',variable)]
	DT[,c('AU.FL','sugar'):=colsplit(variable,'.g',c('AU.FL','sugar'))]
	DT[,sugar:=paste0('g',sugar)]
	if(vnull ==T)DT[,variable:=NULL]
}

DTx<-copy(densDT[!is.na(clone.named)])
DTX<-merge(DTx,clusterframe0,by='clone.named')

DTX[,pairing:=mutant_gene]
WTFrame<-DTX[WT==T]
a<-copy(WTFrame);b<-copy(WTFrame);cx<-copy(WTFrame);aa<-copy(WTFrame);bb<-copy(WTFrame);cxc<-copy(WTFrame)
a[,pairing:='GAL80 GAL4'];b[,pairing:='GAL3 GAL80'];cx[,pairing:='GAL80 GAL4'];aa[,pairing:='GAL3 GAL4'];bb[,pairing:='GAL3 GAL4'];cxc[,pairing:='GAL3 GAL80']
DTX1<-rbind(DTX[WT!=T],a,b,cx,aa,bb,cxc)

singleframe<-DTX[single_lv==T]
a<-copy(singleframe[WT!=T]);b<-copy(singleframe[WT!=T]);cx<-copy(singleframe[WT!=T]);aa<-copy(singleframe[WT!=T]);bb<-copy(singleframe[WT!=T]);cxc<-copy(singleframe[WT!=T])
a[,pairing:='GAL80 GAL4'];b[,pairing:='GAL3 GAL80'];cx[,pairing:='GAL80 GAL4'];aa[,pairing:='GAL3 GAL4'];bb[,pairing:='GAL3 GAL4'];cxc[,pairing:='GAL3 GAL80']
DTX2<-rbind(DTX1[single_lv!=T],a,b,cx,aa,bb,cxc)

singleframe2<-merge(phenDT[single_lv==T],clusterframe0,by='clone.named')
singleframe2[grepl('GAL4',mutant_gene),gene:='GAL4']
singleframe2[grepl('GAL3',mutant_gene),gene:='GAL3']
singleframe2[grepl('GAL80',mutant_gene),gene:='GAL80']
summ<-singleframe2[,summary.func(phenotypic.index),by=c('mutants','cluster_Simple','gene')]
summ[,var1:=mean(mean),by=cluster_Simple]
G3<-c(c('GAL3.WT','GAL3.delta'),notin(gsub('GAL3.WT GAL80.WT GAL4.WT','GAL3.WT',summ[gene=='GAL3'|mutants=='GAL3.WT GAL80.WT GAL4.WT'][order(var1,mean)]$mutants),c('GAL3.WT','GAL3.delta')))
G80<-c(c('GAL80.WT','GAL80.delta'),notin(gsub('GAL3.WT GAL80.WT GAL4.WT','GAL80.WT',summ[gene=='GAL80'|mutants=='GAL3.WT GAL80.WT GAL4.WT'][order(var1,mean)]$mutants),c('GAL80.WT','GAL80.delta')))
G4<-c(c('GAL4.WT','GAL4.delta'),notin(gsub('GAL3.WT GAL80.WT GAL4.WT','GAL4.WT',summ[gene=='GAL4'|mutants=='GAL3.WT GAL80.WT GAL4.WT'][order(var1,mean)]$mutants),c('GAL4.WT','GAL4.delta')))


DTX2[,aa:=factor(aa,levels=G4)]
DTX2[,bb:=factor(bb,levels=G3)]
DTX2[,cc:=factor(cc,levels=G80)]


mDT0<-melt(DTX2[,c(dens.cols,genos,'cluster_Simple'),with=F],id=c(genos,'cluster_Simple'))#;mDT<-mDT[order(sugar,AU.FL,decreasing=T)]

mDT<-mDT0[,list(density_mean=mean(value,na.rm=T),density_sd=sd(value,na.rm=T)),by=c(genos,'cluster_Simple','variable')]
denscolconv(mDT)
mDT[,AU.FL:=as.numeric(gsub('_','',AU.FL))]
mDT[density_mean>=max(ylimits),density_mean:=max(ylimits)]

G4<-c('GAL4.WT', 'GAL4-L868*','GAL4-L868C','GAL4-L868E','GAL4-L868F','GAL4-L868G','GAL4-L868K','GAL4-L868P','GAL4-L868S','GAL4.delta','GAL4.30') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80.delta','GAL80S-1','GAL80S-2','GAL80S-0','GAL80S-V352E')
G3<-c('GAL3.WT')
mDT1<-mDT[aa%in%G4&cc%in%G80&bb=='GAL3.WT'|(cc=='GAL80.WT'&aa=='GAL4.WT')&pairing=='GAL80 GAL4']
mDT1[,c('x','y'):=list(2,0.2)]
sample(LETTERS,4)
pKUQF<-ExpDensLinePlot(mDT1)+facet_grid(cc~aa)+geom_text(data=mDT1[,mean(density_mean),by=c('cluster_Simple','aa','cc','x','y')],aes(x=x,y=y,label=cluster_Simple),size=2,col='grey70')

if(writeplot==T){
	w<-6.8;h<-4.5
	ggsave(paste0(figout,'1707xx-181211-pKUQF-SubsetForPresentations.png'), pKUQF,width=w,height=h)
	
}


G4<-c('GAL4.WT','GAL4-L868G') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80S-1')
G3<-c('GAL3.WT')
mDT1<-mDT[aa%in%G4&cc%in%G80&bb=='GAL3.WT']
mDT1[,c('x','y'):=list(2,0.2)]
sample(LETTERS,4)
pAGCO<-ExpDensLinePlot(mDT1)+facet_grid(cc~aa)+geom_text(data=mDT1[,mean(density_mean),by=c('cluster_Simple','aa','cc','x','y')],aes(x=x,y=y,label=cluster_Simple),size=2,col='grey70')

if(writeplot==T){
	w<-2.35;h<-1.88
	ggsave(paste0(figout,'1707xx-181211-pAGCO-SubsetForPresentations.png'), pAGCO,width=w,height=h)
	
}


G4<-c('GAL4.WT', 'GAL4-L868E','GAL4-L868G','GAL4-L868P','GAL4.delta') # G4<-'GAL4-L868P'
G80<-c('GAL80.WT','GAL80.delta','GAL80S-1','GAL80.07')
G3<-c('GAL3.WT')
mDT1<-mDT[aa%in%G4&cc%in%G80&bb=='GAL3.WT'|(cc=='GAL80.WT'&aa=='GAL4.WT')&pairing=='GAL80 GAL4']
mDT1[,c('x','y'):=list(2,0.2)]
sample(LETTERS,4)
pPSBD<-ExpDensLinePlot(mDT1)+facet_grid(cc~aa)
ExpDensLinePlot(mDT1)+facet_grid(cc~aa)+geom_text(data=mDT1[,mean(density_mean),by=c('cluster_Simple','aa','cc','x','y')],aes(x=x,y=y,label=cluster_Simple),size=2,col='grey70')
if(writeplot==T){
	w<-3.75;h<-2.8
	ggsave(paste0(figout,'1707xx-181211-pPSBD-SubsetForPresentations.png'), pPSBD,width=w,height=h)
	
}


mDT<-mDT0[,list(density_mean=mean(value,na.rm=T),density_sd=sd(value,na.rm=T)),by=c('cluster_Simple','variable')]
denscolconv(mDT)
mDT[,AU.FL:=as.numeric(gsub('_','',AU.FL))]
mDT[density_mean>=max(ylimits),density_mean:=max(ylimits)]
mDT1<-mDT
mDT1[,c('x','y'):=list(2,0.2)]
sample(LETTERS,4)
pWJGO<-ExpDensLinePlot(mDT1)+facet_grid(~cluster_Simple)+geom_text(data=mDT1[,mean(density_mean),by=c('cluster_Simple','aa','cc','x','y')],aes(x=x,y=y,label=cluster_Simple),size=2,col='grey70')

if(writeplot==T){
	w<-3.9;h<-1.45
	ggsave(paste0(figout,'1707xx-181211-pWJGO-FiveClustersForPresentations.png'), pWJGO,width=w,height=h)
	
}


