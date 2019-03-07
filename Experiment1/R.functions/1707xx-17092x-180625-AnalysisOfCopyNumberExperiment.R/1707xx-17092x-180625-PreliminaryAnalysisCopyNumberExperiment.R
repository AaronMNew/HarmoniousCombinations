
# Scripts for dealing with flow cytometry data.
 # biocLite("flowCore")
# biocLite("flowViz")
# biocLite("ggcyto")
# install.packages("diptest")
source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'
source(source.code.path.file)
load.all() # this function loads all the standard packages you use in this script

head.dir<-"/Users/anew/Google Drive/002_lehnerlabData/1707xx-BigMixNMatch"
date<-'1707xx' # beginning date clones were measured
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'/figure.output/'
datout<-'/output.data/'
setwd(head.dir)

br<-seq(-0.5,2.7,length=64)

###################################################
### custom functions for analysis
###################################################
source(source.code.path.file)

t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 
{
    if( equal.variance==FALSE ) 
    {
        se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        # welch-satterthwaite df
        df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    } else
    {
        # pooled standard deviation, scaled by the sample sizes
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
        df <- n1+n2-2
    }      
    t <- (m1-m2-m0)/se 
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
    names(dat) <- c("diff.of.means", "std.err", "t", "p.value")
    return(dat) 
}

t.test3 <- function(m1,m2,se,df,m0=0){
# m1, m2: the sample means
# s1, s2: the sample standard deviations
# se: the standard error pre-calculated
# df: the degrees of freedom, pre-calculated
    t <- (m1-m2-m0)/se 
    dat <- list(m1-m2, se, t, pt(-abs(t),df))    
    names(dat) <- list("diff.of.means", "std.err", "t", "p.value")
    return(dat) 
}
	
merge_recurse2<-function(qDT=NULL,dt.list,by.list,names.list=NULL){
	x<-qDT
	if(is.null(qDT))x<-dt.list[[1]]
	if(length(dt.list)!=(length(by.list)))stop('list lengths look funny')
	i<-2
	while(i <= length(dt.list)){
		x<-merge(x,dt.list[[i]],by=by.list[[i]],all=T)
		i<-i+1
	}
	x
}

grepincols<-function(DT,variables){
	x<-colnames(DT)
	unlist(lapply(variables,function(xv)x[grep(xv,x)]))
}

summary.func.all<-function(DAT,variablelist,bylist){
	summarylist<-lapply(1:length(variablelist),function(i){
		DTx<-DAT[,c(variablelist[i],bylist),with=F]
		setnames(DTx,variablelist[i],'V1')
		DTx[,summary.func(V1,paste0(variablelist[i])),by=bylist]
	})
	bylist2<-lapply(1:length(variablelist),function(x)return(bylist))
	
	return(merge_recurse2(dt.list=summarylist,by.list=bylist2))
}

splitByVarlist<-function(DT,variablelist,bylist){
	out<-lapply(1:length(variablelist),function(i){
		DT[,grepincols(DT,c(bylist,variablelist[i])),with=F]
	})
	names(out)<-variablelist
	out
}

TukeyInteractionFlip<-function(df){
	# df must have A and B columns
	DTemp1<-na.exclude(data.table(df,data.table(colsplit(rownames(df),'-',c('A','B')))[,lapply(.SD,as.character)]))
	uniqfac<-unique(c(DTemp1 $B, DTemp1 $A))
	temp<-rbind(t(combn(uniqfac,2)),t(combn(uniqfac,2))[,c(2,1)])
	fullDT<-data.table(temp)[,lapply(.SD, as.character)]
	colnames(fullDT)<-c('A','B')
	DTemp1a<-copy(DTemp1)
	DTemp2<-copy(DTemp1)
	
	DTemp2[,c('A','B'):=list(B,A)]
	DTemp3<-rbind(DTemp1a,DTemp2)
	test<-merge(DTemp3,fullDT,by=c('B','A'))
	uniqtab1<-data.table(A=uniqfac,A.arb=(1:length(uniqfac)))
	uniqtab2<-data.table(B=uniqfac,B.arb=(1:length(uniqfac)))
	DTemp0<-merge(uniqtab2,merge(test,uniqtab1,by='A'),by='B')
	DTemp0[,temp:=apply(data.frame(A,B),1,function(x){
		a<-max(x)
		b<-min(x)
		paste0(a,b)
	})]
	DTemp4<-DTemp0[!duplicated(temp)]
	DTemp4[,temp:=NULL]
	
	asdf<-merge(DTemp4[,c('A','B','diff')],DTemp1[,c('A','B','diff')],by=c('A','B'),all=T)
	asdf[is.na(diff.y),flipflag:=-1]
	asdf[!is.na(diff.y),flipflag:=1]
	asdf[,diff.flipped:=diff.x*flipflag]
	DTemp<-merge(asdf[,c('A','B','diff.flipped')],DTemp4,by=c('A','B'))
	
	return(DTemp)
	}

varexplained<-function(y1,yexp){
	1-sum((y1-yexp)^2,na.rm=T)/sum(((y1-mean(y1,na.rm=T))^2),na.rm=T)
}


###################################################
### read in data, merge with layout and do preliminary calculations of growth rate etc.
###################################################
# FCS files are processed in 1707xx-17092x-FullAnalysis-PrettyGood
output<-fread(paste(date,'-170927-170928-reanalyzed.tab',sep=''),header=T,sep='\t')
layout<-na.exclude(fread('1707xx-17092x-LayoutByHand.txt',header=T,sep='\t'))

layout2<-fread('1707xx-layout.txt',header=T,sep='\t')

merged<-cbind(layout2[order(samp_id)][tx.date== '1709x1'],layout,output)
merged2<-merged[,!duplicated(colnames(merged)),with=F]
strainbackgroundSummaries<-fread('1707xx-17092x-StrainBackroundSummaries.txt',header=T,sep='\t')
merged2 <-merge(merged2,strainbackgroundSummaries,by='genotype.strain')


nrow(merged2)
merged2$genotype<-paste(merged2$GAL3.allele,merged2$GAL80.allele,merged2$GAL4.allele)
merged2$replicate<-paste(merged2$genotype,merged2$measurement.date)
merged2$clone.named<-paste(merged2$GAL3.clone.named,merged2$GAL80.clone.named,merged2$GAL4.clone.named)

old.cols<-colnames(merged2)[c(grep("\\]",colnames(merged2)))]
dens.cols<-paste('dens',round(br[1:(length(br)-4)],2),sep='.')
setnames(merged2,old= old.cols,new= dens.cols)
dens.norm<-data.table(t(apply(merged2[, dens.cols,with=FALSE],1,row.normalization)))
colnames(dens.norm)<-dens.cols
merged2[,dens.cols]<-dens.norm
# for cell density per ml, cells.per.sec * 1/µl per second sampling rate * 1000 µl / ml

merged2$cell.per.ml<-merged2$cells.per.sec*(1/as.numeric(merged2$sample_flow_rate))*1000
DT.glu<-merged2[glu==1]
DT.gal<-merged2[gal==1]


# to couple direct measurements with each other (glucose > galactose) i need to create an id variable that will allow joining glu with gal measurements
DT.glu$next.day<-as.numeric(DT.glu$measurement.date)+1
DT.gal$next.day<-NA
DT.glu$id<-paste(DT.glu$source.plate,DT.glu$rc,DT.glu$next.day,DT.glu$biol.rep)
DT.gal$id<-gsub('\n','',paste(DT.gal$source.plate,DT.gal$rc,DT.gal$measurement.date,DT.gal$biol.rep))



# DT.gal[!id %in% DT.glu$id]
# DT.glu[!id %in% DT.gal$id]
# DT.glu[grep('170928',id)]
DT0<-merged2[tx.date=='1709x1']
clone.nameds<-unique(DT0$clone.named)
dev.new()
par(mfrow=c(2,4))
ylimit<-c(0,55000)
for(i in 1:length(clone.nameds)){
	test<-grepcols(DT0[genotype.strain!='KV447'&clone.named==clone.nameds[i]&gal=='0'],c('yfp.mean','copies'))[,1:4]	
	boxplot(as.numeric(yfp.mean)~.,data=test,main=clone.nameds[i],horizontal=T,las=1,ylim=ylimit)
}

dev.new()
par(mfrow=c(2,4))
ylimit<-c(0,55000)
for(i in 1:length(clone.nameds)){
	test<-grepcols(DT0[genotype.strain!='KV447'&clone.named==clone.nameds[i]&gal=='1'],c('yfp.mean','copies'))[,1:4]	
	boxplot(as.numeric(yfp.mean)~.,data=test,main=clone.nameds[i],horizontal=T,las=1,ylim=ylimit)
}

# i re-merge them here and then calculate generations, growth rate and phenotypic index

#variables definition
dil.factor<-9/150
hours <-12

setnames(DT.glu,colnames(DT.glu),paste0(colnames(DT.glu),'.glu'))
setnames(DT.gal,colnames(DT.gal),paste0(colnames(DT.gal),'.gal'))
setnames(DT.glu,'id.glu','id')
setnames(DT.gal,'id.gal','id')
DT.merge<-merge(DT.glu,DT.gal,by='id',all=TRUE)
DT.merge[is.na(DT.merge $samp_id.glu)]
DT.merge[grep('170928',id)]

DT.merge$generations<-log((DT.merge$cell.per.ml.gal/((DT.merge$cell.per.ml.glu)*dil.factor)),base=2)
DT.merge$growth.rate<-log((DT.merge$cell.per.ml.gal/((DT.merge$cell.per.ml.glu)*dil.factor)))/hours
DT.merge$phenotypic.index<-DT.merge$fracon.glu+DT.merge$fracon.gal

plasmid.genotype.summary<-as.vector(sapply(DT.merge$clone.named.glu,function(x){
	xa<-strsplit(x,' ')[[1]]
	paste(as.numeric(grepl('WT',xa)),collapse=' ')
}))
DT.merge$plasmid.genotype.summary<-plasmid.genotype.summary
grepcols(DT.merge,'cell.per')
summary.func(DT.merge$cell.per.ml.glu)
hist(log10(DT.merge$cell.per.ml.gal),50)
DT.merge.copynumber<-DT.merge

# there is an obvious contamination i discovered below during analysis that was not flagged during the recording process. 
# remove this from the data.table
DT.merge[,gr.temp:=paste0(growth.rate)]
a<-DT.merge[gr.temp=='0.167383942360909']$samp_id.glu
DT.merge[samp_id.glu==a]<-NA
DT.merge$gr.temp<-NULL

qDTx<-DT.merge[genotype.strain.glu!='KV447']
gencopies<-colnames(qDTx)[grep(glob2rx('copies*genome.glu'), colnames(qDTx))]
plascopies<-colnames(qDTx)[grep(glob2rx('copies*plasmid.glu'), colnames(qDTx))]
setnames(qDTx, gencopies,gsub('.glu','', gencopies))
setnames(qDTx, plascopies,gsub('.glu','', plascopies))
gencopies1<-gsub('.glu','', gencopies)
plascopies1<-gsub('.glu','', plascopies)

qDTx[,plascopies:=apply(qDTx[,plascopies1,with=F],1,function(x)paste(x[1],x[2],x[3],collapse=' '))]
qDTx[,gencopies:=apply(qDTx[,gencopies1,with=F],1,function(x)paste(x[1],x[2],x[3],collapse=' '))]
qDTx $plasgeno<-apply(qDTx[,c(plascopies1,gencopies1),with=F],1,function(x){paste(x,collapse=' ')})
qDTx $genoplas<-apply(qDTx[,c(gencopies1,plascopies1),with=F],1,function(x){paste(x,collapse=' ')})
qDTx $binarycount<-apply(qDTx[,c(plascopies1,gencopies1),with=F],1,function(x){x<-as.numeric(x);paste(as.numeric(as.logical(x[1:3]+x[4:6])),collapse=' ')})
qDTx $totalcopies<-apply(qDTx[,c(plascopies1,gencopies1),with=F],1,function(x){x<-as.numeric(x);paste(x[1:3]+x[4:6],collapse=' ')})
genobin<-qDTx[,c(grep(glob2rx('copies*genome'),colnames(qDTx)),grep(glob2rx('copies*plasmid'),colnames(qDTx))),with=F]
copies.GAL3<-apply(grepcols(genobin,'GAL3'),1,function(x){sum(as.numeric(x))})
copies.GAL80<-apply(grepcols(genobin,'GAL80'),1,function(x){sum(as.numeric(x))})
copies.GAL4<-apply(grepcols(genobin,'GAL4'),1,function(x){sum(as.numeric(x))})
DTcopies<-data.table(apply(data.frame(copies.GAL3,copies.GAL80,copies.GAL4),2,as.factor))
DTbinary<-data.table(apply(data.frame(copies.GAL3,copies.GAL80,copies.GAL4),2,function(x){
	as.factor(as.numeric(x>0))
	}))
setnames(DTbinary,colnames(DTbinary),gsub('copies','binary',colnames(DTbinary)))
qDT01<-data.table(qDTx,DTcopies,DTbinary)
qdf<-data.table(source.plate.glu =unique(qDT01 $source.plate.glu),tx.rep=letters[1:length(unique(qDT01$source.plate.glu))])
qDT<-merge(qdf,qDT01,by='source.plate.glu')


backg<-mean(qDT[copies.GAL4==0]$growth.rate)
backgsd<-sd(qDT[copies.GAL4==0]$growth.rate)/sqrt(length(qDT[copies.GAL4==0]$growth.rate))



qDT[growth.rate>0.15&yfp.mean.gal<20000]<-NA
qDT[, genome.genotype.summary:=background.genotype.summary.2.gal]
qplot(qDT $growth.rate, qDT $yfp.mean.gal)
quickDT<-data.table(qDT[,c('growth.rate','yfp.mean.gal','samp_id.gal','samp_id.glu'),with=F])
mod1<-lm(growth.rate~plasgeno,qDT)
summary(mod1)
summary(aov(mod1))
credible.growth.rate.background<-qDT[copies.GAL4==0,summary.func(growth.rate)][,mean+sd*1.96]

qplot(mod1$fitted.values,qDT$growth.rate)+geom_abline(intercept=0,slope=1)+geom_vline(xintercept=credible.growth.rate.background,col='grey70')+geom_hline(yintercept=credible.growth.rate.background,col='grey70')
accuracy(mod1)




qDT1<-data.table(grepcols(qDT,'.clone.named.glu'),qDT[,'growth.rate',with=F])
mod1<-lm(growth.rate~.,qDT1)
summary(mod1)

qDT2<-qDT[,c('background.genotype.summary.2.gal','plasmid.genotype.summary'),with=F]
qtrans<-data.table(t(data.frame(apply(qDT2,1,function(x){
	c(strsplit(x[1],' ')[[1]],strsplit(x[2],' ')[[1]])
}))))
colnames(qtrans)<-c(paste(c('GAL3','GAL80','GAL4'),c('.genome'),sep=''),paste(c('GAL3','GAL80','GAL4'),c('.plasmid'),sep=''))
qDT2a<-data.table(qDT[,'growth.rate',with=F],apply(data.frame(qtrans),2,as.factor))
mod2<-lm(growth.rate~.^3,qDT2a)
summary(mod2)
TukeyHSD(aov(mod2))
dens.cols.new<-colnames(DT.merge)[grep('dens',colnames(DT.merge))]


###################################################
### technical variation of experiment
###################################################
DT0<-copy(qDT)
DT0[,dumm:=1]
grepincols(DT0,'tx.rep')
phenotypes<-c('dumm','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal','growth.rate','phenotypic.index')
batchfac1<-c("biol.rep.glu",'plasgeno')
batchSumm<-summary.func.all(DT0,phenotypes,batchfac1)
b<-batchSumm[,!grepincols(batchSumm,'dumm'),with=F]
a<-melt(b,id=batchfac1)
a[,c('phenotype','stat'):=colsplit(variable,'_',c('phenotype','stat'))][,variable:=NULL]
biolRep<-dcast(a,formula(plasgeno+phenotype~stat+biol.rep.glu),value.var='value')
pRep<-ggplot(biolRep,aes(mean_1,mean_2,xmin=lower_1,xmax=upper_1,ymin=lower_2,ymax=upper_2))+geom_errorbarh()+geom_errorbar()+geom_point()+facet_wrap(~phenotype,scales='free')+geom_abline()

batchfac2<-c('plasgeno','biol.rep.glu', 'tx.rep','samp_id.glu')
phenotypes<-c('growth.rate','phenotypic.index','yfp.mean.glu','yfp.mean.gal','fracon.glu','fracon.gal')
batchSumm<-DT0[,lapply(.SD,mean),.SDcols=phenotypes,by='plasgeno']
DTm<-merge(batchSumm,DT0[,c(batchfac2,phenotypes),with=F],by='plasgeno',suffixes=c('_x','_y'))
mDT<-melt(DTm,id=batchfac2)
mDT[,c('phenotype','var'):=colsplit(variable,'_',c('phenotype','var'))][, variable:=NULL]
cDT<-dcast(mDT,formula(plasgeno+phenotype+ tx.rep +biol.rep.glu+samp_id.glu~var),value.var='value')
cDT[,varexp:=round(varexplained(x,y),3),by='phenotype']
cDT[,varexplained:=paste0('variance explained: ',varexp)]
cDT[,phenotype:=factor(phenotype,levels=phenotypes)]
pLZF_Rep<-ggplot(cDT,aes(x,y))+geom_point(aes(col= tx.rep),shape=21)+geom_abline()+facet_wrap(~phenotype+varexplained,scales='free')+xlab('within-genotype means')+ylab('independent observations')+
	theme(axis.text.x = element_text(size = 10, angle=30,hjust=1))

w<-10;h<-w*0.7
ggsave(paste0(head.dir,figout,'1707xx-17092x-160622-pLZF_Rep-technical_reprodcibility.png'), pLZF_Rep,width=w,height=h,limitsize=F)



# a few assemblies gave weird measurements for fracon.gal; have a quick look
mod1<-lm(y~phenotype*(plasgeno+ tx.rep +biol.rep.glu),cDT)
summary(mod1)
asdf<-DT0[,summary.func(fracon.gal),by=plasgeno]
qplot(asdf$mean,asdf$CoV)
mod1<-lm(log(CoV)~log(mean),asdf)
asdf[,outlier:=exp(mod1$residuals)>15]
asdf[outlier==T]
ggplot(asdf,aes((mean),(CoV),col=outlier))+geom_point()
hist(mod1$residuals)
qplot(asdf$mean,asdf$CoV)


###################################################
### basic phenotypic overview plotting
###################################################

##################
# single-gene effects expression density
##################

# expression density plots
coregenotypes<-c('0 0 0 1 1 1', # WT
	'0 0 0 1 1 0', # GAL4.delta
	'0 0 0 0 1 1', # GAL3.delta
	'0 0 0 1 0 1' # GAL80.delta
	)
DTemp<-qDT[,c('plasgeno',dens.cols.new),with=F][plasgeno%in%coregenotypes]
DTemp$plasgeno <-factor(DTemp$plasgeno,levels=coregenotypes)
a<-melt(DTemp);setnames(a,'value','density')
a$variable<-gsub('dens.','',a$variable)
b<-cbind(a,colsplit(a$variable,'\\.g',c('AU.FL.VAL','sugar')))
b$sugar<-paste0('g',b$sugar)
cx<-b[,mean(density),by=c('plasgeno','sugar','AU.FL.VAL')];setnames(cx,'V1','density')
mdf<-merge(b,cx,by=c('plasgeno','sugar','AU.FL.VAL'))
mdf$plasgeno <-factor(mdf$plasgeno,levels=coregenotypes)
paste(sample(LETTERS,3),collapse='')
xlabel2<-'pseudo-log10 fluorescence A.U.'
ylabel2<-'density'

p00IZF<-ggplot(mdf[!is.na(plasgeno)],aes(x=AU.FL.VAL,y=density.y))+
	geom_line(aes(colour=sugar),size=1)+
	scale_y_continuous(breaks=c(0,0.1,0.2,0.3))+
	facet_wrap(~plasgeno,ncol=4)+
	scale_colour_manual(values=c('cornflowerblue','indianred'))+coord_cartesian(ylim=c(0,0.4))+xlab(xlabel2)+ylab(ylabel2)+theme(legend.position = 'none')

##################
# change in cell density of few genotypes with simple dot plots at the beginning and end
##################

DTemp<-qDT[!is.na(plasgeno),c('plasgeno', 'plascopies','gencopies','cell.per.ml.glu','cell.per.ml.gal'),with=F][plasgeno%in%coregenotypes]
DTemp[,denisty.gal:=cell.per.ml.gal*1/dil.factor]
summgal<-DTemp[,summary.func(denisty.gal,'gal'),by=plasgeno]
#summgal$`hours elapsed`<-12
summglu<-DTemp[,summary.func(cell.per.ml.glu,'glu'),by=plasgeno]
#summglu$`hours elapsed`<-0
summ<-merge(summgal,summglu,by='plasgeno')
mDT1<-melt(summ)
mDT<-data.table(mDT1,colsplit(mDT1$variable,'_',c('sugar','stat')))
mDTlower<-mDT[stat=='lower']
mDTupper<-mDT[stat=='upper']
mDT0<-merge(merge(mDT[stat=='mean'],mDTlower,by=c('plasgeno','sugar')),mDTupper,by=c('plasgeno','sugar'))
tmp<-data.table(`hours elapsed`=c(0,12),sugar=c('glu','gal'))
mDT01<-merge(tmp,mDT0,by='sugar')
mDT01$plasgeno<-factor(mDT01$plasgeno,levels=coregenotypes)
p01FCP<-ggplot(mDT01,aes(x=`hours elapsed`,y=value.x,ymax=value,ymin=value.y,colour=sugar,group=plasgeno))+geom_errorbar(width=1)+geom_point()+ geom_line(colour='grey70',size=1)+scale_colour_manual(values=c('cornflowerblue','indianred'))+facet_wrap(~plasgeno,ncol=4)+ylab('cell density')+theme(strip.text.x=element_blank())+theme(legend.position = c(0.9, 0.3))


pBLANK<-ggplot()

ROWS<-1000;COLS<-ROWS
lmatblank<-data.frame(matrix(nrow=ROWS,ncol=COLS))
N<-0.028
a<-0.55
lmatblank[1:(a*ROWS),1:(N*COLS)]<-1
lmatblank[1:(a*ROWS),((N*COLS)+1):COLS]<-2
lmatblank[(a*ROWS+1):(ROWS),1:COLS]<-3

lmat<-as.matrix(lmatblank)           
pOPF<-list(pBLANK,p00IZF,p01FCP,layout_matrix=lmat)
 
do.call(grid.arrange, pOPF)

h<-4;w<-h*2.2

ggsave(paste0(head.dir,figout,'1707xx-17092x-pOPF-Basic3GenotypePhenotypicOverviewForIntro.pdf'),do.call(grid.arrange, pOPF),width=w,height=h,limitsize=F)
ggsave(paste0(head.dir,figout,'1707xx-17092x-pOPF-Basic3GenotypePhenotypicOverviewForIntro.png'),do.call(grid.arrange, pOPF),width=w,height=h,limitsize=F)

##################
# growth rate heatmap + expression distributions for all clones
##################

DTemp<-qDT[!is.na(plasgeno),c('plasgeno', 'plascopies','gencopies','growth.rate'),with=F]
collim<-signif(as.vector(c(quantile(DTemp $growth.rate,0.005,na.rm=T),quantile(DTemp $growth.rate,0.995,na.rm=T))),3)
DTemp[growth.rate<collim[1]]$growth.rate <-collim[1]
DTemp[growth.rate>collim[2]]$growth.rate <-collim[2]
brks<-signif(seq(collim[1],collim[2],length=5),3)
labs<-round(brks,3)
gencopyunique<-unique(DTemp$gencopies)
DTemp$gencopies<-factor(DTemp$gencopies,levels=gencopyunique[order(gencopyunique,decreasing=TRUE)])
ylabel1<-'plasmid copies [GAL3 GAL80 GAL4]'
xlabel1<-'genome copies [GAL3 GAL80 GAL4]'

p1HIU<-ggplot(DTemp,aes(x=plascopies,gencopies))+ geom_tile(aes(fill=growth.rate))+ 
	scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),guide = guide_legend(title = bquote('growth rate'~h^-1)),breaks=brks, labels=labs,limits = collim) +
	xlab(xlabel1) + ylab(ylabel1) +
	scale_x_discrete(position='top')+scale_y_discrete(position='right')+
	theme(axis.text.x = element_text(size = 10, angle=0))+
	theme(axis.text.y = element_text(size = 10))

# expression density plots
DTemp<-qDT[,c('plasmid.genotype.summary','genome.genotype.summary',dens.cols.new),with=F]
a<-melt(DTemp,id=c('plasmid.genotype.summary','genome.genotype.summary'));setnames(a,'value','density')
a$variable<-gsub('dens.','',a$variable)
b<-cbind(a,colsplit(a$variable,'\\.g',c('AU.FL.VAL','sugar')))
b$sugar<-paste0('g',b$sugar)
cx<-b[,mean(density),by=c('plasmid.genotype.summary','genome.genotype.summary','sugar','AU.FL.VAL')];setnames(cx,'V1','density')
mdf<-merge(b,cx,by=c('plasmid.genotype.summary','genome.genotype.summary','sugar','AU.FL.VAL'))
paste(sample(LETTERS,3),collapse='')
xlabel2<-'pseudo-log10 fluorescence A.U.'
ylabel2<-'density'
p2UWK<-ggplot(mdf[!is.na(plasmid.genotype.summary)],aes(x=AU.FL.VAL,y=density.y))+
	geom_point(aes(x=AU.FL.VAL,y=density.x,colour=sugar),size=0.1,alpha=0.2)+
	geom_line(aes(colour=sugar),size=0.3)+
	scale_y_continuous(breaks=c(0,0.08,0.16))+
	facet_grid(plasmid.genotype.summary~genome.genotype.summary)+
	scale_colour_manual(values=c('cornflowerblue','indianred'))+coord_cartesian(ylim=c(0,0.2))+xlab(xlabel2)+ylab(ylabel2)


pBLANK<-ggplot()

ROWS<-1000;COLS<-ROWS
lmatblank<-data.frame(matrix(nrow=ROWS,ncol=COLS))
a<-0.5
X<-0.04
Y<-0.91
lmatblank[1:(a*ROWS),seq(1,COLS*X,by=1)]<-1
lmatblank[1:(a*ROWS),((COLS*X+1):COLS)]<-2
lmatblank[(a*ROWS+1):ROWS,seq(1,COLS*Y,by=1)]<-3
lmatblank[(a*ROWS+1):ROWS,seq((COLS*Y+1),COLS,by=1)]<-1

lmat<-as.matrix(lmatblank)            

pUKB<-list(pBLANK,p1HIU, p2UWK,layout_matrix=lmat)
#pUKB<-list(p1HIU, p2UWK,ncol=1)
# do.call(grid.arrange, pUKB)
w<-15;h<-w*1.2

ggsave(paste0(head.dir,figout,'1707xx-17092x-pUKB-CopyNumberPhenotypeOverview.png'),do.call(grid.arrange, pUKB),width=w,height=h,limitsize=F)
ggsave(paste0(head.dir,figout,'1707xx-17092x-pUKB-CopyNumberPhenotypeOverview.pdf'),do.call(grid.arrange, pUKB),width=w,height=h,limitsize=F)


##################
# dotplots and boxplots of phenotypes
##################
DTtemp<-DT.merge[genotype.strain.glu!='KV447']
autofl<-DT.merge[genotype.strain.glu=='KV447']
autofluorescent.range.gal<-autofl$yfp.95th.gal+abs(autofl$yfp.5th.gal)
autofluorescent.range.glu<-autofl$yfp.95th.glu+abs(autofl$yfp.5th.glu)
autofluorescent.fracon.gal<-autofl$fracon.gal
autofluorescent.fracon.glu<-autofl$fracon.gal
autofluorescent.PI<-autofl$phenotypic.index

credible.growth.rate.background<-qDT[copies.GAL4==0,summary.func(growth.rate)][,mean+sd*1.96]
credible.yfp.diff.gal<-mean(autofluorescent.range.gal)+1.96*sd(autofluorescent.range.gal)
credible.yfp.diff.glu<-mean(autofluorescent.range.glu)+1.96*sd(autofluorescent.range.glu)
credible.fracon.diff.gal<-mean(autofluorescent.fracon.gal)+1.96*sd(autofluorescent.fracon.gal)
credible.fracon.diff.glu<-mean(autofluorescent.fracon.glu)+1.96*sd(autofluorescent.fracon.glu)
credible.PI.diff<-mean(autofluorescent.PI)+1.96*sd(autofluorescent.PI)



genobin<-DTtemp[,c('copies.GAL3.genome.glu','copies.GAL80.genome.glu','copies.GAL4.genome.glu','copies.GAL3.plasmid.glu','copies.GAL80.plasmid.glu','copies.GAL4.plasmid.glu'),with=F]
setnames(genobin,colnames(genobin),gsub('copies.','',gsub('.glu','',colnames(genobin))))
genobin[,genomecopies:=(apply(genobin[,grepincols(genobin,'genome'),with=F],1,function(x){paste(x,collapse=' ')}))]
genobin[,plasmidcopies:=(apply(genobin[,grepincols(genobin,'plasmid'),with=F],1,function(x){paste(x,collapse=' ')}))]


# flag where plasmid and genome genotypes are identical and leave them out
genobin[,plasgeno:=apply(data.table(genobin$plasmidcopies,genobin$genomecopies),1,function(x){paste(x, collapse=' ')})]
genobin[,genoplas:=apply(data.table(genobin$genomecopies,genobin$plasmidcopies ),1,function(x){paste(x, collapse=' ')})]
genobin[plasgeno!=genoplas]
genobin$totalcopies<-apply(grepcols(grepcols(DTtemp,'total.copies'),'glu'),1,function(x){paste(x,collapse=' ')})

phenotypes<-c('fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal','growth.rate','phenotypic.index')
DT<-data.table(DTtemp[,phenotypes,with=F],genobin)
DT[,GAL80GAL4UNLINKED:=as.factor((plasmidcopies=='0 0 1'&genomecopies=='0 1 0') |(plasmidcopies=='0 1 0'&genomecopies=='0 0 1') | (plasmidcopies=='0 1 1'&genomecopies=='0 1 0') |(plasmidcopies=='0 1 0'&genomecopies=='0 1 1') |  (plasmidcopies=='0 1 1'&genomecopies=='0 0 1') |(plasmidcopies=='0 1 0'&genomecopies=='0 1 1'))]

DAT<-DT[plasgeno!=genoplas]
DAT[order(totalcopies,plasmidcopies,genomecopies)][1:10]

phenotypes<-c('fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal','growth.rate','phenotypic.index')
variablelist<-phenotypes
bylist<-c('totalcopies')

Summarized<-summary.func.all(DAT,variablelist=phenotypes,bylist=bylist)


mPvGmerge1<-melt(Summarized)
mPvGmerge<-data.table(mPvGmerge1,colsplit(mPvGmerge1 $variable,'_',c('phenotype','statistic')))
castDAT<-dcast(mPvGmerge,formula(totalcopies+phenotype~statistic),valuevar='value')
phenofactors<-c('growth.rate','phenotypic.index','fracon.gal','fracon.glu','yfp.mean.gal','yfp.mean.glu')
castDAT[,phenotype:=factor(phenotype,levels=phenofactors)]
castDAT[,se:=sd*sqrt(N)]
castDAT[,ActivatorsTo80:=factor(sapply(totalcopies,function(x){
	ax<-as.numeric(strsplit(x,'\\ ')[[1]])
	ax[1]+ax[3]-ax[2]
}))]
genofactors<-castDAT[phenotype=='growth.rate'][order(mean)]$totalcopies
castDAT[,totalcopies:=factor(totalcopies,levels=genofactors)]

# This one isn't very good
p1a<-ggplot(castDAT,aes(x=mean,y=totalcopies,xmin=mean-se*1.96,xmax=mean+se*1.96))+
	geom_errorbarh(col='grey50')+geom_point(colour='black',size=3)+geom_point(aes(col= ActivatorsTo80),size=2)+
	geom_vline(xintercept=0)+ theme(axis.text.x = element_text(size = 10, angle=30,vjust = 0.5))+
	facet_wrap(~phenotype,ncol=length(unique(phenotypes)),scales='free_x')+
	xlab('mean phenotype')+ylab('genotype - total copies [GAL3 GAL80 GAL4]')+
	theme(panel.spacing=unit(1.5,'lines'), axis.text.x =element_text(size=10))



mPvGmerge1<-melt(DAT[,c('totalcopies',phenotypes),with=F],id='totalcopies')
mPvGmerge1[, ActivatorsToG80:=factor(sapply(totalcopies,function(x){
	ax<-as.numeric(strsplit(x,'\\ ')[[1]])
	ax[1]+ax[3]-ax[2]
}))]
mPvGmerge1[, `sum(activator copies) - sum(repressor copies)`:=ActivatorsToG80]

genofactors<-castDAT[phenotype=='growth.rate'][order(mean)]$totalcopies
mPvGmerge1[,totalcopies:=factor(totalcopies,levels=genofactors)]
phenofactors<-c('growth.rate','phenotypic.index','fracon.gal','fracon.glu','yfp.mean.gal','yfp.mean.glu')
mPvGmerge1[,variable:=factor(variable,levels=phenofactors)]
### NEXT ADD IN THE BINARY GENOTYPE

pXQMa<-ggplot(mPvGmerge1,aes(x=totalcopies,y=value))+
	geom_boxplot(aes(colour=`sum(activator copies) - sum(repressor copies)`))+coord_flip()+
	theme(axis.text.x = element_text(size = 10, angle=30,vjust = 0.5))+
	facet_wrap(~variable,ncol=length((phenotypes)),scales='free_x')+
	ylab('phenotype value')+xlab('GENOTYPE\ntotal copies [GAL3 GAL80 GAL4]')+scale_colour_manual(values= brewer.pal(5, 'Set1'),guide=guide_legend(title.position='top'))+
	theme(panel.spacing=unit(1.5,'lines'), axis.text.x =element_text(size=10),legend.position='bottom')
pXQMb<-ggplot(mPvGmerge1,aes(x=totalcopies,y=value))+
	geom_boxplot(aes(colour=`sum(activator copies) - sum(repressor copies)`))+coord_flip()+
	theme(axis.text.x = element_text(size = 10, angle=30,vjust = 0.5))+
	facet_wrap(~variable,ncol=(length((phenotypes))/2),scales='free_x')+
	ylab('phenotype value')+xlab('GENOTYPE\ntotal copies [GAL3 GAL80 GAL4]')+scale_colour_manual(values= brewer.pal(5, 'Set1'),guide=guide_legend(title.position='top'))+
	theme(panel.spacing=unit(1.5,'lines'), axis.text.x =element_text(size=10),legend.position='bottom')


h<-8;w<-h*2
ggsave(paste0(head.dir,figout,'1707xx-17092x-pXQMb-PhentotypeOverviewAsFunctionOfTotalCopies.png'), pXQMb,width=w,height=h)



##################
# plasmid vs genome copies
##################
phenotypes<-c('fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal','growth.rate','phenotypic.index')

variablelist<-phenotypes
bylist<-c('totalcopies','plasmidcopies','genomecopies','plasgeno','genoplas')
Summarized<-summary.func.all(DAT,variablelist=phenotypes,bylist=bylist)

# check that your contaminant was removed successfully above
# DAT[plasgeno==x[growth.rate_CoV>0.5]$plasgeno]
# DAT[,gr.temp:=paste0(growth.rate)]
# DAT[gr.temp=='0.167383942360909']

# function to correct for unreasonably low errors

x<-splitByVarlist(Summarized,variablelist,bylist)$growth.rate

DATX<-copy(x)
qplot(DATX$fracon.gal_mean,DATX$fracon.gal_CoV)#qplot(DATX$growth.rate_mean,DATX$growth.rate_CoV)
varcols<-c('growth.rate_mean','growth.rate_CoV','growth.rate_sd')#varcols<-c('fracon.gal_mean','fracon.gal_CoV','fracon.gal_sd')
varname<-'growth.rate'
ErrorCorrection<-function(DT,varcols=NULL,varname=NULL){
	# idea with this function is to take an arbitrary dataset with given mean, sd and CoV, flag samples with low error, correct that error
	# it assumes that CoV and mean are linearly anticorrelated in log space, therefore this transformation is used to flag samples that
	# have very low CoV given their mean 
	
	# regarding input: varcols must have the names of the columns for mean, CoV and sd respectively as a c() object
	# the DT is an object
	DATX<-copy(DT)
	varcols1<-c('mean','CoV','sd')
	if(!is.null(varcols))setnames(DATX,varcols,varcols1)
	qplot(DATX$mean,DATX$CoV)

	mod1<-lm(log(CoV)~log(mean), DATX[,varcols1,with=F])

	# we're flagging very low error values, but we don't want to penalize ourselves by 
	# a full fit to the logged data, so we take the lower bound of the standard error estimate
	# of the effect for the log of the mean
	lower<-0#as.data.frame(summary(mod1)$coefficients)$`Std. Error`[2]
	resid<-mod1$residuals
	fit<-exp((mod1$fitted.values-lower))
	
	# set a cutoff of 1 standard deviation below the expectation for "a sample with very low measurement error"
	# this cutoff, and the manipulation of the distribution of errors i.e. to exclude extreme outliers, could be much more sophisticated. this seems to work for now.
	cutoff<-(-sd(resid)*0.7)
	flag<-resid<cutoff
	ggplot(data.table(DATX,fit,flag),aes(x=mean,y=CoV,col=flag))+geom_point()+geom_line(aes(x=mean,y=fit))
	temp<-data.table(fit,resid,flag,DATX[,varcols1,with=F])
	temp[,sd.fit:=fit*mean]
	temp[flag==FALSE,sd.new:=sd]
	temp[flag==TRUE,sd.new:=sd.fit]
	temp[,check:=sd==sd.new]
	out<-DATX[,sd.temp:=temp$sd.new]
	out[, stdev.corrected:=sd!=sd.temp]
	out[,sd:=sd.temp]
	out[,sd.temp:=NULL]
	out[,CoV:=sd/mean]
	if(!is.null(varname))setnames(out,c(varcols1,'stdev.corrected'),paste0(varname,"_",c(varcols1,'stdev.corrected')))
	data.frame(out)
}

x<-splitByVarlist(Summarized,variablelist,bylist)$growth.rate
plasmid.vs.genome.comparisons<-lapply(splitByVarlist(Summarized,variablelist,bylist),function(x){
	x<-data.table(x)
	DATX<-x[,grepincols(x,c(bylist,'_mean','_sd','_N','_CoV')),with=F]
	variable<-strsplit(grepincols(x,c('mean')),'_')[[1]][1]
	DT<-data.table(ErrorCorrection(DATX,varcols=grepincols(DATX,c('_mean','_CoV','_sd')),varname=variable))
	DT$mean<-DT[,grepincols(DT,'_mean'),with=F]
	DT$sd<-DT[,grepincols(DT,'_sd'),with=F]
	DT$N<-DT[,grepincols(DT,'_N'),with=F]

	matches<-match(DT $plasgeno, DT $genoplas)
	matchDT<-data.table(plas=1:nrow(DT),gen=matches)
	matchDT[,dupflag:=apply(data.frame(plas,gen),1,function(x)paste0(min(x),max(x)))]
	
	diff<-DT[,'mean',with=F]-DT[matches,'mean',with=F]
	se<-sqrt((DT[,'sd',with=F]/sqrt(DT[,'N',with=F]))^2+(DT[matches,'sd',with=F]/sqrt(DT[matches,'N',with=F]))^2)
	t<-(diff/se)
	df<-DT[,'N',with=F]+DT[matches,'N',with=F]-2
	p<-apply(data.frame(t,df),1,function(xa){
		pt(-abs(xa[1]),xa[2])
		})
	if(grepl('yfp.mean.glu',variable)){
		vartemp<-credible.yfp.diff.glu
		p[abs(diff)<vartemp]<-0.5
		diff[abs(diff)<vartemp]<-0
	}
	if(grepl('yfp.mean.gal',variable)){
		vartemp<-credible.yfp.diff.gal
		p[abs(diff)<vartemp]<-0.5
		diff[abs(diff)<vartemp]<-0
	}

	if(grepl('fracon.gal',variable)){
		vartemp<-credible.fracon.diff.gal
		p[abs(diff)<vartemp]<-0.5
		diff[abs(diff)<vartemp]<-0
	}

	if(grepl('fracon.glu',variable)){
		vartemp<-credible.fracon.diff.glu
		p[abs(diff)<vartemp]<-0.5
		diff[abs(diff)<vartemp]<-0
	}

	if(grepl('phenotypic.index',variable)){
		vartemp<-credible.PI.diff
		p[abs(diff)<vartemp]<-0.5
		diff[abs(diff)<vartemp]<-0
	}

	if(grepl('growth.rate',variable)){
		vartemp<-credible.growth.rate.background
		p[DT$mean<vartemp|DT[matches]$mean<vartemp]<-0.5
		diff[DT$mean<vartemp|DT[matches]$mean<vartemp]<-0
	}
	
	x$diff<-diff
	x$se<-se
	x$df<-df
	x$p.value<-p*2
	x$p.adj<-p.adjust(p,'bonferroni')
	x$t<-t
	renamecols<-c('diff','se','df','p.value','p.adj','t')
	
	
	setnames(x,renamecols, paste0(variable,'_diffs_',renamecols))
	out1<-data.table(x,matchDT)[!duplicated(dupflag)]
	out<-out1[,!colnames(matchDT),with=F]
	out
})

# check that this worked correctly
# this does not work for all samples when you do the error correction...
DTtest<-plasmid.vs.genome.comparisons$growth.rate
plasgenoA<-'1 1 1 1 0 1' 
plasgenoB<-'1 0 1 1 1 1' 
A<-DAT[plasgeno== plasgenoA,'growth.rate']
B<-DAT[plasgeno== plasgenoB,'growth.rate']
test<-t.test(A,B,var.equal=TRUE,alternative='two.sided')
data.table(paste0(c(test$statistic,test$parameter,test$p.value)), 
	paste0(DTtest[plasgeno==plasgenoA|plasgeno==plasgenoB,c('growth.rate_diffs_t','growth.rate_diffs_df','growth.rate_diffs_p.value')]))

bylistbylist<-lapply(1:length(plasmid.vs.genome.comparisons),function(i)return(bylist))
mPvGmerge1<-melt(merge_recurse2(dt.list=plasmid.vs.genome.comparisons,by.list=bylistbylist))
mPvGmerge<-data.table(mPvGmerge1,colsplit(mPvGmerge1 $variable,'_',c('phenotype','statistic')))
castDAT<-dcast(mPvGmerge,formula(plasgeno+phenotype~statistic),valuevar='value')
castDAT[,p.adj.all:=p.adjust(diffs_p.value,'fdr')]
sigcutoff<-0.15
castDAT[,sig:= diffs_p.adj <sigcutoff]
phenofactors<-c('growth.rate','phenotypic.index','fracon.gal','fracon.glu','yfp.mean.gal','yfp.mean.glu')
castDAT[,phenotype:=factor(phenotype,levels=phenofactors)]

p<-ggplot(castDAT,aes(x=diffs_diff,y=plasgeno,col=-(log10(diffs_p.adj)),xmin=diffs_diff-diffs_se*1.96,xmax=diffs_diff+diffs_se*1.96))+geom_point()+geom_errorbarh()+
	scale_colour_gradientn(colours=rainbow(5))+geom_vline(xintercept=0)+ theme(axis.text.x = element_text(size = 10, angle=30,vjust = 0.5))+
	facet_wrap(~phenotype,ncol=length(unique(phenotypes)),scales='free_x')

p<-ggplot(castDAT,aes(x=diffs_diff,y=plasgeno,xmin=diffs_diff-diffs_se*1.96,xmax=diffs_diff+diffs_se*1.96))+geom_point(aes(col=as.factor(sig)))+geom_errorbarh(col='grey50')+
	geom_vline(xintercept=0)+ theme(axis.text.x = element_text(size = 10, angle=30,vjust = 0.5))+scale_colour_manual(values=c('cornflowerblue','indianred'))+
	facet_wrap(~phenotype,ncol=length(unique(phenotypes)),scales='free_x')


castDAT2<-merge(castDAT,castDAT[phenotype=='growth.rate',c('plasgeno', 'sig'),with=F],by='plasgeno')
p<-ggplot(castDAT2[sig.y==T],aes(x=diffs_diff,y=plasgeno,xmin=diffs_diff-diffs_se*1.96,xmax=diffs_diff+diffs_se*1.96))+geom_point(aes(col=as.factor(sig.x)))+geom_errorbarh(col='grey50')+
	geom_vline(xintercept=0)+ theme(axis.text.x = element_text(size = 10, angle=30,vjust = 0.5))+scale_colour_manual(values=c('cornflowerblue','indianred'))+
	facet_wrap(~phenotype,ncol=length(unique(phenotypes)),scales='free_x')



### make all growth rate differences greater than 0 to help make more sense of the genotypes and corresponding underlying expression phenotypes; requires some massaging

grflag<-castDAT[phenotype=='growth.rate',c('diffs_diff','plasgeno'),with=F]
grflag[, genotype:=lapply(1:nrow(grflag),function(i){
	x<-c(diffs_diff[i],plasgeno[i])
	a<-x[2]
	b<-1
	if(as.numeric(x[1])<0){
		a<-paste0(c(strsplit(a,'\\ ')[[1]][4:6],strsplit(a,'\\ ')[[1]][1:3]),collapse=' ')
		b<-(-1)
	}
	a
})][,genotype:=as.factor(unlist(genotype))]
grflag[, flipfac:=lapply(1:nrow(grflag),function(i){
	x<-c(diffs_diff[i],plasgeno[i])
	a<-x[2]
	b<-1
	if(as.numeric(x[1])<0){
		a<-paste0(c(strsplit(a,'\\ ')[[1]][4:6],strsplit(a,'\\ ')[[1]][1:3]),collapse=' ')
		b<-(-1)
	}
	b
})]

castDAT3<-merge(castDAT,grflag[,c('plasgeno','genotype','flipfac')],by='plasgeno')
castDAT3[,diffs_flip:=diffs_diff*as.numeric(flipfac)]

p<-ggplot(castDAT3,aes(x=diffs_flip,y=genotype,col=-(log10(diffs_p.adj)),xmin=diffs_flip-diffs_se*1.96,xmax=diffs_flip+diffs_se*1.96))+geom_point()+geom_errorbarh()+
	scale_colour_gradientn(colours=rainbow(5))+geom_vline(xintercept=0)+ theme(axis.text.x = element_text(size = 10, angle=30,vjust = 0.5))+
	facet_wrap(~phenotype,ncol=length(unique(phenotypes)),scales='free_x')

p<-ggplot(castDAT3,aes(x=diffs_flip,y=genotype,xmin=diffs_flip-diffs_se*1.96,xmax=diffs_flip+diffs_se*1.96))+geom_point(aes(col=as.factor(sig)))+geom_errorbarh(col='grey50')+
	geom_vline(xintercept=0)+ theme(axis.text.x = element_text(size = 10, angle=30,vjust = 0.5))+scale_colour_manual(values=c('cornflowerblue','indianred'))+
	facet_wrap(~phenotype,ncol=length(unique(phenotypes)),scales='free_x')

# Polish up dataset with new levels to significance, genotypic summary information
castDAT4<-merge(castDAT3, castDAT3[phenotype=='growth.rate',c('genotype', 'sig'),with=F],by='genotype')

castDAT4[,totalcopies:=lapply(genotype,function(x){
	a<-as.numeric(strsplit(as.character(x),'\\ ')[[1]])
	paste0(a[1:3]+a[4:6],collapse=' ')
})]
x<-castDAT4$genotype[1]
castDAT4[,ActivatorsToG80:=lapply(genotype,function(x){
	a<-as.numeric(strsplit(as.character(x),'\\ ')[[1]])
	sum(a[c(1,3,4,6)])-sum(a[c(2,5)])
})]
castDAT4[,GAL4.binary:=lapply(genotype,function(x){
	a<-as.numeric(strsplit(as.character(x),'\\ ')[[1]])
	sum(a[c(3,6)])>0
})]
castDAT4[,significant:=factor(apply(data.frame(unlist(GAL4.binary),sig.x),1,function(x){
	a<-x[2]
	if(x[1]!=TRUE)a<-'GAL4 deletion'
	a
}),levels=c('FALSE','TRUE','GAL4 deletion'))]
castDAT4$genotype<-droplevels(factor(castDAT4$genotype,levels=unique(castDAT3$genotype)))
castDAT4[,genotype2:=paste0(ActivatorsToG80,'  ||  ',totalcopies,'  ||  ',genotype)]
genoorder1<-c(paste0(castDAT4[GAL4.binary==F &phenotype=='growth.rate',c('genotype','diffs_flip'),with=F][order(diffs_flip)]$genotype),
	paste0(castDAT4[GAL4.binary==T &phenotype=='growth.rate',c('genotype','diffs_flip'),with=F][order(diffs_flip)]$genotype))
castDAT4$genotype<-factor(castDAT4$genotype,levels=genoorder1)
genoorder2<-c(paste0(castDAT4[GAL4.binary==FALSE&phenotype=='growth.rate',c('genotype2','diffs_flip'),with=F][order(diffs_flip)]$genotype2),
	paste0(castDAT4[phenotype=='growth.rate'&GAL4.binary==T,c('genotype2','diffs_flip'),with=F][order(diffs_flip)]$genotype2))
castDAT4$genotype2<-factor(castDAT4$genotype2,levels=genoorder2)
castDAT4[,`FDR < 0.15`:=sig.x]

# make pretty plots
p1a<-ggplot(castDAT4,aes(x=diffs_flip,y=genotype,xmin=diffs_flip-diffs_se*1.96,xmax=diffs_flip+diffs_se*1.96))+
	geom_errorbarh(col='grey50')+geom_point(colour='black',size=3)+geom_point(aes(col=`FDR < 0.15`),size=2)+
	geom_vline(xintercept=0)+ theme(axis.text.x = element_text(size = 10, angle=30,vjust = 0.5))+scale_colour_manual(values=c('cornflowerblue','indianred','grey50'))+
	facet_wrap(~phenotype,ncol=length(unique(phenotypes)),scales='free_x')+
	xlab('difference in phenotypic means (given plasmid:genome genotype compared to the complementary genome:plasmid genotype)')+ylab('genotype')+
	theme(panel.spacing=unit(1.5,'lines'), axis.text.x =element_text(size=10))
#&ActivatorsToG80>=0
p1b<-ggplot(castDAT4,aes(x=diffs_flip,y=genotype2,xmin=diffs_flip-diffs_se*1.96,xmax=diffs_flip+diffs_se*1.96))+
	geom_errorbarh(col='grey50')+geom_point(colour='black',size=3)+geom_point(aes(col= `FDR < 0.15`),size=2)+#geom_point(aes(col= significant),size=1)+
	geom_vline(xintercept=0)+ theme(axis.text.x = element_text(size = 10, angle=30,vjust = 0.5))+scale_colour_manual(values=c('cornflowerblue','indianred'))+	
	#scale_colour_manual(values=c('cornflowerblue','indianred','grey50'))+
	facet_wrap(~phenotype,ncol=length(unique(phenotypes)),scales='free_x')+
	xlab('difference in phenotypic means of given [plasmid:genome] - complementary [genome:plasmid] genotypes')+ylab('genotype')+
	theme(panel.spacing=unit(1.5,'lines'), axis.text.x =element_text(size=10))+ggtitle('all genotypes where [plasmid:genome] differs from [genome:plasmid]')+theme(plot.title = element_text(face = "plain"))

castDAT4[,binarycopies:=sapply(totalcopies,function(x){
	xa<-strsplit(x,'\\ ')[[1]]
	paste0(sapply(xa,function(boo){as.numeric(boo>0)}),collapse=' ')
})]
castDAT4a<-castDAT4[binarycopies=='1 1 1']
p1c<-ggplot(castDAT4a,aes(x=diffs_flip,y=genotype2,xmin=diffs_flip-diffs_se*1.96,xmax=diffs_flip+diffs_se*1.96))+
	geom_errorbarh(col='grey50')+geom_point(colour='black',size=3)+geom_point(aes(col= `FDR < 0.15`),size=2)+#geom_point(aes(col= significant),size=1)+
	geom_vline(xintercept=0)+ theme(axis.text.x = element_text(size = 10, angle=30,vjust = 0.5))+scale_colour_manual(values=c('cornflowerblue','indianred'))+	
	#scale_colour_manual(values=c('cornflowerblue','indianred','grey50'))+
	facet_wrap(~phenotype,ncol=length(unique(phenotypes)),scales='free_x')+
	xlab('difference in phenotypic means of given [plasmid:genome] - complementary [genome:plasmid] genotypes')+ylab('genotype')+
	theme(panel.spacing=unit(1.5,'lines'), axis.text.x =element_text(size=10))+ggtitle('inducible genotypes [GAL3.WT GAL80.WT GAL4.WT] where [plasmid:genome] differs from [genome:plasmid]')+theme(plot.title = element_text(face = "plain"))


paste(sample(LETTERS,3),collapse='')


pBLANK<-ggplot()
ROWS<-100;COLS<-ROWS
lmatblank<-data.frame(matrix(nrow=ROWS,ncol=COLS))
a<-0.6
b<-0.4

lmatblank[1:(a*ROWS),1:COLS]<-1
lmatblank[(a*ROWS+1):ROWS,1:COLS]<-2

lmat<-as.matrix(lmatblank)            
pUWS<-list(p1b,p1c,layout_matrix=lmat)
do.call(grid.arrange,pUWS)
h<-10;w<-h*2.2
ggsave(paste0(head.dir,figout,'1707xx-17092x-pUWS-PlasmidVsGenomeCopies2.png'),do.call(grid.arrange, pUWS),width=w,height=h)
ggsave(paste0(head.dir,figout,'1707xx-17092x-pUWS-PlasmidVsGenomeCopies2.pdf'),do.call(grid.arrange, pUWS),width=w,height=h)





#### OLD

means<-DAT[,lapply(.SD,mean,na.rm=T),by=c('totalcopies','plasmidcopies','genomecopies','plasgeno'),.SDcols=phenotypes][order(totalcopies,plasmidcopies,genomecopies)]
sds<-DAT[,lapply(.SD,sd,na.rm=T),by=c('totalcopies','plasmidcopies','genomecopies','plasgeno'),.SDcols=phenotypes][order(totalcopies,plasmidcopies,genomecopies)]
a<-melt(means[seq(from=1,to=nrow(means),by=2)])[order(totalcopies)]
b<-melt(means[seq(from=2,to=nrow(means),by=2)])[order(totalcopies)]
asd<-melt(sds[seq(from=1,to=nrow(means),by=2)])[order(totalcopies)]
bsd<-melt(sds[seq(from=2,to=nrow(means),by=2)])[order(totalcopies)]

df<-data.table(plasgeno=factor(a$plasgeno),totalcopies=a$totalcopies,variable=a$variable,copy.in.genome=a$value,copy.in.plasmid=b$value,gsd=asd$value,psd=bsd$value)[order(plasgeno,totalcopies,variable)]
df$diff<-df$copy.in.plasmid-df$copy.in.genome
df$N<-8
df$differr<-sqrt(((df$gsd/sqrt(df$N))^2+(df$psd/sqrt(df$N))^2))
df$upperdiff<-df$diff+df$differr
df$lowerdiff<-df$diff-df$differr

# DO TTESTS 

numericcols<-c('copy.in.genome','copy.in.plasmid','gsd','psd','N','N')

ttests<-data.table(t(apply(df[,numericcols,with=F],1,function(x){
	names(x)<-c('m1','m2','s1','s2','n1','n2')
	with(as.list(x),{t.test2(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)})
})))

df1<-cbind(df,ttests)
df1$p.adj<-p.adjust(df1$p.value,'bonferroni')
df1$significant1<-df1$p.adj<0.05
df1$significant2<-df1$p.adj<0.01
df1$significant3 <-factor(apply(data.table(df1$significant1,df1$significant2),1,sum))
dumdf<-data.table(significant=c('NS','p<0.05','p<0.01'), significant3 =c('0','1','2'))
df2<-merge(df1,dumdf,by='significant3')[order(plasgeno)]
df2$significant<-factor(df2$significant)


DTz<-data.table(DAT[,lapply(.SD,zscore),.SDcols=phenotypes],DAT[,c('totalcopies','plasmidcopies','genomecopies','plasgeno'),with=F])
means<-DTz[,lapply(.SD,mean,na.rm=T),by=c('totalcopies','plasmidcopies','genomecopies','plasgeno'),.SDcols=phenotypes][order(totalcopies,plasmidcopies,genomecopies)]
sds<-DTz[,lapply(.SD,sd,na.rm=T),by=c('totalcopies','plasmidcopies','genomecopies','plasgeno'),.SDcols=phenotypes][order(totalcopies,plasmidcopies,genomecopies)]
Ns<-DTz[,lapply(.SD,length),by=c('totalcopies','plasmidcopies','genomecopies','plasgeno'),.SDcols=phenotypes][order(totalcopies,plasmidcopies,genomecopies)]
a<-melt(means[seq(from=1,to=nrow(means),by=2)])[order(totalcopies)]
b<-melt(means[seq(from=2,to=nrow(means),by=2)])[order(totalcopies)]
asd<-melt(sds[seq(from=1,to=nrow(means),by=2)])[order(totalcopies)]
bsd<-melt(sds[seq(from=2,to=nrow(means),by=2)])[order(totalcopies)]

dfz<-data.table(plasgeno =factor(a$plasgeno),totalcopies=a$totalcopies,variable=a$variable,copy.in.genome=a$value,copy.in.plasmid=b$value,gsd=asd$value,psd=bsd$value)[order(plasgeno,totalcopies,variable)]
dfz $diff<-dfz $copy.in.plasmid-dfz $copy.in.genome
dfz $N<-8
dfz $differr<-sqrt(((dfz $gsd/sqrt(dfz $N))^2+(dfz $psd/sqrt(dfz $N))^2))
dfz $upperdiff<-dfz $diff+ dfz $differr
dfz $lowerdiff<-dfz $diff-dfz $differr
dfz$meanz<-apply(data.table(dfz$copy.in.plasmid,dfz$copy.in.genome),1,mean)
dfz[,col1:=lapply(.SD,function(x){ecdf(x)(x)}),by='variable',.SDcols=c('meanz')]
dfz[,col2:=lapply(.SD,function(x){ecdf(x)(x)}),by='variable',.SDcols=c('diff')]


DT<-merge(df2,dfz,by=c('plasgeno','totalcopies','variable'))
alph<-0.5
p1<-ggplot(DT,aes(x=copy.in.genome.x,y=copy.in.plasmid.x))+
	geom_errorbar(aes(ymax = copy.in.plasmid.x+psd.x, ymin = copy.in.plasmid.x-psd.x),alpha=alph)+
	geom_errorbarh(aes(xmax = copy.in.genome.x+gsd.x, xmin = copy.in.genome.x-gsd.x),alpha=alph)+
	geom_point(aes(colour=log10(p.adj)))+
	geom_abline(slope=1,intercept=0,col='grey50')+
	scale_colour_gradientn(colours=rainbow(5))+
	facet_wrap(~variable,ncol=5,scales='free')+theme(panel.spacing=unit(1.5,'lines'), axis.text.x =element_text(size=10))+xlab('value when genotype is in genome')+ylab('value when\ngenotype is\nin plasmid')

p2<-ggplot(DT,aes(y=diff.x,x=plasgeno))+
	geom_hline(yintercept=0,col='grey50',alpha=0.5)+
	geom_errorbar(aes(ymax=upperdiff.x,ymin=lowerdiff.x))+
	geom_point(aes(colour= col1),size=2)+
	scale_colour_gradientn(colours=rainbow(5))+coord_flip()+
	facet_grid(~variable,scales='free')

p3<-ggplot(DT,aes(y=diff.x,x=plasgeno))+
	geom_hline(yintercept=0,col='grey50',alpha=0.5)+
	geom_point(colour='black',size=3)+ geom_point(colour='white',size=2.5)+
	geom_errorbar(aes(ymax=upperdiff.x,ymin=lowerdiff.x))+
	geom_point(aes(colour= log10(p.adj),shape=significant),size=1.5)+
	scale_colour_gradientn(colours=rainbow(5))+coord_flip()+
	facet_grid(~variable,scales='free')+theme(panel.spacing=unit(1.5,'lines'), axis.text.x =element_text(size=10))+ylab('difference in means [plasmid genotype] - isogenic [genome genotype]')+xlab('genotype')

pBLANK<-ggplot()
ROWS<-100;COLS<-ROWS
lmatblank<-data.frame(matrix(nrow=ROWS,ncol=COLS))
a<-0.01
b<-0.3
c<-1
lmatblank[1:(a*ROWS),1:COLS]<-1
lmatblank[(a*ROWS+1):(b*ROWS),1:COLS]<-2
lmatblank[(b*ROWS+1):(c*ROWS),1:COLS]<-3

lmat<-as.matrix(lmatblank)            
pVNF<-list(pBLANK,p1,p3,layout_matrix=lmat)
do.call(grid.arrange, pVNF)

w<-14;h=w*1
ggsave(paste0(head.dir,figout,'1707xx-17092x-pVNF-PlasmidVsGenomeCopies.pdf'),do.call(grid.arrange,pVNF),width=w,height=h)
ggsave(paste0(head.dir,figout,'1707xx-17092x-pVNF-PlasmidVsGenomeCopies.pdf'),do.call(grid.arrange,pVNF),width=w,height=h)

paste(sample(LETTERS,3),collapse='')






