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
rm(letters)

applyPaste<-function(x,collapse=''){
	apply(x,1,function(x)paste0(x,collapse=collapse))
}

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
	i<-1
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

radial_plot_space<-function(DT,genos){
	genos1<-DT[,genos,with=F]
	genos2<-lapply(1:length(genos),function(i){
		paste0(genos[1:i],collapse='')
	})
	
	genos3<-lapply(1:length(genos),function(i){
		apply(DT[,genos[1:i],with=F],1,function(x)paste0(x,collapse=''))
	})
	genos4<-data.table(dummy=1:nrow(DT))
	for(i in 1:length(genos3)){
		genos4<-cbind(genos4,genos3[[i]])
	};genos4[,dummy:=NULL]
	colnames(genos4)<-paste0('x',unlist(genos2))
	
	radiusID<-data.table(apply(genos4,2,function(goo){
		
		out<-c(TRUE)
		for(i in 2:length(goo)){
			out<-c(out,(goo[i]!=goo[i-1]))
		}
		out
	}))
	
	nums<-1:(length(genos)+1)
	colnames(radiusID)<-paste0('r_',nums[-1])
	radiusID[,r_1:=c(TRUE,rep(FALSE,(nrow(radiusID)-1)))]
	radiusID<-(radiusID[,order(colnames(radiusID)),with=F])
	
	i<-1
	radius<-sapply (1:length(nums),function(i){
		x<-unlist(radiusID[,colnames(radiusID)[i],with=F])
		b<-gsub('TRUE',nums[i],x)
		as.numeric(gsub('FALSE',NA,b))
	})
	colnames(radius)<-colnames(radiusID)
	x<-radius[,1]
	starts<-apply(radius,2,function(x){
		a<-seq(0,(1-1/length(na.exclude(x))),length=length(na.exclude(x)))*2*pi
		x[!is.na(x)]<-a
		x
	})
	
	x<-starts[,2]
	ends<-apply(starts,2,function(x){
		b<-na.exclude(x)
		a<-pi*2
		if(length(b)>1)	a<-c(b[2:length(b)],pi*2)
		x[!is.na(x)]<-a
		x
	
	})
	
	colnames(starts)<-paste0('start_',nums)
	colnames(ends)<-paste0('end_',nums)
	DT1<-data.table(genos1,radius,starts,ends)
	merge(DT,DT1,by=genos)	
	
}
	
walsh_matrix = function(order){
    m = 1
    for(i in 1:order){
        m = rbind(  cbind(m,m),
                    cbind(m,m*-1)
        )
    }
    m
}

relative_matrix = function(order){
    m = 1
    for(i in 1:order){
        m = rbind(  cbind(m,0*m),
                    cbind(-1*m,m)
        )
    }
    m
}

weighting_matrix = function(order){
    if (order < 1) {
        m = 1
    } else {
        m = 1
        zm = 0
        for(i in 1:order){
            m = rbind(  cbind(0.5*m,zm),
                        cbind(zm,m*-1))
            zm = matrix(0,nrow(m),ncol(m))
        }
    }
    m
}


varexplained<-function(y1,yexp, DT=NULL){
	if(!is.null(DT)){df<-data.frame(DT);y1<-df[,1];yexp=df[,2]
		1-sum((y1-yexp)^2,na.rm=T)/sum(((y1-mean(y1,na.rm=T))^2),na.rm=T)
	}else(1-sum((y1-yexp)^2,na.rm=T)/sum(((y1-mean(y1,na.rm=T))^2),na.rm=T))
	
	}

logconv<-function(DT,x='x',backg=0){
	DT<-data.table(data.frame(DT))
	DT[,lx:=log(abs(x-backg))]
	data.table(data.frame(DT[,!'x',with=F]))
}
x<-'observed_mean'
DT<-copy(DATcur1)

logconv2<-function(DT,x='x',backg=0,backgsd =0,genos,flip=F){
	if(x!='x')setnames(DT,x,'x')
	DT<-data.table(data.frame(DT))
	credmin<-backg+backsd
	DT[,corx:=x]
	DT[corx<credmin,corx:=backg+10^-16]
	DT[,subx:=abs(x-backg)]
	DT[,flip:=(x-backg)/abs(x-backg)]
	DT[,lx:=log(subx)]
	DT[,lx:=log(abs(x-backg))]
	a<-data.table(data.frame(DT[,c('lx',genos),with=F]))
	if(flip==T){a<-data.table(data.frame(DT[,c('lx','flip',genos),with=F]))}
	a
}
expMod<-function(DT,mod,backg=0,genos,nameFit=NULL){
	DT<-data.table(data.frame(DT[,mean(x,na.rm=T),by=genos][,genos,with=F]))
	pred<-predict(mod,DT,se.fit=T)
	DT[,c('fit_mean','fit_se'):=list(exp(pred$fit)+backg,exp(pred$fit)*pred$se.fit)]
	DT[,c('fit_upper','fit_lower'):=list(fit_mean+fit_se,fit_mean-fit_se)]
	if(!is.null(nameFit))setnames(DT,c('fit_mean','fit_se','fit_upper','fit_lower'),paste0(nameFit,'__',c('fit_mean','fit_se','fit_upper','fit_lower')))
	data.table(data.frame(DT))
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
qDT<-data.table(qDTx,DTcopies,DTbinary)


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

varexplained<-function(y1,yexp, DT=NULL){
	if(!is.null(DT)){df<-data.frame(DT);y1<-df[,1];yexp=df[,2]
		1-sum((y1-yexp)^2,na.rm=T)/sum(((y1-mean(y1,na.rm=T))^2),na.rm=T)
	}else(1-sum((y1-yexp)^2,na.rm=T)/sum(((y1-mean(y1,na.rm=T))^2),na.rm=T))
	
	}

totvar<-varexplained(qDT$growth.rate,mod1$fitted.values)
qplot(mod1$fitted.values,qDT$growth.rate)+geom_abline(intercept=0,slope=1)+geom_vline(xintercept=credible.growth.rate.background,col='grey70')+
	geom_hline(yintercept=credible.growth.rate.background,col='grey70')+xlab('within-genotype means')+ylab('observed across N=8 observations')+
	xlim(c(0,0.25))+ylim(c(0,0.25))+ggtitle('overview of plasmid vs genome copies expt')+
	annotate('text',x=0.07,y=0.20,label=paste0('total variation explained =\n',round(totvar,3)))


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
# TukeyHSD(aov(mod2))
dens.cols.new<-colnames(DT.merge)[grep('dens',colnames(DT.merge))]

###################################################
### Define genotypic ensembles
###################################################

phenotype<-'growth.rate'
backg<-mean(qDT[copies.GAL4==0]$growth.rate)
backgsd<-sd(qDT[copies.GAL4==0]$growth.rate)

	### binary genotypes
	genocols<-c('binarycount','binary.GAL4','binary.GAL3','binary.GAL80')
	genosalt<-c('plasgeno')
	DAT1<-qDT[,c(phenotype,genocols),with=F]
	N<-length(genocols)-1
	GENOS<-paste0(LETTERS[1:N])
	genos<-paste0(letters[1:N])
	GENOSx<-quote(mget(GENOS))
	genosx<-quote(mget(genos))
	setnames(DAT1,colnames(DAT1),c('x','genotypei',GENOS))
	dum<-DAT1[,GENOS,with=F][,lapply(.SD,function(x)as.factor(1-as.numeric(x)))];colnames(dum)<-genos
	DAT2<-data.table(DAT1,dum,qDT[,genosalt,with=F])
	DATbinary<-copy(DAT2)
	bingeno<-genos
	
	### total copies
	DAT0 <-qDT[,c(phenotype,'totalcopies','plasgeno'),with=F]
	DAT0$arb<-1:nrow(DAT0)
	a<-(t(unlist(sapply(DAT0 $totalcopies,function(x)strsplit(x,'\\ ')[[1]]))))
	b<-data.table(a);b$totalcopies<-rownames(a);b$arb<-1:nrow(b)
	clonenamedcols<-c('GAL3','GAL80','GAL4')
	colnames(b)<-c(clonenamedcols,'totalcopies','arb')
	N<-sum(apply(b[,clonenamedcols,with=F],2,function(x)length(unique(x))-1))
	DAT0x<-merge(DAT0,b,by=c('totalcopies','arb'))
	genocols<-c('totalcopies','GAL4','GAL3','GAL80')
	genosalt<-c('plasgeno')
	DAT1<-DAT0x[,c(phenotype,genocols,genosalt),with=F]
	N0<-length(genocols)-1
	GENOS<-paste0(LETTERS[1: N0])
	genos<-paste0(letters[1: N0])
	GENOSx<-quote(mget(GENOS))
	genosx<-quote(mget(genos))
	setnames(DAT1,colnames(DAT1),c('x','genotypei',GENOS,genosalt))
	dum<-DAT1[,GENOS,with=F][,lapply(.SD,function(x)as.factor(2-as.numeric(x)))];colnames(dum)<-genos
	DAT2<-data.table(DAT1,dum)
	DATtotalcopies<-copy(DAT2)
	totcopgenos<-genos

	### full genotypes
	genocols<-c('plasgeno',"copies.GAL4.plasmid","copies.GAL4.genome","copies.GAL3.plasmid","copies.GAL3.genome","copies.GAL80.plasmid","copies.GAL80.genome")
	genosalt<-c('plasgeno')
	DAT1 <-qDT[,unique(c(phenotype,genosalt,genocols)),with=F]
	N<-length(genocols)-1
	GENOS<-paste0(LETTERS[1:N])
	genos<-paste0(letters[1:N])
	GENOSx<-quote(mget(GENOS))
	genosx<-quote(mget(genos))
	setnames(DAT1,colnames(DAT1),c('x','genotypei',GENOS))
	N0<-sum(DAT1[,apply(data.frame(eval(GENOSx)),2,function(x){length(unique(x))-1})])
	dum<-DAT1[,GENOS,with=F][,lapply(.SD,function(x)as.factor(1-as.numeric(x)))];colnames(dum)<-genos
	DAT2<-data.table(DAT1,dum)
	DAT2[,plasgeno:=genotypei]
	backg=DAT2[A==0&B==0,mean(x)];backsd=DAT2[A==0&B==0,sd(x)]
	DATfullgenotypes<-copy(DAT2)
	fullgenos<-genos

DATbinary[,genotype_ensemble:=as.character('binary')]
DATbinary[,genotypei:=applyPaste(data.frame(A,B,C),' ')]
DATtotalcopies[,genotype_ensemble:=as.character('total copies')]
DATtotalcopies[, genotypei:=applyPaste(data.frame(A,B,C),' ')]
DATfullgenotypes[,genotype_ensemble:=as.character('full')]
DATfullgenotypes[,genotypei:=applyPaste(data.frame(A,B,C,D,E,F),' ')]


genolist<-list(bingeno,totcopgenos,fullgenos)
dtlist<-list(DATbinary,DATtotalcopies,DATfullgenotypes)
phenotype<-'growth.rate'
backg<-mean(qDT[copies.GAL4==0]$growth.rate)
backgsd<-sd(qDT[copies.GAL4==0]$growth.rate)
i<-1
nameFits<-data.table(order=1:6,interaction_order=c('first','second','third','fourth','fifth','sixth'))
genotype_ensembles<-c('binary','total copies','full')


###################################################
### Prediction of growth rate phenotype from genotype
###################################################

genos<-list(bingeno,totcopgenos,fullgenos)
dtlist<-list(DATbinary,DATtotalcopies,DATfullgenotypes)
phenotype<-'growth.rate'
backg<-mean(qDT[copies.GAL4==0]$growth.rate)
backgsd<-sd(qDT[copies.GAL4==0]$growth.rate)
i<-1
nameFits<-data.table(order=1:6,interaction_order=c('first','second','third','fourth','fifth','sixth'))
genotype_ensembles<-c('binary','total copies','full')
predictions<-data.table(ldply(lapply(1:3,function(i){
	genocur<-copy(genos[[i]])
	DATcur<-copy(dtlist[[i]])
	DATcur0<-copy(DATcur[,c('plasgeno',genocur,'x'),with=F])	
	DATcur1<-DATcur0[,summary.func(x,'observed'),by=c('plasgeno',genocur)]
	
	lDAT0<-logconv2(DATcur0,x='x',backg=backg,backgsd=backgsd,genos=genocur)
	preds<-lapply(1:length(genocur),function(i){
		ord<-nameFits$order[i]
		nameFit<-nameFits$interaction_order[i]
		form<-paste0('lx~.')
		if(ord>1)form<-paste0('lx~.^',ord)
		expMod(mod=lm(as.formula(form),data=lDAT0),DT=DATcur,backg=backg,nameFit=nameFit,genos=genocur)
	})
	PREDmerge<-merge_recurse(preds,c(genocur))
	mDT<-melt(PREDmerge,id=genocur)
	mDT[,c('interaction_order','stat'):=colsplit(variable,'__',c('interaction_order','stat'))][,variable:=NULL]
	castform<-paste0(paste0(c(genocur,'interaction_order'),collapse='+'),'~stat')
	cDT<-dcast(mDT,as.formula(castform),value.var='value')
	DATmerge<-merge(cDT, DATcur1[,c('plasgeno',genocur,grepincols(DATcur1,c('_mean','_se','_upper','_lower'))),with=F],by=c(genocur),allow.cartesian=T)
#	logconv2(copy(DATcur1),x='observed_mean',backg=backg,genos=genocur,flip=T)
	DATmerge[,epistasis:=observed_mean-fit_mean]
	DATmerge[,genotype_ensemble:= genotype_ensembles[i]]
	DATmerge[,varexplained:=round(varexplained(observed_mean,fit_mean),3),by='interaction_order']
	DATmerge[,rsq:=round(summary(lm(observed_mean~fit_mean))$r.squared,3),by='interaction_order']
})))

predictions[,genotype_ensemble:=factor(genotype_ensemble,levels=genotype_ensembles)]
predictions[, interaction_order:=factor(interaction_order,levels= nameFits$interaction_order)]



p1<-ggplot(predictions,aes(fit_mean,observed_mean,xmax=fit_upper,xmin=fit_lower,ymax=observed_upper,ymin=observed_lower))+geom_errorbar()+geom_errorbarh()+
	geom_point(aes(colour=interaction_order),shape=21)+xlab(bquote('predicted growth rate'~h^-1))+ylab(bquote('observed growth rate'~h^-1))+
	geom_abline(col='grey70')+facet_grid(interaction_order~ genotype_ensemble,scales='free')+ theme(legend.position='none')+
	theme_minimal()+ggtitle('Observed vs predicted across all genotype ensembles and levels of interaction')
	
p2<-ggplot(predictions,aes(epistasis,colour= interaction_order))+stat_ecdf(alpha=0.5,size=1)+ylab('cumulative proportion')+ xlab(expression(epsilon))+
	facet_grid(~genotype_ensemble)+theme_minimal()+ggtitle('Distribution of epistasis')
ggplot(predictions,aes(as.numeric(interaction_order),rsq))+geom_line()+geom_point()+ylab('R-squared')+ xlab('interaction order')+
	facet_grid(~genotype_ensemble)+theme_minimal()
do.call(grid.arrange,list(p1,p2))

########
### proportion variance explained from predictions
########
predsplit<-split(predictions,by=c('interaction_order','genotype_ensemble'))
predsplit1<-predsplit[unlist(lapply(predsplit,function(x)nrow(x)>1))]
DTx<-predsplit1[[1]]
models<-data.table(ldply(lapply(predsplit1,function(DTx){
	DTx0<-merge(DATfullgenotypes,DTx[,c('plasgeno','fit_mean'),with=F],by='plasgeno')
	mod<-lm(x~fit_mean, DTx0)
	summod<-summary(mod)
	coefdf<-as.data.frame(summod$coefficients)
	coefDT<-data.table(coefdf, coef=rownames(coefdf),coef1=gsub('2','',gsub('1','',rownames(coefdf))))
	coefDT[,analysis:='lm_coef']
	coefDT[,varexp2:=NA]
	aovmod<-anova(mod)
	aovdf<-as.data.frame(aovmod)
	aovDT<-data.table(aovdf,coef=rownames(aovdf))
	#aovDT[,genotype_ensemble:= genotype_ensembles[i]]
	aovDT[,analysis:='anova']	
	aovDT[,coef1:=coef]
	aovDT[,varexp2:=varexplained(DTx0$x, DTx0$fit_mean)]
	meltID<-c('coef','coef1','analysis')
	out<-rbind(melt(coefDT,meltID),melt(aovDT,meltID))
	setnames(out,'variable','stat')
	out[,coef:=gsub('\\(','',coef)]
	out[,coef:=gsub('\\)','',coef)]
	out[,coef:=gsub('Intercept','int',coef)]
	out[,coef:=gsub('Residuals','resid',coef)]
	out[,interaction_order:=unique(DTx$interaction_order)]
	out[,genotype_ensemble:=unique(DTx$genotype_ensemble)]
	out	
})))

models[,.id:=NULL]

models[,allele_effect:=!coef%in%c('resid','int')]
models[allele_effect==T,coef_order:=as.character(unlist(sapply(coef,function(x)length(strsplit(x,'\\:')[[1]]))))]
models[coef =='int',coef_order:='int']
models[coef =='resid',coef_order:='resid']
coef_levels<-c('resid',6:1,'int')
models[,coef_order:=factor(coef_order,levels=coef_levels)]
models[,genotype_ensemble:=factor(genotype_ensemble,levels=genotype_ensembles)]


ssq<-models[analysis=='anova'&stat=='Sum Sq']

ssq[,varexp:=value/sum(value),by=c('genotype_ensemble','interaction_order')]
tech_var<-ssq[,lapply(.SD,length),by=c('genotype_ensemble','interaction_order')]
tech_var[,coef:=as.character(coef)]
tech_var[,coef:='technical']
tech_var[,coef_order:=as.character(coef)]
tech_var[,coef_order:='technical']

tech_var[, varexp:=as.numeric(value)]
tech_var[, varexp:=min(1-ssq$varexp)]
ssq1<-rbind(ssq,tech_var)
ssq1[coef=='resid', varexp:= varexp-min(1-ssq$varexp)]
varframe<-data.table(coef_order=unique(ssq1$coef_order),`variation source`=factor(c('genetic','residual','technical'),levels=c('technical','residual','genetic')))
ssq2<-merge(ssq1, varframe,by='coef_order')
p1<-ggplot(ssq2,aes(x=factor(interaction_order),y=varexp,fill= `variation source`))+geom_bar(stat='identity')+facet_grid(~genotype_ensemble,scales='free')+
	scale_fill_manual(values=c('grey50','cornflowerblue','indianred'))+geom_hline(yintercept=1)+
	ylab('proportion variance explained')+xlab('interaction order considered in linear model')+theme(axis.text.x=element_text(angle=30,hjust=1))
p1+	geom_text(aes(x=factor(interaction_order),y=1-varexp,label = paste0(round(varexp,2)*100,"%")), size=4)
h<-3;w<-h*2.5
sample(letters,4)
ggsave(paste0('.',figout,'180627-pPTIG-VarianceExplainedDifferentGeneticInteractionsLMs.png'),p1,height=h,width=w)


ssq_1<-na.exclude(models[stat=='varexp2'])
ssq_2<-na.exclude(models[stat=='varexp2'])
ssq_3<-na.exclude(models[stat=='varexp2'])
ssq_2[,value:=1-value-min(1-ssq$varexp]
ssq_3[,value:=min(1-ssq$varexp)]
ssq_1[,`variation source`:='genetic']
ssq_2[,`variation source`:='residual']
ssq_3[,`variation source`:='technical']
ssq_0<-rbind(ssq_1,ssq_2,ssq_3)
p11<-ggplot(ssq_0,aes(x=factor(interaction_order),y=value,fill= `variation source`))+geom_bar(stat='identity')+facet_grid(~genotype_ensemble,scales='free')+
	scale_fill_manual(values=c('grey50','cornflowerblue','indianred'))+geom_hline(yintercept=1)+
	ylab('proportion variance explained')+xlab('interaction order considered in linear model')+theme(axis.text.x=element_text(angle=30,hjust=1))
p11 +	geom_text(aes(x=factor(interaction_order),y=1-value,label = paste0(round(value,2)*100,"%")), size=4)


########
### pathway epistasis case for genotypic ensemble binary
########
mergecols<-c('plasgeno','g3g80',letters[1:3])
DATcur0 <-copy(DATbinary)
DATcur0[,g3g80:=factor(apply(data.frame(a,b,c),1,function(x)as.numeric(x[1]==0&x[2]==1&x[3]==1)))]
lDAT<-logconv2(DATcur0,x='x',backg=backg,backgsd=backgsd,genos= letters[1:3])
mod<-lm(lx~a+b+b*c,lDAT)
predDT<-expMod(DATcur0,mod=mod,backg=backg,genos=letters[1:3])
predMod<-merge(predDT, DATcur0,by=letters[1:3])
mod2<-lm(x~fit_mean, predMod)
ssq<-as.data.frame(anova(mod2))$`Sum Sq`
((ssq/sum(ssq))[1])
varexplained(predMod$x,predMod$fit_mean)
varexplained(predMod$x,predMod$fit_mean*mod2$coefficients[2]+mod2$coefficients[1])

#qplot(predMod$x,predMod$fit_mean)+coord_flip()
DTm<-merge(predMod[,summary.func(x,'observed'),by= mergecols], predDT,by=letters[1:3])
DTg3g80<-copy(DTm)
ggplot(DTm,aes(fit_mean,observed_mean,xmax=fit_upper,xmin=fit_lower,ymax=observed_upper,ymin=observed_lower,col=g3g80))+geom_errorbar(col='grey70')+geom_errorbarh(col='grey70')+geom_point()+
	geom_abline(col='grey70')

DATcur0 <-copy(DATbinary)
DATcur0[,g3g80:=factor(apply(data.frame(a,b,c),1,function(x)as.numeric(x[1]==0&x[2]==1&x[3]==1)))]
lDAT<-logconv2(DATcur0,x='x',backg=backg,backgsd=backgsd,genos= letters[1:3])
mod<-lm(lx~a+b+c,lDAT)
predDT<-expMod(DATcur0,mod=mod,backg=backg,genos=letters[1:3])
predMod<-merge(predDT, DATcur0,by=letters[1:3])
varexplained(predMod$x,predMod$fit_mean)
mod2<-lm(x~fit_mean, predMod)
ssq<-as.data.frame(anova(mod2))$`Sum Sq`
((ssq/sum(ssq))[1])
#qplot(predMod$x,predMod$fit_mean)+coord_flip()
DTm<-merge(predMod[,summary.func(x,'observed'),by=mergecols], predDT,by=letters[1:3])
DTsingles<-copy(DTm)
ggplot(DTm,aes(fit_mean,observed_mean,xmax=fit_upper,xmin=fit_lower,ymax=observed_upper,ymin=observed_lower,col=g3g80))+geom_errorbar(col='grey70')+geom_errorbarh(col='grey70')+geom_point()+
	geom_abline(col='grey70')

mDT<-melt(merge(DTsingles,DTg3g80,by= mergecols,suffixes=c('.singles','.singles + GAL3:GAL80 epistatic term')),id= mergecols)
mDT[,c('stat','mod'):=colsplit(variable,'\\.',c('stat','var'))][,variable:=NULL]
cDT<-cast(mDT,formula(plasgeno+g3g80+mod~stat),value.var='value')
paste0(sample(LETTERS,4),collapse='')
pKCGL <- ggplot(cDT,aes(fit_mean,observed_mean,xmax=fit_upper,xmin=fit_lower,ymax=observed_upper,ymin=observed_lower,col=g3g80))+geom_errorbar(col='grey70')+geom_errorbarh(col='grey70')+geom_point()+
	geom_abline(col='grey70')+facet_wrap(~mod)+scale_colour_manual(values=c('cornflowerblue','indianred'))+
	ylab(bquote('growth rate across all genotypes'~h^-1))+	xlab(bquote('predicted growth rate based on binary ensemble'~h^-1))+
	NULL
w<-8;h<-w*0.5
ggsave(paste0('.',figout,'180627-pKCGL-SingleEffectsVsGAL80GAL3EpistaticEffect-BinaryEnsemble.png'),pKCGL,height=h,width=w)


########
### explain residual variation in WT
########

DATcur0 <-copy(DATbinary[,!'genotype_ensemble',with=F])
DATcur1<-copy(DATtotalcopies[,!'genotype_ensemble',with=F])
rncols<-c(LETTERS[1:3],letters[1:3],'genotypei')
setnames(DATcur1,rncols,paste0(rncols,'_TC'))
genosi<-paste0(letters[1:3],'_TC')
DTM<-merge(DATcur0,DATcur1,by=c('x','plasgeno'))
WT<-DTM[genotypei=='1 1 1']
lWT<-logconv(WT[,c('x',genosi),with=F],backg=backg)
lWT[,paste0(genosi):=lapply(.SD,as.factor),.SDcols=genosi]
mod1<-lm(lx~.,lWT)
mod2<-lm(lx~.^2,lWT)
mod3<-lm(lx~.^3,lWT)
em1<-expMod(WT,mod1,backg=backg,genos=genosi,nameFit='first order')
em2<-expMod(WT,mod2,backg=backg,genos=genosi,nameFit='second order')
em3<-expMod(WT,mod3,backg=backg,genos=genosi,nameFit='third order')
DTm<-merge_recurse(list(em1,em2,em3),genosi)
mDT<-melt(DTm,id=genosi)
mDT[,c('order','stat'):=colsplit(variable,'__',c('order','stat'))];mDT[,variable:=NULL]
cDT<-cast(mDT,formula(a_TC+b_TC+c_TC + order ~ stat),value.var='value')

WTsumm<-WT[,summary.func(x),by=c(genosi,'plasgeno')]
WTm<-data.table(merge(cDT,merge(WT,WTsumm,by=c(genosi,'plasgeno')),by=genosi))
WTm[,varexp_full:=varexplained(x,mean),by='order']
Vg<-unique(WTm$varexp_full)
WTm[,varexp_tc:=varexplained(x,fit_mean),by='order']
WTm[,Vg_Cn:=varexp_tc/Vg]
WTm[,order:=factor(gsub('\\.','\\ ',order),levels=c('first order','second order','third order'))]
WTm[,facet:=applyPaste(data.frame(order,round(varexp_full,2),round(Vg_Cn,2)),collapse='\n')]
WTm[,rel_ActRep:=(apply(data.frame(a_TC,b_TC,c_TC),1,function(xa){
	x<-as.numeric(xa)
	3-(1+(x[1]+x[2]-x[3]))
}))]


WTm[,`Activators - Repressors`:=factor(rel_ActRep,levels=c(0,1,2,3))]
pIWTG_ALT<-ggplot(WTm,aes(fit_mean,mean,ymax=upper,ymin=lower,xmax=fit_upper,xmin=fit_lower))+geom_abline(col='grey70')+geom_errorbar()+geom_errorbarh()+geom_point(aes(col=`Activators - Repressors`))+
	ylab(expression(paste('all measurements of binary [111] ',mu~h^-1)))+xlab(expression(paste('prediction lm (',mu,' ~ GAL4 + GAL3 + GAL80 copy number)')))+
	facet_grid(~facet)
w<-10;h<-w*0.4

ggsave(paste0('.',figout,'180702-pIWTG_ALT-EffectOfCopyNumberScatterplot-TotalCopyEnsemble.png'), pIWTG_ALT,height=h,width=w)


########
### explain residual variation in WT: activators - repressors
########
qDT[,wt:= factor(apply(colsplit(binarycount,'\\ ',c('a','b','c')),1,function(xa){
	x<-as.numeric(xa)
	as.logical(x[1]>0&x[2]>0&x[3]>0)}))]

DATcur0 <-copy(DATbinary)
DATcur0[,g3g80:=factor(apply(data.frame(a,b,c),1,function(x)as.numeric(x[1]==0&x[2]==1&x[3]==1)))]

DATcur1 <-copy(DATtotalcopies)
DATcur1[,wt:= factor(apply(colsplit(genotypei,'\\ ',c('a','b','c')),1,function(xa){
	x<-as.numeric(xa)
	as.logical(x[1]>0&x[2]>0&x[3]>0)}))]
DATcur1[wt==T,g3g80_rel:=factor(apply(data.frame(a,b,c),1,function(xa){
	x<-as.numeric(xa)
	x[1]+x[2]-x[3]
	}))]
DATcur1[is.na(g3g80_rel), g3g80_rel:=factor('G4delta')]
DATcur1[,c('a1','b1','c1'):=list(a,b,c)]
DATcur1[,c('a','b','c'):=factor(apply(data.frame(a1,b1,c1),1,function(xa){
	x<-as.numeric(xa)
	as.numeric(x[1]>0&x[2]>0&x[3]>0)}))]

genos<-c(letters[1:3],'g3g80_rel','wt')
mergecols<-genos
lDAT<-logconv2(DATcur1,x='x',backg=backg,backgsd=backgsd,genos=genos)
mod<-lm(lx~g3g80_rel,droplevels(lDAT[wt==T]))
predDT<-expMod(DATcur1[wt==T],mod=mod,backg=backg,genos=c(LETTERS[1:3],genos))
predMod<-merge(predDT, DATcur1,by=genos)
summary(lm(x~ g3g80_rel,predMod))
summ1<-DATcur1[wt==T,summary.func(x,'growth.rate'),by=c(LETTERS[1:3],'plasgeno')]
summ2<-qDT[,summary.func(fracon.gal,'fracon.gal'),by=c('plasgeno')]
predum<-qDT[wt==T,c('plasgeno','fracon.gal','growth.rate'),with=F]
setnames(predum,'growth.rate','x')
summDT<-merge(merge(DATcur1[wt==T],merge(summ1,summ2,by='plasgeno'),by=c(LETTERS[1:3],'plasgeno'),all.x=T), predum,by=c('x','plasgeno'),all.x=T)
a<-varexplained(predMod$x,predMod$fit_mean)
b<-varexplained(summDT$x, summDT $mean)
a/b
DTm<-predMod[,summary.func(x),by= g3g80_rel]
DTsplit<-split(data.table(g3g80_rel=predMod$g3g80_rel,genotypei=predMod$genotypei),by='g3g80_rel')
x<-DTsplit[[2]]
labelTab<-data.table(ldply(lapply(DTsplit,function(x){
	out<-data.table(genotypei=unique(x$genotypei))
	a<-as.numeric(unique(x$g3g80_rel))
	b<-4-a
	out[, g3g80_rel:= unique(x$g3g80_rel)]
	out[, `Activators - Repressors`:=paste0(b,'\n\n',paste0(paste0('[',unique(x$genotypei),']',collapse='\n'),collapse='\n'),collapse='')]
	out
})))
labelTab[,.id:=NULL]
levs<-unique(labelTab$`Activators - Repressors`)[4:1]
labelTab[,`Activators - Repressors`:=factor(`Activators - Repressors`,levels=c(levs))]


DTm<-merge(summDT,labelTab,by=c('g3g80_rel','genotypei'))
DTm[,`total copy genotype`:=genotypei]
mDT<-melt(DTm[,c('g3g80_rel','genotypei','total copy genotype','x','plasgeno','fracon.gal','Activators - Repressors'),with=F],id=c('total copy genotype','g3g80_rel','genotypei','plasgeno','Activators - Repressors'))
mDT[,variable:=gsub('x','growth rate',variable)]
mDT[,variable:=gsub('fracon.gal','fraction ON in galactose',variable)]
paste0(sample(LETTERS,4),collapse='')
pIWTG<-ggplot(mDT,aes(`Activators - Repressors`,value))+geom_boxplot()+geom_jitter(width=0.1,aes(col=`total copy genotype`),alpha=0.7,shape=21)+geom_smooth(method='lm')+
	ylab('all measurements of binary [111]')+xlab('Activators - Repressors')+facet_grid(variable~.,scales='free',switch='y')
w<-7;h<-w
ggsave(paste0('.',figout,'180627-pIWTG-EffectOfCopyNumberBoxplot-TotalCopyEnsemble.png'), pIWTG,height=h,width=w)



########
### jacknifing above for "cross validation"
########

jacknife<-function(frac){
	pred<-data.table(ldply(lapply(1:3,function(i){
		genocur<-copy(genos[[i]])
		DATcur<-copy(dtlist[[i]])
		DATcur[,arb:=factor(1:nrow(dtlist[[i]]))]
		DATcur0_train<-sample_frac(copy(DATcur[,c('plasgeno',genocur,'x','arb'),with=F]),replace=F,size= frac)
		DATcur0_test<-DATcur[!DATcur$arb%in%DATcur0_train$arb, summary.func(x,'observed'),by=c('plasgeno',genocur)]
		
		lDAT0_train<-logconv2(DATcur0_train,x='x',backg=backg,backgsd=backgsd,genos=genocur)
		preds<-lapply(1:length(genocur),function(i){
			ord<-nameFits$order[i]
			nameFit<-nameFits$interaction_order[i]
			form<-paste0('lx~.')
			if(ord>1)form<-paste0('lx~.^',ord)
			expMod(mod=lm(as.formula(form),data= lDAT0_train),DT=DATcur,backg=backg,nameFit=nameFit,genos=genocur)
		})
		PREDmerge<-merge_recurse(preds,c(genocur))
		mDT<-melt(PREDmerge,id=genocur)
		mDT[,c('interaction_order','stat'):=colsplit(variable,'__',c('interaction_order','stat'))][,variable:=NULL]
		castform<-paste0(paste0(c(genocur,'interaction_order'),collapse='+'),'~stat')
		cDT<-dcast(mDT,as.formula(castform),value.var='value')
		DATmerge<-merge(cDT, DATcur0_test[,c('plasgeno',genocur,grepincols(DATcur1,c('_mean','_se','_upper','_lower'))),with=F],by=c(genocur),allow.cartesian=T)
	#	logconv2(copy(DATcur1),x='observed_mean',backg=backg,genos=genocur,flip=T)
		DATmerge[,epistasis:=observed_mean-fit_mean]
		DATmerge[,genotype_ensemble:= genotype_ensembles[i]]
		DATmerge[,varexplained:=round(varexplained(observed_mean,fit_mean),3),by='interaction_order']
		DATmerge[,rsq:=round(summary(lm(observed_mean~fit_mean))$r.squared,3),by='interaction_order']
	})))
	c(pred[,unique(rsq),by=mergecols]$V1,pred[,unique(varexplained),by=mergecols]$V1)
}


predTemplate<-data.table(predictions[,unique(rsq),by=mergecols][,mergecols,with=F]);predTemplate<-rbind(predTemplate,predTemplate)
predTemplate[,stat:=c(rep('rsq',12),rep('varexp',12))]
fracs<-seq(0.1,0.9,by=0.1)

# # jacknifing here -- takes about 10 min
# jacknifed<-lapply(fracs,function(frac){
	# a<-data.table(replicate(100,jacknife(frac))	)
	# a[,frac:=1-frac]
	# a
	# })

# save(jacknifed,file=paste0('.',datout,'180625-JacknifeCrossValidation-100Reps.rData'))

# just load waht you ran before
load(paste0('.',datout,'180625-JacknifeCrossValidation-100Reps.rData'))
jacknifed

jckSumm1<-data.table(ldply(lapply(jacknifed,function(DT){
	a<-data.table(predTemplate,t(apply(DT,1,function(x)c(mean(x,na.rm=T),sd(x,na.rm=T)))),frac=DT$frac)
	setnames(a,c('V1','V2'),c('value_mean','value_sd'))
	})))
jckSumm1[,value_upper:=value_mean+value_sd]
jckSumm1[,value_lower:=value_mean-value_sd]
jckSumm1[value_upper>1,value_upper:=1]
jckSumm1[value_lower<0,value_lower:=0]

jckSumm2<-data.table(ldply(lapply(jacknifed,function(DT){
	a<-data.table(predTemplate,DT)
	melt(a,id=c('genotype_ensemble','interaction_order','stat','frac'))
	})))
jckSumm2[,val_cor:=value]
jckSumm2[val_cor<0,val_cor:=0]
p3 <-ggplot(jckSumm2[stat=='rsq'],aes(factor(frac),value))+geom_boxplot()+facet_grid(interaction_order~ genotype_ensemble)+ylim(c(0,1))+theme(axis.text.x=element_text(angle=30,hjust=1))+
	xlab('proportion left out from training')+ylab('R-squared for predicted vs observed')+theme_minimal()+ggtitle('N-fold cross validation: R-squared obs vs predicted')
p4<-ggplot(jckSumm2[stat=='varexp'],aes(factor(frac), val_cor))+geom_boxplot()+ facet_grid(interaction_order~ genotype_ensemble)+ylim(c(0,1))+theme(axis.text.x=element_text(angle=30,hjust=1))+
	xlab('proportion left out from training')+ylab('proportion variance explained for predicted vs observed')+theme_minimal()+ggtitle('N-fold cross validation: proportion variance explained obs vs predicted')
p5<-ggplot(jckSumm1[stat=='rsq'],aes(frac,value_mean,ymax=value_upper,ymin=value_lower))+geom_point(col='grey70')+geom_errorbar(col='grey70')+geom_line()+
	facet_grid(interaction_order~ genotype_ensemble)+geom_hline(yintercept=1,col='indianred',linetype='dotted')+geom_hline(yintercept=0,col='cornflowerblue',linetype='dotted')+
	xlab('proportion left out from training')+ylab('R-squared for predicted vs observed')+theme_minimal()+ggtitle('N-fold cross validation: R-squared obs vs predicted')
sample(letters,3)
newdir<-paste0('.',figout,'1707xx-17092x-180625-pQEI-JacknifeAllEnsemblesAndLevels')
if(!dir.exists(newdir))system(paste0('mkdir ', newdir))
setwd(newdir)
w<-12;h<-w*0.7
whs<-c(1,0.5,0.7,0.7)
plotlist<-list(p1,p2,p3,p4,p5);i<-1
lapply(seq(length(plotlist)),function(i){
	ggsave(paste0('180625-pQEI-Jacknifing-plot-',i,'.png'), plotlist[[i]],height=h*whs[i],width=w*whs[i])
})
setwd(head.dir)



###################################################
### proportion variance explained: just looking at log values
###################################################

logconv<-function(DT,x='x',backg=0){
	DT<-data.table(data.frame(DT))
	DT[,lx:=log(abs(x-backg))]
	data.table(data.frame(DT[,!'x',with=F]))
}
x<-'observed_mean'
DT<-copy(DATcur1)

logconv2<-function(DT,x='x',backg=0,backgsd =0,genos,flip=F){
	if(x!='x')setnames(DT,x,'x')
	DT<-data.table(data.frame(DT))
	credmin<-backg+backsd
	DT[,corx:=x]
	DT[corx<credmin,corx:=backg+10^-16]
	DT[,subx:=abs(x-backg)]
	DT[,flip:=(x-backg)/abs(x-backg)]
	DT[,lx:=log2(subx)]
	DT[,lx:=log2(abs(x-backg))]
	a<-data.table(data.frame(DT[,c('lx',genos),with=F]))
	if(flip==T){a<-data.table(data.frame(DT[,c('lx','flip',genos),with=F]))}
	a
}
VarExpMod<-function(mod){
	
	DT<-data.table(data.frame(DT[,mean(x,na.rm=T),by=genos][,genos,with=F]))
	pred<-predict(mod,DT,se.fit=T)
	DT[,c('fit_mean','fit_se'):=list(exp(pred$fit)+backg,exp(pred$fit)*pred$se.fit)]
	DT[,c('fit_upper','fit_lower'):=list(fit_mean+fit_se,fit_mean-fit_se)]
	if(!is.null(nameFit))setnames(DT,c('fit_mean','fit_se','fit_upper','fit_lower'),paste0(nameFit,'__',c('fit_mean','fit_se','fit_upper','fit_lower')))
	data.table(data.frame(DT))
}

genos<-list(bingeno,totcopgenos,fullgenos)
dtlist<-list(DATbinary,DATtotalcopies,DATfullgenotypes)
phenotype<-'growth.rate'
backg<-mean(qDT[copies.GAL4==0]$growth.rate)
backgsd<-sd(qDT[copies.GAL4==0]$growth.rate)
i<-1
nameFits<-data.table(order=1:6,interaction_order=c('first','second','third','fourth','fifth','sixth'))
genotype_ensembles<-c('binary','total copies','full')
models<-data.table(ldply(lapply(1:3,function(i){
	genocur<-copy(genos[[i]])
	DATcur<-copy(dtlist[[i]])
	DATcur0<-copy(DATcur[,c('plasgeno',genocur,'x'),with=F])	
	DATcur1<-DATcur0[,summary.func(x,'observed'),by=c('plasgeno',genocur)]
	
	lDAT0<-logconv2(DATcur0,x='x',backg=backg,backgsd=backgsd,genos=genocur)
	coefs<-data.table(ldply(lapply(1:length(genocur),function(i){
		ord<-nameFits$order[i]
		nameFit<-nameFits$interaction_order[i]
		form<-paste0('lx~.')
		if(ord>1)form<-paste0('lx~.^',ord)
		mod<-lm(as.formula(form),data=lDAT0)
		summod<-summary(mod)
		coefdf<-as.data.frame(summod$coefficients)
		coefDT<-data.table(coefdf, coef=rownames(coefdf),coef1=gsub('2','',gsub('1','',rownames(coefdf))))
		coefDT[, interaction_order:= ord]	
		coefDT[,analysis:='lm_coef']	
		aovmod<-anova(mod)
		aovdf<-as.data.frame(aovmod)
		aovDT<-data.table(aovdf,coef=rownames(aovdf))
		#aovDT[,genotype_ensemble:= genotype_ensembles[i]]
		aovDT[, interaction_order:= ord]
		aovDT[,analysis:='anova']	
		aovDT[,coef1:=coef]
		meltID<-c('coef','coef1','interaction_order','analysis')
		out<-rbind(melt(coefDT,meltID),melt(aovDT,meltID))
		setnames(out,'variable','stat')
		out[,coef:=gsub('\\(','',coef)]
		out[,coef:=gsub('\\)','',coef)]
		out[,coef:=gsub('Intercept','int',coef)]
		out[,coef:=gsub('Residuals','resid',coef)]
		out
	})))
	coefs[,genotype_ensemble:= genotype_ensembles[i]]
	coefs[order(genotype_ensemble,interaction_order,analysis,coef)]
})))
models[,allele_effect:=!coef%in%c('resid','int')]
models[allele_effect==T,coef_order:=as.character(unlist(sapply(coef,function(x)length(strsplit(x,'\\:')[[1]]))))]
models[coef =='int',coef_order:='int']
models[coef =='resid',coef_order:='resid']
coef_levels<-c('resid',6:1,'int')
models[,coef_order:=factor(coef_order,levels=coef_levels)]
models[,genotype_ensemble:=factor(genotype_ensemble,levels=genotype_ensembles)]

ssq<-models[analysis=='anova'&stat=='Sum Sq']
ssq[,varexp:=value/sum(value),by=c('genotype_ensemble','interaction_order')]

ggplot(ssq,aes(x=factor(interaction_order),y=varexp,fill= coef_order))+geom_bar(stat='identity')+facet_grid(~genotype_ensemble,scales='free')


lm_coefs<-models[analysis=='lm_coef']
cDT<-dcast(lm_coefs,formula(coef+interaction_order+analysis+ genotype_ensemble+ coef_order~stat),value.var='value')
coef_levels<-cDT[,unique(analysis),by=c('coef','coef_order')][order(coef_order,decreasing=F)]$coef
cDT[,coef_levels:=factor(coef,levels=coef_levels)]
cDT[,sig:=`Pr(>|t|)`<0.05]
ggplot(cDT,aes(coef_levels,Estimate))+geom_point(aes(colour=sig))+facet_grid(interaction_order~genotype_ensemble)+coord_flip()

cDT[interaction_order=='2'&genotype_ensemble=='full'][order(`Pr(>|t|)`)]


###################################################
### Mapping phenotype onto radial fitness landscapes
###################################################
# custom function for making the labels
# make a label data frame
labeldfun<-function(DT,rotation_factor=0){
	mergesumm<-copy(DT)
	logiframe<-as.matrix(mergesumm[,grepincols(mergesumm,'start'),with=F][,lapply(.SD,function(x)as.numeric(as.logical(x+10)))][,-1])
	labels<-melt(data.table(mergesumm[,lapply(.SD,as.numeric),.SDcols=c(LETTERS[1:3])]*logiframe,arb=factor(1:nrow(mergesumm))),id='arb')[,!c('arb','variable')]
	setnames(labels,'value','label')
	rs<-melt(data.table(mergesumm[,grepincols(mergesumm,'r_'),with=F][,-1]*((logiframe)), arb=factor(1:nrow(mergesumm))),id='arb')[,!c('arb','variable')]
	setnames(rs,'value','r')
	thetas<-melt(data.table(mergesumm[,grepincols(mergesumm,'start_'),with=F][,-1]*logiframe, arb=factor(1:nrow(mergesumm))),id='arb')[,!c('arb')]
	setnames(thetas,'value','theta')
	genotypes<-(melt(data.table(mergesumm[,'genotypei',with=F],logiframe, arb=factor(1:nrow(mergesumm)))))
	labelDT<-na.exclude(data.table(labels,rs,thetas,genotypes[,c('arb','genotypei'),with=F]))
	labelDT[,mod_theta1:=(1/(.N)*pi),by='variable']
	lsplit<-split(labelDT,by='variable')
	levs<-1:length(lsplit)
	df<-ldply(lapply(levs,function(i){
		a<-0
		if(i>1){	a<-unique(lsplit[[i]]$mod_theta1)}
		if(i>2){	a<-unique(lsplit[[i]]$mod_theta1)+ unique(lsplit[[i-1]]$mod_theta1)}
		a
	}))
	df$variable<-paste0('start_',(1:nrow(df)+1))
	labelDT1<-merge(labelDT,df,by='variable',all=T)
	labelDT1[,x:=(r-0.5)*cos(-(theta)+V1+ rotation_factor)]
	labelDT1[,y:=(r-0.5)*sin(-(theta)+V1+ rotation_factor)]	
	return(data.frame(labelDT1))
}


library(ggforce)
pRAD <- ggplot() + theme_no_axes() + coord_fixed()

# get data together

genolist<-list(bingeno,totcopgenos,fullgenos)
dtlist<-list(DATbinary,DATtotalcopies,DATfullgenotypes)
phenotype<-'growth.rate'
backg<-mean(qDT[copies.GAL4==0]$growth.rate)
backgsd<-sd(qDT[copies.GAL4==0]$growth.rate)
i<-1
nameFits<-data.table(order=1:6,interaction_order=c('first','second','third','fourth','fifth','sixth'))
genotype_ensembles<-c('binary','total copies','full')
predictions[,genotype_ensemble:=as.character(genotype_ensemble)]
bin<-dtlist[[1]][,unique(a),by=c(letters[1:3],LETTERS[1:3])][,!'V1']
tot<-dtlist[[2]][,unique(a),by=c(letters[1:3],LETTERS[1:3])][,!'V1']
full<-dtlist[[3]][,unique(a),by=c(letters[1:6],LETTERS[1:6])][,!'V1']
bin[,genotype_ensemble:=as.character('binary')]
tot[,genotype_ensemble:=as.character('total copies')]
full[,genotype_ensemble:=as.character('full')]
bin[,genotypei:=applyPaste(data.frame(A,B,C),' ')]
tot[,genotype_ensemble:=as.character('total copies')]
tot[,genotypei:=applyPaste(data.frame(A,B,C),' ')]
full[,genotype_ensemble:=as.character('full')]
full[,genotypei:=applyPaste(data.frame(A,B,C,D,E,F),' ')]

binpred<-merge(predictions,bin,by=c(letters[1:3],'genotype_ensemble'))
totpred<-merge(predictions,tot,by=c(letters[1:3],'genotype_ensemble'))
fullpred<-merge(predictions,full,by=c(letters[1:6],'genotype_ensemble'))

keepcols<-c('genotypei',LETTERS[1:3],)

# make plots for binary
DT<-copy(binpred[interaction_order=='first',mean(observed_mean),by=c('genotypei',LETTERS[1:3],letters[1:3])])[order(a,b,c)]
setnames(DT,'V1','mean')
genos<-c('A','B','C')
mergesumm<-radial_plot_space(DT,genos)[order(a,b,c)]
dummysumm<-data.table(label1=c('WT','GAL4 binary','GAL3 binary','GAL80 binary'))
dummysumm[,c('x1','y1'):=list(c(0,0,0,0),c(0,1.7,2.7,3.7))]
rotation_factor<-0

labelDT1<-labeldfun(mergesumm,rotation_factor=rotation_factor)


phenolim<-c(0.035,0.22)
pBinary<-pRAD + 
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_4, start = start_4, end = end_4, fill = mean),col='grey80', data = mergesumm[r_4==4])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_3, start = start_3, end = end_3, fill = mean),col='grey80', data = mergesumm[r_3==3])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_2, start = start_2, end = end_2, fill = mean),col='grey80', data = mergesumm[r_2==2])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_1, start = start_1, end = end_1, fill = mean),col='grey80', data = mergesumm[r_1==1])+
	geom_circle(aes(x0=0,y0=0,r=1),col='grey80')+
	geom_text(data=dummysumm,aes(x=x1,y=y1+0.1,label= label1))+
	geom_text(data=labelDT1,aes(x=x,y=y,label= genotypei))+
	scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),name='observed growth rate',limits= phenolim)+theme(panel.border=element_blank())
pBinary

# make plots for total copies
DT<-copy(totpred[interaction_order=='first',mean(observed_mean),by=c('genotypei',LETTERS[1:3],letters[1:3])])[order(a,b,c)]
setnames(DT,'V1','mean')
genos<-c('A','B','C')
mergesumm<-radial_plot_space(DT,genos)[order(a,b,c)]
dummysumm<-data.table(label1=c('WT','GAL4 copies','GAL3 copies','GAL80 copies'))
dummysumm[,c('x1','y1'):=list(c(0,0,0,0),c(0,1.7,2.7,3.7))]
rotation_factor<-1/3.2*pi

labelDT1<-labeldfun(mergesumm,rotation_factor=rotation_factor)

pRAD <- ggplot() + theme_no_axes() + coord_fixed()

phenolim<-c(0.035,0.22)
pTotCop<-pRAD + 
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_4, start = start_4, end = end_4, fill = mean),col='grey80', data = mergesumm[r_4==4])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_3, start = start_3, end = end_3, fill = mean),col='grey80', data = mergesumm[r_3==3])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_2, start = start_2, end = end_2, fill = mean),col='grey80', data = mergesumm[r_2==2])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_1, start = start_1, end = end_1, fill = mean),col='grey80', data = mergesumm[r_1==1])+
	geom_circle(aes(x0=0,y0=0,r=1),col='grey80')+
	geom_text(data=dummysumm,aes(x=x1,y=y1+0.1,label= label1))+
	geom_text(data=labelDT1,aes(x=x,y=y,label= genotypei))+
	scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),name='observed growth rate',limits= phenolim)+theme(panel.border=element_blank())
pTotCop+ggtitle('phenotypic variation at total copy number ensemble')
dev.new()
# make plots for full genotype landscape
DT<-copy(fullpred[interaction_order=='first',mean(observed_mean),by=c('genotypei',LETTERS[1:6],letters[1:6])])[order(a,b,c,d,e,f)]
setnames(DT,'V1','mean')
genos<-c('A','B','C','D','E','F')
mergesumm<-radial_plot_space(DT,genos)[order(a,b,c,d,e,f)]
dummysumm<-data.table(label1=c('WT','GAL4 plasmid','GAL4 genome','GAL3 plasmid','GAL3 genome','GAL80 plasmid','GAL80 genome'))
dummysumm[,c('x1','y1'):=list(rep(0,7),c(0,seq(1:6)+0.2))]
rotation_factor<-0

#labelDT1<-labeldfun(mergesumm,rotation_factor=rotation_factor)


phenolim<-c(0.035,0.22)
pFullLandscape<-pRAD + 
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_7, start = start_7, end = end_7, fill = mean),col='grey80', data = mergesumm[r_7==7])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_6, start = start_6, end = end_6, fill = mean),col='grey80', data = mergesumm[r_6==6])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_5, start = start_5, end = end_5, fill = mean),col='grey80', data = mergesumm[r_5==5])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_4, start = start_4, end = end_4, fill = mean),col='grey80', data = mergesumm[r_4==4])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_3, start = start_3, end = end_3, fill = mean),col='grey80', data = mergesumm[r_3==3])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_2, start = start_2, end = end_2, fill = mean),col='grey80', data = mergesumm[r_2==2])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_1, start = start_1, end = end_1, fill = mean),col='grey80', data = mergesumm[r_1==1])+
	geom_circle(aes(x0=0,y0=0,r=1),col='grey80')+
	geom_text(data=dummysumm,aes(x=x1,y=y1+0.1,label= label1))+
#	geom_text(data=labelDT1,aes(x=x,y=y,label= genotypei))+
	scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),name='observed growth rate',limits= phenolim)+theme(panel.border=element_blank())
pFullLandscape


sample(LETTERS,4)
plotlist<-list(pBinary,pTotCop,pFullLandscape,ncol=3)
newdir<-paste0('.',figout,'180627-pWTZV-RadialInteractionPlots')
if(!dir.exists(newdir))system(paste0('mkdir ', newdir))
setwd(newdir)
w<-12;h<-w*0.7
whs<-c(0.5,1,3)
lapply(seq(length(plotlist)),function(i){
	ggsave(paste0('180627-RadialInteractionPlots-',i,'.png'), plotlist[[i]],height=h*whs[i],width=w*whs[i])
})
setwd(head.dir)


do.call(grid.arrange,)

###################################################
### Mapping how ensemble genotypes predict 
###################################################
# custom function for making the labels
# make a label data frame
labeldfun2<-function(DT,rotation_factor=0){
	mergesumm<-copy(DT)
	logiframe<-as.matrix(mergesumm[,grepincols(mergesumm,'start'),with=F][,lapply(.SD,function(x)as.numeric(as.logical(x+10)))][,-1])
	labels<-melt(data.table(mergesumm[,lapply(.SD,as.numeric),.SDcols=c(LETTERS[1:3])]*logiframe,arb=factor(1:nrow(mergesumm))),id='arb')[,!c('arb','variable')]
	setnames(labels,'value','label')
	rs<-melt(data.table(mergesumm[,grepincols(mergesumm,'r_'),with=F][,-1]*((logiframe)), arb=factor(1:nrow(mergesumm))),id='arb')[,!c('arb','variable')]
	setnames(rs,'value','r')
	thetas<-melt(data.table(mergesumm[,grepincols(mergesumm,'start_'),with=F][,-1]*logiframe, arb=factor(1:nrow(mergesumm))),id='arb')[,!c('arb')]
	setnames(thetas,'value','theta')
	sigs<-melt(data.table(mergesumm$sig*logiframe, arb=factor(1:nrow(mergesumm))),id='arb')[,!c('variable','arb')]
	setnames(sigs,'value','sig')
	genotypes<-(melt(data.table(mergesumm[,'genotypei',with=F],logiframe, arb=factor(1:nrow(mergesumm)))))
	labelDT<-(data.table(labels,rs,thetas,sigs,genotypes[,c('arb','genotypei'),with=F]))
	labelDT[,mod_theta1:=(1/(.N)*pi),by='variable']
	lsplit<-split(labelDT,by='variable')
	levs<-1:length(lsplit)
	df<-ldply(lapply(levs,function(i){
		a<-0
		if(i>1){	a<-unique(lsplit[[i]]$mod_theta1)}
		if(i>2){	a<-unique(lsplit[[i]]$mod_theta1)+ unique(lsplit[[i-1]]$mod_theta1)}
		a
	}))
	df$variable<-paste0('start_',(1:nrow(df)+1))
	labelDT1<-merge(labelDT,df,by='variable',all=T)
	labelDT1[,x:=(r-0.5)*cos(-(theta)+V1+ rotation_factor)]
	labelDT1[,y:=(r-0.5)*sin(-(theta)+V1+ rotation_factor)]	
	return(data.frame(labelDT1[!is.na(x)]))
}


library(ggforce)
pRAD <- ggplot() + theme_no_axes() + coord_fixed()

# get data together

genolist<-list(bingeno,totcopgenos,fullgenos)
dtlist<-list(DATbinary,DATtotalcopies,DATfullgenotypes)
phenotype<-'growth.rate'
backg<-mean(qDT[copies.GAL4==0]$growth.rate)
backgsd<-sd(qDT[copies.GAL4==0]$growth.rate)
i<-1
nameFits<-data.table(order=1:6,interaction_order=c('first','second','third','fourth','fifth','sixth'))
genotype_ensembles<-c('binary','total copies','full')
predictions[,genotype_ensemble:=as.character(genotype_ensemble)]
bin<-dtlist[[1]][,unique(a),by=c(letters[1:3],LETTERS[1:3])][,!'V1']
tot<-dtlist[[2]][,unique(a),by=c(letters[1:3],LETTERS[1:3])][,!'V1']
full<-dtlist[[3]][,unique(a),by=c(letters[1:6],LETTERS[1:6])][,!'V1']
bin[,genotype_ensemble:=as.character('binary')]
tot[,genotype_ensemble:=as.character('total copies')]
full[,genotype_ensemble:=as.character('full')]
bin[,genotypei:=applyPaste(data.frame(A,B,C),' ')]
tot[,genotype_ensemble:=as.character('total copies')]
tot[,genotypei:=applyPaste(data.frame(A,B,C),' ')]
full[,genotype_ensemble:=as.character('full')]
full[,genotypei:=applyPaste(data.frame(A,B,C,D,E,F),' ')]

bycols1<-c(letters[1:3],LETTERS[1:3],'genotypei')
bycols2<-c(letters[1:6],LETTERS[1:6],'genotypei')
dtlist<-list(DATbinary,DATtotalcopies,DATfullgenotypes)

DATbinary[,genotype_ensemble:=as.character('binary')]
DATtotalcopies
full2<-merge(dtlist[[3]],full,by=bycols2)

# get data together
DT1<-DATbinary[,summary.func(x,'binary'),by=bycols1][order(a,b,c)]
DT2<-DATtotalcopies[,summary.func(x,'totalcopy'),by=bycols1][order(a,b,c)]
DT3<-DATfullgenotypes[,summary.func(x,'full'),by=bycols2][order(a,b,c,d,e,f)]

DT1[,bin_tot:=genotypei]
DT1[,bin_full:=genotypei]
DT2[,bin_tot:=apply(data.frame(A,B,C),1,function(x){
	paste0(as.numeric(as.logical(as.numeric(x))),collapse=' ')
	})]
DT2[,tot_full:=genotypei]
attach(DT3);detach(DT3)
DT3[,tot_full:=apply(data.frame(A,B,C,D,E,F),1,function(x){
	xa<-as.numeric(x)
	paste(sum(xa[1:2]),sum(xa[3:4]),sum(xa[5:6]),collapse=' ')
})]
DT3[,bin_full:=apply(data.frame(A,B,C,D,E,F),1,function(x){
	xa<-as.numeric(x)
	paste(as.numeric(as.logical(sum(xa[1:2]))),as.numeric(as.logical(sum(xa[3:4]))),as.numeric(as.logical(sum(xa[5:6]))),collapse=' ')
	})]


# compare total copies to binary

bin_tot<-merge(DT1[,c('bin_tot',grepincols(DT1,'binary')),with=F],DT2,by='bin_tot')[order(a,b,c)]
bin_tot[,diff_mean:= totalcopy_mean-binary_mean]
bin_tot[,diff_se:=sqrt(totalcopy_se^2+binary_se^2)]
bin_tot[,diff_df:=binary_N+totalcopy_N-1]
attach(bin_tot);pval<-data.table(t(apply(data.frame(totalcopy_mean,binary_mean,totalcopy_sd,binary_sd,totalcopy_N,binary_N),1,function(x){
	test<-t.test2(m1=x[1],m2=x[2],s1=x[3],s2=x[4],n1=x[5],n2=x[6])
})))$p.value;detach(bin_tot)
bin_tot$diff_p.val<-pval
bin_tot[,p.adj:=p.adjust(diff_p.val)]
bin_tot[,sig:=p.adj<0.05]

bin_tot[,jitter_mean:=jitter(binary_mean,factor=20)]
bin_tot[,jitter_upper:=jitter_mean-binary_mean+binary_upper]
bin_tot[,jitter_lower:=jitter_mean-binary_mean+binary_lower]
bin_tot[,significant:=sig]

pBinVsTot1<-ggplot(bin_tot,aes(x= jitter_mean,xmax= jitter_upper,xmin= jitter_lower,y=totalcopy_mean,ymax=totalcopy_upper,ymin=totalcopy_lower,col= significant))+geom_abline()+geom_errorbar()+geom_errorbarh()+geom_point()+
	ylab(bquote('growth rate of total copy genotypes'~h^-1))+	xlab(bquote('jittered growth rate of binary genotypes'~h^-1))+
	scale_colour_manual(values=c('cornflowerblue','indianred'))

genos<-c('A','B','C')
mergesumm<-radial_plot_space(bin_tot[order(a,b,c)],genos)[order(a,b,c)]

dummysumm<-data.table(label1=c('WT','GAL4','GAL','GAL80'))
dummysumm[,c('x1','y1'):=list(c(0,0,0,0),c(0,1.7,2.7,3.7))]
rotation_factor<-pi/2.5

labelDT1<-data.table(labeldfun2(mergesumm,rotation_factor=rotation_factor))
labelDT1[, significant:=as.logical(sig)]
range(mergesumm$diff)
difflim<-c(-0.024,0.024)
diffcols<-c('cornflowerblue','white','indianred')
pBinVsTot<-pRAD + 
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_4, start = start_4, end = end_4, fill = diff_mean),col='grey80', data = mergesumm[r_4==4])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_3, start = start_3, end = end_3, fill = diff_mean),col='grey80', data = mergesumm[r_3==3])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_2, start = start_2, end = end_2, fill = diff_mean),col='grey80', data = mergesumm[r_2==2])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_1, start = start_1, end = end_1, fill = diff_mean),col='grey80', data = mergesumm[r_1==1])+
	geom_circle(aes(x0=0,y0=0,r=1),col='grey80')+
	geom_text(data=dummysumm,aes(x=x1,y=y1+0.1,label= label1))+
	geom_text(data=labelDT1,aes(x=x,y=y,label= genotypei))+
	geom_point(data=labelDT1[sig==T],aes(x=x,y=y-0.2),colour='black',size=2)+geom_point(data=labelDT1[sig==T],aes(x=x,y=y-0.2,colour=significant))+
	scale_fill_gradientn(colours= diffcols,name='deviation binary copy phenotype',limits= difflim)+theme(panel.border=element_blank())
pBinVsTot

pBinVsTot


sample(LETTERS,4)
plotlist<-list(pBinVsTot)
newdir<-paste0('.',figout,'180627-pWTZV-RadialInteractionPlots')
if(!dir.exists(newdir))system(paste0('mkdir ', newdir))
setwd(newdir)
w<-12;h<-w*0.7
whs<-c(0.5,1,3)
h<-8;w<-h*1.2
ggsave(paste0('180627-RadialInteractionPlots-4.png'), pBinVsTot,height=h,width=w)
h<-4;w<-h*1.2
ggsave(paste0('180627-RadialInteractionPlots-5.png'), pBinVsTot1,height=h,width=w)
setwd(head.dir)

# compare full data set to total copies
tot_full<-merge(DT2[,c('tot_full',grepincols(DT2,'totalcopy')),with=F],DT3,by='tot_full')[order(a,b,c,d,e,f)]
tot_full[,diff_mean:= full_mean-totalcopy_mean]
tot_full[,diff_se:=sqrt(totalcopy_se^2+full_se^2)]
tot_full[,diff_df:=full_N+totalcopy_N-1]
attach(tot_full);pval<-data.table(t(apply(data.frame(totalcopy_mean,full_mean,totalcopy_sd,full_sd,totalcopy_N,full_N),1,function(x){
	test<-t.test2(m1=x[1],m2=x[2],s1=x[3],s2=x[4],n1=x[5],n2=x[6])
})))$p.value;detach(tot_full)
tot_full $diff_p.val<-pval
tot_full[,p.adj:=p.adjust(diff_p.val)]
tot_full[,sig:=p.adj<0.05]

ggplot(tot_full,aes(x=full_mean,xmax=full_upper,xmin=full_lower,y=totalcopy_mean,ymax=totalcopy_upper,ymin=totalcopy_lower,col=sig))+geom_errorbar()+geom_errorbarh()+geom_point()+geom_abline()+coord_flip()
tot_full[sig==T]
genos<-c('A','B','C','D','E','F')
mergesumm<-radial_plot_space(tot_full[order(a,b,c,d,e,f)],genos)[order(a,b,c,d,e,f)]

dummysumm<-data.table(label1=c('WT','GAL4 plasmid','GAL4 genome','GAL3 plasmid','GAL3 genome','GAL80 plasmid','GAL80 genome'))
dummysumm[,c('x1','y1'):=list(rep(0,7),c(0,seq(1:6)+0.2))]
rotation_factor<-pi/2.2


labelDT1<-data.table(labeldfun2(mergesumm,rotation_factor=rotation_factor))
labelDT1[, significant:=as.logical(sig)]
range(mergesumm$diff_mean)
difflim<-c(-0.035,0.035)
diffcols<-c('cornflowerblue','white','indianred')
pTotVsFull<-pRAD + 
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_7, start = start_7, end = end_7, fill = diff_mean),col='grey80', data = mergesumm[r_7==7])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_6, start = start_6, end = end_6, fill = diff_mean),col='grey80', data = mergesumm[r_6==6])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_5, start = start_5, end = end_5, fill = diff_mean),col='grey80', data = mergesumm[r_5==5])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_4, start = start_4, end = end_4, fill = diff_mean),col='grey80', data = mergesumm[r_4==4])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_3, start = start_3, end = end_3, fill = diff_mean),col='grey80', data = mergesumm[r_3==3])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_2, start = start_2, end = end_2, fill = diff_mean),col='grey80', data = mergesumm[r_2==2])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_1, start = start_1, end = end_1, fill = diff_mean),col='grey80', data = mergesumm[r_1==1])+
	geom_circle(aes(x0=0,y0=0,r=1),col='grey80')+
	geom_text(data=dummysumm,aes(x=x1,y=y1+0.1,label= label1))+
#	geom_text(data=labelDT1,aes(x=x,y=y,label= genotypei))+
	geom_point(data=labelDT1[sig==T],aes(x=x,y=y-0.2),colour='black',size=2)+geom_point(data=labelDT1[sig==T],aes(x=x,y=y-0.2,colour=significant))+
	scale_fill_gradientn(colours=diffcols,name='deviation total copy phenotype',limits= difflim)+theme(panel.border=element_blank())
pTotVsFull

# compare full data set to total copies
bin_full<-merge(DT1[,c('bin_full',grepincols(DT1,'binary')),with=F],DT3,by='bin_full')[order(a,b,c,d,e,f)]
bin_full[,diff_mean:= full_mean-binary_mean]
bin_full[,diff_se:=sqrt(binary_se^2+full_se^2)]
bin_full[,diff_df:=full_N+binary_N-1]
attach(bin_full);pval<-data.table(t(apply(data.frame(binary_mean,full_mean,binary_sd,full_sd,binary_N,full_N),1,function(x){
	test<-t.test2(m1=x[1],m2=x[2],s1=x[3],s2=x[4],n1=x[5],n2=x[6])
})))$p.value;detach(bin_full)
bin_full $diff_p.val<-pval
bin_full[,p.adj:=p.adjust(diff_p.val)]
bin_full[,sig:=p.adj<0.05]

ggplot(bin_full,aes(x=full_mean,xmax=full_upper,xmin=full_lower,y=binary_mean,ymax=binary_upper,ymin=binary_lower,col=sig))+geom_errorbar()+geom_errorbarh()+geom_point()+geom_abline()+coord_flip()
bin_full[sig==T]
genos<-c('A','B','C','D','E','F')
mergesumm<-radial_plot_space(bin_full[order(a,b,c,d,e,f)],genos)[order(a,b,c,d,e,f)]

dummysumm<-data.table(label1=c('WT','GAL4 plasmid','GAL4 genome','GAL3 plasmid','GAL3 genome','GAL80 plasmid','GAL80 genome'))
dummysumm[,c('x1','y1'):=list(rep(0,7),c(0,seq(1:6)+0.2))]


rotation_factor<-pi/2.2
labelDT1<-data.table(labeldfun2(mergesumm,rotation_factor=rotation_factor))
labelDT1[, significant:=as.logical(sig)]
range(mergesumm$diff_mean)
difflim<-c(-0.06,0.06)
diffcols<-c('cornflowerblue','white','indianred')
pBinVsFull<-pRAD + 
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_7, start = start_7, end = end_7, fill = diff_mean),col='grey80', data = mergesumm[r_7==7])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_6, start = start_6, end = end_6, fill = diff_mean),col='grey80', data = mergesumm[r_6==6])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_5, start = start_5, end = end_5, fill = diff_mean),col='grey80', data = mergesumm[r_5==5])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_4, start = start_4, end = end_4, fill = diff_mean),col='grey80', data = mergesumm[r_4==4])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_3, start = start_3, end = end_3, fill = diff_mean),col='grey80', data = mergesumm[r_3==3])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_2, start = start_2, end = end_2, fill = diff_mean),col='grey80', data = mergesumm[r_2==2])+
	geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r =r_1, start = start_1, end = end_1, fill = diff_mean),col='grey80', data = mergesumm[r_1==1])+
	geom_circle(aes(x0=0,y0=0,r=1),col='grey80')+
	geom_text(data=dummysumm,aes(x=x1,y=y1+0.1,label= label1))+
#	geom_text(data=labelDT1,aes(x=x,y=y,label= genotypei))+
	geom_point(data=labelDT1[sig==T],aes(x=x,y=y-0.2),colour='black',size=2)+geom_point(data=labelDT1[sig==T],aes(x=x,y=y-0.2,colour=significant))+
	scale_fill_gradientn(colours=diffcols,name='deviation from binary phenotype',limits= difflim)+theme(panel.border=element_blank())
pBinVsFull


