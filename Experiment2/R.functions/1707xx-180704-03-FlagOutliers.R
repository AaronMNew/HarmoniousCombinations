# Scripts for flagging outliers.
source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'

head.dir<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis'
layout.dir<-'./layout.sample.data'
date<-'1707xx' # beginning date clones were measured
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/'
outputdata<-'./output.data/'
setwd(head.dir)

br<-seq(-0.5,2.7,length=61)
library(tidyr)

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
	
merge_recurse2<-function(qDT=NULL,dt.list,by.list,names.list=NULL,all_lv=T){
	x<-qDT
	if(is.null(qDT))x<-dt.list[[1]]
	if(length(dt.list)!=(length(by.list)))stop('list lengths look funny')
	i<-1
	while(i <= length(dt.list)){
		x<-merge(x,dt.list[[i]],by=by.list[[i]],all= all_lv)
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

plateplot<-function(DT,phenot,collim=NULL){
	# DT must have r, c and phenotype columns, measurement.plate.gal and plate.gal
	x<-copy(DT)
	x[,row:=factor(r,levels=LETTERS[8:1])]
	x[,column:=factor(c,levels=1:12)]
	setnames(x,phenot,'pheno')
	ggplot(x,aes(column,row))+geom_tile(aes(fill=abs(pheno)))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),na.value='transparent',limits=collim)+
		theme(legend.title=element_blank())+scale_colour_manual(values=c('black'),na.value='transparent')+
		facet_grid(measurement.date.gal~plate.gal)+ggtitle(paste0(phenot," ",unique(x$source.plate)))	
}

notcols<-function(DT,x,grepl=T){
	xa<-x
	if(grepl==T)xa<-unlist(lapply(x,function(xa)colnames(DT)[grepl(xa,colnames(DT))]))
	colnames(DT)[!colnames(DT)%in%xa]
}
castform<-function(leftside, rightside){
	left<-paste0(leftside,collapse='+')
	right<-paste0(rightside,collapse='+')
	paste0(left,'~',right)
}
snames<-names(summary.func(1:10))
# rside<-c('single')
# lside<-notcols(mDT,c(rside,snames))
# form<-castform(lside, rside)
# cDT<-cast(mDT,as.formula(form),value='mean')


###################################################
### load data
###################################################
DT0<-fread(paste0(outputdata,'1707xx-180704-preprocessed_data.txt'),sep='\t',header=T)

###################################################
### set variables and add labels
###################################################
snames<-names(summary.func(1:10))

DT0[,aa:=GAL4.clone.named]
DT0[,bb:=GAL3.clone.named]
DT0[,cc:=GAL80.clone.named]
DT0[,abc:=clone.named]

DT0[,G4.single:=cc=='GAL80.WT'&bb=='GAL3.WT']
DT0[,G3.single:=cc=='GAL80.WT'&aa=='GAL4.WT']
DT0[,G80.single:=bb=='GAL3.WT'&aa=='GAL4.WT']

DT0[,WT:=aa=='GAL4.WT'&bb=='GAL3.WT'&cc=='GAL80.WT']
DT0[,GAL3vGAL80:=aa=='GAL4.WT']
DT0[,GAL4vGAL80:=bb=='GAL3.WT']
DT0[,GAL3vGAL4:=cc=='GAL80.WT']

DT0[,G4vG3.double:=GAL3vGAL4==T&GAL3vGAL80==F]
DT0[,G4vG80.double:=GAL3vGAL4==T&GAL3vGAL80==F]

DT0[,replicate.allele:=1:.N,by=c('allele.named')]
DT0[,replicate.clone.named:=1:.N,by=c('clone.named')]
table(DT0[,'replicate.allele',with=F])
table(DT0[,'replicate.clone.named',with=F])

phenocols<-c('yfp.mean','mean.on.fraction','var.on.fraction','fracon','dip.stat')
galcols<-c(c('growth.rate','phenotypic.index'),paste(phenocols,'.gal',sep=''))
glucols<-paste(phenocols,'.glu',sep='')
coordcols<-colnames(DT0)[ grep('Dim',colnames(DT0))]
dependents<-c(galcols,glucols,coordcols)
genocols<-c('genotype','clone.named','GAL4.allele', 'GAL80.allele', 'GAL3.allele')
clone.named.set<-colnames(DT0)[grep('clone.name',colnames(DT0))]
groupings<-c('GAL3.clone.named','GAL80.clone.named','GAL4.clone.named','clone.named')
GALgenes<-c('GAL3','GAL80','GAL4')
dens.cols<-grepincols(DT0,'dens')

G3alleles<-na.exclude(unique(DT0$GAL3.clone.named))
G80alleles<-na.exclude(unique(DT0$GAL80.clone.named))
G4alleles<-na.exclude(unique(DT0$GAL4.clone.named))
G3bin<-data.table(t(ldply(lapply(G3alleles,function(x){
	grepl(x,DT0$GAL3.clone.named)
}))))
table(apply(G3bin,1,sum))
G80bin<-data.table(t(ldply(lapply(G80alleles,function(x){
	grepl(x,DT0$GAL80.clone.named)
}))))
table(apply(G80bin,1,sum))

G4bin<-data.table(t(ldply(lapply(G4alleles,function(x){
	grepl(x,DT0$GAL4.clone.named)
}))))
table(apply(G4bin,1,sum))

binaryalleles<-data.table(G3bin,G80bin,G4bin)
colnames(binaryalleles)<-c(G3alleles,G80alleles,G4alleles)

batchcols<-c('samp_id.glu','samp_id.gal',
	'biol.rep.gal',
	'replicate.gal')

batchcols1<-c('replicate.glu', # 
	'number.mutants.glu', # this scores the number of mutations present compared to WT (0,1, or 2)
	'id',# id joins the glucose sample with the same one the next day that was measured in galactose (date term is the date gal measurement was made)
	'samp_id.glu',
	'measurement.date.glu',
	'measurement.date.gal',
	'media.batch.date.glu',
	'source.plate.glu',
	'plate.name.in.source.fcs.glu',
	'plate.glu',
	'plate.gal',
	'biol.rep.glu',
	'onemut.glu',
	'twomut.glu',
	'samp_id.gal',
	'biol.rep.gal',
	'replicate.gal',
	'replicate.allele',
	'replicate.clone.named')
	
batchcols2<-c('measurement.date.glu',
	'GAL3.clone.named',
	'GAL80.clone.named',
	'GAL4.clone.named',
	'GAL3.allele',
	'GAL80.allele',
	'GAL4.allele',
	'genotype',
	'measurement.date.gal',
	'source.plate',
	'clone.named',
	'tx.date'
	)

###################################################
### a bit of code for plots to check how you're doing as you go along
###################################################

# check how you're doing
allele.plot<-function(DT0,plot=T){
	bycols<-c('allele.named','WT', 'G3.single','G80.single','G4.single')
	gr<-DT0[G3.single==T|G80.single==T|G4.single,summary.func(growth.rate,'GR'),by=bycols]
	pi<-DT0[G3.single==T|G80.single==T|G4.single,summary.func(phenotypic.index,'PI'),by=bycols]
	DTm<-merge(gr,pi,by=bycols)
	mDT1<-melt(DTm,id=c(colnames(DTm)[!colnames(DTm)%in%c('G3.single','G80.single','G4.single')]))
	mDT<-mDT1[value==T]
	p<-ggplot(mDT,aes(PI_mean,GR_mean,xmin=PI_lower,xmax=PI_upper,ymax=GR_upper,ymin=GR_lower))+geom_errorbar()+geom_errorbarh()+geom_point()+
		geom_point(data=DTm[WT==T],aes(PI_mean,GR_mean),col='red',size=3,shape=21,stroke=1)+facet_grid(~variable)+theme_light()
	if(plot==F){return(mDT)}
	if(plot==T){return(p)}
}
clone.plot<-function(DT0,plot=T){
	bycols<-c('clone.named','WT', 'G3.single','G80.single','G4.single')
	gr<-DT0[G3.single==T|G80.single==T|G4.single,summary.func(growth.rate,'GR'),by=bycols]
	pi<-DT0[G3.single==T|G80.single==T|G4.single,summary.func(phenotypic.index,'PI'),by=bycols]
	DTm<-merge(gr,pi,by=bycols)
	mDT1<-melt(DTm,id=c(colnames(DTm)[!colnames(DTm)%in%c('G3.single','G80.single','G4.single')]))
	mDT<-mDT1[value==T]
	p<-ggplot(mDT,aes(PI_mean,GR_mean,xmin=PI_lower,xmax=PI_upper,ymax=GR_upper,ymin=GR_lower))+geom_errorbar()+geom_errorbarh()+geom_point()+
		geom_point(data=DTm[WT==T],aes(PI_mean,GR_mean),col='red',size=3,shape=21,stroke=1)+facet_grid(~variable)+theme_light()
	if(plot==F){return(mDT)}
	if(plot==T){return(p)}
}
pTechrep<-function(DT0,genos='clone.named',plot=T){
	batchcols<-c('plate.gal','plate.glu',
		'biol.rep.glu','biol.rep.gal','samp_id.glu','samp_id.gal',
		'replicate.allele',
		'replicate.clone.named','r','c','source.plate')
	qDT<-copy(DT0)
	groupings<-c(genos,'GAL3vGAL80','GAL4vGAL80','GAL3vGAL4')
	phenocols<-c('growth.rate','phenotypic.index','yfp.mean.gal','yfp.mean.glu','fracon.gal','fracon.glu')
	phenotype<-'growth.rate'
	qSub<-qDT[,c(phenotype,groupings,batchcols),with=F]
	setnames(qSub,phenotype,'x')
	qSub[,mean:=mean(x),by=genos]
	mDT<-melt(qSub,id=c('x','mean',genos,batchcols))
	ggplot(mDT[value==T],aes(mean,x))+geom_point()+geom_abline()+facet_wrap(~variable)
	
	qSub<-qDT[,c(phenocols,groupings,batchcols),with=F]
	qSub[,dumm:=1]
	qSumm1<-summary.func.all(qSub,c('dumm',phenocols),bylist= genos)
	qSumm<-qSumm1[,!c(grepincols(qSumm1,'dumm')),with=F];qSub[,dumm:=NULL]
	mSub1<-melt(qSub,id=c(groupings,batchcols))
	setnames(mSub1,c('variable','value'),c('phenotype','observed'))
	mSub<-melt(mSub1,id=c('phenotype','observed', genos,batchcols))[value==T]
	mSub[,value:=NULL]
	setnames(mSub,'variable','pairing')
	
	mSumm1<-melt(qSumm,id= genos)
	mSumm1[,c('phenotype','stat'):=colsplit(variable,'_',c('phenotype','stat'))];mSumm1[,variable:=NULL]
	form<-paste0(paste(genos,'phenotype',sep='+'),'~stat')
	mSumm<-dcast(mSumm1,as.formula(form),value.var='value')
	mDT<-merge(mSumm,mSub,by=c(genos,'phenotype'))
	mDT[,diff:=observed-mean]
	mDT[,diff_scale:=scale(diff),by=c('phenotype')]
	mDT[,diff_3SD:=abs(diff_scale)>3]
	mDT[,observed_scale:=scale(observed),by=c('phenotype')]
	mDT[,mean_scale:=scale(mean),by=c('phenotype')]
	mDT[,var_exp_phen:=varexplained(observed,mean),by=c('phenotype')]
	mDT[,var_exp_pair:=varexplained(observed,mean),by=c('pairing')]
	mDT[,var_exp_fac_phen:=paste0(round(var_exp_phen,2))]
	mDT[,var_exp_fac_pair:=paste0(round(var_exp_pair,2))]
	mDT[,pairing2:=paste0(pairing,'\n', var_exp_fac_pair)]
	mDT[,phenotype2:=paste0(phenotype,'\n', var_exp_fac_phen)]
	pTechRep<-ggplot(mDT,aes(mean_scale,observed_scale,col= diff_3SD))+geom_point(shape=21,alpha=0.2)+geom_abline()+facet_grid(pairing2~phenotype2,scales='free')+
		ylab('observed measurements')+xlab('within-genotype mean of measurement')
	if(plot==T){return(pTechRep)}
	if(plot==F){return(mDT)}
}



###################################################
### look at technical variation by and allele.named sets
###################################################

p1<-clone.plot(DT0=DT0,plot=T)
p2<-allele.plot(DT0,plot=T)
do.call(grid.arrange,list(p1,p2))
pTechrep(DT0,'clone.named',plot=T)
pTechrep(DT0,'allele.named',plot=T)



########
### look at all plates looking for systematic inoculation problems
########

DTemp<-copy(DT0)
mu_gal<-mu_glu<-0
mDT<-as.data.table(pTechrep(DT0,'allele.named',plot=F))
problemclones<-unique(mDT[diff_3SD==T]$allele.named)
DTemp[,problemclone:=allele.named%in%problemclones]
DTemp[problemclone==F,problemclone:=NA]
DTemp[problemclone==T]
pSplit<-split(DTemp,by='source.plate')
w<-14.2;h<-8.2

# ONLY RUN ONCE ...
# newdir<-paste0(figout,'1707xx-AllPlatesPhenotypes')
# if(!dir.exists(newdir))system(paste0('mkdir ', newdir))
# setwd(newdir)
# lapply(seq(length(pSplit)),function(i){
	# x<-copy(pSplit[[i]])
	# x[,log10.yfp.mean.gal:=log10(yfp.mean.gal-mu_gal+1000)-3]
	# x[,log10.yfp.mean.glu:=log10(yfp.mean.gal-mu_glu+1000)-3]
	# p1<-plateplot(x,'fracon.gal',collim=c(0,1))+geom_point(aes(colour=problemclone))
	# p2<-plateplot(x,'fracon.glu',collim=c(0,1))+geom_point(aes(colour=problemclone))
	# p3<-plateplot(x,'phenotypic.index',collim=c(0,2))+geom_point(aes(colour=problemclone))
	# p4<-plateplot(x,'log10.yfp.mean.gal',collim=log10(c(1,10^4)))+geom_point(aes(colour=problemclone))
	# p5<-plateplot(x,'log10.yfp.mean.glu',collim=log10(c(1,10^4)))+geom_point(aes(colour=problemclone))
	# p6<-plateplot(x,'growth.rate',collim=c(0,0.3))+geom_point(aes(colour=problemclone))
	# ggsave(paste0(date,'-180704-AllPlates-',unique(x$source.plate),'.png'),do.call(grid.arrange,list(p1,p4,p2,p5,p3,p6,ncol=2)),height=h,width=w)
# })
# setwd(head.dir)

# Plate_023 source plate 170823-Plate_015 row A was mis-inoculated

# source.plate Plate_025 grew differentially, likely due to some delay in the morning.
# plateplot(DT0[source.plate=='Plate_025'],'growth.rate')
qM<-merge(DT0[source.plate=='Plate_025'&plate.gal=='170708-Plate_001'],DT0[source.plate=='Plate_025'&plate.gal=='170708-Plate_002'],by=c('r','c'))
# qplot(qM$growth.rate.x,qM[,growth.rate.x-growth.rate.y])
corfac<-mean(qM[growth.rate.x>0.1,growth.rate.x-growth.rate.y])
DTemp<-copy(DT0)
DTemp[source.plate=='Plate_025'&plate.gal=='170708-Plate_002'&growth.rate>0.1,growth.rate:=growth.rate+corfac]
# plateplot(DTemp[source.plate=='Plate_025'],'growth.rate')
# correct DT0
DT0[source.plate=='Plate_025'&plate.gal=='170708-Plate_002'&growth.rate>0.1,growth.rate:=growth.rate+corfac]


# Source.plate Plate_034 is missing most phenotypes. look more closely there.
# Souce plate Plate_035 C11 should be masked. mis-transformed.
# Souce plate Plate_037 H1 should be masked. mis-transformed or contaminant.
# Plate_041 plate.gal 170718-Plate_013 had a problem during acquisition. Delete.
# Plate_041 plate.gal 170720-Plate_021 B11 B12 are contaminated
# only one measurement of source.plate Plate_047 for most phenotypes because unpaired, but then sometimes there are two measurements in the "NA" category
# and Plate_048 and Plate_063
# Plate_060 G6 is contaminant and should be masked
# Fracon.glu measurement of F7 in plate.gal 170825-Plate_011 should be masked; some kind of contaminant. Just going to remove the

#plateplot(DT0[source.plate=='Plate_023'],'growth.rate')
DT0<-DT0[!(plate.gal=='170823-Plate_015'&r=='A')]
DT0<-DT0[!(source.plate=='Plate_035'&rc=='C11')]
DT0<-DT0[!(source.plate=='Plate_037'&rc=='H1')]
#plateplot(DT0[plate.gal=='170718-Plate_013'|plate.gal=='170718-Plate_014'],'growth.rate')
#plateplot(DT0[source.plate=='Plate_041'],'growth.rate')
DT0<-DT0[!(plate.gal=='170718-Plate_013')]
#plateplot(DT0[plate.gal=='170720-Plate_021'&!(rc=='B11'|rc=='B12')],'growth.rate')
DT0<-DT0[!(plate.gal=='170720-Plate_021'&(rc=='B11'|rc=='B12'))]
DT0<-DT0[!(source.plate=='Plate_060'&rc=='G6')]
DT0[plate.gal=='170825-Plate_011'&rc=='F7',fracon.glu:=NA]
DT0[plate.gal=='170825-Plate_011'&rc=='F7',phenotypic.index:=NA]

# sp plate Plate_057 plate.gal 170825-Plate_011 is mis-inoculated
test<-copy(mDT[pairing=='GAL3vGAL4'&phenotype=='phenotypic.index'])
test[,diff:= mean_scale-observed_scale]
test[,diff.score:=abs(diff)>0.3]
pTechRep<-ggplot(test,aes(mean_scale,observed_scale,col= diff.score))+geom_point(shape=21,alpha=0.2)+geom_abline()+facet_grid(pairing2~phenotype2,scales='free')+
	ylab('observed measurements')+xlab('within-genotype mean of measurement')
# pTechRep
test[diff.score==T]
# plateplot(DT0[source.plate=='Plate_057'],'fracon.glu')
# plateplot(DT0[source.plate=='Plate_057'&!(plate.gal=='170825-Plate_011'&(r=='G'|r=='H'))],'fracon.glu')
DT0<-DT0[!(plate.gal=='170825-Plate_011'&(r=='G'|r=='H'))]


# sp Plate_003 170718-Plate_003 row F was mis-inoculated or 170718-Plate_004 row B was inoculated into 170718-Plate_004 row F.
# plateplot(DT0[source.plate=='Plate_003'],'growth.rate')
DT0[plate.gal=='170718-Plate_004',grepincols(DT0,'clone'),with=F]
DT0[plate.gal=='170718-Plate_004'&r=='F',grepincols(DT0,'clone'),with=F]
# hist(DT0[cc=='GAL80.06'&aa=='GAL4.WT']$phenotypic.index)
# hist(DT0[cc=='GAL80.05'&aa=='GAL4.WT']$phenotypic.index)
G3s<-DT0[plate.gal=='170718-Plate_004'&r=='F',grepincols(DT0,'clone'),with=F]$GAL3.clone.named
qSumm<-DT0[G80.single==T,summary.func(phenotypic.index),by='cc']
qSumm[,const:=mean>1.4]
ccs<-qSumm[const==T]$cc
# ggplot(DT0[bb%in%G3s&aa=='GAL4.WT'&!(cc%in%ccs),summary.func(growth.rate),by=c('clone.named','bb')],aes(bb,mean))+geom_boxplot()
# all of the GAL3's in that row are dead  in WT GAL4 and GAL80 backgrounds and GAL80.05 and GAL80.06 behave like WT, so this was an mis-inoculation of 170718-Plate_003 row F
DTemp<-DT0[!(plate.gal=='170718-Plate_003'&r=='F')]
# plateplot(DTemp[source.plate=='Plate_003'],'growth.rate')
DT0<-DT0[!(plate.gal=='170718-Plate_003'&r=='F')]


# Growth profiles of sp Plate_006 plates 170706-Plate_005 170706-Plate_006 look different enough to have been inoculation errors. Especially rows G and H of 170706-Plate_006 look weird.
	# for this comparison will be tricky. First remind yourself of what the GAL4 GAL3 GAL80 variants are. 170706-Plate_005's rows CD look very much like rows GH
# plateplot(DT0[source.plate=='Plate_006'],'growth.rate')
# table(DT0[source.plate=='Plate_006',grepincols(DT0,'clone'),with=F])
DT0[plate.gal=='170718-Plate_004'&r=='F',grepincols(DT0,'clone'),with=F]
# hist(DT0[cc=='GAL80.11'&aa=='GAL4.WT']$phenotypic.index)
# hist(DT0[cc=='GAL80.12'&aa=='GAL4.WT']$phenotypic.index)
# no notes on how long the samples sat at room temperature. it seems that the E,F, G and H were mis-inoculated nad the plate over-grew in general. 
# last possibility is that it was inoculated twice.
DTemp<-copy(DT0[plate.gal=='170706-Plate_006'])
DTemp[,cell.per.ml.gal:=cell.per.ml.gal/2]
DTemp[,growth.rate:=log(cell.per.ml.gal/(cell.per.ml.glu*9/150))/12]
DTemp[,plate.gal:='reduced.density.by.half']
# plateplot(rbind(DTemp,DT0[source.plate=='Plate_006']),'growth.rate')
# this doesn't seem likely.
# my decision is to remove those rows and keep the rest
# plateplot(DT0[source.plate=='Plate_006'&!(plate.gal=='170706-Plate_006'&(r=='E'|r=='F'|r=='G'|r=='H'))],'growth.rate')
DT0<-DT0[!(plate.gal=='170706-Plate_006'&(r=='E'|r=='F'|r=='G'|r=='H'))]
# plateplot(DT0[source.plate=='Plate_006'],'growth.rate')

# similar situation with SP Plate_008 rows A and B have very different growth profiles in the different measurement plates. Also 170706-Plate_009 D11 and DT12 are contaminated
# plateplot(DT0[source.plate=='Plate_008'],'growth.rate')
# table(DT0[source.plate=='Plate_008',grepincols(DT0,'clone'),with=F])
DT0[plate.gal=='170706-Plate_009'&(r=='A'|r=='B'),grepincols(DT0,'clone'),with=F]
G3s<-DT0[plate.gal=='170706-Plate_009'&(r=='A'|r=='B'),grepincols(DT0,'clone'),with=F]$GAL3.clone.named

qSumm<-DT0[G80.single==T,summary.func(phenotypic.index),by='cc']
qSumm[,const:=mean>1.4]
ccs<-qSumm[const==T]$cc
# ggplot(DT0[bb%in%G3s&aa=='GAL4.WT'&!(cc%in%ccs),summary.func(growth.rate),by=c('clone.named','bb')],aes(bb,mean))+geom_boxplot()+coord_flip()

cDT<-dcast(DT0[source.plate=='Plate_008',summary.func(growth.rate,'growth.rate'),by=c('plate.gal','allele.named')],formula(allele.named~plate.gal),value.var= 'growth.rate_mean')
colnames(cDT)<-c('allele.named','p1','p2')
# qplot(cDT[,p1],cDT[,p2])+geom_abline()
corfac<-cDT[,mean(p1-p2)]
cDT[,pCor:=p2+corfac]
# qplot(unlist(cDT[,p1]),unlist(cDT[,3]))+geom_abline()
# qplot(cDT[,p1],cDT[,pCor])+geom_abline()
DTx<-copy(DT0[source.plate=='Plate_008'])
DTx[,gr.cor:=growth.rate]
# table(DTx$plate.gal)
DTx[plate.gal=='170706-Plate_010',gr.cor:=growth.rate+corfac]
# plateplot(DTx,'gr.cor')

# they are roughly correlated with no obvoius mis-inoculations, also as measured by phenotypic index. I imagine these was ain inoculation problem (not enough cells in a few rows) and that the plates sat longer at room temperature
# this is evident because dense measurements have higher YFP.gal expression, suggesting that these are un-matured fluorophores.
# I think they're just crap measurements. Will change the growth rate to the corrected one to correct for systematic biases

DT0[source.plate=='Plate_008'&plate.gal=='170706-Plate_010',growth.rate:=growth.rate+corfac]
# clone.plot(DT0)

# now this plate has one particular GAL80 alelle that lookst like crap.
# look again find out source.plate Plate_008 has problem with GAL80.15 outlying GAL80 sample. this grew very fast and might be contaminant
DTal<-allele.plot(DT0,plot=F)
DTal[variable=='G80.single'&GR_upper>0.3]
DT0[,arb:=factor(1:nrow(DT0))]
DTx<-na.exclude(DT0[,c('arb','growth.rate','median.fsc.gal','median.ssc.gal')])

# countaminants could be samples with very different shape characteristics given known growth-rate shape relationhips
# pairs.pears(DT0[,c('growth.rate','median.fsc.gal','median.ssc.gal')])
# make model of that relationship and calculated deviations from it
mod1<-lm(growth.rate~median.fsc.gal+median.ssc.gal,DTx)
summary(mod1)
DTx[,size_shape_residual:=growth.rate-mod1$fitted.values]
DT0x<-merge(DTx,DT0,by=c('arb','growth.rate','median.fsc.gal','median.ssc.gal'))
# hist(DT0[cc=='GAL80.15']$growth.rate)
rc<-DT0[cc=='GAL80.15'&growth.rate>0.29]$rc
plat<-DT0[cc=='GAL80.15'&growth.rate>0.29]$source.plate
# look at these variables across the plates.
# plateplot(DT0[source.plate==plat],'median.fsc.gal')
# plateplot(DT0[source.plate==plat],'median.fsc.glu',collim=c(2000,6000))
# plateplot(DT0[source.plate==plat],'median.fsc.gal',collim=c(5000,11000))
# key plots are here
# plateplot(DT0x[source.plate==plat],'size_shape_residual')
# dev.new();plateplot(DT0x[source.plate==plat],'growth.rate')
# the suspected contaminated wells actually have very low residuals of expected growth rate given size and shape measurements
# it appears their cell shape characteristics match their measured growth rate variables. Something weird happened with inoculation then.


# Plate_009 growth profiles look different between samples. Middle removed sample might have been contaminated; seems this way because wells usrrpounding it are high growth
# plateplot(DT0[source.plate=='Plate_009'],'growth.rate')
cDT<-dcast(DT0[source.plate=='Plate_009',summary.func(growth.rate,'growth.rate'),by=c('plate.gal','allele.named')],formula(allele.named~plate.gal),value.var= 'growth.rate_mean')
colnames(cDT)<-c('allele.named','p1','p2')
# qplot(cDT[,p1],cDT[,p2])+geom_abline()
corfac<-cDT[,mean(p1-p2)]
cDT[,pCor:=p2+corfac]
# qplot(unlist(cDT[,p1]),unlist(cDT[,3]))+geom_abline()
# qplot(cDT[,p1],cDT[,pCor])+geom_abline()
# inoculations all look fine, and systematic bias between plates not obvoius. just leaving it as it is: crap data


# Plate_010 has very different growth profiles. seems that 170706-Plate_013 row C might have been inoculated twice, but not across all wells.
# plateplot(DT0[source.plate=='Plate_010'],'growth.rate')
temp<-copy(DT0[source.plate=='Plate_010'])
temp[,gr.cor:=growth.rate]
# divide row C by two to see if that changes things.
temp[plate.gal=='170706-Plate_013'&r=='C',gr.cor:=growth.rate/2]
# plateplot(temp,'gr.cor')
cDT<-dcast(temp[source.plate=='Plate_010',summary.func(gr.cor,'growth.rate'),by=c('plate.gal','allele.named')],formula(allele.named~plate.gal),value.var= 'growth.rate_mean')
colnames(cDT)<-c('allele.named','p1','p2')
# qplot(cDT[,p1],cDT[,p2])+geom_abline()
test<-melt(DT0[source.plate=='Plate_010'&r=='C',c('cell.per.ml.gal','plate.gal','c')],id=c('plate.gal','c'))
test2<-(dcast(test,formula(c~plate.gal),value.var='value'))
colnames(test2)<-c('c','plate1','plate2')
test2[,diff:=plate1/plate2]
# looks like it was just a crap plate :/

# no fluorescence but high growth
plateplot(DT0[source.plate==unique(DT0[clone.named=='GAL3.42 GAL80.10 GAL4.WT']$source.plate)],'growth.rate')
# this is a contaminant
DT0<-DT0[!clone.named=='GAL3.42 GAL80.10 GAL4.WT']

# discovered 181025
problemclones<-c('GAL3.13 GAL80.WT GAL4-L868G','GAL3.13 GAL80.WT GAL4-L868S')
problemplates<-unique(DT0[clone.named%in%problemclones]$plate.gal)
problemplates1<-unique(DT0[clone.named%in%problemclones]$plate.glu)
DT0[plate.gal%in%problemplates&clone.named%in%problemclones]
DT0[plate.gal%in%&clone.named%in%problemclones]
plateplot(DT0[plate.gal%in%problemplates],'yfp.mean.gal')
plateplot(DT0[plate.gal%in%c(problemplates,"170712-Plate_011","170712-Plate_012","170712-Plate_015","170712-Plate_016")],'yfp.mean.gal')
# Iinoculateion problems with D1 and D2 for plates "170712-Plate_013" "170712-Plate_014"
# mask those

DTemp<-DT0[!samp_id.gal%in%c("170712-Plate_013.D1","170712-Plate_014.D1","170712-Plate_013.D2","170712-Plate_014.D2")]
plateplot(DTemp[plate.gal%in%c(problemplates,"170712-Plate_011","170712-Plate_012","170712-Plate_015","170712-Plate_016")],'yfp.mean.gal')
DT0<-DTemp

# discovered 181030
# these two clones must have had a problem with transformation. complete outliers
DT0[clone.named%in%c("GAL3.WT GAL80.01 GAL4.18","GAL3.WT GAL80.01 GAL4.19")]<-NA


########
### make another ptechrep
########
DT0[,clone.named.old:=clone.named]
DT0[,clone.named:=applyPaste(data.frame(GAL3.clone.named,GAL80.clone.named,GAL4.clone.named),collapse=' ')]
pTechrep(DT0=DT0,genos='allele.named')
pTechrep(DT0=DT0,genos='clone.named')
# still problems with clone.named samples but for the most part we're good to go.
clone.plot(DT0,plot=T)

########
### write all modifications to new table
########

# annotate by hand hwich columns are worth keeping
write.table(data.table(colnames(DT0)),paste0(layout.dir,'/1707xx-DT0cols.txt'),row.names=F)
keepcolDT<-fread(paste0(layout.dir,'/1707xx-180710-ColumnsToKeep.txt'),sep='\t',header=T)
keep<-keepcolDT[keepcols==T]$cols

write.table(DT0[,keep,with=F],paste0(outputdata,'1707xx-180704-OutliersRemoved.txt'),sep='\t',row.names=F)
DT0<-fread(paste0(outputdata,'1707xx-180704-OutliersRemoved.txt'),sep='\t',header=T)



