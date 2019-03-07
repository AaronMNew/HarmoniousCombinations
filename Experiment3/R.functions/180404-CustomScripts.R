###################################################
### custom functions for analysis
###################################################

applyPaste<-function(x,collapse=''){
	apply(x,1,function(x)paste0(x,collapse=collapse))
}

t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,list=F,equal.variance=FALSE)
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
    if(list==T)dat<-list(m1-m2, se, t, 2*pt(-abs(t),df))
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
merge_recurse3<-function(qDT=NULL,dt.list,by.list,names.list=NULL){
	x<-qDT
	if(is.null(qDT))x<-dt.list[[1]]
	if(length(dt.list)==1)stop('you only have one data.table in the list')
	if(length(dt.list)!=(length(by.list)))stop('list lengths look funny')
	i<-1
	if(is.null(qDT))i<-2
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
	
	return(merge_recurse3(qDT=NULL,dt.list=summarylist,by.list=bylist2))
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

###################################################
### extra code
###################################################

notcols<-function(DT,x,grepl=T){
	xa<-x
	if(grepl==T)xa<-unlist(lapply(x,function(xa)colnames(DT)[grepl(xa,colnames(DT))]))
	colnames(DT)[!colnames(DT)%in%xa]
}

notin<-function(x,notx,grepl=T){
	out<-x[!x%in%notx]
	if(grepl==T){
		notcols<-unlist(lapply(notx,function(xa)x[grepl(xa,x)]))
		out<-x[!x%in%notcols]
		}
	out
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

pTechrep2<-function(DT,geno='clone.named',phenosh=NULL,plot=T){
	qDT<-copy(DT)
	groupings<-c(geno,'GAL3 v GAL80','GAL4 v GAL80','GAL3 v GAL4')
	phenocols<-c('growth.rate','phenotypic.index','yfp.mean.gal','yfp.mean.glu','fracon.gal','fracon.glu')
	if(is.null(phenosh)) phenocols <-phenosh <-c('GR','PI','fONglu','fONgal','ymglu','ymgal')
	if(!is.null(phenosh))phenocols<-phenosh
	
	qSub<-qDT[,c(phenocols,groupings),with=F]
	qSumm<-summary.func.all(qSub,c(phenocols),bylist= geno)
	mSub1<-melt(qSub,id=c(groupings))
	setnames(mSub1,c('variable','value'),c('phenotype','observed'))
	mSub<-melt(mSub1,id=c('phenotype','observed', geno))[value==T]
	mSub[,value:=NULL]
	setnames(mSub,'variable','pairing')
	
	mSumm1<-melt(qSumm,id= geno)
	mSumm1[,c('phenotype','stat'):=colsplit(variable,'_',c('phenotype','stat'))];mSumm1[,variable:=NULL]
	form<-paste0(paste(geno,'phenotype',sep='+'),'~stat')
	mSumm<-dcast(mSumm1,as.formula(form),value.var='value')
	mDT<-merge(mSub,mSumm,by=c(geno,'phenotype'))
	# check that all phenotypes are in there
	(identical(table(mDT$phenotype),table(mSub$phenotype)))
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
	pTechRep<-ggplot(mDT,aes(mean_scale,observed_scale))+geom_point(size=0.1,alpha=0.2,col='grey70')+
		geom_abline()+facet_grid(pairing2~phenotype2,scales='free')+
		ylab('observed measurements')+xlab('within-genotype mean of measurement')
	if(plot==T){return(pTechRep)}
	if(plot==F){return(mDT)}
}


subsetList<-function(lDT,commonCols,return.list,newnames=F,collapsefac='_'){
	lapply(1:length(lDT),function(i){
		x<-data.table(copy(lDT[[i]]))
		xa<-x[,c(return.list[[i]],commonCols),with=F]
		if(newnames==T)setnames(xa,commonCols,paste0(return.list[[i]],collapsefac,commonCols))
		xa
	})
}

CloneNamedSummary<-function(DT,allele=T){
	boo<-data.table(t(sapply(unique(DT$clone.named),function(x){
		xa<-strsplit(x,'\\ ')[[1]]
		lv<-grepl('WT',xa)
		mutcount<-3-sum(as.numeric(lv))
		single<-FALSE
		if(sum(as.numeric(lv))==2) single <-TRUE
		double<-FALSE
		if(sum(as.numeric(lv))==1) double <-TRUE
		WT<-FALSE
		if(sum(as.numeric(lv))==3)WT<-TRUE
		if(WT==T)single<-T;double<-T
		mut<-x
		if(WT!=T)mut<-paste0(c(xa[!lv]),collapse=' ')
		wt<-paste0(c(xa[lv]),collapse=' ')
		G3<-!lv[1]
		G80<-!lv[2]
		G4<-!lv[3]
		G3vG80<-as.logical(prod(G3|G80|WT))
		G3vG4<-as.logical(prod(G3|G4|WT))
		G4vG80<-as.logical(prod(G4|G80|WT))
		single_locus<-NA
		aa<-xa[3]
		bb<-xa[1]
		cc<-xa[2]
		if(single==T&WT!=T)single_locus<-paste0(strsplit(mut,'')[[1]][1:4],collapse='')
		if(grepl('8',single_locus)) single_locus<-paste0(strsplit(mut,'')[[1]][1:5],collapse='')
		list(x,mutcount,mut,wt,single,double,WT,WT,G3,G80,G4, G3vG80,G3vG4, G4vG80,single_locus,aa,bb,cc)
	})))
	if(allele==T){
		boop<-data.table(t(sapply(unique(DT$allele.named),function(x){
		xa<-strsplit(x,'\\ ')[[1]]
		c(allele.named=x,aa1=xa[3],bb1=xa[1],cc1=xa[2])
			})))
		}
	
	colnames(boo)<-c('clone.named','mutcount','mutants','WTs','single_lv','double_lv','WT','WT_lv','G3_mut_lv','G80_mut_lv','G4_mut_lv','GAL3 v GAL80','GAL3 v GAL4','GAL4 v GAL80','single_locus','aa','bb','cc')
	boo2<-boo[,lapply(.SD,unlist)]
	boo2[,G3.single:= (G3_mut_lv==T& mutcount==1)|WT==T]
	boo2[,G4.single:= (G4_mut_lv==1& mutcount==1)|WT==T]
	boo2[,G80.single:= (G80_mut_lv==1& mutcount==1)|WT==T]
	boo2[,mutant_gene:=unlist(sapply(mutants,function(x){
		xa<-strsplit(x,'\\ ')[[1]]
		bbbas<-sapply(xa,function(xaa){
			out<-paste0(strsplit(xaa,'')[[1]][1:4],collapse='')
			if(grepl('8',out))out<-paste0(out,'0',collapse='')
			out
		})
		paste0(bbbas,collapse=' ')
	}))]

	out<-merge(copy(DT),boo2,by='clone.named')
	if(allele==T)out<-merge(out,boop,by='allele.named') 
	return(out)
}


