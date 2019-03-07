
###################################################
### calculate expectations for binary using regression of log transformed and output predicted values
###################################################

applyPaste<-function(x,collapse=''){
	apply(x,1,function(x)paste0(x,collapse=collapse))
}
modelAnyLevel<-function(DT,phenotypex=NULL,genos,interactionlevel=1,backg=0,backsd=0,sig=0.05,genosalt=NULL){
#	DT<-copy(DAT2); interactionlevel =2;backg=0; phenotypex=NULL; backsd=0
	if(!is.null(phenotypex))setnames(DT,phenotypex,'x')
	genosx<-quote(mget(genos))
	N<-length(genos)
	indep<-paste('(',paste(genos,collapse='+'),')')
	ff<-as.formula(paste('log(abs(x-backg)) ~',indep))
	if(interactionlevel>1)ff<-as.formula(paste('log(abs(x-backg)) ~',paste(indep,'^', interactionlevel,collapse='')))
	credmin=backg+backsd
	mergecols<-genos
	if(!is.null(genosalt))mergecols<-c(genos,genosalt)
	mod1<-lm(ff,DT)
	summ<-DT[,summary.func(x),by=mergecols ]
	summ[,order:=apply(data.frame(eval(genosx)),1,function(x)sum(as.numeric(x)))]
	pred<-predict(mod1,summ,se.fit=T)
	summ[,muhat_mean:=exp(pred$fit)+backg]
	summ[,muhat_se:=(exp(pred$fit)*pred$se.fit)]	
	summ[,muhat_df:=pred$df]	
	summ[,muhat_resiscale:=pred$residual.scale]
	summ$interactionlevel<-interactionlevel
	summ[,epistasis_mean:=mean-muhat_mean]
	summ[,epistasis_se:=sqrt(se^2+ muhat_se^2)]
	summ[,epistasis_p.val:=t.test3(m1=mean,m2= muhat_mean,se=epistasis_se,df=(muhat_df-1))$p.value]
	summ[,epistasis_p.adj:=p.adjust(epistasis_p.val,'fdr')]
	summ[,epistasis_sig:=epistasis_p.adj<=0.05]
	summ[mean<credmin&muhat_mean<credmin,epistasis_sig:=FALSE]
	return(data.table(summ))
}

phenotype<-'growth.rate'
phenotypes<-c('growth.rate','phenotypic.index','fracon.glu','fracon.gal','yfp.mean.glu','yfp.mean.gal')


DTx0<-data.table(ldply(lapply(1:length(phenotypes),function(i){
	phenotype<-phenotypes[i]
	
	### binary
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
	x1<-data.table(ldply(lapply(1:N,function(interactionlevel){
		modelAnyLevel(DAT2,genos=genos,interactionlevel=interactionlevel,backg=DAT2[A==0,mean(x)], backsd =DAT2[A==0,sd(x)],genosalt=genosalt)
	})))
	x1[,var:='binary genotypes']
	
	p1<-ggplot(x1,aes(muhat_mean,mean,xmin=muhat_mean-muhat_se*1.96,xmax=muhat_mean+muhat_se*1.96,ymin=mean-se*1.96,ymax=mean+se*1.96))+geom_abline()+
			geom_errorbar()+geom_errorbarh()+geom_point(aes(colour=epistasis_sig))+facet_grid(~interactionlevel)
	
	
	### total copies
	DAT0 <-qDT[,c(phenotype,'totalcopies','plasgeno'),with=F]
	DAT0$arb<-1:nrow(DAT0)
	a<-(t(unlist(sapply(DAT0 $totalcopies,function(x)strsplit(x,'\\ ')[[1]]))))
	b<-data.table(a);b$totalcopies<-rownames(a);b$arb<-1:nrow(b)
	clonenamedcols<-c('GAL3','GAL80','GAL4')
	colnames(b)<-c(clonenamedcols,'totalcopies','arb')
	N0<-sum(apply(b[,clonenamedcols,with=F],2,function(x)length(unique(x))-1))
	DAT0x<-merge(DAT0,b,by=c('totalcopies','arb'))
	
	genocols<-c('totalcopies','GAL4','GAL3','GAL80')
	genosalt<-c('plasgeno')
	DAT1<-DAT0x[,c(phenotype,genocols,genosalt),with=F]
	N<-length(genocols)-1
	GENOS<-paste0(LETTERS[1:N])
	genos<-paste0(letters[1:N])
	GENOSx<-quote(mget(GENOS))
	genosx<-quote(mget(genos))
	setnames(DAT1,colnames(DAT1),c('x','genotypei',GENOS,genosalt))
	dum<-DAT1[,GENOS,with=F][,lapply(.SD,function(x)as.factor(x))];colnames(dum)<-genos
	DAT2<-data.table(DAT1,dum)
	x2<-data.table(ldply(lapply(1:N0,function(interactionlevel){
		modelAnyLevel(DAT2,genos=genos,interactionlevel=interactionlevel,backg=DAT2[A==0,mean(x)], backsd =DAT2[A==0,sd(x)],genosalt=genosalt)
	})))
	x2[,var:='totalcopies genotypes']
	
	p2<-ggplot(x2,aes(muhat_mean,mean,xmin=muhat_mean-muhat_se*1.96,xmax=muhat_mean+muhat_se*1.96,ymin=mean-se*1.96,ymax=mean+se*1.96))+geom_abline()+
			geom_errorbar()+geom_errorbarh()+geom_point(aes(colour=epistasis_sig))+facet_grid(~interactionlevel)
	
	
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
	x3<-data.table(ldply(lapply(1:N0,function(interactionlevel){
		modelAnyLevel(DAT2,genos=genos,interactionlevel=interactionlevel,backg=DAT2[A==0&B==0,mean(x)],backsd=DAT2[A==0&B==0,sd(x)],sig=0.1,genosalt='genotypei')
	})))
	setnames(x3,'genotypei',genosalt)
	x3[,var:='all genotypes']
	p3<-ggplot(x3,aes(muhat_mean,mean,xmin=muhat_mean-muhat_se*1.96,xmax=muhat_mean+muhat_se*1.96,ymin=mean-se*1.96,ymax=mean+se*1.96))+geom_abline()+
			geom_errorbar()+geom_errorbarh()+geom_point(aes(colour=epistasis_sig))+facet_grid(~interactionlevel,scales='free')
	
	meltcols<-c('plasgeno','var','interactionlevel',grepincols(x1,c('se','mean','epistasis')))
	xdum<-x1[,meltcols,with=F]
	DTx<-rbind(x1[,meltcols,with=F],x2[,meltcols,with=F],x3[,meltcols,with=F])
	DTx[,intlev:=paste0('interaction order: ',interactionlevel)]
	DTx[,var:=factor(var,levels= unique(DTx$var))]
	DTx$phenotype<-phenotype

	# ## NOT RUN
	# # example of plot to make from each DTx
	# alph<-1
	# p4<-ggplot(DTx,aes(muhat_mean,mean,xmin=muhat_mean-muhat_se*1.96,xmax=muhat_mean+muhat_se*1.96,ymin=mean-se*1.96,ymax=mean+se*1.96))+geom_abline()+
		# geom_errorbar(alpha=alph,colour='grey80')+geom_errorbarh(alpha=alph,colour='grey80')+geom_point(aes(colour=epistasis_sig),shape=21,size=3,stroke=1)+
		# scale_colour_manual(values=c('cornflowerblue','indianred'))+
		# xlab(paste0('prediction based on linear model of log-transformed ',phenotype,' values'))+ylab(paste0('measured ',phenotype,' for all plasmid and genome genotypes'))+
		# facet_grid(var~ intlev,scales='free')		
	# h<-10;w<-h*2.7
	# ggsave(paste0('.',figout,'1707xx-17092x-180518-',phenotype,'-AllOrdersOfInteractionsAcrossGenotypeSimplifications.png'),p4,height=h,width=w)

	DTx
	
})))




DTx0[,obs_scaled_mean:=zscore(mean),by='phenotype']
DTx0[,exp_scaled_mean:=zscore(muhat_mean),by='phenotype']
DTx0[,obs_scaled_se:= obs_scaled_mean/mean*se,by='phenotype']
DTx0[,exp_scaled_se:=exp_scaled_mean/muhat_mean*muhat_se,by='phenotype']
DTx0[,phenotype:=factor(phenotype,levels=phenotypes)]

alph<-1
p4<-ggplot(DTx0[phenotype=='growth.rate'],aes(exp_scaled_mean, obs_scaled_mean,xmin= exp_scaled_mean-exp_scaled_se*1.96,xmax= exp_scaled_mean + exp_scaled_se*1.96,ymin= obs_scaled_mean-obs_scaled_se*1.96,ymax= obs_scaled_mean + obs_scaled_se*1.96))+geom_abline()+
	geom_errorbar(alpha=alph,colour='grey80')+geom_errorbarh(alpha=alph,colour='grey80')+geom_point(aes(colour=epistasis_sig),shape=21,size=3,stroke=1)+
	scale_colour_manual(values=c('cornflowerblue','indianred'))+ggtitle(phenotype)+
	xlab(paste0('scaled prediction based on linear model of log-transformed phenotype values'))+ylab(paste0('scaled phenotype, measured for all plasmid and genome genotypes'))+
	facet_grid(phenotype+var~intlev,scales='free')		
l<-6
h<-l*5;w<-l*2.7
ggsave(paste0('.',figout,'1707xx-17092x-180518-GrowthRateScaled-AllOrdersOfInteractionsAcrossGenotypeSimplifications.png'),p4,height=h,width=w)


p5<-ggplot(DTx0,aes(exp_scaled_mean, obs_scaled_mean,xmin= exp_scaled_mean-exp_scaled_se*1.96,xmax= exp_scaled_mean + exp_scaled_se*1.96,ymin= obs_scaled_mean-obs_scaled_se*1.96,ymax= obs_scaled_mean + obs_scaled_se*1.96))+geom_abline()+
	geom_errorbar(alpha=alph,colour='grey80')+geom_errorbarh(alpha=alph,colour='grey80')+geom_point(aes(colour=epistasis_sig),shape=21,size=3,stroke=1)+
	scale_colour_manual(values=c('cornflowerblue','indianred'))+ggtitle(phenotype)+
	xlab(paste0('scaled prediction based on linear model of log-transformed phenotype values'))+ylab(paste0('scaled phenotype, measured for all plasmid and genome genotypes'))+
	facet_grid(phenotype+var~intlev,scales='free')		
l<-6
h<-l*5;w<-l*2.7
ggsave(paste0('.',figout,'1707xx-17092x-180518-AllPhenotypesScaled-AllOrdersOfInteractionsAcrossGenotypeSimplifications.png'),p5,height=h,width=w)


p6a<-ggplot(DTx0[phenotype=='growth.rate'],aes(muhat_mean,mean,xmin=muhat_mean-muhat_se*1.96,xmax=muhat_mean+muhat_se*1.96,ymin=mean-se*1.96,ymax=mean+se*1.96))+geom_abline()+
	geom_errorbar(alpha=alph,colour='grey80')+geom_errorbarh(alpha=alph,colour='grey80')+geom_point(aes(colour=epistasis_sig),shape=21,size=3,stroke=1)+
	scale_colour_manual(values=c('cornflowerblue','indianred'))+ggtitle(phenotype)+
	xlab(paste0('scaled prediction based on linear model of log-transformed phenotype values'))+ylab(paste0('scaled phenotype, measured for all plasmid and genome genotypes'))+
	facet_grid( phenotype+var~intlev ,scales='free')		
l<-6
w<-l*5;h<-l*2.7

p6<-ggplot(DTx0,aes(muhat_mean,mean,xmin=muhat_mean-muhat_se*1.96,xmax=muhat_mean+muhat_se*1.96,ymin=mean-se*1.96,ymax=mean+se*1.96))+geom_abline()+
	geom_errorbar(alpha=alph,colour='grey80')+geom_errorbarh(alpha=alph,colour='grey80')+geom_point(aes(colour=epistasis_sig),shape=21,size=3,stroke=1)+
	scale_colour_manual(values=c('cornflowerblue','indianred'))+ggtitle(phenotype)+
	xlab(paste0('scaled prediction based on linear model of log-transformed phenotype values'))+ylab(paste0('scaled phenotype, measured for all plasmid and genome genotypes'))+
	facet_grid(intlev~ phenotype+var ,scales='free')		
l<-6
w<-l*5;h<-l*2.7
ggsave(paste0('.',figout,'1707xx-17092x-180518-AllPhenotypes-AllOrdersOfInteractionsAcrossGenotypeSimplifications.png'),p6,height=h,width=w)

###################################################
### use regsubsets to identify the best parameters for linear regression
###################################################
library(leaps)

### binary genotypes

phenotype<-'growth.rate'
phenotypes<-c('growth.rate','phenotypic.index','fracon.gal','yfp.mean.glu','yfp.mean.gal')
i<-1
phenotype<-phenotypes[i]
### binary
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


DATfull<-copy(DAT2)
credmin<-backg+backsd*1.96
DATfull[,corx:=x]
DATfull[corx<credmin,corx:=backg+0.0000001]
DATfull[,subx:=abs(x-backg)]
DATfull[,lx:=log(subx)]

DATfull[subx<0,]

mergecols<-c(genos,'plasgeno')

best.subset<-regsubsets(lx~a*b*c,DATfull,really.big=T,nvmax=8)
#best.subset<-regsubsets(lx~A*B*C*D*E*F,DATfull,really.big=T)

best.summary<-summary(best.subset)

get_model_formula <- function(id, object, outcome){
  # get models data
  models <- summary(object)$which[id,-1]
  # Get outcome variable
  #form <- as.formula(object$call[[2]])
  #outcome <- all.vars(form)[1]
  # Get model predictors
  predictors <- names(which(models == TRUE))
  predictors <- paste(predictors, collapse = "+")
  # Build model formula
  form<-(gsub('0','',gsub('1','',paste0(outcome, "~", predictors))))
  as.formula(form)
}
Ni<-nrow(best.summary$which)
forms<-lapply(1:Ni,function(i){
	get_model_formula(i, best.subset, 'lx')
})
forms[[Ni+1]]<-as.formula(paste0('lx ~ a + b + c + b:c'))

mods<-lapply(forms,function(form){
	mod1<-lm(form,DATfull)
	list(summary(mod1),form)
})

i<-1
modframes<-lapply(forms,function(form){
	print(form)
#	form<-get_model_formula(i[[1]], best.subset, 'lx')
	mod1<-lm(form,DATfull)
	summ<-DATfull[,summary.func(exp(lx)+backg),by=mergecols]
	summ[,lx:=mean]
	pred<-predict(mod1,summ, se.fit=T)
	summ[,c('fit_mean','fit_se'):=list(exp(pred$fit)+backg,(exp(pred$fit)*pred$se.fit))]
	summ[,fit_df:=pred$df]
	summ[,epistasis_mean:=mean-fit_mean]
	summ[,epistasis_se:=sqrt(se^2+ fit_se^2)]
	summ[,epistasis_p.val:=t.test3(m1=mean,m2= fit_mean,se=epistasis_se,df=(fit_df-1))$p.value]
	summ[,epistasis_p.adj:=p.adjust(epistasis_p.val,'fdr')]
	summ[,epistasis_sig:=epistasis_p.adj<=0.05]
	summ[,modtit:=capture.output(form)[[1]]]
	return(summ)
})
modF<-data.table(ldply(modframes))

ggplot(modF,aes(x=fit_mean,y=mean,ymax=upper,ymin=lower,xmin=fit_mean-fit_se,xmax=fit_mean+fit_se))+geom_abline(col='grey80')+geom_errorbar()+geom_errorbarh()+geom_point(aes(colour=epistasis_sig))+
	facet_wrap(~ modtit)



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
DATfull<-copy(DAT2)
credmin<-backg+backsd*1.96
DATfull[,corx:=x]
DATfull[corx<credmin,corx:=backg+0.0000001]
DATfull[,subx:=abs(x-backg)]
DATfull[,lx:=log(subx)]

DATfull[subx<0,]

mergecols<-c(genos,'plasgeno')

best.subset<-regsubsets(lx~a*b*c*d*e*f,DATfull,really.big=T)
#best.subset<-regsubsets(lx~A*B*C*D*E*F,DATfull,really.big=T)

plot(best.subset)
best.summary<-summary(best.subset)

get_model_formula <- function(id, object, outcome){
  # get models data
  models <- summary(object)$which[id,-1]
  # Get outcome variable
  #form <- as.formula(object$call[[2]])
  #outcome <- all.vars(form)[1]
  # Get model predictors
  predictors <- names(which(models == TRUE))
  predictors <- paste(predictors, collapse = "+")
  # Build model formula
  form<-(gsub('0','',gsub('1','',paste0(outcome, "~", predictors))))
  as.formula(form)
}
Ni<-nrow(best.summary$which)
forms<-lapply(1:Ni,function(i){
	get_model_formula(i, best.subset, 'lx')
})
forms[[Ni+1]]<-as.formula(paste0('lx ~ a:b + c:d + c:d:e:f'))
#forms[[Ni+1]]<-as.formula(paste0('lx ~ A:B + C:D + C:D:E:F'))

mods<-lapply(forms,function(form){
	mod1<-lm(form,DATfull)
	list(summary(mod1),form)
})

i<-1
modframes<-lapply(forms,function(form){
	print(form)
#	form<-get_model_formula(i[[1]], best.subset, 'lx')
	mod1<-lm(form,DATfull)
	summ<-DATfull[,summary.func(x),by=mergecols]
	summ[,lx:=mean]
	pred<-predict(mod1,summ, se.fit=T)
	summ[,c('fit_mean','fit_se'):=list(exp(pred$fit)+backg,(exp(pred$fit)*pred$se.fit))]
	summ[,fit_df:=pred$df]
	summ[,epistasis_mean:=mean-fit_mean]
	summ[,epistasis_se:=sqrt(se^2+ fit_se^2)]
	summ[,epistasis_p.val:=t.test3(m1=mean,m2= fit_mean,se=epistasis_se,df=(fit_df-1))$p.value]
	summ[,epistasis_p.adj:=p.adjust(epistasis_p.val,'fdr')]
	summ[,epistasis_sig:=epistasis_p.adj<=0.05]
	summ[,modtit:=capture.output(form)[[1]]]
	return(summ)
})
modF<-data.table(ldply(modframes))

ggplot(modF,aes(x=fit_mean,y=mean,ymax=upper,ymin=lower,xmin=fit_mean-fit_se,xmax=fit_mean+fit_se))+geom_abline(col='grey80')+geom_errorbar()+geom_errorbarh()+geom_point(aes(colour=epistasis_sig))+
	facet_wrap(~ modtit)


modframe<-data.frame(model=ldply(lapply(mods,function(x)paste0(x[[2]])))$V3,rsq=ldply(lapply(mods,function(x)x[[1]]$r.squared)))


###################################################
### exhaustively plot searching all models for binary genotypes to show fit as function of different (interaction) parameters
###################################################
library(leaps)

### binary genotypes

phenotype<-'growth.rate'
phenotypes<-c('growth.rate','phenotypic.index','fracon.gal','yfp.mean.glu','yfp.mean.gal')
i<-1
phenotype<-phenotypes[i]
### binary
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


DATfull<-copy(DAT2)
credmin<-backg+backsd*1.96
DATfull[,corx:=x]
DATfull[corx<credmin,corx:=backg+0.0000001]
DATfull[,subx:=abs(x-backg)]
DATfull[,lx:=log(subx)]

summ<-DATfull[,summary.func(subx),by=c(genos)][order(a,b,c)]
mmat<-(model.matrix(mean~a*b*c,summ))
mmatlv<-as.data.frame(matrix(as.logical(mmat),ncol=length(colnames(mmat)))[-1,-1])
# modvar<-gsub('1','',colnames(mmat)[-1])
modvar<-c('a','b','c','a*b','b*c','a*c','a*b*c')


list1<-(sapply(1:length(modvar),function(x)combn(modvar,x)))

forms<-unlist(lapply(list1,function(x){
	a<-t(x)
	b<-apply(a,1,function(x)paste0(x,collapse='+'))
	c(paste0('lx~',b))
	}))
formlength<-unlist(lapply(list1,function(x){
	a<-(t(x))
	(rep(ncol(a),nrow(a)))
	}))
x<-forms[127]
formorders<-lapply(forms,function(x){
	a<-strsplit(x,'\\~')[[1]][2]
	b<-strsplit(a,'\\+')[[1]]
	df<-t(data.frame(table(sapply(b,function(xa){
		d<-1
		if(grepl('\\*',xa))d<-length(strsplit(xa,'\\*')[[1]])
		if(grepl('\\:',xa))d<-length(strsplit(xa,'\\:')[[1]])
		d
	}))))
	out<-data.table(t(df[-1,]))
	colnames(out)<-paste0('order_',df[1,])
	out
	})
formframe<-cbind(modtit =forms,formlength=formlength,rbindlist(formorders,fill=T))

i<-1
modframes<-lapply(1:nrow(formframe),function(i){
	form<-paste0(formframe[i,1])
	print(form)
#	form<-get_model_formula(i[[1]], best.subset, 'lx')
	mod1<-lm(form,DATfull)
	summ<-DATfull[,summary.func(x),by=mergecols]
	summ[,lx:=mean]
	pred<-predict(mod1,summ, se.fit=T)
	summ[,c('fit_mean','fit_se'):=list(exp(pred$fit)+backg,(exp(pred$fit)*pred$se.fit))]
	summ[,fit_df:=pred$df]
	summ[,epistasis_mean:=mean-fit_mean]
	summ[,epistasis_se:=sqrt(se^2+ fit_se^2)]
	summ[,epistasis_p.val:=t.test3(m1=mean,m2= fit_mean,se=epistasis_se,df=(fit_df-1))$p.value]
	summ[,epistasis_p.adj:=p.adjust(epistasis_p.val,'fdr')]
	summ[,epistasis_sig:=epistasis_p.adj<=0.05]
	#summ[,modtit:=capture.output(form)[[1]]]
	cbind(summ,formframe[i])
	
})
moduniqs<-lapply(modframes,function(x){
	paste0(as.character(round(x$fit_mean,6)),collapse='')
})
uniqlv<-!duplicated(moduniqs)
modF<-data.table(ldply(modframes[uniqlv]))
modF[,order:=apply(data.frame(a,b,c),1,function(x)sum(as.numeric(x)))]
modF[,unique(order),by='modtit']
applyPaste<-function(x,collapse=''){
	apply(x,1,function(x)paste0(x,collapse=collapse))
}


modF[,modtit:=factor(modtit,levels=unique(modF$modtit))]
modF[,modtitnull:=as.character(modtit)]

modF[,modDim:=unlist(lapply((modtitnull),function(x){
	dim(model.matrix(as.formula(x),DATfull))[2]
	})),by='modtit']
modF[,VarExp:=round(varexplained(mean,fit_mean),2),by='modtit']
modF[,Rsq:=round(summary(lm(mean~fit_mean))$r.squared,2),by='modtit']
modF[,modtit:=applyPaste(data.frame(modtit,modDim,Rsq,VarExp),collapse='\n')]
modtitLevels1<-modF[,lapply(.SD,unique),by='modtit',.SDcols=c('modtitnull','modDim','Rsq','VarExp')][order(modDim,Rsq,VarExp,modtitnull)]$modtit
modtitLevels2<-modF[,lapply(.SD,unique),by='modtit',.SDcols=c('modtitnull','modDim','Rsq','VarExp')][order(modDim,Rsq,VarExp,modtitnull)]$modtitnull
modF[,modtit:=factor(modtit,levels=modtitLevels1)]
modF[,modtitnull:=factor(modtitnull,levels=modtitLevels2)]
p0<-ggplot(modF,aes(x=fit_mean,y=mean,ymax=upper,ymin=lower,xmin=fit_mean-fit_se,xmax=fit_mean+fit_se))+geom_abline(col='grey80')+geom_errorbar()+geom_errorbarh()+
	geom_point(aes(colour=epistasis_sig),shape=21,stroke=0.5)+scale_colour_manual(values=c('cornflowerblue','indianred'))+
	facet_grid(order~ modtit)+ggtitle('predicted vs observed for all binary genotype linear models stratified by order of genotype\nfacet titles are the model, the number of terms in the model, the rsq and the variance explained')+theme(plot.title=element_text(face='plain'))

mergecols<-c('a','b','c','plasgeno','modtit','modtitnull','formlength',grepincols(modF,'order_'))
meanEpis1<-modF[,(sqrt(mean(epistasis_mean^2))),by=mergecols]
meanEpis1[,abc:=applyPaste(data.frame(a,b,c))]
meanEpis1[,RMSE:=V1]
meanEpis2<-modF[,epistasis_mean,by=mergecols]
meanEpis<-merge(meanEpis1,meanEpis2,by=mergecols)

p1<-ggplot(meanEpis,aes(x=abc,y=modtitnull))+geom_tile(aes(fill=RMSE))+scale_fill_gradientn(colours=c('white','cornflowerblue','yellow','indianred'))+
	theme(axis.text.y=element_text(size=8),plot.title=element_text(face='plain'))+
	xlab('genotype [GAL4 GAL3 GAL80] where 0 = WT and 1 = deletion mutation\n011 is the GAL3.delta GAL80.delta GAL4.WT strain, which displays pathway epistasis')+
	ylab('models, increasing by complexity and explanatory power')+
	ggtitle('Root mean squared error for all possible linear models for binary genotype set stratified by binary genotype')

p2<-ggplot(meanEpis,aes(x=plasgeno,y=modtitnull))+geom_tile(aes(fill=RMSE))+scale_fill_gradientn(colours=c('white','cornflowerblue','yellow','indianred'))+
	theme(axis.text.y=element_text(size=8),plot.title=element_text(face='plain'),axis.text.x=element_text(angle=30,size=10,hjust=1))+
	xlab('genotype [plasmid[GAL3 GAL80 GAL4]genome[GAL3 GAL80 GAL4]] where 1 = WT and 0 = deletion mutation\n001001 001000 and 000001 are the GAL3.delta GAL80.delta GAL4.WT strain, which displays pathway epistasis')+
	ylab('models, increasing by complexity and explanatory power')+
	ggtitle('Root mean squared error for all possible linear models for binary genotype set stratified by full genotype')

p3<-ggplot(meanEpis,aes(x=abc,y=modtitnull))+geom_tile(aes(fill= epistasis_mean))+scale_fill_gradientn(colours=c('cornflowerblue','white','indianred'),limits=c(-0.2,0.2))+
	theme(axis.text.y=element_text(size=8),plot.title=element_text(face='plain'))+
	xlab('genotype [GAL4 GAL3 GAL80] where 0 = WT and 1 = deletion mutation\n011 is the GAL3.delta GAL80.delta GAL4.WT strain, which displays pathway epistasis')+
	ylab('models, increasing by complexity and explanatory power')+
	ggtitle('Epistasis for all possible linear models for binary genotype set stratified by binary genotype')

p4<-ggplot(meanEpis,aes(x=plasgeno,y=modtitnull))+geom_tile(aes(fill=epistasis_mean))+scale_fill_gradientn(colours=c('cornflowerblue','white','indianred'),limits=c(-0.2,0.2))+
	theme(axis.text.y=element_text(size=8),plot.title=element_text(face='plain'),axis.text.x=element_text(angle=30,size=10,hjust=1))+
	xlab('genotype [plasmid[GAL3 GAL80 GAL4]genome[GAL3 GAL80 GAL4]] where 1 = WT and 0 = deletion mutation\n001001 001000 and 000001 are the GAL3.delta GAL80.delta GAL4.WT strain, which displays pathway epistasis')+
	ylab('models, increasing by complexity and explanatory power')+
	ggtitle('Epistasis for all possible linear models for binary genotype set stratified by full genotype')

setwd(head.dir)
newdir<-paste0('.',figout,'1707xx-17092x-180606-AllLinearModelsBinaryGenotypePlots')
if(!dir.exists(newdir))system(paste0('mkdir ', newdir))
setwd(newdir)
h<-10;w<-h*2
plotlist<-list(p0,p1,p2,p3,p4)
lapply(seq(length(plotlist)),function(i){
	ggsave(paste0('180606-plot-',i,'.pdf'), plotlist[[i]],height=h,width=w)
})
setwd(head.dir)