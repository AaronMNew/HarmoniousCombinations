
all.ttest<-function(DTx,pheno='pheno',bylist,plot=F,points=T,padj='fdr',sig.cutoff=0.05){
	# DTx<-copy(phenDT)
	# pheno <- 'growth.rate'
	# bylist<- c('aa','bb','cc','dd')
	DT<-copy(DTx)
	setnames(DT,pheno,'pheno')
	summ<-DT[,summary.func(pheno),by=bylist]
	df<-data.frame(summ[,bylist,with=F])
	summ[,Var1:=applyPaste(df,' ')]
	summ[,Var2:=applyPaste(df,' ')]
	DT1<-data.table(expand.grid(summ$Var1,summ$Var2))
	DT2<-merge(summ[,!'Var2'],merge(DT1,summ[,!'Var1'],by='Var2'),by='Var1')
	temp<-data.frame(dcast(DT2,formula=as.formula('Var1~Var2'),value.var='mean.x'))[,-1]
	tri_lv_frame<-data.table(melt(lower.tri(as.matrix(temp))))
	uniqfac<-unique(c(summ $Var1, summ $Var2));uniqfac<-uniqfac[order(uniqfac)]
	DT3<-merge(data.table(X2=1:length(uniqfac),Var2=uniqfac),merge(data.table(X1=1:length(uniqfac),Var1=uniqfac), DT2,by='Var1'),by='Var2')
	DT4<-merge(DT3,tri_lv_frame,by=c('X1','X2'))[value==T][,!c('X1','X2','value')]
	DT4[,diff:=mean.x-mean.y]
	DT4[,df:=(N.x+N.y)-2]
	DT4[,se:=sqrt(sd.x^2+sd.y^2)/sqrt(df+2)]
	DT5 <-DT4[,.(m1=diff,m2=0,se=se,df=df)]
	DT6 <-data.table(DT4, DT5[,t.test3(m1,m2,se,df),by=seq_len(nrow(DT5))][,!c('diff.of.means','std.err','seq_len')])
	DT6[,p.adj:=p.adjust(p.value,padj)]
	DT6[,sig:=p.adj<=sig.cutoff]

	out<-DT6
	if(plot==F)return(data.table(out))else{
		lims<-c(-max(abs(range(out $diff))),max(abs(range(out $diff))))
	
		p<-ggplot(out,aes(Var1,Var2,fill=diff))+geom_tile()+
			scale_fill_gradientn(colours=c('cornflowerblue','white','indianred'),limits=lims)+
			theme(axis.text.x=element_text(angle=30,hjust=01))+
			theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))
		fdrname<-paste0('p adj < ', sig.cutoff)
		if(points==T)p<-p+geom_point(aes(col=sig))+scale_colour_manual(values=c('orange1','darkgreen'))+labs(value= paste0('ttest ',fdrname))
		return(list(data.table(out),p))
	}

}



all.ttest<-function(DTx,pheno='pheno',bylist,plot=F,points=T,padj='fdr',sig.cutoff=0.05){
	# DTx<-copy(phenDT)
	# pheno <- 'growth.rate'
	# bylist<- c('aa','bb','cc','dd')
	DT<-copy(DTx)
	setnames(DT,pheno,'pheno')
	summ<-DT[,summary.func(pheno),by=bylist]
	df<-data.frame(summ[,bylist,with=F])
	summ[,Var1:=applyPaste(df,' ')]
	summ[,Var2:=applyPaste(df,' ')]
	DT1<-data.table(expand.grid(summ$Var1,summ$Var2))
	DT2<-merge(summ[,!'Var2'],(merge(DT1,summ[,!'Var1'],by='Var2')),by='Var1')
	dum1<-data.table(Var1=summ$Var1,Var2=summ$Var2,X1=1:nrow(summ),X2=1:nrow(summ))
	DT3<-merge(dum1[,!c('Var1','X1')],(merge(DT2,dum1[,!c('X2','Var2')],by=c('Var1')) ),by='Var2')
	temp<-data.frame(dcast(DT3,formula=as.formula('Var1~Var2'),value.var='mean.x'))
	tri_lv_frame<-data.table(melt(lower.tri(as.matrix(temp[,-1]))))
	DT4<-merge(DT3,tri_lv_frame,by=c('X1','X2'))[value==T][,!c('X1','X2','value')]	
	
	DT4[,diff:=mean.x-mean.y]
	DT4[,df:=(N.x+N.y)-2]
	DT4[,se:=sqrt(sd.x^2+sd.y^2)/sqrt(df+2)]
	DT5<-DT4[,.(m1=diff,m2=0,se=se,df=df)]
	DT6<-data.table(DT4,DT5[,t.test3(m1,m2,se,df),by=seq_len(nrow(DT5))][,!c('diff.of.means','std.err','seq_len')])

	DT7<-copy(DT6)
	DT7[,c('Var1','Var2','diff','t'):=list(Var2,Var1,-diff,-t)]
	DT8<-rbind(DT7,DT6)

	# take lower triangle of matrix to avoid redundant information. this expects that melt returns X1 and X2 as arbitrary column names
	tri_lv_frame<-melt(lower.tri(as.matrix(dcast(DT8,formula=as.formula('Var1~Var2'),value.var='N.x'))))
	# make a dummy frame to merge the lower triangle with later
	uniqfac<-unique(c(DT8 $Var1, DT8 $Var2));uniqfac<-uniqfac[order(uniqfac)]
	dum<-merge(data.table(X2=1:length(uniqfac),Var2=uniqfac),merge(data.table(X1=1:length(uniqfac),Var1=uniqfac),DT8,by='Var1'),by='Var2')
	# merge logical vector frame with data and drop out the extra columns used for merging and filtering triangle
	DT9<-merge(dum,tri_lv_frame,by=c('X1','X2'))[value==T][,!c('X1','X2','value')]
	DT9[,p.adj:=p.adjust(p.value,padj)]
	DT9[,sig:=p.adj<=sig.cutoff]

	if(plot==F)return(data.table(DT6))else{
		lims<-c(-max(abs(range(DT6$diff))),max(abs(range(DT6$diff))))
	
		p<-ggplot(DT9,aes(Var1,Var2,fill=diff))+geom_tile()+
			scale_fill_gradientn(colours=c('cornflowerblue','white','indianred'),limits=lims)+
			theme(axis.text.x=element_text(angle=30,hjust=01))+
			theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))
		fdrname<-paste0('p adj < ', sig.cutoff)
		if(points==T)p<-p+geom_point(aes(col=sig))+scale_colour_manual(values=c('orange1','darkgreen'))+labs(value= paste0('ttest ',fdrname))
		return(list(data.table(DT6),p))
	}

}


test<-all.ttest(phenDT,pheno='growth.rate',bylist=c('aa','bb','cc','dd'),plot=T,points=T)
test2<-test[[1]][aa.x=='GAL4.WT'&aa.y=='GAL4.WT'&cc.x=='GAL80.07'&cc.y=='GAL80.07'&dd.x!='HIS5.Sch_pom'&dd.y!='HIS5.Sch_pom']
test3<-test[[1]][aa.x=='GAL4.WT'&aa.y=='GAL4.WT'&cc.x=='GAL80.07'&cc.y=='GAL80.07'&dd.x!='HIS5.Sch_pom'&dd.y!='HIS5.Sch_pom'&dd.x!='GALK.Can_alb'&dd.y!='GALK.Can_alb']
duplicated(abs(test[[1]][aa.x=='GAL4.WT'&aa.y=='GAL4.WT'&cc.x=='GAL80.07'&cc.y=='GAL80.07'&dd.x!='HIS5.Sch_pomb'&dd.y!='HIS5.Sch_pomb']$diff))







TukeyInteractionFlip2<-function(DTx,plot=T,points=T,tukeyPoint=T,degfreedom=2,fdr=0.05){
	# improved over first implementation to make more sense and give a nice ggplot plot every time.
	# points = T means you put points over the tiles for significance
	# tukeyPoint = T means the colour of the point is from the p adj from the Tukey test
	fdrname<-paste0('p adj < ',fdr)
	# need to get all pairs and reciprocal differences together
	DTemp1<-na.exclude(data.table(DTx,data.table(colsplit(rownames(DTx),'-',c('A','B')))[,lapply(.SD,as.character)]))
	DTemp2<-copy(DTemp1)
	DTemp2[,c('A','B','diff','lwr','upr'):=list(B,A,-diff,-lwr,-upr)]
	DT<-rbind(DTemp1,DTemp2)
	setnames(DT,'diff','diff A-B')
	DT[,sig:=`p adj`<0.05]
	
	# take lower triangle of matrix to avoid redundant information. this expects that melt returns X1 and X2 as arbitrary column names
	tri_lv_frame<-melt(lower.tri(as.matrix(dcast(DT,formula=as.formula('A~B'),value.var='sig'))))
	# make a dummy frame to merge the lower triangle with later
	uniqfac<-unique(c(DT $B, DT $A));uniqfac<-uniqfac[order(uniqfac)]
	dum<-merge(data.table(X2=1:length(uniqfac),B=uniqfac),merge(data.table(X1=1:length(uniqfac),A=uniqfac),DT,by='A'),by='B')
	# merge logical vector frame with data and drop out the extra columns used for merging and filtering triangle
	out1<-merge(dum,tri_lv_frame,by=c('X1','X2'))[value==T][,!c('X1','X2','value')]
	out1[,sig1:=sig]
	setnames(out1,'sig1',paste0('Tukey ',fdrname))
	ttests<-data.table(t(out1[,apply(data.table(m1= abs(`diff A-B`),m2=0,se=abs(`diff A-B`-lwr),df= degfreedom),1,function(x)t.test3(m1=x[1],m2=x[2],se=x[3],df=x[4]))]))
	ttests[,p.adj.fdr:=p.adjust(V4,'fdr')]
	ttests[,c('t','ttest.pval'):=list(as.numeric(V3,V4))]
	ttests[,sig2:=p.adj.fdr<=0.05]
	out<-data.table(out1,ttests[,.(t, ttest.pval, p.adj.fdr,sig2)])
	setnames(ttests,'sig2',paste0('ttest ',fdrname))
	# calculate limits for plotting colours
	lims<-c(-max(abs(range(out $`diff A-B`))),max(abs(range(out$`diff A-B`))))
	# make plot
	p<-ggplot(out,aes(A,B,fill=`diff A-B`))+geom_tile()+
		scale_fill_gradientn(colours=c('cornflowerblue','white','indianred'),limits=lims)+
		theme(axis.text.x=element_text(angle=30,hjust=01))+
		theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))
	# output
	if(points==T& tukeyPoint ==T)p<-p+geom_point(aes(col=sig))+scale_colour_manual(values=c('orange1','darkgreen'))+labs(value= paste0('Tukey ',fdrname))
	if(points==T& tukeyPoint ==F)p<-p+geom_point(aes(col=sig2))+scale_colour_manual(values=c('orange1','darkgreen'))+labs(color= paste0('ttest ',fdrname))
	if(plot==T)return(p)
	if(plot==F)return(out[,!'sig'])
}
