

# Scripts for dealing with flow cytometry data.
source.code.path.file<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/CustomFunctions.R'

head.dir<-'/Users/anew/Google Drive/002_lehnerlabData/1707xx-RepeatAnalysis'
layout.dir<-'./layout.sample.data'
date<-'1707xx' # beginning date clones were measured
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/'
outputdata<-'./output.data/'
setwd(head.dir)

br<-seq(-0.5,2.7,length=61)

###################################################
### custom functions for analysis
###################################################
source(source.code.path.file)
load.all() # this function loads all the standard packages you use in this experiment
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
	x[,row:=factor(r,levels=unique(x$r)[8:1])]
	x[,column:=factor(c,levels=1:12)]
	setnames(x,phenot,'pheno')
	ggplot(x,aes(column,row))+geom_tile(aes(fill=abs(pheno)))+scale_fill_gradientn(colours=c('cornflowerblue','yellow','indianred'),na.value='transparent',limits=collim)+
		theme(legend.title=element_blank())+scale_colour_manual(values=c('black'),na.value='transparent')+
		facet_grid(measurement.date.gal~plate.gal)+ggtitle(paste0(phenot," ",unique(x$source.plate)))	
}





###################################################
### load data
###################################################

setwd(head.dir)
# analyzed data have already been written file and here are fread output written the first time around.
output<-fread(paste0(outputdata,'1707xx-180704-analyzed.tab'),header=TRUE,sep='\t')
nullsamps.summary<-fread(paste0(outputdata,'1707xx-180704-null.samps.tab'),header=TRUE,sep='\t')
# layout is a table with headers and at least a "samp_id" column either generated in the layout.maker script or by hand in excel.
layout<-fread(paste0(layout.dir,'/',date,'-180704-layout.txt',sep=''),header=TRUE,sep='\t')

if(!is.null(nullsamps.summary)){
	#create a dummy data frame for null samples generated and add to output
	nullsamps.layout <-as.data.frame(t(apply(nullsamps.summary,1,function(x){
		c(rep(NA,(ncol(output)-1)),x)
	})))
	colnames(nullsamps.layout)<-colnames(output)
	output.dummy<-rbind(output,nullsamps.layout,make.row.names=FALSE)

	output.final<-data.table(t(apply(output.dummy[,1:(ncol(output.dummy)-1)],1,as.numeric)),output.dummy$samp_id)
	colnames(output.final)<-colnames(output)
	}else{output.final<-output}

output.final[,samp_id:=as.character(samp_id)]
layout[,samp_id:=as.character(samp_id)]
dens.cols<-as.character(colnames(output.final)[c(grep("\\]",colnames(output.final)))])
stripfun<-function(x)mean(as.numeric(strsplit(gsub('\\]','',gsub('\\(','',x)),',')[[1]]))
new.names<-as.character(paste0('dens__',as.numeric(as.character(unlist(sapply(dens.cols,stripfun))))))
setnames(output.final,dens.cols,new.names)
dens.cols<-new.names
merged<-merge(output.final,layout,by='samp_id')



###################################################
### add a few things: compute number of generations + phenotypic index
###################################################
dens.norm<-data.table(t(apply(merged[, dens.cols,with=FALSE],1,function(x)round(row.normalization(x),3))))
merged[,dens.cols]<-dens.norm

merged[,c('measurement.date','samp.plate.rc'):=colsplit(samp_id,'-',c('measurement.date','samp.plate.rc'))]
merged[,allele.named:=applyPaste(data.frame(GAL3.allele,GAL80.allele,GAL4.allele),' ')]
merged[,replicate.day:=applyPaste(data.frame(allele.named,measurement.date),' ')]
merged[,clone.named:=applyPaste(data.frame(GAL3.clone.named,GAL80.clone.named,GAL4.clone.named),' ')]

# for cell density per ml, cells.per.sec * 1/µl per second sampling rate * 1000 µl / ml
dim(merged)

merged[,cell.per.ml:=cells.per.sec*(1/as.numeric(sample_flow_rate))*1000]
DT.glu<-merged[glu==1]
DT.gal<-merged[gal==1]

# to couple direct measurements with each other (glucose > galactose) i need to create an id variable that will allow joining glu with gal measurements
DT.glu[,next.day:=as.numeric(measurement.date)+1]
DT.gal[,next.day:=NA]
DT.glu[,id:=applyPaste(data.frame(allele.named,next.day,plate.name.in.source.fcs),' ')]
DT.gal[,id:=applyPaste(data.frame(allele.named,measurement.date,plate.name.in.source.fcs),' ')]
DT.gal[,gal:=NULL]
DT.gal[,glu:=NULL]
DT.glu[,gal:=NULL]
DT.glu[,glu:=NULL]

#variables definition
dil.factor<-9/150
hours <-12
galphenos<-colnames(output)
# i re-merge them here and then calculate generations, growth rate and phenotypic index
phenos1<-c(colnames(output.final),"record_sample_volume","mixing_volume","mixing_speed","number_mixes","wash_volume","tube.name", "plate.name.in.source.fcs",'keep.plate', "biol.rep", "mask.scoring",
	'replicate.day','cell.per.ml','next.day','measurement.date','plate','comments')
setnames(DT.gal,phenos1,paste0(phenos1,'.gal'))
setnames(DT.glu,phenos1,paste0(phenos1,'.glu'))
DT.merge1<-merge(DT.glu,DT.gal,by=c('id'),suffixes=c('','___y'),all=TRUE)
DT.merge<-DT.merge1[,!grepincols(DT.merge1,'___'),with=F]
DT.merge[,generations:=log((cell.per.ml.gal/((cell.per.ml.glu)*dil.factor)),base=2)]
DT.merge[,growth.rate:=log((cell.per.ml.gal/((cell.per.ml.glu)*dil.factor)))/hours]
DT.merge[,phenotypic.index:=fracon.glu+fracon.gal]
mDT<-melt(DT.merge[,c('generations','growth.rate','phenotypic.index')])
ggplot(mDT,aes(value))+geom_histogram()+facet_grid(~variable,scales='free')




###################################################
### deal with unpaired glu > gal observations
###################################################

DT0x<-copy(DT.merge)
DT0x[,comments.gal:=paste0(comments.gal)]
DT0x[,comments.glu:=paste0(comments.glu)]
# arose from defective 'id' variable calculation from a repeat set of expts that happened in August
gal_unmerged<-DT0x[is.na(plate.glu)]# DT0x[is.na(plate.glu),c(grepincols(DT0x,'.gal'),'comments'),with=F]
# OK this GAL_unmerged is unpaired because the glucose plate was dropped.
glu_unmerged<-DT0x[is.na(plate.gal)]# DT0x[is.na(plate.gal),c(grepincols(DT0x,'.glu'),'comments'),with=F]
# OK 170707-Plate_014's pair 170708-Plate_014 was not included because The source folder and .fcs files for this plate have been lost.
# 170717-Plate_010 (glu plate) is unpaired in 24 wells because of problems with they cytometer's acquisition of the GAL samples in rows D and E
# 170824-Plate_001's gal plate counterpart was dropped
# 170824-Plate_008's gal counterpart was mis-inoculated
DT0x[is.na(plate.glu),comments.gal:=paste0('unmerged because glu plate was dropped')]
DT0x[plate.glu=='170707-Plate_014',comments.gal:='The source folder and .fcs files for this plate have been lost.']
table(DT0x[plate.glu=='170717-Plate_010',comments.gal])
table(DT0x[plate.glu=='170717-Plate_010',comments.gal])

########
### strangely duplicated some plates and samples
########
# DTdedup<-DT0x[,lapply(.SD,unique),by=c('samp_id.glu','samp_id.gal')]
DT0x[samp_id.gal%in% DT0x[duplicated(samp_id.gal)]$samp_id.gal]
DT0x[samp_id.glu%in% DT0x[duplicated(samp_id.glu)]$samp_id.glu]

DT0x[,vartemp:=applyPaste(data.frame(samp_id.glu,samp_id.gal))]

DT0x<-DT0x[!duplicated(vartemp)]

write.table(DT0x,paste0(outputdata,'1707xx-180704-preprocessed_data.txt'),sep='\t',row.names=F)
DT0<-fread(paste0(outputdata,'1707xx-180704-preprocessed_data.txt'),sep='\t',header=T)

