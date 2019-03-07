
library(plyr)
library(XML)
library(data.table)


head.dir<-"/Users/anew/Google Drive/002_lehnerlabData/180404-VarGALK.With.G80xG4_2"
date<-'180404' # beginning date clones were measured
plates.dir<-paste(head.dir,'/ExperimentalPlates',sep='')
figout<-'./figure.output/'
outputdata<-'./output.data/'
layoutdir<-'./layout.sample.data/'
setwd(head.dir)




# # # code below parses a plate
# data <- xmlParse('Plate_014.xml')
# xml_data <- xmlToList(data)

# # data for A1
# x<-xml_data$experiment$tray$specimen[3]
# test<-NULL

# names(x$tube$loader_user_info)
# # below returns a table with useful information about the experiment
# test<-c("sample_flow_rate","sample_volume","record_sample_volume","mixing_volume","mixing_speed","number_mixes", "wash_volume",'tube.name')
# for(i in 3:87){
	# x<-xml_data$experiment$tray$specimen[i]
	# a<-unlist(c(as.vector(x$tube$loader_user_info),as.vector(x$tube$.attrs)))
	# test<-rbind(test,a)
# }


# extend the individual plate parsing above into entire experiments
template.directory<-'/Templates'
xml.dir.index<-list.files(paste(head.dir,template.directory,sep=''))
xmls<-xml.dir.index[grep('.xml',xml.dir.index)]
output1<-lapply(xmls,function(x){
	xab<-paste(head.dir,template.directory,'/',x,sep='')
	data2 <- xmlParse(xab) 
	xml_data2 <- xmlToList(data2)
	#as.data.frame(names(xml_data2$experiment)[9:(length(xml_data2$experiment)-3)]) # check you have identified the "trays" which are the plates
	test<-c("sample_flow_rate","sample_volume","record_sample_volume","mixing_volume","mixing_speed","number_mixes", "wash_volume",'tube.name','plate.num')
	plates<-xml_data2$experiment[9:(length(xml_data2$experiment)-3)]
	for (ix in 1:length(plates)){
		# define plate
		plate<-plates[ix]
		# define the range of tubes included in this plate
		len<-(length(plate$tray$specimen)-1)
		# get each well of the experiment into a specimens object you cycle through extracting information
		specimens<-plate$tray$specimen[3:len]
		print(length(specimens))
		for(i in 1:length(specimens)){
			# create tube-specific object
			x<-specimens[i]
			# the numbers below for date.plate are tailored for the naming of the different plates 
			# (in this case, a 6-digit date, a dash, then the standard Plate_XXX)
			date.plate<-substr(xab,(nchar(xab)-22),(nchar(xab)-23+6))
			# pull out the information you want
			a<-unlist(c(as.vector(x$tube$loader_user_info),as.vector(x$tube$.attrs), as.vector(plate$tray$.attrs[1]),date.plate=date.plate))
			# append it onto the test dummy variable
			test<-rbind(test,a)
		}
	}
	rownames(test)<-NULL
	output<-as.data.frame(test[2:nrow(test),])

	colnames(output)<-test[1,]
	nchar(xab)
	
	
	output$samp_id<-paste(date.plate, '-',output$plate.num,'.',output$tube.name,sep='')
	output
})
lapply(output1,function(d)nrow(d)/96)
lapply(output1,function(d)table(d$plate.num))

sampling.data <- data.table(ldply(output1, data.table))
sampling.data$sample_flow_rate.1<-NULL
setwd(head.dir)
write.table(sampling.data,paste0(layoutdir,date,'-SamplingData.txt'),row.names=FALSE,sep='\t')
fread(paste0(layoutdir,date,'-SamplingData.txt'),header=T,sep='\t')




