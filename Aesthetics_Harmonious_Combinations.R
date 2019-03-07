# aesthetics.dir<-'/Users/anew/Google Drive/002_lehnerlabData/RCustomFuctions/181119-Aesthetics_Harmonious_Combinations.R'
# source(aesthetics.dir)

###################################################
### custom ggplot themes for harmonious combinations paper
###################################################

factorDefault<-c('orange1','steelblue','darkgreen','red','chartreuse3')
factorDefault2<-scale_colour_brewer(palette = "Set1")
factorTwoSugar<-c('orangered','deepskyblue')
heatmapColRange<-c('cornflowerblue','yellow','indianred')
heatmapColBounded<-c('white','cornflowerblue','yellow','indianred')
heatmapColCentered<-c('cornflowerblue','white','indianred')

cust_font_facet<-theme(strip.text.x = element_text(size = 6),strip.text.y = element_text(size = 6))
cust_font_axis<-theme(axis.text = element_text(size = 6))
cust_font_axislab<-theme(axis.title = element_text(size = 8))
cust_font_legend<-theme(legend.text=element_text(size=6),legend.title=element_text(size=8))

cust_font<-cust_font_facet+cust_font_axis+cust_font_axislab+ cust_font_legend

cust_panel<-theme(panel.background=element_blank(),panel.border=element_blank())
cust_strip<-theme(strip.background=element_rect(fill='grey90',colour='transparent'))
cust_axis_arrow<-theme(axis.line=element_line(arrow=arrow(angle=20,length=unit(0.15,'cm')),
	size=0.3,colour='grey50'),
	axis.ticks=element_blank())
cust_axis_noarrow<-theme(axis.line=element_line(size=0.3,colour='grey50'),axis.ticks=element_blank())

theme_NewPub<-theme_classic()+cust_strip+ cust_axis_arrow + cust_font+ cust_panel
theme_NewPub_noarrow<-theme_classic()+cust_strip+ cust_axis_noarrow + cust_font+ cust_panel

# booasdfawfgb<-data.table(x=1:10,y=1:10,fac=factor(rep(c('A','B'),5)))

# ggplot(booasdfawfgb,aes(x,y))+geom_point()+geom_line(aes(col=fac))+
	# theme_NewPub +
	# scale_colour_manual(values=factorDefault)+
	# facet_grid(~fac)

