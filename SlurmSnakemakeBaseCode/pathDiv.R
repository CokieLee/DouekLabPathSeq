args <- commandArgs(TRUE)
a<-args[1]
b<-args[2]
c<-args[3]
d<-args[4]
e<-args[5]
pathdiv=function(a,b,c,d,e)
{
	#a is output file from PathSeq
	#b is diversity output file name, ex. name.csv
	#c is distance matrix output file name, ex. name.csv
	#d is tree output name.pdf
	#e is heatmap distance output name.pdf
	
	Sys.setenv(libopen = "/data/vrc_his/douek_lab/projects/PathSeq/current_scripts/rlib/vegan/libs/libopenblas.so.0")

	library(vegan, lib.loc="/data/vrc_his/douek_lab/projects/PathSeq/current_scripts/rlib/")
	library(gplots, lib.loc="/data/vrc_his/douek_lab/projects/PathSeq/current_scripts/rlib/")
	hmcols<-colorRampPalette(c("white","green","green4","violet"))(100)
	
	getwd()
	
	data <- read.csv(a, row.names=1, check.names=F)
	data <- t(data)
	shannon_entropy=diversity(data,"shannon")
	richness=specnumber(data)
	evenness=shannon_entropy/log(richness)
	out=data.frame(shannon_entropy,richness,evenness)
	write.csv(out, file=b, quote=F)
	
	my.horn <- vegdist(data, method="horn")
	my.distance<-as.dist(my.horn)
	my.dist_matrix<-as.matrix(my.distance)
	write.csv(my.dist_matrix, c, quote=F)
	
	my.tree <- hclust(my.horn, method="complete")
	pdf(d, width=20, height=16.66)
	plot(my.tree)
	dev.off()
	
	pdf(e, width=20, height=20)
	heatmap.2(my.dist_matrix,col= hmcols,scale="none", key=T, keysize=1.5,density.info="none", trace="none",cexCol=1,cexRow=1,dendrogram="both",na.rm=T,symm=T,margins=c(10,10), Rowv=T)
	dev.off()
}
pathdiv(a,b,c,d,e)