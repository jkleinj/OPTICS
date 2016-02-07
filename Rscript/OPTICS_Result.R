######################################################################
#
# Analyse OPTICS results (xyz example)
#
#######################################################################
# load packages
library(scatterplot3d)
# Load data
d = read.table("input.dat")
o = read.table("output.dat", header = T)
cl = read.table("cluster.dat", header = T)
#######################################################################
# name columns
names(d) = c('x','y','z')
# generate reachability plot
png(file = 'RDplot.png', width = 1000, height = 800)
layout(matrix(1:2, nrow=2, byrow=TRUE), heights = c(0.8,0.2))
barplot(o[,'RD'], border = NA)
plot.new()
for(clid in c(1:dim(cl)[1])){
  y = 1 - (cl[clid, 'parent'] + 1) / (max(cl[,'parent']) + 1)
  x1 = (cl[clid, 'start'] + 1) / dim(d)[1]
  x2 = (cl[clid, 'end'] + 1) / dim(d)[1]  
  segments(x1,y,x2,y, lwd = 3, col = 'blue')
}
dev.off()

# create cluster ordered table
od = d[o$dataId + 1,]
for(i in c(1:dim(cl)[1])){
  od[,paste('cl', i,sep = '')] = rep(0)
}
# populate cluster table
for(clid in c(1:dim(cl)[1])){
  clname = paste('cl', clid, sep = '')
  start = cl[clid, 'start'] + 1
  end = cl[clid, 'end'] + 1
  od[c(start:end), clname] = 1
}
# plot root cluster division
nRootCl = dim(cl[cl$parent == -1, ])[1]
od$col = rep(1)
for(i in c(1:nRootCl)){
  clname = paste('cl', i, sep = '')
  od$col = od$col + od[,clname]  * i
}
png(file = 'S3plot.png', width = 1000, height = 800)
s3 = scatterplot3d(od[,c(1:3)], color = od$col, pch = 16)
dev.off()
