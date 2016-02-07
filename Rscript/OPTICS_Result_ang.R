######################################################################
#
# Analyse OPTICS results (ang example)
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
names(d) = c("protein","aa","pos","phi1","phi2","theta","ss","a1phi1","a1psi","a2phi1","a2psi","a3phi1","a3psi","a4phi1","a4psi")
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
od[,'lastCl'] = rep(0)
# populate cluster table
for(clid in c(1:dim(cl)[1])){
  clname = paste('cl', clid, sep = '')
  start = cl[clid, 'start'] + 1
  end = cl[clid, 'end'] + 1
  od[c(start:end), clname] = 1
  od[c(start:end), 'lastCl'] = clid
}
# plot root cluster division
nRootCl = dim(cl[cl$parent == -1, ])[1]
od$col = rep(1)
for(i in c(1:nRootCl)){
  clname = paste('cl', i, sep = '')
  od$col = od$col + od[,clname]  * i
}
png(file = 'S3plot.png', width = 2000, height = 800)
par(mfrow = c(1,2))
s3 = scatterplot3d(od[,c('phi1', 'phi2', 'theta')], color = od$col, pch = 16)
s3 = scatterplot3d(od[,c('phi1', 'phi2', 'theta')], color = od$lastCl, pch = od$lastCl)
dev.off()
