library(mclust)
args = commandArgs(TRUE)

#Prefix of current experiment
currexpt = args[1]
#cellfloor = No. reads to be considered a cell
#Either specify a number or "mclust" to calculate automatically
cellfloor = args[2]

report2 = read.table(paste0(currexpt,".report.txt"),header=T)
bkgd.ind = grep("bkgd",report2$Tag)
if(length(bkgd.ind) > 0){
	nobkgdmat = report2[-bkgd.ind,]
} else {
	nobkgdmat = report2
}
cellcall = Mclust(data.frame(log10(nobkgdmat$Total)),G=2)
if(cellfloor == "mclust"){
  cellfloor = min(nobkgdmat[which(cellcall$classification == 2 & cellcall$uncertainty < 0.05),2])
}else{
  cellfloor = as.numeric(cellfloor)
}
#cellfloor = 1000

bloom_collision <- function(n1, n2, n12){
  #from https://github.com/jbloomlab/multiplet_freq
  n <- n1 * n2 / n12
  mu1 <- -log((n - n1) / n)
  mu2 <- -log((n - n2) / n)
  mu <- mu1 + mu2
  return (1 - mu * exp(-mu) / (1 - exp(-mu)))
}

subsetmat = nobkgdmat[which(nobkgdmat$Total >= cellfloor),]
subsetmice = subsetmat[which(subsetmat[,4]/subsetmat[,2] >= 0.9),]
subsethumans = subsetmat[which(subsetmat[,3]/subsetmat[,2] >= 0.9),]
outtable = cbind(as.character(rownames(subsetmat)),as.character(subsetmat[,1]))
outtable = cbind(outtable,rep("Collision",times=dim(outtable)[1]))
outtable[which(subsetmat[,3]/subsetmat[,2] >= 0.9),3] = "Human"
outtable[which(subsetmat[,4]/subsetmat[,2] >= 0.9),3] = "Mouse"
write.table(outtable,paste0(currexpt,".readdepth.cells.indextable.txt"),row.names=F,col.names=F,sep="\t",quote=F)

pdf(paste0(currexpt,".results.pdf"),height=8,width=8)
par(mfrow=c(2,2))
subsamples = levels(nobkgdmat$Tag)

for(i in 1:length(subsamples)){
  if(subsamples[i] == "bkgd"){next}
  write(subsamples[i], stdout())
  currind = grep(subsamples[i],	nobkgdmat$Tag)
  currsubind = grep(subsamples[i], subsetmat$Tag)
  currsub = subset(nobkgdmat,Tag == subsamples[i])
  currsubcells = which(currsub$Total >= cellfloor)
  if(length(currsubcells) == 0){next}
  currcellcall = Mclust(data.frame(log10(currsub$Total)),G=2)
  currcellfloor = min(currsub[which(currcellcall$classification == 2 & currcellcall$uncertainty < 0.05),2])
  hist(log10(currsub$Total)[currsub$Total > 4],breaks=60,col="mediumseagreen",main=subsamples[i],
       xlab="Number of Reads (log10)",las=1,xlim=c(0,7))
  abline(v=log10(cellfloor),lty="dashed",  lwd=2,col="blue")
  abline(v=log10(currcellfloor),lty="dashed",  lwd=2,col="red")
  legend("topright",legend = c(paste0("global ",cellfloor),paste0("current ",currcellfloor)),fill=c("blue","red"))
  boxplot(list(subset(subsethumans,Tag == subsamples[i])[,5]/subset(subsethumans,Tag == subsamples[i])[,3],subset(subsetmice,Tag == subsamples[i])[,6]/subset(subsetmice,Tag == subsamples[i])[,4]),
        ylim=c(0,1),ylab="Fraction of Reads Mapping to DHS",col="dodgerblue2",lwd=2,pch=20,las=1)
  abline(h=0.5,lwd=2,lty="dashed")
  boxplot(list(subset(subsethumans,Tag == subsamples[i])[,5]/subset(subsethumans,Tag == subsamples[i])[,3],subset(subsetmice,Tag == subsamples[i])[,6]/subset(subsetmice,Tag == subsamples[i])[,4]),
        col="dodgerblue2",lwd=2,pch=20,add=T,las=1)
}
dev.off()

pdf(paste0(currexpt,".results.hists.pdf"),height=12,width=12)
par(mfrow=c(2,2))
for(i in 1:length(subsamples)){
  if(subsamples[i] == "bkgd"){next}
  currind = grep(subsamples[i], nobkgdmat$Tag)
  currsubind = grep(subsamples[i], subsetmat$Tag)
  currsub = subset(nobkgdmat,Tag == subsamples[i])
  currsubcells = which(currsub$Total >= cellfloor)
  if(length(currsubcells) == 0){next}
  hist(log10(currsub$Total)[currsub$Total > 4],breaks=60,col="mediumseagreen",main=subsamples[i],
       xlab="Number of Reads (log10)",las=1,xlim=c(0,7))
  abline(v=log10(cellfloor),lty="dashed",  lwd=2)
  legend("topright",c(paste0("Total Reads: ",sum(currsub$Total)),
                    paste0("\n Total Reads (cells only): ",sum(currsub[currsubcells,2])),
                    paste0("\n Total Barcodes: ",length(currsub$Total)),
                    paste0("\n Number of Cells: ",length(subsetmat$Total[currsubcells])),
                    paste0("\n Median Reads/Cell: ",median(currsub[currsubcells,2])),
                    paste0("\n Range of Reads/Cell: ",min(currsub[currsubcells,2])," - ",max(currsub[currsubcells,2]))),bty="n")
  mousecounts = length(which(currsub[currsubcells,4]/currsub[currsubcells,2] >= 0.9))
  humancounts = length(which(currsub[currsubcells,3]/currsub[currsubcells,2] >= 0.9))
  totalcounts = length(currsubcells)
  humanfrac = humancounts/(humancounts + mousecounts)
  mousefrac = mousecounts/(humancounts + mousecounts)
  collisioninflation = ((humanfrac^2) + (mousefrac^2))/(2*humanfrac*mousefrac)
  plot(currsub[currsubcells,3],currsub[currsubcells,4],pch=20,xlab="Human reads",
       ylab="Mouse reads")
  legend("topright",c(paste0("Total Cells: ",length(currsubcells)),
                      paste0("Human Cells: ",humancounts),
                      paste0("Mouse Cells: ",mousecounts),
                      paste0("Observed Collision Rate: ",signif(1-(humancounts + mousecounts)/length(currsubcells),4)),
                      paste0("Calculated Collision Rate: ",signif(2*collisioninflation*(1-(humancounts + mousecounts)/length(currsubcells)),4)),
                      paste0("Bloom Collision Rate: ",signif(bloom_collision(totalcounts-mousecounts,totalcounts-humancounts,totalcounts-(humancounts+mousecounts)),4))))
}

hist(log10(nobkgdmat$Total)[nobkgdmat$Total > 4],breaks=60,col="mediumseagreen",main="Overall",
	xlab="Number of Reads (log10)",las=1,xlim=c(0,7))
abline(v=log10(cellfloor),lty="dashed",  lwd=2)
legend("topright",c(paste0("Total Reads: ",sum(nobkgdmat$Total)),
	paste0("Total Reads (cells only): ",sum(subsetmat$Total)),
	paste0("\n Total Barcodes: ",length(nobkgdmat$Total)),
	paste0("\n Number of Cells: ",length(subsetmat$Total)),
	paste0("\n Median Reads/Cell: ",median(subsetmat[,2])),
	paste0("\n Range of Reads/Cell: ",min(subsetmat[,2])," - ",max(subsetmat[,2]))),bty="n")
mousecounts = length(which(subsetmat[,4]/subsetmat[,2] >= 0.9))
humancounts = length(which(subsetmat[,3]/subsetmat[,2] >= 0.9))
totalcounts = nrow(subsetmat)
humanfrac = humancounts/(humancounts + mousecounts)
mousefrac = mousecounts/(humancounts + mousecounts)
collisioninflation = ((humanfrac^2) + (mousefrac^2))/(2*humanfrac*mousefrac)
plot(subsetmat[,3],subsetmat[,4],pch=20,xlab="Human reads",
     ylab="Mouse reads")
legend("topright",c(paste0("Total Cells: ",totalcounts),
                    paste0("Human Cells: ",humancounts),
                    paste0("Mouse Cells: ",mousecounts),
                    paste0("Observed Collision Rate: ",signif(1-(humancounts + mousecounts)/totalcounts,4)),
                    paste0("Calculated Collision Rate: ",signif(2*collisioninflation*(1-(humancounts + mousecounts)/totalcounts),4)),
                    paste0("Bloom Collision Rate: ",signif(bloom_collision(totalcounts-mousecounts,totalcounts-humancounts,totalcounts-(humancounts+mousecounts)),4))))
dev.off()
