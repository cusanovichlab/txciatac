args = commandArgs(TRUE)
report = read.table(args[1],header=T)
indextable = args[2]
output = args[3]
if(indextable != "None"){
  indices = read.table(indextable)
}

#This code is directly adapted from the picard source code: http://picard-tools.sourcearchive.com/documentation/1.25-1/classnet_1_1sf_1_1picard_1_1sam_1_1DuplicationMetrics_a0a76ff106a01cb46966cbdf0d778843.html#a0a76ff106a01cb46966cbdf0d778843
fer = function(x,c,n){
  return(c/x - 1 + exp(-n/x))
}

picardcomp = function(c = 50,n = 100,m = 1, M = 10){
  if(c >= n || fer(m*c,c,n) <= 0){
    print("Invalid inputs!")
    break
  }
  while(fer(M*c,c,n) >= 0){
    M = M*10
  }
  
  for(i in 1:40){
    r = (m+M)/2
    u = fer(r*c,c,n)
    if(u == 0){
      break
    }else{
      if(u > 0){
        m = r
      }else{
        if(u < 0){
          M = r
        }
      }
    }
  }
  return(c*(m+M)/2)
}

deduped = which(report[,2] > report[,3])
totalfrags = apply(report[deduped,2:3],1,function(x){picardcomp(x[2],x[1])})
dereport = cbind(report[deduped,],totalfrags)
write.table(dereport,paste0(output,".complexity_report.txt"),row.names=F,sep="\t",quote=F)

if(indextable != "None"){
  commonind = intersect(dereport[,1],indices[,1])
  cleanindices = indices[match(commonind,indices[,1]),]
  cleanreport = dereport[match(commonind,dereport[,1]),]
}

png(paste0(output,".dedup_plot.png"), width = 600, height = 1500)
par(mfrow= c(3, 1))
plot((report[,2]-report[,3])/report[,2],log10(report[,3]),xlab="Fraction PCR Duplicate",ylab="Log10(Unique Reads)",pch=20)
plot(report[deduped,3]/totalfrags,log10(report[deduped,3]),xlab="Fraction Unique Molecules Observed",ylab="Log10(Unique Reads)",pch=20)
if(indextable != "None"){
  groupmedians = aggregate(cleanreport$totalfrags, by=list(cleanindices[,2]), FUN=median)
  par(mar=c(12,3,2,2))
  par(cex.axis=0.5)
  if(nrow(groupmedians) == 1){
    boxplot(log10(cleanreport$totalfrags) ~ cleanindices[,2],xlab=paste0(groupmedians[,1],"\n(Median = ",round(groupmedians[,2],1),")"),las=2,cex.names=0.5)
  }else{
    boxplot(log10(cleanreport$totalfrags) ~ cleanindices[,2],names=paste0(groupmedians[,1],"\n(Median = ",round(groupmedians[,2],1),")"),las=2,cex.names=0.5)
  }
}
dev.off()
