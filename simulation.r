#SimulationCensoredEPD
rm(list =ls())
Packages<-c("ggplot2","parallel","R.utils","ReIns","grid","plotly")
for(i in 1:length(Packages)){
  if(all(installed.packages()[,1]!=Packages[i])){
    install.packages(Packages[i])}
  require(Packages[i], character.only = TRUE)
}

#Run subroutines for all estimators
source("cSubroutines.r")

#---------------------------------------------------
SimFun<-function(sim){ 
  set.seed(sim+34)
  u1<-runif(n); u2<-runif(n)
  # X<-(-log(u1))^(-gamma1)
  # Y<-(-log(u2))^(-1/beta)
  X<-(eta1*((1-u1)^(-1/lambda1)-1))^(1/tau1)
  Y<-(eta2*((1-u2)^(-1/lambda2)-1))^(1/tau2)
  Z<-cbind(X,Y)
  Z<-apply(Z,1,min)
  Z1<-order(Z)
  Z<-sort(Z)
  del<-(X<=Y)
  del1<-del[Z1]
  K<-bb:(k.max-1)
  Worms2<-CR_Worms(Z,(del1))$H
  rhos<-seq(-3,-0.5,by=0.1)
  LeurgBR_pen<-BRW(Z,(1-del1),L_kn=Worms2,w=1,pen=T,rho=rhos)$gamma[K,]
  LeurgBR<-BRW(Z,(1-del1),L_kn=Worms2,pen=F,rho=rhos)$gamma[K,]
  ro<-c(-2,-1.5,-1,-0.5)
  ind1<-which(rhos%in%ro)
  var_rho<-apply(LeurgBR,2,var,na.rm=T)
  ind2<-which.min(var_rho)
  var_rhoE<-apply(LeurgBR_pen,2,var,na.rm=T)
  ind2E<-which.min(var_rhoE)
  return(c(c(LeurgBR[,c(ind1,ind2)]),c(LeurgBR_pen[,c(ind1,ind2E)])))
}

 nsimul<-1000
 n<-500
 ######################
 eta1<-10
 tau1<-2
 lambda1<-2
 ########################
 eta2<-10
 tau2<-5
 lambda2<-2
 alpha<-tau1*lambda1
 gamma1<-1/alpha
 bb<-1#round(0.05*n,0)
 k.max<-n-1#round(n/2,0)
 nn<-length((k.max):bb)
 rhos<-c(-2,-1.5,-1,-0.5,"MinVar")

#Initiate Variables and matrices to be split
#----------------------------------------------
t.<-(2*length(rhos))
ML.part<-matrix(0,t.*(k.max-bb),nsimul)
SlvGrid<- matrix(seq(1:nsimul),1,nsimul)
#Run the SimFun in Parallel
#----------------------------------------------

library(parallel)
#************Checks what machine you're running on
if(Sys.info()[['nodename']]=="node0407"||Sys.info()[['nodename']]=="node0408"){#Remote Desktop
  nodelist<-ifelse(nsimul<64,nsimul,64)
}else if(Sys.info()[['sysname']]=="Windows"){#Normal Computer (leave 1 CPU free)
  nodelist <-ifelse(nsimul<12,nsimul,12)
} else{#Then you must be running on a Linux Cluster (MAC OS not allowed, sorry)
  nodelist <- unlist(c(read.table('nodelist.txt',sep="+")))}
cl <- makeCluster(nodelist)
clusterEvalQ(cl,{c(library(ReIns))})
clusterExport(cl, ls())
t<-system.time({ML.parts <- tryCatch({parSapplyLB(cl,SlvGrid,SimFun)},warning=function(w){print("Surpressed")})})
stopCluster(cl)

#Organize the returned variables coloumnwise
#-------------------------------------------------
K <- bb:(k.max-1)
ML.part <- ML.parts
bias <- matrix((rowMeans((ML.part-gamma1),na.rm=T)),(k.max-bb))
mseMeans <- matrix((rowMeans((ML.part-gamma1)^2,na.rm=T)),(k.max-bb))
alliEstimates <- data.frame(bias,mseMeans)

library(ggplot2)
library(grid)
library(reshape2)
EVI_pen<-data.frame(K,alliEstimates[1:(t./2)]);EVI_W<-data.frame(K,alliEstimates[(1+t./2):(t.)])
MSE_pen<-data.frame(K,alliEstimates[(t.+1):(3/2*t.)]);MSE_W<-data.frame(K,alliEstimates[(3/2*t.+1):(2*t.)])
EVI_pen_long<-melt(EVI_pen,id="K",variable.name="Estimators",value.name="Gamma");EVI_W_long<-melt(EVI_W,id="K",variable.name="Estimators",value.name="Gamma")
MSE_pen_long<-melt(MSE_pen,id="K",variable.name="Estimators",value.name="MSE");MSE_W_long<-melt(MSE_W,id="K",variable.name="Estimators",value.name="MSE")

MSE_pen_long$MSE<-sqrt(MSE_pen_long$MSE);MSE_W_long$MSE<-sqrt(MSE_W_long$MSE)

.scales_W<-.scales_Wp<-list()

LisRhos<-rhos#c(rhos[1],"Fraga")
for(j in 1:length(rhos)){
  .scales_W[[j]]<-substitute(paste(hat(bold(gamma)["1,k"])^"(s,W)","(",rho,"=",rho.x,")"), list(rho.x = LisRhos[j]))
  .scales_Wp[[j]]<-substitute(paste(hat(bold(gamma)["1,k"])^"*(BR,W)","(",rho,"=",rho.x,")"), list(rho.x = LisRhos[j])) 
}

cbbPalette <- c("#003300", "#0072B2", "#D55E00", "#CC79A7","#999999","#660000","#FF0000", "#56B4E9", "#009E73")

eviPlot_W<-ggplot(data=EVI_W_long,aes(x=K,y=Gamma,color=Estimators,linetype=Estimators,size=Estimators))+geom_line()+
  scale_linetype_manual("",labels=c(.scales_W),values=c(4,5,1,3,1,2,1))+
  scale_colour_manual("",labels=c(.scales_W),values=cbbPalette)+
  scale_size_manual("", labels=c(.scales_W), values=1+c(1.3,1.4,1.5,1.6,1.7,1.8,1.9))+
  xlab("K") + ylab("Bias") + ggtitle("Bias")+
  geom_hline(aes(yintercept=0))


eviPlot_Wp<-ggplot(data=EVI_pen_long,aes(x=K,y=Gamma,color=Estimators,linetype=Estimators,size=Estimators))+geom_line()+
  scale_linetype_manual("",labels=c(.scales_Wp),values=c(4,5,1,3,1,2,1))+
  scale_colour_manual("",labels=c(.scales_Wp),values=cbbPalette)+
  scale_size_manual("", labels=c(.scales_Wp), values=1+c(1.6,1.7,1.8,1.9,1.3,1.4,1.5))+
  xlab("K") + ylab("Bias") + ggtitle("Bias")+
  geom_hline(aes(yintercept=0))


msePlot_W<-ggplot(MSE_W_long,aes(x=K,y=MSE,color=Estimators,linetype=Estimators,size=Estimators))+geom_line()+
  scale_linetype_manual("",labels=c(.scales_W),values=c(4,5,1,3,1,2,1))+
  scale_colour_manual("",labels=c(.scales_W),values=cbbPalette)+
  scale_size_manual("", labels=c(.scales_W), values=1+c(1.6,1.7,1.8,1.9,1.3,1.4,1.5))+
  xlab("K") + ylab("RMSE") + ggtitle("Root Mean square error")+
  geom_hline(aes(yintercept=0))

msePlot_Wp<-ggplot(MSE_pen_long,aes(x=K,y=MSE,color=Estimators,linetype=Estimators,size=Estimators))+geom_line()+
  scale_linetype_manual("",labels=c(.scales_Wp),values=c(4,5,1,3,1,2,1))+
  scale_colour_manual("",labels=c(.scales_Wp),values=cbbPalette)+
  scale_size_manual("", labels=c(.scales_Wp), values=1+c(1.6,1.7,1.8,1.9,1.3,1.4,1.5))+
  xlab("K") + ylab("RMSE") + ggtitle("Root Mean square error")+
  geom_hline(aes(yintercept=0))


EHfeatures<-theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background=element_rect(fill=NA,colour='black',size=1,linetype='solid'))+
  theme(axis.line=element_line(colour="black",size=2))+       
  theme(legend.title = element_text(face = "bold",size=20),
        axis.text.y = element_text(angle = 0, hjust = 1, face="plain",size=15,color="black"),
        axis.text.x = element_text(angle = 0, hjust = 1, face="plain",size=15,color="black"),
        plot.title = element_text(hjust=0.5,lineheight=10, face="bold", color="black", size=40))+
  theme(axis.text=element_text(size=20),  
        axis.title=element_text(size=20,face="plain"))+
  theme(legend.text=element_text(size=20))+ theme(plot.margin=unit(c(0.5, 3, 0.5, 0.5), units="line"),
                                                  legend.text.align = 0,legend.key.width=unit(3,"line"),
                                                  legend.key.size = unit(4.5, "line"))

#Trim Coordinates
print((eviPlot_W+EHfeatures+scale_y_continuous(limit=c(-1,0.5))));
print((msePlot_W+EHfeatures+scale_y_continuous(limit=c(0,1))));

print((eviPlot_Wp+EHfeatures+scale_y_continuous(limit=c(-1,0.5))))
print((msePlot_Wp+EHfeatures+scale_y_continuous(limit=c(0,1))))


EHfeatures<-theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background=element_rect(fill=NA,colour='black',size=1,linetype='solid'))+
  theme(axis.line=element_line(colour="black",size=2))+       
  theme(legend.title = element_text(face = "bold",size=30),
        axis.text.y = element_text(angle = 0, hjust = 1, face="plain",size=30,color="black"),
        axis.text.x = element_text(angle = 0, hjust = 1, face="plain",size=30,color="black"),
        plot.title = element_text(hjust=0.5,lineheight=10, face="bold", color="black", size=45))+
  theme(axis.text=element_text(size=35),  
        axis.title=element_text(size=35,face="plain"))+
  theme(legend.text=element_text(size=35))+ theme(plot.margin=unit(c(0.6, 5, 2, 1), units="line"),
                                                  legend.text.align = 0,legend.key.width=unit(4,"line"),
                                                  legend.key.size = unit(4.5, "line"))


# pdf(file = paste("rhoBias_W.pdf",sep=""), width=18, height=10, paper="special",compress=T, family="Times",pointsize=20)
# print((eviPlot_W+EHfeatures+scale_y_continuous(limit=c(-1,0.5))));
# dev.off()
# 
# pdf(file = paste("rhoBias_Ws.pdf",sep=""), width=18, height=10, paper="special",compress=T, family="Times",pointsize = 20)
# print((eviPlot_Wp+EHfeatures+scale_y_continuous(limit=c(-1,0.5))))
# dev.off()
# 
# 
# pdf(file = paste("rhoRMSE_W.pdf",sep=""), width=18, height=10, paper="special",compress=T, family="Times",pointsize = 20)
# print((msePlot_W+EHfeatures+scale_y_continuous(limit=c(0,1))));
# dev.off()
# 
# pdf(file = paste("rhoRMSE_Ws.pdf",sep=""), width=18, height=10, paper="special",compress=T, family="Times",pointsize = 20)
# print((msePlot_Wp+EHfeatures+scale_y_continuous(limit=c(0,1))))
# dev.off()

#save.image(file="rhos.RData")
