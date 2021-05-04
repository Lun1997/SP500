library(splines)
library(ggplot2)
library(tidyr)
str1=function(string){
  index=gsub("-","/",string) 
  return(index)
}
str2=function(string){
  index=gregexpr("-",string)[[1]][1]+1
  index2=substring(string,index,nchar(string))
  index3=sub("-","/",index2) 
  return(index3)
}
#data====
data=read.csv("建模資料.csv")
date=data$X[-1] #日期對應VFINX
data$VFINX[1:nrow(data)-1]=data$VFINX[2:nrow(data)]
data=data[1:nrow(data)-1,]
data=data[,-1]
head(data)
#NA:沒有遺失值====
sum(is.na(data))

#B-spline====
knot=function(x,lowbond,upperbond){
  if(x==0){return(NULL)}
  else{return(seq(lowbond,upperbond,by=(upperbond-lowbond)/(x+1))[c(-1,-(x+2))])}
}

#spline
order=4
var_num=ncol(data)-1
Y=data$VFINX

#超參數取得:節點組合
#交叉驗證:LOOCV
kn=15:30
knotset=expand.grid(kn,kn,kn)
#迴圈====
# ptm <- proc.time()
# RSSCV=c()
# for(N in 1:nrow(knotset)){
#   basis=list()
#   #基底矩陣
#   for(var in 1:var_num){
#     x=data[,var]
#     low=min(x)
#     up=max(x)
#     B=bs(x,deg=order-1, knots=knot(knotset[N,var],low,up), Boundary.knots=c(low,up),intercept=F)
#     basis[[var]]=B
#   }
#   #基底調整
#   basismix=do.call(cbind,basis)
#   u=apply(basismix,2,mean)
#   for(i in 1:ncol(basismix)){
#     basismix[,i]=basismix[,i]+u[i]
#   }
#   basismix=cbind(rep(1,nrow(data)),basismix,Y)
#   colnames(basismix)=c(paste("V",0:(ncol(basismix)-2),sep=""),"y")
#   basismix=as.data.frame(basismix)
#   #model
#   model=lm(y~.-1,data=basismix)
#   hii.v=lm.influence(model)$hat
#   rsscv=mean((model$resid/(1-hii.v))^2)
#   RSSCV[N]=rsscv
# }
# proc.time() - ptm
#平行運算====
library(parallel)
myCPU=detectCores()
cl=makeCluster(myCPU-1)
clusterEvalQ(cl,c(library(splines)))
clusterExport(cl,c("var_num","data","order","knot","knotset","Y","bs"))


parll=function(N){
basis=list()
#基底矩陣
for(var in 1:var_num){
  x=data[,var]
  low=min(x)
  up=max(x)
  B=bs(x,deg=order-1, knots=knot(knotset[N,var],low,up), Boundary.knots=c(low,up),intercept=F)
  basis[[var]]=B
}
#基底調整
basismix=do.call(cbind,basis)
u=apply(basismix,2,mean)
for(i in 1:ncol(basismix)){
  basismix[,i]=basismix[,i]+u[i]
}
basismix=cbind(rep(1,nrow(data)),basismix,Y)
colnames(basismix)=c(paste("V",0:(ncol(basismix)-2),sep=""),"y")
basismix=as.data.frame(basismix)
#model
model=lm(y~.-1,data=basismix)
hii.v=lm.influence(model)$hat
rsscv=mean((model$resid/(1-hii.v))^2)
return(rsscv)
}
ptm <- proc.time()
RSSCV=parSapply(cl,1:nrow(knotset),parll)
proc.time() - ptm
stopCluster(cl)
#====
bestset=which.min(RSSCV)
bestknot=knotset[bestset,]
bestknot=as.vector(as.matrix(bestknot))
#最佳節點組合
bestknot
#用最佳節點執行B-spline近似====
basis=list()
for(var in 1:var_num){
  x=data[,var]
  low=min(x)
  up=max(x)
  B=bs(x,deg=order-1, knots=knot(bestknot[var],low,up), Boundary.knots=c(low,up),intercept=F)
  basis[[var]]=B
}
#基底調整
basismix=do.call(cbind,basis)
u=apply(basismix,2,mean)
for(i in 1:ncol(basismix)){
  basismix[,i]=basismix[,i]+u[i]
}
#basismix=apply(basismix,2,Scale)
basismix=cbind(rep(1,nrow(data)),basismix,Y)
colnames(basismix)=c(paste("V",0:(ncol(basismix)-2),sep=""),"y")
basismix=as.data.frame(basismix)
#model
model=lm(y~.-1,data=basismix)
RSS=sum((model$residuals)^2)
AIC=2*(ncol(basismix)-1) + nrow(basismix)*log(RSS/nrow(basismix))
#plot
date=str1(date)

plots=gather(data.frame(True=Y,Fit=model$fitted.values),
             key="keys",value="values",
             True,Fit)
plots$Date=1:length(Y)

P1=ggplot(plots,aes(x=Date,y=values,group=keys,lty=keys,col=keys))+
   geom_line()+
   scale_x_continuous(breaks=c(seq(1,length(Y),500),length(Y)),
                      
                      labels=date[c(seq(1,length(Y),500),length(Y))])+
  
   labs(x="day",y="VFINX",title="VFINX",
       subtitle="1990.01.03~2019.09.27 訓練資料擬合")+
  
   theme(plot.title=element_text(hjust = 0.5,face="bold"))+
   theme(plot.subtitle=element_text(hjust = 0.5,face="bold"))+
   scale_colour_manual(name="",breaks=c("True","Fit"),
                       values = c("Blue","red"),
                       labels=c("實際值","擬合值"))+
   scale_linetype_manual(name="",breaks=c("True","Fit"),
                        values = c(2,1),
                        labels=c("實際值","擬合值"))+
   theme(legend.position="top")

P1+annotate("text",
            label=paste("RSS =",RSS,"\n","AIC =",AIC),
            x = 1000, y =200, 
            size =4, 
            colour = "Black")

#測試====
test=read.csv("預測區間.csv")
date_test=test$X[-1] #日期對應VFINX
test$VFINX[1:nrow(test)-1]=test$VFINX[2:nrow(test)]
test=test[1:nrow(test)-1,]
test=test[,-1]
head(test)
Y_test=test$VFINX

basis_test=list()
for(var in 1:var_num){
  x=test[,var]
  low=min(x)
  up=max(x)
  B=bs(x,deg=order-1, knots=knot(bestknot[var],low,up), Boundary.knots=c(low,up),intercept=F)
  basis_test[[var]]=B
}
#基底調整
basismix_test=do.call(cbind,basis_test)
#basismix_test=apply(basismix_test,2,Scale)
u=apply(basismix_test[,1:ncol(basismix_test)],2,mean)
for(i in 1:ncol(basismix_test)){
  basismix_test[,i]=basismix_test[,i]+u[i]
}
basismix_test=cbind(rep(1,nrow(test)),basismix_test,Y_test)
colnames(basismix_test)=c(paste("V",0:(ncol(basismix_test)-2),sep=""),"y")


basismix_test=as.data.frame(basismix_test)


pred=predict(model,basismix_test)
RSME=sqrt(mean((pred-Y_test)^2))
RSME
sd(Y_test)

#畫圖====
date_test=str2(date_test)
#長資料畫法

plotdata=gather(data.frame(True=Y_test,Pred=pred),key="keys",value="values",
           True,Pred)
plotdata$Date=1:length(Y_test)

P2=ggplot(plotdata,aes(x=Date,y=values,group=keys,lty=keys,col=keys))+
   geom_line()+
   scale_x_continuous(breaks=c(seq(1,length(Y_test),15),length(Y_test)),
                     labels=date_test[c(seq(1,length(Y_test),15),length(Y_test))])+
   labs(x="day",y="VFINX",title="VFINX",
        subtitle="2019.09.30~2020.09.30 測試資料預測"
        )+
   theme(plot.title=element_text(hjust = 0.5,face="bold"))+
   theme(plot.subtitle=element_text(hjust = 0.5,face="bold"))+
   scale_colour_manual(name="",breaks=c("True","Pred"),
                       values = c("black","red"),
                       labels=c("實際值","預測值"))+
   scale_linetype_manual(name="",breaks=c("True","Pred"),
                         values = c(1,1),
                         labels=c("實際值","預測值"))+
   theme(legend.position="top")+
   geom_vline(xintercept=which(date_test=="12/31"),lty=2,col="blue")

P2+annotate("text",
            label=paste("RSME =",RSME,"\n","SD =",sd(Y_test)),
            x = 220, y =150, 
            size =4, 
            colour = "Blue")+
  annotate("text",
           label="2019",
           x = 50, y =350, 
           size =4, 
           colour = "Blue")+
  annotate("text",
           label="2020",
           x = 80, y =350, 
           size =4, 
           colour = "Blue")
  
