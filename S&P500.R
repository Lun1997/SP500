library(splines)
#data====
data=read.csv("測試資料.csv")
head(data)
date=data$Date[-1]
data$spx[1:nrow(data)-1]=data$spx[2:nrow(data)]
data=data[1:nrow(data)-1,]
data=data[,-1]
head(data)
#NA====
sum(is.na(data))
na_col=which(is.na(data))%/%nrow(data) +1
na_row=which(is.na(data))%%nrow(data)
if(na_row==0){na_row=nrow(data)}
plot(data[,na_col],type="l")
points(c(na_row-1,na_row+1),
       data[c(na_row-1,na_row+1),na_col],
       pch=16,col="blue",cex=0.6)
abline(v=na_row,lty=2,col=2)
#遺漏值位在相對平穩的區段，用前後值取平均補值
data[na_row,na_col]=mean(data[na_row-1,na_col],
                         data[na_row+1,na_col])

#B-spline====
knot=function(x,lowbond,upperbond){
  if(x==0){return(NULL)}
  else{return(seq(lowbond,upperbond,by=(upperbond-lowbond)/(x+1))[c(-1,-(x+2))])}
}

#Quadratic spline
order=3 
var_num=ncol(data)-1
Y=data$spx

#嘗試節點組合
kn=0:10
knotset=expand.grid(kn,kn,kn,kn,kn)
RSSCV=c()
for(N in 1:nrow(knotset)){
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
  basismix=cbind(rep(1,nrow(data)),basismix,Y)
  colnames(basismix)=c(paste("V",0:(ncol(basismix)-2),sep=""),"y")
  basismix=as.data.frame(basismix)
  #model
  model=lm(y~.-1,data=basismix)
  hii.v=lm.influence(model)$hat
  rsscv=mean((model$resid/(1-hii.v))^2)
  RSSCV[N]=rsscv
}

min(RSSCV)
bestset=which.min(RSSCV)
bestknot=knotset[bestset,]
bestknot=as.vector(as.matrix(bestknot))
#最佳節點組合
bestknot

#用最佳節點執行B-spline近似

for(var in 1:var_num){
  x=data[,var]
  low=min(x)
  up=max(x)
  B=bs(x,deg=order-1, knots=knot(bestknot[var],low,up), Boundary.knots=c(low,up),intercept=F)
  basis[[var]]=B
}
#基底調整
basismix=do.call(cbind,basis)
basismix=cbind(rep(1,nrow(data)),basismix,Y)
colnames(basismix)=c(paste("V",0:(ncol(basismix)-2),sep=""),"y")
basismix=as.data.frame(basismix)
#model
model=lm(y~.-1,data=basismix)
#LOOCV====
library(DAAG)
LOOCV=cv.lm(basismix,model,m=nrow(basismix),plotit = F)
head(LOOCV)
library(ggplot2)
date=as.character(date)
str=function(string){
  mon=substring(string,6,7)
  day=substring(string,9,10)
  return(paste(mon,"/",day,sep=""))
}
date=str(date)

#視覺化====
plotdata=data.frame(Date=1:length(Y),true=Y,Pred=LOOCV$cvpred)
ggplot(data=plotdata,aes(x=Date))+
  geom_line(aes(y=true,col="實際值"))+
  geom_line(aes(y=Pred,col="預測值"))+
  scale_x_continuous(breaks=c(seq(1,length(Y),9),length(Y)),
                     labels=date[c(seq(1,length(Y),9),length(Y))])+
  labs(x="day",y="S&P500",title="S&P500",
       subtitle="2020.01~2020.10",
       caption="leave one out prediction")+
  theme(plot.title=element_text(hjust = 0.5,face="bold"))+
  theme(plot.subtitle=element_text(hjust = 0.5,face="bold"))+
  scale_colour_manual(name = "",
                      values = c("Blue","red"))

#模型診斷====
#模型殘差
sqrt(mean((LOOCV$y-LOOCV$Predicted)^2))
#預測誤差
sqrt(mean((LOOCV$y-LOOCV$cvpred)^2))