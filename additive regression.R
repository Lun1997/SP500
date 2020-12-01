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
colnames(data)[na_col]
date[na_row]
library(ggplot2)
plot=data.frame(day=1:nrow(data),rom=data[,na_col])
ggplot(data=plot)+geom_line(aes(x=day,y=rom))+
  geom_vline(xintercept=na_row,col=2,lty=2)+
  geom_point(x=na_row-1,y=data[na_row-1,na_col],col="Blue",cex=1)+
  geom_point(x=na_row+1,y=data[na_row+1,na_col],col="Blue",cex=1)+
  scale_x_continuous(breaks=na_row,labels=date[na_row])

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


bestset=which.min(RSSCV)
bestknot=knotset[bestset,]
bestknot=as.vector(as.matrix(bestknot))
#最佳節點組合
bestknot

#用最佳節點執行B-spline近似
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
basismix=cbind(rep(1,nrow(data)),basismix,Y)
colnames(basismix)=c(paste("V",0:(ncol(basismix)-2),sep=""),"y")
basismix=as.data.frame(basismix)
#model
model=lm(y~.-1,data=basismix)
#LOOCV====
library(DAAG)
LOOCV=cv.lm(basismix,model,m=nrow(basismix),plotit = F)
head(LOOCV)


#視覺化====
library(ggplot2)
date=as.character(date)
str=function(string){
  mon=substring(string,6,7)
  day=substring(string,9,10)
  return(paste(mon,"/",day,sep=""))
}
date=str(date)
plotdata=data.frame(Date=1:length(Y),true=Y,Pred=LOOCV$cvpred,fit=LOOCV$Predicted)
P=ggplot(data=plotdata,aes(x=Date))+
   geom_line(aes(y=fit,col="擬合值"))+
   geom_line(aes(y=true,col="實際值"),lty=2)+
   geom_line(aes(y=Pred,col="LOOCV預測值"))+
   scale_x_continuous(breaks=c(seq(1,length(Y),9),length(Y)),
                      labels=date[c(seq(1,length(Y),9),length(Y))])+
   labs(x="day",y="S&P500",title="S&P500",
        subtitle="2020.01~2020.10",
        caption="leave one out prediction")+
   theme(plot.title=element_text(hjust = 0.5,face="bold"))+
   theme(plot.subtitle=element_text(hjust = 0.5,face="bold"))+
   scale_colour_manual(name="",breaks=c("實際值","擬合值","LOOCV預測值"),
                       values = c("black","blue","red"))
   
   
   
  

#統計====
a1=paste("模型殘差     :",sqrt(mean((LOOCV$y-LOOCV$Predicted)^2)))
a2=paste("預測誤差     :",sqrt(mean((LOOCV$y-LOOCV$cvpred)^2)))
a3=paste("S&P500標準差 :",sd(Y))
paste(a1,"\n",a2,"\n",a3)
P+geom_text(x=160,y=5000,
            aes(label=paste(a1,"\n",a2,"\n",a3)),
            size=4,family="Times New Roman")

  

