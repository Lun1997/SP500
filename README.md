# Model : additive regression(nonlinear) 

* 目標:VFINX調整收盤價
* 訓練資料:1990.01.02~2019.09.27
* feacture:VIX VBMFX SP500

資料型態:7494個營業日 x 3種金融指標

| Date     | VIX         | S&P500      | VBMFX       | VFINX       |
|----------|-------------|-------------|-------------|-------------|
| 1990/1/2 | 17.23999977 | 386.1600037 | 2.03178525  | 17.95616531 |
| 1990/1/3 | 18.19000053 | 385.1700134 | 2.027475595 | 17.90896797 |
| 1990/1/4 | 19.21999931 | 382.019989  | 2.029630184 | 17.76213264 |
| 1990/1/5 | 20.11000061 | 378.2999878 | 2.027475595 | 17.58907127 |
| 1990/1/8 | 20.26000023 | 380.0400085 | 2.025319576 | 17.67298698 |

#### 今日指數預測明日目標
```R
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
```

```R
#資料讀取&整理====
#data====
data=read.csv("建模資料.csv")
date=data$X[-1] #日期對應VFINX，用前一天的特徵預測當日 VFINX
data$VFINX[1:nrow(data)-1]=data$VFINX[2:nrow(data)]
data=data[1:nrow(data)-1,]
data=data[,-1]
head(data)
#NA:沒有遺失值====
sum(is.na(data))
```
```
[1] 0 
```

#### 配置可加模型，基於B-spline基底近似估計

模型調參:用leave one out交叉驗證選出最佳參數組合
```R
knot=function(x,lowbond,upperbond){
  if(x==0){return(NULL)}
  else{return(seq(lowbond,upperbond,by=(upperbond-lowbond)/(x+1))[c(-1,-(x+2))])}
}
```
```R
#cubic spline
order=4
var_num=ncol(data)-1
Y=data$VFINX
```
#### 超參數:B-spline節點組合
```R
#試行節點組合，以LOOCV取最佳組合
kn=15:30
knotset=expand.grid(kn,kn,kn)
```
#### 交叉驗證(平行運算)
```R
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
```

#### 用最佳參數配置模型
```R
bestset=which.min(RSSCV)
bestknot=knotset[bestset,]
bestknot=as.vector(as.matrix(bestknot))
```
```R
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
u=apply(basismix,2,mean)
for(i in 1:ncol(basismix)){
  basismix[,i]=basismix[,i]+u[i]
}
basismix=cbind(rep(1,nrow(data)),basismix,Y)
colnames(basismix)=c(paste("V",0:(ncol(basismix)-2),sep=""),"y")
basismix=as.data.frame(basismix)
```
```R
#model
model=lm(y~.-1,data=basismix)
RSS=sum((model$residuals)^2)
AIC=2*(ncol(basismix)-1) + nrow(basismix)*log(RSS/nrow(basismix))
```
#### 模型擬合情形
```R
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
```
![image](https://github.com/Lun1997/SP500/blob/main/VFINX%20%E6%93%AC%E5%90%88.jpg)

#### 預測
```R
#測試
test=read.csv("預測區間.csv")
date_test=as.character(test$X[-1]) #日期對應VFINX
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
u=apply(basismix_test[,1:ncol(basismix_test)],2,mean)
for(i in 1:ncol(basismix_test)){
  basismix_test[,i]=basismix_test[,i]+u[i]
}
basismix_test=cbind(rep(1,nrow(test)),basismix_test,Y_test)
colnames(basismix_test)=c(paste("V",0:(ncol(basismix_test)-2),sep=""),"y")
basismix_test=as.data.frame(basismix_test)
```
#### 預測效果

```R
pred=predict(model,basismix_test)
RMSE=sqrt(mean((pred-Y_test)^2))
RMSE 
sd(Y_test)

#plot
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
            label=paste("RMSE =",RMSE,"\n","SD =",sd(Y_test)),
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
```
![image](https://github.com/Lun1997/SP500/blob/main/VFINX%20%E9%A0%90%E6%B8%AC.jpg)

#### 結論
* 訓練過程overfitting，預測效果有限
* 此模型在實務應用上參數量會快速膨脹導致運算效果降低



