library(jsonlite)
library(httr)
library(xts)
library(fitdistrplus)
library(lubridate)
library(forecast)
#----Config File Location

Config.json=("log/Config.json")

#---Query Functions

#Obtiene xts de datos puntuales (servers a5/a6)
GetDBSerie<-function(Id,Start,End,Agg='none',Fill='yes',configFile=Config.json){
  data=list()  
  config=fromJSON(configFile)
  URI=paste0(config$api$url,"/obs/puntual/series/",Id,"/observaciones?timestart=",as.Date(Start),"&timeend=",as.Date(End)+1)
  request=GET(URI,add_headers("Authorization"=config$api$token),accept_json())
  x=fromJSON(content(request,"text"))
  data=xts(x=x$valor,order.by=as.Date(x$timestart)) #transforma consulta a objeto ts
  if(Agg!='none'){
    data=AggTimeSeries(data,Agg)
  }
  if(Fill=='yes'){
    return(na.approx(data))
  }
  else{
    return(data)
  }
}

#Obtiene xts de datos areales (servers a5/a6)
getDBArealSerie<-function(Id,Start,End,Agg='none',Fill='yes',configFile=Config.json){
  data=list()  
  config=fromJSON(configFile)
  URI=paste0(config$api$url,"/obs/areal/series/",Id,"?timestart=",as.Date(Start),"&timeend=",as.Date(End)+1)
  request=GET(URI,add_headers("Authorization"=config$api$token),accept_json())
  x=fromJSON(content(request,"text"))
  data=xts(x=x$observaciones$valor,order.by=as.Date(x$observaciones$timestart)) #transforma consulta a objeto ts
  if(Agg!='none'){
    data=AggTimeSeries(data,Agg)
  }
  if(Fill=='yes'){
    return(na.approx(data))
  }
  else{
    return(data)
  }
}

#Obtiene xts de previsiones puntuales (servers a5/a6)
GetDBRun<-function(StId,CalId,RunId,Start=Sys.Date()+1,End=Sys.Date()+15,Agg='none',Fill='yes',configFile=Config.json,VarId=4,Qualy='main'){
  data=list()  
  config=fromJSON(configFile)
  URI=paste0(config$api$url,"/sim/calibrados/",CalId,"/corridas/last?estacion_id=",StId,"&var_id=",VarId,"&qualifier=",Qualy,"&includeProno=true",sep="")
  request=GET(URI,add_headers("Authorization"=config$api$token),accept_json())
  x=fromJSON(content(request,"text"))
  data=xts(x=as.numeric(x$series$pronosticos[[1]][,3]),order.by=as.Date(x$series$pronosticos[[1]][,1])) #transforma consulta a objeto ts
  if(Agg!='none'){
    data=AggTimeSeries(data,Agg)
  }
  if(Fill=='yes'){
    return(na.approx(data))
  }
  else{
    return(data)
  }
}

#Obtiene xts de conjunto de previsiones (servers a5/a6)
getDBEnsembleRun<-function(siteId=2,modelId=c(65,69,32,319,320,321),VarId=c(4,4,4,4,4,4)){
  i=1
  runs=xts()
  for(id in modelId[2:length(modelId)]){
    v=GetDBRun(siteId,id,VarId[i])
    colnames(v)=c(paste0('model_',id))
    runs=cbind(runs,v)
    i=i+1
    }  
  return(runs)
}

#Postprocesa xts de conjunto de previsiones a fin de obtener tendencia central y espacio probable
get_forecasts_space<-function(xts_forecasts,probs_thresholds=c(0.01,0.99)){
  forecast_space=c()
  for(i in seq(1,dim(xts_forecasts)[1])){
    q=as.numeric(quantile(xts_forecasts[i,],probs=probs_thresholds))
    avg=mean(c(xts_forecasts[i,]))
    median=median(c(xts_forecasts[i,]))
    central=mean(c(median,avg))
    values=c(q,central)
    forecast_space=rbind(forecast_space,values)   
  }
  forecast_space=xts(order.by=index(xts_forecasts),forecast_space)
  names=c()
  for(i in seq(1,length(probs_thresholds))){
    names[i]=c(paste0("p_",probs_thresholds[i]))
  }
  names[dim(forecast_space)[2]]=paste("main")
  colnames(forecast_space)=names
  return(forecast_space)
}

#---Analysis Functions

#Rutina para creación de función distribución sobre la base de n reservorios lineales con tiempo promedio de residencia igual a k. Integra por trapecio, a subpasos de cálculo de longitud dt.
CreateGamma<-function(k=1,n=1,m=10*n*k,dt=0.01){
  t=0
  Q=c()
  while(t<=m){
    u=1/(k*gamma(n))*(t/k)^(n-1)*exp(-t/k)*1/2
    for(j in seq(t+dt,t+1-dt,dt)){
      u=u+1/(k*gamma(n))*(j/k)^(n-1)*exp(-j/k)
    }
    u=u+1/(k*gamma(n))*((t+1)/k)^(n-1)*exp(-(t+1)/k)*1/2
    t=t+1
    Q[t]=u
  }
  return(Q/sum(Q))
}


#Obtiene la curva de régimen intra-anual para un xts de datos puntuales (servers a5/a6)
getMeanRegime=function(seriesId,refStart='1991-01-01',refEnd='2020-12-31'){
  regime=list()
  regime$Id=seriesId
  regime$refStart=as.Date(refStart)
  regime$refEnd=as.Date(refEnd)
  serie=GetDBSerie(seriesId,regime$refStart,regime$refEnd,Agg='Monthly')
  regime$reg=array(aggregate(serie,cycle(as.yearmon(index(serie))),mean))
  return(regime)
}

#Extended deficit analysis (McMahon, 2007)

  #obtención de xts EDA para xts puntual (servers a5/a6)
extendedDeficitAnalysis<-function(serieId,thresHold,start='2010-01-01',end=Sys.Date(),f=10^6,z0=0){
  serie=GetDBSerie(serieId,start,end,Agg='Daily')
  z=c()
  z[1]=z0
  for(i in seq(2,length(serie))){
    z[i]=min(c((as.numeric(serie[i])-thresHold)*86400+z[i-1],z0))
  }
  return(xts(order.by=index(serie),z/f))
}

#obtención de xts EDA a partir de xts custom
extendedDeficitAnalysisbySerie<-function(serie,thresHold,f=10^6,z0=0){
  z=c()
  z[1]=z0
  for(i in seq(2,length(serie))){
    z[i]=min(c((as.numeric(serie[i])-thresHold)*86400+z[i-1],z0))
  }
  return(xts(order.by=index(serie),z/f))
}

#Obtención de Curva Almacenamiento-Descarga sobre la base de análisis de xts EDA
getVolumenReleasebyDefecitiAnalysis<-function(serie,thresHolds,f=10^6,z0=0,k=1.68,p=.99){
  r=c()
  for(threshold in thresHolds){
    defSerie=extendedDeficitAnalysisbySerie(serie,threshold,f,z0)
    sam=peakdet(as.numeric(-1*defSerie),delta = k*sd(defSerie))
    sam=-1*defSerie[sam$maxtab$pos]
    q=quantile(sam,probs=c(p))
    r=rbind(r,c(threshold,q))
  }
  colnames(r)=c('release[m³/s]','volume_deficit[hm³]')
  return(r)
}

#Anàlisis de correlación cruzada entre 2 xts
getCrossCor<-function(upSerie,dSerie,maxLag=5,ini=1){
  r=c()
  for(i in seq(ini,maxLag)){
    j=i-ini+1
    uL=upSerie[i:c(length(upSerie))-i]
    index(uL)=index(uL)+i
    r[j]=cor(cbind(uL,dSerie)[c(i+1):dim(dSerie)[1]],use='complete.obs')[2,1]
  }
  return(r)
}


#--- Modelos de propagaciòn lineal

#A. Estáticos
#Obtención de modelo lineal predictivo a partir de anàlisis de correlación cruzada entre 2 xts (lag and lineal fit)
getLagAndLinealFit<-function(upSerie,downSerie){
  lag=getCrossCor(upSerie,downSerie)
  index(upSerie)=index(upSerie)+which.max(lag)
  tsFrame=cbind(upSerie,downSerie)
  colnames(tsFrame)=c('upSerie','downSerie')
  model=summary(lm(tsFrame$downSerie~tsFrame$upSerie))
  return(model)
}

#B. Dinámicos
#Rutina para la generación una matriz de convolución. Requiere vector x de datos de hidrograma de entrada y el vector Disfun función de distribución
GetPulsesMatrix<-function(x=data,DisFun=u){
  m=length(DisFun) #ordenadas de función distribución
  n=length(x) #ordenadas de hidrograma de entrada
  Q=matrix(nrow=sum(n,m),ncol=m)
  k=1
  for(col in 1:m){
    for(row in 1:sum(n,m)){
      if(col>row){
        Q[row,col]=0
      }
      else
      {
        if(row>=n+k){
          Q[row,col]=0
        }
        else{
          Q[row,col]=x[row-k+1]
        }
        
      }
    }
    k=k+1
  }
  return(Q)
}

#Rutina de trànsito lineal utilizando kernel obtenido mediante gamma(k,n)
LinearTransit<-function(ts,k,n,U){
  if(missing(U)){
    U=CreateGamma(k,n) #crea función de distrubución (tránsito por N reservorios lineales con tiempo promedio de residencia igual a K)
  }
  Response=list() #declara objeto respuesta (lista, general)está 
  PulsesMatrix=GetPulsesMatrix(x=as.matrix(ts)[1:length(ts)],DisFun=U) #crea matriz de pulsos
  Response$val=PulsesMatrix%*%U #realiza convolución de pulsos
  Response=xts(x=Response$val[1:length(ts)],order.by=index(ts)) #transforma response a objeto ts
  return(na.fill(Response,0))
  # na.fill(return(window(Response,start=as.Date(strsplit(x = Interval,split="/")[[1]][1]),end=as.Date(strsplit(x = Interval,split="/")[[1]][2]))),0)
}


#Modelo Lag And Route para tránsito de señales (Desplazamiento rígido + efecto de almacenamiento)
lagAndRoute<-function(upSerie,lag=0,k=0.01,n=1,warmUp=30){
  route=upSerie
  index(route)=index(route)+lag
  route=LinearTransit(route,k,n)
  return(route[warmUp:dim(route)[1]])
}

#Modelo Lag And Route para tránsito de señales con ajuste de sesgo por regresión lineal
lagAndRouteBiasAdj<-function(upSerie,downSerie,lag=0,k=0.01,n=1,warmUp=30){
  route=lagAndRoute(upSerie,lag,k,n,warmUp)
  series=cbind(downSerie,route)
  series=series[warmUp:dim(series)[1]]
  colnames(series)=c('obs','sim')
  model=summary(lm(series$obs~series$sim))
  message(paste0("Utilizando ajuste de sesgo por regresión lineal R²=",model$r.squared))
  series$sim=model$coef[1]+model$coef[2]*series$sim
  return(series)
}

