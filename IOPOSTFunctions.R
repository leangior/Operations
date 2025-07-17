library(jsonlite)
library(httr)
library(xts)
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
  for(id in modelId[2:length(modelsId)]){
    v=GetDBRun(siteId,id,VarId[i])
    colnames(v)=c(paste0('model_',id))
    runs=cbind(runs,v)
    i=i+1
    }  
  return(runs)
}

#Analysis Functions

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

#Extended déficit analysis (McMahon, 2007)

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

#Obtención de modelo lineal predictivo a partir de anàlisis de correlación cruzada entre 2 xts (lag and lineal fit)
getLagAndLinealFit<-function(upSerie,downSerie){
  lag=getCrossCor(upSerie,downSerie)
  index(upSerie)=index(upSerie)+which.max(lag)
  tsFrame=cbind(upSerie,downSerie)
  colnames(tsFrame)=c('upSerie','downSerie')
  model=summary(lm(tsFrame$downSerie~tsFrame$upSerie))
  return(model)
}
