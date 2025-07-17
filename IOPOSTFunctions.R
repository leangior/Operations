library(jsonlite)
library(httr)
library(xts)
#----Config File Location
Config.json=("log/Config.json")

#---Query Functions
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
