require(chron)
require(reshape2)
require(shiny)

serverAgg <- function(input,output){
  
  allDomains <- reactive({
    filenames  <- list.files(input$group)
    
    varsaggSpace <- character(0)
    if(input$space.Agg=="byDomain"){
      varsaggSpace = "domainID"
    }else{
      if(input$space.Agg=="bySite"){
        varsaggSpace = c("domainID","siteID")
      }else{
        varsaggSpace = c("domainID","siteID","plotID")
      }
    }
    
    varsaggTime <- character(0)
    if(input$time.Agg=="aggAll"){
      varsaggTime = NULL
    }else{
      if(input$time.Agg=="byBout"){
        varsaggTime = c("yr","month","Bout")
      }else{
        if(input$time.Agg=="byMonthYear"){
          varsaggTime = c("yr","month")
        }else{
          varsaggTime = c("yr")
        }
      }
    }
    if(input$group=="smallMammal"){
      otherFile = "capturedata"
      fieldFiles = "perplotnight"
      uniqueid = c("plotID","eventID","date")
      mixVars = c("samplingEffort")
    }else{
      if(input$group=="carabid"){
        otherFile = "IDandpinning"
        fieldFiles = "fielddata"
        uniqueid = c("plotID","trapID","collectDate")
        mixVars = c("daysOfTrapping","decimalLatitude","decimalLongitude","boutNumber")
      }else{
        if(input$typeFiles=="none"){
          stop("choose 'file type for plant presence' different to 'none'")
        }else{
          otherFile = input$typeFiles
        }
        fieldFiles = "Variables"
        uniqueid = c("plotID","subplotID","date")
        mixVars = NULL#c("percentCover")  
        #if(input$typeFiles=='1m2Data'){
        #  mixVars = "percentCover"
        #}else{
        #  mixVars = NULL#c("percentCover")  
        #}
      }
    }
    currentFiles   <- filenames[grep(otherFile,filenames)]  
    fielddataFiles   <- filenames[grep(fieldFiles,filenames)]  
    
    fileData <- matrix( unlist(strsplit(currentFiles,'[.]')),length(currentFiles),byrow=T)
    siteData <- matrix( unlist(strsplit(fielddataFiles,'[.]')),length(fielddataFiles),byrow=T)
    
    domain  <- fileData[,2]
    site  <- fileData[,3]
    ndomain  <- length(unique(domain))
    nsite  <- length(site)
    
    dataAllDomains <- dataAllField <- list()
    
    for(k in 1:nsite){
      kfile  <- paste(input$group,'/',currentFiles[k],sep='')
      kdat <- read.csv(kfile,header=T) 
      dataAllDomains[[site[k]]] <- kdat
      
      ffile <-  paste(input$group,'/',fielddataFiles[k],sep='')
      fdat <- read.csv(ffile,header=T) 
      dataAllField[[site[k]]] <- fdat[,-which(names(fdat)=="remarks")]
    }
    
    dataAllDomains <- data.frame(do.call(rbind,dataAllDomains),stringsAsFactors = F)
    dataAllField <- do.call(rbind,dataAllField)[,c(uniqueid,mixVars)]
    
    
    dataAllDomains[dataAllDomains==""] = NA
    dateVar <- ifelse(otherFile!='IDandpinning',"date","collectDate")
    formatTime <- "y-m-d"
    if(otherFile%in%c('capturedata','IDandpinning')){
      dataAllDomains$individualID <- as.character(dataAllDomains$individualID)
      dataAllDomains = dataAllDomains[!is.na(dataAllDomains$individualID),]    
    }else{
      dataAllDomains = dataAllDomains[!is.na(dataAllDomains$eventID),]
    }
    dateValuesDom = chron(as.character(dataAllDomains[,dateVar]),format=formatTime)
    dataAllDomains$yr <- as.character(years(dateValuesDom))
    dataAllDomains$month <- as.numeric(months(dateValuesDom))
    if(input$group=="smallMammal"){
      eventIDunique <- as.character(unique(dataAllDomains$eventID))
      ndays <- 0
      indDom <- indField <- logical(0)
      for(id in eventIDunique){
        indDom <- (dataAllDomains$eventID==id)
        indField <- (dataAllField$eventID==id)
        ndays =as.numeric(max(dateValuesDom[indDom])-min(dateValuesDom[indDom]))
        dataAllField$samplingEffort[indField]=dataAllField$samplingEffort[indField]*ndays
      }
    }
    #dateValuesField = chron(as.character(dataAllField[,dateVar]),format=formatTime)
    #dataAllField$yr <- as.character(years(dateValuesField))
    #dataAllField$month <- as.numeric(months(dateValuesField))
    
    dataAllDomains$taxonID <- as.character(dataAllDomains$taxonID)
    dataAllDomains$taxonID[is.na(dataAllDomains$taxonID)] <- "notID"
    
    if(otherFile=="IDandpinning"){
      mm <- as.character(dataAllDomains$month)
      mm[nchar(mm)==1] <- paste0("0",mm[nchar(mm)==1],sep="")
      dataAllDomains$Bout <- with(dataAllDomains,paste(yr,mm,sep="."))  
    }else{
      dataAllDomains$Bout <- unlist(lapply(1:nrow(dataAllDomains),function(rr){
        bout=strsplit(as.character(dataAllDomains$eventID)[rr],split=paste0(dataAllDomains$siteID[rr],"."))[[1]][2]
        bout = strsplit(bout,split=".",fixed=T)[[1]]
        bout= bout[(length(bout)-1):length(bout)]
        if(nchar(bout[2])==1){bout[2]=paste0("0",bout[2])}
        bout = paste0(bout,collapse=".")
        return(bout)
      }))  
    }
    #-----------
    colKeep <- c(names(dataAllDomains),mixVars)
    if(input$group!="plantPresenceCover"){
      dataAllDomains <- merge(x=dataAllDomains, y=dataAllField, by = uniqueid,all.x = T,all.y = F)
      dataAllDomains <- dataAllDomains[!duplicated(dataAllDomains),]  
    }
    list(data=dataAllDomains,varsaggTime=varsaggTime,varsaggSpace=varsaggSpace,mixVars=mixVars)
  })
  
  aggregatedData <- reactive({
    dataAllDomains <- allDomains()$data
    varsaggTime <-  allDomains()$varsaggTime
    varsaggSpace <- allDomains()$varsaggSpace
    mixVars <- allDomains()$mixVars
    
    aggVars <- c(varsaggSpace,varsaggTime)
    
    if(input$group=="carabid"){
      #if(input$time.Agg=="byBout"){
      #  aggVars[aggVars=="Bout"]="boutNumber"
      #  mixVars <- mixVars[mixVars!="boutNumber"]
      #}
      
      aggData <- dcast(dataAllDomains,as.formula(paste(paste(aggVars,collapse="+"),"~ taxonID")),
                       value.var="taxonID",fun.aggregate=length,fill=0)
      
      mixVars <- unique(c(mixVars,"decimalLatitude","decimalLongitude"))
      dataAllDomains$daysOfTrapping <- dataAllDomains$daysOfTrapping*4
      otheraggData <- dcast(melt(dataAllDomains[,c(aggVars,mixVars)],
                                 id.vars=aggVars),as.formula(paste(paste(aggVars,collapse="+"),"~ variable")),
                            fun.aggregate=median,na.rm=T,fill=0)
      otheraggData$daysOfTrapping[otheraggData$daysOfTrapping==0] = NA
    }else{
      #aggData <- dcast(dataAllDomains,as.formula(paste(paste(aggVars,collapse="+"),"~ taxonID")),
      #                 value.var="taxonID",fun.aggregate=length,fill=0)
      #mixVars <- unique(c(mixVars,"decimalLatitude","decimalLongitude"))
      #otheraggData <- dcast(melt(dataAllDomains[,c(aggVars,mixVars)],id.vars=aggVars),
      #                      as.formula(paste(paste(aggVars,collapse="+"),"~ variable")),
      #                      fun.aggregate=median,na.rm=T,fill=0)  
      if((input$group=="plantPresenceCover")&(input$typeFiles=='1m2Data')){
        aggData <- dcast(dataAllDomains,as.formula(paste(paste(aggVars,collapse="+"),"~ taxonID")),
                         value.var="percentCover",fun.aggregate=mean,na.rm=T,fill=0)
        mixVars <- unique(c(mixVars,"decimalLatitude","decimalLongitude"))
        otheraggData <- dcast(melt(dataAllDomains[,c(aggVars,mixVars)],id.vars=aggVars),
                              as.formula(paste(paste(aggVars,collapse="+"),"~ variable")),
                              fun.aggregate=median,na.rm=T,fill=0)  
      }else{
       aggData <- dcast(dataAllDomains,as.formula(paste(paste(aggVars,collapse="+"),"~ taxonID")),
                        value.var="taxonID",fun.aggregate=length,fill=0)
       mixVars <- unique(c(mixVars,"decimalLatitude","decimalLongitude"))
       otheraggData <- dcast(melt(dataAllDomains[,c(aggVars,mixVars)],id.vars=aggVars),
                             as.formula(paste(paste(aggVars,collapse="+"),"~ variable")),
                             fun.aggregate=median,na.rm=T,fill=0)  
      }
    }
    
    aggData <- merge(y=aggData, x=otheraggData, all.x = T,all.y = T)
    #get rid of subsetting
    aggData
  })
  
  summTable <- reactive({
    dataAllDomains <- allDomains()$data
    
    summ<- data.frame(c("Data type:","temporal scale:",
                        "spatial scale:","number of domains:",
                        "number of sites:","number of plots:",
                        "number of species:"),
                      c(input$group,input$time.Agg,input$space.Agg,
                        length(unique(dataAllDomains$domainID)),
                        length(unique(dataAllDomains$siteID)),
                        length(unique(dataAllDomains$plotID)),
                        length(unique(dataAllDomains$taxonID))),row.names = NULL)
    names(summ) = c("","")
    summ
    })
  
  output$aggData <- renderDataTable({aggregatedData()[,1:input$numcols]},
                                    options = list(lengthMenu = c(5, 30, 50), pageLength = 10,
                                                   searching=F))
  
  output$downloadData <- downloadHandler(
    filename <- function(){ paste0(paste(input$group,input$space.Agg,input$time.Agg,sep="_"),".csv")},
    content = function(file){
      write.csv(aggregatedData(), file)
    }
  )
  
  output$summary <- renderDataTable({summTable()},
                                    , options = list(lengthMenu = c(5, 30, 50), pageLength = 10, 
                                                     searching = FALSE))
  
}

