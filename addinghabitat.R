##############################################################
##############################################################
## Script title: adding habitat function
## Repository title: Testing connectivity indicators in R
##
## Author: Jenny A. Hodgson
## 
## Copyright: (c) Jenny A. Hodgson 2020
##  
## Email: jenny.hodgson@liverpool.ac.uk; github@condatis.org.uk
## 
## Preferred attribution: "derived from" Hodgson (2020) 
## 'Testing connectivity indicators in R' www.github.com/condatis/indicators2020, 
## DOI: [quote DOI given on github page after release] 
## 
##  License: Open Government License v.3 (see license file shared in the repository)
## 
## Disclaimer:The Information is licensed 'as is' and the Information Provider
## and/or Licensor excludes all representations, warranties, obligations and 
## liabilities in relation to the Information to the maximum extent permitted by law. 
## 
## IP NOTE: Any IP introduced in this script is background IP introduced
## by Jenny Hodgson, within the meaning of IP clauses in 
## CWI-Defra-CEH Aligned Research Agreement,Subcontract PROJECT NUMBER 07111,
## Inteum Number 7123, between UK Centre for Ecology & Hydrology and
## University of Liverpool
## 
## Purpose of script: create new habitat in the landscape
##                    according to 3 spatial arrangements;
##                    aggregated, random and "stepping stones"
##
## Warning: although much of this code has potential to be re-used 
## in projects other than the one for which this script was intended,
## it was not written with broad re-usability in mind. Only
## minimal comments are included for guidance 
##############################################################
##############################################################

addhabitat3ways<- function(options){
  
  # load libraries
  library(raster)
  library(spdep)
  library(rgeos)
  library(rgdal)
  
 
  #extract the options
  regionfname<- options$regionfname # label to use for study region
  hablabel<- options$hablabel # label to use for study habitat
  habtext<- options$habtext #e.g. "semi.natural.grassland",#what the focal habitat is called in the sizes lookup table
  datapath<- options$datapath #file path for writing outputs, e.g. getwd()
  findfhabitatreg<- options$findfhabitatreg # file location for baseline habitat raster
  findhk500<- options$findhk500 #  object name for short-distance Hanski connectivity calculated on baseline habitat raster by function 'hanski'
  findhk2500<- options$findhk2500 # object name for longt-distance Hanski connectivity calculated on baseline habitat raster by function 'hanski'
  fhabitatbasics<- options$fhabitatbasics  # list of landscape stats as created by 'makehabitat' function
   addinglevels<- options$addinglevels # vector of proportions of habitat to add, e.g. c(0.01,0.05)
    resn<-options$resn # resolution of all rasters
  sizesfname<-  options$sizesfname #location of csv file of patch size distributions, e.g. "sizes.csv"#sizesfname should include any necessary path
  
  #load data about starting landscape
  hk500<-get(findhk500)
  hk2500<-get(findhk2500)
  fhabitatreg <- raster(findfhabitatreg)

  #load patch area distribution
  sizes<-read.csv(sizesfname)
  
  cat(paste0("\n", 'Data loaded', "\n", Sys.time(), "\n"))
  
  ###set max sample size 
  
  propnonhab<- (1-fhabitatbasics$proparea)*(fhabitatbasics$regarea/
                                              (ncell(fhabitatreg)*resn*resn/1000/1000))
  
  
  ###create distributions of patch diameter and patch area
  didist<- rep(sqrt(sizes$npix),times= round(sizes[,habtext]*1000))
  ardist<- didist^2*resn*resn/1000/1000
  
  #max sample determined as 10* that needed at the max level of addition desired, given the mean area of each patch and loss of pixels randomly overlapping habitat  
  maxsamp<- round(10* fhabitatbasics$totarea/mean(ardist) * max(addinglevels) /propnonhab )
  
  
  chosen<- sampleRandom(fhabitatreg, size=maxsamp, na.rm=F,  
                        cells=TRUE, rowcol=FALSE, xy=T, sp=FALSE, asRaster=FALSE,)

  chosen<-as.data.frame(chosen)
  
  #remove from sample those on habitat or outside study area (NA)
  nhchosen<-nhchosen<-chosen[!is.na(chosen[,4]) & chosen[,4]==0,]
  
  rm(chosen)
  
  #assign areas to candidates
  nhchosen$di<- base::sample(didist, dim(nhchosen)[1], replace =T)
  
  nhchosen$patchid <- 2:(1+ dim(nhchosen)[1]) #the existing habitat will have id 1
  
  
  #add cells to create bigger patches (around cell centres)
  
  morecells<-nhchosen
  centrecells<-nhchosen
  
  maxring<- (sqrt(max(sizes$npix))-1)/2
  
  for(i in c(1:maxring)){#this will be the no. of cells to count away in each direction
    
    side<-centrecells[(centrecells$di -1)/2 >= i,]
    
    if(dim(side)[1]>0){
      sideall<- side
      sideall$x<- sideall$x - (i*resn)
      sideall$y<- sideall$y - (i*resn)
      
      for(j in (-i+1):i){
        
        #left side
        side2<- side
        side2$x<- side2$x - (i*resn)
        side2$y<- side2$y + (j*resn)
        
        sideall<- rbind(sideall,side2)
        
        #right side
        side2<- side
        side2$y<- side2$y - (i*resn)
        side2$x<- side2$x + (j*resn)
        
        sideall<- rbind(sideall,side2)
        
        #right side
        side2<- side
        side2$x<- side2$x + (i*resn)
        side2$y<- side2$y + (j*resn)
        
        sideall<- rbind(sideall,side2)
        
        #top side
        side2<- side
        side2$y<- side2$y + (i*resn)
        side2$x<- side2$x + (j*resn)
        
        if (j<i){#this to prevent last, top-right cell being duplicated
          sideall<- rbind(sideall,side2)
        }#end if
      }#end j loop: filling 4 sides of square
      morecells<-rbind(morecells,sideall)
    }}#end if and i loop
  
  cat(paste0("\n", 'Cells added', "\n", Sys.time(), "\n"))
  
  #candidates not overlapping each other, favour lowest cell number for consistency
  
  morecells2<-aggregate(morecells[,c(1,4)],by=morecells[,c("x","y")],FUN=min)

  morecells2<-merge(morecells2,morecells,all.x=T,all.y=F)
  
  #candidates not overlapping existing habitat
  
  nonhab <- subs(x=fhabitatreg,y=data.frame(id=1, v=NA), subsWithNA=F)
  
  morecells<-rasterize(x=morecells2[,c("x","y")],y=nonhab, field=as.matrix(morecells2$patchid),update=TRUE,updateValue='!NA')

  #turn back to points for selection
  
  morecells2<-rasterToPoints(morecells,fun=function(x){x>0})

  
  #aggregate by id
  
  nhchosen<-nhchosen[nhchosen$patchid %in% morecells2[,"layer"],] #so,only centres of those patches that have some area left
  
  morecells2<-as.data.frame(morecells2)
  names(morecells2)[3]<-"patchid"
  
  nhchosen$arealeft<- as.vector(table(morecells2$patchid)[as.character(nhchosen$patchid)])  *resn*resn/1000/1000
  
  #get connectivity values - here where hanski rasters get used
  
  nhchosen$hk2500<- extract(x=hk2500$conr,y=nhchosen[,c("x","y")], method = "bilinear")

  nhchosen$hk500<- extract(x=hk500$conr,y=nhchosen[,c("x","y")], method = "bilinear")

  
  nhchosen$priority1<- (nhchosen$hk500 + nhchosen$hk2500)/2
  
  nhchosen$priority2<- (1+ nhchosen$hk2500 - nhchosen$hk500)^2
  
  nhchosen$priority2<- nhchosen$priority2 - min(nhchosen$priority2,na.rm=T)#somewhat arbitrary but there is no intuitive baseline to choose for this metric
  
  nhchosen$priority1[is.na(nhchosen$priority1 )]<-0 #there should not be NAs once edge effects are fixed
  nhchosen$priority2[is.na(nhchosen$priority2 )]<-0 
  
  #order candidates in 3 ways and place thresholds
  addingtargets<- fhabitatbasics$totarea* addinglevels
  maxartarget<- fhabitatbasics$totarea* max(addinglevels,na.rm=T)
  
  quant1<- maxartarget/mean(nhchosen$arealeft,na.rm=T)
  
  ##safer no.of patches needed based on sd
  
  asd<- sd(ardist)
  
  sol1<- 2.95*asd + sqrt( (2.95*asd)^2 + 4*mean(nhchosen$arealeft,na.rm=T)^2*quant1)
  sol1<- (sol1/2/mean(nhchosen$arealeft,na.rm=T))^2
  
  
  pick1<- sample.int(dim(nhchosen)[1], size = round(sol1), replace = FALSE, prob = nhchosen$priority1)
  
  pick2<- sample.int(dim(nhchosen)[1], size = round(sol1), replace = FALSE, prob = nhchosen$priority2)
  
  ##add buffer to area target to prevent consistently undershooting
  addingtargets<- addingtargets+ mean(nhchosen$arealeft,na.rm=T)/2
  
  #create data frame to return, to check basic stats
  verify<-expand.grid(add=addinglevels,method=c("random","hihanski","ridges"))
  
  ###merge different thresholded versions with original habitat
  
  #first the random picking
  nhchosen<-nhchosen[order(nhchosen$patchid),]
  
  nhchosen$cumarea<-cumsum(nhchosen$arealeft)
  
  nhchosen$cumarea<-cut(nhchosen$cumarea,breaks=c(0,addingtargets,max(nhchosen$cumarea)))

  
  for(i in 1:length(addingtargets)){
    
    pickedids<- nhchosen$patchid[as.numeric(nhchosen$cumarea)<=i]
    
    morecells<-rasterize(x=morecells2[morecells2$patchid %in% pickedids ,c("x","y")],y=fhabitatreg, field=1,update=TRUE,updateValue='all')

    writeRaster(morecells, filename=paste0(datapath, "/", paste(regionfname,hablabel,"random.plus",addinglevels[i]*100,"tif",sep=".")),overwrite=TRUE)
 
    verify$nadded[verify$add==addinglevels[i] & verify$method=="random"] <- length(pickedids)
    
    verify$areaadded[verify$add==addinglevels[i] & verify$method=="random"] <- sum(nhchosen$arealeft[as.numeric(nhchosen$cumarea)<=i])
  }
  
  cat(paste0("\n", 'Random created', "\n", Sys.time(), "\n"))
  
  #next picking highest connectivity (priority1)
  
  nhchosen1<-nhchosen[pick1,]
  
  nhchosen1$cumarea<-cumsum(nhchosen1$arealeft)
  
  nhchosen1$cumarea<-cut(nhchosen1$cumarea,breaks=c(0,addingtargets,max(nhchosen1$cumarea)))
 
  
  for(i in 1:length(addingtargets)){
    
    pickedids<- nhchosen1$patchid[as.numeric(nhchosen1$cumarea)<=i]
    
    morecells<-rasterize(x=morecells2[morecells2$patchid %in% pickedids ,c("x","y")],y=fhabitatreg, field=1,update=TRUE,updateValue='all')

    writeRaster(morecells, filename=paste0(datapath, "/", paste(regionfname,hablabel,"hihanski.plus",addinglevels[i]*100,"tif",sep=".")),overwrite=TRUE)

    verify$nadded[verify$add==addinglevels[i] & verify$method=="hihanski"] <- length(pickedids)
    
    verify$areaadded[verify$add==addinglevels[i] & verify$method=="hihanski"] <- sum(nhchosen1$arealeft[as.numeric(nhchosen1$cumarea)<=i])
  }
  
  cat(paste0("\n", 'High connectivity created', "\n", Sys.time(), "\n"))  
  
  #next picking ridges (priority2)
  
  nhchosen1<-nhchosen[pick2,]
  
  nhchosen1$cumarea<-cumsum(nhchosen1$arealeft)
  
  nhchosen1$cumarea<-cut(nhchosen1$cumarea,breaks=c(0,addingtargets,max(nhchosen1$cumarea)))
 
  for(i in 1:length(addingtargets)){
    
    pickedids<- nhchosen1$patchid[as.numeric(nhchosen1$cumarea)<=i]
    
    morecells<-rasterize(x=morecells2[morecells2$patchid %in% pickedids ,c("x","y")],y=fhabitatreg, field=1,update=TRUE,updateValue='all')

    writeRaster(morecells, filename=paste0(datapath, "/", paste(regionfname,hablabel,"ridges.plus",addinglevels[i]*100,"tif",sep=".")),overwrite=TRUE)
  
    verify$nadded[verify$add==addinglevels[i] & verify$method=="ridges"] <- length(pickedids)
    
    verify$areaadded[verify$add==addinglevels[i] & verify$method=="ridges"] <- sum(nhchosen1$arealeft[as.numeric(nhchosen1$cumarea)<=i])
  }
  
  cat(paste0("\n", 'Stepping stones created', "\n", Sys.time(), "\n"))
  
  return(verify) ##return this to check what happened
}#end function
