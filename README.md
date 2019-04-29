# Film-cyclic-voltammetry-analysis
A code in Rscript that allows cyclic voltammograms showing film-based redox signals to be analysed

#This code reads in xy files. These files should be systematically named
#The scan rates used should be part of the file names. List the scan rates in the vector below

Scanrates = c(1000,750,500,400,300,200,100,50,25)

#Define a vector of colours to be used in the graphs
Colours=c("Black","blue","red","purple","brown","orange","grey53","green","pink","yellow")

#If desired, the following function can be used to filter out systematic noise (i.e. mains frequency)
Smoothed=function(data){
  x=as.vector(data[,1])
  y=as.vector(data[,2])
  datatoSmooth = smooth.spline(x, y, spar=0.5)
  predictdata=predict(datatoSmooth, newdata=data.frame(x=x.new))
  
  return(predictdata)
}

#These functions find the half height coordinates of any extracted peaks
findhalfheightback=function(data){
  peak=min(data[,2], na.rm=TRUE)
  
  halfheight=peak/2
  halfheights=c()
  for(i in 1:length(data[,1])){
    
    if (round(data[i,2], digits = 1) == round(halfheight, digits = 1)){ 
      halfheights=c(halfheights,data[i,1])
    }
    
  }
  
  widthathalfheight=abs(max(halfheights)-min(halfheights))
  
  return(widthathalfheight)}


findhalfheightforward=function(data){
  peak=max(data[,2], na.rm=TRUE)
  
  halfheight=peak/2
  halfheights=c()
  for(i in 1:length(data[,1])){
    
    if (round(data[i,2], digits = 1) == round(halfheight, digits = 1)){ 
      halfheights=c(halfheights,data[i,1])
    }
    
  }
  
  widthathalfheight=abs(max(halfheights)-min(halfheights))
  
  return(widthathalfheight)}

#This function finds the position of the extracted anodic peak
Peakanodic=function(data){
  peakanoidic=max(data[,2], na.rm=TRUE)
  for(i in 1:length(data[,1])){
    if (data[i,2] == peakanoidic){
      centre=data[i,1]}
    
  }
  return(centre)
}

#This function finds the position of the extracted cathodic peak
Peakcathoic=function(data){
  peakcathoic=min(data[,2], na.rm=TRUE)
  for(i in 1:length(data[,1])){
    if (data[i,2] == peakcathoic){
      centre=data[i,1]}
    
  }
  return(centre)
}

#These functions split the xy dataset into half so that the cathodic and anodic sweeps can be analysed separately
Chopf = function(xy){
  a = length(xy[,1])
  b = length(xy[,1])/2
  
  forward = xy[-c(b:a),]
  return(forward)
}

Chopb = function(xy){
  a = length(xy[,1])
  b = length(xy[,1])/2
  
  back = xy[-c(1:b),]
  return(back)
}

#These functions are used to remove parts of the cyclic voltammogram that aren't of use for peak extraction
Clip = function(xy,clipfrommax,clipfrommin){
  x = abs(xy[2,1]-xy[1,1])
  y=length(xy[,1])
  z=xy[y,1]
  q=xy[1,1]
  datapointstoremove1=round(clipfrommax/x)
  Rowstoremove1 = c(y + 1 - 1:datapointstoremove1)
  xy = xy[-Rowstoremove1,]
  
  datapointstoremove2=round(clipfrommin/x)
  Rowstoremove2 = c(1:datapointstoremove2)
  xy = xy[-Rowstoremove2,]
  
  return(xy)
}

Clip2 = function(xy,Fitbaselineusing1,Fitbaselineusing2){
  
  y=length(xy[,1])
  dist=c()
  
  for(i in 1:y-1){
    x = abs(xy[i+1,1]-xy[i,1])
    dist=c(dist,x)}
  
  x=mean(dist)
  z=xy[y,1]
  q=xy[1,1]
  g = xy[y,1]-Fitbaselineusing1
  datapointstokeep1=round(g/x)
  Rowstokeep1 = c(y + 1 - 1:datapointstokeep1)
  
  f = Fitbaselineusing2-xy[1,1]
  datapointstokeep2=round(f/x)
  Rowstokeep2 = c(1:datapointstokeep2)
  xy = xy[c(Rowstokeep2,Rowstokeep1),]
  
  return(xy)
}

#This function reads in the xy files and processes it contain only 2 columns in the desired units (it would need modification for different xy data file layouts)
Process = function(targetfile){
  
  t <- paste(targetfile, sep = "")
  
  x <- read.table(t, skip=1, header=TRUE 
  )
  
  rownames(x) <- 1:nrow(x)
  x=x[,-1]
  x=x[,-3]
  x[,2]=x[,2]*1000000
  
  colnames(x) <- c("Voltage (V)", " Current (µA)")
  
  return (x)
}

#This function is used to integrate isolated signals and equate it to surface coverage of redox active units, where n is number of electrons tranferred in the redox couple
Findcoverage=function(Isolated_peak,electrodearea,scanrate,n){
  
  Isolated_peak[ "New1" ] <- NA
  Isolated_peak[ "New2" ] <- NA
  Isolated_peak[ "New3" ] <- NA
  
  colnames(Isolated_peak) <- c("Voltage (V)", " Current (µA)", "E2-E1", "Current2+Current1", "Trapesium")
  
  for (i in 1:length(Isolated_peak[,1])){
    Isolated_peak[i+1,3]=Isolated_peak[i+1,1] - Isolated_peak[i,1]
    Isolated_peak[i+1,4]=Isolated_peak[i+1,2] + Isolated_peak[i,2]
  }
  
  Isolated_peak=Isolated_peak[-c(1),]
  
  for (i in 1:length(Isolated_peak[,1])){
    Isolated_peak[i,5]=0.5*Isolated_peak[i,3]*Isolated_peak[i,4]
  }
  
  Integral_of_peak=sum(Isolated_peak$Trapesium,na.rm=TRUE)
  
    x=Integral_of_peak/1000000
    y=scanrate/1000
    z=x/y
    g=z*1000000000000/96485.3329
    j=g/n
    A=j/electrodearea
    
    return(A)
}
  
#The code below generates a list of the names of the systematic xy data files, i.e. a file could be called "1000_3" if it is scan 3 at 1000 mV s^-1
Namelist=c()

for(i in 1:length(Scanrates)){
Scanrate = Scanrates[i]
Combinedname = paste(Scanrate,"_3", sep="")
Namelist = c(Namelist,Combinedname)
}

#The code below is a long loop that shows the oxidative sweeps of the CV and allows you to isolate the oxidative faradic peaks via clipping and polynomical fitting
Funclist = list()
Extractedpeaksforward = list()

for(i in 1:length(Namelist)){
  Funclist[[i]] = Process(Namelist[i])
}

for(i in 1:length(Scanrates)){

message(Scanrates[i])  
Func1 = Funclist[[i]]
Forward=Chopf(Func1)
plot(Forward)
lines(Smoothed(Forward), col = "red")


Clipfrommax1forward <- readline(prompt="Enter Clipfrommax1forward: ")
Clipfrommax1forward = as.numeric(Clipfrommax1forward)
Clipfrommin1forward <- readline(prompt="Enter Clipfrommin1forward: ")
Clipfrommin1forward = as.numeric(Clipfrommin1forward)

Forward = Clip(Forward,Clipfrommax1forward,Clipfrommin1forward)
plot(Forward)
lines(Smoothed(Forward), col = "red")

Fitbaselineusing1forward <- readline(prompt="Enter Fitbaselineusing1forward: ")
Fitbaselineusing1forward = as.numeric(Fitbaselineusing1forward)
Fitbaselineusing2forward <- readline(prompt="Enter Fitbaselineusing2forward: ")
Fitbaselineusing2forward = as.numeric(Fitbaselineusing2forward)


Baseline = Clip2(Forward,Fitbaselineusing1forward, Fitbaselineusing2forward)
Voltagevalues= as.vector(as.matrix(Forward[,1])) 
Voltage=Baseline[,1]
Current = Baseline[,2]

Polyorder <- readline(prompt="Enter Polyorder: ")
Polyorder = as.numeric(Polyorder)

poly=lm(Current ~ poly(Voltage, Polyorder, raw=TRUE))
  
Baselineprediction=predict(poly, list(Voltage=Voltagevalues))
Baselinetoplot = data.frame(Voltagevalues,Baselineprediction)
  
  plot(Forward)
  lines(Smoothed(Forward), col = "red")
  lines(Baseline, col="blue")
  lines(Baselinetoplot, col="red")
  
  GraphName= paste("Ox",Scanrates[i],".png", sep="")
  
  png(file = GraphName)
  plot(Forward)
  lines(Baseline, col="blue")
  lines(Baselinetoplot, col="red")
  dev.off()
  
p <- readline(prompt="Press 1 to continue, 0 to go back: ") 

for(n in 1:100) {if(p == 1) {break}
else
Polyorder <- readline(prompt="Enter Polyorder: ")
Polyorder = as.numeric(Polyorder)

poly=lm(Current ~ poly(Voltage, Polyorder, raw=TRUE))

Baselineprediction=predict(poly, list(Voltage=Voltagevalues))
Baselinetoplot = data.frame(Voltagevalues,Baselineprediction)

plot(Forward)
lines(Baseline, col="blue")
lines(Baselinetoplot, col="red")

png(file = GraphName)
plot(Forward)
lines(Baseline, col="blue")
lines(Baselinetoplot, col="red")
dev.off()

p <- readline(prompt="Press 1 to continue, 0 to go back: ") 
}


Extractedsignalcurrent=as.vector(as.matrix(Forward[,2]-Baselinetoplot[,2]))
Extractedsignalforward=data.frame(Voltagevalues,Extractedsignalcurrent)

x=as.vector(Extractedsignalforward[,1])
y=as.vector(Extractedsignalforward[,2])

smoothingSpline = smooth.spline(x, y, spar=0.5)
predictvalues=predict(smoothingSpline, newdata=data.frame(x=x.new))
plot(x,y, type="l")
lines(predictvalues, col="red")

#message(Scanrates[i])
Extractedpeaksforward[[i]] = predictvalues
  
}
#The code below is a long loop that shows the reductive sweeps of the CV and allows you to isolate the reductive faradic peaks via clipping and polynomical fitting
Extractedpeaksbackwards = list()

for(i in 1:length(Namelist)){
  Funclist[[i]] = Process(Namelist[i])
}

for(i in 1:length(Scanrates)){


Func1 = Process(Namelist[i])
Backward = Chopb(Func1)
message(Scanrates[i])  
plot(Backward)
lines(Smoothed(Backward), col = "red")

Clipfrommax1back <- readline(prompt="Enter Clipfrommax1back: ")
Clipfrommax1back = as.numeric(Clipfrommax1back)
Clipfrommin1back <- readline(prompt="Enter Clipfrommin1back: ")
Clipfrommin1back = as.numeric(Clipfrommin1back)

Backward=Backward[order(nrow(Backward):1),]
Backward = Clip(Backward,Clipfrommax1back,Clipfrommin1back)
plot(Backward)
lines(Smoothed(Backward), col = "red")

Fitbaselineusing1back <- readline(prompt="Enter Fitbaselineusing1back: ")
Fitbaselineusing1back = as.numeric(Fitbaselineusing1back)
Fitbaselineusing2back <- readline(prompt="Enter Fitbaselineusing2back: ")
Fitbaselineusing2back = as.numeric(Fitbaselineusing2back)

Baseline = Clip2(Backward,Fitbaselineusing1back, Fitbaselineusing2back)
Voltagevalues= as.vector(as.matrix(Backward[,1])) 
Voltage=Baseline[,1]
Current = Baseline[,2]

Polyorder <- readline(prompt="Enter Polyorder: ")
Polyorder = as.numeric(Polyorder)

poly=lm(Current ~ poly(Voltage, Polyorder, raw=TRUE))

Baselineprediction=predict(poly, list(Voltage=Voltagevalues))
Baselinetoplot = data.frame(Voltagevalues,Baselineprediction)

plot(Backward)
lines(Smoothed(Backward), col = "red")
lines(Baseline, col="blue")
lines(Baselinetoplot, col="red")

  GraphName1= paste("Red",Scanrates[i],".png", sep="")

  png(file = GraphName1)
  plot(Backward)
  lines(Baseline, col="blue")
  lines(Baselinetoplot, col="red")
  dev.off()

p <- readline(prompt="Press 1 to continue, 0 to go back: ") 

for(n in 1:100) {if(p == 1) {break}
  else
    Polyorder <- readline(prompt="Enter Polyorder: ")
  Polyorder = as.numeric(Polyorder)
  
  poly=lm(Current ~ poly(Voltage, Polyorder, raw=TRUE))
  
  Baselineprediction=predict(poly, list(Voltage=Voltagevalues))
  Baselinetoplot = data.frame(Voltagevalues,Baselineprediction)
  
  plot(Backward)
  lines(Baseline, col="blue")
  lines(Baselinetoplot, col="red")
  
  png(file = GraphName1)
  plot(Backward)
  lines(Baseline, col="blue")
  lines(Baselinetoplot, col="red")
  dev.off()
  
  p <- readline(prompt="Press 1 to continue, 0 to go back: ") 
}


Extractedsignalcurrent=as.vector(as.matrix(Backward[,2]-Baselinetoplot[,2]))
Extractedsignalbackward=data.frame(Voltagevalues,Extractedsignalcurrent)

x=as.vector(Extractedsignalbackward[,1])
y=as.vector(Extractedsignalbackward[,2])

z = smooth.spline(x, y, spar=0.5)
predictvalues=predict(z, newdata=data.frame(x=x.new))
plot(x,y, type="l")
lines(predictvalues, col="red")

#message(Scanrates[i])
Extractedpeaksbackwards[[i]] = predictvalues

}

#The lines of code below analyse the data based on the extracted peaks

Dataframenames=c()
for (i in 1:length(Scanrates)){
  Dataframenames=c(Dataframenames,paste(Scanrates[i],"mVs-1", sep=""))
}

peakmaximaforwards=c()
for(i in 1:length(Extractedpeaksforward)){
  x=(as.data.frame(Extractedpeaksforward[i]))
  y=max(x[,2],na.rm=TRUE)
  peakmaximaforwards=c(peakmaximaforwards, y)
}
peakmaximabackwards=c()
for(i in 1:length(Extractedpeaksbackwards)){
  x=(as.data.frame(Extractedpeaksbackwards[i]))
  y=-min(x[,2],na.rm=TRUE)
  peakmaximabackwards=c(peakmaximabackwards, y)
}
Coverageforward=c()
Coveragebackwards=c()

for(i in 1:length(Scanrates)){
  Coverageforward=c(as.numeric(Findcoverage(as.data.frame(Extractedpeaksforward[[i]]),0.0706858347,Scanrates[i],1)),Coverageforward)
  Coveragebackwards=c(-as.numeric(Findcoverage(as.data.frame(Extractedpeaksbackwards[[i]]),0.0706858347,Scanrates[i],1)),Coveragebackwards)
}
Coverageforward=as.vector(Coverageforward)
Coveragebackwards=as.vector(Coveragebackwards)


Combinedcoverages=c(Coverageforward,Coveragebackwards)
Combinedcoverages_df=as.data.frame(Combinedcoverages)
Combinedcoverages_df[2]="Combined"
colnames(Combinedcoverages_df)=c("Coverage","Type")
Coverageforward_df=as.data.frame(Coverageforward)
Coverageforward_df[2]="Oxidative"
colnames(Coverageforward_df)=c("Coverage","Type")
Coveragebackwards_df=as.data.frame(Coveragebackwards)
Coveragebackwards_df[2]="Reductive"
colnames(Coveragebackwards_df)=c("Coverage","Type")


Coveragestoplot=rbind(Coverageforward_df, Combinedcoverages_df,Coveragebackwards_df)

anodicpeaks=c()
for(i in 1:length(Scanrates)){
  anodicpeak=Peakanodic(as.data.frame(Extractedpeaksforward[[i]]))
  anodicpeaks=c(anodicpeaks,anodicpeak)
}


cathodicpeaks=c()
for(i in 1:length(Scanrates)){
  cathodicpeak=Peakcathoic(as.data.frame(Extractedpeaksbackwards[[i]]))
  cathodicpeaks=c(cathodicpeaks,cathodicpeak)
}


peakseparationsr=anodicpeaks-cathodicpeaks
peakseparationsr=peakseparationsr*1000

peakmids=cathodicpeaks+anodicpeaks
peakmids=peakmids/2
peakmids=peakmids[-length(peakmids-1)]
peakmids=peakmids[-1]

Widthsathalfheightforward=c()
for(i in 1:length(Scanrates)){
  widthforward=findhalfheightforward(as.data.frame(Extractedpeaksforward[[i]]))
  Widthsathalfheightforward=c(Widthsathalfheightforward,widthforward)
}

Widthsathalfheightback=c()
for(i in 1:length(Scanrates)){
  widthback=findhalfheightback(as.data.frame(Extractedpeaksbackwards[[i]]))
  Widthsathalfheightback=c(Widthsathalfheightback,widthback)
}
Combinedwidthsathalfheight=c(Widthsathalfheightforward,Widthsathalfheightback)

Combinedwidthsathalfheight_df=as.data.frame(Combinedwidthsathalfheight)
Combinedwidthsathalfheight_df[2]="Combined"
colnames(Combinedwidthsathalfheight_df)=c("PWHH","Type")
Widthsathalfheightforward_df=as.data.frame(Widthsathalfheightforward)
Widthsathalfheightforward_df[2]="Oxidative"
colnames(Widthsathalfheightforward_df)=c("PWHH","Type")
Widthsathalfheightback_df=as.data.frame(Widthsathalfheightback)
Widthsathalfheightback_df[2]="Reductive"
colnames(Widthsathalfheightback_df)=c("PWHH","Type")


PWHHtoplot=rbind(Widthsathalfheightforward_df, Combinedwidthsathalfheight_df,Widthsathalfheightback_df)

rootScanrates=Scanrates^0.5

#A summary plot of your analysis is generated
par(mfrow=c(2,3))
plot(Extractedpeaksforward[[1]], type = "l", lwd=2,
     xlab="Voltage (V)", ylab="Current (µA)", ylim=c(-10,8), xlim=c(-0.3,0.9))
for(i in 2:length(Extractedpeaksforward)){
  lines(Extractedpeaksforward[[i]], col=Colours[i], lwd=2)
}
for(i in 1:length(Extractedpeaksbackwards)){
  lines(Extractedpeaksbackwards[[i]], col=Colours[i], lwd=2)
}

plot(Scanrates,peakmaximaforwards, col="red", ylab="Peak current (µA)", xlab="Scanrate (mV / s)", ylim=c(0,10), xlim=c(0,1000))
points(Scanrates,peakmaximabackwards, col="blue")
abline(lm(peakmaximaforwards~Scanrates), col="red")
abline(lm(peakmaximabackwards~Scanrates), col="blue")
legend(2, 40, legend=c("Oxidative","Reductive"),
       col=c("red","blue"), lwd=1, cex=0.9,bty="n")

plot(rootScanrates,peakmaximaforwards, col="red", ylab="Peak current (µA)", xlab="Scanrate^0.5 (mV / s)^0.5", ylim=c(0,10),xlim=c(0,35))
points(rootScanrates,peakmaximabackwards, col="blue")
abline(lm(peakmaximaforwards~rootScanrates), col="red")
abline(lm(peakmaximabackwards~rootScanrates), col="blue")
legend(2, 40, legend=c("Oxidative","Reductive"),
       col=c("red","blue"), lwd=1, cex=0.9,bty="n")

plot(Scanrates/1000,peakseparationsr, xlab="Scanrate (V / s)", ylab = "Peak separation (mV)",ylim=c(0,600))

boxplot(Coverage~Type, data=Coveragestoplot, ylab="Estimated coverage (pmol /cm2)", col=c("purple","Red","Blue"),ylim=c(100,500)) 

boxplot(peakmids, col="green", ylab="Peak midpoint (V)", ylim=c(0.25, 0.4))

#The code below can be used export xy data for the extracted peaks if desired
for (h in 1:length(Scanrates)){
Oxfile = Scanrates[h]
Redfile = Scanrates[h]

ExportOx=as.data.frame(Extractedpeaksforward[h])
colnames(ExportOx) <- c("Voltage (V)", " Current (µA)")
ExportRed=as.data.frame(Extractedpeaksbackwards[h])
colnames(ExportRed) <- c("Voltage (V)", " Current (µA)")
write.table(ExportOx,Oxfile)
 write.table(ExportRed,Redfile)}
