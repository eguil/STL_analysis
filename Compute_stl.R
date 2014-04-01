#!/usr/bin/env Rscript
#
# ----------------------------------------------------------------------------------------
#  Procedure to compute STL decomposition of monthly serie of 2D (and soon 3D) fields
#  Called by scr_stl.sh (to be translated in python as an exercise)
# ----------------------------------------------------------------------------------------
#  E. Guilyardi (2014) with help from J. Servonnat
#
# Examples using sh call script
# 
# scr_stl.sh -v sosstsst -s 13 -f 7 -t 373,1344 -l sosst/SST/DegC -o ../STL_out HadISST1_1m_187001_201212_reg1m.nc 
# scr_stl.sh -v sozotaux -s 13 -f 7 -l sozot/Taux/N/m2 -o ../STL_out ERA40_1m_195801_200112_grid_U.nc 
# scr_stl.sh -v thetao -s13 -f7 -l votemp/thetao/DegC,m -t1,48 -o ../STL_out ORAS4_1m_195801_201112_gridT_thetao.nc 
#  Hadisst 1870_2012: 1900-2012 -> init=373 / countt=1344
#
# -- Packages needed
msg.trap <- capture.output( suppressMessages(library(ncdf4))) # netcDF 4
msg.trap <- capture.output( suppressMessages(library(car))) # misc
msg.trap <- capture.output( suppressMessages(library(zoo))) # misc
msg.trap <- capture.output( suppressMessages(library(doMC))) # parallel library
msg.trap <- capture.output( suppressMessages(library(matrixStats))) # matrix stats

host=Sys.info()[['nodename']]
if ( host == 'crunchy.llnl.gov' ) {
    source ('/work/guilyardi/STL_analysis/support_nc.R') # custom fonctions
} else {
    source('/Users/ericg/Projets/SF_2013-2014/ENSO_AC_div/STL_analysis/support_nc.R') # custom fonctions
}
# -- Arguments
#
# <file> <var> <iniyear> <init_idx> <count> <STL seasonal window> <STL trend window>
#./Compute_stl.R /Users/ericg/Desktop/Data/database/HadISST1_1m_187001_201212_reg1m.nc sosstsst 1900 373 1344 13 241
# needs to be executed in <dir> of <file>

# Recover arguments

args=commandArgs(trailingOnly = T)

file=args[1]
varname=args[2]
iniyear=as.numeric(args[3])
init=as.numeric(args[4])
countt=as.numeric(args[5])
seasonalwin=as.numeric(args[6])
trendwin=as.numeric(args[7])
debug=as.numeric(args[8])
varout=unlist(strsplit(args[9],"/"))

if ( debug == 1 ) {
    print(" Arguments:")
    print(args)
}


# ------------------------------------------------------------- #
# -- Open file and retrieve data and dimensions              -- #
# ------------------------------------------------------------- #

# -- On ouvre le fichier netcdf avec la commande open.ncdf
nc = nc_open(file)
# -- En interactif, tu peux regarder ce qu'il y a dans l'objet nc (en tapant juste nc dans le prompt, et enter)
#nc

# -- Tu peux aussi parcourir le fichier avec la commande names() pour savoir ce qu'il y a dans l'objet nc
if ( debug == 1 ) {
    cat(" Contents of file ", file)
    cat(names(nc),"\n")
    cat(names(nc$var),"\n")
    cat(names(nc$dim),"\n")
}
# 2Dxtime or 3Dxtime data ?
dep_name=sub('.*\"(.*)\"','\\1',capture.output(z=get.dim(nc, "depth"))[1])
if (debug == 1) { cat( " dep_name =",dep_name,"\n") }
if ( dep_name != "character(0)" ) {
  # 3D data 
    start <- c(1,1,1,init)
    count <- c(-1,-1,-1,countt)
    sw3d=1
    if (debug == 1) {
#        start <- c(20,60,1,init)
#        count <- c(40,80,10,countt)
#        start <- c(1,20,1,init)
#        count <- c(20,25,10,countt)
    }
    
} else {
  # 2D data
    start <- c(1,1,init)    

    count <- c(-1,-1,countt)
    sw3d=0
}
if ( debug == 1) {
    cat(" 3D data ? =",sw3d,"\n")
    cat(" Start/count =",start, count,"\n")
}

# read data array
dat = ncvar_get(nc,varid=varname, start=start,count=count)
if (debug == 1) { cat( " dim(dat) =",dim(dat),"\n") }


# -- Get dimensions
x=get.dim(nc, "lon")
lon_name=sub('.*\"(.*)\"','\\1',capture.output(x=get.dim(nc, "lon"))[1])
y=get.dim(nc, "lat")
lat_name=sub('.*\"(.*)\"','\\1',capture.output(y=get.dim(nc, "lat"))[1])
nav_lon = ncvar_get(nc,varid=lon_name)
nav_lat = ncvar_get(nc,varid=lat_name)
if ( sw3d == 1 ) {
    z=get.dim(nc, dep_name)
    depth = ncvar_get(nc,varid=dep_name)
}
tim_name=sub('.*\"(.*)\"','\\1',capture.output(t=get.dim(nc, "time"))[1])
time = ncvar_get(nc,varid=tim_name, start=init,count=countt)


if (debug == 1) {
    cat( " Dimension of nav_lon, nav_lat, time =", dim(nav_lon), dim(nav_lat), dim(time),"\n")
    if ( sw3d == 1) {cat( " Dimension of depth =", dim(depth),"\n")}
}
# 1D TAO-like data ?
if ((dim(nav_lon) == 1) & (dim(nav_lat)==1)) {
    sw1d = 1
} else {sw1d = 0}

# collapse to 1Dxtime

datc=dat
if (sw3d == 0) {
    ntimes=dim(dat)[3]
    ngrid=dim(dat)[1]*dim(dat)[2]
} else {
    if (sw1d == 1) {
        ntimes=dim(dat)[2]
        ngrid=dim(dat)[1]
    } else {
        ntimes=dim(dat)[4]
        ngrid=dim(dat)[1]*dim(dat)[2]*dim(dat)[3]
    }
}
if (debug == 1) { cat( " ngrid, ntimes =",ngrid,ntimes,"\n") }
    dim(datc) <- c(ngrid,ntimes)

# field already masked with NA ?
if (length(which(is.na(dat) == TRUE)) == 0) { 

# get missval 
    missval=nc$var[[varname]]$missval

    if (debug == 1) { cat( " missval=", missval,"\n") }
    if (nchar(missval) == 0 ){ stop(" *** missval not found ***") }

# for HadISST... BUG in reading missval
#missval=-missval 
#if (missval < 0 ) dat[which(dat < missval/10)] <- -missval
#missval=-missval # for HadISST...


# index of non masked points (only make STL on these points)
  idx <- which(datc[,1] < missval/10)
  mask <- rep(0,ngrid)
  mask[idx]=1
# make sure all dates have same mask
#for (i in 1:length(idx)) {
#  idxti <- which(datc[idx[i],] > missval/10)
#  if (length(idxti) > 0) mask[idx[i]]=0
#}
  idxm <- which(mask == 0)
  idxp <- which(mask == 1)
  datc[idxm,] <- NA
} else {
  idxm <- which(is.na(datc[,1]) == TRUE)
  idxp <- which(is.na(datc[,1]) == FALSE)
}

# output arrays
datstl1c=datc*NA    # annual harmonic
datstl2c=datc*NA    # trend
datstl3c=datc*NA    # residual (interannual)

# -------------
#  Compute STL
# -------------
mindim=1
maxdim=length(idxp)

if (debug == 1) {
    mindim=1 # for test (tropical pac with reg1m HadISST1)
    maxdim=400 # for test   
}
cat(" Number of valid grid points: ", length(idxp[mindim:maxdim]), "vs." ,ngrid,"(",length(idxp[mindim:maxdim])/ngrid*100,"%) \n")
# Parallel set up    
registerDoMC()
paral=1
cat(" Number of processor used: ", getDoParWorkers(), "\n")
# time process
ptm <- proc.time()
if (paral == 0) {
# loop on non masked points
    # sequential
    for (i in idxp[mindim:maxdim]) {
        tmp=stl(ts(datc[i,], start=c(iniyear,1), frequency=12), s.window=seasonalwin, robust=TRUE, t.window=trendwin, na.action=na.locf)
        annual=tmp$time.series[,"seasonal"]
        trend=tmp$time.series[,"trend"]
        remaind=tmp$time.series[,"remainder"]
        datstl1c[i,]=as.vector(annual)
        datstl2c[i,]=as.vector(trend)
        datstl3c[i,]=as.vector(remaind)
    }
} else {
    # parallel

    stlpar=foreach(i=1:length(idxp[mindim:maxdim]), .combine=cbind,.inorder=TRUE) %dopar% {
        tmp=stl(ts(datc[idxp[mindim+i-1],], start=c(iniyear,1), frequency=12), s.window=seasonalwin, robust=TRUE, t.window=trendwin, na.action=na.locf)
        annual=tmp$time.series[,"seasonal"]
        trend=tmp$time.series[,"trend"]
        remaind=tmp$time.series[,"remainder"]
        cbind(as.vector(annual),as.vector(trend),as.vector(remaind))
    }
    stlpar=aperm(stlpar)
    
    datstl1c[idxp[mindim:maxdim],]=stlpar[seq(1, nrow(stlpar), by=3),]
    datstl2c[idxp[mindim:maxdim],]=stlpar[seq(2, nrow(stlpar), by=3),]
    datstl3c[idxp[mindim:maxdim],]=stlpar[seq(3, nrow(stlpar), by=3),]
}

proc.time() - ptm

# -----------------------------
#  Mean AC and STD computation
# -----------------------------
acmea <- c(ngrid)      # AC cycle mean amplitude
stdam <- c(ngrid)      # std of annual harmonic - mean AC = modulation AC cycle
stdre <- c(ngrid)      # std of remainder (interannual variations)

# compute AC amplitude for each year and average to have mean AC amplitude

nyears=ntimes/12
dim(datstl1c) <- c(dim(datstl1c)[1], 12, nyears)
acmeai=datstl1c[,1,]*NA
dim(acmeai) <- c(dim(datstl1c)[1], nyears)


for (iy in 1:nyears){
  daty=datstl1c[,,iy]
  dim(daty) <- c(dim(datstl1c)[1], 12)
  # stdev
   acmeai[,iy]=rowSds(daty, na.rm = TRUE)
  # max-min
  #acmeai[,iy]=max(daty, na.rm = TRUE)-min(daty, na.rm = TRUE)
  
}
acmea=rowMeans(acmeai, na.rm=TRUE, dims=1)

dim(datstl1c) <- c(dim(datstl1c)[1], ntimes)

# compute std(remain) and stdev (seasonal minus mean seasonal) 

# compute std (remain)
stdre=rowSds(datstl3c, na.rm = TRUE)
# compute mean AC
meanac=foreach(i=1:12, .combine=cbind, .inorder=TRUE) %dopar% { 
    rowMeans(datstl1c[,seq(i,ncol(datstl1c), by=12)])
}
# replicate nyears times
meanacmult=apply(meanac, 2, rep, nyears)
dim(meanacmult) <- dim(datstl1c)
# remove from annual harmonic and take stddev
stdam=rowSds(datstl1c-meanacmult, na.rm = TRUE)

# compute ratio of stddev (remain) to AC amplitude
strat=stdre/acmea
#print(strat[18500])


# re-organise in 3D/2D including time dimension
if ( sw3d == 0 ) {
    dim(datstl1c) <- c(dim(dat)[1],dim(dat)[2],dim(dat)[3])
    dim(datstl2c) <- c(dim(dat)[1],dim(dat)[2],dim(dat)[3])
    dim(datstl3c) <- c(dim(dat)[1],dim(dat)[2],dim(dat)[3])
    
    dim(acmea) <- c(dim(dat)[1],dim(dat)[2],1)
    dim(stdam) <- c(dim(dat)[1],dim(dat)[2],1)
    dim(stdre) <- c(dim(dat)[1],dim(dat)[2],1)
    dim(strat) <- c(dim(dat)[1],dim(dat)[2],1)
} else {
    dim(datstl1c) <- c(dim(dat)[1],dim(dat)[2],dim(dat)[3],dim(dat)[4])
    dim(datstl2c) <- c(dim(dat)[1],dim(dat)[2],dim(dat)[3],dim(dat)[4])
    dim(datstl3c) <- c(dim(dat)[1],dim(dat)[2],dim(dat)[3],dim(dat)[4])
    
    dim(acmea) <- c(dim(dat)[1],dim(dat)[2],dim(dat)[3],1)
    dim(stdam) <- c(dim(dat)[1],dim(dat)[2],dim(dat)[3],1)
    dim(stdre) <- c(dim(dat)[1],dim(dat)[2],dim(dat)[3],1)
    dim(strat) <- c(dim(dat)[1],dim(dat)[2],dim(dat)[3],1)
}

if (debug == 1) {
    cat( " Dim STL arrays=", dim(datstl1c), dim(datstl2c), dim(datstl3c),"\n")
    cat( " Dim Mean STL arrays=", dim(acmea), dim(stdam), dim(stdre), dim(strat),"\n")
}


#plot(datstl1c[20,90,], type="l")
#contour(nav_lon, rev(nav_lat),datstl1c[,,30],nlevels=10)

#(dats <- ts(dat[30,90,], start=c(1980,1), frequency=12))
#plot(stl(dats,"per"))
#rm(list=ls())


outdir=""

# ------------------------------------------------------------- #
# -- NETCDF OUTPUT                                           -- #
# ------------------------------------------------------------- #
# -- On commence par definir les dimensions

TIME=ncdim_def("time_counter","time_counter",time,create_dimvar=TRUE,unlim=TRUE)

# -- On cree un objet varnc qui va accueillir les definitions des dimensions et des variables
varnc=list()
# puis pour nav_lon et nav_lat
if (length(dim(nav_lon)) == 1) {  # 1D dimension arrays
    X=ncdim_def("x",names(nc$dim)[1],x,create_dimvar=TRUE)
    Y=ncdim_def("y",names(nc$dim)[2],y,create_dimvar=TRUE)
    varnc[["nav_lon"]]=ncvar_def("nav_lon","degE",list(X), 1*10^30)
    varnc[["nav_lat"]]=ncvar_def("nav_lat","degN",list(Y)  ,1*10^30)
} else { # 2D dimension arrays
    X=ncdim_def("x",names(nc$dim)[1],nav_lon[,1],create_dimvar=TRUE)
    Y=ncdim_def("y",names(nc$dim)[2],nav_lat[1,],create_dimvar=TRUE)
    varnc[["nav_lon"]]=ncvar_def("nav_lon","degE",list(X,Y),1*10^30)
    varnc[["nav_lat"]]=ncvar_def("nav_lat","degN",list(X,Y),1*10^30)
}
if ( sw3d == 1 ) {
    Z=ncdim_def("z",names(nc$dim)[3],z,create_dimvar=TRUE)
    varnc[["depth"]]=ncvar_def("depth","m",list(Z),1*10^30)
}
# are we writing th 1m file ?
#if (nchar(varout[4]) == 0) {

    cat (" Write file 1m","\n")
# Pour les STL...

    var1=paste(varout[1],"ac",sep="")
    var2=paste(varout[1],"tr",sep="")
    var3=paste(varout[1],"re",sep="")
    if ( sw3d == 0 ) {
        dims=list(X,Y,TIME)
    } else {
        dims=list(X,Y,Z,TIME)
    }
    varnc[[var1]]=ncvar_def(var1,varout[3],dims,1*10^30, longname=paste("STL",varout[2],"seasonal"))
    varnc[[var2]]=ncvar_def(var2,varout[3],dims,1*10^30, longname=paste("STL",varout[2],"trend"))
    varnc[[var3]]=ncvar_def(var3,varout[3],dims,1*10^30, longname=paste("STL",varout[2],"remainder"))

    fileout="Out_stl.nc"

    if ( debug == 1 ) {
        cat (" lon, lat & time dim= ", dim(nav_lon), dim(nav_lat), dim(time),"\n")
        if (sw3d == 1) {cat (" depth dim= ", dim(depth),"\n")}
        cat (" varnc=",names(varnc),"\n")
        cat (" create file...\n")
    }

    ncnew=nc_create(paste(outdir,fileout,sep=""),varnc)

# -- On met les donnees des variables dans le fichier netcdf qu'on vient de creer

    if ( debug == 1 ) { cat (" write file...\n") }
    ncvar_put(ncnew,var1,datstl1c,start=NA,count=dim(datstl1c))
    ncvar_put(ncnew,var2,datstl2c,start=NA,count=dim(datstl2c))
    ncvar_put(ncnew,var3,datstl3c,start=NA,count=dim(datstl3c))
    ncvar_put(ncnew,"nav_lon",nav_lon,start=NA,count=dim(nav_lon))
    ncvar_put(ncnew,"nav_lat",nav_lat,start=NA,count=dim(nav_lat))
    if ( sw3d == 1 ) {ncvar_put(ncnew,"depth",depth,start=NA,count=dim(depth))}
    
# -- Et on ferme le fichier
    nc_close(ncnew)
#}
# ---------------------------------------------------
# -- store stds and ratio in new file
# ---------------------------------------------------
#
cat (" Write file mean","\n")


time <- c(1)
#
#
TIME=ncdim_def("time_counter","time_counter",time,create_dimvar=TRUE,unlim=TRUE) # On rajoute unlim=TRUE pour avoir une dimension UNLIMITED
if ( sw3d == 0 ) {
    dims=list(X,Y,TIME)
} else {
    dims=list(X,Y,Z,TIME)
}

varncs=list()
# Pour les STD + AC amplitude
var1=paste(varout[1],"acd",sep="")
var2=paste(varout[1],"acm",sep="")
var3=paste(varout[1],"dev",sep="")
var4=paste(varout[1],"rat",sep="")
varncs[[var1]]=ncvar_def(var1,varout[3],dims,1*10^30,longname=paste(varout[2],"AC mean amplitude"))
varncs[[var2]]=ncvar_def(var2,varout[3],dims,1*10^30,longname=paste(varout[2],"AC modulation"    ))
varncs[[var3]]=ncvar_def(var3,varout[3],dims,1*10^30,longname=paste(varout[2],"Std Dev"          ))
varncs[[var4]]=ncvar_def(var4," "      ,dims,1*10^30,longname=paste(varout[2],"Std Dev Int-AC ratio"   ))
varncs[["nav_lon"]]=varnc[["nav_lon"]]
varncs[["nav_lat"]]=varnc[["nav_lat"]]
if ( sw3d == 1 ) { varncs[["depth"]]=varnc[["depth"]] }

fileout="Out_std_stl.nc"
ncnews=nc_create(paste(outdir,fileout,sep=""),varncs)

ncvar_put(ncnews,var1,acmea,start=NA,count=dim(acmea))
ncvar_put(ncnews,var2,stdam,start=NA,count=dim(stdam))
ncvar_put(ncnews,var3,stdre,start=NA,count=dim(stdre))
ncvar_put(ncnews,var4,strat,start=NA,count=dim(strat))
ncvar_put(ncnews,"nav_lon",nav_lon,start=NA,count=dim(nav_lon))
ncvar_put(ncnews,"nav_lat",nav_lat,start=NA,count=dim(nav_lat))
if ( sw3d == 1 ) { ncvar_put(ncnews,"depth",depth,start=NA,count=dim(depth)) }

nc_close(ncnews)




# misc

# 360x180x 20 ans =
# MBA: 313.103   3.516 327.142
# cratos: 331.001       3.224     349.186
# 360x180x 40 ans =
# MBA: 452.761  23.186 553.522
# 100 ans: MBA 880.639  307.283 1784.621 

# for 4000 pts: 
#   user  system elapsed 
# 13.498   0.202  15.632
# 12.297   0.333  13.036 
# ORCA2 350 months
# user  system elapsed 
# 52.314   0.562  53.189 

# use multiproc
#maxdim=length(idxp)
#maxdim=4000 # for test
# time process
#ptm <- proc.time()
#stlpar=foreach(i=1:length(idxp[1:maxdim]),.combine=cbind,.inorder=TRUE) %dopar% {
#  tmp=stl(ts(datc[idxp[i],], start=c(1980,1), frequency=12), "per")
#  annual=tmp$time.series[,"seasonal"]
#  trend=tmp$time.series[,"trend"]
#  remain=tmp$time.series[,"remainder"]
#  cbind(annual,trend,remain)
#}

#datstl1c[idx[1:maxdim],]=stlpar[,seq(1,ncol(stlpar), by=3)]
#datstl2c[idx[1:maxdim],]=stlpar[,seq(2,ncol(stlpar), by=3)]
#datstl2c[idx[1:maxdim],]=stlpar[,seq(3,ncol(stlpar), by=3)]
# Stop the clock
#proc.time() - ptm
# for 4000 pts: 
#   user  system elapsed 
# 26.261   4.028  23.622 
# 10.156   3.733  23.621 

  #plot(stl(ts(dat[80,100,], start=c(1980,1), frequency=12), "per", robust=TRUE, t.window=500))
  #indian: plot(stl(ts(dat[250,90,600:1344], start=c(1950,1), frequency=12), s.window=13, t.window=241, na.action=na.locf,robust=TRUE))
  #t.window=npx1.5 a npx2
  #s.window:
   # Using R to perform STL decomposition, s.window controls how rapidly the seasonal component can change. Small values allow more rapid change. Setting the seasonal window to be infinite is equivalent to forcing the seasonal component to be periodic (i.e., identical across years).
  #t.window = trend window





