#
# ====================================================================================
#
# define get.dim function (de Jerome S)

"get.dim"=
function(nc,dimname,varname=NULL,verbose=2){

# L'interet de cette fonction est de recuperer systematiquement la variable qui correspond a une dimension classique, comme la latitude

dimnames=names(nc$dim)
if (is.null(varname)==FALSE){
dimnames=dimnames[nc$var[[varname]]$dimids]
}#end if


if (dimname=="x" | dimname=="longitude" | dimname=="lon" | dimname=="LON" | dimname=="X"){
   names=c("x","X","lon","long","longitude","LON","LONGITUDE","Lon","Longitude")
   dum=c()
   for (l in 1:length(dimnames)){
     dum=c(dum,(any(names == dimnames[l])))
   }#end for l
   nyname=dimnames[which(dum==TRUE)]
   if (verbose==2) print(nyname)
   vals=nc$dim[[nyname]]$vals
}#end if

if (dimname=="lat" | dimname=="latitude" | dimname=="y" | dimname=="Y" | dimname=="LAT"){
   names=c("y","Y","lat","lati","latitude","LAT","LATITUDE","Lat","Latitude")
   dum=c()
   for (l in 1:length(dimnames)){
     dum=c(dum,(any(names == dimnames[l])))
   }#end for l
   nyname=dimnames[which(dum==TRUE)]
   if (verbose==2) print(nyname)
   vals=nc$dim[[nyname]]$vals
}#end if

if (dimname=="z" | dimname=="vertical" | dimname=="levels" | dimname=="level" | dimname=="depth" | dimname=="deptht" ){
   names=c("z","Z","presnivs","deptht","depth","depthw","pressure","lev","level","levels","DEPTH",
           "sigma","sigma0","sigma1","sigma2","sigma3","sigma4","SIGMA","AXSIGMA","axsigma", "DEPTH1_28")
   dum=c()
   for (l in 1:length(dimnames)){
    dum=c(dum,(any(names == dimnames[l])))
   }#end for l
   nyname=dimnames[which(dum==TRUE)]
   if (verbose==2) print(nyname)
   if ( length (nyname) > 0 ) {
       vals=nc$dim[[nyname]]$vals
   } else {
       vals=""
   }
}#end if


if (dimname=="t" | dimname=="time" | dimname=="time_counter" | dimname=="time"){
   names=c("t","T","time","time_counter","TIME_COUNTER","TIME","TMEANSC")
   dum=c()
   for (l in 1:length(dimnames)){
    dum=c(dum,(any(names == dimnames[l])))
   }#end for l
   nyname=dimnames[which(dum==TRUE)]
   if (verbose==2) print(nyname)
   vals=nc$dim[[nyname]]$vals
}#end if

vals

}#end def get.dim
