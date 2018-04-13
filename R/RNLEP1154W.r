RNLEP1154W <- function(
MSlist,
mz,
dmz=5,
ppm=TRUE,
RT,
dRT,
int,
kdTree=FALSE
){
if(length(dmz)==1){dmz<-rep(dmz,length(mz))}
if(length(dRT)==1){dRT<-rep(dRT,length(mz))}
result<-rep("FALSE",length(mz));
options(digits=10)
inter<-interactive()
if(kdTree[1]=="FALSE"){
cat("\n Build kdtree ...");
if(inter){pbar <- txtProgressBar( min = 0, max = length(MSlist[[4]][[2]][,1]), style = 3 )}
rawTree<-.Call("kdtree4",
as.matrix(MSlist[[4]][[2]][,1:3]),
as.integer(inter),
pbar,
PACKAGE="nontarget"
);
if(inter){close(pbar)}
cat("done.")
}
bounds<-matrix(ncol=2,nrow=3,0)
bounds[2,2]<-max(MSlist[[4]][[2]][,2])
cat("\n Screen raw data...");
if(inter){pbar <- txtProgressBar( min = 0, max = length(mz), style = 3 )}
for(i in 1:length(mz)){
if(inter){setTxtProgressBar(pbar,i,title = NULL, label = NULL)}
if(ppm){
bounds[1,1]<-(mz[i]-(mz[i]*dmz[i]/1E6))
bounds[1,2]<-(mz[i]+(mz[i]*dmz[i]/1E6))
}else{
bounds[1,1]<-(mz[i]-dmz)
bounds[1,2]<-(mz[i]+dmz)
}
bounds[2,1]<-int[i]
bounds[3,1]<-(RT[i]-dRT[i])
bounds[3,2]<-(RT[i]+dRT[i])
found<-.Call("search_kdtree",
MSlist[[4]][[2]][,1:3],
rawTree,
bounds,
PACKAGE="nontarget"
)
if(length(found)>0){
result[i]<-"TRUE"
}
}
if(inter){close(pbar)};
cat("done.")
return(result)
}
