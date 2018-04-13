MAKLM3979Y <- function(
profileList,
profpeaks2,
progbar=FALSE
){
if(!profileList[[1]][[3]]){stop("profileList not profiled; aborted.")}
atPOSIX<-profileList[[3]];
sampletype<-profileList[[9]];
sampleID<-profileList[[4]];
keep<-((sampletype=="sample")|(sampletype=="blank"))
atPOSIX<-atPOSIX[keep]
sampletype<-sampletype[keep]
sampleID<-sampleID[keep]
atdate<-c();
attime<-c();
for(i in 1:length(atPOSIX)){
atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
}
attime<-as.difftime(attime);
atdate<-as.Date(atdate, tz="GMT");
ord<-order(as.numeric(atdate),as.numeric(attime),sampleID);
atPOSIXsort<-atPOSIX[ord];
atdate<-atdate[ord];
attime<-attime[ord];
sampleID<-sampleID[ord];
sampletype<-sampletype[ord];
timeset<-matrix(nrow=length(atPOSIX),ncol=5,0);
for(i in 1:length(sampleID)){
if(sampletype[i]=="sample"){
timeset[i,2]<-as.numeric(sampleID[i]);
}
if(sampletype[i]=="blank"){
timeset[i,3]<-as.numeric(sampleID[i]);
}
}
numtime<-(as.numeric(atdate)+as.numeric(attime/(24*60*60)))
colnames(timeset)<-c("above blank?","sampleID","blankID","sample_int","blank_int")
leng_list<-length(profpeaks2[,1])
leng<-length(timeset[,1])
export<-matrix(ncol=(3+(2*leng_list)),nrow=(leng+3),0)
export[4:(leng+3),c(2,3)]<-timeset[,c(2,3)]
export[4:(leng+3),1]<-numtime
if(progbar==TRUE){  prog<-winProgressBar("Extract intensity differences...",min=1,max=leng_list);
setWinProgressBar(prog, 0, title = "Extract intensity differences...", label = NULL);}
for(k in 1:leng_list){
if(progbar==TRUE){ setWinProgressBar(prog, k, title = "Extract intensity differences...", label = NULL) }
ID<-as.numeric(profpeaks2[k,10])
getit<-seq(1,length(profileList[[7]][,4]),1)[profileList[[7]][,4]==ID]
if(length(getit)>1){stop("Debug me - different profiles with same IDs found!")}
timeset[,c(4,5)]<-0;
timeset[,c(4,5)] <-.Call("_enviMass_fill_timeset",
as.numeric(timeset),
as.numeric(profileList[[2]][(profileList[[7]][getit,1]:profileList[[7]][getit,2]),6]),
as.numeric(profileList[[2]][(profileList[[7]][getit,1]:profileList[[7]][getit,2]),2]),
as.integer(length(timeset[,1])),
PACKAGE="enviMass"
)
at<-(3+(2*(k-1)+c(1,2)))
export[4:(leng+3),at]<-timeset[,c(4,5)]
export[1,at[1]]<-profileList[[7]][getit,4][[1]]
export[2,at[1]]<-profileList[[7]][getit,14][[1]]
export[3,at[1]]<-profileList[[7]][getit,15][[1]]
}
if(progbar==TRUE){close(prog)};
return(export)
}
