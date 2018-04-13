LNGXE3418J <- function(
profileList,
profileID
){
mz<-(as.numeric(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),1]))
int<-(as.numeric(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),2]))
plot(mz,int,pch=19,cex=.7,xlab="m/z",ylab="Intensity")
abline(v=mean(mz),lwd=2,lty=2,col="lightblue")
}
