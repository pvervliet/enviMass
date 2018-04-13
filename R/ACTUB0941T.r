ACTUB0941T <- function(Date,increment){
if(!is.character(Date)){stop("Date not a character string.")}
if(!is.numeric(increment)){stop("increment not numeric")}
new_Date<-as.POSIXct(Date)
new_Date<-(new_Date+(increment*60*60*26))
new_Date<-strsplit(as.character(new_Date)," ")[[1]][1]
return(new_Date);
}
