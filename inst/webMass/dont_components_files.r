remo<-list.files(file.path(logfile[[1]],"results","componentization","components"))
if(length(remo)>0){
	for(i in 1:length(remo)){
		file.remove(file.path(logfile[[1]],"results","componentization","components",remo[i])) 
	}
}
