
remo<-list.files(file.path(logfile[[1]],"results","componentization","isotopologues"))
if(length(remo)>0){
	for(i in 1:length(remo)){
		file.remove(file.path(logfile[[1]],"results","componentization","isotopologues",remo[i])) 
	}
}


