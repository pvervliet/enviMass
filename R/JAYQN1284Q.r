JAYQN1284Q <- function(pro_name, pro_dir){
say_path<-"Project path ok"
skip_rest<-FALSE
if(!file.exists(pro_dir)){
say_path<-paste("File or folder ",pro_dir," does not exist!")
skip_rest<-TRUE
}
if(file.access(pro_dir, mode = 0)!=0 & !skip_rest){
say_path<-paste("File or folder ",pro_dir," does not exist!")
}
if(file.access(pro_dir, mode = 2)!=0 & !skip_rest){
say_path<-paste("Not allowed to write into ",pro_dir,". Revise user permissions?")
}
if(file.access(pro_dir, mode = 4)!=0 & !skip_rest){
say_path<-paste("Not allowed to read from ",pro_dir,". Revise user permissions?")
}
a<-try({
dir.create(
file.path(pro_dir,pro_name)
)
},silent=TRUE)
if((a==FALSE) || (class(a)=="try-error" )){
say_path<-paste("Cannot create ",pro_dir,". Revise user permissions, path validity, ... ?")
}else{
if(FALSE){
file.remove(file.path(pro_dir,pro_name))
}
}
return(say_path);
}
