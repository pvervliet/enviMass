GSZEL3799Y <- function(a,b,row_order=FALSE,column_order_a=TRUE,column_order_b=TRUE,get_index=FALSE){
if(dim(a)[2]!=dim(b)[2]){
stop("\n Matrices a and b have different number of columns?")
}
len1<-dim(a)[1]
if(row_order){
a<-t(apply(t(a), 2, sort))
b<-t(apply(t(b), 2, sort))
}
if(column_order_a){
ord1<-do.call(order,as.data.frame(a))
a<-a[ord1,]
}
len2<-dim(b)[1]
if(column_order_b){
ord2<-do.call(order,as.data.frame(b))
b<-b[ord2,]
}
if(FALSE){
found<-rep(0,len1)
dims<-dim(a)[2]
do_dims<-c(1:dim(a)[2])
at<-1
for(i in 1:len1){
for(m in do_dims){
while(b[at,m]<a[i,m]){
at<-(at+1)
if(at>len2){break}
}
if(at>len2){break}
if(b[at,m]>a[i,m]){break}
if(m==dims){found[i]<-at}
}
if(at>len2){break}
}
}
if(TRUE){
results<-rep(0,len1)
found <- .Call("_enviMass_compare",
as.matrix(a),
as.matrix(b),
as.matrix(results),
PACKAGE="enviMass"
)
}
if(!get_index){
if(!column_order_a & !column_order_b){
return(found>0)
}
if(column_order_a & !column_order_b){
found<-found[order(ord1)]
return(found>0)
}
if(column_order_a & column_order_b){
found<-found[order(ord1)]
found<-order(ord2)[found]
return(found>0)
}
}else{
if(!column_order_a & !column_order_b){
return(found)
}
if(column_order_a & !column_order_b){
found<-found[order(ord1)]
return(found)
}
if(column_order_a & column_order_b){
found<-found[order(ord1)]
found<-ord2[found]
return(found)
}
}
}
