init <- function(){
  if(!file.exists(paste(main,"/","work data",sep=""))) dir.create(path=paste(main,"/","work data",sep="")) 
  if(!file.exists(paste(main,"/","tables",sep=""))) dir.create(path=paste(main,"/","tables",sep=""))
  if(!file.exists(paste(main,"/","figures",sep=""))) dir.create(path=paste(main,"/","figures",sep=""))
}
init()