init <- function(){
  if(!file.exists(paste(main,"/","work data",sep=""))) dir.create(path=paste(main,"/","work data",sep="")) 
  if(!file.exists(paste(main,"/","tables",sep=""))) dir.create(path=paste(main,"/","tables",sep=""))
  if(!file.exists(paste(main,"/","figures",sep=""))) dir.create(path=paste(main,"/","figures",sep=""))
  
  subdirs <- list(scripts = paste(main, "scripts", sep = "/"),
                  org.DEUII = paste(main, "original data sets", "DEU2", sep = "/"),
                  leg.cycles = paste(main, "original data sets", "Legislative Cycles", sep = "/"),
                  council.voting = paste(main, "original data sets", "Council Voting Weights", sep = "/"),
                  corrected.EAs = paste(main, "original data sets", "EAsCorrected", sep = "/"),
                  same.country = paste(main, "original data sets", "ConflictOfInterest", sep = "/"),
                  type.data = paste(main, "original data sets", "type of file", sep = "/"),
                  ea.data = paste(main, "original data sets", "Early Agreements", sep = "/"),
                  workdata = paste(main, "work data", sep = "/"),
                  graphs = paste(main, "figures", sep = "/"),
                  tables = paste(main, "tables", sep = "/"))
  return(subdirs)
}