library(readxl)

if(data.all){
  sl.var <- c("GDP", "IP", "ESI", "car_registrations", "PMI", "EUR12","COVID cases","google mobility")
}else{
  sl.var <- c("GDP", "IP", "ESI", "car_registrations", "PMI", "EUR12")
}

if(cn.run=="DE"){
  DE <- readxl::read_excel("raw-data/DE.xlsx", sheet="data_set_clean")
  trans.codes <- as.numeric(DE[1,]); trans.codes[[1]] <- 3
  Y.m <- as.matrix(DE[2:nrow(DE),])
  rownames(Y.m) <- Y.m[,1]
  class(Y.m) <- "numeric"
  
  Y.raw <- Y.m <- Y.m[2:(nrow(Y.m)-1),sl.var]
  Y.raw[,1] <- Y.m[,1] <- 100*Y.raw[,1]
  Y.raw[is.nan(Y.raw)] <- NA
  
  # real time data
  RT_GDP <- read_excel("raw-data/DE_RT_GDP.xlsx", skip = 6)[-1,-2]
  RT_GDP[,-1] <- apply(RT_GDP[,-1],2,as.numeric)
  colnames(RT_GDP) <- c("date",as.character(seq(as.Date("2000-07-01"),as.Date("2020-09-01"),by="month")))
  RT_GDP <- RT_GDP[-nrow(RT_GDP),]
  RT_GDP <- RT_GDP[,c("date",as.character(seq(as.Date("2011-01-01"),as.Date("2020-06-01"),by="month")))]
  RT_GDP[,-1] <- rbind(NA,apply(RT_GDP[,-1],2,function(x){diff(log(x))}))*100
  RT_GDP <- as.matrix(RT_GDP[RT_GDP$date %in% c("Q2-2005","Q3-2005","Q4-2005",paste0("Q",1:4,"-",rep(2006:2019,each=4)),"Q1-2020","Q2-2020"),-1])
  
  RT_IP <- read_excel("raw-data/DE_RT_IP.xlsx", skip = 6)[-1,-2]
  RT_IP[,-1] <- apply(RT_IP[,-1],2,as.numeric)
  colnames(RT_IP) <- c("date",as.character(seq(as.Date("2000-03-01"),as.Date("2020-09-01"),by="month")))
  RT_IP <- RT_IP[-nrow(RT_IP),]
  RT_IP <- RT_IP[,c("date",as.character(seq(as.Date("2011-01-01"),as.Date("2020-06-01"),by="month")))]
  RT_IP[,-1] <- rbind(NA,apply(RT_IP[,-1],2,function(x){diff(log(x))}))*100
  IP.rt <- as.matrix(RT_IP[RT_IP$date %in% as.character(format(seq(as.Date("2005-04-01"),as.Date("2020-06-01"),by="month"),format="%b-%Y")),-1])
  
  GDP.rt <- kronecker(RT_GDP,matrix(1,nrow=3,ncol=1))
  GDP.rt[-seq(3,nrow(GDP.rt),by=3),] <- NA
}else if(cn.run=="ES"){
  ES <- readxl::read_excel("raw-data/ES.xlsx", sheet="data_set_clean")
  trans.codes <- as.numeric(ES[1,]); trans.codes[[1]] <- 3
  Y.m <- as.matrix(ES[2:nrow(ES),])
  rownames(Y.m) <- Y.m[,1]
  class(Y.m) <- "numeric"
  
  Y.raw <- Y.m <- Y.m[2:(nrow(Y.m)-1),sl.var]
  Y.raw[,1] <- Y.m[,1] <- 100*Y.raw[,1]
  
  # real time data
  RT_GDP <- read_excel("raw-data/ES_RT_GDP.xlsx", skip = 6)[-1,-2]
  RT_GDP[,-1] <- apply(RT_GDP[,-1],2,as.numeric)
  colnames(RT_GDP) <- c("date",as.character(seq(as.Date("2000-07-01"),as.Date("2020-09-01"),by="month")))
  RT_GDP <- RT_GDP[-nrow(RT_GDP),]
  RT_GDP <- RT_GDP[,c("date",as.character(seq(as.Date("2011-01-01"),as.Date("2020-06-01"),by="month")))]
  RT_GDP[,-1] <- rbind(NA,apply(RT_GDP[,-1],2,function(x){diff(log(x))}))*100
  RT_GDP <- as.matrix(RT_GDP[RT_GDP$date %in% c("Q2-2005","Q3-2005","Q4-2005",paste0("Q",1:4,"-",rep(2006:2019,each=4)),"Q1-2020","Q2-2020"),-1])
  
  RT_IP <- read_excel("raw-data/ES_RT_IP.xlsx", skip = 6)[-1,-2]
  RT_IP[,-1] <- apply(RT_IP[,-1],2,as.numeric)
  colnames(RT_IP) <- c("date",as.character(seq(as.Date("2000-03-01"),as.Date("2020-09-01"),by="month")))
  RT_IP <- RT_IP[-nrow(RT_IP),]
  RT_IP <- RT_IP[,c("date",as.character(seq(as.Date("2011-01-01"),as.Date("2020-06-01"),by="month")))]
  RT_IP[,-1] <- rbind(NA,apply(RT_IP[,-1],2,function(x){diff(log(x))}))*100
  IP.rt <- as.matrix(RT_IP[RT_IP$date %in% as.character(format(seq(as.Date("2005-04-01"),as.Date("2020-06-01"),by="month"),format="%b-%Y")),-1])
  
  GDP.rt <- kronecker(RT_GDP,matrix(1,nrow=3,ncol=1))
  GDP.rt[-seq(3,nrow(GDP.rt),by=3),] <- NA
}else if(cn.run=="FR"){
  FR <- readxl::read_excel("FR.xlsx", sheet="data_set_clean")
  trans.codes <- as.numeric(FR[1,]); trans.codes[[1]] <- 3
  Y.m <- as.matrix(FR[2:nrow(FR),])
  rownames(Y.m) <- Y.m[,1]
  class(Y.m) <- "numeric"
  
  Y.raw <- Y.m <- Y.m[2:(nrow(Y.m)-1),sl.var]
  Y.raw[,1] <- Y.m[,1] <- 100*Y.raw[,1]
  
  # real time data
  RT_GDP <- read_excel("raw-data/FR_RT_GDP.xlsx", skip = 6)[-1,-2]
  RT_GDP[,-1] <- apply(RT_GDP[,-1],2,as.numeric)
  colnames(RT_GDP) <- c("date",as.character(seq(as.Date("2000-07-01"),as.Date("2020-09-01"),by="month")))
  RT_GDP <- RT_GDP[-nrow(RT_GDP),]
  RT_GDP <- RT_GDP[,c("date",as.character(seq(as.Date("2011-01-01"),as.Date("2020-06-01"),by="month")))]
  RT_GDP[,-1] <- rbind(NA,apply(RT_GDP[,-1],2,function(x){diff(log(x))}))*100
  RT_GDP <- as.matrix(RT_GDP[RT_GDP$date %in% c("Q2-2005","Q3-2005","Q4-2005",paste0("Q",1:4,"-",rep(2006:2019,each=4)),"Q1-2020","Q2-2020"),-1])
  
  RT_IP <- read_excel("raw-data/FR_RT_IP.xlsx", skip = 6)[-1,-2]
  RT_IP[,-1] <- apply(RT_IP[,-1],2,as.numeric)
  colnames(RT_IP) <- c("date",as.character(seq(as.Date("2000-03-01"),as.Date("2020-09-01"),by="month")))
  RT_IP <- RT_IP[-nrow(RT_IP),]
  RT_IP <- RT_IP[,c("date",as.character(seq(as.Date("2011-01-01"),as.Date("2020-06-01"),by="month")))]
  RT_IP[,-1] <- rbind(NA,apply(RT_IP[,-1],2,function(x){diff(log(x))}))*100
  IP.rt <- as.matrix(RT_IP[RT_IP$date %in% as.character(format(seq(as.Date("2005-04-01"),as.Date("2020-06-01"),by="month"),format="%b-%Y")),-1])
  
  GDP.rt <- kronecker(RT_GDP,matrix(1,nrow=3,ncol=1))
  GDP.rt[-seq(3,nrow(GDP.rt),by=3),] <- NA
}else if(cn.run=="IT"){
  IT <- readxl::read_excel("raw-data/IT.xlsx", sheet="data_set_clean")
  trans.codes <- as.numeric(IT[1,]); trans.codes[[1]] <- 3
  Y.m <- as.matrix(IT[2:nrow(IT),])
  rownames(Y.m) <- Y.m[,1]
  class(Y.m) <- "numeric"
  
  Y.raw <- Y.m <- Y.m[2:(nrow(Y.m)-1),sl.var]
  Y.raw[,1] <- Y.m[,1] <- 100*Y.raw[,1]
  
  # real time data
  RT_GDP <- read_excel("raw-data/IT_RT_GDP.xlsx", skip = 6)[-1,-2]
  RT_GDP[,-1] <- apply(RT_GDP[,-1],2,as.numeric)
  colnames(RT_GDP) <- c("date",as.character(seq(as.Date("2000-07-01"),as.Date("2020-09-01"),by="month")))
  RT_GDP <- RT_GDP[-nrow(RT_GDP),]
  RT_GDP <- RT_GDP[,c("date",as.character(seq(as.Date("2011-01-01"),as.Date("2020-06-01"),by="month")))]
  RT_GDP[,-1] <- rbind(NA,apply(RT_GDP[,-1],2,function(x){diff(log(x))}))*100
  RT_GDP <- as.matrix(RT_GDP[RT_GDP$date %in% c("Q2-2005","Q3-2005","Q4-2005",paste0("Q",1:4,"-",rep(2006:2019,each=4)),"Q1-2020","Q2-2020"),-1])
  
  RT_IP <- read_excel("raw-data/IT_RT_IP.xlsx", skip = 6)[-1,-2]
  RT_IP[,-1] <- apply(RT_IP[,-1],2,as.numeric)
  colnames(RT_IP) <- c("date",as.character(seq(as.Date("2000-03-01"),as.Date("2020-09-01"),by="month")))
  RT_IP <- RT_IP[-nrow(RT_IP),]
  RT_IP <- RT_IP[,c("date",as.character(seq(as.Date("2011-01-01"),as.Date("2020-06-01"),by="month")))]
  RT_IP[,-1] <- rbind(NA,apply(RT_IP[,-1],2,function(x){diff(log(x))}))*100
  IP.rt <- as.matrix(RT_IP[RT_IP$date %in% as.character(format(seq(as.Date("2005-04-01"),as.Date("2020-06-01"),by="month"),format="%b-%Y")),-1])
  
  GDP.rt <- kronecker(RT_GDP,matrix(1,nrow=3,ncol=1))
  GDP.rt[-seq(3,nrow(GDP.rt),by=3),] <- NA
}
