library(broom)
library(plyr)
library(dplyr)

#input part
PerM <- 15
setwd("D:/Bedrest_biorhythm/actiheart原始文件")
#print("Please input data frequency, one in ? seconds.")
#PerM <- readline()
#print("Please input data directory请输入文件夹名称")
#datadir <- readline("")
#setwd(datadir)
#file.choose()

#get the number and name of file
data_name <- list.files()
#dir = paste("./dataoutput/",datadir,sep="")
#setwd(dir)
print("File ready to process:")
cat(data_name)
datalength = length(data_name)                                 


#day and light data
statisfunction <- function(datapr){
        datap <- as.numeric(as.character(datapr))
        BPMmean <- mean(datap, na.rm = T)
        BPMsd <- sd(datap, na.rm = T)
        dataprlength = 0
        for (i in 1:length(datap))
        {
                if (is.na(datap[i])== TRUE) dataprlength = dataprlength + 1
        }
        BPMsem <- BPMsd/sqrt(dataprlength)
        BMPlist <- data.frame(BPMmean = BPMmean, BPMsd = BPMsd
                              ,stringsAsFactors=FALSE)
        return(BMPlist)
}        

#Standardization of the output
BPMoutputall <- data.frame(BPMmax = "48H", BPMmin = "48H",BPMmean = "48H",
                           BPMsd = "48H", stringsAsFactors=FALSE)
BPMlightdataall <- data.frame(BPMmean = "light", BPMsd = "light"
                              ,stringsAsFactors=FALSE)
BPMdarkdataall <- data.frame(BPMmean = "dark", BPMsd = "dark"
                             ,stringsAsFactors=FALSE)
BMPperhourlistall <- data.frame(Number = datalength)

#read all file in the dictionary and process
for (j in 1 : datalength) {
        
        #read file and process
        acolnames <- list("a1")
        for(i in 2:48){
                acolnames<- cbind(acolnames,paste("a", i, sep=""))
        }
        CDPcolnames <- c("Time","Activity","BPM","IBI's",acolnames)
        CDPdata <- read.csv(data_name[j],col.names = CDPcolnames, header = F)
        CDPdataA <- CDPdata[1:13,]
        CDPdataB <- CDPdata[16:nrow(CDPdata), ]
        names(CDPdataB)[1:4]<-c("Time","Activity","BPM","IBI's")
        CDPdataB$BPM[CDPdataB$BPM == "0"] <- ""
        CDPdataB$BPM[as.numeric(as.character(CDPdataB$BPM)) > 166] <- ""
        CDPdataB$BPM[as.numeric(as.character(CDPdataB$BPM)) < 43] <- ""
        #data12[strptime("06:50:00","%H:%M:%S")>data12
        #       &data12>strptime("06:30:00","%H:%M:%S")]
        uniformtime <- strptime(CDPdataB$Time,"%H:%M:%S")
        CDPdataBuni <- cbind(CDPdataB,uniformtime)
        BPMlight <- subset(CDPdataBuni, uniformtime>=strptime("06:30:00","%H:%M:%S")
                           &uniformtime<strptime("22:30:00","%H:%M:%S"))
        BPMdark <- subset(CDPdataBuni, uniformtime<strptime("06:30:00","%H:%M:%S")
                          |uniformtime>=strptime("22:30:00","%H:%M:%S"))
        BPMlightdata <- statisfunction(BPMlight$BPM)
        #        BMPlightdata <- cbind("lightdata",BMPlightdata[,-1])
        BPMdarkdata <- statisfunction(BPMdark$BPM)
        #        BPMdarkdata <- cbind("darkdata",BPMdarkdata[,-1])
        
        #Chronos-fit data output
        CDPdataC <- cbind(CDPdataB$Time,CDPdataB$BPM)
        for (l in 1: length(CDPdataB$Time)){
                CDPdataCtime<- unclass(as.POSIXlt(uniformtime[l]))
                CDPdataC[l,1] <- CDPdataCtime$hour + CDPdataCtime$min/60 + CDPdataCtime$sec/3600
        }
        for (g in 1: length(CDPdataB$Time)) {
                if (CDPdataC[g,2] < 3){
                        CDPdataC[g,2] <- 0
                }
        }
        write.table(CDPdataC,
                  file = paste("D:/Bedrest_biorhythm/分析结果/",
                               CDPdataA[2,2], "Chronos.txt", sep=""),
                  row.names = F, col.names = F)
        
        #plot
        hourlist<- as.numeric(strftime(CDPdataBuni$uniformtime[1], format = "%H"))
        nhourlist <- hourlist
        oriminutes <- as.numeric(strftime(CDPdataBuni$uniformtime[1], format = "%M"))
        oriseconds <- as.numeric(strftime(CDPdataBuni$uniformtime[1], format = "%S"))
        oriinteger <- ((60-oriminutes)*60-oriseconds)/PerM-1
        BMPmeanperhour <- mean(as.numeric
                               (as.character(CDPdataBuni$BPM[1:oriinteger])),
                               na.rm = T)
        #calculate the BPM per hour
        for (i in 1:((length(CDPdataBuni$uniformtime)*PerM/3600)-1)) {
                hourlist <- c(hourlist, nhourlist + i)
                NBMPmeanperhour <- mean(as.numeric(as.character(CDPdataBuni$BPM))
                                        [((i-1)*3600/PerM+oriinteger+1):(i*3600/PerM+oriinteger+1)],na.rm = T)
                BMPmeanperhour <- c(BMPmeanperhour, NBMPmeanperhour)
        }
        
        #BMP data processing
        CDPBPM <- as.numeric(as.character(CDPdataB$BPM))
        BPMmax <- max(CDPBPM, na.rm = T)
        BPMmin <- min(CDPBPM, na.rm = T)
        BPMmean <- mean(CDPBPM, na.rm = T)
        BPMsd <- sd(CDPBPM, na.rm = T)
        BPMmedian <- median(CDPBPM, na.rm = T)
        BPMamplitude <- (BPMmax - BPMmin)/2
        BPMnlsdata <- data.frame(x = hourlist, y = BMPmeanperhour)
        BPMnls <- nls(y ~ the1 * sin(((x*pi)/12)+the3)+the2, start = list(the1 = 10, the2=50, the3=5),
                      data = BPMnlsdata, trace = T )
        BPMtmax = as.numeric(tidy(BPMnls)[2,2]) + abs(as.numeric(tidy(BPMnls)[1,2]))
        BPMtmin = as.numeric(tidy(BPMnls)[2,2]) - abs(as.numeric(tidy(BPMnls)[1,2]))
        BPMoutput <- data.frame(BPMmax = BPMmax, BPMmin = BPMmin,BPMmean = BPMmean,
                                BPMsd = BPMsd, stringsAsFactors=FALSE)
        
        
        #calculate the nls, output tmin and tmax
        
        #for (i in 1:length(hourlist)) {
        #        hourlist[i] <- as.character(as.numeric(hourlist[i])%%24)
        #        hourlist[i] <- c(hourlist[i],"day")
        #}
        #output part
        pdf(file = paste(CDPdataA[2,2], "output.pdf", sep=""))
        plot(hourlist, BMPmeanperhour, main = "Heart Rate",
             xlab = "hours", ylab = "BPM",
             col= "black", pch= 19, type = "l")
        BPMnlsy <-as.numeric(tidy(BPMnls)[1,2]) * sin(hourlist * 24/pi) + as.numeric(tidy(BPMnls)[2,2])
        lines(hourlist, BPMnlsy, col = "red")
        dev.off()
        for (g in (1:length(hourlist))){
                hourlist[g] = as.numeric(hourlist[g]) %% 24
        }
        BMPperhourlist <- data.frame(cbind(hourlist,BMPmeanperhour))
        names(BMPperhourlist) <- c(paste(CDPdataA[2,2],"时间",sep=""),
                                   paste(CDPdataA[2,2],"平均BPM",sep=""))
        write.csv(BMPperhourlist
                  ,file =paste("D:/Bedrest_biorhythm/分析结果/",
                               CDPdataA[2,2], "每小时分析.csv", sep=""),row.names = F)
        finaldataoutput <- cbind(c("48H",BPMoutput,"light", BPMlightdata,"dark",BPMdarkdata),
                                 CDPdataA[1:2])
        for (i in 2:13){
                for (k in 1:11) {
                        finaldataoutput[i,k]=""
                }
        }
        write.csv(finaldataoutput,
                  file = paste("D:/Bedrest_biorhythm/分析结果/",
                               CDPdataA[2,2], "线性分析.csv", sep=""),row.names = F)
        
        
        BPMoutputall <- rbind(BPMoutputall,BPMoutput)
        BPMlightdataall <- rbind(BPMlightdataall,BPMlightdata)
        BPMdarkdataall <- rbind(BPMdarkdataall,BPMdarkdata)
        BMPperhourlistall <- rbind.fill(BMPperhourlistall, BMPperhourlist)
        }

nowdate <- Sys.time()
nowdate <- gsub(':', '-', nowdate)

finaldataoutputall <- cbind(BPMoutputall,BPMlightdataall,BPMdarkdataall)
write.csv(finaldataoutputall,
          file = paste("D:/Bedrest_biorhythm/分析结果/",
                nowdate, "线性分析汇总.csv", sep=""),row.names = c("1",data_name))


write.csv(BMPperhourlistall,
          file = paste("D:/Bedrest_biorhythm/分析结果/",
                nowdate, "每小时分析汇总.csv", sep=""),row.names = F)

print("处理完毕")
