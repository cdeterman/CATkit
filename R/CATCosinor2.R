
#' @title Cosinor Analysis 2.0
#' @param TimeCol Column(s) where time is found.  c(1,2) can be used when date is in column 1 and time is in column 2.  
#' Date and time can be in any two columns, i.e., c(6,3). If time is in two columns, the formats should be \%d/\%m/\%y 
#' in the first of the two columns, and \%H:\%M:\%S in the second of the two (Any time format in timeFormat paramater is 
#' ignored in this case.)
#' @param Y Column holding data to be analyzed
#' @param Components Default=1.  Indicates if this is a single or multiple component cosinor analysis, where the number of components is specified (>0).  
#' @param window Data will be convolved with a window:  Hanning, Hamming, Bartlett, Blackman.  Default="noTaper"
#' @param RefDateTime Default=NA  Date used as reference, and subtracted from all data dates, to make the number smaller. 
#' The format is always "\%Y\%m\%d\%H\%M" not implemented---if RefDateTime = NA, use the 1st date of the data as the RefDateTime
#' if RefDateTime = 0, use midnight of the same day as the data starts when using years with leading 0s  (0500 AD) be 
#' sure to include the leading 0, and put quotes around the RefDateTime
#' @param timeFormat RefDateTime and the time in the data file should both use the same date/time format.   
#' Can be "numeric" or "\%Y\%m\%d\%H\%M"??? i.e., "\%Y\%m\%d\%H\%M", or if = "numeric", time column in data file can be 
#' simple numbers (0 - 99999...) if "numeric", data is sorted by time to be sure they are ordered ascending.  
#' First must be smallest , and last largest. Time can also be in two columns (indicate in TimeCol);  
#' timeFormat is ignored when time is in two columns -- the format use is \%d/\%m/\%y in the first of the two 
#' columns, and \%H:\%M:\%S or \%H:\%M in the second of the two
#' @param Units Units (hour, year, week or day) for Interval and Increment arguments, as well as 
#' Period arguments, and any time variable. ToDo:  modify to use units other than hours -- currently only works for hours
#' @param dt When equidistant data, dt indicates the sampling interval.  If dt =0, no periodogram is done.
#' Data is assumed to be equidistant when this is nonzero
#' @export
CATCosinor2 <-
    function (data, TimeCol=1,Y=2, Components=1, window="noTaper", RefDateTime=NA,  
              timeFormat="%Y%m%d%H%M", RangeDateTime=list(Start=NA, End=NA), 
              Units="hours", dt=0, Progressive=list(Interval=0, Increment=0), 
              Period=list(Set=0,Start=0,Increment=1,End=0),
              Debug=FALSE,
              IDcol="fileName", 
              GrpCol = NA) {   
        # Future:  NAs are interpolated;  0s are modified to a very close to 0 value as some functions will not accept 0 data values.
        # Any column selected to be read (Y=column) will be scanned for missing data, and that entire row will be removed.
        #          If you specify two or more columns Y=c(1:3), where ever data is missing in any one column, the entire row will be removed
        #          If this is not acceptable, the program must be run with a single column specified.  Then only rows where data is missing for that column will be removed
        # Progressive
        #   $Interval : length of the time span being analyzed (in Units)  -- multiple spans calculated
        #               If 0, assumes no progession, Interval is set to the full dataset length, and Increment = full data set  
        #  $Increment: number of Days, Wks or Yrs  (uses same unit as set for Interval) to shift forward for each successive Interval analyses
        #               If 0, assumes no progession, Interval is set to the full dataset length, and Increment = full data set 
        # Period
        #       $End :  [only used if $Set=0] Last (and smallest) period to calculate, in units of Days, Wks or Yrs (as set by Units)  EXCLUSIVE
        #               Defaults to 2*dt or 4;  (1 is too small) 0 is invalid  -- default will be used
        #  $Increment : [only used if $Set=0] Increment to starting period, in units of Days, Wks or Yrs (as set by Units)
        #               Defaults to 1;   0 is invalid  -- default will be used
        #     $Start :  [only used if $Set=0] First (and largest) period to calculate, in units of Days, Wks or Yrs (as set by Units);  (Interval/1)
        #               0 is Default: the full time range of the data file is analyzed [in hours?] (RangeDateTime$End-RangeDateTime$Start)= MyData_length; or if progressive, Interval/1;     
        #               Important:  It is normally best if the user sets $Start to a multiple of the largest period of interest  -- Fourier frequencies calculated will be closest to period of interest this way.
        #                           So if one is interested in a weekly period, use a multiple of 1 week as $Start
        #      $Set :  If Set=0, a series of periods are analyzed (spectrum) according to Period$Start, $Increment, $End (or their default value, if not specified)
        #              If not equal to 0, Overrides Period$Start and $Increment, to completely specify an exact set of periods to analyze (as set by Units).  
        #              Can be in the format of a single number, or an array:  c(1,2,3) or c(8,12:24)
        #              When Components=1, each period specified in the vector will be assessed by cosinor independently.
        #              When parameter Components is >1, Period$Set must have a correspondig number of components, which assessed together in a multiple component cosinor.
        #              When 0, only the maximum period, or the maximum period per Interval, from a spectrum is listed on p1 of the graphics file.
        #              Otherwise, all periods are displayed on the graphic   
        # Debug       Turn on when you want to see more diagnostic output for programming debug purposes
        # IDcol       If IDcol="fileName" then use the fileName as the subject ID;  because sometimes the filename is the subject id
        #             Otherwise, this is the data file column number where the subject ID is found, because sometimes the subject id is in the data
        
        #library(season)
        
        #####################################################################################################
        #                                                                                                   #
        #        Set up some initial parameter values, error messages, and key variables                    #
        #                                                                                                   #
        #####################################################################################################
        tz="GMT"
        opar = par()
        # read data;  if no data file name given, prompt the user for the filename
        # if (length(fileName)==0){
        #     fileName <- file.choose()
        # }
        is.wholenumber <-  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
        
        LCM<-Period$Set[1]
        # calculate greatest common divisor and least common multiple; 
        # used later to calculate time in one degree
        if (all(is.wholenumber(Period$Set))){
            tempPeriods<-Period$Set   # no need for a fraction modifier
            fractionModifier<-1
        } else {
            fractionModifier<-100
            tempPeriods<-round(Period$Set, digits=2)*fractionModifier
        }
        minPeriod<-min(Period$Set)  # get the shortest period for use in calculating plotting
        #browser()
        if (Period$Set[1]>0){
            LCM<-max(Period$Set)
            if (Components>1){ 
                "%gcd%" <- function(u, v) {ifelse(u %% v != 0, v %gcd% (u%%v), v)}
                "%lcm%" <- function(u, v) { abs(u*v)/(u %gcd% v)}
                LCM<-1
                
                for (i in 1:Components){
                    LCM<-as.integer(tempPeriods[i]) %lcm% LCM
                }
                LCM<-LCM/fractionModifier
            }
        } else {LCM<-24}      #  default to 24 if spectrum given
        print(LCM)
        
        #  if not yet set, set defaults;  This should be added as an undocumented parameter to override which Graphs print by default
        # if lss1 ..  else if non-progressive2 ... else if gliding spectrum progressive3 ... else if 1-comp progressive4 or multi-comp progressive5
        #  raw data, model, MESOR, Amplitude, Phase, N (# pts), heatmap
        
        
        
        # rename parameters
        StartDate<-RangeDateTime$Start
        EndDate<-RangeDateTime$End
        Interval<-Progressive$Interval
        Increment<-Progressive$Increment
        FreqInc<-Period$Increment
        #Period$Start<-Period$Start
        #Period$End<-Period$End
        oneCycle<-sort(Period$Set,decreasing=TRUE)     # can be a vector of one or more
        Ys_end<-length(Y)
        
        if(!is.numeric(TimeCol)){
            if(TimeCol %in% colnames(data)){
                TimeCol <- which(colnames(data) == TimeCol)    
            }else{
                stop("TimeCol not found in data")
            }
            
        }
        
        Y.orig <- Y
        if(!is.numeric(Y)){
            if(all(Y %in% colnames(data))){
                Y <- which(colnames(data) %in% Y)    
            }else{
                stop("Some 'Y' columns not found in data")
            }
            
        }
        
        paramMsg<-paste("\n  TimeCol=",TimeCol,",  Y=",Y, 
                        "\n --  Periods=",Period["Set"],", Units=",Units,
                        ",  Interval=",format(Interval,nsmall=3), 
                        ",  Increment=",format(Increment,nsmall=3), 
                        "\nPeriod$Start=",format(Period$Start,nsmall=3),
                        ",  FreqInc=",format(FreqInc,nsmall=3), 
                        ",  Period$End=",format(Period$End,nsmall=3),
                        "\nRefDateTime=",RefDateTime, 
                        ", StartDate=",format(StartDate,nsmall=3),
                        ", EndDate=",format(EndDate,nsmall=3),
                        "\nPercent of missing (blank) sample values: unknown")
        
        #  validate and error if too large -- should it be Components*2?  no Feb 26
        CosMatrixDim<-Components*2+1         
        #if (oneCycle[1]>0){
        if (Components>1 && length(oneCycle)!=Components){
            errMsg<-paste("Error:  When parameter Components is >1, Period$Set must have a corresponding number of components.  \nComponents=",Components,"; Period$Set=",Period$Set)
            #closeOutput(file=paste("ERR",fileName),output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
            stop(errMsg)
            # Period$Set will contain multiple possible components to check separately, or as a multi-component cosinor, depending on param Components
            # if Components=1, Period$Set can have any number of components.  If Components>1, length(oneCyce) must equal Components
        }else if (Components==1 && length(oneCycle)!=Components){
            errMsg<-paste("Error:  If Components equals 1, Period$Set must also be only one component.  \nComponents=",Components,"; Period$Set=",Period$Set)
            #closeOutput(file=paste("ERR",fileName),output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
            stop(errMsg)
        }
        #} 
        errMsg<-""            #   this error message gets printed at the beginning of the rtf output file
        #  symbols currently in use for Err values, printed in the matrix
        errKey<-c("*","**","+","++","@","*@","@@@","@@","*@@","*@@@")   
        # messages that correspond to the errKey symbols
        keyMsgs<- vector(mode="character",length = 10)
        keyMsgs[1]<-paste("*  : Error1:  The matrix cannot be inverted for these points and period. This may happen if there are not enough data points for the period selected. Check sampling interval: there should be at least 2p+1 data points per period.\n")
        keyMsgs[2]<-paste("** : Error2:  The matrix cannot be inverted for these points and period. This may happen if there are not enough data points for the period selected. Check sampling interval: there should be at least 2p+1 data points per period.\n")
        keyMsgs[3]<-"+  : Error:  This model does not allow calculation of PR for the individual components. PR is negative, which is invalid.\n"
        keyMsgs[4]<-"++ : Error:  Coefs cannot be calculated\n"        #  this one is likely not needed
        keyMsgs[5]<-"@  : Error:  Requested analysis violates Nyquist.  Time span should be > 2C, where C=#cycles.  Results may be invalid.\n"
        keyMsgs[6]<-"*@ : Warning:  RSS=0.  Not enough data to calculate S.E.  \n"
        #  if there is a period of 8 hours being assessed, there must be 4 data points in 8 hours.
        keyMsgs[7]<-"@@@: Error:  This interval must be at least as long as the trial period (Period$Set).  No analysis performed when Interval<90% of target period.\n"
        keyMsgs[8]<-"@@: Warning:  The interval must be at least as long as the trial period (Period$Set).  Results may be unreliable.\n"
        keyMsgs[9]<-"*@@: ERROR:  RSS/df2 will be infinite because df2<0, so this case was skipped.  (Number of data points <= number of components*2+1)\n"
        keyMsgs[10]<-"*@@@: Warning:  RSS/df2 will be infinite because df2=0, so s.e.s, F and P calculations are skipped.  (Number of data points <= number of components*2+1)\n"
        
        # extract filename from path for display on output 
        # if (.Platform$file.sep=="/"){
        #     m=regexec("[^/]*$",fileName)            #  ./(.$)    pc=[^\\]*$   mac=[^/]*$
        # } else {
        #     m=regexec("[^\\]*$",fileName)
        # }
        # fileName1<-regmatches(fileName,m) 
        # fileLen<-nchar(fileName1)
        # fileName6<-substring(fileName1,1,fileLen-4)
        # send output to file  -- build name of file with path to data f
        BaseTime<-Sys.time()        
        thisTime <- format(BaseTime, "--%d%b%y--%H-%M-%OS")
        # fileName2<-paste(fileName,window,thisTime,functionName,"Cos.txt",sep="")
        # fileName3<-paste(fileName,window,thisTime,functionName,"Cos.rtf",sep="")
        # if (Output$Txt){
        #     sink(fileName2, append=FALSE, split=TRUE)
        # }
        
        
        #####################################################################################################
        #                                                                                                   #
        #        Read data file;  assume tab delimited;  if that fails, read as CSV;  get header;           #
        #                omit missing data and count missing data;                                          #
        #                                                                                                   #
        #####################################################################################################
        
        
        
        # if (is.character(MyDataRead[,-TimeCol[1]])){    # are any columns character, excluding date col??
        #   str(MyDataRead)
        #   message<-paste("ERROR:  Data is being read as character data.  This will cause all data to be read as character.")
        #   errMsg<-paste(errMsg,message)
        #   closeOutput(file=fileName3,output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
        #   stop(message) 
        # }
        # if (IDcol=="fileName"){    #set once here
        #     SubjectID<-fileName6
        # } else {
            # SubjectID<-data[1,IDcol]     #  set for each Loop, within i Loop
        # }
        
        # output <- vector(mode = "list", length = length(Y))
        grp_inc <- ifelse(is.na(GrpCol), 0, length(GrpCol))
        
        grp_ids <- colnames(data)[which(colnames(data) %in% GrpCol)]
        
        # print(grp_inc)
        # print(grp_ids)
        df <- setNames(data.frame(matrix(nrow = 1, ncol = 10 + grp_inc)), 
                 c(grp_ids,
                   c("Interval",
                     "Period",
                     "P",
                     "PR",
                     "Mesor",
                     "M.s.e",
                     "Amp",
                     "A.s.e",
                     "Phi",
                     "P.s.e")))
        df[,grp_ids] <- as.character(df[,grp_ids])
        df[,(ncol(df) - 9):ncol(df)] <- as.numeric(df[,(ncol(df) - 9):ncol(df)])
        
        output <- lapply(Y, function(x) df)
        names(output) <- Y.orig
        
        # print(str(output))
        # stop("stop now")
        
        # make sure drop unused levels
        data <- droplevels(data)
        
        # get levels to iterate over
        subjectIDs <- levels(data[,IDcol])
        
        for(sub in seq_along(subjectIDs)){
            
            if(sub > 1){
                StartDate <- NA
                EndDate <- NA
            }
            
            # print(head(data))
            
            # print(sub)
            
            # print(which(data[,IDcol] == subjectIDs[sub]))
            
            mydata <- data[which(data[,IDcol] == subjectIDs[sub]),]
            mydata <- droplevels(mydata)
            
            print(head(mydata))
            
            keepers<-c()
            for (y in c(TimeCol,Y)){
                # print(y)
                # print(head(mydata[,y]))
                #MyDataRead<-subset(data, all(!is.na(data[,c(TimeCol,Y)])) )     # omit data row if any Time or Y column is NA
                #  MyDataRead<-na.omit(data[,c(TimeCol,Y)])   changes column numbers, so problematic
                is.na(mydata[,y])<-which(mydata[,y]=="")    #   replace spaces with NA
                if (all(is.na(mydata[,y]))){
                    message<-paste("ERROR:  You have selected a column that has no data in it:  column",y)
                    errMsg<-paste(errMsg,message)
                    # closeOutput(file=fileName3,output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
                    stop(message) 
                    
                }
                keepersY<-na.omit(mydata[,y])            #  eliminate NA
                keepers<-union(na.action(keepersY),keepers)     #  Add NA rows to previous NA rows
                #   use inverse of omitted NA indexes to get good values
            }
            if (length(keepers)>0){
                MyDataRead<-mydata[-keepers,]
                message<-"Note:  Missing data was found.  All Y columns values in that row were omitted."
                print(message)
                errMsg=paste(errMsg,message)
            } else {
                MyDataRead<-mydata
            }
            missingData<-(length(data[,TimeCol])-length(MyDataRead[,TimeCol]))/length(data[,TimeCol])    #  rowsData<-length(data[,1])
            
            #MyDataRead<-subset(data, !is.na(data[,Y]))       # omit data row if Y column is NA
            #MyDataRead = na.omit(head(MyDataRead1,n=-2))   #  will drop the very last row in case it had no data
            #####################################################################################################
            #                                                                                                   #
            #        Time column: convert time format to numeric if needed, or to number of Units of time       #
            #                May be in one or two columns;  2 columns must be appended before conversion        #
            #                                                                                                   #
            #####################################################################################################
            
            timeFormatOrig <- timeFormat
            if (timeFormat=="numeric"){
                BaseDate<-format(Sys.Date(), "%Y%m%d%H%M")   # as.POSIXct(strptime(Sys.Date(), "%Y%m%d%H%M",tz=tz))     # midnight of first day of data
                #BaseDate2 <- paste(substr(BaseDate,1,8),"0000",sep="")
                # apply timezone hack to get correct time
                origin <- as.POSIXct('1970-01-01 00:00:00',tz=tz) 
                offset <- as.numeric(origin)
                #browser()
                newOrder<-MyDataRead[order(MyDataRead[,TimeCol], na.last = TRUE, decreasing = FALSE),]
                MyDataRead<-newOrder
                #browser()
                MyDataRead[,TimeCol]<-format(as.POSIXct(as.numeric(as.POSIXct(strptime(BaseDate, "%Y%m%d%H%M",tz=tz)))+(as.numeric(MyDataRead[,TimeCol])*3600),origin=origin,tz),"%Y%m%d%H%M")
                timeFormat<-"%Y%m%d%H%M"
                timeFormatOrig<-"numeric"
                
            } else {   #convert time from 2 columns  (only one format allowed if time is in 2 columns:   17/10/07 \t 19:00:20)
                if (grepl(".",MyDataRead[1,TimeCol[1]],fixed=TRUE)){    # timeFormat not == 'numeric' but time is numeric
                    str(MyDataRead)
                    #browser()
                    message<-paste("ERROR:  The time format does not match timeFormat setting.  Time is in a numeric/decimal format, but timeFormat is not set to 'numeric'.")
                    errMsg<-paste(errMsg,message)
                    # closeOutput(file=fileName3,output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
                    stop(message) 
                }
                if (length(TimeCol)==2){
                    MyDataRead[,TimeCol[1]] <- as.POSIXct(strptime(x=paste(MyDataRead[,TimeCol[1]],MyDataRead[,TimeCol[2]]), format=timeFormat, tz=tz))      #  "%d/%m/%y %H:%M:%S"
                    MyDataRead[,TimeCol[2]] <- format(MyDataRead[,TimeCol[1]],"%Y%m%d%H%M")
                    names(MyDataRead)[TimeCol[2]] <- "Datetime"
                    TimeCol<-TimeCol[2]
                } else {   #convert time from 1 columns
                    if (length(TimeCol)==1){
                        MyDataRead[,TimeCol[1]] <- as.POSIXct(strptime(x=MyDataRead[,TimeCol[1]], format=timeFormat, tz=tz))
                        MyDataRead[,TimeCol[1]] <- format(MyDataRead[,TimeCol[1]],"%Y%m%d%H%M")
                        names(MyDataRead)[TimeCol[1]] <- "Datetime"
                    }}}
            #####################################################################################################
            #                                                                                                   #
            #        RangeDateTime$Start and $End: calculate total time length, StartDate and EndDate           #
            #                May be in one or two columns;  2 columns must be appended before conversion        #
            #                    NA- use the 1st time in the data set                                           #
            #                    0- use Midnight at the start/end of the data                                   #
            #                Calculate start and ending indexes in the data                                     #
            #                print actual hours to be used for analysis, and 1st and last data point in file    #
            #                                                                                                   #
            ##################################################################################################### 
            #browser()
            # str(MyDataRead)                                         #  removed +3600 from *tz=tz))+3600)*       6/14
            AnalyzeLength<-(as.POSIXct(strptime(tail(MyDataRead[,TimeCol],1), "%Y%m%d%H%M",tz=tz))) - as.POSIXct(strptime(MyDataRead[1,TimeCol], "%Y%m%d%H%M",tz=tz)) 
            MyData_HRs<-(AnalyzeLength)*24  # hours in length of actual data
            if (is.na(StartDate)){
                StartDate <- format(as.POSIXct(strptime(MyDataRead[1,TimeCol], "%Y%m%d%H%M",tz=tz)),"%Y%m%d%H%M")
            } else {if (StartDate == 0){
                StartDate <-paste(substr(MyDataRead[1,TimeCol],1,8),"0000",sep="")                #  Use StartTime as start of analysis date-time
            } else StartDate<-format(as.POSIXct(strptime(StartDate, timeFormat,tz=tz)),"%Y%m%d%H%M")     #  needed for when StartDate is specified in a different format--test!!
            } #as.POSIXct(strptime(MyData[1,TimeCol], "%Y%m%d",tz=tz))
            #MyData_length<-length(MyDataRead[,1])
            StartIndex<-which(MyDataRead[,TimeCol]>=StartDate)   #  where is the desired startTime?  (if start point is later than first data point)
            MyData_length<-length(StartIndex)  
            
            print("check ranges")
            print(StartIndex)
            print(EndDate)
            
            if (is.na(EndDate)){                                                     #  removed tz=tz))+3600      6/14
                EndDate <- format(as.POSIXct(strptime(tail(MyDataRead[StartIndex,TimeCol],1), "%Y%m%d%H%M",tz=tz)),"%Y%m%d%H%M") # should these be pasted or posix?  todo
            } else {if (EndDate == 0){
                #EndDate <-paste(substr(as.numeric(MyDataRead[MyData_length,TimeCol])+as.numeric(10000),1,8),"0000",sep="")                #  Use StartTime as start of analysis date-time
                EndDate<-paste(substr(format(as.POSIXct(strptime(MyDataRead[MyData_length,TimeCol], "%Y%m%d%H%M",tz=tz))+86400,"%Y%m%d%H%M") ,1,8),"0000",sep="")
            } else EndDate<-format(as.POSIXct(strptime(EndDate, timeFormat,tz=tz)),"%Y%m%d%H%M")     #  needed for when EndDate is specified in a different format--test!!
            }
            
            # print(head(MyDataRead))
            
            EndIndex<-which(MyDataRead[,TimeCol]<=EndDate)   #  where is the desired endTime? 
            #cat("Check that your data file and End Date date are correct.  Data file begins:",MyDataRead[1,TimeCol],"; End Date is: ",EndDate,"\n")
            
            print("check after")
            print(StartIndex[1])
            print(EndIndex)
            
            MyData<-MyDataRead[StartIndex[1]:tail(EndIndex,n=1),]     #  grab all the data from array after StartTime, before EndTime
            
            #  Analysis will proceed from",StartDate," to ",EndDate,"\n")
            cat("There are ",MyData_HRs," actual hours of data in this file, ",MyDataRead[1,TimeCol]," to ",MyDataRead[tail(EndIndex,n=1),TimeCol],", and ",MyData_length,"data points.\n  %",missingData*100,"of data points are missing.\n")      
            
            print(MyDataRead[1,])
            print(tail(MyDataRead,1))
            #####################################################################################################
            #                                                                                                   #
            #        RefDateTime: get correct reference time;  subtract from data, start tim, and end time      #
            #                convert time to a count in Units, as set by Units parameter                        #
            #                    calculate MyData$time.n, MyData$time.hour and  MyData_length (used in analysis)#                                          #
            #                    print total hours used in analysis and start/end times                         #
            #                                                                                                   #
            ##################################################################################################### 
            
            # convert time
            MyData$time=as.POSIXct(strptime(MyData[,TimeCol], "%Y%m%d%H%M", tz=tz))     #   tz="GMT"
            # convert time to a number (time in seconds, converted into hours)
            # should probably be done with a hard coded time zone, so CDT and DST don't cause problems
            
            if (is.na(RefDateTime)){
                RefTime <-as.POSIXct(strptime(MyData[1,TimeCol], "%Y%m%d%H%M",tz=tz))                #  11/16 chg from timeFormat to "%Y%m%d%H%M": when RefDatTime=NA  Use StartTime as reference date-time
                RefTimeString <- MyData[1,TimeCol]
            } else if (RefDateTime==0){
                if (timeFormatOrig=="numeric"){
                    RefTime<-as.POSIXct(strptime(MyData[1,TimeCol], "%Y%m%d",tz=tz))     # midnight of first day of data   why would it be "%Y-%m-%d"??  changed
                    RefTimeString <- paste(substr(MyData[1,TimeCol],1,10),"0000",sep="")
                } else {
                    RefTime<-as.POSIXct(strptime(MyData[1,TimeCol], "%Y%m%d",tz=tz))     # midnight of first day of data
                    RefTimeString <- paste(substr(MyData[1,TimeCol],1,8),"0000",sep="")
                }
            } else if (RefDateTime<=24 && timeFormatOrig=="numeric"){
                #RefTime<-as.POSIXct(strptime(MyData[1,TimeCol], "%Y%m%d",tz=tz))     # midnight of first day of data   why would it be "%Y-%m-%d"??  changed
                # RefDateTimePart<-round(RefDateTime * 3600, units="min")
                # RefTime <- as.POSIXct(strptime(x=paste(MyData[1,TimeCol],RefDateTimePart), format=timeFormat, tz=tz))      #  "%d/%m/%y %H:%M:%S"
                # RefTime <- format(RefTime,"%Y%m%d%H%M")
                # 
                # BaseDate<-format(Sys.Date(), "%Y%m%d%H%M")         # midnight of first day of data
                # # apply timezone hack to get correct time
                # origin <- as.POSIXct('1970-01-01 00:00:00',tz=tz) 
                # offset <- as.numeric(origin)
                RefTime<-as.POSIXct(as.numeric(as.POSIXct(strptime(BaseDate, "%Y%m%d%H%M",tz=tz)))+(as.numeric(RefDateTime)*3600),origin=origin,tz)
                RefTimeString <- format(RefTime,"%Y%m%d%H%M")    #  ,"%Y%m%d%H%M")
                
                
                
            } else {  #  RefDateTime is a specific date and time
                RefTime<-as.POSIXct(strptime(RefDateTime, timeFormat,tz=tz))
                RefTimeString <- RefDateTime
            }
            
            StartTime<-(as.numeric(as.POSIXct(strptime(StartDate, "%Y%m%d%H%M",tz=tz)))-as.numeric(RefTime))/3600      # convert to hours
            if (is.na(StartTime) ){
                message<-paste("ERROR:  StartTime or StartDate or Date format is invalid: StartTime",StartTime," StartDate: ",StartDate," RefTime:",RefTime)
                errMsg<-paste(errMsg,message)
                # closeOutput(file=fileName3,output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
                stop(message) 
            } 
            EndTime<-(as.numeric(as.POSIXct(strptime(EndDate, "%Y%m%d%H%M",tz=tz)))-as.numeric(RefTime))/3600   # convert to hours
            
            if (is.na(EndTime) ){
                message<-paste("ERROR:  EndDate or Date format is invalid: ",EndTime)
                errMsg<-paste(errMsg,message)
                # closeOutput(file=fileName3,output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
                stop(message) 
            } 
            
            MyData$time.n=as.numeric(MyData$time) - as.numeric(RefTime)
            MyData$time.hour<-MyData$time.n/3600
            MyData_length<-length(MyData[,1])
            drawPgram<-FALSE
            if (dt==0){
                drawPgram<-TRUE
                sumTab<-diff(MyData$time.hour)          #  acts only on integer values, use .n
                #  dt<-sumTab[median(which(sumTab == max(sumTab, na.rm=TRUE)))]
                dt<-median(sumTab, na.rm=TRUE)
            }
            #   if (timeFormatOrig=="numeric"){
            #     EndTime<-max(MyData[StartIndex[1]:tail(EndIndex,n=1),TimeCol])
            #     StartTime<-min(MyData[StartIndex[1]:tail(EndIndex,n=1),TimeCol])
            #   } 
            MyData_hours<-as.numeric(EndTime-StartTime) + dt         #   get number of hours to analyze
            
            cat("The program is using ",MyData_hours," hours for the analysis, starting at hour",StartTime,"and ending at hour",EndTime," relative to RefTime:",RefTimeString,"\n")
            # calculate number of hours in the Interval and Increment (Units*Interval, Units*Increment)
            # this is used for inner loop, progressive analysis of intervals over the full data set
            if (Interval==0 || Increment==0){     #  || oneCycle[1]>0){     #  No progression -- analyze full dataset
                Interval <-MyData_hours
                Increment<-MyData_hours
            } else if (Units=='hours'){
                Interval<-Interval
                Increment<-Increment
            } else if (Units=='days'){
                Interval<-Interval*24
                Increment<-Increment*24
            } else if (Units=='years'){
                Interval<-Interval*24*365
                Increment<-Increment*24*365
            } else if (Units=='weeks'){
                Interval<-24*7*Interval      
                Increment<-24*7*Increment
            } else {
                message<-paste("ERROR:  Units is improperly given: ",Units," Valid values:  hours, days, weeks, years")
                errMsg<-paste(errMsg,message)
                # closeOutput(file=fileName3,output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
                stop(message) 
            }
            if (Period$Start==0){       
                Period$Start <-Interval
            } else if (Period$Start<0){
                message<-"ERROR:  The parameter Period$Start cannot be negative.\n"
                errMsg<-paste(errMsg,message)
                # closeOutput(file=fileName3,output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
                stop(message)
            } else if (Period$Start>(1.25*MyData_hours)){        #  validation of Period$Start   
                message<-paste("ERROR:  The parameter Period$Start",Period$Start, "is much longer than the observation period",MyData_hours,".  Unreliable results.\n")
                errMsg<-paste(errMsg,message)
                print(message)
            } else if (Period$Start>MyData_hours){        #  validation of Period$Start   
                message<-paste("Warning:  The parameter Period$Start",Period$Start, "is slightly longer than the observation period",MyData_hours,".  Results may be unreliable.\n")
                errMsg<-paste(errMsg,message)
                print(message)
            } # end else Period$Start==0
            if (FreqInc==0){
                FreqInc<-1
                paramMsg<-paste("\n  TimeCol=",TimeCol,",  Y=",Y,
                                "\n --  Periods=",Period["Set"],
                                ", Units=",Units, ",  Interval=",format(Interval,nsmall=3), 
                                ",  Increment=",format(Increment,nsmall=3), 
                                "\nPeriod$Start=",format(Period$Start,nsmall=3),
                                ",  FreqInc=",format(FreqInc,nsmall=3),
                                ",  Period$End=",format(Period$End,nsmall=3), 
                                "\nRefDateTime=",RefDateTime, ", StartTime=",format(StartTime,nsmall=3),
                                ", EndTime=",format(EndTime,nsmall=3),
                                "\nPercent of missing (blank) sample values: %",missingData*100,
                                "\n")
            } else if (FreqInc<0){
                message<-"ERROR:  The parameter Period$Increment cannot be negative.\n"
                errMsg<-paste(errMsg,message)
                # closeOutput(file=fileName3,output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
                stop(message)
            } 
            #browser()
            if (Period$End==0){
                if (dt==0){
                    Period$End<-4
                } else if (Period$End<0){
                    message<-"ERROR:  The parameter Period$End cannot be negative.\n"
                    errMsg<-paste(errMsg,message)
                    # closeOutput(file=fileName3,output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
                    stop(message)
                } else if (Period$End>MyData_hours){
                    message<-"ERROR:  The parameter Period$End is larger than the number of hours.  ???\n"
                    errMsg<-paste(errMsg,message)
                    print(message)
                } else Period$End<-3*dt
                if (Period$End>Period$Start){
                    Period$End<-Period$Start/4
                    FreqInc<-.1
                }
            }    #    end Period$End==0
            cat(Period$Start," Period$Start ",Period$End," Period$End ",Period$Start,"Period$Start",Interval,"Interval",Increment,"Increment\n")
            
            paramMsg<-paste("\n  TimeCol=",TimeCol,",  Y=",Y,"\n --  Periods=",Period["Set"],", Units=",Units, ",  Interval=",format(Interval,nsmall=3), ",  Increment=",format(Increment,nsmall=3), "\nPeriod$Start=",format(Period$Start,nsmall=3), ",  FreqInc=",format(FreqInc,nsmall=3), ",  Period$End=",format(Period$End,nsmall=3), "\nRefDateTime=",RefDateTime, ", StartTime=",format(StartTime,nsmall=3),", EndTime=",format(EndTime,nsmall=3),"\nPercent of missing (blank) sample values: %",missingData*100,"\n")
            
            par(mar=c(4,4,1,1)) # set up margins
            
            #  #####    note on indexes and variables:  one set tracks through the ACTUAL times and data points:  #####
            #                EndIdx, StartIdx, thisIdxCnt, MyData_length, TimeIdx
            #  #####    another set tracks through the periods, in hours, being calculated:       ######
            #                MyData_hours; StartSpans, Interval, Increment, Progression_end, StartTime, StartSpanTime, EndTime
            #  #####    These two must, of course, align. See:   TimeIdx<-which(StartSpanTime <= MyData$time.hour) 
            EndIdx<-0
            
            if (Interval > MyData_hours){    #  warn if interval chosen is longer than data
                message<-paste("Warning: Chosen Progressive$Interval (",Interval,") is longer than data file (",MyData_hours,") (MyData_hours).\n")
                print(message)
                errMsg<-paste(errMsg,message)
            } 
            if (Increment> MyData_hours) {
                message<-paste("Warning:  Chosen Progressive$Increment (",Increment,") is longer than data file (",MyData_hours,") (MyData_hours).\n")
                print(message)
                errMsg=paste(errMsg,message)
            }
            # seq(72-84+1,192-50.4,by=8.4)
            if (MyData_hours==Interval){    #StartTime==MyData$time.hour[1] && 
                StartSpans<-seq(from=0,to=MyData_hours-1,by=Increment)        #StartSpans<-seq(1,Interval, by=Increment)
                #} else {StartSpans<-seq(from=(StartTime),to=(StartTime+MyData_hours-Interval),by=Increment)  
            } else {if ((MyData_hours-Interval+Increment)<=0){
                errMsg<-paste(errMsg,"Error:  MyData_hours-Progressive$Interval+Progressive$Increment is < 1:  ",MyData_hours,"-",Interval,"+",Increment,".  Using ",MyData_hours,"as Interval.\n")
                # closeOutput(file=fileName3,output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
                Interval<-MyData_hours     # ????? is this ever run?
            }
                StartSpans<-seq(from=0,to=(MyData_hours-Interval+Increment),by=Increment)     #  one too many (SS), perfect for all (NA/NA)
                #  to=(MyData_hours-1-Interval):  added + Increment because 3varsProg was doing 1 too few spans
                #  to=(MyData_hours-1-Interval+Increment):  #best  2 too few when inc=.5 interval=24  serial section
                #  to=(MyData_hours+1-Increment)         7 too many when inc=.5 interval=24 serial section
                #  to=(MyData_hours-Interval+Increment)  one too many  (should be LESS THAN last span value) when inc=.5 interval=24  serial section
            }   #StartSpans<-seq(1,Interval, by=Increment)
            
            Progression_end<-length(StartSpans)
            #  print Y title, and points Interval
            sumN <- matrix(data=NA,nrow=Ys_end,ncol=Progression_end+1)
            sumLow <- matrix(data=NA,nrow=Ys_end,ncol=Progression_end+1)
            sumHi <- matrix(data=NA,nrow=Ys_end,ncol=Progression_end+1)
            sumMean <- matrix(data=NA,nrow=Ys_end,ncol=Progression_end+1)
            sumMedian <- matrix(data=NA,nrow=Ys_end,ncol=Progression_end+1)
            sumMode <- matrix(data=NA,nrow=Ys_end,ncol=Progression_end+1)
            sumSD <- matrix(data=NA,nrow=Ys_end,ncol=Progression_end+1)
            sumT <- matrix(data=NA,nrow=Ys_end,ncol=Progression_end+1)
            #sumN[Progression_end+1]<-MyData_length
            yNew<-0
            
            for (y in Y){   #  parse each Y columng to get the summary component of each;  goes in final row after all j summaries
                yNew<-yNew+1
                sumN[yNew,Progression_end+1]<-MyData_length
                sumLow[yNew,Progression_end+1]<-min(MyData[,y], na.rm=TRUE)
                sumHi[yNew,Progression_end+1]<-max(MyData[,y], na.rm=TRUE)
                sumMean[yNew,Progression_end+1]<-mean(MyData[,y], na.rm=TRUE)
                sumMedian[yNew,Progression_end+1]<-median(MyData[,y], na.rm=TRUE)   #  could be proper dt if VERY irregular data
                sumSD[yNew,Progression_end+1]<-sqrt(var(MyData[,y],y=NULL))
                sumTab<-tabulate(MyData[,y]*1000)          #  acts only on integer values, ensure integers;   each value is counted in it's cardinal location     
                dTtest<-which(sumTab == max(sumTab, na.rm=TRUE))/1000   #  which data value is most often used? 
                sumT[yNew,Progression_end+1]<-dt
                if (length(dTtest)>1){    # if more than 4 are the same as the max dT, average them to get Mode
                    sumMode[yNew,Progression_end+1] <- mean(dTtest, na.rm=TRUE)
                } else {sumMode[yNew,Progression_end+1]<-dTtest   #  The one used the most is the Mode
                }
                #   sumT      hours+dt=T acknowledges that most times points are binned, so should be more than actual hours in file
            }
            
            Page<-0
            if (oneCycle[1]==0){                    #  if oneCycle != 0, Period$End and FreqInc are not used.
                if (Period$End <= 0 || FreqInc <= 0){
                    message<-paste("ERROR:  Period$End (",Period$End,") and Period$Increment (",FreqInc,") cannot be 0 when Period$Set is 0.\n")
                    errMsg<-paste(errMsg,message)
                    # closeOutput(file=fileName3,output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
                    stop(message)
                }
                if (Period$Start<Period$End){
                    message<-paste("ERROR:  The parameter Period$Start (",Period$Start,") cannot be smaller than Period$End (",Period$End,") when Period$Set is 0.\n")
                    errMsg<-paste(errMsg,message)
                    # closeOutput(file=fileName3,output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
                    stop(message)
                }
                if (FreqInc > Period$Start){    #  is this needed??    6/14
                    message<-paste("ERROR:  Period$Increment (",FreqInc,") cannot be larger than Period$Start (",Period$Start,") (or Progressive$Interval) when Period$Set is 0.\n")
                    errMsg<-paste(errMsg,message)
                    # closeOutput(file=fileName3,output=Output,console=Console,opar=opar,ERR=TRUE,errMsg=errMsg,paramMsg=paramMsg)
                    stop(message)
                } else if (FreqInc>1) {
                    message<-"Warning:  Chosen Increment (>1) is less than optimal.  Increments would optimally be <=1.\n"
                    print(message)
                    #closeOutput(file=fileName3,output=Output,console=Console,opar=opar)
                    errMsg=paste(errMsg,message)
                }
                
                # How many loops?  the integer (truncated) part of 1 + (tau-s/tau-e -1)/delta.
                RowCnt<-floor(1 + ((Period$Start/Period$End)-1)/FreqInc)      #  using ceiling with out 1+ is not equivalent, as non integer result is 1 too small
                RowCntAlt2<-floor(1 + ((1/Period$End)/(1/Period$Start)-1)/FreqInc)    #  for example:  start=52, end=2, inc = 5 should give 6 results;  (52/2)+1 harmonics
                
                if (RowCnt<1){
                    RowCnt<-1                              #  8/11/14  for LS needs to be 1   (was 2)
                }
            } else {if (Components>1){
                RowCnt<-1}
                else {RowCnt <- length(oneCycle)}   #  if oneCycle != 0, reset Period$End and FreqInc to default=1
                print("resetting freq start and inc")
            }
            
            # moved array definitions from HERE
            for (y in 1:Ys_end){
                Ycol<-Y[y]
                #par(mfrow=c(3,1),mar=c(4,4,2,1))   # 6/10/2016
                #  test 6/10/2016  layout(matrix(c(1,2,3), 3, 1, byrow = TRUE), widths=c(1), heights=c(1,3,3))   #   this gets reset in CatWindow, but is needed when j>1       
                paramMsg<-paste("\n  TimeCol=",TimeCol,",  Y=",Y, "\n --  Periods=",Period["Set"],", Units=",Units, ",  Interval=",format(Interval,nsmall=3), ",  Increment=",format(Increment,nsmall=3), "\nPeriod$Start=",format(Period$Start,nsmall=3), ",  FreqInc=",format(FreqInc,nsmall=3), ",  Period$End=",format(Period$End,nsmall=3), "\nRefDateTime=",RefDateTime, ", StartTime=",format(StartTime,nsmall=3),", EndTime=",format(EndTime,nsmall=3),"\nPercent of missing (blank) sample values: %",missingData*100,"\n")
                
                n<-length(MyData[,Ycol])
                if (is.numeric(MyData[n, TimeCol])){
                    startDateP<-as.POSIXct(strptime(x=MyData[1,TimeCol], format="%Y%m%d%H%M", tz=tz))
                    endDateP<-as.POSIXct(strptime(x=MyData[n,TimeCol], format="%Y%m%d%H%M", tz=tz))
                } else {startDateP<-MyData[1,TimeCol]
                endDateP<-MyData[n,TimeCol]}
                
                realLCM<-LCM
                if (minPeriod==0){  #  if this is a LS
                    modelLen<-360
                } else {
                    modelLen<-(LCM/minPeriod)*6       #  model should be long enough that each cycle of the shortest period (minPeriod) has 6 pts (LCM/minPeriod is how many cycles in LCM)
                    if (LCM>MyData_hours){
                        LCM<-MyData_hours
                        modelLen<-(LCM/minPeriod)*24
                    }
                }
                
                if (modelLen<360){      # should be at least 360 points generated in a model
                    modelLen<-360        # how many points to calculate in the model (for valid plotting representation of cycle -- to avoid aliasing)  
                }
                
                #SubjectID<-vector(mode="integer",length = Components)
                M <- matrix(data=NA,nrow=RowCnt, ncol=Progression_end)
                Model_Y <- matrix(data=NA,nrow=RowCnt, ncol=modelLen)
                Model_Y_mag <- matrix(data=NA,nrow=RowCnt, ncol=modelLen)
                #plotModel <- matrix(data=NA,nrow=RowCnt, ncol=modelLen)     # 6/27/2016
                plotModel <- vector(mode="numeric",length = modelLen)
                multi_Model_Y <- matrix(data=NA,nrow=RowCnt, ncol=modelLen)
                multi_Model_mag <- matrix(data=NA,nrow=RowCnt, ncol=modelLen)
                newPR <- array(data=NA,c(RowCnt, Progression_end,Components))
                PHI <- array(data=NA,c(RowCnt, Progression_end,Components))
                PHIr <- array(data=NA,c(RowCnt, Progression_end,Components))
                A <- array(data=NA,c(RowCnt, Progression_end,Components))
                PR <- array(data=NA,c(RowCnt, Progression_end,Components))
                P <- array(data=NA,c(RowCnt, Progression_end,Components))
                testPg <- array(data=NA,c(RowCnt, Progression_end,Components))
                testPb <- array(data=NA,c(RowCnt, Progression_end,Components))
                F <- array(data=NA,c(RowCnt, Progression_end,Components))
                mesor_se <- matrix(data=NA,nrow=RowCnt, ncol=Progression_end)
                phi_se <- array(data=NA,c(RowCnt, Progression_end,Components))
                amp_se <- array(data=NA,c(RowCnt, Progression_end,Components))
                Cycle <- array(data=NA,c(RowCnt, Progression_end,Components))
                Err <- array(data="",c(RowCnt, Progression_end,Components))
                hours <-matrix(data=NA,nrow=RowCnt, ncol=Progression_end)
                nPts<- matrix(data=NA,nrow=RowCnt, ncol=Progression_end)
                sPts<- matrix(data=NA,nrow=RowCnt, ncol=Progression_end)
                time<- matrix(data=NA,nrow=RowCnt, ncol=Progression_end)
                yVar<- matrix(data=NA,nrow=RowCnt, ncol=Progression_end)
                MSS<- vector(mode="integer",length = Components)
                printP<- array(data=NA,c(RowCnt, Progression_end,Components))
                if (Components>1){
                    multi_P<-matrix(data=NA,nrow=RowCnt, ncol=Progression_end)
                    multi_PR<-matrix(data=NA,nrow=RowCnt, ncol=Progression_end)
                    Magnitude<-matrix(data=NA,nrow=RowCnt, ncol=Progression_end)
                }
                Orthophase<-matrix(data=NA,nrow=RowCnt, ncol=Progression_end)
                Bathyphase<-matrix(data=NA,nrow=RowCnt, ncol=Progression_end)
                newData<- matrix(data=NA,nrow=RowCnt, ncol=Progression_end)    # holds filtered or unfiltered data, myData for each Interval
                #----------------------array definitions
                
                Page<-0
                for (j in 1:Progression_end){   #  make a vector with incremental values.  use var(i) for necessary indexes
                    maxAmp<-max2Amp<-max3Amp <- 1  
                    StartSpanTime<-StartTime+StartSpans[j]    # Feb 2014 added back in; 11/4  StartTime+StartSpans[j]-1 # 10/18 added StartTime+ DH;  9/30  MyData$time.hour[1]+StartSpans[j]-1
                    EndSpanTime<-StartSpanTime +Interval          # 6/14 no subtraction; 5/2014 subtract 1 oneSec instead; 9/30  -1***adj to match gc  Add Interval length to the last Interval StartSpans[j]+Interval to get the ending of the Interval
                    # when using 24 hour progressive Interval, if should not go from 8:00 to 8:00 in one Interval, should only go from 8:00 to 7:59:59!  
                    TimeIdx<-which(StartSpanTime <= MyData$time.hour)   #   get indices past the  time of this Interval
                    if (is.na(TimeIdx[1])){
                        StartIdx<-length(MyData$time.hour)     #  use the largest index (last data point)
                    } else StartIdx<-TimeIdx[1]                # index of starting time is FIRST index in array of indices
                    if (EndSpanTime>MyData_hours+StartTime){   # Feb 2014 had to add + StartTime to these 2 lines
                        EndSpanTime<-MyData_hours+StartTime
                    }
                    #if (EndSpan>EndTime){
                    #  print("**************what to do about this???****************")
                    #  EndSpan<-MyData_hours
                    #}
                    if (Progression_end==j){                  #  6/14   last (or only) spans need to be treated differently than all others
                        # needed to round because equal numbers didn't compare as equal due to differing precisions  (116.5333333 <> 116.5333)
                        TimeIdx<-which(round(EndSpanTime,digits=10) >= round(MyData$time.hour,digits=10))   #  should be gte;  get indices UP TO the end of this Interval
                        #  in what case does this need to be gte?   ActivityP1-ZCM-SSend2 serial section might need it to be equal
                    } else TimeIdx<-which(EndSpanTime > MyData$time.hour)    #  6/14
                    EndIdx<-which.max(TimeIdx)               # index of ending time is LAST index in array if indices
                    Interval_hours<-EndSpanTime-StartSpanTime
                    
                    thisIdxCnt<-EndIdx-StartIdx+1   # data points in this Interval  -- used for S11, and SigmaHat
                    Loops<-RowCnt          
                    if (Components==1 && oneCycle[1] > 0){              #  a specified period in oneCycle overrides spectrum calculation  (and Period$Start/End?)
                        Cycles<-oneCycle
                        Loops<-length(Cycles)        
                    } else if (Components>1){
                        Cycles<-oneCycle
                        Loops<-1
                    }
                    #cat("StSpan ",StartSpans[j]," Interval ",Interval," EndSpan ",EndSpan," StrtIdx ",StartIdx, " EndIdx ",EndIdx," thisIdxCnt ",thisIdxCnt," j ",j,"\n")
                    if (window != "noTaper"){          # calculates a window for the Interval
                        # pass exact interval-1!!!!
                        #browser()
                        CatWindowList<-CatWindow(window=window,myData=MyData[StartIdx:EndIdx,Ycol],Time=MyData$time.hour[StartIdx:EndIdx],dataN=Interval_hours,Start=StartSpanTime, debug=Debug)
                        #CatWindowList<-CatWindow(window=window,myData=MyData[,Ycol],Time=MyData$time.hour,dataN=MyData_length,Start=MyData$time.hour[1]+StartSpans[1]-1)    # earlier in CATcosinor
                        #9/13 CatWindowList<-CatWindow(window=window,myData=MyData[,Ycol],Time=MyData$time.hour,dataN=MyData_hours,Start=MyData$time.hour[1]+StartSpans[1]-1) 
                        #CatWindow(window=window,myData=MyData[StartIdx:EndIdx,Ycol],Time=MyData$time.hour[StartIdx:EndIdx],dataN=(EndIdx-StartIdx+1),Start=StartSpanTime)   #later in catcosinor
                        newData<-CatWindowList$data
                    } else { newData<-MyData[StartIdx:EndIdx,Ycol]}
                    
                    #  summary is on the filtered data.
                    sumN[y,j]<-length(newData)  # should be changed to this
                    #sumN<-length(newData[,Ycol])
                    sumLow[y,j]<-min(newData, na.rm=TRUE)
                    sumHi[y,j]<-max(newData, na.rm=TRUE)
                    sumMean[y,j]<-mean(newData, na.rm=TRUE)
                    sumMedian[y,j]<-median(newData, na.rm=TRUE)   
                    
                    sumSD[y,j]<-sqrt(var(newData,y=NULL))
                    sumTab<-tabulate(newData*1000)          #  acts only on integer values, ensure integer;   each value is counted in it's cardinal location     
                    dTtest<-which(sumTab == max(sumTab, na.rm=TRUE))/1000   #  
                    
                    sumT[y,j]<-dt
                    if (length(dTtest)>1){    # if more than 4 are the same as the max dT, use mean instead of median
                        sumMode[y,j] <- mean(dTtest, na.rm=TRUE)
                    } else {sumMode[y,j]<-dTtest   #  ????way to find dt for non-equidistant data
                    }
                    
                    for (i in 1:Loops){
                        # fit sinusoidal model
                        if (oneCycle[1]==0){                  #  if a period has been specified, use present cycle
                            #cycle<-Interval/(i-(StartSpans[j]-1))       #  divisor should start at one each time this loop restarts
                            cycle<-Period$Start/(1+((i-1)*FreqInc))          ####  Interval/(1+((i-1)*FreqInc))  from Interval/i
                        } else cycle<-Cycles[i]
                        
                        #  set subject ID
                        # if (IDcol!="fileName"){    # (this could be taken outside the Y loop and made an outer loop)
                        #     if (!is.numeric(IDcol) ){    # not "fileName" and not numeric is an error!!!! fill with spaces
                        #         SubjectID[i]<-" "
                        #     } else {
                        #         SubjectID[i]<-MyData[StartIdx,IDcol]     #  set for each Loop
                        #     }
                        # }
                        
                        Cycle[i,j,1]<-cycle                    # store for later 
                        yVar[i,j]<-names(MyData)[Y[y]]
                        hours[i,j] <- paste(format(StartSpanTime,nsmall=1)," - ",format(EndSpanTime,nsmall=1),sep="")  # took :59 off because can be decimal hrs;  this shows ALL used
                        cycleCnt<-(EndSpanTime-StartSpanTime)/cycle
                        minPtCnt<-2*cycleCnt+1
                        nPts[i,j] <- EndIdx-StartIdx+1
                        sPts[i,j] <- paste(StartIdx,"-",EndIdx,"\n(",nPts[i,j],")")
                        time[i,j] <- paste(MyData[StartIdx,TimeCol],"-\n",MyData[EndIdx,TimeCol])
                        if (nPts[i,j]<=minPtCnt) {        #  not enough data points for the trial period being analzyed.
                            print(paste(keyMsgs[5],"C=",2*cycleCnt+1,"min pts=",minPtCnt))
                            Err[i,j,]<-paste(Err[i,j,],errKey[5],"min pts=",minPtCnt)   # "@"
                            #next   should not skip calculations  -- calculate as many as possible  jiii
                        }
                        
                        #if (cycle>((MyData$time.hour[EndIdx]+dt)-MyData$time.hour[StartIdx])) {     #  data covers 90% of timeperiod < target period   (8-2-2016 added +dt for when "numeric" hours used -Garth)
                        if (cycle>Interval) {     #  data covers 90% of timeperiod < target period   (8-2-2016 Useing Interval because StartIdx and EndIdx are first and last, not always interval size when units=numeric)
                            if ((.9*cycle)>(MyData$time.hour[EndIdx]-MyData$time.hour[StartIdx])) {    #  Error if much smaller than cycle
                                print(keyMsgs[7])
                                Err[i,j,]<-paste(Err[i,j,],":",errKey[7])   # "@@@"
                                #next
                            } else {  #  Warning if cycle about same size as data Interval
                                print(keyMsgs[8])
                                Err[i,j,]<-paste(Err[i,j,],":",errKey[8])   # "@@"
                            }
                        }
                        
                        # still need to check for gaps
                        CosMatrix<-matrix(data=NA,nrow=CosMatrixDim, ncol=CosMatrixDim)
                        testXmean<-matrix(data=NA,nrow=CosMatrixDim, ncol=Components)
                        YMatrix<-matrix(data=NA,nrow=CosMatrixDim, ncol=1)
                        CosCoMatrix1<-matrix(data=NA,nrow=CosMatrixDim, ncol=thisIdxCnt)
                        CosCoMatrix1[1,]<-rep(x=1,times=thisIdxCnt)         # need to be 0s for CosMatrix<-Sum...
                        
                        for (m in 1:Components){
                            # calculate dTime for use in magnitude:  if all components are multiples of largest, divide one period into modelLen (outside loop)
                            #       otherwise, divide full length of data hours into modelLen;  -->dT     Starttime+dT-->T
                            if (oneCycle[1]>0){    # for all cases except LS
                                if (LCM<=Interval){
                                    dTime<-LCM/modelLen        # length of time units (hours plotted)/#dots [not for each degree of one period]
                                } else {
                                    dTime<-Interval/modelLen
                                }
                                
                                magTime<-vector(mode="integer",length = modelLen)
                                # relative to RefTime, so start at 0
                                magDt<-(modelLen*dTime)              #  this is the #hours from RefTime of the last time
                                magTime<-seq(from=StartSpanTime,by=dTime,to=(StartSpanTime+(magDt-dTime)))       #  changed from 0 to magDt-dTime
                                #magTimeEnd<-RefTime + (magDt*3600) 
                            }  # end oneCycle[1]>0
                            if (Components>1){    #  if components >1, override value of cycle set above.
                                cycle<-Cycles[m]
                                Cycle[i,j,m]<-cycle 
                                #  this is the last time, in time format (convert hrs to sec)
                            }
                            
                            CosCoMatrix1[m*2,]<-cos(MyData$time.hour[StartIdx:EndIdx]*(2*pi)/cycle)     # x1 for  cycle m
                            CosCoMatrix1[m*2+1,]<-sin(MyData$time.hour[StartIdx:EndIdx]*(2*pi)/cycle)   # z1 for  cycle m
                            YMatrix[m*2] <- sum(newData%*%(cos(MyData$time.hour[StartIdx:EndIdx]*(2*pi)/cycle)))     # Y*x1
                            YMatrix[m*2+1] <- sum(newData%*%(sin(MyData$time.hour[StartIdx:EndIdx]*(2*pi)/cycle)))   # Y*z1 
                        }
                        
                        # vector multiplication to build matrix for cosine
                        for (m in 1:CosMatrixDim){      #  columns
                            for (o in 1:CosMatrixDim){    #  rows
                                CosMatrix[o,m]<-sum(CosCoMatrix1[m,] * CosCoMatrix1[o,])
                            }
                        }
                        CosMatrix[1,1]<-thisIdxCnt
                        YMatrix[1] <- sum(newData)
                        
                        mdat<-CosMatrix
                        mdat2<-YMatrix
                        #invert matrix to get the solution
                        if (any(is.na(mdat))) {
                            Err[i,j,]<-paste(Err[i,j,],":",errKey[1])         # "*"     #  Err holds all error symbols for this element;  add this error
                            print(keyMsgs[1])
                            next
                        } else {
                            if (det(mdat)<.00000000001){
                                Err[i,j,]<-paste(Err[i,j,],":",errKey[2])         # "**     #  Err holds all error symbols for this element;  add this error
                                print(keyMsgs[2])
                                next
                            }  else {mdatInv <- solve(mdat)
                            
                            #multiply the inverted matrix by the dependent variable matrix to get vector M, B, G, B2, G2, ...   (x=inv(S)b)
                            coefs <- mdatInv  %*%  mdat2
                            #browser()
                            if (Debug==TRUE || i==Interval){
                                #print(TimeIdx)
                                #print(crsprd)
                                #print(s)
                                #print(model$coefficients)
                                #print(Amp)
                                #cat("data22 ",MyData$time.hour[StartIdx:EndIdx]*(2*pi)/cycle," cycle ",cycle,"\n")
                                #cat("cos22 ",cos(MyData$time.hour[StartIdx:EndIdx]*(2*pi)/cycle)," cycle ",cycle,"\n")
                                #cat("StartIdx ",StartIdx,"  EndIdx ",EndIdx," thisIdxCnt ",thisIdxCnt,"\n")
                                print(mdat)
                                print(mdat2)
                                print(mdatInv)
                                cat("matrix coefs ",coefs,"\n")
                                #browser()         #########   Great Debug location
                            }
                            #calculate MESOR only once, regardless of # of components
                            M[i,j]<-coefs[1]
                            #calculate the Mean of Ys
                            MeanR<-mean(newData, na.rm=TRUE)
                            multi_Model_mag<-multi_Model_Y<-Model_Y <- M[i,j]          # Model_Y = Mesor +  bX + gZ 
                            
                            for (m in 1:Components){
                                beta<-m*2
                                gamma<-m*2+1
                                if (Components>1){    #  if components >1, override value of cycle set above.
                                    cycle<-Cycles[m]
                                    Cycle[i,j,m]<-cycle
                                }
                                #sometimes coefs will be undefined.  skip calculations and move to next cycle
                                if (!is.na(coefs[beta]) && !is.na(coefs[gamma])){
                                    #Calculate Amplitude using B and G   (coefs[beta]=Beta)
                                    A[i,j,m]<-(coefs[beta]^2+coefs[gamma]^2)^0.5
                                    #  cat("A",A[i,j,m],"max",A[maxAmp,j,m],"m",m,"maxAmp",maxAmp)
                                    #preserve the index for the cycle with the maximum amplitude  
                                    #browser()
                                    if (is.na(A[maxAmp,j,m])){   #  this does not capture the maxAmp for multiple components properly
                                        maxAmp<-i
                                    } else if (A[i,j,m]>A[maxAmp,j,m]) {
                                        max3Amp<-max2Amp     # getting top 3 amplitudes
                                        max2Amp<-maxAmp
                                        maxAmp<-i
                                    } else if (A[i,j,m]>A[max2Amp,j,m]) { #  this does not capture the maxAmp for multiple components properly
                                        max3Amp<-max2Amp     # getting top 3 amplitudes
                                        max2Amp<-i
                                    } else if (A[i,j,m]>A[max3Amp,j,m]) {
                                        max3Amp<-i     # getting top 3 amplitudes
                                    }
                                    
                                    # Calculate Model:  Y(t) = M + bX + gZ    where X=cos(t2pi/cycle) and Z = sin(t2pi/cycle)
                                    Model_Y = M[i,j] + (coefs[beta]*cos(2*pi*MyData$time.hour[StartIdx:EndIdx]/cycle)) + (coefs[gamma]*sin(2*pi*MyData$time.hour[StartIdx:EndIdx]/cycle))
                                    if (oneCycle[1]>0){  
                                        #Calculate model of 1 period of longest period for plot:  Y(t) = M + bX + gZ + e(t)where  b = Acosf  and  g = -Asinfand  X = cos(2pt/t)  and Z = sin(2pt/t)
                                        Model_Y_mag[i,] = M[i,j] + (coefs[beta]*cos(2*pi*magTime/cycle)) + (coefs[gamma]*sin(2*pi*magTime/cycle))
                                    }
                                    if (Components>1){      # multiple components model 
                                        multi_Model_Y = multi_Model_Y + (coefs[beta]*cos(2*pi*MyData$time.hour[StartIdx:EndIdx]/cycle)) + (coefs[gamma]*sin(2*pi*MyData$time.hour[StartIdx:EndIdx]/cycle))  
                                        
                                        #Calculate model of 1 period for magnitude:  Y(t) = M + bX + gZ + e(t)where  b = Acosf  and  g = -Asinfand  X = cos(2pt/t)  and Z = sin(2pt/t)
                                        multi_Model_mag = multi_Model_mag + coefs[beta]*cos(2*pi*magTime/Cycles[m]) + coefs[gamma]*sin(2*pi*magTime/Cycles[m])
                                    }
                                    
                                    MYhat_Ymean<- Model_Y-MeanR     #   calculates MSS for single component model
                                    MSS[m]<-sum(MYhat_Ymean^2)        #   calculates MSS for single component model
                                    #  new PR 
                                    yTSS<-newData-MeanR       #  vector Y-Ybar for each Y
                                    TSS_2<-sum((yTSS)^2)       #   sum of Y-Ybar^2
                                    #beta1   *  sum[(Y-Ybar)x1]  * gamma1  *  sum[(Y-Ybar)z1] /TSS^2
                                    #beta2   *  sum[(Y-Ybar)x2]  * gamma2  *  sum[(Y-Ybar)z2] /TSS^2  (one for each component)
                                    newPR[i,j,m]<-100*(coefs[beta]*sum(yTSS*cos(2*pi*MyData$time.hour[StartIdx:EndIdx]/cycle)) + coefs[gamma]*sum(yTSS*sin(2*pi*MyData$time.hour[StartIdx:EndIdx]/cycle)))/TSS_2        #  my Y-Ymean*x1*coef
                                    if (newPR[i,j,m]<0){
                                        Err[i,j,m]<-paste(Err[i,j,m],":",errKey[3])          #  +
                                        print(keyMsgs[3])
                                        #next            should do all calculations as possible;  multiple component PR doesn't get put into PR
                                    }
                                    #calculate phi:  atan(-G/B) + Kpi
                                    ph<-atan(-coefs[gamma]/coefs[beta])
                                    if (coefs[beta]>=0) {phi2<-ph}
                                    if (coefs[beta]<0&coefs[gamma]>=0) {phi2<-ph+pi}
                                    if (coefs[beta]<0&coefs[gamma]<0) {phi2<-ph-pi}
                                    # put in 0 to 2pi range
                                    PHIr[i,j,m]<-phi2
                                    if (phi2<0) { PHIr[i,j,m]<-phi2+(2*pi)} 
                                    if (phi2>(2*pi)) { PHIr[i,j,m]<-phi2-(2*pi)}
                                    # convert to degrees betwen 0 and 360
                                    PHI[i,j,m] <- (PHIr[i,j,m]/(2*pi))*360
                                    if (PHI[i,j,m]>0) {PHI[i,j,m]<-PHI[i,j,m]-360}
                                    
                                } else {  #  end only-if coefs can be calculated  -- this error not needed
                                    Err[i,j,m]<-paste(Err[i,j,m],":",errKey[4])   #  ++
                                    print(keyMsgs[4])
                                    #browser()
                                    next
                                }    #  end only-if coefs can be calculated
                            }   # end coef calculation for each component
                            
                            # sum(Y-meanY)^2 = sum(Y - Yhat)^2 + sum(Yhat - meanY)^2
                            if (Components>1){      # multiple components model 
                                Model_Y <- multi_Model_Y       #  use the full (data derived) model for calculations from here on
                                # calculate MSS for the full model
                                MYhat_Ymean<- Model_Y-MeanR
                                multi_MSS<-sum(MYhat_Ymean^2)
                            }
                            MY_Yhat <- newData- Model_Y     # my Y-Yhat
                            #browser()
                            RSS<-sum(MY_Yhat^2)
                            if (RSS==0) {
                                print(keyMsgs[6])
                                Err[i,j,]<-paste(Err[i,j,],":",errKey[6])         # "@@"
                                next
                            }
                            df1<-Components*2
                            df2<-thisIdxCnt-CosMatrixDim                      # N-2*Components+1
                            if (df2<0){
                                print(paste(keyMsgs[9],"2*Components+1=",CosMatrixDim,"thisIdxCnt=",thisIdxCnt))
                                Err[i,j,]<-paste(Err[i,j,],":",errKey[9],"thisIdxCnt=",thisIdxCnt)   # "*@@"
                                next
                            }else if (df2<1){
                                print(paste(keyMsgs[10],"2*Components+1=",CosMatrixDim,"thisIdxCnt=",thisIdxCnt))
                                Err[i,j,]<-paste(Err[i,j,],":",errKey[10],"thisIdxCnt=",thisIdxCnt)   # "*@@"
                                df2<-NA   #  this will cause se's and F and P to also be NA
                            }
                            
                            if (Components>1){      # multiple components model 
                                multi_F<-(multi_MSS/df1)/(RSS/(df2))       #  df1 = 2, since MSS is for one component;  df2 is for full model
                                multi_P[i,j]<-1-pf(multi_F,df1=df1, df2=df2)              # Fisher-Snedecor (F ) (X2 ) 
                                multi_PR[i,j]<-(multi_MSS/(RSS+multi_MSS))*100
                                cat("F-multi:",multi_F,"\n")
                            } 
                            
                            SigmaHat<-(RSS/(df2))^.5
                            #cat(thisIdxCnt," sHat ",SigmaHat,"\t",PR[i,j],"\t",F[i,j],"\t ",P[i,j],"\t ")
                            
                            mesor_se[i,j] <-SigmaHat * mdatInv[1,1]^.5
                            
                            #   for testing alternate formula for individual P
                            for (m in 1:Components){
                                F[i,j,m]<-(MSS[m]/2)/(RSS/(df2))       #  df1 = 2, since MSS is for one component;  df2 is for full model
                                P[i,j,m]<-1-pf(F[i,j,m],df1=2, df2=df2)              # Fisher-Snedecor (F ) (X2 )
                                if (Debug==TRUE){
                                    cat("F:df1-df2",F[i,j,m],2,df2,"MSS",MSS[m],"RSS",RSS,"\n")
                                }
                                # testing alternate formula for individual P
                                tInd1<-2*m
                                tInd2<-2*m+1
                                beta<-2*m
                                gamma<-2*m+1
                                #   excludes diagonals of inverse matrix, which may be important
                                test<-(mdatInv[tInd2,tInd2]*coefs[beta]^2)+2*mdatInv[tInd1,tInd2]*coefs[beta]*coefs[gamma]+mdatInv[tInd1,tInd1]*coefs[gamma]^2
                                test2<-test/(mdatInv[tInd1,tInd1]*mdatInv[tInd2,tInd2]-mdatInv[tInd1,tInd2]^2)
                                test3<-test2/(2*SigmaHat^2)
                                testPg[i,j,m]<-1-pf(test3,df1=2, df2=df2) 
                                
                                #  alternative calculation for F    p 409 Bingham [33] 
                                testXmean<-(CosMatrix[1,tInd1])/thisIdxCnt    #  X mean
                                testZmean<-(CosMatrix[1,tInd2])/thisIdxCnt    #  Z mean
                                testX2<-sum((CosCoMatrix1[tInd1,] * CosCoMatrix1[1,]-testXmean)^2)   #  X 
                                testZ2<-sum((CosCoMatrix1[tInd2,] * CosCoMatrix1[1,]-testZmean)^2)   #  Z 
                                testXZ<-sum((CosCoMatrix1[tInd1,] * CosCoMatrix1[1,]-testXmean)*(CosCoMatrix1[tInd2,] * CosCoMatrix1[1,]-testZmean))   #  XZ 
                                testF<-(testX2 * coefs[beta]^2 + 2*testXZ * coefs[beta] * coefs[gamma] + testZ2 * coefs[gamma]^2)/(2*SigmaHat^2)
                                testPb[i,j,m]<-1-pf(testF,df1=2, df2=df2)
                                
                                if (Components>1){      # multiple components model 
                                    PR[i,j,m]<-(MSS[m]/(RSS+multi_MSS))*100
                                } else {
                                    PR[i,j,m]<-(MSS[m]/(RSS+MSS[m]))*100
                                }
                                
                                #calculate standard error for mesor, amplitude, phi -- convert se for phi to degrees
                                amp_se[i,j,m] <-SigmaHat * (mdatInv[beta,beta] * cos(PHIr[i,j,m])^2 - (2 * mdatInv[beta,gamma] * cos(PHIr[i,j,m]) * sin(PHIr[i,j,m])) + mdatInv[gamma,gamma] * sin(PHIr[i,j,m])^2)^.5
                                phi_se[i,j,m] <- (SigmaHat * (mdatInv[beta,beta] * sin(PHIr[i,j,m])^2 + (2 * mdatInv[beta,gamma] * cos(PHIr[i,j,m]) * sin(PHIr[i,j,m])) + mdatInv[gamma,gamma] * cos(PHIr[i,j,m])^2)^.5)/A[i,j,m]
                                phi_se[i,j,m]<-phi_se[i,j,m]*180/pi       # Convert to degrees
                                
                                if (Debug==TRUE){
                                    cat(" seM ",mesor_se[i,j]," seA ",amp_se[i,j,m]," sePHI ",phi_se[i,j,m]," \n")
                                    cat("test3",m,"   ",test3,"\n")
                                    cat("testPg",m,"   ",testPg[i,j,m],"\n")
                                    cat("testF",m,"   ",testF,"\n")
                                    cat("testPb",m,"   ",testPb[i,j,m],"\n")
                                    cat(thisIdxCnt,"\t","sigma ",SigmaHat,"\t",PR[i,j,m],"\t",F[i,j,m],"\t ",P[i,j,m],"\t \n")
                                    cat(StartSpans[j],"\t",i,"\t",format(Cycle[i,j,m],width=3), "\t ",M[i,j],"\t", format(mesor_se[i,j],digits=8),"\t",A[i,j,m],"\t",format(amp_se[i,j,m],digits=8),"\t",format(PHI[i,j,m],digits=8),"\t",format(phi_se[i,j,m],digits=8),"\n")
                                }
                                # if (oneCycle[1]>0){   # print this cycle result if  oneCycle parameter set is>0
                                #     #ht<-hdrRow+.5-(i+(j*htVar))+jht-(m-1)
                                #     ht<-hdrRow+.5-(i+(j*htVar))+jht-(m-1)
                                #     #browser()
                                #     #  prints period and estimates for each period in a range
                                #     text(c(-.4,1.2,2,3,3.9,5.1,6.2,7.2,8.2,9.2),c(ht,ht,ht,ht,ht,ht,ht,ht,ht,ht,ht),labels=c(paste(StartSpans[j]," - ",format(StartSpans[j]+Interval-1,nsmall=3)),format(Cycle[i,j,m],nsmall=2),format(P[i,j,m],digits=3,nsmall=3),format(newPR[i,j,m],digits=5),format(M[i,j],digits=8), format(mesor_se[i,j],digits=6),format(A[i,j,m],digits=6),format(amp_se[i,j,m],digits=6),format(PHI[i,j,m],digits=6),format(phi_se[i,j,m],digits=6)),cex=(cexVar*1.2),adj = c(0,0))
                                # }   #   format(F[i,j,m],digits=5),   removed 
                            }  # end  Components loop
                            
                            
                            }  #  end check for matrix validity
                        }   # end matrix cannot be inverted
                    }    #  end i=Loops for loop
                    
                    
                }    #   END for (j in 1:Progression_end){  ,seq(from=magTime[1],by=dTime,to=tail(magTime,1))
                
                tmp <- c(Period = Cycle[1,1,1],
                         P = P[1,1,1],
                         PR = PR[1,1,1],
                         Mesor = M[1,1],
                         M.s.e. = mesor_se[1,1],
                         Amp = A[1,1,1],
                         A.s.e. = amp_se[1,1,1],
                         Phi = PHI[1,1,1],
                         P.s.e = phi_se[1,1,1])
                
                if(!is.na(GrpCol)){
                    grp <- as.character(sapply(mydata[,GrpCol,drop=FALSE], unique))
                    # print(grp)
                    # print(str(output[[y]]))
                    # print(str(output[[y]][sub,1:length(GrpCol)]))
                    output[[y]][sub,1:length(GrpCol)] <- grp 
                }
                
                idx <- (ncol(output[[y]]) - 9)
                output[[y]][sub, idx] <- paste(StartSpans[1], "-", Interval - 1, sep = "")
                output[[y]][sub, (idx + 1):ncol(output[[y]])] <- tmp
                # print("row name assignment?")
                # print(subjectIDs[sub])
                row.names(output[[y]])[sub] <- subjectIDs[sub]
                
                
            }    #    END for (y in 1:Ys_end){
        }
        
        
        
        return(output)
    }
