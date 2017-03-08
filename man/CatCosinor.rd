\name{CATCosinor}
\alias{CATCosinor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cosinor Analysis
        %%  ~~function to do ... ~~
}
\description{ Perform one of various cosinor analyses:  single cosinor, least squares spectrum by cosinor, multiple-component cosinor, progressive single cosinor, progressive least square spectrum, progressive multiple-component cosinor.
              
              (** In time formats such as "_Y_m_d_H_M", the underscore, "_", symbol represents a percent sign - which is not allowed in this help file)
              %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{  
  CATCosinor(TimeCol=1,Y=2, Components=1, window="noTaper", RefDateTime=NA, 
    timeFormat="\%Y\%m\%d\%H\%M", RangeDateTime=list(Start=NA, End=NA), Units="hours", 
    dt=0, Progressive=list(Interval=0, Increment=0), Period=list(Set=0,Start=0,
    Increment=1,End=0),header=FALSE, Skip=0, Colors="BW", Graphics="pdf",Output=list
    (Txt=FALSE,Dat=TRUE,Doc=TRUE,Graphs=FALSE),yLabel="", Console=FALSE,Debug=FALSE,
    IDcol="fileName",fileName=fileName,functionName="") 
  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TimeCol}{Column number where time is found. Time format of the TimeCol is given by timeFormat parameter. Time is read to the minute. Seconds are discarded. Date and Time can be in any column, in any of three formats (seconds are discarded):  
                   1)  numeric (in number of hours from starting time); 
2)  a one column date/time format, where the format is specified by the parameter timeFormat; or
3)  a two column format consisting of Date in the first of the 2 columns specified, and Time in the second:  c(3,6);  where the format is specified by the parameter timeFormat

%%     ~~Describe \code ~~
  }

\item{Y}{The column number(s) of the data to be analyzed.  This is a numeric, either scalar or vector.  May use any valid R expression, such as c(4,5) or a single number.
         %%     ~~Describe \code ~~
}
\item{Components}{Default=1.  Indicates if this is a single or multiple component cosinor analysis, where the number of components is specified (>0).  If doing a single component cosinor, set Components=1.  If doing a multiple components model, set Components equal to the number of frequencies in the model.  
                  %%     ~~Describe \code here~~
}
\item{window}{Valid windowing function to be applied are:  "noTaper","Hanning","Hamming","Bartlett","Blackman"
              %%     ~~Describe \code{ } here~~
}
\item{RefDateTime}{Date used as reference, and subtracted from all data dates, to make the number smaller.  **Must be in the same time zone!!!!!**
                     if RefDateTime = NA, uses the 1st date of the data as the RefDateTime  
                   if RefDateTime = 0, uses midnight of the same day as the data starts 
                   %%     ~~Describe \code{ } here~~
}
\item{timeFormat}{Can be "numeric", or any valid R time conversion specification,  i.e., "\%Y\%m\%d\%H\%M".  See strptime for conversion specifications.
}

\item{RangeDateTime}{specify in the form of a list:  RangeDateTime=list(Start=12, End=0)      
                     $Start:  Analysis is started at this date/time.  May be before 1st data date.
                     if Start = NA, uses the 1st data point as the StartDate
                     if Start =  0, uses Midnight at the start of the 1st date of the data as the StartDate 
                     $End:    Analysis ends at this date/time.  May be after last data date.
                     if End =  NA, use the last data point as the EndDate
                     if End =  0 uses the midnight at the end of the last date of the data as the EndDate   
                     %%     ~~Describe \code{k} here~~
}
\item{Units}{Units (Hour, Year, Week or Day) for use by Interval and Increment arguments, as well as Period arguments (Note: only Hour is currently implemented.)
             %%     ~~Describe \code{yLab} here~~
}
\item{dt}{When equidistant data, dt indicates the sampling interval.  If dt =0, no periodogram is done.  Data is assumed to be equidistant when this is nonzero.
          %%     ~~Describe \code{modulo} here~~
}
\item{Progressive}{specify in the form of a list:  Progressive=list(Interval=0, Increment=0)
                   $Interval: length of the time span being analyzed (in Units) -- multiple spans calculated
                   If 0, assumes no progression, Interval is set to the full dataset length, and Increment = full data set  
                   $Increment: number of Days, Wks or Yrs  (uses same unit as set for Interval) to shift forward for each successive Interval analyses
                   If 0, assumes no progression, Interval is set to the full dataset length, and Increment = full data set 
                   %%     ~~Describe \code{ } here~~
}
\item{Period}{specify in the form of a list:   Period=list(Set=0,Start=0,Increment=1,End=0)
              $Start : [only used if $Set=0]  First (and largest) period to calculate, in units of Days, Wks or Yrs (as set by Units);  (Interval/1).   0 is Default: the full time range of the data file is analyzed [in hours] (RangeDateTime$End-RangeDateTime$Start)= MyData_length; or if progressive, Interval/1;  
              $Increment : [only used if $Set=0] Increment to starting period, in units of Days, Wks or Yrs (as set by Units).  Defaults to 1;   0 is invalid  -- default will be used
              $End : [only used if $Set=0]  Last (and smallest) period to calculate, in units of Days, Wks or Yrs (as set by Units), EXCLUSIVE.  Defaults to 2*dt or 4;  (1 is too small) 0 is invalid  -- default will be used
              $Set :  If Set=0, a series of periods are analyzed (spectrum) according to Period$Start, $Increment, $End (or their default value, if not specified).  If not equal to 0, overrides Period$Start and $Increment, to completely specify a set of periods to analyze (as set by Units), ending at $End.  Can be in the format of a single number, a vector, or R expression:  c(1,2,3) or c(8,12:24) or seq(1,50, by=.75).  When Components=1, each period specified in the vector will be assessed by cosinor independently.  When parameter Components is >1, Period$Set must have a corresponding number of components, which are assessed together in a multiple component cosinor.    When 0, only the largest period, or the largest period per Interval, from a spectrum is listed on page 1 of the graphics output file.  Otherwise, all periods are displayed in the graphic file.    
              %%     ~~Describe \code{ } here~~
}
\item{header}{ T/F to indicate if the file has a header.  Headers are used as variable names.
               %%     ~~Describe \code{Skip} here~~
}
\item{Skip}{ How many lines to skip at the top of the file (before the header). The header will be read next, after Skip lines are skipped (if header=TRUE). 
             %%     ~~Describe \code{ } here~~
}
\item{Colors}{Colors for a heatmap.  "Heat" renders the heat map in colors;  "BW" renders the heatmap in grayscale
              %%     ~~Describe \code{ } here~~
}
\item{Graphics}{The main results of CATCosinor are sent to a graphic file when Console=F.  Default file output type is "pdf".   Possible values: "jpg, pdf, tiff, png".   
                %%     ~~Describe \code{ } here~~
}
                \item{Output}{Specify in the form of a list:  list(Txt=F,Dat=T, Doc=T,Graphs=F).  $Txt=T will capture the console output to a .txt file.  $Dat=T will generate a .txt computer-readable file of tab delimited cosinor results: MESOR, Amplitude, Phase, standard error values, and others. It is useful for subsequent processing (by CAT, or excel, for example).  $Doc=T will generate a nicely-formatted .rtf file, also readable by Word.  $Graphs will enable a set of graphs plotting Data, Model, MESOR, Amplitude, Acrophase over time, or a heatmap.  The exact files generated will vary with the functions performed.  Heatmaps are only created when a progressive analysis of (single component) multiple frequencies is performed. $Graphs=F will disable printing of graphs for faster processing, if you do not need the files.
                %%     ~~Describe \code{ } here~~
                }
\item{yLabel}{Label for the Y axis on data graphs.  If this is left blank (the default="") then the column headers are used for Y axis label.
              %%     ~~Describe \code{yLabel} here~~
  }
                \item{Console}{Default is F.  When Console=T, output will be redirected to the RStudio Console, instead of an output file.  Otherwise, sent to file type indicated in Graphics parameter.
                %%     ~~Describe \code{ } here~~
                }
                \item{Debug}{ Turn on when you want to see more diagnostic output for programming debug purposes.  
                %%     ~~Describe \code{ } here~~
                }
                \item{IDcol}{What to use as a unique identifier:  either "fileName" or a column number.  Default = "fileName"
                %%     ~~Describe \code{ } here~~
                }
                \item{fileName}{The path to a data file to be processed.  If this is blank, CATCosinor will solicit a file from the user by a dialogue box. 
                %%     ~~Describe \code{ } here~~
                }
                \item{functionName}{Brief user-defined name for this run, used to help distinguish it from other runs.  Name is printed in output files.
                %%     ~~Describe \code{ } here~~
                }
}
                \section{Input Data requirements:}{  
                There is no need for data to be equidistant for the Cosinor.  Data columns specified in Y are expected to be numeric.  
                Input data files can be tab delimited, or comma delimited.  There may be multiple columns of data for processing.
                }
                \section{Output Data:}{  
                Sample graphics output file can be found in the Output section on the CAT web site for a sample of a full output file.  All graphs contain the input data filename to clearly identify the data file under analysis, and a timestamp to show the time of analysis.  Each graph also lists the column name being analyzed (or averaged), as well as the starting and ending times of analysis, as they vary slightly from the full data set.
                Graphic output files can be in JPG, PDF, TIFF or PNG.  TIFF and PNG are higher resolution than jpeg and PDF.  Plots in PNG and JPEG formats can easily be converted to many other bitmap formats, and both can be displayed in modern web browsers. The PNG format is lossless and is best for line diagrams and blocks of color. The JPEG format is lossy, but may be useful for image plots, for example. It is most often used in html web pages.  TIFF is a meta-format: the default format written by tiff is lossless and stores RGB (and alpha where appropriate) values uncompressed-such files are widely accepted, which is their main virtue over PNG.
                %%  ~~ If necessary, more details than the description above ~~
                }
                
                \references{http://564394709114639785.weebly.com/running-cat.html
                %% ~put references to the literature/web site here ~
                }
                \author{Cathy Lee Gierke, Germaine Cornelissen-Guillaume, John A Lindgren, Ruth Ann Helget
                
                Maintainer: \email{ <leegi001@umn.edu>}
                }
                %%  ~~who you are~~
                
                %% ~Make other sections like Warning with \section{Warning }{....} ~
                
                \examples{
                #Data is systolic and dialstolic blood pressure; and heart rate at 
                #        30 minute intervals.
                #---------------------  Vignette0    multi-component cosinor.  
                # multiple-components cosinor analysis of 24, 12 and 8 hours.  Columns 
                #       4,6,7 are processed.

                #  Normally you would use these lines to read a file for use in CATkit 
                #        (use filePath format for your OS)
                # filePath<-"~/file/path/Installing CAT/Vignette0"          # use for mac
                # filePath<-"c:\\file\\path\\Installing CAT\\Vignette0"     # use for PC
                # fileName<-file.path(filePath,'CLGi001.dat')
                
                # this line is used instead of the above since it is part of a package.
                fileName<-system.file("extdata", "CLGi001.dat", package = "CATkit")
                #  
                CATCosinor(TimeCol=2,Y=c(4,6,7), Components=3, window="noTaper", RefDateTime
                  ="201302030000",  timeFormat="\%Y\%m\%d\%H\%M",RangeDateTime =list(Start=0, 
                  End=0), Units='hours', dt=0, Progressive=list(Interval=0, Increment=0),
                  Period=list(Set=c(24,12,8),Start=0,Increment=1,End=0),header=FALSE, 
                  Skip=0,Colors="BW",Graphics="pdf",Output=list(Txt=FALSE,Dat=TRUE,
                  Doc=TRUE,Graphs=TRUE),Console=FALSE,Debug=FALSE,IDcol="fileName", 
                  fileName=fileName,functionName='Vignette0')
                
                #---------------------   Vignette1     Least Square spectrum
                # 
                #  Normally you would use these lines to read a file for use in CATkit 
                #     (use filePath format for your OS)
                # filePath<-"~/file/path/Installing CAT/Vignette1"          # use for mac
                # filePath<-"c:\\file\\path\\Installing CAT\\Vignette1"     # use for PC
                # fileName<-file.path(filePath,'Signal10-20.txt')
                
                # this line is used instead of the above since it is part of a package.
                fileName<-system.file("extdata", "Signal10-20.txt", package = "CATkit")
                #
                CATCosinor(TimeCol=1,Y=2, Components=1, window="noTaper", RefDateTime="0",
                timeFormat="\%Y\%m\%d\%H\%M", RangeDateTime=list(Start=0, End=0), 
                Units='hours', dt=0, Progressive=list(Interval=0, Increment=0), Period=
                list(Set=0,Start=144,Increment=1,End=4),header=FALSE, Skip=0, Colors="BW",
                Graphics="pdf",Output=list(Txt=FALSE,Dat=TRUE,Doc=TRUE,Graphs=TRUE),
                Console=FALSE,Debug=FALSE,fileName=fileName,functionName='Vignette1') 
                
                # ---------------------   Vignette2     progressive least squares spectrum
                #
                #  Normally you would use these lines to read a file for use in CATkit 
                #       (use filePath format for your OS)
                # filePath<-"~/file/path/Installing CAT/Vignette2"          # use for mac
                # filePath<-"c:\\file\\path\\Installing CAT\\Vignette2"     # use for PC
                # fileName<-file.path(filePath,'FWedited.txt')
                #
                # this line is used instead of the above since it is part of a package.
                fileName<-system.file("extdata", "FWedited.txt", package = "CATkit")
                #
                CATCosinor(TimeCol=1,Y=2, Components=1, window="noTaper",RefDateTime 
                ="199210192152",  timeFormat="\%Y\%m\%d\%H\%M", RangeDateTime= list(Start
                ="199210200000", End="199211300000"), Units='hours', dt=0, Progressive=list
                (Interval=168, Increment=24), Period=list(Set=0, Start=168,Increment=.5,
                End=9.5),header=FALSE, Skip=0, Colors="BW", Graphics="pdf",Output=list
                (Txt=FALSE,Dat=TRUE,Doc=FALSE,Graphs=TRUE),Console=FALSE,Debug=FALSE,
                fileName=fileName,functionName='Vignette2')
                }
                
                
                % Add one or more standard keywords, see file 'KEYWORDS' in the
                % R documentation directory.
                \keyword{ ~ Cosinor }
                \keyword{ ~ multi-component }
                \keyword{ ~ single-component }
                \keyword{ ~ regression }
                \keyword{ ~ least squares }
                \keyword{ ~ spectrum }
                \keyword{ ~ chronomics }