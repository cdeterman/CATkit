
library(CATkit)

test_that("pnorm works", {
    # print(cpp_pnorm(1))
    
    print(getwd())
    # save(mydata, file = 'tests/testthat/test_data.rda')
    load('test_data.rda')
    
    # mydata <- read.table("clipboard", sep = ",")
    # mydata <- as.matrix(mydata[, 1:4])
    
    # print(mydata)
    
    # .Fortran("test", mydata, nrow(mydata), ncol(mydata))
    # print(mydata)
    
    # print("try c++ again")
    # 
    # cpp_test(mydata)
    # print(mydata)
    
    # str(mydata)
    # gsub('^.*?-', '', mydata$V5[1])
    # mydata$V6 <- gsub('^.*?-', '', mydata$V5)
    # split(mydata, mydata[,6])
    
    mylist <- lapply(split(mydata, mydata[,6]), function(x) pmc(as.matrix(x[,1:4]), perd = 24))
    
    print(str(mylist))
    
    out2 <- do.call('rbind', mylist)
    
    # str(out)
    print(out2)
    
    # print(pmc(mydata, perd = 24))
})
