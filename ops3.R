#########################################################################
### Online Portfolio Selection -- Gradient Descend ######################
#########################################################################

setwd('/Users/lizhuo/Documents/UVA_Statistics/6003 Optimization/project')
library(quantmod)

symbols <- c('MMM', 'AXP', 'AAPL', 'BA', 'CAT', 'CVX', 'CSCO', 'KO', 'DIS', 'XOM', 'GE',
            'GS', 'HD', 'IBM', 'INTC', 'JNJ', 'JPM', 'MCD', 'MRK', 'MSFT', 'NKE', 'PFE', 'PG',
            'TRV', 'UTX', 'UNH', 'VZ', 'V', 'WMT')

for (i in 1:length(symbols)) getSymbols(symbols[i], from = '2007-01-03', to = '2017-11-30')


# Building weekly returns for each of the stocks
data <-sapply(symbols, function(x) Ad(to.weekly(get(x))))
data <- Reduce(cbind, data)
data_returns <- apply(data, 2, function(x) diff(log(x))) #log returns
colnames(data_returns) <- symbols
data_returns[is.na(data_returns)] <- 0 # VISA hasnt negotiations between 2007 and 2008
#write.csv(data_returns, "return.csv")

# try to get more stocks
dat <- read.csv("stock.csv")
Names <-as.character(unique(colnames(dat[,3:491])))
Names <- Names[Names !="AES"]
Names <- Names[Names !="BHF"]
Names <- Names[Names !="ESS"]
Names <- Names[Names !="MMC"]
Names <- Names[Names !="UA"]
for (i in 1:length(Names)) getSymbols(Names[i], from = '2007-01-03', to = '2017-11-30')
df <-sapply(Names, function(x) Ad(to.weekly(get(x))))
df <- Reduce(cbind, df)
df_day <-sapply(Names, function(x) Ad(to.daily(get(x))))
df_month <-sapply(Names, function(x) Ad(to.monthly(get(x))))
df <- Reduce(cbind, df)
df_day <- Reduce(cbind, df_day)
df_month <- Reduce(cbind, df_month)
nalist <- which(is.na(df[1,]))
df <- df[,-which(is.na(df[1,]))]
df <- df[rowSums(is.na(df)) == 0 , ]
sum(is.na(df))
df_returns <- apply(df, 2, function(x) diff(log(x))) #log returns
colnames(df_returns) <- Names[-nalist]
df_returns[is.na(df_returns)] <- 0 # VISA hasnt negotiations between 2007 and 2008
#write.csv(df_returns, "df_return.csv")
nalist <- which(is.na(df_day[1,]))
df_day <- df_day[,-which(is.na(df_day[1,]))]
df_day <- df_day[rowSums(is.na(df_day)) == 0 , ]
sum(is.na(df_day))
day_returns <- apply(df_day, 2, function(x) diff(log(x))) #log returns
colnames(day_returns) <- Names[-nalist]
day_returns[is.na(day_returns)] <- 0 # VISA hasnt negotiations between 2007 and 2008

nalist <- which(is.na(df_month[1,]))
df_month <- df_month[,-which(is.na(df_month[1,]))]
df_month <- df_month[rowSums(is.na(df_month)) == 0 , ]
sum(is.na(df_month))
month_returns <- apply(df_month, 2, function(x) diff(log(x))) #log returns
colnames(month_returns) <- Names[-nalist]
month_returns[is.na(month_returns)] <- 0 # VISA hasnt negotiations between 2007 and 2008

## random sample 50 from sp500
set.seed(9)
ind <- sample(dim(df_returns)[2], 50)
df50 <- df_returns[,ind]
ind <- sample(dim(df_returns)[2], 45)
df45 <- df_returns[,ind]
ind <- sample(dim(df_returns)[2], 40)
df40 <- df_returns[,ind]
ind <- sample(dim(df_returns)[2], 35)
df35 <- df_returns[,ind]
ind <- sample(dim(df_returns)[2], 30)
df30 <- df_returns[,ind]
ind <- sample(dim(df_returns)[2], 25)
df25 <- df_returns[,ind]
ind <- sample(dim(df_returns)[2], 20)
df20 <- df_returns[,ind]
ind <- sample(dim(df_returns)[2], 15)
df15 <- df_returns[,ind]
ind <- sample(dim(df_returns)[2], 10)
df10 <- df_returns[,ind]
ind <- sample(dim(df_returns)[2], 5)
df5 <- df_returns[,ind]

library(quadprogXT)

OGD <- function(base, eta) {
  
  # Gradient of Regret Function
  gradient = function(b, p, r) b + r/(p%*%r)
  
  # Projection onto viable Set
  proj = function(p) {
    
    Dmat = diag(length(p))
    Amat = cbind(diag(rep(1, length(p))), -1)
    bvec = c(rep(0, length(p)), -1)
    
    fit = solveQPXT(Dmat = Dmat, dvec = p, Amat = Amat, bvec = bvec)
    
    return(fit$solution)
  }
  
  T = nrow(base)
  N = ncol(base)
  
  r = as.matrix(base) + 1 # this is because the algo doesnt work directly with log returns
  p = matrix(0, nrow = N, ncol = T); p[,1] = 1/N # initial portfolio
  b = matrix(0, nrow = N, ncol = T); b[,1] = 0
  
  for (i in 2:T) {
    b[,i] = gradient(b[,i-1], p[,i-1], r[i-1,]) # calculating gradient
    p.aux = p[,i-1] + eta*b[,i] # what we would like to play
    p[,i] = proj(p.aux) # projection in the viable set
  }
  
  return(list('portfolio' = p,'gradient' = b))
}

# testing two etas
#portfolio1 <- OGD(base = data_returns, eta = 1/100)
#portfolio2 <- OGD(base = data_returns, eta = 1/1000)

portfolio3 <- OGD(base = df_returns, eta = 1/1000)
#portfolio4 <- OGD(base = df_returns, eta = 1/10000)

## test against portfolio size
portfolio5 <- OGD(base = df5, eta=1/1000)
portfolio10 <- OGD(base = df10, eta=1/1000)
portfolio15 <- OGD(base = df15, eta=1/1000)
portfolio20 <- OGD(base = df20, eta=1/1000)
portfolio25 <- OGD(base = df25, eta=1/1000)
portfolio30 <- OGD(base = df30, eta=1/1000)
portfolio35 <- OGD(base = df35, eta=1/1000)
portfolio40 <- OGD(base = df40, eta=1/1000)
portfolio45 <- OGD(base = df45, eta=1/1000)
portfolio50 <- OGD(base = df50, eta=1/1000)


compound_return <- function(portfolio, r){
  
  return_OGD <- c(); return_OGD[1] = portfolio$portfolio[,1]%*%r[1,]
  portfolio_value <- c(); portfolio_value[1] = 1 + portfolio$portfolio[,1]%*%r[1,]
  
  for (i in 2:nrow(r)){
    return_OGD[i] = portfolio$portfolio[,i]%*%r[i,]
    portfolio_value[i] = portfolio_value[i-1]*(1 + return_OGD[i])
  }  
  
  return(list('value' = portfolio_value, 'return' = return_OGD))
}

# Our portfolios
#portfolio_value1 <- compound_return(portfolio1, data_returns)
#portfolio_value2 <- compound_return(portfolio2, data_returns)
portfolio_value3 <- compound_return(portfolio3, df_returns)
#portfolio_value4 <- compound_return(portfolio4, df_returns)


portfoliod <- OGD(base = day_returns, eta = 1/1000)
portfoliom <- OGD(base = month_returns, eta = 1/1000)

portfoliod_value <- compound_return(portfoliod, day_returns)
portfoliom_value <- compound_return(portfoliom, month_returns)

n <- length(portfolio_value3$return)
#####################################################################################
### test against trading frequency
w <- seq(1, 561)
m <- seq(1,129)
jpeg('freq.jpg')
plot(portfoliod_value$value, type='l', ylim=c(0.4, 2.5), lwd=2, col="coral3",
     xlab="Trading days", ylab="Portfolio Value")
lines(w*7, portfolio_value3$value, type='l', col="brown4", lwd=2)
lines(m*21, portfoliom_value$value, type='l', col="bisque3", lwd=2)
legend("topleft",legend = c("Daily", "Weekly", "Monthly"), 
       col = c("coral3","brown4","bisque3"), lty = 1, cex=0.8)
dev.off()

# APY
portfolio_APYd <- ((cumprod(1+portfoliod_value$return)[2716])^(1/2716))^252-1
portfolio_APY3 <- ((cumprod(1+portfolio_value3$return)[561])^(1/561))^52-1
portfolio_APYm <- ((cumprod(1+portfoliom_value$return)[129])^(1/129))^12-1


# Sharpe Ratio
Sharpe_ratio <- (portfolio_APY3-0.0074)/sd(portfolio_value3$return)
Sharpe_ratiod <- (portfolio_APYd-0.0074)/sd(portfoliod_value$return)
Sharpe_ratiom <- (portfolio_APYm-0.0074)/sd(portfoliom_value$return)



## test against portfolio size
portfolio_value5 <- compound_return(portfolio5, df5)
portfolio_value10 <- compound_return(portfolio10, df10)
portfolio_value15 <- compound_return(portfolio15, df15)
portfolio_value20 <- compound_return(portfolio20, df20)
portfolio_value25 <- compound_return(portfolio25, df25)
portfolio_value30 <- compound_return(portfolio30, df30)
portfolio_value35 <- compound_return(portfolio35, df35)
portfolio_value40 <- compound_return(portfolio40, df40)
portfolio_value45 <- compound_return(portfolio45, df45)
portfolio_value50 <- compound_return(portfolio50, df50)

n <- length(portfolio_value3$return)
## final portfolio value of different size
size <- c(portfolio_value5$value[564], portfolio_value10$value[564], portfolio_value15$value[564],
          portfolio_value20$value[564], portfolio_value25$value[564], portfolio_value30$value[564],
          portfolio_value35$value[564], portfolio_value40$value[564], portfolio_value45$value[564], 
          portfolio_value50$value[564])

size


#########################################################################
## Sample size
files <- list()
for (i in 1:10){
  # sample stocks
  ind50 <- sample(dim(df_returns)[2], 50)
  assign(paste("df50",i,sep="_"), (df_returns[,ind50]))
  ind45 <- sample(dim(df_returns)[2], 45)
  assign(paste("df45",i,sep="_"), (df_returns[,ind45]))
  ind40 <- sample(dim(df_returns)[2], 40)
  assign(paste("df40",i,sep="_"), (df_returns[,ind40]))
  ind35 <- sample(dim(df_returns)[2], 35)
  assign(paste("df35",i,sep="_"), (df_returns[,ind35]))
  ind30 <- sample(dim(df_returns)[2], 30)
  assign(paste("df30",i,sep="_"), (df_returns[,ind30]))
  ind25 <- sample(dim(df_returns)[2], 25)
  assign(paste("df25",i,sep="_"), (df_returns[,ind25]))
  ind20 <- sample(dim(df_returns)[2], 20)
  assign(paste("df20",i,sep="_"), (df_returns[,ind20]))
  ind15 <- sample(dim(df_returns)[2], 15)
  assign(paste("df15",i,sep="_"), (df_returns[,ind15]))
  ind10 <- sample(dim(df_returns)[2], 10)
  assign(paste("df10",i,sep="_"), (df_returns[,ind10]))
  ind5 <- sample(dim(df_returns)[2], 5)
  assign(paste("df5",i,sep="_"), (df_returns[,ind5]))
  #files <- c(files, (paste("df50",i,sep="_")))
  #files <- c(files, (paste("df45",i,sep="_")))
  #files <- c(files, (paste("df40",i,sep="_")))
  #files <- c(files, (paste("df35",i,sep="_")))
  #files <- c(files, (paste("df30",i,sep="_")))
  #files <- c(files, (paste("df25",i,sep="_")))
  #files <- c(files, (paste("df20",i,sep="_")))
  #files <- c(files, (paste("df15",i,sep="_")))
  #files <- c(files, (paste("df10",i,sep="_")))
  #files <- c(files, (paste("df5",i,sep="_")))
  
}

i <- 0
for (df in list(df5_1, df5_2, df5_3, df5_4, df5_5, df5_6, df5_7, df5_8, df5_9, df5_10)){
  assign(paste("p5",i,sep="_"), compound_return(OGD(base = df, eta=1/1000), df))
  i <- i+1
}

i <- 0
for (df in list(df10_1, df10_2, df10_3, df10_4, df10_5, df10_6, df10_7, df10_8, df10_9, df10_10)){
  assign(paste("p10",i,sep="_"), compound_return(OGD(base = df, eta=1/1000), df))
  i <- i+1
}





compound_return(portfolio5, df5)

#########################################################################
## get index
getSymbols('^DJI', src = 'yahoo', from = '2007-01-03', to = '2017-11-30')
DJIA_returns = as.numeric(cumprod(weeklyReturn(DJI) + 1))

getSymbols('^GSPC', src = 'yahoo', from = '2007-01-03', to = '2017-11-30')
SP = weeklyReturn(GSPC)
SP_returns = as.numeric(cumprod(weeklyReturn(GSPC) + 1))

SP_day = dailyReturn(GSPC)
SP_returnsday = as.numeric(cumprod(dailyReturn(GSPC) + 1))

SP_month = monthlyReturn(GSPC)
SP_returnsmonth = as.numeric(cumprod(monthlyReturn(GSPC) + 1))
#########################################################################
########################## Simple equal weight portfolio
portfolio_ew <- apply(df_returns, 1, mean)
portfolio_ew_value <- cumprod(1+portfolio_ew)

portfolio_ewd <- apply(day_returns, 1, mean)
portfolio_ew_valued <- cumprod(1+portfolio_ewd)

portfolio_ewm <- apply(month_returns, 1, mean)
portfolio_ew_valuem <- cumprod(1+portfolio_ewm)

plot(portfolio_ew_value, type="l")
plot(portfolio_ew, type="l")
# APY
portfolio_ew_APY <- ((cumprod(1+portfolio_ew)[n])^(1/n))^52-1
portfolio_ew_APYd <- ((cumprod(1+portfolio_ewd)[length(portfoliod_value$return)])^(1/length(portfoliod_value$return)))^252-1
portfolio_ew_APYm <- ((cumprod(1+portfolio_ewm)[length(portfoliom_value$return)])^(1/length(portfoliom_value$return)))^12-1

# change weekly return into APY
#portfolio_APY1 <- (1+portfolio_value1$return)^52-1
#portfolio_APY2 <- (1+portfolio_value2$return)^52-1
portfolio_APY3 <- ((cumprod(1+portfolio_value3$return)[n])^(1/n))^52-1
#portfolio_APY4 <- (1+portfolio_value4$return)^52-1

portfolio_APY5 <- ((cumprod(1+portfolio_value5$return)[564])^(1/564))^52-1
portfolio_APY10 <- ((cumprod(1+portfolio_value10$return)[564])^(1/564))^52-1
portfolio_APY15 <- ((cumprod(1+portfolio_value15$return)[564])^(1/564))^52-1
portfolio_APY20 <- ((cumprod(1+portfolio_value20$return)[564])^(1/564))^52-1
portfolio_APY25 <- ((cumprod(1+portfolio_value25$return)[564])^(1/564))^52-1
portfolio_APY30 <- ((cumprod(1+portfolio_value30$return)[564])^(1/564))^52-1
portfolio_APY35 <- ((cumprod(1+portfolio_value35$return)[564])^(1/564))^52-1
portfolio_APY40 <- ((cumprod(1+portfolio_value40$return)[564])^(1/564))^52-1
portfolio_APY45 <- ((cumprod(1+portfolio_value45$return)[564])^(1/564))^52-1
portfolio_APY50 <- ((cumprod(1+portfolio_value50$return)[564])^(1/564))^52-1



# get SD of returns
#portfolio_sd1 <- sd(portfolio_value1$return)
#portfolio_sd2 <- sd(portfolio_value2$return)
portfolio_sd3 <- sd(portfolio_value3$return)
#portfolio_sd4 <- sd(portfolio_value4$return)

portfolio_sd5 <- sd(portfolio_value5$return)
portfolio_sd10 <- sd(portfolio_value10$return)
portfolio_sd15 <- sd(portfolio_value15$return)
portfolio_sd20 <- sd(portfolio_value20$return)
portfolio_sd25 <- sd(portfolio_value25$return)
portfolio_sd30 <- sd(portfolio_value30$return)
portfolio_sd35 <- sd(portfolio_value35$return)
portfolio_sd40 <- sd(portfolio_value40$return)
portfolio_sd45 <- sd(portfolio_value45$return)
portfolio_sd50 <- sd(portfolio_value50$return)

## the T bill return from 2007 to 2016, average is about 0.0074
## compute sharpe ratio
#Sharpe_ratio1 <- (portfolio_meanAPY1-0.05028)/portfolio_sd1
#Sharpe_ratio2 <- (portfolio_meanAPY2-0.05028)/portfolio_sd2
Sharpe_ratio3 <- (portfolio_APY3-0.0074)/portfolio_sd3
#Sharpe_ratio4 <- (portfolio_meanAPY4-0.05028)/portfolio_sd4
SP_APY <- ((SP_returns[570])^(1/570))^52-1
SP_APYd <- ((SP_returnsday[length(SP_returnsday)])^(1/length(SP_returnsday)))^252-1
SP_APYm <- ((SP_returnsday[length(SP_returnsmonth)])^(1/length(SP_returnsmonth)))^12-1
Sharpe_ratioSP <- (SP_APY-0.0074)/sd(SP)
Sharpe_ratioSPd <- (SP_APYd-0.0074)/sd(SP_day)
Sharpe_ratioSPm <- (SP_APYm-0.0074)/sd(SP_month)

Sharpe_ratio5 <- (portfolio_APY5-0.0074)/portfolio_sd5
Sharpe_ratio10 <- (portfolio_APY10-0.0074)/portfolio_sd10
Sharpe_ratio15 <- (portfolio_APY15-0.0074)/portfolio_sd15
Sharpe_ratio20 <- (portfolio_APY20-0.0074)/portfolio_sd20
Sharpe_ratio25 <- (portfolio_APY25-0.0074)/portfolio_sd25
Sharpe_ratio30 <- (portfolio_APY30-0.0074)/portfolio_sd30
Sharpe_ratio35 <- (portfolio_APY35-0.0074)/portfolio_sd35
Sharpe_ratio40 <- (portfolio_APY40-0.0074)/portfolio_sd40
Sharpe_ratio45 <- (portfolio_APY45-0.0074)/portfolio_sd45
Sharpe_ratio50 <- (portfolio_APY50-0.0074)/portfolio_sd50

Sharpe_ratioEW <- (portfolio_ew_APY-0.0074)/sd(portfolio_ew)
Sharpe_ratioEWd <- (portfolio_ew_APYd-0.0074)/sd(portfolio_ewd)
Sharpe_ratioEWm <- (portfolio_ew_APYm-0.0074)/sd(portfolio_ewm)
######################### Sharpe ratio table
SR_table <- data.frame(c(Sharpe_ratioSP,Sharpe_ratioEW,Sharpe_ratio3,Sharpe_ratio5,Sharpe_ratio10,Sharpe_ratio15,
                         Sharpe_ratio20, Sharpe_ratio25, Sharpe_ratio30, Sharpe_ratio35,
                         Sharpe_ratio40, Sharpe_ratio45, Sharpe_ratio50))
row.names(SR_table) <- c("SP500 Index","Equal Weight Portfolio","All SP500 Stocks","Size 5","Size 10","Size 15", "Size 20",
                         "Size 25", "Size 30", "Size 35", "Size 40", "Size 45", "Size 50")
colnames(SR_table) <- "Sharpe ratio"
library(gridExtra)
grid.table(SR_table)


####################### Fig 2 mean APY of different size
mean_APY_size <- c(portfolio_APY5, portfolio_APY10, portfolio_APY15,
                   portfolio_APY20, portfolio_APY25, portfolio_APY30,
                   portfolio_APY35, portfolio_APY40, portfolio_APY45,
                   portfolio_APY50)
x <- seq(5, 50, 5)
plot(x,mean_APY_size, type = "b",lty=2,xlab="Number of Stocks", ylab="mean APY")
lines(x, mean_APY_size)

dev.off()
###################### Fig 5 mean APY of 1 2 3 4
#mean_APY_1234 <- c(portfolio_meanAPY1,portfolio_meanAPY2,portfolio_meanAPY3,portfolio_meanAPY4)
#barplot(mean_APY_1234, names.arg=c("date_return,eta=0.01","date_return,eta=0.001","df_return,eta=0.01","df_return,eta=0.001")
#        ,ylab = "mean APY", xlab = "Portfolio", border = "coral3", col="bisque3")



plot(x, size, type='l', lty=2, xlab="Number of Stocks", "")
par(mfrow=c(1,2))
plot(portfolio_value3$return, type='l')
plot(portfolio_value3$value, type='l')



# Individual stocks
stocks_returns = apply(data_returns + 1, 2, cumprod)

# Our portfolios
portfolio_value1 <- compound_return(portfolio1, data_returns)
portfolio_value2 <- compound_return(portfolio2, data_returns)

plot(SP_returns, type="l")

par(mfrow=c(1,2))
plot(portfolio_value3$return, type='l')

################ Fig 6 trading days
jpeg('p6.jpg')
plot(portfolio_value3$value, type='l', ylab="Portfolio Value", xlab="Trading weeks", ylim=c(0.4, 3))
lines(portfolio_value10$value, type='l',  col="bisque3")
lines(SP_returns, type='l', col="brown4")
lines(portfolio_value5$value, type='l', col="coral3")
lines(portfolio_value50$value, type='l', col="chocolate2")
legend("topleft",legend = c("All SP500 Stock", "Portfolio with size 5"
                  ,"SP500 Index Return","Portfolio with size 5", "Portfolio of size 50"), col = c("black","bisque3",
                  "brown4","coral3", "chocolate2"), lty = 1, cex=0.8)
dev.off()


risk_analysis <- function(return) {
  
  report = matrix(0, ncol = 2, nrow = 2)
  
  report[1,] = quantile(return, probs = c(0.01, 0.05))
  report[2,] = c(mean(return[which(return = report[1,1])]), mean(return[which(return = report[1,2])]))
  
  rownames(report) = c('VaR', 'CVaR')
  colnames(report) = c('1%', '5%')
  
  return(round(report, 3))
}

report1 <- risk_analysis(portfolio_value1$return)
report2 <- risk_analysis(portfolio_value2$return)
report_DJIA <- risk_analysis(weeklyReturn(DJI))



###################### Portfolio by variance
Var <- apply(day_returns, 2, var)
Var <- sort(Var)
Var_low <- Var[1:200]
Var_high <- Var[201:400]

ind_low <- sample(Var_low, 20)
ind_high <- sample(Var_high, 20)

names_low <- names(Var_low)
names_high <- names(Var_high)

return_low <- day_returns[ind_low]
return_high <- day_returns[ind_high]

return_low <- subset(day_returns, select=names_low)
return_high <- subset(day_returns, select=names_high)

portfolio_low <- OGD(base = return_low, eta = 1/1000)
portfolio_value_low <- compound_return(portfolio_low, return_low)

portfolio_high <- OGD(base = return_high, eta = 1/1000)
portfolio_value_high <- compound_return(portfolio_high, return_high)

# APY
portfolio_APYlow <- ((cumprod(1+portfolio_value_low$return)[2716])^(1/2716))^252-1
portfolio_APYhigh <- ((cumprod(1+portfolio_value_high$return)[2716])^(1/2718))^252-1


(portfolio_APYhigh-0.0074)/sd(portfolio_value_high$return)
(portfolio_APYlow-0.0074)/sd(portfolio_value_low$return)




### final sample
#[1] 185 426 379 341 155 406 253 148 429 269 208 140 427 311 225 251
#[17] 325  31 337 384 175  36 386 133   3  39 260 203  50 259 404 286
#[33] 275 115 391 282 264  82 303 120 322  44 183 400 368 257 139 310
#[49] 394 288













