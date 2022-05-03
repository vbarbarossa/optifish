library(dplyr); library(tidyr); library(lubridate)

################### import data 
nyse <- data.table::fread('data/archive/prices.csv')

nyse <- nyse %>% 
  mutate(date = ymd(date)) %>% 
  filter(year(date) == 2015,
         month(date) %in% c(1:3))

head(nyse)

securities <- data.table::fread("data/archive/securities.csv")
securities <- securities %>% 
  select(`Ticker symbol`, Security) %>% 
  rename(symbol = `Ticker symbol`,
         name = Security)

nyse <- nyse %>% 
  left_join(securities, by = c("symbol")) %>% 
  select(date, symbol, name, everything())

head(nyse)

set.seed(13)
selected_stock <- sample(nyse$symbol, 30)

nyse <- nyse %>% 
  filter(symbol %in% selected_stock)
head(nyse)

#################### calculate returns
nyse <- nyse %>%   
  rename(price = close) %>% 
  select(date, symbol, name, price) %>% 
  group_by(symbol, name) %>%
  mutate(price_prev = lag(price),
         returns = (price - price_prev)/price_prev) %>% 
  slice(-1) %>% 
  ungroup()

head(nyse)

mean_stock <- nyse %>% 
  group_by(symbol) %>% 
  summarise(mean = mean(returns)) %>% 
  arrange(desc(mean))

head(mean_stock)

rf <- 0.04/100

########### covariance matrix between portfolios
nyse_wide <- nyse %>%
  pivot_wider(id_cols = c(date, symbol), names_from = symbol, values_from = returns) %>% 
  select(-date)

# Create Excess Return
for (i in 1:n_distinct(nyse$symbol)) {
  nyse_wide[,i]<- nyse_wide[,i] - as.numeric(mean_stock[i,2])
}

head(nyse_wide)

nyse_cov <- cov(x = nyse_wide)
############# define fitness function
fitness <- function(x){
  # Assign weight for each stocks
  weight_stock <- x
  
  # Calculate the total returns
  f1 <- numeric()
  for (i in 1:n_distinct(nyse$symbol)) {
    f1[i] <- weight_stock[i]*mean_stock$mean[i]
  }
  mean_return <- sum(f1) - 1e9 * (round(sum(weight_stock),10)-1)^2 
  
  # Calculate the total risk
  f2 <- numeric()
  for (i in 1:n_distinct(nyse$symbol)) {
    f3 <- numeric()
    
    for (j in 1:n_distinct(nyse$symbol)) {
      f3[j] <- weight_stock[i]*weight_stock[j]*nyse_cov[i,j]
    }
    f2[i] <- sum(f3)
  }
  risk <- sum(f2) + 1e9 * (round(sum(weight_stock),10)-1)^2
  
  # Calculate the number of assets
  card <- length(weight_stock[weight_stock > 0]) 
  
  return(c(-mean_return, risk, card))
}


set.seed(123)

n_asset <- n_distinct(nyse$symbol)

library(nsga2R)

finance_optim <- nsga2R(fn = fitness, varNo = n_asset, objDim = 3, generations = 1000,
                        mprob = 0.2, popSize = 200, cprob = 0.8,
                        lowerBounds = rep(0, n_asset), upperBounds = rep(1, n_asset))



########## our dataset
dams <- data.frame(
  id = 1:30,
  IC = rnorm(30,1000,100)
) %>%
  mutate(vol = IC*3.1*rnorm(30,3,0.5))

fitness <- function(x){
  
  decision <- round(x,0)
  
  totIC <- sum(dams$IC*decision)
  totVol <- sum(dams$vol*decision)
  
  # here plug in connectivity
  # totCI <- a number based on the config of the dams
  
  return(c(-totIC,totVol))
  
}

n <- nrow(dams)

# optim <- nsga2R(fn = fitness, varNo = n, objDim = 2, generations = 100,
#                 mprob = 0.2, popSize = 20, cprob = 0.8,
#                 lowerBounds = rep(0, n), upperBounds = rep(1, n))
# 
# plot(optim$objectives)

# <<<<<<<<<<<<<<<<<< much faster with mco, no printing time
library(mco)
optim <- nsga2(fn = fitness, idim = n, odim = 2, generations = 1000,
                mprob = 0.2, popsize = 100, cprob = 0.8,
                lower.bounds = rep(0, n), upper.bounds = rep(1, n))
plot(optim)

