### Two-way fixed effects models in air pollution epidemiology
### a proposed framework for model specifications and recommendations for future studies
### Yiqun Ma, UCSD, 08/26/2025

library(dplyr)
library(fixest)
library(splines)

### (1) Main Analysis
## load data
data.all <- readRDS("CA_weekly_ZCTA_data_for_analysis.RDS")
FE.table <- read.csv("FE_combination.csv")

## prepare models
FEs <- FE.table$combination

table <- data.frame(model = c(paste0("Model ",1:length(FEs))),
                    category = FE.table$category,
                    FE = FEs)

formulas <- c(paste0("respiratory ~ PM25_24h + ns(Tmean, df=5)|",
                     FEs))

## run each model
for (i in 1:length(formulas)){
  print(i)
  data <- data.all
  
  formula.i <- formulas[i]
  
  model = feglm(as.formula(formula.i), 
                data = data,
                offset = log(data$pop.tot),
                vcov = "iid",
                family = "quasipoisson")
  
  coef <- model$coefficients[1]
  SE <- model$se[1]
  P <- summary(model)$coeftable[1,4]
  FE.size <- sum(model$fixef_sizes)
  
  table[i,"coef"] <- coef
  table[i,"SE"] <- SE
  table[i,"P"] <- P
  table[i,"FE.size"] <- FE.size
  
}

table <- table %>%
  mutate(pct = (exp(coef) - 1) * 100,
         pct_lower = (exp(coef - 1.96*SE) - 1) * 100,
         pct_upper = (exp(coef + 1.96*SE) - 1) * 100)


### (2) Permutation tests
## Temporal permutation
n <- 1000
for (i in 1:length(formulas)){
  print(i)
  formula.i <- formulas[i]
  model.name <- table[i,"model"]

  simulation.coefs <- foreach(
    j = 1:n, 
    .combine = 'c',
    .packages = c("fixest","dplyr","splines")
  ) %dopar% {
    data <- data.all
    
    data <- data %>%
      group_by(ZCTA) %>%
      mutate(PM25_24h = sample(PM25_24h))
    
    model = feglm(as.formula(formula.i), 
                  data = data,
                  offset = log(data$pop.tot),
                  vcov = "iid",
                  family = "quasipoisson")
    
    model$coefficients[1]
  }
  
  temporal.table <- data.frame(n=1:n,
                              coef = simulation.coefs)
  
  saveRDS(spatial.table, paste0("Randomization_temporal_",model.name,".RDS"))
  
}


## Spatial permutation
n <- 1000
for (i in 1:length(formulas)){
  print(i)
  formula.i <- formulas[i]
  model.name <- table[i,"model"]

  simulation.coefs <- foreach(
    j = 1:n, 
    .combine = 'c',
    .packages = c("fixest","dplyr","splines")
  ) %dopar% {
    data <- data.all
    
    data <- data %>%
      group_by(week) %>%
      mutate(PM25_24h = sample(PM25_24h))
    
    model = feglm(as.formula(formula.i), 
                  data = data,
                  offset = log(data$pop.tot),
                  vcov = "iid",
                  family = "quasipoisson")
    
    model$coefficients[1]
  }
  
  spatial.table <- data.frame(n=1:n,
                              coef = simulation.coefs)
  
  saveRDS(spatial.table, paste0("Randomization_spatial_",model.name,".RDS"))
  
}


### (3) Different SE calculations
VCOVs <- c( "iid", "hetero","cluster by ZCTA", "cluster by county", "twoway")

table.SE <- data.frame(VCOV = VCOVs,
                    model = "Model 38",
                    FE = "week + COUNTY^year + ZCTA^month")

### (2) run each model
for (i in 1:length(VCOVs)){
  print(VCOVs[i])
  data <- data.all
  
  if(VCOVs[i] == "cluster by ZCTA"){
    model = feglm(respiratory ~ PM25_24h + ns(Tmean, df=5)| week + COUNTY^year+ZCTA^month, 
                  data = data,
                  offset = log(data$pop.tot),
                  cluster = data$ZCTA,
                  family = "quasipoisson")
  }else if (VCOVs[i] == "cluster by county"){
    model = feglm(respiratory ~ PM25_24h + ns(Tmean, df=5)| week + COUNTY^year+ZCTA^month, 
                  data = data,
                  offset = log(data$pop.tot),
                  cluster = data$COUNTY,
                  family = "quasipoisson")
  }else if (VCOVs[i] == "twoway"){
    model = feglm(respiratory ~ PM25_24h + ns(Tmean, df=5)| week + COUNTY^year+ZCTA^month, 
                  data = data,
                  offset = log(data$pop.tot),
                  cluster = c("week", "ZCTA"),
                  family = "quasipoisson")
  }else{
    model = feglm(respiratory ~ PM25_24h + ns(Tmean, df=5)| week + COUNTY^year+ZCTA^month, 
                  data = data,
                  offset = log(data$pop.tot),
                  vcov = VCOVs[i],
                  family = "quasipoisson")
  }
  
  coef <- model$coefficients[1]
  SE <- model$se[1]
  P <- summary(model)$coeftable.SE[1,4]
  
  table.SE[i,"coef"] <- coef
  table.SE[i,"SE"] <- SE
  table.SE[i,"P"] <- P
}
