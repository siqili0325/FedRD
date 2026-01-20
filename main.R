rm(list = ls())
library(ggplot2)
library(tidyverse)
library(Metrics)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("fun.R")

#----- Set ---------------------------------------------------
beta1.true <- 1
beta2.true <- 0.5
beta3.true <- 0.5
n1 <- 100
n2 <- 100
n3 <- 100
n4 <- 100
n5 <- 100
l1 <- 0.02 
u1 <- 1.28

set.seed(42)
# simulation
repli <- 500
p <- 3
n_site <- 5 
est <- matrix(0, nrow = repli, ncol = p*(4+n_site))
se <- matrix(0, nrow = repli, ncol = p*(4+n_site))
timerec <- matrix(0, nrow = repli, ncol = (4+n_site))

for(i in 1:repli){
  site_list <- c("data1", "data2", "data3", "data4", "data5")
  Lambda0 <- function(t){return(t^2)}
  data1 <- get_data(n=n1, l1=0.02, u1=1.28)
  Lambda0 <- function(t){return(t^3)}
  data2 <- get_data(n=n2, l1=0.02, u1=1.28)
  Lambda0 <- function(t){return(t^4)}
  data3 <- get_data(n=n3, l1=0.02, u1=1.28)
  Lambda0 <- function(t){return(log(1+t)+t^3)}
  data4 <- get_data(n=n4, l1=0.02, u1=1.28)
  Lambda0 <- function(t){return(log(1+t)+t^4)}
  data5 <- get_data(n=n5, l1=0.02, u1=1.28)
  
  ##===========
  ## Fit model
  ##===========
  ## ========== Federated estimator -unstratified model ============================================================
  ## First, calculate summary estimator X_bar
  time0 <- Sys.time()
  sumY <- c(data1$label_time, data2$label_time, data3$label_time, data4$label_time, data5$label_time)
  y_sort<- sort(sumY)
  X1_bar <- Xk_bar(yobs = data1$label_time, X=subset(data1, select = -c(label_time, delta)), y_sort )
  X2_bar <- Xk_bar(yobs = data2$label_time, X=subset(data2, select = -c(label_time, delta)), y_sort )
  X3_bar <- Xk_bar(yobs = data3$label_time, X=subset(data3, select = -c(label_time, delta)), y_sort )
  X4_bar <- Xk_bar(yobs = data4$label_time, X=subset(data4, select = -c(label_time, delta)), y_sort )
  X5_bar <- Xk_bar(yobs = data5$label_time, X=subset(data5, select = -c(label_time, delta)), y_sort )
  X_bar_fenzi <- (X1_bar$Xk_bar_fenzi + X2_bar$Xk_bar_fenzi + X3_bar$Xk_bar_fenzi + X4_bar$Xk_bar_fenzi + X5_bar$Xk_bar_fenzi) #N*p
  X_bar_fenmu <- (X1_bar$Xk_bar_fenmu + X2_bar$Xk_bar_fenmu + X3_bar$Xk_bar_fenmu + X4_bar$Xk_bar_fenmu + X5_bar$Xk_bar_fenmu) #N*1
  X_bar <- X_bar_fenzi/X_bar_fenmu #N*p
  
  ## Then, calculate V for each site
  Vu1 <- V_k_un(yobs = data1$label_time, delta = data1$delta, X=subset(data1, select = -c(label_time, delta)), y_sort, X_bar)
  Vu2 <- V_k_un(yobs = data2$label_time, delta = data2$delta, X=subset(data2, select = -c(label_time, delta)), y_sort, X_bar)
  Vu3 <- V_k_un(yobs = data3$label_time, delta = data3$delta, X=subset(data3, select = -c(label_time, delta)), y_sort, X_bar)
  Vu4 <- V_k_un(yobs = data4$label_time, delta = data4$delta, X=subset(data4, select = -c(label_time, delta)), y_sort, X_bar)
  Vu5 <- V_k_un(yobs = data5$label_time, delta = data5$delta, X=subset(data5, select = -c(label_time, delta)), y_sort, X_bar)
  ## --------- Save the information ---------------------------------------------
  write_site_V(sitename = "data1", infor = Vu1)
  write_site_V(sitename = "data2", infor = Vu2)
  write_site_V(sitename = "data3", infor = Vu3)
  write_site_V(sitename = "data4", infor = Vu4)
  write_site_V(sitename = "data5", infor = Vu5)
  
  ## --------- Read information and conduct Federate estimate -------------------
  Fed_unest <- get_Fedbeta(site_list=c("data1", "data2", "data3", "data4", "data5"), var_p = p )
  est[i, 1:p] <- Fed_unest$Est
  se[i, 1:p] <- Fed_unest$SE
  time1 <- Sys.time()
  timerec[i,1] <- time1-time0
  
  ## ========== Federated estimator -Stratified model ============================================================  
  time0 <- Sys.time()
  ## --------- Calculate information of each site -------------------------------
  # Note: data and X must be data.frame
  V1 <- V_k(yobs = data1$label_time, delta = data1$delta,  X=subset(data1, select = -c(label_time, delta))  )
  V2 <- V_k(yobs = data2$label_time, delta = data2$delta,  X=subset(data2, select = -c(label_time, delta)) )
  V3 <- V_k(yobs = data3$label_time, delta = data3$delta,  X=subset(data3, select = -c(label_time, delta)) )
  V4 <- V_k(yobs = data4$label_time, delta = data4$delta,  X=subset(data4, select = -c(label_time, delta)) )
  V5 <- V_k(yobs = data5$label_time, delta = data5$delta,  X=subset(data5, select = -c(label_time, delta)) )
  
  ## --------- Save the information ---------------------------------------------
  write_site_V(sitename = "data1", infor = V1)
  write_site_V(sitename = "data2", infor = V2)
  write_site_V(sitename = "data3", infor = V3)
  write_site_V(sitename = "data4", infor = V4)
  write_site_V(sitename = "data5", infor = V5)
  
  ## --------- Read information and conduct Federate estimate -------------------
  Fed_est <- get_Fedbeta(site_list=c("data1", "data2", "data3", "data4", "data5"), var_p = p )
  
  est[i, (p+1):(2*p)] <- Fed_est$Est
  se[i, (p+1):(2*p)] <- Fed_est$SE
  time1 <- Sys.time()
  timerec[i,2] <- time1-time0
  
  ## ========== Local estimator =================================================================
  ## --------- local site1 for beta-------------------------------------------------------------
  time0 <-time00<- Sys.time()
  site1_est <- AR.fit(data1$label_time, data1$delta, X=subset(data1, select = -c(label_time, delta)) )$coef
  est[i, (2*p+1):(3*p)] <- site1_est$Est
  se[i, (2*p+1):(3*p)] <- site1_est$SE
  time1 <- Sys.time()
  timerec[i,3] <- time1-time0
  
  ## --------- local site2-------------------------------------------------------------
  time0 <- Sys.time()
  site2_est <- AR.fit(data2$label_time, data2$delta, X=subset(data2, select = -c(label_time, delta)) )$coef
  est[i, (3*p+1):(4*p)] <- site2_est$Est
  se[i, (3*p+1):(4*p)] <- site2_est$SE
  time1 <- Sys.time()
  timerec[i,4] <- time1-time0
  
  ## --------- local site3-------------------------------------------------------------
  time0 <- Sys.time()
  site3_est <- AR.fit(data3$label_time, data3$delta, X=subset(data3, select = -c(label_time, delta)) )$coef
  est[i, (4*p+1):(5*p)] <- site3_est$Est
  se[i, (4*p+1):(5*p)] <- site3_est$SE
  time1 <- Sys.time()
  timerec[i,5] <- time1-time0
  
  ## --------- local site4-------------------------------------------------------------
  time0 <- Sys.time()
  site4_est <- AR.fit(data4$label_time, data4$delta, X=subset(data4, select = -c(label_time, delta)) )$coef
  est[i, (5*p+1):(6*p)] <- site4_est$Est
  se[i, (5*p+1):(6*p)] <- site4_est$SE
  time1 <- Sys.time()
  timerec[i,6] <- time1-time0
  
  ## --------- local site5-------------------------------------------------------------
  time0 <- Sys.time()
  site5_est <- AR.fit(data5$label_time, data5$delta, X=subset(data5, select = -c(label_time, delta)) )$coef
  est[i, (6*p+1):(7*p)] <- site5_est$Est
  se[i, (6*p+1):(7*p)] <- site5_est$SE
  time1 <- Sys.time()
  timerec[i,7] <- time1-time0
  
  ## ========= Ens estimator ==================================================================
  # est[i, (7*p+1):(8*p)] <- (site1_est$Est + site2_est$Est + site3_est$Est + site4_est$Est + site5_est$Est)/length(site_list)
  est[i, (7*p+1):(8*p)] <- (site1_est$Est/site1_est$SE + site2_est$Est/site2_est$SE + site3_est$Est/site3_est$SE
                            + site4_est$Est/site4_est$SE + site5_est$Est/site5_est$SE)/(1/site1_est$SE 
                                                                                        + 1/site2_est$SE 
                                                                                        + 1/site3_est$SE 
                                                                                        + 1/site4_est$SE 
                                                                                        + 1/site5_est$SE)
   se[i, (7*p+1):(8*p)] <- (site1_est$SE + site2_est$SE + site3_est$SE + site4_est$SE + site5_est$SE)/(length(site_list)*2)
  time1 <- Sys.time()
  timerec[i,8] <- time1-time00
  ## ========= Pooled estimator ==================================================================
  time0 <- Sys.time()
  data_sum <- rbind(data1, data2, data3, data4, data5)
  Pool_est <- AR.fit(data_sum$label_time, data_sum$delta, X= subset(data_sum, select=-c(label_time, delta)) )$coef
  est[i, (8*p+1):(9*p)] <- Pool_est$Est 
  se[i, (8*p+1):(9*p)] <- Pool_est$SE
  time1 <- Sys.time()
  timerec[i,9] <- time1-time0
  
  print(i)
}

## =======================================
## Summary
## =======================================
est_ave <- apply(est, 2, mean)
bet_true <- rep(c(beta1.true, beta2.true, beta3.true), 4+n_site)
bias <- est_ave - bet_true
ssd <- apply(est, 2, sd)
see <- apply(se, 2, mean)
timemean <- apply(timerec, 2, mean)
timerun <- rep(timemean, each=3)

cp <- mse<- rep(0,p*(4+n_site))
for (ii in 1:(p*(4+n_site)) ){
  mse[ii] <- mse(est[,ii], bet_true[ii])
  if (ii %in% c(1,4,7,10,13,16,19, 22,25)) {
    cp[ii] <- mean(as.matrix(beta1.true <= est[,ii] + 1.96*se[,ii] & 
                               beta1.true >= est[,ii] - 1.96*se[,ii]))
  }else if (ii %in% c(2, 5, 8, 11, 14, 17, 20,23, 26)) {
    cp[ii] <- mean(as.matrix(beta2.true <= est[,ii] + 1.96*se[,ii] & 
                               beta2.true >= est[,ii] - 1.96*se[,ii]))
  }
  else {
    cp[ii] <- mean(as.matrix(beta3.true <= est[,ii] + 1.96*se[,ii] & 
                               beta3.true >= est[,ii] - 1.96*se[,ii]))
  }
}

results <- cbind(bias,ssd,see,cp,mse,timerun)
rownames(results)=c("FedRD(U)_b1","FedRD(U)_b2","FedRD(U)_b3", 
                    "FedRD(S)_b1","FedRD(S)_b2","FedRD(S)_b3", 
                    "Site1_b1", "Site1_b2","Site1_b3",
                    "Site2_b1","Site2_b2","Site2_b3", 
                    "Site3_b1","Site3_b2","Site3_b3",
                    "Site4_b1","Site4_b2","Site4_b3",
                    "Site5_b1","Site5_b2","Site5_b3",
                    "Ens_b1", "Ens_b2", "Ens_b3",
                    "Pool_b1", "Pool_b2", "Pool_b3")
print(results)

## ---- Plot ------------------------------------

results_est <- data.frame(
  Estimator = as.vector(est),
  Methods = c(rep("FedRD-U", 3*repli),rep("FedRD-S", 3*repli), rep("Local1", 3*repli),  rep("Local2", 3*repli), rep("Local3", 3*repli),
              rep("Local4", 3*repli), rep("Local5", 3*repli),
              rep("Meta",3*repli), rep("Pooled", 3*repli)),
  Coef = rep(c(rep("beta1",repli), rep("beta2",repli), rep("beta3",repli)), 4+n_site))

results_est$Methods <- factor(results_est$Methods, levels = c('Pooled', 'FedRD-U', 'FedRD-S', 'Meta', 'Local1', 'Local2', 'Local3', 'Local4', 'Local5'))
mycolors <- c("#EBEBEB","#F1BFB5", "#FBE7C0","#CCECE6","#CAD4E7","#E8D8EA", "#beaed4", "#899CCB",  "#899CCB", "#788AB6")

library(ggplot2)

lines_data <- data.frame(y = c(1, 0.5, 0.5),
                         Coef = c("beta1", "beta2", "beta3")) %>% 
  mutate(Coef = factor(Coef,
                       levels = c("beta1", "beta2", "beta3"),
                       labels = c(expression(beta[1]), 
                                  expression(beta[2]),
                                  expression(beta[3]))
  ))

boxplot <- results_est %>%
  mutate(Coef = factor(Coef,
                       levels = c("beta1", "beta2", "beta3"),
                       labels = c(expression(beta[1]), 
                                  expression(beta[2]),
                                  expression(beta[3]))
  )) %>% ggplot( aes(Methods, Estimator, fill = Methods)) + 
  geom_boxplot() + 
  facet_wrap(~ Coef, ncol = 3, labeller = label_parsed) + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 1.5, color = "blue")+
  geom_hline(data = lines_data, aes(yintercept = y),linetype = "dashed", color = "red")+
  scale_y_continuous(breaks = seq(-4, 4, by = 0.5)) +
  labs(x = "") + 
  scale_fill_manual(values = mycolors) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 70, vjust = 0.1, hjust = 0.1),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 16),      # Legend title "Method"
        legend.text = element_text(size = 16),       # Legend labels
        legend.key.size = unit(1, "cm")
  ) 

ggsave("FedRD.pdf", plot=boxplot, width= 20, height = 6.5)

