### Packages ####
install.packages("pacman")
pacman::p_load(nnet,gmodels,psych,MASS,VGAM,Hmisc,pROC,clusterPower,ordinal, DescTools,R.utils,mvtnorm,robustbase,rms)

install.packages("mice")
library(mice)
library(ggplot2)
library(tidyr)

#Set working directory

setwd("/Users/vc/Documents/Studie/MSc Epi/Research/Code")

####################### DATA GENERATING MECHANISM ##############################

#### Seed ####
set.seed(1234567)

#### 1. Multinominal Logistic Regression (MLR) ####

simulate_data_k3_mlr <- function( n, # sample size
                                  mu_x, # mean x
                                  cprev, # prevalence categories in c
                                  beta_int,
                                  beta_c2,
                                  beta_c3,
                                  beta_x
){
  # sample c multinomial
  c <- sample( 1:3,
               n,
               replace = TRUE,
               prob = cprev)
  
  # sample x with dependency with c
  x <- rnorm( n,
              mean = 0,
              sd = 1) + ( (c == 1) * mu_x[1] +
                            ( c == 2) * mu_x[2] +
                            ( c == 3) * mu_x[3])
  
  # create data frame
  datax <- as.data.frame( cbind( x,c)) #data frame
  names( datax) <- c( "x","c") #namen van kolomen
  datax$c <- as.factor( datax$c) #as factor nodig gehad
  
  
  # store probabilities of c given x --> misschien niet relevant in onze studie?
  pdf1 <- dnorm( x = datax[ , "x"],
                 mean = mu_x[1],
                 sd = 1)
  pdf2 <- dnorm( x = datax[ , "x"],
                 mean = mu_x[2],
                 sd = 1)
  pdf3 <- dnorm( x = datax[ , "x"],
                 mean = mu_x[3],
                 sd = 1)
  
  #true probability for categories
  datax$ptrue1 <- ( cprev[1]*pdf1) / ( cprev[1] * pdf1 + cprev[2] * pdf2 + cprev[3] * pdf3)
  datax$ptrue2 <- ( cprev[2]*pdf2) / ( cprev[1] * pdf1 + cprev[2] * pdf2 + cprev[3] * pdf3)
  datax$ptrue3 <- ( cprev[3]*pdf3) / ( cprev[1] * pdf1 + cprev[2] * pdf2 + cprev[3] * pdf3)
  
  # generating y
  datax$y <- rnorm( n = n,
                    mean = beta_int +
                      beta_c2 * (datax$c %in% "2") + #2 and 3 as factors
                      beta_c3 * (datax$c %in% "3") +
                      beta_x * datax$x,
                    sd = 1)
  return(datax)
}

# Paramters for mlr

data_mlr.n <- simulate_data_k3_mlr( n = 100000,
                                  mu_x = c(0.0, 0.40, 0.80), 
                                  cprev = c(0.50, 0.35, 0.15),
                                  beta_int = 1,
                                  beta_c2 = 0.5,
                                  beta_c3 = 2,
                                  beta_x = 1.5)

# Save and load file
save(data_mlr.n, file = "data_mlr.n.RData") #full_data is before ampute

load("data_mlr.n.RData") #load full mlr data

### FREQUENCY GRAPH ###
# Convert 'c' column to factor
data_mlr.n$c <- as.factor(data_mlr.n$c)

# Create a data frame with the counts of each value in 'c'
count_data.n <- data.frame(table(data_mlr.n$c))

# Change the labels of 'Var1' to "C1", "C2", and "C3"
count_data.n$Var1 <- factor(count_data.n$Var1, levels = 1:3, labels = c("C1", "C2", "C3"))

# Create the histogram with text labels for frequency
ggplot(count_data.n, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "white") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3, color = "black") +
  labs(title = "Frequency of C from MLR",
       x = "Categories",
       y = "Frequency")



# fitting of mlr model to check whether DGM works
VGAM::vglm( c ~ x,
            family = VGAM::multinomial( refLevel = "1"),
            data = data_mlr,n) 

#glm fit for our data set
glm( y ~ c + x, data = data_mlr.n)


#Distribution of C check
frequency_data_mlr.n <- table(data_mlr.n$c)
percentage_mlr.n <- prop.table(frequency_data_mlr.n) * 100
percentage_mlr.n

################################################################################
################################################################################

#### Cumulative Logit ####

simulate_data_k3_cl <- function( n, # sample size
                                 mu_x, # gemiddelde x, 1 waarde
                                 cprev, # prevalentie van categorien in c
                                 alpha_c2,
                                 alpha_c3,
                                 gamma_c2,
                                 gamma_c3,
                                 beta_int,
                                 beta_c2,
                                 beta_c3,
                                 beta_x
){
  x <- rnorm( n, mu_x, sd = 1.21) #simulatie x
  lp_c2 <- alpha_c2 - gamma_c2*x #linear predictor categorie 2
  lp_c3 <- alpha_c3 - gamma_c3*x
  
  ptrue1 <- exp( lp_c2)/( 1 + exp( lp_c2)) #probability categorie
  ptrue2 <- exp( lp_c3)/( 1 + exp( lp_c3)) - ptrue1
  ptrue3 <- 1 - ptrue1 - ptrue2
  
  c <- matrix( data = NA,
               nrow = length(lp_c2),
               ncol = 1)
  
  for (i in 1:length(lp_c2)){
    c[ i] = which( grepl( 1,
                          rmultinom( 1, 1, c( ptrue1[i], ptrue2[i], ptrue3[i]))))
  }
  
  datax <- as.data.frame(cbind(x, c, ptrue1, ptrue2, ptrue3))
  names( datax)=c("x","c","ptrue1","ptrue2","ptrue3")
  datax$c <- as.factor( datax$c)
  
  # nu nog y toevoegen #betas hier
  datax$y <- beta_int +
    beta_c2 * (datax$c %in% "2") + 
    beta_c3 * (datax$c %in% "3") +
    beta_x * datax$x + rnorm( n = n,
                    mean = 0,
                    sd = 1)
  
  return(datax)
}

data_cl <- simulate_data_k3_cl( n = 100000,
                                mu_x = 0.4,
                                alpha_c2 = -0.18,
                                alpha_c3 = 1.55,
                                gamma_c2 = -0.55,
                                gamma_c3 = -0.55,
                                beta_int = 1,
                                beta_c2 = 0.5,
                                beta_c3 = 2,
                                beta_x = 1.5)



# Save and load file
save(data_cl, file = "data_cl.RData") #data before ampute

load("full_data_cl.RData") #load full cl data


### FREQUENCY GRAPH ###
# Convert 'c' column to factor
data_cl$c <- as.factor(data_cl$c)

# Create a data frame with the counts of each value in 'c'
count_data <- data.frame(table(data_cl$c))

# Change the labels of 'Var1' to "C1", "C2", and "C3"
count_data$Var1 <- factor(count_data$Var1, levels = 1:3, labels = c("C1", "C2", "C3"))

# Create the histogram with text labels for frequency
ggplot(count_data, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "white") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3, color = "black") +
  labs(title = "Frequency of C from CL",
       x = "Categories",
       y = "Frequency")



# Fitting for our data set to check correctness of DGM
VGAM::vglm( c ~ x,
            family = "cumulative",
            data = data_cl)

#glm for the data set CL
glm( y ~ c + x, data = data_cl)

vglm( c ~ x, family=multinomial(refLevel = "1"), data=data_cl) 




