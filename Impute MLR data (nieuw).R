############## mlr.n DATA  IMPUTATION ##############
#Packages
install.packages("kableExtra")
install.packages("reshape2")
install.packages("xtable")
install.packages("flextable")
install.packages("tidyverse")
library(mice)
library(ggplot2)
library(ggmice)
library(kableExtra)
library(reshape2)
library(flextable) #package <- voor export naar word
library(tidyverse)
library(stringr) 
library(ggpubr)

#Seed + file locations
set.seed(1234567)

setwd("/Users/vc/Documents/Studie/MSc Epi/Research/Code")

#to test methods: polr, pmm, polyreg

################################################################################
################################################################################
# Data loading, sample size, and subset creation
#Data loading
load("df_missing_mlr.n.RData")

load("data_mlr.n.RData")

frequency_c_mlr.n <- table(data_mlr.n$c)

#Frequency c in amputed full data
frequency_c_full_mlr.n <- table(df_missing_mlr.n$c)


#sample size, for turning the frequencies into %
ss <- 100000

#subsetting data to show only y, x, and c in amputed data set
sub.df_missing_mlr.n <- subset(df_missing_mlr.n, select = c("y", "c", "x"))


################################################################################
################################################################################

#### FULL DATA ####
##### PMM #####
md.pattern(sub.df_missing_mlr.n)

df_missing_mlr.n$c <- as.factor(df_missing_mlr.n$c)

# Perform multiple imputation using the "pmm" method (cystem time is time elapsed)
imputed_data_mlr.n_full_pmm <- mice(sub.df_missing_mlr.n, m = 15, maxit = 20, method = 'pmm')
system.time(imputed_data_mlr.n_full_pmm <- mice(sub.df_missing_mlr.n, m = 15, maxit = 20, method = 'pmm'))

#speed t = 13,5secs

# Save the imputed data set
save(imputed_data_mlr.n_full_pmm, file = "imputed_data_mlr.n_full_pmm.RData")

#convergence check
mlr.n.full.pmm.con.plot <- plot(imputed_data_mlr.n_full_pmm, y = "c", layout = c(1,2))

# Start the analyses
MIpmm.mlr.n.full <- with(data = imputed_data_mlr.n_full_pmm , glm(y ~ x + c))
MIpmm.mlr.n.full

#Rubins rule
pooledEst.pmm.mlr.n.full <- pool(MIpmm.mlr.n.full)

# Create a data frame with the coefficients and standard errors
coeff_se_full_data_mlr.n_pmm <- data.frame(
  Paramters = pooledEst.pmm.mlr.n.full$pooled["term"],
  Estimate = pooledEst.pmm.mlr.n.full$pooled["estimate"],
  Std_Error = sqrt(pooledEst.pmm.mlr.n.full$pooled["t"]) #t = total variance
)

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.pmm.full <- coeff_se_full_data_mlr.n_pmm[coeff_se_full_data_mlr.n_pmm$term %in% c("c"), ]

# View the imputed dataset
completed_data_mlr.n_full_pmm <- complete(imputed_data_mlr.n_full_pmm)

md.pattern(completed_data_mlr.n_full_pmm)

#frquency moet nog

frequency_c_mlr.n_imputed_pmm <- table(completed_data_mlr.n_full_pmm$c)

full_data_imputed_difference_pmm <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_mlr.n/ss)*100,
  Imputed <- (frequency_c_mlr.n_imputed_pmm/ss)*100
)

frequency_c_mlr.n
frequency_c_full_mlr.n
################################################################################
### POLYREG ###
md.pattern(sub.df_missing_mlr.n)

# Convert 'c' to a factor variable 
sub.df_missing_mlr.n$c <- as.factor(sub.df_missing_mlr.n$c)


# Perform multiple imputation using the "polyreg" method
imputed_data_mlr.n_full_polyreg <- mice(sub.df_missing_mlr.n, m = 15, maxit = 20, method = 'polyreg')

system.time(imputed_data_mlr.n_full_polyreg <- mice(sub.df_missing_mlr.n, m = 15, maxit = 20, method = 'polyreg'))

#time proces: t = 136 sec

# Save the data set
save(imputed_data_mlr.n_full_polyreg , file = "imputed_data_mlr.n_full_polyreg .RData")

# View the imputed dataset
completed_data_mlr.n_full_polyreg <- complete(imputed_data_mlr.n_full_polyreg)

md.pattern(completed_data_mlr.n_full_polyreg)

frequency_c_mlr.n_imputed_polyreg <- table(completed_data_mlr.n_full_polyreg$c)
frequency_c_mlr.n_imputed_polyreg


#convergence check
mlr.n.full.polyreg.con.plot <-plot(imputed_data_mlr.n_full_polyreg, y = "c", layout = c(1,2))

# Start the analyses
MIpolyreg.mlr.n.full <- with(data = imputed_data_mlr.n_full_polyreg , glm(y ~ x + c))
MIpolyreg.mlr.n.full

#Rubins rule
pooledEst.polyreg.mlr.n.full <- pool(MIpolyreg.mlr.n.full)

# Create a data frame with the coefficients and standard errors
coeff_se_full_data_mlr.n_polyreg <- data.frame(
  Coefficients = c("", "", "b2", "b3"),
  #Paramters = pooledEst.polyreg.mlr.n.full$pooled["term"],
  Estimate = pooledEst.polyreg.mlr.n.full$pooled["estimate"],
  Std_Error = sqrt(pooledEst.polyreg.mlr.n.full$pooled["t"]) #t = total variance
)

#change column names
colnames(coeff_se_full_data_mlr.n_polyreg) <- c("Parameter", "Estimate", "Std_Error")

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.polyreg.full <- data.frame(coeff_se_full_data_mlr.n_polyreg [coeff_se_full_data_mlr.n_polyreg $Coefficients  %in% c("b2", "b3"), ])


################################################################################
### POLR ###

md.pattern(df_missing_mlr.n)


# Convert 'c' to a factor variable 
sub.df_missing_mlr.n$c <- as.factor(sub.df_missing_mlr.n$c)

# Perform multiple imputation using the "polr" method
imputed_data_mlr.n_full_polr <- mice(sub.df_missing_mlr.n, m = 15, maxit = 20, method = 'polr')

system.time(syimputed_data_mlr.n_full_polr <- mice(sub.df_missing_mlr.n, m = 15, maxit = 20, method = 'polr'))

#time elapsed polyr: t = 284 seconds
# Save the data set
save(syimputed_data_mlr.n_full_polr , file = "imputed_data_mlr.n_full_polr .RData")

# View the imputed dataset
completed_data_mlr.n_full_polr <- complete(syimputed_data_mlr.n_full_polr)

md.pattern(completed_data_mlr.n_full_polr)

frequency_c_mlr.n_imputed_polr <- table(completed_data_mlr.n_full_polr$c)
frequency_c_mlr.n_imputed_polr


#convergence check
mlr.n.full.polr.con.plot <- plot(syimputed_data_mlr.n_full_polr, y = "c", layout = c(1,2))

# Start the analyses
MIpolr.mlr.n.full <- with(data = syimputed_data_mlr.n_full_polr , glm(y ~ x + c))
MIpolr.mlr.n.full

#Rubins rule
pooledEst.polr.mlr.n.full <- pool(MIpolr.mlr.n.full)


# Create a data frame with the coefficients and standard errors
coeff_se_full_data_mlr.n_polr <- data.frame(
  Paramters = pooledEst.polr.mlr.n.full$pooled["term"],
  Estimate = pooledEst.polr.mlr.n.full$pooled["estimate"],
  Std_Error = sqrt(pooledEst.polr.mlr.n.full$pooled["t"]) #t = total variance
)

#change column names
colnames(coeff_se_full_data_mlr.n_polr) <- c("Parameter", "Estimate", "Std_Error")

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.polr.full <- coeff_se_full_data_mlr.n_polr [coeff_se_full_data_mlr.n_polr $Parameter %in% c("c2", "c3"), ]

# View the imputed dataset
completed_data_mlr.n_full_polr <- complete(syimputed_data_mlr.n_full_polr)

md.pattern(completed_data_mlr.n_full_polr)

#frquency moet nog

frequency_c_mlr.n_imputed_polr <- table(completed_data_mlr.n_full_polr$c)

full_data_imputed_difference_polr <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_mlr.n/ss)*100,
  Imputed <- (frequency_c_mlr.n_imputed_polr/ss)*100
)


################################################################################
#GRAPHS

### FREQUENCY
# Function to create the bar chart
create_bar_chart_mlr.n <- function(data, title) {
  ggplot(data, aes(x = Category)) +
    geom_bar(aes(y = Before, fill = "Before"), stat = "identity", position = "dodge") +
    geom_bar(aes(y = Imputed, fill = "Imputed"), stat = "identity", position = "dodge") +
    labs(title = title,
         x = "Category",
         y = "") +
    scale_fill_manual(values = c("Missing" = "skyblue", "Imputed" = "red")) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Specify percent scale for y-axis
    theme_minimal()
}

#, "Imputed" = "#1f77b4
# Create the plots
frq.plot1.full <- create_bar_chart_mlr.n(full_data_imputed_difference_pmm, "Imputed with PMM")
frq.plot2.full <- create_bar_chart_mlr.n(full_data_imputed_difference_polyreg, "Imputed with POLYREG")
frq.plot3.full <- create_bar_chart_mlr.n(full_data_imputed_difference_polr, "Imputed with POLR")

#grid the frequency plots
ggarrange(frq.plot1.full, frq.plot2.full, frq.plot3.full, ncol=3, common.legend = TRUE, legend="bottom")

ggarrange(mlr.n.full.polr.con.plot, mlr.n.full.polyreg.con.plot, mlr.n.full.pmm.con.plot, 
          common.legend = TRUE, legend="bottom")


# Table coef combined
comb.mlr.n.polyreg.polr.full <- rbind(c2c3.coef.se.polyreg.full, c2c3.coef.se.polr.full)
comb.mlr.n.polyreg.polr.full
a <-as.data.frame(c2c3.coef.se.polyreg.full)
b <- as.data.frame(c2c3.coef.se.polr.full)

flextable(data = list(a,b) , col_keys = "beta coef", "estimate", "std.error")

################################################################################
################################################################################
### mlr.n X and C scenario ###
#load mlr.n x and c scenario file

load("df_missing_mlr.n_x_c.RData")

frequency_xc_mlr.n <- table(df_missing_mlr.n_x_c$c)

#CCA

cca_xc.mlr.S1 <- glm(y ~ x + c, data = df_missing_mlr.n_x_c)

summary(cca_xc.mlr.S1)

#Frequency c in amputed full data
frequency_xc_missing_mlr.n <- (table(df_missing_mlr.n_x_c$c)/ss)*100

##### PMM #####
md.pattern(df_missing_mlr.n_x_c)

# Convert 'c' to a factor variable 
df_missing_mlr.n_x_c$c <- as.factor(df_missing_mlr.n_x_c$c)


# Perform multiple imputation using the "pmm" method
imputed_data_mlr.n_xc_pmm <- mice(df_missing_mlr.n_x_c, m = 15, maxit = 20, method = 'pmm')
system.time(imputed_data_mlr.n_xc_pmm <- mice(df_missing_mlr.n_x_c, m = 15, maxit = 20, method = 'pmm'))
#time elapsed: t = 14sec
# Save the data set
save(imputed_data_mlr.n_xc_pmm, file = "imputed_data_mlr.n_xc_pmm.RData")

#convergence check
mlr.n.xc.pmm.con.plot <- plot(imputed_data_mlr.n_xc_pmm , y = "c", layout = c(1,2))
mlr.n.xc.pmm.con.plot

# Start the analyses
MIpmm.mlr.n.xc <- with(data = imputed_data_mlr.n_xc_pmm , glm(y ~ x + c))
MIpmm.mlr.n.xc

#Rubins rule
pooledEst.pmm.mlr.n.xc <- pool(MIpmm.mlr.n.xc)


# Create a data frame with the coefficients and standard errors
coeff_se_xc_data_mlr.n_pmm <- data.frame(
  Paramters = pooledEst.pmm.mlr.n.xc$pooled["term"],
  Estimate = pooledEst.pmm.mlr.n.xc$pooled["estimate"],
  Std_Error = sqrt(pooledEst.pmm.mlr.n.xc$pooled["t"]) #t = total variance
)

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.pmm.xc <- coeff_se_xc_data_mlr.n_pmm[coeff_se_xc_data_mlr.n_pmm$term %in% c("c2", "c3"), ]

# View the imputed dataset
completed_data_mlr.n_xc_pmm <- complete(imputed_data_mlr.n_xc_pmm)

md.pattern(completed_data_mlr.n_xc_pmm)

#frquency moet nog

frequency_xc_mlr.n_imputed_pmm <- table(completed_data_mlr.n_xc_pmm$c)

xc_data_imputed_difference_pmm <- data.frame(
  Category = c("1", "2", "3"),
  Missing <- ((frequency_c_mlr.n - frequency_c_full_mlr.n)/ss)*100,
  #After <- (frequency_c_mlr.n_imputed_pmm/ss)*100,
  Imputed <- (( frequency_xc_mlr.n_imputed_pmm - frequency_c_mlr.n)/ss)*100
)

#GRAPH: side-by-side bar chart full data vs pmm imputed data

# Create tables for frequency

before.fq.mlr.n <- (frequency_c_mlr.n/ss)*100

fq.full.data.mlr.n <- data.frame(
  Category = c("1", "2", "3"),
  Value = before.fq.mlr.n 
)

fq.pmm.xc.imputed.mlr.n <- (frequency_xc_mlr.n_imputed_pmm/ss)*100

fq.xc.pmm.data.mlr.n  <- data.frame(
  Category = c("1", "2", "3"),
  Value = fq.pmm.xc.imputed.mlr.n
)



fq.xc.amputed.data.mlr.n <- data.frame(
  Category = c("1", "2", "3"),
  Value = frequency_xc_missing_mlr.n
)

# Add an identifier column to distinguish between the tables
fq.full.data.mlr.n$Table <- "1. Full data set"
fq.xc.pmm.data.mlr.n$Table <- "2. Imputed data set"
fq.xc.amputed.data.mlr.n$Table <- "3. Amputed data set"

# Combine the tables into one data frame
combined_data.full.and.xc.pmm <- rbind(fq.full.data.mlr.n, fq.xc.pmm.data.mlr.n, fq.xc.amputed.data.mlr.n)

# Create the side-by-side bar plot
plot1.mlr.n.xc.pmm <- ggplot(combined_data.full.and.xc.pmm, aes(x = Category, y = Value.Freq, fill = Table)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(title = "1. Imputed with pmm method",
       x = "Category",
       y = "Frequency") +
  scale_fill_manual(values = c("1. Full data set"= "skyblue", "2. Imputed data set" = "lightgreen", "3. Amputed data set" = "#1f77b4")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal()

plot1.mlr.n.xc.pmm 




################################################################################
### POLYREG ###

# Convert 'c' to a factor variable 
md.pattern(df_missing_mlr.n_x_c)

df_missing_mlr.n_x_c$c <- as.factor(df_missing_mlr.n_x_c$c)

# Perform multiple imputation using the "polyreg" method
imputed_data_mlr.n_xc_polyreg <- mice(df_missing_mlr.n_x_c, m = 15, maxit = 20, method = 'polyreg')

system.time(imputed_data_mlr.n_xc_polyreg <- mice(df_missing_mlr.n_x_c, m = 15, maxit = 20, method = 'polyreg'))
#time elapsed: t = 137
# Save the data set
save(imputed_data_mlr.n_xc_polyreg, file = "imputed_data_mlr.n_xc_polyreg.RData")

# View the imputed dataset
completed_data_mlr.n_xc_polyreg <- complete(imputed_data_mlr.n_xc_polyreg)

md.pattern(completed_data_mlr.n_xc_polyreg)

frequency_xc_mlr.n_imputed_polyreg <- table(completed_data_mlr.n_xc_polyreg$c)
frequency_xc_mlr.n_imputed_polyreg


#convergence check
mlr.n.xc.polyreg.con.plot <-plot(imputed_data_mlr.n_xc_polyreg, y = "c", layout = c(1,2))

# Start the analyses
MIpolyreg.mlr.n.xc <- with(data = imputed_data_mlr.n_xc_polyreg , glm(y ~ x + c))
MIpolyreg.mlr.n.xc

#Rubins rule
pooledEst.polyreg.mlr.n.xc <- pool(MIpolyreg.mlr.n.xc)

# Create a data frame with the coefficients and standard errors
coeff_se_xc_data_mlr.n_polyreg <- data.frame(
  Paramters = pooledEst.polyreg.mlr.n.xc$pooled["term"],
  Estimate = pooledEst.polyreg.mlr.n.xc$pooled["estimate"],
  Std_Error = sqrt(pooledEst.polyreg.mlr.n.xc$pooled["t"]) #t = total variance
)

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.polyreg.xc <- coeff_se_xc_data_mlr.n_polyreg[coeff_se_xc_data_mlr.n_polyreg$term %in% c("c2", "c3"), ]

# View the imputed dataset
completed_data_mlr.n_xc_polyreg <- complete(imputed_data_mlr.n_xc_polyreg)

md.pattern(completed_data_mlr.n_xc_polyreg)

#frquency moet nog

frequency_xc_mlr.n_imputed_polyreg <- table(completed_data_mlr.n_xc_polyreg$c)

xc_data_imputed_difference_polyreg <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_full_mlr.n/ss)*100,
  After <- (frequency_xc_mlr.n_imputed_polyreg/ss)*100,
  Imputed <- (( frequency_xc_mlr.n_imputed_polyreg - frequency_c_mlr.n)/ss)*100
)


#GRAPH: side-by-side bar chart full data vs polyreg imputed data

# Create tables for frequency

before.fq.mlr.n <- (frequency_c_mlr.n/ss)*100

fq.full.data.mlr.n <- data.frame(
  Category = c("1", "2", "3"),
  Value = before.fq.mlr.n 
)

fq.polyreg.xc.imputed.mlr.n <- (frequency_xc_mlr.n_imputed_polyreg/ss)*100

fq.xc.polyreg.data.mlr.n  <- data.frame(
  Category = c("1", "2", "3"),
  Value = fq.polyreg.xc.imputed.mlr.n
)


# Add an identifier column to distinguish between the tables
fq.full.data.mlr.n$Table <- "1. Full data set"
fq.xc.polyreg.data.mlr.n$Table <- "2. Imputed data set"
fq.xc.amputed.data.mlr.n$Table <- "3. Amputed data set"

# Combine the tables into one data frame
combined_data.full.and.xc.polyreg <- rbind(fq.full.data.mlr.n, fq.xc.polyreg.data.mlr.n, fq.xc.amputed.data.mlr.n)

# Create the side-by-side bar plot
plot1.mlr.n.xc.polyreg <- ggplot(combined_data.full.and.xc.polyreg, aes(x = Category, y = Value.Freq, fill = Table)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(title = "2. Impution with polyreg method",
       x = "Category",
       y = "Frequency") +
  scale_fill_manual(values = c("1. Full data set"= "skyblue", "2. Imputed data set" = "lightgreen", "3. Amputed data set" = "#1f77b4")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal()

plot1.mlr.n.xc.polyreg

################################################################################
### POLR ###

md.pattern(df_missing_mlr.n)

# Convert 'c' to a factor variable 
md.pattern(df_missing_mlr.n_x_c)

df_missing_mlr.n_x_c$c <- as.factor(df_missing_mlr.n_x_c$c)

# Perform multiple imputation using the "polr" method
imputed_data_mlr.n_xc_polr <- mice(df_missing_mlr.n_x_c, m = 15, maxit = 20, method = 'polr')

system.time(imputed_data_mlr.n_xc_polr <- mice(df_missing_mlr.n_x_c, m = 15, maxit = 20, method = 'polr'))
#time elapsed: t = 295
# Save the data set
save(imputed_data_mlr.n_xc_polr, file = "imputed_data_mlr.n_xc_polr.RData")

# View the imputed dataset
completed_data_mlr.n_xc_polr <- complete(imputed_data_mlr.n_xc_polr)

md.pattern(completed_data_mlr.n_xc_polr)

frequency_xc_mlr.n_imputed_polr <- table(completed_data_mlr.n_xc_polr$c)
frequency_xc_mlr.n_imputed_polr


#convergence check
mlr.n.xc.polr.con.plot <- plot(imputed_data_mlr.n_xc_polr, y = "c", layout = c(1,2))

# Start the analyses
MIpolr.mlr.n.xc <- with(data = imputed_data_mlr.n_xc_polr , glm(y ~ x + c))
MIpolr.mlr.n.xc

#Rubins rule
pooledEst.polr.mlr.n.xc <- pool(MIpolr.mlr.n.xc)

# Create a data frame with the coefficients and standard errors
coeff_se_xc_data_mlr.n_polr <- data.frame(
  Paramters = pooledEst.polr.mlr.n.xc$pooled["term"],
  Estimate = pooledEst.polr.mlr.n.xc$pooled["estimate"],
  Std_Error = sqrt(pooledEst.polr.mlr.n.xc$pooled["t"]) #t = total variance
)

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.polr.xc <- coeff_se_xc_data_mlr.n_polr[coeff_se_xc_data_mlr.n_polr$term %in% c("c2", "c3"), ]

# View the imputed dataset
completed_data_mlr.n_xc_polr <- complete(imputed_data_mlr.n_xc_polr)

md.pattern(completed_data_mlr.n_xc_polr)

#frquency moet nog

frequency_xc_mlr.n_imputed_polr <- table(completed_data_mlr.n_xc_polr$c)

xc_data_imputed_difference_polr <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_full_mlr.n/ss)*100,
  After <- (frequency_xc_mlr.n_imputed_polr/ss)*100,
  Imputed <- (( frequency_xc_mlr.n_imputed_polr - frequency_c_mlr.n)/ss)*100
)


#GRAPH: side-by-side bar chart full data vs polyreg imputed data

# Create tables for frequency

before.fq.mlr.n <- (frequency_c_mlr.n/ss)*100

fq.full.data.mlr.n <- data.frame(
  Category = c("1", "2", "3"),
  Value = before.fq.mlr.n 
)

fq.polr.xc.imputed.mlr.n <- (frequency_xc_mlr.n_imputed_polr/ss)*100

fq.xc.polr.data.mlr.n  <- data.frame(
  Category = c("1", "2", "3"),
  Value = fq.polr.xc.imputed.mlr.n
)



# Add an identifier column to distinguish between the tables
fq.full.data.mlr.n$Table <- "1. Full data set"
fq.xc.polr.data.mlr.n$Table <- "2. Imputed data set"
fq.xc.amputed.data.mlr.n$Table <- "3. Amputed data set"

# Combine the tables into one data frame
combined_data.full.and.xc.polr <- rbind(fq.full.data.mlr.n, fq.xc.polr.data.mlr.n, fq.xc.amputed.data.mlr.n)

# Create the side-by-side bar plot
plot1.mlr.n.xc.polr <- ggplot(combined_data.full.and.xc.polr, aes(x = Category, y = Value.Freq, fill = Table)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(title = "3. Imputed with polr method",
       x = "Category",
       y = "Frequency") +
  scale_fill_manual(values = c("1. Full data set"= "skyblue", "2. Imputed data set" = "lightgreen", "3. Amputed data set" = "#1f77b4")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal()

plot1.mlr.n.xc.polr


################################################################################
#GRAPHS

### FREQUENCY

#grid the frequency plots
ggarrange(plot1.mlr.n.xc.pmm, plot1.mlr.n.xc.polyreg, plot1.mlr.n.xc.polr, 
          nrow = 3, common.legend = TRUE, legend="bottom")




################################################################################
################################################################################
### mlr.n X, Y and C scenario ###
#load mlr.n x, y c scenario file

load("df_missing_mlr.n_xyc.RData")

frequency_xyc_mlr.n <- table(df_missing_mlr.n_xyc$c)


cca_xyc.mlr.S1 <- glm(y ~ x + c, data = df_missing_mlr.n_xyc)

summary(cca_xyc.mlr.S1)

#Frequency c in amputed full data
frequency_xyc_missing_mlr.n <- (table(df_missing_mlr.n_xyc$c)/ss)*100

##### PMM #####
md.pattern(df_missing_mlr.n_xyc)

# Convert 'c' to a factor variable 
df_missing_mlr.n_xyc$c <- as.factor(df_missing_mlr.n_xyc$c)

# Perform multiple imputation using the "pmm" method
imputed_data_mlr.n_xyc_pmm <- mice(df_missing_mlr.n_xyc, m = 15, maxit = 20, method = 'pmm')

system.time(imputed_data_mlr.n_xyc_pmm <- mice(df_missing_mlr.n_xyc, m = 15, maxit = 20, method = 'pmm'))
#time elapsed: t = 13 sec
# Save the data set
save(imputed_data_mlr.n_xyc_pmm, file = "imputed_data_mlr.n_xyc_pmm.RData")

#convergence check
mlr.n.xyc.pmm.con.plot <- plot(imputed_data_mlr.n_xyc_pmm , y = "c", layout = c(1,2))

# Start the analyses
MIpmm.mlr.n.xyc <- with(data = imputed_data_mlr.n_xyc_pmm , glm(y ~ x + c))
MIpmm.mlr.n.xyc

#Rubins rule
pooledEst.pmm.mlr.n.xyc <- pool(MIpmm.mlr.n.xyc)

# Create a data frame with the coefficients and standard errors
coeff_se_xyc_data_mlr.n_pmm <- data.frame(
  Paramters = pooledEst.pmm.mlr.n.xyc$pooled["term"],
  Estimate =  pooledEst.pmm.mlr.n.xyc$pooled["estimate"],
  Std_Error = sqrt(pooledEst.pmm.mlr.n.xyc$pooled["t"]) #t = total variance
)

c2c3.coef.se.pmm.xyc <- coeff_se_xyc_data_mlr.n_pmm[coeff_se_xyc_data_mlr.n_pmm$term %in% c("c2", "c3"), ]

# View the imputed dataset
completed_data_mlr.n_xyc_pmm <- complete(imputed_data_mlr.n_xyc_pmm)

md.pattern(completed_data_mlr.n_xyc_pmm)

#frquency moet nog

frequency_xyc_mlr.n_imputed_pmm <- table(completed_data_mlr.n_xyc_pmm$c)

xyc_data_imputed_difference_pmm <- data.frame(
  Category = c("1", "2", "3"),
  Missing <- ((frequency_c_mlr.n - frequency_c_full_mlr.n)/ss)*100,
  #After <- (frequency_c_mlr.n_imputed_pmm/ss)*100,
  Imputed <- (( frequency_xyc_mlr.n_imputed_pmm - frequency_c_mlr.n)/ss)*100
)

#GRAPH: side-by-side bar chart full data vs polyreg imputed data

# Create tables for frequency

before.fq.mlr.n <- (frequency_c_mlr.n/ss)*100

fq.full.data.mlr.n <- data.frame(
  Category = c("1", "2", "3"),
  Value = before.fq.mlr.n 
)

fq.pmm.xyc.imputed.mlr.n <- (frequency_xyc_mlr.n_imputed_pmm/ss)*100

fq.xyc.pmm.data.mlr.n  <- data.frame(
  Category = c("1", "2", "3"),
  Value = fq.pmm.xyc.imputed.mlr.n
)

fq.xyc.amputed.full.data.mlr.n <- data.frame(
  Category = c("1", "2", "3"),
  Value = frequency_xyc_missing_mlr.n
)

# Add an identifier column to distinguish between the tables
fq.full.data.mlr.n$Table <- "1. Full data set"
fq.xyc.pmm.data.mlr.n$Table <- "2. Imputed data set"
fq.xyc.amputed.full.data.mlr.n$Table <- "3. Amputed data set"

# Combine the tables into one data frame
combined_data.full.and.xyc.pmm <- rbind(fq.full.data.mlr.n, fq.xyc.pmm.data.mlr.n, fq.xyc.amputed.full.data.mlr.n)

# Create the side-by-side bar plot
plot1.mlr.n.xyc.pmm <- ggplot(combined_data.full.and.xyc.pmm, aes(x = Category, y = Value.Freq, fill = Table)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(title = "1. Impution with pmm method",
       x = "Category",
       y = "Frequency") +
  scale_fill_manual(values = c("1. Full data set"= "skyblue", "2. Imputed data set" = "lightgreen", "3. Amputed data set" = "#1f77b4")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal()

plot1.mlr.n.xyc.pmm 



################################################################################
### POLYREG ###

# Convert 'c' to a factor variable 
md.pattern(df_missing_mlr.n_xyc)

df_missing_mlr.n_xyc$c <- as.factor(df_missing_mlr.n_xyc$c)

# Perform multiple imputation using the "polyreg" method
imputed_data_mlr.n_xyc_polyreg <- mice(df_missing_mlr.n_xyc, m = 15, maxit = 20, method = 'polyreg')

system.time(imputed_data_mlr.n_xyc_polyreg <- mice(df_missing_mlr.n_xyc, m = 15, maxit = 20, method = 'polyreg'))
#time elapsed: t = 130secs
# Save the data set
save(imputed_data_mlr.n_xyc_polyreg, file = "imputed_data_mlr.n_xyc_polyreg.RData")

# View the imputed dataset
completed_data_mlr.n_xyc_polyreg <- complete(imputed_data_mlr.n_xyc_polyreg)

md.pattern(completed_data_mlr.n_xyc_polyreg)

frequency_xyc_mlr.n_imputed_polyreg <- table(completed_data_mlr.n_xyc_polyreg$c)
frequency_xyc_mlr.n_imputed_polyreg


#convergence check
mlr.n.xyc.polyreg.con.plot <-plot(imputed_data_mlr.n_xyc_polyreg, y = "c", layout = c(1,2))

# Start the analyses
MIpolyreg.mlr.n.xyc <- with(data = imputed_data_mlr.n_xyc_polyreg , glm(y ~ x + c))
MIpolyreg.mlr.n.xyc

#Rubins rule
pooledEst.polyreg.mlr.n.xyc <- pool(MIpolyreg.mlr.n.xyc)


# Create a data frame with the coefficients and standard errors
coeff_se_xyc_data_mlr.n_polyreg <- data.frame(
  Paramters = pooledEst.polyreg.mlr.n.xyc$pooled["term"],
  Estimate = pooledEst.polyreg.mlr.n.xyc$pooled["estimate"],
  Std_Error = sqrt(pooledEst.polyreg.mlr.n.xyc$pooled["t"]) #t = total variance
)

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.polyreg.xyc <- coeff_se_xyc_data_mlr.n_polyreg[coeff_se_xyc_data_mlr.n_polyreg$term %in% c("c2", "c3"), ]


# View the imputed dataset
completed_data_mlr.n_xyc_polyreg <- complete(imputed_data_mlr.n_xyc_polyreg)

md.pattern(completed_data_mlr.n_xyc_polyreg)

#frquency moet nog

frequency_xyc_mlr.n_imputed_polyreg <- table(completed_data_mlr.n_xyc_polyreg$c)

xyc_data_imputed_difference_polyreg <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_full_mlr.n/ss)*100,
  After <- (frequency_xyc_mlr.n_imputed_polyreg/ss)*100,
  Imputed <- (( frequency_xyc_mlr.n_imputed_polyreg - frequency_c_mlr.n)/ss)*100
)


#GRAPH: side-by-side bar chart full data vs polyreg imputed data

# Create tables for frequency

before.fq.mlr.n <- (frequency_c_mlr.n/ss)*100

fq.full.data.mlr.n <- data.frame(
  Category = c("1", "2", "3"),
  Value = before.fq.mlr.n 
)

fq.polyreg.xyc.imputed.mlr.n <- (frequency_xyc_mlr.n_imputed_polyreg/ss)*100

fq.xyc.polyreg.data.mlr.n  <- data.frame(
  Category = c("1", "2", "3"),
  Value = fq.polyreg.xyc.imputed.mlr.n
)

# Add an identifier column to distinguish between the tables
fq.full.data.mlr.n$Table <- "1. Full data set"
fq.xyc.polyreg.data.mlr.n$Table <- "2. Imputed data set"
fq.xyc.amputed.full.data.mlr.n$Table <- "3. Amputed data set"

# Combine the tables into one data frame
combined_data.full.and.xc.polyreg <- rbind(fq.full.data.mlr.n, fq.xyc.polyreg.data.mlr.n, fq.xyc.amputed.full.data.mlr.n)

# Create the side-by-side bar plot
plot1.mlr.n.xyc.polyreg <- ggplot(combined_data.full.and.xc.polyreg, aes(x = Category, y = Value.Freq, fill = Table)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(title = "2. Imputed with polyreg method",
       x = "Category",
       y = "Frequency") +
  scale_fill_manual(values = c("1. Full data set"= "skyblue", "2. Imputed data set" = "lightgreen", "3. Amputed data set" = "#1f77b4")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal()

plot1.mlr.n.xyc.polyreg

################################################################################
### POLR ###

md.pattern(df_missing_mlr.n)

# Convert 'c' to a factor variable 
md.pattern(df_missing_mlr.n_xyc)

df_missing_mlr.n_xyc$c <- as.factor(df_missing_mlr.n_xyc$c)

# Perform multiple imputation using the "polr" method
imputed_data_mlr.n_xyc_polr <- mice(df_missing_mlr.n_xyc, m = 15, maxit = 20, method = 'polr')

system.time(imputed_data_mlr.n_xyc_polr <- mice(df_missing_mlr.n_xyc, m = 15, maxit = 20, method = 'polr'))
#time elapsed: t = 288 secs
# Save the data set
save(imputed_data_mlr.n_xyc_polr, file = "imputed_data_mlr.n_xyc_polr.RData")

# View the imputed dataset
completed_data_mlr.n_xyc_polr <- complete(imputed_data_mlr.n_xyc_polr)

md.pattern(completed_data_mlr.n_xyc_polr)

frequency_xyc_mlr.n_imputed_polr <- table(completed_data_mlr.n_xyc_polr$c)
frequency_xyc_mlr.n_imputed_polr


#convergence check
mlr.n.xyc.polr.con.plot <- plot(imputed_data_mlr.n_xyc_polr, y = "c", layout = c(1,2))

# Start the analyses
MIpolr.mlr.n.xyc <- with(data = imputed_data_mlr.n_xyc_polr , glm(y ~ x + c))
MIpolr.mlr.n.xyc

#Rubins rule
pooledEst.polr.mlr.n.xyc <- pool(MIpolr.mlr.n.xyc)


# Create a data frame with the coefficients and standard errors
coeff_se_xyc_data_mlr.n_polr <- data.frame(
  Paramters = pooledEst.polr.mlr.n.xyc$pooled["term"],
  Estimate = pooledEst.polr.mlr.n.xyc$pooled["estimate"],
  Std_Error = sqrt(pooledEst.polr.mlr.n.xyc$pooled["t"]) #t = total variance
)

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.polr.xyc <- coeff_se_xyc_data_mlr.n_polr[coeff_se_xyc_data_mlr.n_polr$term %in% c("c2", "c3"), ]

# View the imputed dataset
completed_data_mlr.n_xyc_polr <- complete(imputed_data_mlr.n_xyc_polr)

md.pattern(completed_data_mlr.n_xyc_polr)

#frquency moet nog

frequency_xyc_mlr.n_imputed_polr <- table(completed_data_mlr.n_xyc_polr$c)

xyc_data_imputed_difference_polr <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_full_mlr.n/ss)*100,
  After <- (frequency_xyc_mlr.n_imputed_polr/ss)*100,
  Imputed <- (( frequency_xyc_mlr.n_imputed_polr - frequency_c_mlr.n)/ss)*100
)

#GRAPH: side-by-side bar chart full data vs polyreg imputed data

# Create tables for frequency

before.fq.mlr.n <- (frequency_c_mlr.n/ss)*100

fq.full.data.mlr.n <- data.frame(
  Category = c("1", "2", "3"),
  Value = before.fq.mlr.n 
)

fq.polr.xyc.imputed.mlr.n <- (frequency_xyc_mlr.n_imputed_polr/ss)*100

fq.xyc.polr.data.mlr.n  <- data.frame(
  Category = c("1", "2", "3"),
  Value = fq.polr.xyc.imputed.mlr.n
)

# Add an identifier column to distinguish between the tables
fq.full.data.mlr.n$Table <- "1. Full data set"
fq.xyc.polr.data.mlr.n$Table <- "2. Imputed data set"
fq.xyc.amputed.full.data.mlr.n$Table <- "3. Amputed data set"

# Combine the tables into one data frame
combined_data.full.and.xc.polr <- rbind(fq.full.data.mlr.n, fq.xyc.polr.data.mlr.n, fq.xyc.amputed.full.data.mlr.n)

# Create the side-by-side bar plot
plot1.mlr.n.xyc.polr <- ggplot(combined_data.full.and.xc.polr, aes(x = Category, y = Value.Freq, fill = Table)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(title = "3. Imputed with polr method",
       x = "Category",
       y = "Frequency") +
  scale_fill_manual(values = c("1. Full data set"= "skyblue", "2. Imputed data set" = "lightgreen", "3. Amputed data set" = "#1f77b4")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal()

plot1.mlr.n.xyc.polr

################################################################################
#GRAPHS

### FREQUENCY
#grid the frequency plots
ggarrange(plot1.mlr.n.xyc.pmm, plot1.mlr.n.xyc.polyreg, plot1.mlr.n.xyc.polr, 
          nrow = 3, common.legend = TRUE, legend="bottom")


###### PLOT all coef in grid
create_coefficient_point_plot.xyc <- function(data, title) {
  ggplot(data, aes(x = term, y = estimate)) +
    geom_point() +
    labs(
      title = title,
      x = "Paramter",
      y = "Estimate"
    )
}


# Plots:
coef.plot1.xyc <- create_coefficient_point_plot(coeff_se_xyc_data_mlr.n_pmm, "PMM")
coef.plot2.xyc <- create_coefficient_point_plot(coeff_se_xyc_data_mlr.n_polyreg, "Polyreg")
coef.plot3.xyc <- create_coefficient_point_plot(coeff_se_xyc_data_mlr.n_polr, "POLR")

# Arrange all the graphs in a 2x2 grid
grid.arrange(coef.plot1.xyc, coef.plot2.xyc, coef.plot3.xyc, ncol = 3)




### FREQUENCY
#grid the frequency plots
mlr.xyc.fqr.plot <- ggarrange(plot1.mlr.n.xyc.pmm, plot1.mlr.n.xyc.polyreg, plot1.mlr.n.xyc.polr, 
                             nrow = 3, common.legend = TRUE, legend="bottom")

annotate_figure(mlr.xyc.fqr.plot, top = text_grob("Multinominal Logistic Regression Model Scenario 2", 
                                                 color = "blue", face = "bold", size = 14))

mlr.xc.fqr.plot <- ggarrange(plot1.mlr.n.xc.pmm, plot1.mlr.n.xc.polyreg, plot1.mlr.n.xc.polr, 
                            nrow = 3, common.legend = TRUE, legend="bottom")

annotate_figure(mlr.xc.fqr.plot, top = text_grob("Multinominal Logistic Regression Model Scenario 1", 
                                                color = "blue", face = "bold", size = 14))


mlr.n.xc.pmm.con.plot
mlr.n.xc.polyreg.con.plot
mlr.n.xc.polr.con.plot
mlr.n.xyc.pmm.con.plot
mlr.n.xyc.polyreg.con.plot
mlr.n.xyc.polr.con.plot
                              


