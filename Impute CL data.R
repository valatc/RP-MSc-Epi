############## CL DATA  IMPUTATION ##############
#Packages
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
load("df_missing_cl.RData")

load("data_cl.RData")

frequency_c_cl <- table(data_cl$c)

#Frequency c in amputed full data
frequency_c_full_cl <- table(df_missing_cl$c)


#sample size, for turning the frequencies into %
ss <- 100000

#subsetting data to show only y, x, and c in amputed data set
sub.df_missing_cl <- subset(df_missing_cl, select = c("y", "c", "x"))

#CCA of amputed data

#S1

load()

################################################################################
################################################################################

#### FULL DATA ####
##### PMM #####
md.pattern(sub.df_missing_cl)

df_missing_cl$c <- as.factor(df_missing_cl$c)

# Perform multiple imputation using the "pmm" method
imputed_data_cl_full_pmm <- mice(sub.df_missing_cl, m = 15, maxit = 20, method = 'pmm')


# Save the imputed data set
save(imputed_data_cl_full_pmm, file = "imputed_data_cl_full_pmm.RData")

laod("imputed_data_cl_full_pmm.RData")
#convergence check
cl.full.pmm.con.plot <- plot(imputed_data_cl_full_pmm, y = "c", layout = c(1,2))

# Start the analyses
MIpmm.cl.full <- with(data = imputed_data_cl_full_pmm , glm(y ~ x + c))
MIpmm.cl.full

#Rubins rule
pooledEst.pmm.cl.full <- pool(MIpmm.cl.full)

# Create a data frame with the coefficients and standard errors
coeff_se_full_data_cl_pmm <- data.frame(
  Paramters = pooledEst.pmm.cl.full$pooled["term"],
  Estimate = pooledEst.pmm.cl.full$pooled["estimate"],
  Std_Error = sqrt(pooledEst.pmm.cl.full$pooled["t"]) #t = total variance
)

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.pmm.full <- coeff_se_full_data_cl_pmm[coeff_se_full_data_cl_pmm$term %in% c("c"), ]

# View the imputed dataset
completed_data_cl_full_pmm <- complete(imputed_data_cl_full_pmm)

md.pattern(completed_data_cl_full_pmm)

#frquency moet nog

frequency_c_cl_imputed_pmm <- table(completed_data_cl_full_pmm$c)

full_data_imputed_difference_pmm <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_cl/ss)*100,
  Imputed <- (frequency_c_cl_imputed_pmm/ss)*100
)

frequency_c_cl
frequency_c_full_cl
################################################################################
### POLYREG ###
md.pattern(sub.df_missing_cl)

# Convert 'c' to a factor variable 
sub.df_missing_cl$c <- as.factor(sub.df_missing_cl$c)


# Perform multiple imputation using the "polyreg" method
imputed_data_cl_full_polyreg <- mice(sub.df_missing_cl, m = 15, maxit = 20, method = 'polyreg')

# Save the data set
save(imputed_data_cl_full_polyreg , file = "imputed_data_cl_full_polyreg .RData")

# View the imputed dataset
completed_data_cl_full_polyreg <- complete(imputed_data_cl_full_polyreg)

md.pattern(completed_data_cl_full_polyreg)

frequency_c_cl_imputed_polyreg <- table(completed_data_cl_full_polyreg$c)
frequency_c_cl_imputed_polyreg


#convergence check
cl.full.polyreg.con.plot <-plot(imputed_data_cl_full_polyreg, y = "c", layout = c(1,2))

# Start the analyses
MIpolyreg.cl.full <- with(data = imputed_data_cl_full_polyreg , glm(y ~ x + c))
MIpolyreg.cl.full

#Rubins rule
pooledEst.polyreg.cl.full <- pool(MIpolyreg.cl.full)

# Create a data frame with the coefficients and standard errors
coeff_se_full_data_cl_polyreg <- data.frame(
  Coefficients = c("", "", "b2", "b3"),
  #Paramters = pooledEst.polyreg.cl.full$pooled["term"],
  Estimate = pooledEst.polyreg.cl.full$pooled["estimate"],
  Std_Error = sqrt(pooledEst.polyreg.cl.full$pooled["t"]) #t = total variance
)

#change column names
colnames(coeff_se_full_data_cl_polyreg) <- c("Parameter", "Estimate", "Std_Error")

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.polyreg.full <- data.frame(coeff_se_full_data_cl_polyreg [coeff_se_full_data_cl_polyreg $Coefficients  %in% c("b2", "b3"), ])


#frequency
f.before <- (frequency_c_cl/ss)*100
f.imputed_polyreg <- (frequency_c_cl_imputed_pmm/ss)*100

coef.f.cl.polyreg <- data.frame(c2c3.coef.se.polyreg.full)

comb.polyreg <- rbind(f.before, f.imputed_polyreg)



################################################################################
### POLR ###

md.pattern(df_missing_cl)


# Convert 'c' to a factor variable 
sub.df_missing_cl$c <- as.factor(sub.df_missing_cl$c)

# Perform multiple imputation using the "polr" method
imputed_data_cl_full_polr <- mice(sub.df_missing_cl, m = 15, maxit = 20, method = 'polr')

# Save the data set
save(imputed_data_cl_full_polr , file = "imputed_data_cl_full_polr .RData")

# View the imputed dataset
completed_data_cl_full_polr <- complete(imputed_data_cl_full_polr)

md.pattern(completed_data_cl_full_polr)

frequency_c_cl_imputed_polr <- table(completed_data_cl_full_polr$c)
frequency_c_cl_imputed_polr


#convergence check
cl.full.polr.con.plot <- plot(imputed_data_cl_full_polr, y = "c", layout = c(1,2))

# Start the analyses
MIpolr.cl.full <- with(data = imputed_data_cl_full_polr , glm(y ~ x + c))
MIpolr.cl.full

#Rubins rule
pooledEst.polr.cl.full <- pool(MIpolr.cl.full)


# Create a data frame with the coefficients and standard errors
coeff_se_full_data_cl_polr <- data.frame(
  Paramters = pooledEst.polr.cl.full$pooled["term"],
  Estimate = pooledEst.polr.cl.full$pooled["estimate"],
  Std_Error = sqrt(pooledEst.polr.cl.full$pooled["t"]) #t = total variance
)

#change column names
colnames(coeff_se_full_data_cl_polr) <- c("Parameter", "Estimate", "Std_Error")

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.polr.full <- coeff_se_full_data_cl_polr [coeff_se_full_data_cl_polr $Parameter %in% c("c2", "c3"), ]

# View the imputed dataset
completed_data_cl_full_polr <- complete(imputed_data_cl_full_polr)

md.pattern(completed_data_cl_full_polr)

#frquency moet nog

frequency_c_cl_imputed_polr <- table(completed_data_cl_full_polr$c)

full_data_imputed_difference_polr <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_cl/ss)*100,
  Imputed <- (frequency_c_cl_imputed_polr/ss)*100
)


################################################################################
#GRAPHS

### FREQUENCY
# Function to create the bar chart
create_bar_chart_cl <- function(data, title) {
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
frq.plot1.full <- create_bar_chart_cl(full_data_imputed_difference_pmm, "Imputed with PMM")
frq.plot2.full <- create_bar_chart_cl(full_data_imputed_difference_polyreg, "Imputed with POLYREG")
frq.plot3.full <- create_bar_chart_cl(full_data_imputed_difference_polr, "Imputed with POLR")

#grid the frequency plots
ggarrange(frq.plot1.full, frq.plot2.full, frq.plot3.full, ncol=3, common.legend = TRUE, legend="bottom")

ggarrange(cl.full.polr.con.plot, cl.full.polyreg.con.plot, cl.full.pmm.con.plot, 
          common.legend = TRUE, legend="bottom")


# Table coef combined
comb.cl.polyreg.polr.full <- rbind(c2c3.coef.se.polyreg.full, c2c3.coef.se.polr.full)
comb.cl.polyreg.polr.full
a <-as.data.frame(c2c3.coef.se.polyreg.full)
b <- as.data.frame(c2c3.coef.se.polr.full)

flextable(data = list(a,b) , col_keys = "beta coef", "estimate", "std.error")

################################################################################
################################################################################
### CL X and C scenario ###
#load cl x and c scenario file

load("df_missing_cl_x_c.RData")


#CCA

cca_xc.cl.S1 <- glm(y ~ x + c, data = df_missing_cl_x_c)

summary(cca_xc.cl.S1)

##### PMM #####
# Convert 'c' to a factor variable 
md.pattern(df_missing_cl_x_c)

df_missing_cl_x_c$c <- as.factor(df_missing_cl_x_c$c)

# Perform multiple imputation using the "pmm" method
imputed_data_cl_xc_pmm <- mice(df_missing_cl_x_c, m = 15, maxit = 20, method = 'pmm')

# Save the data set
save(imputed_data_cl_xc_pmm, file = "imputed_data_cl_xc_pmm.RData")
load("imputed_data_cl_xc_pmm.RData")

frequency_xc_cl <- table(df_missing_cl_x_c$c)

#Frequency c in amputed full data
frequency_xc_missing_cl <- (table(df_missing_cl_x_c$c)/ss)*100


#convergence check
cl.xc.pmm.con.plot <- plot(imputed_data_cl_xc_pmm , y = "c", layout = c(1,2))

# Start the analyses
MIpmm.cl.xc <- with(data = imputed_data_cl_xc_pmm , glm(y ~ x + c))
MIpmm.cl.xc

#Rubins rule
pooledEst.pmm.cl.xc <- pool(MIpmm.cl.xc)


# Create a data frame with the coefficients and standard errors
coeff_se_xc_data_cl_pmm <- data.frame(
  Paramters = pooledEst.pmm.cl.xc$pooled["term"],
  Estimate = pooledEst.pmm.cl.xc$pooled["estimate"],
  Std_Error = sqrt(pooledEst.pmm.cl.xc$pooled["t"]) #t = total variance
)

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.pmm.xc <- coeff_se_xc_data_cl_pmm[coeff_se_xc_data_cl_pmm$term %in% c("c2", "c3"), ]

# View the imputed dataset
completed_data_cl_xc_pmm <- complete(imputed_data_cl_xc_pmm)

md.pattern(completed_data_cl_xc_pmm)

#frquency moet nog

frequency_xc_cl_imputed_pmm <- table(completed_data_cl_xc_pmm$c)

xc_data_imputed_difference_pmm <- data.frame(
  Category = c("1", "2", "3"),
  Missing <- ((frequency_c_cl - frequency_c_full_cl)/ss)*100,
  #After <- (frequency_c_cl_imputed_pmm/ss)*100,
  Imputed <- (( frequency_xc_cl_imputed_pmm - frequency_c_cl)/ss)*100
)

#GRAPH: side-by-side bar chart full data vs pmm imputed data

# Create tables for frequency

before.fq.cl <- (frequency_c_cl/ss)*100

fq.full.data.cl <- data.frame(
  Category = c("1", "2", "3"),
  Value = before.fq.cl 
)

fq.pmm.xc.imputed.cl <- (frequency_xc_cl_imputed_pmm/ss)*100

fq.xc.pmm.data.cl  <- data.frame(
  Category = c("1", "2", "3"),
  Value = fq.pmm.xc.imputed.cl
)


fq.xc.amputed.data.cl <- data.frame(
  Category = c("1", "2", "3"),
  Value = frequency_xc_missing_cl
)

# Add an identifier column to distinguish between the tables
fq.full.data.cl$Table <- "1. Full data set"
fq.xc.pmm.data.cl$Table <- "2. Imputed data set"
fq.xc.amputed.data.cl$Table <- "3. Amputed data set"

# Combine the tables into one data frame
combined_cl.data.full.and.xc.pmm <- rbind(fq.full.data.cl, fq.xc.pmm.data.cl, fq.xc.amputed.data.cl)

# Create the side-by-side bar plot
plot1.cl.xc.pmm <- ggplot(combined_cl.data.full.and.xc.pmm, aes(x = Category, y = Value.Freq, fill = Table)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(title = "1. Imputed with pmm method",
       x = "Category",
       y = "Frequency") +
  scale_fill_manual(values = c("1. Full data set"= "skyblue", "2. Imputed data set" = "lightgreen", "3. Amputed data set" = "#1f77b4")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal()

plot1.cl.xc.pmm





################################################################################
### POLYREG ###

# Convert 'c' to a factor variable 
md.pattern(df_missing_cl_x_c)

df_missing_cl_x_c$c <- as.factor(df_missing_cl_x_c$c)

# Perform multiple imputation using the "polyreg" method
imputed_data_cl_xc_polyreg <- mice(df_missing_cl_x_c, m = 15, maxit = 20, method = 'polyreg')

# Save the data set
save(imputed_data_cl_xc_polyreg, file = "imputed_data_cl_xc_polyreg.RData")

load("imputed_data_cl_xc_polyreg.RData")

# View the imputed dataset
completed_data_cl_xc_polyreg <- complete(imputed_data_cl_xc_polyreg)

md.pattern(completed_data_cl_xc_polyreg)

frequency_xc_cl_imputed_polyreg <- table(completed_data_cl_xc_polyreg$c)
frequency_xc_cl_imputed_polyreg


#convergence check
cl.xc.polyreg.con.plot <-plot(imputed_data_cl_xc_polyreg, y = "c", layout = c(1,2))

# Start the analyses
MIpolyreg.cl.xc <- with(data = imputed_data_cl_xc_polyreg , glm(y ~ x + c))
MIpolyreg.cl.xc

#Rubins rule
pooledEst.polyreg.cl.xc <- pool(MIpolyreg.cl.xc)

# Create a data frame with the coefficients and standard errors
coeff_se_xc_data_cl_polyreg <- data.frame(
  Paramters = pooledEst.polyreg.cl.xc$pooled["term"],
  Estimate = pooledEst.polyreg.cl.xc$pooled["estimate"],
  Std_Error = sqrt(pooledEst.polyreg.cl.xc$pooled["t"]) #t = total variance
)

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.polyreg.xc <- coeff_se_xc_data_cl_polyreg[coeff_se_xc_data_cl_polyreg$term %in% c("c2", "c3"), ]

# View the imputed dataset
completed_data_cl_xc_polyreg <- complete(imputed_data_cl_xc_polyreg)

md.pattern(completed_data_cl_xc_polyreg)

#frquency moet nog

frequency_xc_cl_imputed_polyreg <- table(completed_data_cl_xc_polyreg$c)

xc_data_imputed_difference_polyreg <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_full_cl/ss)*100,
  After <- (frequency_xc_cl_imputed_polyreg/ss)*100,
  Imputed <- (( frequency_xc_cl_imputed_polyreg - frequency_c_cl)/ss)*100
)


#GRAPH: side-by-side bar chart full data vs polyreg imputed data

# Create tables for frequency

before.fq.cl <- (frequency_c_cl/ss)*100

fq.full.data.cl <- data.frame(
  Category = c("1", "2", "3"),
  Value = before.fq.cl 
)

fq.polyreg.xc.imputed.cl <- (frequency_xc_cl_imputed_polyreg/ss)*100

fq.xc.polyreg.data.cl  <- data.frame(
  Category = c("1", "2", "3"),
  Value = fq.polyreg.xc.imputed.cl
)


fq.xyc.amputed.data.cl <- data.frame(
  Category = c("1", "2", "3"),
  Value = frequency_xyc_missing_cl
)


# Add an identifier column to distinguish between the tables
fq.full.data.cl$Table <- "1. Full data set"
fq.xc.polyreg.data.cl$Table <- "2. Imputed data set"
fq.xc.amputed.data.cl$Table <- "3. Amputed data set"

# Combine the tables into one data frame
combined_cl.data.full.and.xc.polyreg <- rbind(fq.full.data.cl, fq.xc.polyreg.data.cl, fq.xc.amputed.data.cl)

# Create the side-by-side bar plot
plot1.cl.xc.polyreg <- ggplot(combined_cl.data.full.and.xc.polyreg, aes(x = Category, y = Value.Freq, fill = Table)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(title = "2. Imputed with polyreg method",
       x = "Category",
       y = "Frequency") +
  scale_fill_manual(values = c("1. Full data set"= "skyblue", "2. Imputed data set" = "lightgreen", "3. Amputed data set" = "#1f77b4")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal()

plot1.cl.xc.polyreg

################################################################################
### POLR ###

md.pattern(df_missing_cl)

# Convert 'c' to a factor variable 
md.pattern(df_missing_cl_x_c)

df_missing_cl_x_c$c <- as.factor(df_missing_cl_x_c$c)

# Perform multiple imputation using the "polr" method
imputed_data_cl_xc_polr <- mice(df_missing_cl_x_c, m = 15, maxit = 20, method = 'polr')

# Save the data set
save(imputed_data_cl_xc_polr, file = "imputed_data_cl_xc_polr.RData")

load("imputed_data_cl_xc_polr.RData")
# View the imputed dataset
completed_data_cl_xc_polr <- complete(imputed_data_cl_xc_polr)

md.pattern(completed_data_cl_xc_polr)

frequency_xc_cl_imputed_polr <- table(completed_data_cl_xc_polr$c)
frequency_xc_cl_imputed_polr


#convergence check
cl.xc.polr.con.plot <- plot(imputed_data_cl_xc_polr, y = "c", layout = c(1,2))

# Start the analyses
MIpolr.cl.xc <- with(data = imputed_data_cl_xc_polr , glm(y ~ x + c))
MIpolr.cl.xc

#Rubins rule
pooledEst.polr.cl.xc <- pool(MIpolr.cl.xc)

# Create a data frame with the coefficients and standard errors
coeff_se_xc_data_cl_polr <- data.frame(
  Paramters = pooledEst.polr.cl.xc$pooled["term"],
  Estimate = pooledEst.polr.cl.xc$pooled["estimate"],
  Std_Error = sqrt(pooledEst.polr.cl.xc$pooled["t"]) #t = total variance
)

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.polr.xc <- coeff_se_xc_data_cl_polr[coeff_se_xc_data_cl_polr$term %in% c("c2", "c3"), ]

# View the imputed dataset
completed_data_cl_xc_polr <- complete(imputed_data_cl_xc_polr)

md.pattern(completed_data_cl_xc_polr)

#frquency moet nog

frequency_xc_cl_imputed_polr <- table(completed_data_cl_xc_polr$c)

xc_data_imputed_difference_polr <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_full_cl/ss)*100,
  After <- (frequency_xc_cl_imputed_polr/ss)*100,
  Imputed <- (( frequency_xc_cl_imputed_polr - frequency_c_cl)/ss)*100
)


#GRAPH: side-by-side bar chart full data vs polyreg imputed data

# Create tables for frequency

before.fq.cl <- (frequency_c_cl/ss)*100

fq.full.data.cl <- data.frame(
  Category = c("1", "2", "3"),
  Value = before.fq.cl 
)

fq.polr.xc.imputed.cl <- (frequency_xc_cl_imputed_polr/ss)*100

fq.xc.polr.data.cl  <- data.frame(
  Category = c("1", "2", "3"),
  Value = fq.polr.xc.imputed.cl
)

# Add an identifier column to distinguish between the tables
fq.full.data.cl$Table <- "1. Full data set"
fq.xc.polr.data.cl$Table <- "2. Imputed data set"
fq.xc.amputed.data.cl$Table <- "3. Amputed data set"

# Combine the tables into one data frame
combined_cl.data.full.and.xc.polr <- rbind(fq.full.data.cl, fq.xc.polr.data.cl, fq.xc.amputed.data.cl)

# Create the side-by-side bar plot
plot1.cl.xc.polr <- ggplot(combined_cl.data.full.and.xc.polr, aes(x = Category, y = Value.Freq, fill = Table)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(title = "3. Imputed with polr method",
       x = "Category",
       y = "Frequency") +
  scale_fill_manual(values = c("1. Full data set"= "skyblue", "2. Imputed data set" = "lightgreen", "3. Amputed data set" = "#1f77b4")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal()

plot1.cl.xc.polr

################################################################################
#GRAPHS

### FREQUENCY

#grid the frequency plots
grids1 <- ggarrange(plot1.cl.xc.pmm, plot1.cl.xc.polyreg, plot1.cl.xc.polr, 
          nrow = 3, common.legend = TRUE, legend="bottom")


###### PLOT all coef in grid
create_coefficient_point_plot.xc <- function(data, title) {
  ggplot(data, aes(x = term, y = estimate)) +
    geom_point() +
    labs(
      title = title,
      x = "Paramter",
      y = "Estimate"
    )
}


# Plots:
coef.plot1.xc <- create_coefficient_point_plot(coeff_se_xc_data_cl_pmm, "PMM")
coef.plot2.xc <- create_coefficient_point_plot(coeff_se_xc_data_cl_polyreg, "Polyreg")
coef.plot3.xc <- create_coefficient_point_plot(coeff_se_xc_data_cl_polr, "POLR")

# Arrange all the graphs in a 2x2 grid
grid.arrange(coef.plot1.xc, coef.plot2.xc, coef.plot3.xc, ncol = 3)




################################################################################
################################################################################
### cl X, Y and C scenario ###
#load cl x, y c scenario file

load("df_missing_cl_xyc.RData")

frequency_xyc_missing_cl <- (table(df_missing_cl_xyc$c)/ss)*100

#CCA

cca_xyc.cl.S2 <- glm(y ~ x + c, data = df_missing_cl_xyc)

summary(cca_xyc.cl.S2)
##### PMM #####
# Convert 'c' to a factor variable 
md.pattern(df_missing_cl_xyc)

df_missing_cl_xyc$c <- as.factor(df_missing_cl_xyc$c)

# Perform multiple imputation using the "pmm" method
imputed_data_cl_xyc_pmm <- mice(df_missing_cl_xyc, m = 15, maxit = 20, method = 'pmm')

# Save the data set
save(imputed_data_cl_xyc_pmm, file = "imputed_data_cl_xyc_pmm.RData")

laod("imputed_data_cl_xyc_pmm.RData")

#convergence check
cl.xyc.pmm.con.plot <- plot(imputed_data_cl_xyc_pmm , y = "c", layout = c(1,2))

# Start the analyses
MIpmm.cl.xyc <- with(data = imputed_data_cl_xyc_pmm , glm(y ~ x + c))
MIpmm.cl.xyc

#Rubins rule
pooledEst.pmm.cl.xyc <- pool(MIpmm.cl.xyc)

# Create a data frame with the coefficients and standard errors
coeff_se_xyc_data_cl_pmm <- data.frame(
  Paramters = pooledEst.pmm.cl.xyc$pooled["term"],
  Estimate =  pooledEst.pmm.cl.xyc$pooled["estimate"],
  Std_Error = sqrt(pooledEst.pmm.cl.xyc$pooled["t"]) #t = total variance
)

c2c3.coef.se.pmm.xyc <- coeff_se_xyc_data_cl_pmm[coeff_se_xyc_data_cl_pmm$term %in% c("c2", "c3"), ]

# View the imputed dataset
completed_data_cl_xyc_pmm <- complete(imputed_data_cl_xyc_pmm)

md.pattern(completed_data_cl_xyc_pmm)

#frquency moet nog

frequency_xyc_cl_imputed_pmm <- table(completed_data_cl_xyc_pmm$c)

xyc_data_imputed_difference_pmm <- data.frame(
  Category = c("1", "2", "3"),
  Missing <- ((frequency_c_cl - frequency_c_full_cl)/ss)*100,
  #After <- (frequency_c_cl_imputed_pmm/ss)*100,
  Imputed <- (( frequency_xyc_cl_imputed_pmm - frequency_c_cl)/ss)*100
)

#GRAPH: side-by-side bar chart full data vs polyreg imputed data

# Create tables for frequency

before.fq.cl <- (frequency_c_cl/ss)*100

fq.full.data.cl <- data.frame(
  Category = c("1", "2", "3"),
  Value = before.fq.cl 
)

fq.pmm.xyc.imputed.cl <- (frequency_xyc_cl_imputed_pmm/ss)*100

fq.xyc.pmm.data.cl  <- data.frame(
  Category = c("1", "2", "3"),
  Value = fq.pmm.xyc.imputed.cl
)

fq.xyc.amputed.data.cl <- data.frame(
  Category = c("1", "2", "3"),
  Value = frequency_xyc_missing_cl
)

# Add an identifier column to distinguish between the tables
fq.full.data.cl$Table <- "1. Full data set"
fq.xyc.pmm.data.cl$Table <- "2. Imputed data set"
fq.xyc.amputed.data.cl$Table <- "3. Amputed data set"

# Combine the tables into one data frame
combined_cl.data.full.and.xyc.pmm <- rbind(fq.full.data.cl, fq.xyc.pmm.data.cl, fq.xyc.amputed.data.cl)

# Create the side-by-side bar plot
plot1.cl.xyc.pmm <- ggplot(combined_cl.data.full.and.xyc.pmm, aes(x = Category, y = Value.Freq, fill = Table)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(title = "1. Imputed with pmm method",
       x = "Category",
       y = "Frequency") +
  scale_fill_manual(values = c("1. Full data set"= "skyblue", "2. Imputed data set" = "lightgreen", "3. Amputed data set" = "#1f77b4")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal()

plot1.cl.xyc.pmm



################################################################################
### POLYREG ###

# Convert 'c' to a factor variable 
md.pattern(df_missing_cl_xyc)

df_missing_cl_xyc$c <- as.factor(df_missing_cl_xyc$c)

# Perform multiple imputation using the "polyreg" method
imputed_data_cl_xyc_polyreg <- mice(df_missing_cl_xyc, m = 15, maxit = 20, method = 'polyreg')

# Save the data set
save(imputed_data_cl_xyc_polyreg, file = "imputed_data_cl_xyc_polyreg.RData")

load("imputed_data_cl_xyc_polyreg.RData")
# View the imputed dataset
completed_data_cl_xyc_polyreg <- complete(imputed_data_cl_xyc_polyreg)

md.pattern(completed_data_cl_xyc_polyreg)

frequency_xyc_cl_imputed_polyreg <- table(completed_data_cl_xyc_polyreg$c)
frequency_xyc_cl_imputed_polyreg


#convergence check
cl.xyc.polyreg.con.plot <-plot(imputed_data_cl_xyc_polyreg, y = "c", layout = c(1,2))

# Start the analyses
MIpolyreg.cl.xyc <- with(data = imputed_data_cl_xyc_polyreg , glm(y ~ x + c))
MIpolyreg.cl.xyc

#Rubins rule
pooledEst.polyreg.cl.xyc <- pool(MIpolyreg.cl.xyc)


# Create a data frame with the coefficients and standard errors
coeff_se_xyc_data_cl_polyreg <- data.frame(
  Paramters = pooledEst.polyreg.cl.xyc$pooled["term"],
  Estimate = pooledEst.polyreg.cl.xyc$pooled["estimate"],
  Std_Error = sqrt(pooledEst.polyreg.cl.xyc$pooled["t"]) #t = total variance
)

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.polyreg.xyc <- coeff_se_xyc_data_cl_polyreg[coeff_se_xyc_data_cl_polyreg$term %in% c("c2", "c3"), ]
c2c3.coef.se.polyreg.xyc

# View the imputed dataset
completed_data_cl_xyc_polyreg <- complete(imputed_data_cl_xyc_polyreg)

md.pattern(completed_data_cl_xyc_polyreg)

#frquency moet nog

frequency_xyc_cl_imputed_polyreg <- table(completed_data_cl_xyc_polyreg$c)

xyc_data_imputed_difference_polyreg <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_full_cl/ss)*100,
  After <- (frequency_xyc_cl_imputed_polyreg/ss)*100,
  Imputed <- (( frequency_xyc_cl_imputed_polyreg - frequency_c_cl)/ss)*100
)


#GRAPH: side-by-side bar chart full data vs polyreg imputed data

# Create tables for frequency

before.fq.cl <- (frequency_c_cl/ss)*100

fq.full.data.cl <- data.frame(
  Category = c("1", "2", "3"),
  Value = before.fq.cl 
)

fq.polyreg.xyc.imputed.cl <- (frequency_xyc_cl_imputed_polyreg/ss)*100

fq.xyc.polyreg.data.cl  <- data.frame(
  Category = c("1", "2", "3"),
  Value = fq.polyreg.xyc.imputed.cl
)

# Add an identifier column to distinguish between the tables
fq.full.data.cl$Table <- "1. Full data set"
fq.xyc.polyreg.data.cl$Table <- "2. Imputed data set"
fq.xyc.amputed.data.cl$Table <- "3. Amputed data set"

# Combine the tables into one data frame
combined_cl.data.full.and.xyc.polyreg <- rbind(fq.full.data.cl, fq.xyc.polyreg.data.cl, fq.xc.amputed.data.cl)

# Create the side-by-side bar plot
plot1.cl.xyc.polyreg <- ggplot(combined_cl.data.full.and.xyc.polyreg, aes(x = Category, y = Value.Freq, fill = Table)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(title = "2. Imputed with polyreg method",
       x = "Category",
       y = "Frequency") +
  scale_fill_manual(values = c("1. Full data set"= "skyblue", "2. Imputed data set" = "lightgreen", "3. Amputed data set" = "#1f77b4")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal()

plot1.cl.xyc.polyreg

################################################################################
### POLR ###

md.pattern(df_missing_cl)

# Convert 'c' to a factor variable 
md.pattern(df_missing_cl_xyc)

df_missing_cl_xyc$c <- as.factor(df_missing_cl_xyc$c)

# Perform multiple imputation using the "polr" method
imputed_data_cl_xyc_polr <- mice(df_missing_cl_xyc, m = 15, maxit = 20, method = 'polr')

# Save the data set
save(imputed_data_cl_xyc_polr, file = "imputed_data_cl_xyc_polr.RData")

load("imputed_data_cl_xyc_polr.RData")

# View the imputed dataset
completed_data_cl_xyc_polr <- complete(imputed_data_cl_xyc_polr)

md.pattern(completed_data_cl_xyc_polr)

frequency_xyc_cl_imputed_polr <- table(completed_data_cl_xyc_polr$c)
frequency_xyc_cl_imputed_polr


#convergence check
cl.xyc.polr.con.plot <- plot(imputed_data_cl_xyc_polr, y = "c", layout = c(1,2))

# Start the analyses
MIpolr.cl.xyc <- with(data = imputed_data_cl_xyc_polr , glm(y ~ x + c))
MIpolr.cl.xyc

#Rubins rule
pooledEst.polr.cl.xyc <- pool(MIpolr.cl.xyc)


# Create a data frame with the coefficients and standard errors
coeff_se_xyc_data_cl_polr <- data.frame(
  Paramters = pooledEst.polr.cl.xyc$pooled["term"],
  Estimate = pooledEst.polr.cl.xyc$pooled["estimate"],
  Std_Error = sqrt(pooledEst.polr.cl.xyc$pooled["t"]) #t = total variance
)

# Filter and display only rows where term is "c2" or "c3"
c2c3.coef.se.polr.xyc <- coeff_se_xyc_data_cl_polr[coeff_se_xyc_data_cl_polr$term %in% c("c2", "c3"), ]

# View the imputed dataset
completed_data_cl_xyc_polr <- complete(imputed_data_cl_xyc_polr)

md.pattern(completed_data_cl_xyc_polr)

#frquency moet nog

frequency_xyc_cl_imputed_polr <- table(completed_data_cl_xyc_polr$c)

xyc_data_imputed_difference_polr <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_full_cl/ss)*100,
  After <- (frequency_xyc_cl_imputed_polr/ss)*100,
  Imputed <- (( frequency_xyc_cl_imputed_polr - frequency_c_cl)/ss)*100
)

#GRAPH: side-by-side bar chart full data vs polyreg imputed data

# Create tables for frequency

before.fq.cl <- (frequency_c_cl/ss)*100

fq.full.data.cl <- data.frame(
  Category = c("1", "2", "3"),
  Value = before.fq.cl 
)

fq.polr.xyc.imputed.cl <- (frequency_xyc_cl_imputed_polr/ss)*100

fq.xyc.polr.data.cl  <- data.frame(
  Category = c("1", "2", "3"),
  Value = fq.polr.xyc.imputed.cl
)

# Add an identifier column to distinguish between the tables
fq.full.data.cl$Table <- "1. Full data set"
fq.xyc.polr.data.cl$Table <- "2. Imputed data set"
fq.xyc.amputed.data.cl$Table <- "3. Amputed data set"

# Combine the tables into one data frame
combined_cl.data.full.and.xyc.polr <- rbind(fq.full.data.cl, fq.xyc.polr.data.cl, fq.xyc.amputed.data.cl)

# Create the side-by-side bar plot
plot1.cl.xyc.polr <- ggplot(combined_cl.data.full.and.xyc.polr, aes(x = Category, y = Value.Freq, fill = Table)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(title = "3. Imputed with polr method",
       x = "Category",
       y = "Frequency") +
  scale_fill_manual(values = c("1. Full data set"= "skyblue", "2. Imputed data set" = "lightgreen", "3. Amputed data set" = "#1f77b4")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal()

plot1.cl.xyc.polr
################################################################################
#GRAPHS

### FREQUENCY
#grid the frequency plots
cl.xyc.fqr.plot <- ggarrange(plot1.cl.xyc.pmm, plot1.cl.xyc.polyreg, plot1.cl.xyc.polr, 
          nrow = 3, common.legend = TRUE, legend="bottom")

annotate_figure(cl.xyc.fqr.plot, top = text_grob("Cumulative Logit Model Scenario 2", 
                                      color = "blue", face = "bold", size = 14))

cl.xc.fqr.plot <- ggarrange(plot1.cl.xc.pmm, plot1.cl.xc.polyreg, plot1.cl.xc.polr, 
                    nrow = 3, common.legend = TRUE, legend="bottom")

annotate_figure(cl.xc.fqr.plot, top = text_grob("Cumulative Logit Model Scenario 1", 
                                                 color = "blue", face = "bold", size = 14))

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
coef.plot1.xyc <- create_coefficient_point_plot(coeff_se_xyc_data_cl_pmm, "PMM")
coef.plot2.xyc <- create_coefficient_point_plot(coeff_se_xyc_data_cl_polyreg, "Polyreg")
coef.plot3.xyc <- create_coefficient_point_plot(coeff_se_xyc_data_cl_polr, "POLR")

# Arrange all the graphs in a 2x2 grid
grid.arrange(coef.plot1.xyc, coef.plot2.xyc, coef.plot3.xyc, ncol = 3)

cl.xc.pmm.con.plot
cl.xc.polyreg.con.plot
cl.xc.polr.con.plot
cl.xyc.pmm.con.plot
cl.xyc.polyreg.con.plot
cl.xyc.polr.con.plot


