
#Packages
install.packages("mice")
install.packages("ggmice")
install.packages("ggpubr")
library(ggpubr)
library(mice)
library(ggplot2)
library(ggmice)
library(gridExtra)


#set seed + working directory 
set.seed(1234567)

setwd("/Users/vc/Documents/Studie/MSc Epi/Research/Code")

################################################################################
################################################################################

# mlr.n ampute 

load("data_mlr.n.RData")

# Save the Y variable from the  mlr.n original data frame to add later 
y_mlr.n <- data_mlr.n$y

#Frequency c in full data prior to ampute
frequency_c_mlr.n <- table(data_mlr.n$c)
frequency_c_mlr.n

#sample size, for turning the frequencies into %

ss <- 100000

#subsetting data to show only y, x, and c
sub.data_mlr.n <- subset(data_mlr.n, select = c("y", "c", "x"))

model_full_mlr.n_data <- glm(y ~ x + c, data = sub.data_mlr.n)

# Extract the coefficients and standard errors
coef_full_data_mlr.n <- coef(model_full_mlr.n_data)
se_full_data_mlr.n <- summary(model_full_mlr.n_data)$coef[, "Std. Error"]

# Create a data frame with the coefficients and standard errors
coeff_se_full_data_mlr.n <- data.frame(
  Coefficient = names(coef_full_data_mlr.n),
  Estimate = coef_full_data_mlr.n,
  Std_Error = se_full_data_mlr.n
)

summary(model_full_mlr.n_data)


# Create a point diagram of the coefficients, no se
plot.coef.full.data.mlr.n <- ggplot(coeff_se_full_data_mlr.n, aes(x = Coefficient, y = Estimate)) +
  geom_point() +
  labs(
    title = "Coefficients of Linear Regression Model",
    x = "Coefficient",
    y = "Estimate"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################################################################################
# Full data

# Introduce missingness using ampute full data
missing_mlr.n <- ampute(
  data_mlr.n,
  prop = 0.15,
  patterns = c(1, 0, 1, 1, 1, 1),
  mech = "MAR",
  weights = NULL,
  std = TRUE,
  cont = TRUE,
  type = NULL,
  odds = NULL,
  bycases = TRUE,
  run = TRUE
)

#Make amputed data frame and save as RData file
df_missing_mlr.n <- as.data.frame(missing_mlr.n$amp)

save(df_missing_mlr.n, file = "df_missing_mlr.n.RData")

load("df_missing_mlr.n.RData")

#Show only data y, x, and c

sub.df_missing_mlr.n <- subset(df_missing_mlr.n, select = c("y", "c", "x"))
  
#Missing pattern
plot_pattern(sub.df_missing_mlr.n)

#frequency of missingness
freq_missing_mlr.n <- table(sub.df_missing_mlr.n$c)

#Table to calculte the difference prior ampute vs post ampute
full_data_difference_c <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_mlr.n/ss)*100,
  After <- (freq_missing_mlr.n/ss)*100,
  Amputed <- ((frequency_c_mlr.n - freq_missing_mlr.n)/ss)*100
)

#regression coefficients
model_full_mlr.n_aa <- glm(y ~ x + c, data = sub.df_missing_mlr.n)

summary(model_full_mlr.n_aa)

# Extract the coefficients and standard errors
coef_full_data_mlr.n_aa <- coef(model_full_mlr.n_aa )
se_full_data_mlr.n_aa <- summary(model_full_mlr.n_aa )$coef[, "Std. Error"]

# Create a data frame with the coefficients and standard errors
coeff_se_full_data_mlr.n_aa <- data.frame(
  Coefficient = names(coef_full_data_mlr.n_aa),
  Estimate = coef_full_data_mlr.n_aa,
  Std_Error = se_full_data_mlr.n_aa
)



################################################################################
# x and c only
# Subset the dataset to include only "x" and "c" variables
data_mlr.n_x_c <- data_mlr.n[, c("x", "c")]
data_mlr.n_x_c$c <- factor(data_mlr.n_x_c$c) 

# Convert the subset into a matrix and c to categorical
df_mlr.n_x_c <- as.data.frame(data_mlr.n_x_c)
df_mlr.n_x_c$c <- factor(df_mlr.n_x_c$c)

# Introduce missingness using ampute x and c
missing_mlr.n_x_c <- ampute(
  data_mlr.n_x_c,
  prop = 0.15,
  patterns = c(1, 0),
  mech = "MAR",
  weights = NULL,
  std = TRUE,
  cont = TRUE,
  type = NULL,
  odds = NULL,
  bycases = TRUE,
  run = TRUE
)

#maak c categorical 
#Save data
df_missing_mlr.n_x_c <- as.data.frame(missing_mlr.n_x_c$amp)

# Add the Y variable back to the data frame 
df_missing_mlr.n_x_c$y <- y_mlr.n

#save as file + load to check
save(df_missing_mlr.n_x_c, file = "df_missing_mlr.n_x_c.RData")

load("df_missing_mlr.n_x_c.RData")

#pattern missingness
plot_pattern(df_missing_mlr.n_x_c)

#point graph
ggmice(df_missing_mlr.n_x_c, aes(x, c)) +
  geom_point()

freq_missing_mlr.n_xc <- table(df_missing_mlr.n_x_c$c)

#Table to calculte the difference prior ampute vs post ampute
full_data_difference_cx <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_mlr.n/ss)*100,
  After <- (freq_missing_mlr.n/ss)*100,
  Amputed <- ((frequency_c_mlr.n - freq_missing_mlr.n_xc)/ss)*100
)

#regression coefficients
model_mlr.n_aa_xc <- glm(y ~ x + c, data = df_missing_mlr.n_x_c)

summary(model_mlr.n_aa_xc)

# Extract the coefficients and standard errors
coef_data_mlr.n_aa_xc <- coef(model_mlr.n_aa_xc)
se_data_mlr.n_aa_xc <- summary(model_mlr.n_aa_xc)$coef[, "Std. Error"]

# Create a data frame with the coefficients and standard errors
coeff_se_data_mlr.n_xc <- data.frame(
  Coefficient = names(coef_data_mlr.n_aa_xc),
  Estimate = coef_data_mlr.n_aa_xc,
  Std_Error = se_data_mlr.n_aa_xc
)

#plot coefficients
plot.coef.xc.data.mlr.n <- ggplot(coeff_se_data_mlr.n_xc, aes(x = Coefficient, y = Estimate)) +
  geom_point() +
  labs(
    title = "Coefficients of Linear Regression Model",
    x = "Coefficient",
    y = "Estimate"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################################################################################
# x, y, and c only
# Subset the dataset to include only "x", "c", and "y" variables
data_mlr.n_xyc <- data_mlr.n[, c("x", "c", "y")]
data_mlr.n_xyc$c <- factor(data_mlr.n_xyc$c)


# Introduce missingness using ampute x, c, and y
missing_mlr.n_xyc <- ampute(
  data_mlr.n_xyc,
  prop = 0.15,
  patterns = c(1, 0, 1),
  mech = "MAR",
  weights = NULL,
  std = TRUE,
  cont = TRUE,
  type = NULL,
  odds = NULL,
  bycases = TRUE,
  run = TRUE
)

#Save as RData
df_missing_mlr.n_xyc <- as.data.frame(missing_mlr.n_xyc$amp)
save(df_missing_mlr.n_xyc, file = "df_missing_mlr.n_xyc.RData")

load("df_missing_mlr.n_xyc.RData")
#point graph
ggmice(df_missing_mlr.n_xyc, aes(x, c)) +
  geom_point()


freq_missing_mlr.n_xyc <- table(df_missing_mlr.n_xyc$c)

#Table to calculte the difference prior ampute vs post ampute
full_data_difference_xyc <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_mlr.n/ss)*100,
  After <- (freq_missing_mlr.n/ss)*100,
  Amputed <- ((frequency_c_mlr.n - freq_missing_mlr.n_xyc)/ss)*100
)

#regression coefficients

model_mlr.n_aa_xyc <- glm(y ~ x + c, data = df_missing_mlr.n_xyc)
summary(model_mlr.n_aa_xyc)

# Extract the coefficients and standard errors
coef_data_mlr.n_aa_xyc <- coef(model_mlr.n_aa_xyc)
se_data_mlr.n_aa_xyc <- summary(model_mlr.n_aa_xyc)$coef[, "Std. Error"]

# Create a data frame with the coefficients and standard errors
coeff_se_data_mlr.n_xyc <- data.frame(
  Coefficient = names(coef_data_mlr.n_aa_xyc),
  Estimate = coef_data_mlr.n_aa_xyc,
  Std_Error = se_data_mlr.n_aa_xyc
)

#plot coefficients
plot.coef.xyc.data.mlr.n <- ggplot(coeff_se_data_mlr.n_xyc, aes(x = Coefficient, y = Estimate)) +
  geom_point() +
  labs(
    title = "Coefficients of Linear Regression Model",
    x = "Coefficient",
    y = "Estimate"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###### PLOT all coef in grid

# Create point graphs for each set of coefficients
plot_full_data_mlr.n_aa <- ggplot(coeff_se_full_data_mlr.n, aes(x = Coefficient, y = Estimate)) +
  geom_point() +
  labs(
    title = "Full data after ampute",
    x = "Coefficient",
    y = "Estimate"
  )

plot_mlr.n_aa_xc <- ggplot(coeff_se_data_mlr.n_xc, aes(x = Coefficient, y = Estimate)) +
  geom_point() +
  labs(
    title = "Data after ampute X and C",
    x = "Coefficient",
    y = "Estimate"
  )

plot_mlr.n_aa_xyc <- ggplot(coeff_se_data_mlr.n_xyc, aes(x = Coefficient, y = Estimate)) +
  geom_point() +
  labs(
    title = "Data after ampute X, Y and C",
    x = "Coefficient",
    y = "Estimate"
  )

# Create a point graph for the new coefficients
plot_full_data_mlr.n_before <- ggplot(coeff_se_full_data_mlr.n, aes(x = Coefficient, y = Estimate)) +
  geom_point() +
  labs(
    title = "Coefficients Before ampute",
    x = "Coefficient",
    y = "Estimate"
  )


# Arrange all the graphs in a 2x2 grid
grid.arrange(plot_full_data_mlr.n_before, plot_full_data_mlr.n_aa, plot_mlr.n_aa_xc, plot_mlr.n_aa_xyc, ncol = 2)

# Arrange two scenarios only the graphs in a 2x2 grid
grid.arrange( plot_mlr.n_aa_xc, plot_mlr.n_aa_xyc, ncol = 2)

### ALL three data plots in one grid

# Function to create the bar chart
create_bar_chart_mlr.n <- function(data, title) {
  ggplot(data, aes(x = Category)) +
    geom_bar(aes(y = Before, fill = "Before"), stat = "identity", position = "dodge") +
    geom_bar(aes(y = After, fill = "After"), stat = "identity", position = "dodge") +
    geom_bar(aes(x = as.numeric(Category) + 0.2, y = Amputed, fill = "Amputed"), stat = "identity", position = "dodge") +
    labs(title = title,
         x = "Category",
         y = "") +
    scale_fill_manual(values = c("Before" = "skyblue", "After" = "lightblue", "Amputed" = "#1f77b4")) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Specify percent scale for y-axis
    theme_minimal()
}


# Create the plots
plot1 <- create_bar_chart_mlr.n(full_data_difference_c, "Amputation for full data")
plot2 <- create_bar_chart_mlr.n(full_data_difference_cx, "Amputation under MLR S1")
plot3 <- create_bar_chart_mlr.n(full_data_difference_xyc, "Amputation under MLR S2")


#grid the frequency plots
ggarrange(plot1, plot2, plot3, ncol=3, common.legend = TRUE, legend="bottom")

#grid the frequency plots
ggarrange(plot2, plot3, ncol=2, common.legend = TRUE, legend="bottom")

#grid coefficents graphs
ggarrange(plot_mlr.n_aa_xc, plot_mlr.n_aa_xyc, ncol=2, common.legend = TRUE, legend="bottom")
