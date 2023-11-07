#Packages
install.packages("mice")
install.packages("ggmice")
install.packages("ggpubr")
library(ggpubr)
library(mice)
library(ggplot2)
library(ggmice)
library(gridExtra)


set.seed(1234567)

setwd("/Users/vc/Documents/Studie/MSc Epi/Research/Code")


################################################################################
################################################################################
#CL ampute
load("/Users/vc/Documents/Studie/MSc Epi/Research/Code/data_cl.RData")

# Save the Y variable from the  CL original data frame to add later 
y_cl <- data_cl$y

#Frequency c in full data prior to ampute
frequency_c_cl <- table(data_cl$c)
frequency_c_cl

#set sammple size for frequency calculation in %
ss <- 100000

#subsetting data to show only y, x, and c
sub.data_cl <- subset(data_cl, select = c("y", "c", "x"))

model_full_cl_data <- glm(y ~ x + c, data = sub.data_cl)

# Extract the coefficients and standard errors
coef_full_data_cl <- coef(model_full_cl_data)
se_full_data_cl <- summary(model_full_cl_data)$coef[, "Std. Error"]

# Create a data frame with the coefficients and standard errors
coeff_se_full_data_cl <- data.frame(
  Coefficient = names(coef_full_data_cl),
  Estimate = coef_full_data_cl,
  Std_Error = se_full_data_cl
)

summary(model_full_cl_data)


# Create a point diagram of the coefficients, no se
plot.coef.full.data.cl <- ggplot(coeff_se_full_data_cl, aes(x = Coefficient, y = Estimate)) +
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
missing_cl <- ampute(
  data_cl,
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
df_missing_cl <- as.data.frame(missing_cl$amp)

save(df_missing_cl, file = "df_missing_cl.RData")

load("df_missing_cl.RData")
#Show only data y, x, and c

sub.df_missing_cl <- subset(df_missing_cl, select = c("y", "c", "x"))

#Missing pattern
plot_pattern(sub.df_missing_cl)

#frequency of missingness
freq_missing_cl <- table(sub.df_missing_cl$c)
freq_missing_cl #automatisering

#Table to calculte the difference prior ampute vs post ampute
full_data_difference_cl <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_cl/ss)*100,
  After <- (freq_missing_cl/ss)*100,
  Amputed <- ((frequency_c_cl - freq_missing_cl)/ss)*100
)

#regression coefficients
model_full_cl_aa <- glm(y ~ x + c, data = sub.df_missing_cl)

summary(model_full_cl_aa)

# Extract the coefficients and standard errors
coef_full_data_cl_aa <- coef(model_full_cl_aa )
se_full_data_cl_aa <- summary(model_full_cl_aa )$coef[, "Std. Error"]

# Create a data frame with the coefficients and standard errors
coeff_se_full_data_cl_aa <- data.frame(
  Coefficient = names(coef_full_data_cl_aa),
  Estimate = coef_full_data_cl_aa,
  Std_Error = se_full_data_cl_aa
)

# Create a point diagram of the coefficients, no se
plot.coef.full.data.cl.aa <- ggplot(coeff_se_full_data_cl_aa, aes(x = Coefficient, y = Estimate)) +
  geom_point() +
  labs(
    title = "Coefficients of Linear Regression Model",
    x = "Coefficient",
    y = "Estimate"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



################################################################################

# x and c only

# Subset the dataset to include only "x" and "c" variables
data_cl_x_c <- data_cl[, c("x", "c")]
data_cl_x_c$c <- factor(data_cl_x_c$c) 

# Introduce missingness using ampute 
missing_cl_x_c <- ampute(
  data_cl_x_c,
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

#Make amputed data frame and save as RData file
df_missing_cl_x_c <- as.data.frame(missing_cl_x_c$amp)

# Add the Y variable back to the data frame 
df_missing_cl_x_c$y <- y_cl

#save and check
save(df_missing_cl_x_c, file = "df_missing_cl_x_c.RData")

load("df_missing_cl_x_c.RData")

#Missing pattern
plot_pattern(df_missing_cl_x_c)

#frequency of missingness
freq_missing_cl_x_c <- table(df_missing_cl_x_c$c)
freq_missing_cl_x_c

#Table to calculte the difference prior ampute vs post ampute
missing_data_difference_cl_xc <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_cl/ss)*100,
  After <- (freq_missing_cl_x_c/ss)*100,
  Amputed <- ((frequency_c_cl - freq_missing_cl_x_c)/ss)*100
)


#regression coefficients
model_cl_aa_xc <- glm(y ~ x + c, data = df_missing_cl_x_c)

summary(model_cl_aa_xc)

# Extract the coefficients and standard errors
coef_data_cl_aa_xc <- coef(model_cl_aa_xc)
se_data_cl_aa_xc <- summary(model_cl_aa_xc)$coef[, "Std. Error"]

# Create a data frame with the coefficients and standard errors
coeff_se_data_cl_xc <- data.frame(
  Coefficient = names(coef_data_cl_aa_xc),
  Estimate = coef_data_cl_aa_xc,
  Std_Error = se_data_cl_aa_xc
)

#plot coefficients
plot.coef.xc.data.cl <- ggplot(coeff_se_data_cl_xc, aes(x = Coefficient, y = Estimate)) +
  geom_point() +
  labs(
    title = "Coefficients of Linear Regression Model",
    x = "Coefficient",
    y = "Estimate"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################################################################################

# x, y,  and c only

# Subset the dataset to include only "x" and "c" variables
data_cl_xyc <- data_cl[, c("x", "c", "y")]
data_cl_xyc$c <- factor(data_cl_xyc$c) 

# Introduce missingness using ampute 
missing_cl_xyc <- ampute(
  data_cl_xyc,
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

#Make amputed data frame and save as RData file
df_missing_cl_xyc <- as.data.frame(missing_cl_xyc$amp)

save(df_missing_cl_xyc, file = "df_missing_cl_xyc.RData")

load("df_missing_cl_xyc.RData")
#Missing pattern
plot_pattern(df_missing_cl_xyc)

#frequency of missingness
freq_missing_cl_xyc <- table(df_missing_cl_xyc$c)
freq_missing_cl_xyc


#Table to calculte the difference prior ampute vs post ampute
missing_data_difference_cl_xyc <- data.frame(
  Category = c("1", "2", "3"),
  Before <- (frequency_c_cl/ss)*100,
  After <- (freq_missing_cl_xyc/ss)*100,
  Amputed <- ((frequency_c_cl - freq_missing_cl_xyc)/ss)*100
)


#regression coefficients
model_cl_aa_xyc <- glm(y ~ x + c, data = df_missing_cl_xyc)

summary(model_cl_aa_xyc)

# Extract the coefficients and standard errors
coef_data_cl_aa_xyc <- coef(model_cl_aa_xyc)
se_data_cl_aa_xyc <- summary(model_cl_aa_xyc)$coef[, "Std. Error"]

# Create a data frame with the coefficients and standard errors
coeff_se_data_cl_xyc <- data.frame(
  Coefficient = names(coef_data_cl_aa_xyc),
  Estimate = coef_data_cl_aa_xyc,
  Std_Error = se_data_cl_aa_xyc
)

#plot coefficients
plot.coef.xyc.data.cl <- ggplot(coeff_se_data_cl_xyc, aes(x = Coefficient, y = Estimate)) +
  geom_point() +
  labs(
    title = "Coefficients of Linear Regression Model",
    x = "Coefficient",
    y = "Estimate"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### All plots for CL in one grid

create_bar_chart_cl <- function(data, title) {
  ggplot(data, aes(x = Category)) +
    geom_bar(aes(y = Before, fill = "Before"), stat = "identity", position = "dodge") +
    geom_bar(aes(y = After, fill = "After"), stat = "identity", position = "dodge") +
    geom_bar(aes(x = as.numeric(Category) + 0.2, y = Amputed, fill = "Amputed"), stat = "identity", position = "dodge") +
    labs(title = title,
         x = "Category",
         y = "Value") +
    scale_fill_manual(values = c("Before" = "skyblue", "After" = "lightblue", "Amputed" = "#1f77b4")) +
    theme_minimal()
}


# Create the plots
plot4.cl <- create_bar_chart_cl(full_data_difference_cl, "Amputation for Full data")
plot5.cl <- create_bar_chart_cl(missing_data_difference_cl_xc, "Amputation for C and X relation")
plot6.cl <- create_bar_chart_cl(missing_data_difference_cl_xyc, "Amputation for C, X, and Y relation")

# Arrange the plots
grid.arrange(plot4.cl, plot5.cl, plot6.cl, ncol = 3)

ggarrange(plot4.cl, plot5.cl, plot6.cl, ncol=3, common.legend = TRUE, legend="bottom")



### ALL three data plots in one grid

# Function to create the bar chart
create_bar_chart_cl <- function(data, title) {
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
plot1.cl <- create_bar_chart_cl(full_data_difference_c, "Amputation for Full data")
plot2.cl <- create_bar_chart_cl(full_data_difference_cx, "Amputation under CL S1")
plot3.cl <- create_bar_chart_cl(full_data_difference_xyc, "Amputation under CL S2")



#grid coefficents graphs
ggarrange(plot.coef.xc.data.cl, plot.coef.xyc.data.cl, ncol=2, common.legend = TRUE, legend="bottom")

ggarrange(plot2.cl, plot3.cl, ncol=2, common.legend = TRUE, legend="bottom")
