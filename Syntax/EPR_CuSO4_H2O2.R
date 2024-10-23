# Load libraries
library(ggplot2)
library(pracma)  # For numerical integration
library(gridExtra)  # For arranging plots in a panel
library(readxl)  
library(ggpmisc)

# 2Cu2++H2O2→ 2Cu+ + 2H+ + O2
# 2H2O2→ 2H2O + O2
# 2Cu2++H2O2+2OH−→2Cu+ + O2 + 2H2O

### Load data and pre process ----

# Read the EPR data, skipping lines if necessary
epr_data <- read.table("Data/CUSO4_H2O2/CUSO4_H2O2_7.TXT", header = FALSE, skip = 2)

# Assign meaningful column names
colnames(epr_data) <- c("Index", "Magnetic_Field", "Signal_Intensity")

# Check the structure of the data
str(epr_data)

# Plotting the raw EPR spectrum
ggplot(epr_data, aes(x = Magnetic_Field, y = Signal_Intensity)) +
  geom_line() +
  labs(title = "EPR Spectrum", x = "Magnetic Field (Gauss)", y = "Signal Intensity")


### Smoothing data ---- 

# "If you have a big enough range of data that is mostly 0 with a flat-ish baseline  _Stan"

# Step 1: Smooth the signal using a moving average filter (or other methods)
epr_data$Smoothed_Signal <- movavg(epr_data$Signal_Intensity, n = 5, type = "s")

# Step 2: Subtract the median to adjust the baseline
median_value <- median(epr_data$Smoothed_Signal)
epr_data$Baseline_Adjusted_Signal <- epr_data$Smoothed_Signal - median_value

# Step 3: Plot the adjusted signal
ggplot(epr_data, aes(x = Magnetic_Field, y = Baseline_Adjusted_Signal)) +
  geom_line() +
  labs(title = "Baseline-Adjusted EPR Spectrum", x = "Magnetic Field (Gauss)", y = "Signal Intensity (Adjusted)")

# Step 4: Optional: Focus on a specific range of the magnetic field 
subset_data <- subset(epr_data, Magnetic_Field >= 2800 & Magnetic_Field <= 3500)

# Plot the subset of the data
ggplot(subset_data, aes(x = Magnetic_Field, y = Baseline_Adjusted_Signal)) +
  geom_line() +
  labs(title = "Subset of Baseline-Adjusted EPR Spectrum", x = "Magnetic Field (Gauss)", y = "Signal Intensity (Adjusted)")



### Calculate the g-Factor for Each Peak ----

# Step 1: Subset the data to the magnetic field range of 2800 to 3500 Gauss
subset_data <- subset(epr_data, Magnetic_Field >= 2800 & Magnetic_Field <= 3500)

# Step 2: Find the maxima and minima points in the corrected signal
max_index <- which.max(subset_data$Baseline_Adjusted_Signal)
min_index <- which.min(subset_data$Baseline_Adjusted_Signal)

# Magnetic field values at the maxima and minima signal intensities
max_field <- subset_data$Magnetic_Field[max_index]
min_field <- subset_data$Magnetic_Field[min_index]

# Step 3: Mark the maxima and minima on the plot
ggplot(subset_data, aes(x = Magnetic_Field, y = Baseline_Adjusted_Signal)) +
  geom_line() +
  geom_point(aes(x = max_field, y = subset_data$Baseline_Adjusted_Signal[max_index]), color = "blue", size = 3) +
  geom_point(aes(x = min_field, y = subset_data$Baseline_Adjusted_Signal[min_index]), color = "red", size = 3) +
  labs(title = "Subset of Baseline-Adjusted EPR Spectrum with Maxima and Minima",
       x = "Magnetic Field (Gauss)", 
       y = "Signal Intensity (Adjusted)") +
  annotate("text", x = max_field, y = subset_data$Baseline_Adjusted_Signal[max_index], label = paste("Max:", round(max_field, 2)), vjust = -1) +
  annotate("text", x = min_field, y = subset_data$Baseline_Adjusted_Signal[min_index], label = paste("Min:", round(min_field, 2)), vjust = 1)

# Step 4: Calculate the g-Factor for Each Peak
# Midpoint of the magnetic field values
midpoint_field <- (max_field + min_field) / 2

# Convert the midpoint magnetic field from Gauss to Tesla
midpoint_field_T <- midpoint_field * 1e-4  # 1 G = 1e-4 T

# Constants for g-factor calculation
microwave_frequency <- 9.647667e+09  # 9.647667 GHz in Hz
h <- 6.626e-34  # Planck's constant in J·s
mu_B <- 9.274e-24  # Bohr magneton in J/T

# Step 5: Calculate the g-factor using the midpoint magnetic field
g_factor <- h * microwave_frequency / (mu_B * midpoint_field_T)

# Output the calculated g-factor
cat("Calculated g-factor:", g_factor, "\n")



### Calculating the absorbance value ----

# Step 1: Subset the data to the magnetic field range of 2800 to 3500 Gauss
subset_data <- subset(epr_data, Magnetic_Field >= 2800 & Magnetic_Field <= 3500)

# Step 2: Calculate the anti-derivative (cumulative integral) for absorbance
# We use the cumulative trapezoidal integration for this purpose
subset_data$Absorbance <- cumtrapz(subset_data$Magnetic_Field, subset_data$Baseline_Adjusted_Signal)

# Step 3: Plot the Absorbance curve
ggplot(subset_data, aes(x = Magnetic_Field, y = Absorbance)) +
  geom_line() +
  labs(title = "Absorbance Curve (Cumulative Integral)", x = "Magnetic Field (Gauss)", y = "Absorbance")

# Step 4: Calculate the total absorbance (area under the curve)
# Using the trapezoidal integration for the absorbance curve
total_absorbance <- trapz(subset_data$Magnetic_Field, subset_data$Absorbance)

# Output the total absorbance value
cat("Total Absorbance (Area Under the Curve):", total_absorbance, "\n")



### Plotting ----

# Step 1: Subset the data to the magnetic field range of 2800 to 3500 Gauss
subset_data <- subset(epr_data, Magnetic_Field >= 2800 & Magnetic_Field <= 3500)

# Step 2: Calculate the anti-derivative (cumulative integral) for absorbance
subset_data$Absorbance <- cumtrapz(subset_data$Magnetic_Field, subset_data$Baseline_Adjusted_Signal)

# Step 3: Create the EPR spectrum plot (top plot)
epr_plot <- ggplot(subset_data, aes(x = Magnetic_Field, y = Baseline_Adjusted_Signal)) +
  geom_line() +
  labs(title = "EPR Spectrum", x = "Magnetic Field (Gauss)", y = "Signal Intensity (Adjusted)") +
  theme(panel.grid.major = element_line(linetype = "dotted", colour = "black"))

# Step 4: Create the Absorbance curve plot (bottom plot)
# Shade the area under the curve with yellow
absorbance_plot <- ggplot(subset_data, aes(x = Magnetic_Field, y = Absorbance)) +
  geom_line() +
  geom_ribbon(aes(ymin = 0, ymax = Absorbance), fill = "yellow", alpha = 0.5) +
  labs(title = "Absorbance Curve (Cumulative Integral)", x = "Magnetic Field (Gauss)", y = "Absorbance") +
  theme(panel.grid.major = element_line(linetype = "dotted", colour = "black"))

# Step 5: Arrange both plots in one panel (EPR on top, Absorbance below)
grid.arrange(epr_plot, absorbance_plot, ncol = 1)





### Summary CuSo4_H2O2 ----

# Read the Excel file (assuming data is in the first sheet)
data2 <- read_excel("Data/CUSO4_T2/Summary_CUSO4.xlsx", sheet = "CuSo4_H2O2")

head(data)

# Create a scatter plot of Absorbance (AUC) vs Conc. (mM)
ggplot(data2, aes(x = `Conc. (mM)`, y = `Absorbance (AUC)`)) +
  geom_point() +
  geom_line() +
  labs(x = "Conc. (mM)", y = "EPR | Absorbance (AUC)", title = "Absorbance vs Conc.") +
  theme_minimal()


#### correlation and calculate the R-squared (R²) ----

# Fit a linear regression model: Absorbance (AUC) vs Concentration (mM)
linear_model <- lm(`Absorbance (AUC)` ~ `Conc. (mM)`, data = data2)

# Get the summary of the linear model to extract R-squared value
summary_model <- summary(linear_model)
r_squared <- summary_model$r.squared

# Print the R-squared value
print(paste("R-squared value:", r_squared))

# Plot the data with the fitted regression line
ggplot(data2, aes(x = `Conc. (mM)`, y = `Absorbance (AUC)`)) +
  geom_point() +                                    # Scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear regression line
  labs(x = "Conc. (mM)", y = "EPR | Absorbance (AUC)", title = "Absorbance vs Conc.") +
  theme(panel.grid.major = element_line(linetype = "dotted", colour = "black"))




# Create the base plot
plot <- ggplot(data2, aes(x = `Conc. (mM)`, y = `Absorbance (AUC)`)) +
  geom_point() +                                    # Scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear regression line
  labs(x = "Conc. (mM)", y = "EPR | Absorbance (AUC)", title = "Absorbance vs Conc.") +
  theme(panel.grid.major = element_line(linetype = "dotted", colour = "black"))

# Add the R-squared value as an annotation on the plot
plot + annotate("text", x = max(data2$`mols`), y = min(data2$`Absorbance (AUC)`),
                label = paste("R² =", round(r_squared, 4)), hjust = 1, vjust = 0)

# Add the R-squared value as an annotation around the coordinates x = 150, y = 20000
plot + annotate("text", 
                x = 100, 
                y = 20000,
                label = paste("R² =", round(r_squared, 4)),
                hjust = 0, vjust = 0)





### Summary CuSO4 &  CuSO4+H2O2 ----

# Read the Excel file (assuming data is in the first sheet)
data3 <- read_excel("Data/CUSO4_T2/Summary_CUSO4.xlsx", sheet = "Combined")
head(data3)

# Split data by Type
data3_CuSO4 <- subset(data3, Type == "CuSO4")
data3_CuSO4_H2O2 <- subset(data3, Type == "CuSO4+H2O2")

# Calculate R-squared values
lm_CuSO4 <- lm(`Absorbance (AUC)` ~ `Conc. (mM)`, data = data3_CuSO4)
lm_CuSO4_H2O2 <- lm(`Absorbance (AUC)` ~ `Conc. (mM)`, data = data3_CuSO4_H2O2)

r_squared_CuSO4 <- summary(lm_CuSO4)$r.squared
r_squared_CuSO4_H2O2 <- summary(lm_CuSO4_H2O2)$r.squared




# Create the plot
ggplot(data3, aes(x = `Conc. (mM)`, y = `Absorbance (AUC)`, color = Type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Absorbance (AUC) vs. Concentration (mM) for CuSO4 and CuSO4+H2O2",
       x = "Concentration (mM)", y = "Absorbance (AUC)") +
  theme_minimal()






# Add R Sq Plots
ggplot(data = data3, aes(x = `Conc. (mM)`, y = `Absorbance (AUC)`, color = Type)) +
  geom_point() + 
  stat_poly_line() + 
  stat_poly_eq() +
  labs(title = "Absorbance (AUC) vs. Concentration (mM) for CuSO₄ and CuSO₄+H₂O₂",
       x = "Concentration (mM)", y = "Absorbance (AUC)") +
  theme_minimal()




# Fit linear models for both CuSO4 and CuSO4+H2O2
lm_CuSO4 <- lm(`Absorbance (AUC)` ~ `Conc. (mM)`, data = subset(data3, Type == "CuSO4"))
lm_CuSO4_H2O2 <- lm(`Absorbance (AUC)` ~ `Conc. (mM)`, data = subset(data3, Type == "CuSO4+H2O2"))

# Extract slopes (coefficients)
slope_CuSO4 <- coef(lm_CuSO4)[2]  # The second coefficient is the slope
slope_CuSO4_H2O2 <- coef(lm_CuSO4_H2O2)[2]  # The second coefficient is the slope

# Print the slopes
slope_CuSO4
slope_CuSO4_H2O2

# Create slope labels
slope_label_CuSO4 <- paste("CuSO4 Slope = ", round(slope_CuSO4, 3))
slope_label_CuSO4_H2O2 <- paste("CuSO4+H2O2 Slope = ", round(slope_CuSO4_H2O2, 3))

# Plot with slopes
ggplot(data3, aes(x = `Conc. (mM)`, y = `Absorbance (AUC)`, color = Type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(x = 75, y = 250000, label = slope_label_CuSO4, color = "blue") +
  geom_text(x = 75, y = 150000, label = slope_label_CuSO4_H2O2, color = "green") +
  labs(title = "Absorbance (AUC) vs. Concentration (mM) for CuSO4 and CuSO4+H2O2",
       x = "Concentration (mM)", y = "Absorbance (AUC)") +
  theme_minimal()



### References ---- 

# https://webhome.auburn.edu/~duinedu/epr/2_pracaspects.pdf
