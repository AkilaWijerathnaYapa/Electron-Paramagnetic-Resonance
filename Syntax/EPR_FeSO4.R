# Load libraries
library(ggplot2)
library(pracma)

# Read the EPR data, skipping lines if necessary
epr_data <- read.table("Data/FESO4_T2/FESO4_D0.TXT", header = FALSE, skip = 2)

# Assign meaningful column names
colnames(epr_data) <- c("Index", "Magnetic_Field", "Signal_Intensity")

# Check the structure of the data
str(epr_data)

# Plotting the raw EPR spectrum
ggplot(epr_data, aes(x = Magnetic_Field, y = Signal_Intensity)) +
  geom_line() +
  labs(title = "EPR Spectrum", x = "Magnetic Field (Gauss)", y = "Signal Intensity")


# Plotting the raw EPR spectrum with x-axis labels every 250 units and rotated 90 degrees
ggplot(epr_data, aes(x = Magnetic_Field, y = Signal_Intensity)) +
  geom_line() +
  labs(title = "EPR Spectrum", x = "Magnetic Field (Gauss)", y = "Signal Intensity") +
  scale_x_continuous(breaks = seq(min(epr_data$Magnetic_Field), max(epr_data$Magnetic_Field), by = 250)) +  # Set x-axis breaks every 250 units
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels 90 degrees



#################################################################
# below sybtax are incorrect for thiese peaks / 

# mark maxima nad minima


###--- Calculate the g-Factor for Each Peak

# Step 1: Find the maximum and minimum points in the corrected signal
max_index <- which.max(epr_data$Signal_Intensity)
min_index <- which.min(epr_data$Signal_Intensity)

# Magnetic field values at the maximum and minimum signal intensities
max_field <- epr_data$Magnetic_Field[max_index]
min_field <- epr_data$Magnetic_Field[min_index]

# Step 2: Calculate the midpoint of the magnetic field values
midpoint_field <- (max_field + min_field) / 2

# Step 3: Convert the midpoint magnetic field from Gauss to Tesla
midpoint_field_T <- midpoint_field * 1e-4  # 1 G = 1e-4 T

# Constants for g-factor calculation
microwave_frequency <- 9.647667e+09 # 9.648322 GHz in Hz
h <- 6.626e-34  # Planck's constant in J·s
mu_B <- 9.274e-24  # Bohr magneton in J/T

# Step 4: Calculate the g-factor using the midpoint magnetic field
g_factor <- h * microwave_frequency / (mu_B * midpoint_field_T)


# Ref https://webhome.auburn.edu/~duinedu/epr/2_pracaspects.pdf



###--- anti-derivation

# Step 1: Define the range of interest
start_field <- 2897.1972
end_field <- 3897.1972

# Step 2: Subset the data based on the defined range
roi_peak <- subset(epr_data, Magnetic_Field >= start_field & Magnetic_Field <= end_field)

# Step 3: Apply a smoothing filter to the signal (optional but recommended)
spline_smooth <- smooth.spline(roi_peak$Magnetic_Field, roi_peak$Signal_Intensity, spar=0.7)
smoothed_signal <- predict(spline_smooth, roi_peak$Magnetic_Field)$y

# Step 4: Implement cumulative trapezoidal rule manually for cumulative integration
cumulative_absorption <- cumsum((smoothed_signal[-1] + smoothed_signal[-length(smoothed_signal)]) / 2 * diff(roi_peak$Magnetic_Field))

# Step 5: Add the first zero value (so that it starts at zero)
cumulative_absorption <- c(0, cumulative_absorption)

# Step 6: Baseline correction - shift the curve to start from zero
cumulative_absorption <- cumulative_absorption - cumulative_absorption[1]

# Step 7: Plot the reconstructed absorption spectrum
plot(roi_peak$Magnetic_Field, cumulative_absorption, type = 'l', col = 'blue',
     xlab = 'Magnetic Field (G)', ylab = 'Absorbance (a.u.)',
     main = 'Reconstructed Absorption Spectrum (Cumulative Trapezoidal Rule)')

# Step 8: Calculate the total area under the curve (AUC) using the trapezoidal rule
auc <- trapz(roi_peak$Magnetic_Field, cumulative_absorption)
print(paste("Area under the curve (AUC):", auc))


# Step 9: Highlight the area under the curve in yellow
polygon(c(roi_peak$Magnetic_Field, rev(roi_peak$Magnetic_Field)),
        c(rep(0, length(cumulative_absorption)), rev(cumulative_absorption)),
        col = rgb(1, 1, 0, 0.5), border = NA)  # Yellow area





