# Read the EPR data, skipping the initial lines that don't contain data
#epr_data <- read.table("cuso.txt", header = FALSE, skip = 2)

# # Assign meaningful column names
# colnames(epr_data) <- c("Index", "Field.G.", "Intensity")
# 
# # Display the structure of the dataframe to verify the import
# str(epr_data)
# 
# # Plotting the EPR spectrum
# plot(epr_data$Field.G., epr_data$Intensity, type = "l",
#      xlab = "Magnetic Field (G)", ylab = "Signal Intensity",
#      main = "EPR Spectrum")
# 
# # Perform linear baseline correction
# baseline_model <- lm(Intensity ~ Field.G., data = epr_data)
# epr_data$CorrectedIntensity <- epr_data$Intensity - predict(baseline_model)
# 
# # Plot the baseline-corrected spectrum
# plot(epr_data$Field.G., epr_data$CorrectedIntensity, type = "l",
#      xlab = "Magnetic Field (G)", ylab = "Corrected Signal Intensity",
#      main = "Baseline-Corrected EPR Spectrum")
# 
# 
# #---- Calculate the Area Under the Curve (Integration)
# 
# # Load the pracma library for numerical integration
# library(pracma)
# 
# # Calculate the area under the curve
# area <- trapz(epr_data$Field.G., epr_data$CorrectedIntensity)
# 
# # Print the calculated area
# print(paste("Area under the EPR curve:", area))


# Load necessary libraries
library(ggplot2)
library(signal)

# Read the EPR data, skipping lines if necessary
epr_data <- read.table("Data/CUSO4_T2/CUSO4_D0.TXT", header = FALSE, skip = 2)

# Assign meaningful column names
colnames(epr_data) <- c("Index", "Magnetic_Field", "Signal_Intensity")

# Check the structure of the data
str(epr_data)

# Plotting the raw EPR spectrum
ggplot(epr_data, aes(x = Magnetic_Field, y = Signal_Intensity)) +
  geom_line() +
  labs(title = "EPR Spectrum", x = "Magnetic Field (Gauss)", y = "Signal Intensity")


# Baseline correction by subtracting the minimum signal intensity
epr_data$Corrected_Signal <- epr_data$Signal_Intensity - min(epr_data$Signal_Intensity)

# Plot the baseline-corrected spectrum
ggplot(epr_data, aes(x = Magnetic_Field, y = Corrected_Signal)) +
  geom_line() +
  labs(title = "Corrected EPR Spectrum", x = "Magnetic Field (Gauss)", y = "Corrected Signal Intensity")



library(pracma)

# Detect peaks in the corrected signal
# The 'minpeakheight' parameter is set to filter out smaller peaks

# findpeaks detects local maxima in the corrected signal.
# minpeakheight is used to filter out smaller peaks below a specified threshold (10% of the maximum signal intensity in this case).
# peak_positions and peak_intensities store the magnetic field values and corresponding intensities of detected peaks.

peaks <- findpeaks(epr_data$Corrected_Signal, minpeakheight = 0.1 * max(epr_data$Corrected_Signal))

# Extract the magnetic field values corresponding to the peaks
peak_positions <- epr_data$Magnetic_Field[peaks[, 2]]  # Column 2 gives the indices of peaks
peak_intensities <- peaks[, 1]  # Column 1 gives the peak values

# Plot the corrected spectrum with peak positions marked
ggplot(epr_data, aes(x = Magnetic_Field, y = Corrected_Signal)) +
  geom_line() +
  geom_point(data = data.frame(Magnetic_Field = peak_positions, Corrected_Signal = peak_intensities),
             aes(x = Magnetic_Field, y = Corrected_Signal), color = "red") +
  labs(title = "EPR Spectrum with Peak Positions", x = "Magnetic Field (Gauss)", y = "Corrected Signal Intensity")


# Calculate the g-Factor for Each Peak

# Step 1: Find the maximum and minimum points in the corrected signal
max_index <- which.max(epr_data$Corrected_Signal)
min_index <- which.min(epr_data$Corrected_Signal)

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

# Assuming 'epr_data$Corrected_Signal' contains the first derivative of the absorption spectrum
library(pracma)

# Integrate the first derivative to reconstruct the absorption spectrum
absorption_spectrum <- cumtrapz(epr_data$Magnetic_Field, epr_data$Corrected_Signal)

# Plot the reconstructed absorption spectrum
plot(epr_data$Magnetic_Field, absorption_spectrum, type = "l",
     xlab = "Magnetic Field (G)", ylab = "Absorption Signal",
     main = "Reconstructed Absorption Spectrum")


# Calculate the area under the reconstructed absorption spectrum
area_under_absorption <- trapz(epr_data$Magnetic_Field, absorption_spectrum)

# Output the area (proportional to radical concentration)
print(paste("Integrated Signal (Area under the absorption curve):", area_under_absorption))
