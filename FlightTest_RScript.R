######### Flight Test Analysis ##########

#author : Olivier JOSSE

#This code is designed to study the flight ability of drosophila during the flight test. 

# Before using this code, the images were cropped to keep only the sticky sheet. During analysis, they will undergo treatments such as conversion to grayscale and thresholding to facilitate segmentation of points representing drosophila. Each high-intensity point corresponds to a drosophila, with coordinates (x, y) on the image. It is therefore essential to manually adjust and test the threshold for each condition, and visually check that the detected points correspond to the drosophila before performing statistical analysis of their distribution.

#/!\ Things to ADAPT CAREFULLY depending on your study and what you do:
#   1) Path of the resized images
#   2) Rotation angle depending on the image
#   3) Threshold depending on the image according to the number of flies used (it should be roughly correct in total_count)
#   4) ground_count1 and ground_count2 which are the number of flies per condition that could not fly during the test and ended up stuck in ethanol on the ground.

######### Detection of drosophila in the flight test (3 conditions) ##########

# Load the imager package
install.packages("imager")
library(imager)

# Load the images
img1_path <- "E:/25.09.03 FlightTest Sparc Bnl RNAi Injured Females/Nouveaux fichiers R/250918_FlightTest_CTRLuninj_2_Cropped.jpg"
img1 <- load.image(img1_path)
img2_path <- "E:/25.09.03 FlightTest Sparc Bnl RNAi Injured Females/Nouveaux fichiers R/250918_FlightTest_CTRLinj_1_Cropped.jpg"
img2 <- load.image(img2_path)
img3_path <- "E:/25.09.03 FlightTest Sparc Bnl RNAi Injured Females/Nouveaux fichiers R/250918_FlightTest_BNLinj_1_Cropped.jpg"
img3 <- load.image(img3_path)

# Convert the image to grayscale
gray_img1 <- grayscale(img1)
gray_img2 <- grayscale(img2)
gray_img3 <- grayscale(img3)

# Display the image 
plot(gray_img1, main = "Grayscale Image")
plot(gray_img2, main = "Grayscale Image")
plot(gray_img3, main = "Grayscale Image")

# Define rotation angle
angle1 <- 0 # degrees
angle2 <- 0  # degrees
angle3 <- 0  # degrees

# Apply rotation
img1_rotated <- imrotate(img1, angle = angle1)
img2_rotated <- imrotate(img2, angle = angle2)
img3_rotated <- imrotate(img3, angle = angle3)

# Convert to grayscale (after rotation)
gray_img1_rotated <- grayscale(img1_rotated)
gray_img2_rotated <- grayscale(img2_rotated)
gray_img3_rotated <- grayscale(img3_rotated)

# Display rotated grayscale image
plot(gray_img1_rotated, main = paste("Rotated Image by", angle1, "degrees"))
plot(gray_img2_rotated, main = paste("Rotated Image by", angle2, "degrees"))
plot(gray_img3_rotated, main = paste("Rotated Image by", angle3, "degrees"))

# Apply threshold to binarize the image
par(mfrow = c(1,5))
plot(gray_img1_rotated > 0.3, main = "Threshold 0.3")
plot(gray_img1_rotated > 0.35, main = "Threshold 0.35")
plot(gray_img1_rotated > 0.4, main = "Threshold 0.4")
plot(gray_img1_rotated > 0.45, main = "Threshold 0.45")
plot(gray_img1_rotated > 0.5, main = "Threshold 0.5")
threshold_value1 <- 0.2  #ATTENTION Adjust this value according to your needs - total_count after this must match the number of flies launched in the test. Otherwise, adjust the threshold.
binary_img1 <- gray_img1_rotated > threshold_value1

par(mfrow = c(1,5))
plot(gray_img2_rotated > 0.3, main = "Threshold 0.3")
plot(gray_img2_rotated > 0.35, main = "Threshold 0.35")
plot(gray_img2_rotated > 0.4, main = "Threshold 0.4")
plot(gray_img2_rotated > 0.45, main = "Threshold 0.45")
plot(gray_img2_rotated > 0.5, main = "Threshold 0.5")
threshold_value2 <- 0.2  #ATTENTION Adjust this value according to your needs
binary_img2 <- gray_img2_rotated > threshold_value2

par(mfrow = c(1,5))
plot(gray_img3_rotated > 0.3, main = "Threshold 0.3")
plot(gray_img3_rotated > 0.35, main = "Threshold 0.35")
plot(gray_img3_rotated > 0.4, main = "Threshold 0.4")
plot(gray_img3_rotated > 0.45, main = "Threshold 0.45")
plot(gray_img3_rotated > 0.5, main = "Threshold 0.5")
threshold_value3 <- 0.25  #ATTENTION Adjust this value according to your needs
binary_img3 <- gray_img3_rotated > threshold_value3

# Identify objects (particles)
labeled_img1 <- label(binary_img1)
labeled_img2 <- label(binary_img2)
labeled_img3 <- label(binary_img3)

# Get centroid coordinates for each object
centroids1 <- list()
for (i in unique(labeled_img1)[-1]) {
  obj_pixels <- which(labeled_img1 == i, arr.ind = TRUE)
  centroid <- colMeans(obj_pixels)
  centroids1[[i]] <- centroid
}

centroids2 <- list()
for (i in unique(labeled_img2)[-1]) {
  obj_pixels <- which(labeled_img2 == i, arr.ind = TRUE)
  centroid <- colMeans(obj_pixels)
  centroids2[[i]] <- centroid
}

centroids3 <- list()
for (i in unique(labeled_img3)[-1]) {
  obj_pixels <- which(labeled_img3 == i, arr.ind = TRUE)
  centroid <- colMeans(obj_pixels)
  centroids3[[i]] <- centroid
}

# Convert centroids to data frame
centroids_df1 <- as.data.frame(do.call(rbind, centroids1))
colnames(centroids_df1) <- c("xcentroid", "ycentroid")

centroids_df2 <- as.data.frame(do.call(rbind, centroids2))
colnames(centroids_df2) <- c("xcentroid", "ycentroid")

centroids_df3 <- as.data.frame(do.call(rbind, centroids3))
colnames(centroids_df3) <- c("xcentroid", "ycentroid")

# Divide the image into three parts (top, middle, bottom)
img1_height <- dim(gray_img1_rotated)[2]  # height (Y axis)
img1_width  <- dim(gray_img1_rotated)[1]  # width (X axis)
third_height1 <- img1_height / 3

img2_height <- dim(gray_img2_rotated)[2]
img2_width  <- dim(gray_img2_rotated)[1]
third_height2 <- img2_height / 3

img3_height <- dim(gray_img3_rotated)[2]
img3_width  <- dim(gray_img3_rotated)[1]
third_height3 <- img3_height / 3

# Define limits for each part
top_part1 <- centroids_df1[centroids_df1$ycentroid < third_height1, ]
middle_part1 <- centroids_df1[centroids_df1$ycentroid >= third_height1 & centroids_df1$ycentroid < 2 * third_height1, ]
bottom_part1 <- centroids_df1[centroids_df1$ycentroid >= 2 * third_height1, ]

top_part2 <- centroids_df2[centroids_df2$ycentroid < third_height2, ]
middle_part2 <- centroids_df2[centroids_df2$ycentroid >= third_height2 & centroids_df2$ycentroid < 2 * third_height2, ]
bottom_part2 <- centroids_df2[centroids_df2$ycentroid >= 2 * third_height2, ]

top_part3 <- centroids_df3[centroids_df3$ycentroid < third_height3, ]
middle_part3 <- centroids_df3[centroids_df3$ycentroid >= third_height3 & centroids_df3$ycentroid < 2 * third_height3, ]
bottom_part3 <- centroids_df3[centroids_df3$ycentroid >= 2 * third_height3, ]

# Count points in each part
top_count1 <- nrow(top_part1)
middle_count1 <- nrow(middle_part1)
bottom_count1 <- nrow(bottom_part1)
ground_count1 <- 0 ##SET THE NUMBER OF FLIES IN ETHANOL BATH

top_count2 <- nrow(top_part2)
middle_count2 <- nrow(middle_part2)
bottom_count2 <- nrow(bottom_part2)
ground_count2 <- 12 ##SET THE NUMBER OF FLIES IN ETHANOL BATH

top_count3 <- nrow(top_part3)
middle_count3 <- nrow(middle_part3)
bottom_count3 <- nrow(bottom_part3)
ground_count3 <- 13 ##SET THE NUMBER OF FLIES IN ETHANOL BATH

# Calculate percentages for each part
total_count1 <- nrow(centroids_df1) + ground_count1
top_percentage1 <- (top_count1 / total_count1) * 100
middle_percentage1 <- (middle_count1 / total_count1) * 100
bottom_percentage1 <- (bottom_count1 / total_count1) * 100
ground_percentage1 <- (ground_count1 / total_count1) * 100

total_count2 <- nrow(centroids_df2) + ground_count2
top_percentage2 <- (top_count2 / total_count2) * 100
middle_percentage2 <- (middle_count2 / total_count2) * 100
bottom_percentage2 <- (bottom_count2 / total_count2) * 100
ground_percentage2 <- (ground_count2 / total_count2) * 100

total_count3 <- nrow(centroids_df3) + ground_count3
top_percentage3 <- (top_count3 / total_count3) * 100
middle_percentage3 <- (middle_count3 / total_count3) * 100
bottom_percentage3 <- (bottom_count3 / total_count3) * 100
ground_percentage3 <- (ground_count3 / total_count3) * 100

# Display results
total_count1
cat("Number of points in top part:", top_count1, "\n")
cat("Number of points in middle part:", middle_count1, "\n")
cat("Number of points in bottom part:", bottom_count1, "\n")
cat("Number of points in ground part:", ground_count1, "\n")
cat("Percentage in top part:", top_percentage1, "%\n")
cat("Percentage in middle part:", middle_percentage1, "%\n")
cat("Percentage in bottom part:", bottom_percentage1, "%\n")
cat("Percentage in ground part:", ground_percentage1, "%\n")

total_count2
cat("Number of points in top part:", top_count2, "\n")
cat("Number of points in middle part:", middle_count2, "\n")
cat("Number of points in bottom part:", bottom_count2, "\n")
cat("Number of points in ground part:", ground_count2, "\n")
cat("Percentage in top part:", top_percentage2, "%\n")
cat("Percentage in middle part:", middle_percentage2, "%\n")
cat("Percentage in bottom part:", bottom_percentage2, "%\n")
cat("Percentage in ground part:", ground_percentage2, "%\n")

total_count3
cat("Number of points in top part:", top_count3, "\n")
cat("Number of points in middle part:", middle_count3, "\n")
cat("Number of points in bottom part:", bottom_count3, "\n")
cat("Number of points in ground part:", ground_count3, "\n")
cat("Percentage in top part:", top_percentage3, "%\n")
cat("Percentage in middle part:", middle_percentage3, "%\n")
cat("Percentage in bottom part:", bottom_percentage3, "%\n")
cat("Percentage in ground part:", ground_percentage3, "%\n")

# Plot the zones on the image
plot(gray_img1_rotated, main = "Grayscale Image with Zones")
abline(h = third_height1, col = "green", lwd = 2)
abline(h = 2 * third_height1, col = "green", lwd = 2)
points(centroids_df1, col = "blue", pch = 19, cex = 0.2)

plot(gray_img2_rotated, main = "Grayscale Image with Zones")
abline(h = third_height2, col = "green", lwd = 2)
abline(h = 2 * third_height2, col = "green", lwd = 2)
points(centroids_df2, col = "blue", pch = 19, cex = 0.2)

plot(gray_img3_rotated, main = "Grayscale Image with Zones")
abline(h = third_height3, col = "green", lwd = 2)
abline(h = 2 * third_height3, col = "green", lwd = 2)
points(centroids_df3, col = "blue", pch = 19, cex = 0.2)

######### Create a simplified plot with original image or without for both conditions #########

par(mfrow = c(1,1)) #Use c(1,2) for two columns of images or c(2,1) for two rows

### Condition 1

# Plot with background image and simplified zones
plot(NA, xlim = c(0, img1_width), ylim = c(img1_height, 0), 
     xlab = "", ylab = "", main = "Image + Simplified Zones Condition 1", asp = 1) # Create an empty plot with correct dimensions and origin at top-left
rasterImage(as.raster(gray_img1_rotated), 
            xleft = 0, ybottom = img1_height, 
            xright = img1_width, ytop = 0) # Display the rotated image as background
abline(h = third_height1, col = "green", lwd = 2) # Draw separation lines
abline(h = 2 * third_height1, col = "green", lwd = 2)
points(centroids_df1$xcentroid, centroids_df1$ycentroid, col = "blue", pch = 19, cex = 0.2) # Plot points

# Simplified plot without background image
plot(NA, xlim = c(0, img1_width), ylim = c(img1_height, 0), 
     xlab = "", ylab = "", main = "Simplified Image Condition 1", asp = 1)
abline(h = third_height1, col = "green", lwd = 2)
abline(h = 2 * third_height1, col = "green", lwd = 2)
points(centroids_df1$xcentroid, centroids_df1$ycentroid, col = "blue", pch = 19, cex = 0.2)


## Save both plots (simplified plot + image plot) in high resolution (300 dpi or more)

# Set output file
png("original_image_plot_Cond1.png", width = 2000, height = 2000, res = 300)
# Repeat the plot code here
plot(NA, xlim = c(0, img1_width), ylim = c(img1_height, 0), 
     xlab = "", ylab = "", main = "Image Condition 1", asp = 1)
rasterImage(as.raster(gray_img1_rotated), 
            xleft = 0, ybottom = img1_height, 
            xright = img1_width, ytop = 0)
abline(h = third_height1, col = "green", lwd = 2)
abline(h = 2 * third_height1, col = "green", lwd = 2)
dev.off()

# Set output file
png("simplified_with_original_image_Cond1.png", width = 2000, height = 2000, res = 300)
plot(NA, xlim = c(0, img1_width), ylim = c(img1_height, 0), 
     xlab = "", ylab = "", main = "Image + Simplified Points Condition 1", asp = 1)
rasterImage(as.raster(gray_img1_rotated), 
            xleft = 0, ybottom = img1_height, 
            xright = img1_width, ytop = 0)
abline(h = third_height1, col = "green", lwd = 2)
abline(h = 2 * third_height1, col = "green", lwd = 2)
points(centroids_df1$xcentroid, centroids_df1$ycentroid, col = "blue", pch = 19, cex = 0.2)
dev.off()

# Set output file
png("simplified_plot_Cond1.png", width = 2000, height = 2000, res = 300)
plot(NA, xlim = c(0, img1_width), ylim = c(img1_height, 0), 
     xlab = "", ylab = "", main = "Simplified Image Condition 1", asp = 1)
abline(h = third_height1, col = "green", lwd = 2)
abline(h = 2 * third_height1, col = "green", lwd = 2)
points(centroids_df1$xcentroid, centroids_df1$ycentroid, col = "blue", pch = 19, cex = 0.2)
dev.off()

# Set output file
png("Image+Simplified_Plots_separate_or_together_Cond1.png", width = 2000, height = 2000, res = 300)
par(mfrow = c(1,3)) # To show all 3 plots together
plot(NA, xlim = c(0, img1_width), ylim = c(img1_height, 0), 
     xlab = "", ylab = "", main = "Image Condition 1", asp = 1)
rasterImage(as.raster(gray_img1_rotated), 
            xleft = 0, ybottom = img1_height, 
            xright = img1_width, ytop = 0)
abline(h = third_height1, col = "green", lwd = 2)
abline(h = 2 * third_height1, col = "green", lwd = 2)
plot(NA, xlim = c(0, img1_width), ylim = c(img1_height, 0), 
     xlab = "", ylab = "", main = "Image + Simplified Points Condition 1", asp = 1)
rasterImage(as.raster(gray_img1_rotated), 
            xleft = 0, ybottom = img1_height, 
            xright = img1_width, ytop = 0)
abline(h = third_height1, col = "green", lwd = 2)
abline(h = 2 * third_height1, col = "green", lwd = 2)
points(centroids_df1$xcentroid, centroids_df1$ycentroid, col = "blue", pch = 19, cex = 0.2)
plot(NA, xlim = c(0, img1_width), ylim = c(img1_height, 0), 
     xlab = "", ylab = "", main = "Simplified Image Condition 1", asp = 1)
abline(h = third_height1, col = "green", lwd = 2)
abline(h = 2 * third_height1, col = "green", lwd = 2)
points(centroids_df1$xcentroid, centroids_df1$ycentroid, col = "blue", pch = 19, cex = 0.2)
dev.off()

### Condition 2

# Plot with background image and simplified zones
plot(NA, xlim = c(0, img2_width), ylim = c(img2_height, 0), 
     xlab = "", ylab = "", main = "Image + Simplified Zones Condition 2", asp = 1)
rasterImage(as.raster(gray_img2_rotated), 
            xleft = 0, ybottom = img2_height, 
            xright = img2_width, ytop = 0)
abline(h = third_height2, col = "green", lwd = 2)
abline(h = 2 * third_height2, col = "green", lwd = 2)
points(centroids_df2$xcentroid, centroids_df2$ycentroid, col = "blue", pch = 19, cex = 0.2)

# Simplified plot without background image
plot(NA, xlim = c(0, img2_width), ylim = c(img2_height, 0), 
     xlab = "", ylab = "", main = "Simplified Image Condition 2", asp = 1)
abline(h = third_height2, col = "green", lwd = 2)
abline(h = 2 * third_height2, col = "green", lwd = 2)
points(centroids_df2$xcentroid, centroids_df2$ycentroid, col = "blue", pch = 19, cex = 0.2)

## Save both plots in high resolution

png("original_image_plot_Cond2.png", width = 2000, height = 2000, res = 300)
plot(NA, xlim = c(0, img2_width), ylim = c(img2_height, 0), 
     xlab = "", ylab = "", main = "Image Condition 2", asp = 1)
rasterImage(as.raster(gray_img2_rotated), 
            xleft = 0, ybottom = img2_height, 
            xright = img2_width, ytop = 0)
abline(h = third_height2, col = "green", lwd = 2)
abline(h = 2 * third_height2, col = "green", lwd = 2)
dev.off()

png("simplified_with_original_image_Cond2.png", width = 2000, height = 2000, res = 300)
plot(NA, xlim = c(0, img2_width), ylim = c(img2_height, 0), 
     xlab = "", ylab = "", main = "Image + Simplified Points Condition 2", asp = 1)
rasterImage(as.raster(gray_img2_rotated), 
            xleft = 0, ybottom = img2_height, 
            xright = img2_width, ytop = 0)
abline(h = third_height2, col = "green", lwd = 2)
abline(h = 2 * third_height2, col = "green", lwd = 2)
points(centroids_df2$xcentroid, centroids_df2$ycentroid, col = "blue", pch = 19, cex = 0.2)
dev.off()

png("simplified_plot_Cond2.png", width = 2000, height = 2000, res = 300)
plot(NA, xlim = c(0, img2_width), ylim = c(img2_height, 0), 
     xlab = "", ylab = "", main = "Simplified Image Condition 2", asp = 1)
abline(h = third_height2, col = "green", lwd = 2)
abline(h = 2 * third_height2, col = "green", lwd = 2)
points(centroids_df2$xcentroid, centroids_df2$ycentroid, col = "blue", pch = 19, cex = 0.2)
dev.off()

png("Image+Simplified_Plots_separate_or_together_Cond2.png", width = 2000, height = 2000, res = 300)
par(mfrow = c(1,3))
plot(NA, xlim = c(0, img2_width), ylim = c(img2_height, 0), 
     xlab = "", ylab = "", main = "Image Condition 2", asp = 1)
rasterImage(as.raster(gray_img2_rotated), 
            xleft = 0, ybottom = img2_height, 
            xright = img2_width, ytop = 0)
abline(h = third_height2, col = "green", lwd = 2)
abline(h = 2 * third_height2, col = "green", lwd = 2)
plot(NA, xlim = c(0, img2_width), ylim = c(img2_height, 0), 
     xlab = "", ylab = "", main = "Image + Simplified Points Condition 2", asp = 1)
rasterImage(as.raster(gray_img2_rotated), 
            xleft = 0, ybottom = img2_height, 
            xright = img2_width, ytop = 0)
abline(h = third_height2, col = "green", lwd = 2)
abline(h = 2 * third_height2, col = "green", lwd = 2)
points(centroids_df2$xcentroid, centroids_df2$ycentroid, col = "blue", pch = 19, cex = 0.2)
plot(NA, xlim = c(0, img2_width), ylim = c(img2_height, 0), 
     xlab = "", ylab = "", main = "Simplified Image Condition 2", asp = 1)
abline(h = third_height2, col = "green", lwd = 2)
abline(h = 2 * third_height2, col = "green", lwd = 2)
points(centroids_df2$xcentroid, centroids_df2$ycentroid, col = "blue", pch = 19, cex = 0.2)
dev.off()

### Condition 3

# Plot with background image and simplified zones
plot(NA, xlim = c(0, img3_width), ylim = c(img3_height, 0), 
     xlab = "", ylab = "", main = "Image + Simplified Zones Condition 3", asp = 1)
rasterImage(as.raster(gray_img3_rotated), 
            xleft = 0, ybottom = img3_height, 
            xright = img3_width, ytop = 0)
abline(h = third_height3, col = "green", lwd = 2)
abline(h = 2 * third_height3, col = "green", lwd = 2)
points(centroids_df3$xcentroid, centroids_df3$ycentroid, col = "blue", pch = 19, cex = 0.2)

# Simplified plot without background image
plot(NA, xlim = c(0, img3_width), ylim = c(img3_height, 0), 
     xlab = "", ylab = "", main = "Simplified Image Condition 3", asp = 1)
abline(h = third_height3, col = "green", lwd = 2)
abline(h = 2 * third_height3, col = "green", lwd = 2)
points(centroids_df3$xcentroid, centroids_df3$ycentroid, col = "blue", pch = 19, cex = 0.2)

## Save both plots in high resolution

png("original_image_plot_Cond3.png", width = 2000, height = 2000, res = 300)
plot(NA, xlim = c(0, img3_width), ylim = c(img3_height, 0), 
     xlab = "", ylab = "", main = "Image Condition 3", asp = 1)
rasterImage(as.raster(gray_img3_rotated), 
            xleft = 0, ybottom = img3_height, 
            xright = img3_width, ytop = 0)
abline(h = third_height3, col = "green", lwd = 2)
abline(h = 2 * third_height3, col = "green", lwd = 2)
dev.off()

png("simplified_with_original_image_Cond3.png", width = 2000, height = 2000, res = 300)
plot(NA, xlim = c(0, img3_width), ylim = c(img3_height, 0), 
     xlab = "", ylab = "", main = "Image + Simplified Points Condition 3", asp = 1)
rasterImage(as.raster(gray_img3_rotated), 
            xleft = 0, ybottom = img3_height, 
            xright = img3_width, ytop = 0)
abline(h = third_height3, col = "green", lwd = 2)
abline(h = 2 * third_height3, col = "green", lwd = 2)
points(centroids_df3$xcentroid, centroids_df3$ycentroid, col = "blue", pch = 19, cex = 0.2)
dev.off()

png("simplified_plot_Cond3.png", width = 2000, height = 2000, res = 300)
plot(NA, xlim = c(0, img3_width), ylim = c(img3_height, 0), 
     xlab = "", ylab = "", main = "Simplified Image Condition 3", asp = 1)
abline(h = third_height3, col = "green", lwd = 2)
abline(h = 2 * third_height3, col = "green", lwd = 2)
points(centroids_df3$xcentroid, centroids_df3$ycentroid, col = "blue", pch = 19, cex = 0.2)
dev.off()

png("Image+Simplified_Plots_separate_or_together_Cond3.png", width = 2000, height = 2000, res = 300)
par(mfrow = c(1,3))
plot(NA, xlim = c(0, img3_width), ylim = c(img3_height, 0), 
     xlab = "", ylab = "", main = "Image Condition 3", asp = 1)
rasterImage(as.raster(gray_img3_rotated), 
            xleft = 0, ybottom = img3_height, 
            xright = img3_width, ytop = 0)
abline(h = third_height3, col = "green", lwd = 2)
abline(h = 2 * third_height3, col = "green", lwd = 2)
plot(NA, xlim = c(0, img3_width), ylim = c(img3_height, 0), 
     xlab = "", ylab = "", main = "Image + Simplified Points Condition 3", asp = 1)
rasterImage(as.raster(gray_img3_rotated), 
            xleft = 0, ybottom = img3_height, 
            xright = img3_width, ytop = 0)
abline(h = third_height3, col = "green", lwd = 2)
abline(h = 2 * third_height3, col = "green", lwd = 2)
points(centroids_df3$xcentroid, centroids_df3$ycentroid, col = "blue", pch = 19, cex = 0.2)
plot(NA, xlim = c(0, img3_width), ylim = c(img3_height, 0), 
     xlab = "", ylab = "", main = "Simplified Image Condition 3", asp = 1)
abline(h = third_height3, col = "green", lwd = 2)
abline(h = 2 * third_height3, col = "green", lwd = 2)
points(centroids_df3$xcentroid, centroids_df3$ycentroid, col = "blue", pch = 19, cex = 0.2)
dev.off()

######### Create a stacked % bar plot and comparison between the three conditions #########

# Create a dataset for the comparative plot
graph_data <- data.frame(
  Section = rep(c("Haut", "Milieu", "Bas", "Sol"), 3),
  Pourcentage = c(top_percentage1, middle_percentage1, bottom_percentage1, ground_percentage1,
                  top_percentage2, middle_percentage2, bottom_percentage2, ground_percentage2,
                  top_percentage3, middle_percentage3, bottom_percentage3, ground_percentage3),
  Condition = rep(c("Control Uninjured", "Control Injured", "BnlRNAi Injured"), each = 4)
)

# Force the order of conditions in the plot (Uninjured on the left, Injured on the right)
graph_data$Condition <- factor(graph_data$Condition,
                               levels = c("Control Uninjured", "Control Injured", "BnlRNAi Injured"))

# Reorder the sections to have "Haut" at the top, "Milieu" in the middle, "Bas" at the bottom
graph_data$Section <- factor(graph_data$Section,
                             levels = c("Haut", "Milieu", "Bas", "Sol"),
                             labels = c("Top", "Middle", "Bottom", "Ground"))

# Use ggplot2 to create a 100% stacked bar plot
library(ggplot2)

ggplot(graph_data, aes(x = Condition, y = Pourcentage, fill = Section)) +
  geom_bar(stat = "identity", position = "fill") +  # position="fill" for 100% total
  labs(title = "Drosophila Distribution by Zone and Condition", x = "Condition", y = "Percentage", fill = "Section") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +  # Display percentages
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.ticks.x = element_blank())

### Statistical tests (Chi2) to see differences between percentages ###

## Cond1 vs Cond2

# Build a contingency table with counts
contingency_table12 <- matrix(
  c(top_count1, middle_count1, bottom_count1, ground_count1,
    top_count2, middle_count2, bottom_count2, ground_count2),
  nrow = 2,
  byrow = TRUE
)
colnames(contingency_table12) <- c("Haut", "Milieu", "Bas", "Sol")
rownames(contingency_table12) <- c("Control Uninjured", "Control Injured")

# Display the table
print("Contingency Table:")
print(contingency_table12)

# Chi-squared test of independence
chi2_result12 <- chisq.test(contingency_table12)

# Test result
print(chi2_result12)

# If any count < 5, Fisher's test is more appropriate:
if(any(contingency_table12 < 5)) {
  fisher_result12 <- fisher.test(contingency_table12)
  print("Fisher's Test (low counts):")
  print(fisher_result12)
}

## Cond2 vs Cond3

# Build a contingency table with counts
contingency_table23 <- matrix(
  c(top_count2, middle_count2, bottom_count2, ground_count2,
    top_count3, middle_count3, bottom_count3, ground_count3),
  nrow = 2, byrow = TRUE
)
colnames(contingency_table23) <- c("Haut", "Milieu", "Bas", "Sol")
rownames(contingency_table23) <- c("Control Injured", "BnlRNAi Injured")

# Display the table
print("Contingency Table:")
print(contingency_table23)

# Chi-squared test of independence
chi2_result23 <- chisq.test(contingency_table23)

# Test result
print(chi2_result23)

# If any count < 5, Fisher's test is more appropriate:
if(any(contingency_table23 < 5)) {
  fisher_result23 <- fisher.test(contingency_table23)
  print("Fisher's Test (low counts):")
  print(fisher_result23)
}

### Add statistical results to ggplot ###

# Chi2 test results
p_val12 <- chi2_result12$p.value
p_val23 <- chi2_result23$p.value

# Convert p-value to annotation
label_p <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "n.s."
}

# Create the graph with comparison lines for Cond1/2 and Cond2/3
annot_data <- data.frame(
  x = c(1.5, 2.5),                           # positions of significance between conditions
  y = c(1.02, 1.02),                         # height of the annotations
  label = c(label_p(p_val12), label_p(p_val23)))
annot_data$vjust <- ifelse(annot_data$label == "ns", -0.4, 0) # Add a specific offset column: ns slightly higher # ↑ only moves the "ns"
ggplot(graph_data, aes(x = Condition, y = Pourcentage, fill = Section)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Répartition des Drosophiles par Zone et Condition", y = "Pourcentage", fill = "Section") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.ticks.x = element_blank()) +
  # Annotation cond1 vs cond2 Horizontal and vertical lines
  geom_segment(aes(x = 1, xend = 2, y = 1.02, yend = 1.02)) +
  geom_segment(aes(x = 1, xend = 1, y = 1.00, yend = 1.02)) +
  geom_segment(aes(x = 2, xend = 2, y = 1.00, yend = 1.02)) +
  # Annotation cond2 vs cond3
  geom_segment(aes(x = 2, xend = 3, y = 1.02, yend = 1.02)) +
  geom_segment(aes(x = 2, xend = 2, y = 1.00, yend = 1.02)) +
  geom_segment(aes(x = 3, xend = 3, y = 1.00, yend = 1.02)) +
  geom_text(
    data = annot_data,
    aes(x = x, y = y, label = label, vjust = vjust),
    inherit.aes = FALSE,   
    size = 18)             

### Save the two graphs (simplified graph + Image and simplified graph) in high resolution (300 or more) ###

# Define output file
png("Comparaison de la répartition des drosophiles dans les tubes.png", width = 5000, height = 7000, res = 600)
# Data for significance annotations and positions
annot_data <- data.frame(
  x = c(1.5, 2.5),                           # positions of significance between conditions
  y = c(1.02, 1.02),                         # height of the annotations
  label = c(label_p(p_val12), label_p(p_val23)))
annot_data$vjust <- ifelse(annot_data$label == "ns", -0.4, 0) # Add a specific offset column: ns slightly higher # ↑ only moves the "ns"
# Reproduce your graph code here
ggplot(graph_data, aes(x = Condition, y = Pourcentage, fill = Section)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Répartition des Drosophiles par Zone et Condition", y = "Distribution of flies across tube sections", fill = "Section") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 18), # Legend angle and size
    axis.text.y  = element_text(size = 14), # enlarged Y-axis numbers
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),  # shifted Y-axis title
    axis.ticks.x = element_blank(),
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.key.size = unit(1.2, "cm"),   # enlarge legend squares
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  # Annotation cond1 vs cond2
  geom_segment(aes(x = 1, xend = 2, y = 1.02, yend = 1.02)) +
  geom_segment(aes(x = 1, xend = 1, y = 1.00, yend = 1.02)) +
  geom_segment(aes(x = 2, xend = 2, y = 1.00, yend = 1.02)) +
  # Annotation cond2 vs cond3
  geom_segment(aes(x = 2, xend = 3, y = 1.02, yend = 1.02)) +
  geom_segment(aes(x = 2, xend = 2, y = 1.00, yend = 1.02)) +
  geom_segment(aes(x = 3, xend = 3, y = 1.00, yend = 1.02)) +
  geom_text(
    data = annot_data,
    aes(x = x, y = y, label = label, vjust = vjust),
    inherit.aes = FALSE,   # very important!
    size = 18)             # uniform size between stars and ns        
# Close the graphics device to save the file
dev.off()

######### Create a stacked % bar plot and comparison between the three conditions with Bottom and Ground combined #########

# Create a dataset for the comparative plot
graph_data <- data.frame(
  Section = rep(c("Haut", "Milieu", "Bas", "Sol"), 3),
  Pourcentage = c(top_percentage1, middle_percentage1, bottomground_percentage1,
                  top_percentage2, middle_percentage2, bottomground_percentage2,
                  top_percentage3, middle_percentage3, bottomground_percentage3),
  Condition = rep(c("Control Uninjured", "Control Injured", "BnlRNAi Injured"), each = 4)
)

# Force the order of conditions in the plot (Uninjured on the left, Injured on the right)
graph_data$Condition <- factor(graph_data$Condition,
                               levels = c("Control Uninjured", "Control Injured", "BnlRNAi Injured"))

# Reorder the sections to have "Haut" at the top, "Milieu" in the middle, "Bas" at the bottom
graph_data$Section <- factor(graph_data$Section,
                             levels = c("Haut", "Milieu", "Bas", "Sol"),
                             labels = c("Top", "Middle", "Bottom", "Ground"))

# Use ggplot2 to create a 100% stacked bar plot
library(ggplot2)

ggplot(graph_data, aes(x = Condition, y = Pourcentage, fill = Section)) +
  geom_bar(stat = "identity", position = "fill") +  # position="fill" for 100% total
  labs(title = "Drosophila Distribution by Zone and Condition", x = "Condition", y = "Percentage", fill = "Section") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +  # Display percentages
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.ticks.x = element_blank())

### Statistical tests (Chi2) to see differences between percentages ###

## Cond1 vs Cond2

# Build a contingency table with counts
contingency_table12 <- matrix(
  c(top_count1, middle_count1, bottomground_count1,
    top_count2, middle_count2, bottomground_count2),
  nrow = 2,
  byrow = TRUE
)
colnames(contingency_table12) <- c("Haut", "Milieu", "Bas", "Sol")
rownames(contingency_table12) <- c("Control Uninjured", "Control Injured")

# Display the table
print("Contingency Table:")
print(contingency_table12)

# Chi-squared test of independence
chi2_result12 <- chisq.test(contingency_table12)

# Test result
print(chi2_result12)

# If any count < 5, Fisher's test is more appropriate:
if(any(contingency_table12 < 5)) {
  fisher_result12 <- fisher.test(contingency_table12)
  print("Fisher's Test (low counts):")
  print(fisher_result12)
}

## Cond2 vs Cond3

# Build a contingency table with counts
contingency_table23 <- matrix(
  c(top_count2, middle_count2, bottomground_count2,
    top_count3, middle_count3, bottomground_count3),
  nrow = 2, byrow = TRUE
)
colnames(contingency_table23) <- c("Haut", "Milieu", "Bas", "Sol")
rownames(contingency_table23) <- c("Control Injured", "BnlRNAi Injured")

# Display the table
print("Contingency Table:")
print(contingency_table23)

# Chi-squared test of independence
chi2_result23 <- chisq.test(contingency_table23)

# Test result
print(chi2_result23)

# If any count < 5, Fisher's test is more appropriate:
if(any(contingency_table23 < 5)) {
  fisher_result23 <- fisher.test(contingency_table23)
  print("Fisher's Test (low counts):")
  print(fisher_result23)
}

### Add statistical results to ggplot ###

# Chi2 test results
p_val12 <- chi2_result12$p.value
p_val23 <- chi2_result23$p.value

# Convert p-value to annotation
label_p <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "n.s."
}

# Create the graph with comparison lines for Cond1/2 and Cond2/3
annot_data <- data.frame(
  x = c(1.5, 2.5),                           # positions of significance between conditions
  y = c(1.02, 1.02),                         # height of the annotations
  label = c(label_p(p_val12), label_p(p_val23)))
annot_data$vjust <- ifelse(annot_data$label == "ns", -0.4, 0) # Add a specific offset column: ns slightly higher # ↑ only moves the "ns"
ggplot(graph_data, aes(x = Condition, y = Pourcentage, fill = Section)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Répartition des Drosophiles par Zone et Condition", y = "Pourcentage", fill = "Section") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.ticks.x = element_blank()) +
  # Annotation cond1 vs cond2 Horizontal and vertical lines
  geom_segment(aes(x = 1, xend = 2, y = 1.02, yend = 1.02)) +
  geom_segment(aes(x = 1, xend = 1, y = 1.00, yend = 1.02)) +
  geom_segment(aes(x = 2, xend = 2, y = 1.00, yend = 1.02)) +
  annotate("text", x = 1.5, y = 1.05, label = label_p(p_val12), size = 14) +
  # Annotation cond2 vs cond3
  geom_segment(aes(x = 2, xend = 3, y = 1.02, yend = 1.02)) +
  geom_segment(aes(x = 2, xend = 2, y = 1.00, yend = 1.02)) +
  geom_segment(aes(x = 3, xend = 3, y = 1.00, yend = 1.02)) +
  annotate("text", x = 2.5, y = 1.05, label = label_p(p_val23), size = 14) # Significance annotation

### Save the two graphs (simplified graph + Image and simplified graph) in high resolution (300 or more) ###

# Define output file
png("Comparaison de la répartition des drosophiles dans les tubes.png", width = 5000, height = 7000, res = 600)
# Data for significance annotations and positions
annot_data <- data.frame(
  x = c(1.5, 2.5),                           # positions of significance between conditions
  y = c(1.02, 1.02),                         # height of the annotations
  label = c(label_p(p_val12), label_p(p_val23)))
annot_data$vjust <- ifelse(annot_data$label == "ns", -0.4, 0) # Add a specific offset column: ns slightly higher # ↑ only moves the "ns"
# Reproduce your graph code here
ggplot(graph_data, aes(x = Condition, y = Pourcentage, fill = Section)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Répartition des Drosophiles par Zone et Condition", y = "Distribution of flies across tube sections", fill = "Section") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 18), # Legend angle and size
    axis.text.y  = element_text(size = 14), # enlarged Y-axis numbers
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),  # shifted Y-axis title
    axis.ticks.x = element_blank(),
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.key.size = unit(1.2, "cm"),   # enlarge legend squares
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  # Annotation cond1 vs cond2
  geom_segment(aes(x = 1, xend = 2, y = 1.02, yend = 1.02)) +
  geom_segment(aes(x = 1, xend = 1, y = 1.00, yend = 1.02)) +
  geom_segment(aes(x = 2, xend = 2, y = 1.00, yend = 1.02)) +
  # Annotation cond2 vs cond3
  geom_segment(aes(x = 2, xend = 3, y = 1.02, yend = 1.02)) +
  geom_segment(aes(x = 2, xend = 2, y = 1.00, yend = 1.02)) +
  geom_segment(aes(x = 3, xend = 3, y = 1.00, yend = 1.02)) +
  geom_text(
    data = annot_data,
    aes(x = x, y = y, label = label, vjust = vjust),
    inherit.aes = FALSE,   # very important!
    size = 18)             # uniform size between stars and ns        
# Close the graphics device to save the file
dev.off()


