# Drosophila Flight Test Analysis

## Description

This R script was developed to **analyze the flight ability of Drosophila** using images from a **flight test**.  
Each image corresponds to an experimental condition and shows a sticky sheet where flies landed after being dropped.

The code performs:
- Image loading, rotation, and preprocessing (grayscale conversion, thresholding).
- Automatic detection of individual flies.
- Division of the image into **three vertical zones (top, middle, bottom)** and a **ground group** (flies that fell into ethanol).
- Calculation of fly distribution percentages per zone.
- Statistical comparison between conditions.
- Automatic export of all figures and results.

---

## Author

**Olivier JOSSE**

---

## Requirements

Make sure you have **R (≥ 4.0)** and the following packages installed:

```{r}
install.packages(c("imager", "ggplot2"))
```

---

## Input Data

Images should:
   - Be cropped to include only the sticky sheet (not the full tube).
   - Be placed in the file paths specified in the script (img1_path, img2_path, img3_path).
   - Be in .jpg format (modifiable).

Each image corresponds to one experimental condition:
   - Control Uninjured
   - Control Injured
   - RNAi Injured

---

## Main Steps of the Script

1) Image Loading
   - Paths for each image are defined manually.
   - Each image is loaded with imager::load.image().

2) Preprocessing
   - Conversion to grayscale.
   - Rotation by a user-defined angle.
   - Manual thresholding to binarize the image. Adjust threshold_valueX to match the expected number of flies.

3) Segmentation
   - Object detection (each fly is identified).
   - Extraction of centroid coordinates for each detected fly.

4) Spatial Analysis
   - Division of the image into 3 equal vertical zones.
   - Counting of flies in each zone.
   - Manual input of ground flies (ground_countX).

5) Calculations and Results
   - Computation of percentages per zone.
   - Console summary of results.
   - Export of:
      . Original and simplified plots per condition.
      . Comparative plots across conditions.

6) Statistical Analysis
   - Chi-square and Fisher’s exact tests for pairwise condition comparisons.
   - Automatic annotation of p-values on the barplots.

---

## Output

The script automatically generates several .png figures:
   - original_image_plot_CondX.png → raw image with zone overlays
   - simplified_plot_CondX.png → schematic representation
   - Image+Simplified_Plots_separate_or_together_CondX.png → combined visualization
   - Comparison_of_Drosophila_Distribution_in_Tubes.png → overall comparison plot

All figures are saved in high resolution (300 dpi).


---

## Parameters to adjust

imgX_path          -->   File path for each image 
angleX             -->   Rotation angle for alignment 
threshold_valueX   -->   Threshold value for fly detection 
ground_countX      -->   Manual count of flies on the ground


---

## Statistical Output 

Two main comparisons are performed:
   - Condition 1 vs Condition 2
   - Condition 2 vs Condition 3

Both Chi-square and Fisher’s exact tests are used.
Significance levels are automatically added to the plots (n.s., *, **, ***).


---

## Output Files 

File type	Description : 
    .png	                    -->   Figures (zones, comparisons, simplified)
    .RData / .csv (optional)	-->   Data summaries (if added later)


---

## How to Run

1. Open R or RStudio
2. Copy and paste the script
3. Update image paths and parameters
4. Run each code section


