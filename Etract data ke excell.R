library(terra)
library(sf)
library(readxl)
library(RStoolbox)
library(dplyr)
library(ggplot2)
library(tidyr)
library(writexl)

# ==============================================================================
# 1. LOAD DATA LAPANGAN & CITRA UTAMA (LANDSAT)
# ==============================================================================

# Load data polygon plot lapangan
plots_poly <- vect("plot_20x20_azam.shp")

# Load data raster Landsat-8/9 sebagai acuan utama (base grid)
img_raw <- rast("B FOREST B2_B3_B4_B8_B11_B12.tif")

# Ubah data DN (0-10000) menjadi Reflektans (0-1)
img <- img_raw / 10000
names(img) <- c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2")

# Samakan proyeksi SHP mengikuti proyeksi citra Landsat
if (crs(plots_poly) != crs(img)) {
  cat("Menyelaraskan CRS plot lapangan...\n")
  plots_poly <- project(plots_poly, crs(img))
}


# ==============================================================================
# 2. LOAD & RESAMPLE DATA TAMBAHAN (SENTINEL-1 & COPERNICUS BIOMASS)
# ==============================================================================

# --- A. SENTINEL-1 (VV & VH) ---
s1_raw <- rast("VV_copernicus_2025_v3.tif") 
names(s1_raw) <- c("VV", "VH")

if (crs(s1_raw) != crs(img)) {
  s1_raw <- project(s1_raw, crs(img))
}
# Resample Sentinel-1 mengikuti grid Landsat (Bilinear cocok untuk radar/reflektans)
s1_resampled <- resample(s1_raw, img, method = "bilinear")


# --- B. COPERNICUS BIOMASS (300m) ---
biomass_raster_raw <- rast("AGB_Carbon_Converted_tonCHa_300m.tif")
names(biomass_raster_raw) <- "Biomass_NASA_300m"

if (crs(biomass_raster_raw) != crs(img)) {
  cat("Menyelaraskan CRS raster biomassa...\n")
  biomass_raster_raw <- project(biomass_raster_raw, crs(img))
}
# Resample Biomassa 300m mengikuti grid Landsat menggunakan metode 'near'
biomass_raster_resampled <- resample(biomass_raster_raw, img, method = "near")


# ==============================================================================
# 3. KALKULASI REMOTE SENSING METRICS
# ==============================================================================

# Indeks Dasar
ndvi <- (img$NIR - img$Red) / (img$NIR + img$Red)
names(ndvi) <- "NDVI"

savi <- ((img$NIR - img$Red) / (img$NIR + img$Red + 0.5)) * 1.5
names(savi) <- "SAVI"

gndvi <- (img$NIR - img$Green) / (img$NIR + img$Green)
names(gndvi) <- "GNDVI"

evi <- 2.5 * ((img$NIR - img$Red) / (img$NIR + 6 * img$Red - 7.5 * img$Blue + 1))
names(evi) <- "EVI"

evi2 <- 2.5 * ((img$NIR - img$Red) / (img$NIR + 2.4 * img$Red + 1))
names(evi2) <- "EVI2"

msavi <- (2 * img$NIR + 1 - sqrt((2 * img$NIR + 1)^2 - 8 * (img$NIR - img$Red))) / 2
names(msavi) <- "MSAVI"

arvi <- (img$NIR - (2 * img$Red - img$Blue)) / (img$NIR + (2 * img$Red - img$Blue))
names(arvi) <- "ARVI"

lai <- 3.618 * ndvi - 0.118
names(lai) <- "LAI"

# Kombinasi Rata-Rata
mean_vi <- (ndvi + savi + evi + gndvi) / 4
names(mean_vi) <- "Mean_VI"

mean_ndvi_savi <- (ndvi + savi) / 2
names(mean_ndvi_savi) <- "Mean_NDVI_SAVI"

mean_gndvi_savi <- (gndvi + savi) / 2
names(mean_gndvi_savi) <- "Mean_GNDVI_SAVI"

mean_evi_savi <- (evi + savi) / 2
names(mean_evi_savi) <- "Mean_EVI_SAVI"

mean_gndvi_swir2 <- (gndvi + img$SWIR2) / 2
names(mean_gndvi_swir2) <- "Mean_GNDVI_SWIR2"

# PCA & Composite
mean_composite <- mean(img)
names(mean_composite) <- "Mean_Composite"

pca_model <- rasterPCA(img, nComp = 3)
pca_bands <- pca_model$map
names(pca_bands) <- c("PC1", "PC2", "PC3")


# ==============================================================================
# 4. PENGGABUNGAN SEMUA LAYER (STACKING) & EKSTRAKSI DATA
# ==============================================================================

# Menggabungkan seluruh data ke dalam satu SpatRaster (Semua sudah beresolusi 30m)
all_metrics <- c(
  img, 
  ndvi, savi, gndvi, evi, evi2, msavi, arvi, lai, 
  mean_composite, mean_vi, mean_ndvi_savi, mean_gndvi_savi, mean_evi_savi, mean_gndvi_swir2, 
  pca_bands,
  s1_resampled,
  biomass_raster_resampled
)

# Cek nama seluruh layer untuk memastikan tidak ada duplikat
print("Daftar layer yang siap diekstrak:")
print(names(all_metrics))

# Ekstrak nilai rata-rata pixel raster yang masuk ke dalam polygon plot lapangan
extracted_values <- terra::extract(all_metrics, plots_poly, fun = mean, na.rm = TRUE)

# Konversi properti SHP ke dataframe dan gabungkan dengan hasil ekstraksi
data_plot_df <- as.data.frame(plots_poly)
final_data <- cbind(data_plot_df, extracted_values[, -1])
final_data_clean <- na.omit(final_data)


# ==============================================================================
# 5. PENYIAPAN FORMAT WIDE & EXPORT EXCEL
# ==============================================================================

# Menyeleksi kolom respons dan seluruh prediktor (dari Blue hingga Biomass_Raster_300m)
data_wide_clean <- final_data_clean %>% 
  select(biomass, Blue:Biomass_NASA_300m)

print("Tampilan atas data wide:")
head(data_wide_clean)

# Export langsung ke format Excel
nama_file_excel <- "Data_Wide_Clean_Lengkap_azam.xlsx"
write_xlsx(data_wide_clean, path = nama_file_excel)

cat("Selesai! Data berhasil disimpan di file:", nama_file_excel, "\n")
