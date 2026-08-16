rm(list = ls())

# Standalone generation of main-manuscript Figures 1--3
# =============================================================================
# Run from the repository root with:
#   source("data analysis/exploratory_figures_seer.R")
#
# This script uses observed SEER counts, expected counts, and county covariates;
# it does not read posterior draws or run MCMC. It creates exactly three PDFs:
#   data analysis/exploratory_figures/
#     figure_1_seer_cancer_sir_maps.pdf
#     figure_2_seer_morans_i_by_distance_band.pdf
#     figure_3_seer_covariate_maps.pdf

options(stringsAsFactors = FALSE, tigris_use_cache = TRUE)

if (!all(file.exists(c("README.md", "data analysis", "data analysis/data")))) {
  stop("Run this script from the repository root.", call. = FALSE)
}

required_packages <- c(
  "ggplot2", "ggpubr", "mapproj", "maps", "RColorBrewer", "readr",
  "sf", "tigris"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    ". Run install.R first.", call. = FALSE
  )
}

data_dir <- file.path("data analysis", "data")
output_dir <- file.path("data analysis", "exploratory_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  sir = file.path(data_dir, "SIR_adjusted.csv"),
  covariates = file.path(data_dir, "covariates.csv"),
  smoking = file.path(data_dir, "smoking.csv")
)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop(
    "Missing required input file(s):\n",
    paste(" -", missing_inputs, collapse = "\n"),
    call. = FALSE
  )
}

normalize_county <- function(x) {
  x <- gsub("\u00a0", " ", as.character(x), fixed = TRUE)
  x <- sub("^CA:\\s*", "", x)
  x <- sub("\\s+(?:County|Registry)(?:\\s+\\([0-9]+\\))?$", "", x)
  tolower(trimws(x))
}

extract_fips <- function(x) {
  x <- as.character(x)
  matched <- grepl("\\([0-9]{5}\\)", x)
  result <- rep(NA_character_, length(x))
  result[matched] <- sub(".*\\(([0-9]{5})\\).*", "\\1", x[matched])
  result
}

parse_encoded_percentage <- function(x, variable_name) {
  result <- suppressWarnings(as.numeric(as.character(x)) / 100)
  if (anyNA(result)) {
    stop("Could not parse ", variable_name, " as percentages.", call. = FALSE)
  }
  result
}

sir <- suppressMessages(readr::read_csv(
  input_files[["sir"]], show_col_types = FALSE, progress = FALSE
))
covariates <- readr::read_csv(
  input_files[["covariates"]], show_col_types = FALSE, progress = FALSE
)
smoking <- utils::read.csv(
  input_files[["smoking"]], check.names = FALSE, stringsAsFactors = FALSE
)

required_sir_columns <- c(
  "State.county", "standard_ratio", "Site.code"
)
if (!all(required_sir_columns %in% names(sir))) {
  stop(
    "SIR_adjusted.csv is missing required columns: ",
    paste(setdiff(required_sir_columns, names(sir)), collapse = ", "),
    call. = FALSE
  )
}
required_covariate_columns <- c(
  "State_county", "V_Persons_age_65_ACS_2012_2016",
  "VFamiliesbelowpovertyACS201220"
)
if (!all(required_covariate_columns %in% names(covariates))) {
  stop(
    "covariates.csv is missing required columns: ",
    paste(setdiff(required_covariate_columns, names(covariates)), collapse = ", "),
    call. = FALSE
  )
}
if (ncol(smoking) < 2L) {
  stop("smoking.csv must contain county and smoking-rate columns.",
       call. = FALSE)
}

sir$fips <- extract_fips(sir$State.county)
sir$county_key <- normalize_county(sir$State.county)
sir$standard_ratio <- suppressWarnings(as.numeric(sir$standard_ratio))
if (anyNA(sir$fips) || anyNA(sir$standard_ratio)) {
  stop("Invalid county FIPS code or standardized incidence ratio.",
       call. = FALSE)
}

site_names <- c(
  lung = "Lung and Bronchus",
  esophageal = "Esophagus",
  larynx = "Larynx",
  colorectal = "Colon and Rectum"
)
site_data <- lapply(site_names, function(site_name) {
  result <- sir[sir$Site.code == site_name, , drop = FALSE]
  if (nrow(result) != 58L || anyDuplicated(result$fips)) {
    stop(
      "Expected 58 unique California counties for site '", site_name, "'.",
      call. = FALSE
    )
  }
  result
})

california_covariates <- covariates[
  substr(covariates$State_county, 1L, 2L) == "CA", , drop = FALSE
]
california_covariates$fips <- extract_fips(
  california_covariates$State_county
)
california_covariates$county_key <- normalize_county(
  california_covariates$State_county
)
if (nrow(california_covariates) != 58L ||
    anyNA(california_covariates$fips) ||
    anyDuplicated(california_covariates$fips)) {
  stop("Expected 58 unique California rows in covariates.csv.",
       call. = FALSE)
}
california_covariates$over_65_rate <- parse_encoded_percentage(
  california_covariates$V_Persons_age_65_ACS_2012_2016,
  "population-over-65 rate"
)
california_covariates$poverty_rate <- parse_encoded_percentage(
  california_covariates$VFamiliesbelowpovertyACS201220,
  "family-poverty rate"
)

smoking_county <- normalize_county(smoking[[1L]])
smoking_text <- gsub("\u00a0", "", as.character(smoking[[2L]]), fixed = TRUE)
smoking_rate <- suppressWarnings(as.numeric(sub("%.*$", "", smoking_text)))
if (length(smoking_rate) != 58L || anyNA(smoking_rate) ||
    anyDuplicated(smoking_county)) {
  stop("Expected 58 unique, numeric county smoking rates.", call. = FALSE)
}
smoking_data <- data.frame(
  county_key = smoking_county,
  smoking_rate = smoking_rate,
  stringsAsFactors = FALSE
)

# Fixed-year Census cartographic boundary geometry for Figures 1 and 3.
ca_counties_sf <- tryCatch(
  tigris::counties(state = "CA", cb = TRUE, year = 2016, class = "sf"),
  error = function(error) {
    stop(
      "Could not obtain the fixed 2016 California county geometry: ",
      conditionMessage(error),
      ". Check the internet connection or the tigris cache.",
      call. = FALSE
    )
  }
)
ca_counties_sf <- sf::st_make_valid(ca_counties_sf)
ca_counties_sf$GEOID <- as.character(ca_counties_sf$GEOID)
ca_counties_sf$county_key <- normalize_county(ca_counties_sf$NAME)
if (nrow(ca_counties_sf) != 58L || anyDuplicated(ca_counties_sf$GEOID)) {
  stop("The California geometry does not contain 58 unique counties.",
       call. = FALSE)
}

for (site in names(site_data)) {
  matched <- match(ca_counties_sf$GEOID, site_data[[site]]$fips)
  if (anyNA(matched)) {
    stop("Unmatched counties for the ", site, " SIR map.", call. = FALSE)
  }
  ca_counties_sf[[paste0("sir_", site)]] <-
    site_data[[site]]$standard_ratio[matched]
}
covariate_match <- match(
  ca_counties_sf$GEOID, california_covariates$fips
)
smoking_match <- match(ca_counties_sf$county_key, smoking_data$county_key)
if (anyNA(covariate_match) || anyNA(smoking_match)) {
  stop("At least one county could not be matched to its covariates.",
       call. = FALSE)
}
ca_counties_sf$smoking_rate <- smoking_data$smoking_rate[smoking_match]
ca_counties_sf$over_65_rate <-
  california_covariates$over_65_rate[covariate_match]
ca_counties_sf$poverty_rate <-
  california_covariates$poverty_rate[covariate_match]

palette <- RColorBrewer::brewer.pal(9L, "YlOrRd")
map_theme <- ggplot2::theme(
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  panel.background = ggplot2::element_blank(),
  axis.line = ggplot2::element_blank(),
  axis.title = ggplot2::element_blank(),
  axis.text = ggplot2::element_blank(),
  axis.ticks = ggplot2::element_blank(),
  plot.title = ggplot2::element_text(
    size = 20, face = "bold", hjust = 0.5
  ),
  legend.text = ggplot2::element_text(size = 12),
  legend.key.size = grid::unit(0.6, "cm"),
  legend.title = ggplot2::element_blank()
)

county_map <- function(fill_variable, title) {
  ggplot2::ggplot(ca_counties_sf) +
    ggplot2::geom_sf(
      ggplot2::aes(fill = .data[[fill_variable]]),
      color = "black", linewidth = 0.20, alpha = 0.6
    ) +
    ggplot2::scale_fill_gradientn(colours = palette, na.value = "grey90") +
    ggplot2::ggtitle(title) +
    map_theme
}

# RESULT: Main Figure 1 - observed cancer SIR maps ----------------------------
figure_1 <- ggpubr::ggarrange(
  county_map("sir_lung", "Lung"),
  county_map("sir_esophageal", "Esophageal"),
  county_map("sir_larynx", "Larynx"),
  county_map("sir_colorectal", "Colorectal"),
  nrow = 1L, ncol = 4L
)
figure_1_file <- file.path(
  output_dir, "figure_1_seer_cancer_sir_maps.pdf"
)
ggplot2::ggsave(
  figure_1_file, figure_1, width = 16, height = 5, units = "in",
  device = grDevices::cairo_pdf
)

# RESULT: Main Figure 2 - Moran's I by distance band --------------------------
# Reproduce the original distance-band calculation using the maps package's
# California county centroids and an Albers projection.
ca_map <- maps::map(
  "county", "california", fill = TRUE, plot = FALSE
)
ca_map_sf <- sf::st_as_sf(ca_map)
sf::st_crs(ca_map_sf) <- NA
ca_map_sf$county_key <- normalize_county(
  sub("^california,", "", ca_map_sf$ID)
)
if (nrow(ca_map_sf) != 58L || anyDuplicated(ca_map_sf$county_key)) {
  stop("The maps package did not return 58 unique California counties.",
       call. = FALSE)
}

map_coordinates <- suppressWarnings(
  sf::st_coordinates(sf::st_centroid(sf::st_geometry(ca_map_sf)))
)
latitude_range <- round(stats::quantile(
  map_coordinates[, 2L], c(0.25, 0.75)
))
albers <- mapproj::mapproject(
  map_coordinates[, 1L], map_coordinates[, 2L],
  projection = "albers", param = latitude_range
)
distance_matrix <- as.matrix(stats::dist(cbind(albers$x, albers$y)))

moran_i <- function(y, adjacency) {
  centered <- y - mean(y)
  weight_sum <- sum(adjacency)
  if (weight_sum <= 0) {
    return(NA_real_)
  }
  length(y) * sum(adjacency * tcrossprod(centered)) /
    (weight_sum * sum(centered^2))
}

moran_rows <- list()
for (site_index in seq_along(site_data)) {
  site <- names(site_data)[[site_index]]
  county_match <- match(
    ca_map_sf$county_key, site_data[[site]]$county_key
  )
  if (anyNA(county_match)) {
    stop("Unmatched counties in the Moran's-I calculation for ", site, ".",
         call. = FALSE)
  }
  sir_values <- site_data[[site]]$standard_ratio[county_match]
  moran_values <- vapply(seq_len(11L), function(lag) {
    adjacency <- 1L * (
      distance_matrix <= 0.01 * lag &
        distance_matrix > 0.01 * (lag - 1L)
    )
    diag(adjacency) <- 0L
    moran_i(sir_values, adjacency)
  }, numeric(1L))
  moran_rows[[site_index]] <- data.frame(
    cancer = c(
      lung = "Lung", esophageal = "Esophageal", larynx = "Larynx",
      colorectal = "Colorectum"
    )[[site]],
    distance_band = seq_len(11L),
    moran_i = moran_values,
    stringsAsFactors = FALSE
  )
}
moran_data <- do.call(rbind, moran_rows)
moran_data$cancer <- factor(
  moran_data$cancer,
  levels = c("Lung", "Esophageal", "Larynx", "Colorectum")
)
figure_2 <- ggplot2::ggplot(
  moran_data, ggplot2::aes(distance_band, moran_i)
) +
  ggplot2::geom_point(size = 3) +
  ggplot2::ylab("Moran's I") +
  ggplot2::xlab("Distance band") +
  ggplot2::facet_wrap(~cancer, nrow = 1L, ncol = 4L) +
  ggplot2::theme_bw() +
  ggplot2::scale_x_continuous(breaks = seq_len(11L)) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 16),
    axis.text = ggplot2::element_text(size = 11),
    strip.text = ggplot2::element_text(size = 13)
  )
figure_2_file <- file.path(
  output_dir, "figure_2_seer_morans_i_by_distance_band.pdf"
)
ggplot2::ggsave(
  figure_2_file, figure_2, width = 14, height = 4.5, units = "in",
  device = grDevices::cairo_pdf
)

# RESULT: Main Figure 3 - observed county covariate maps ----------------------
figure_3 <- ggpubr::ggarrange(
  county_map("smoking_rate", "Smoking rate"),
  county_map("over_65_rate", "Over 65 rate"),
  county_map("poverty_rate", "Below poverty rate"),
  nrow = 1L, ncol = 3L
)
figure_3_file <- file.path(
  output_dir, "figure_3_seer_covariate_maps.pdf"
)
ggplot2::ggsave(
  figure_3_file, figure_3, width = 13, height = 5, units = "in",
  device = grDevices::cairo_pdf
)

figure_files <- c(figure_1_file, figure_2_file, figure_3_file)
if (!all(file.exists(figure_files))) {
  stop("One or more exploratory figures were not created.", call. = FALSE)
}

cat("SEER exploratory figures created:\n")
cat(
  paste0(
    "  - ",
    vapply(
      figure_files,
      normalizePath,
      character(1L),
      winslash = "/",
      mustWork = TRUE
    )
  ),
  sep = "\n"
)
cat("\n")
