Ecological Drivers of Stickleback Morphological Variation
================
John Berini
2026-06-30

### Bring in our data and load the packages

``` r
# Core Data Science & Visualization
library(tidyverse)   # Automatically loads dplyr, ggplot2, forcats, and lubridate
library(tidytext)

# Specialized Machine Learning & Statistics
library(MixRF)
library(randomForest)
library(caret)
library(ranger)      
library(car)         
library(vegan)       
library(lmerTest)    
library(Hmisc)       
library(ggh4x)

# GitHub Safety Guard: Ensure local directories exist to prevent saving errors
if(!dir.exists("figures")) dir.create("figures")

# 1. Load Primary Datasets
sticks <- read.csv("data/stickle_morph.csv", header=TRUE) 
names(sticks)[names(sticks) == "Population"] <- "lake_code"
sticks$lake_code <- trimws(tolower(sticks$lake_code))

covs <- read.csv("data/rf_dat.csv", header=TRUE) 
covs$lake_code <- trimws(tolower(covs$lake_code))

sonde <- read.csv("data/VancouverSondeData2023.csv", header=TRUE)
zoop_raw <- read.csv("data/ZoopIDDataUpdated.csv", header=TRUE)
lake_mapping <- read.csv("data/lake_id_mapping.csv", header=TRUE)

# 2. Build the lookup dictionary
clean_lookup <- lake_mapping %>%
  mutate(
    join_key = trimws(tolower(original_value)),
    LakeID   = trimws(tolower(LakeID))
  ) %>%
  select(join_key, LakeID) %>%
  distinct()

# Map 'sticks' and 'covs' through the lookup dictionary
sticks <- sticks %>%
  mutate(join_key = lake_code) %>%
  left_join(clean_lookup, by = "join_key") %>%
  mutate(lake_code = coalesce(LakeID, lake_code)) %>%
  select(-join_key, -LakeID)

covs <- covs %>%
  mutate(join_key = lake_code) %>%
  left_join(clean_lookup, by = "join_key") %>%
  mutate(lake_code = coalesce(LakeID, lake_code)) %>%
  select(-join_key, -LakeID)
```

\#bridge and join datasets

``` r
sonde_summary <- sonde %>%
  mutate(join_key = trimws(tolower(Site))) %>%
  left_join(clean_lookup, by = "join_key") %>%
  group_by(LakeID, Date) %>%
  summarise(
    Mean_Temp_C      = mean(TemperatureCelcius, na.rm = TRUE),
    SD_Temp_C        = sd(TemperatureCelcius, na.rm = TRUE),
    Mean_DO_Percent  = mean(PercentDissolvedOxygen, na.rm = TRUE),
    SD_DO_Percent    = sd(PercentDissolvedOxygen, na.rm = TRUE),
    Mean_DO_mgL      = mean(DissolvedOxygenmgL, na.rm = TRUE),
    SD_DO_mgL        = sd(DissolvedOxygenmgL, na.rm = TRUE),
    Mean_SPC_uScm    = mean(SPCuSpercm, na.rm = TRUE),
    SD_SPC_uScm      = sd(SPCuSpercm, na.rm = TRUE),
    Mean_pH          = mean(pH, na.rm = TRUE),
    SD_pH            = sd(pH, na.rm = TRUE),
    Mean_Phyco_RFU   = mean(PhycocyaninRFU, na.rm = TRUE),
    SD_Phyco_RFU     = sd(PhycocyaninRFU, na.rm = TRUE),
    Mean_Chl_RFU     = mean(ChlorophyllRFU, na.rm = TRUE),
    SD_Chl_RFU       = sd(ChlorophyllRFU, na.rm = TRUE),
    Mean_fDOM_RFU    = mean(fDOMRFU, na.rm = TRUE),
    SD_fDOM_RFU      = sd(fDOMRFU, na.rm = TRUE),
    .groups = "drop"
  )

zoop_standardized <- zoop_raw %>%
  mutate(join_key = trimws(tolower(LakeName))) %>%
  left_join(clean_lookup, by = "join_key") %>%
  select(-join_key)

master_environmental_file <- sonde_summary %>%
  rename(lake_code = LakeID) %>%
  group_by(lake_code) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

zoop_taxa_cols <- intersect(
  c("Calanoids", "Cyclopoids", "Nauplii", "Daphnia", "Bosmina", "Chydorus", 
    "Diaphanosoma", "Holopedium", "Chaoborus", "WaterMites", "Brachionus.Keratella", 
    "ConochilusColonies", "Asplanchna", "Ciliate", "Amoeboid", "Rotifers..other.", "Polyphemidae"),
  colnames(zoop_standardized)
)

zoop_scores <- data.frame(lake_code = character(), Zoop_Community_PC1 = numeric(), Zoop_Community_PC2 = numeric(), stringsAsFactors = FALSE)

zoop_pca_data <- zoop_standardized %>%
  filter(!is.na(LakeID) & LakeID != "") %>% 
  mutate(across(all_of(zoop_taxa_cols), ~ as.numeric(as.character(.x)))) %>%
  group_by(LakeID) %>%
  summarise(across(all_of(zoop_taxa_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
```

    ## Warning: There were 2 warnings in `mutate()`.
    ## The first warning was:
    ## ℹ In argument: `across(all_of(zoop_taxa_cols), ~as.numeric(as.character(.x)))`.
    ## Caused by warning:
    ## ! NAs introduced by coercion
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.

``` r
if(nrow(zoop_pca_data) > 0) {
  zoop_pca_df <- as.data.frame(zoop_pca_data)
  rownames(zoop_pca_df) <- zoop_pca_df$LakeID
  zoop_pca_df$LakeID <- NULL
  
  zoop_pca_df[is.na(zoop_pca_df)] <- 0
  row_totals <- rowSums(zoop_pca_df)
  zoop_pca_df <- zoop_pca_df[row_totals > 0, ]
  
  if(nrow(zoop_pca_df) > 0) {
    zoop_pca <- vegan::rda(vegan::decostand(zoop_pca_df, method = "hellinger"))
    
    # Extract and format variance explained by the zooplankton axes themselves
    zoop_importance <- summary(zoop_pca)$cont$importance
    zoop_variance_df <- data.frame(
      Axis = c("Zoop_PC1", "Zoop_PC2"),
      Eigenvalue = zoop_importance["Eigenvalue", c(1, 2)],
      Proportion_Explained = zoop_importance["Proportion Explained", c(1, 2)],
      Cumulative_Explained = zoop_importance["Cumulative Proportion", c(1, 2)]
    )
    write.csv(zoop_variance_df, "zooplankton_pca_variance_summary.csv", row.names = FALSE)
    
    extracted_scores <- as.data.frame(vegan::scores(zoop_pca, choices = c(1, 2), display = "sites"))
    colnames(extracted_scores) <- c("Zoop_Community_PC1", "Zoop_Community_PC2")
    
    zoop_scores <- extracted_scores
    zoop_scores$lake_code <- rownames(extracted_scores)
    rownames(zoop_scores) <- NULL
  }
}
```

covariate clearning and imputation log

``` r
combined_covs <- covs %>%
  left_join(master_environmental_file, by = "lake_code") %>%
  left_join(zoop_scores, by = "lake_code")

if("bedrock" %in% colnames(combined_covs)) combined_covs$bedrock <- as.factor(combined_covs$bedrock)
if("watershed" %in% colnames(combined_covs)) combined_covs$watershed <- as.factor(combined_covs$watershed)

# Count total numeric cells across all columns to find the denominator
total_numeric_cells <- combined_covs %>% 
  select(where(is.numeric)) %>% 
  is.na() %>% 
  length()

# Count how many cells are missing (NAs)
total_imputed_values <- combined_covs %>% 
  select(where(is.numeric)) %>% 
  is.na() %>% 
  sum()

# Calculate percent missingness
pct_imputed <- (total_imputed_values / total_numeric_cells) * 100

# Print a summary statement
message(paste0("\n>>> [Reviewer Metric]: Imputed ", total_imputed_values, 
               " missing numeric cells out of ", total_numeric_cells, 
               " total cells (", round(pct_imputed, 2), "% of the numeric matrix).\n"))
```

    ## 
    ## >>> [Reviewer Metric]: Imputed 144 missing numeric cells out of 14006 total cells (1.03% of the numeric matrix).

``` r
# Save imputation log to directory
imputation_log <- data.frame(
  Metric = c("Total_Numeric_Cells", "Total_Imputed_Cells", "Percent_Imputed"),
  Value  = c(total_numeric_cells, total_imputed_values, round(pct_imputed, 4)),
  stringsAsFactors = FALSE
)
write.csv(imputation_log, "metrics_imputation_summary.csv", row.names = FALSE)

# 1. Base lake-level imputation
combined_covs_clean <- combined_covs %>%
  mutate(across(where(is.numeric), ~ {
    col_num <- as.numeric(.x)
    replace_na(col_num, mean(col_num, na.rm = TRUE))
  })) %>%
  mutate(across(where(is.factor), ~ fct_na_value_to_level(.x, level = "Missing")))

# 2. Expand out clean covariates to individual fish observations
sticks_rf_dat <- sticks %>%
  left_join(combined_covs_clean, by = "lake_code") %>%
  drop_na(Standard.Length:Fin.Width)

# Isolate response matrix (morphology traits)
sticks_morph <- sticks_rf_dat %>% select(Standard.Length:Fin.Width)

# 3. Build predictor matrix for Random Forest
rf_dat <- sticks_rf_dat %>%
  select(-all_of(names(sticks_morph)), -lake_code, -any_of("lake_name"), -any_of(zoop_taxa_cols)) %>%
  select(where(~ is.numeric(.) | is.factor(.)))

# 4. Omit complex multi-level factor variables that bias RF 
rf_dat <- rf_dat %>%
  select(where(~ !is.factor(.x) || length(levels(.x)) <= 3))

# 5. Cleanup final predictor matrix
rf_dat <- rf_dat %>%
  select(where(~ any(!is.na(.)))) %>%
  mutate(across(where(is.numeric), ~ replace_na(as.numeric(.x), mean(as.numeric(.x), na.rm = TRUE)))) %>%
  mutate(across(where(is.factor), ~ fct_na_value_to_level(.x, level = "Missing")))

# 6. Drop columns with zero variance
rf_dat <- rf_dat %>%
  select(where(~ !is.numeric(.x) || (!is.na(sd(.x, na.rm = TRUE)) && sd(.x, na.rm = TRUE) != 0)))
```

Random forest variable importance loop

``` r
summary_list <- list()
set.seed(42) 

forest_df <- cbind(rf_dat, sticks_morph)
predictor_names <- colnames(rf_dat)

for (t_name in names(sticks_morph)) {
  rf_fit <- ranger::ranger(
    formula = as.formula(paste(t_name, "~ .")),
    data = forest_df[, c(predictor_names, t_name)],
    importance = "permutation", 
    num.trees = 1000,
    mtry = ceiling(length(predictor_names) / 3)
  )
  
  raw_imp <- rf_fit$variable.importance
  raw_imp[raw_imp < 0] <- 0 
  
  rel_importance = if(sum(raw_imp) > 0) (raw_imp / sum(raw_imp)) * 100 else rep(0, length(raw_imp))
  
  summary_list[[t_name]] <- data.frame(
    trait = t_name,
    variable = names(raw_imp),
    mean_rel_importance = as.numeric(rel_importance),
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, summary_list)
```

Variable Selection & Linear Model Filtering

``` r
n_vars <- ncol(rf_dat)
imp_threshold <- 100 / n_vars

top_vars_dynamic <- summary_table %>%
  group_by(trait) %>%
  filter(mean_rel_importance > imp_threshold) %>%
  slice_max(order_by = mean_rel_importance, n = 10, with_ties = FALSE) %>% 
  ungroup()

final_plot_data <- list()

if(nrow(top_vars_dynamic) > 0) {
  for(i in seq_len(nrow(top_vars_dynamic))) {
    t_name <- top_vars_dynamic$trait[i]
    v_name <- top_vars_dynamic$variable[i]
     
    if(!v_name %in% colnames(rf_dat)) next
     
    # Clean temporary dataset for this specific iteration
    loop_df <- data.frame(
      y_val = sticks_rf_dat[[t_name]],
      x_val = rf_dat[[v_name]],
      stringsAsFactors = FALSE
    ) %>% 
      drop_na(y_val, x_val)
    
    # Data volume check --> Ensure we have data to fit a line
    if(nrow(loop_df) < 5) next
    
    if(is.factor(loop_df$x_val)) {
      fit <- lm(y_val ~ x_val, data = loop_df)
      anova_summary <- anova(fit)
      stats <- as.data.frame(summary(fit)$coefficients)
      
      p_val   <- anova_summary["x_val", "Pr(>F)"]
      est_val <- stats[2, "Estimate"]
      std_err <- stats[2, "White-level Std. Error"] # adjusted dynamically downstream if required
      std_err <- stats[2, "Std. Error"]
      t_val   <- stats[2, "t value"]
    } else {
      loop_df$x_scaled <- scale(as.numeric(loop_df$x_val))
      fit <- lm(y_val ~ x_scaled, data = loop_df)
      stats <- as.data.frame(summary(fit)$coefficients)
      
      p_val   <- stats[2, "Pr(>|t|)"]
      est_val <- stats[2, "Estimate"]
      std_err <- stats[2, "Std. Error"]
      t_val   <- stats[2, "t value"]
    }
     
    if(!is.na(p_val) && p_val < 0.05) { 
      stats_row <- data.frame(
        estimate  = est_val,
        std_error = std_err,
        t_value   = t_val,
        p_value   = p_val,
        variable  = v_name,
        trait     = t_name,
        is_factor = is.factor(loop_df$x_val),
        stringsAsFactors = FALSE
      )
      
      stats_row$sig_tier = case_when(
        stats_row$p_value < 0.001 ~ "***",
        stats_row$p_value < 0.01  ~ "**",
        stats_row$p_value < 0.05  ~ "*"
      )
      
      stats_row$category <- case_when(
        grepl("^Mean_|^SD_", v_name) ~ "Water Column (Sonde)",
        grepl("Zoop_Community_", v_name) ~ "Zooplankton Food Supply",
        TRUE                          ~ "Landscape & Climate"
      )
      
      final_plot_data[[length(final_plot_data) + 1]] <- stats_row
    }
  }
}

# Compile the plot dataset
raw_plot_df <- if(length(final_plot_data) > 0) do.call(rbind, final_plot_data) else data.frame()

all_traits <- names(sticks_morph)
sig_traits = if(nrow(raw_plot_df) > 0) unique(raw_plot_df$trait) else c()
missing_traits <- setdiff(all_traits, sig_traits) 

if(length(missing_traits) > 0) {
  empty_df <- data.frame(
    estimate = 0, std_error = 0, t_value = 0, p_value = 1,
    variable = "(No significant covariates)", trait = missing_traits, 
    sig_tier = "", category = "None", is_factor = FALSE
  )
  plot_df_final <- rbind(raw_plot_df, empty_df)
} else {
  plot_df_final <- raw_plot_df
}

plot_df_final <- plot_df_final %>%
  mutate(
    variable = as.character(variable),
    plot_sort_weight = if_else(is.na(estimate), 0, estimate)
  )
```

Visualize significant environmental drivers

``` r
p_coeff <- ggplot(plot_df_final, aes(y = tidytext::reorder_within(variable, plot_sort_weight, trait), 
                                     x = estimate, 
                                     color = category)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray60", linewidth = 0.6) +
  
  geom_errorbar(data = subset(plot_df_final, variable != "(No significant covariates)"),
                aes(xmin = estimate - (1.96 * std_error), 
                    xmax = estimate + (1.96 * std_error)), 
                width = 0.3, linewidth = 0.6) +
  
  geom_point(data = subset(plot_df_final, variable != "(No significant covariates)"), 
             size = 1.8) +
  
  geom_text(data = subset(plot_df_final, variable != "(No significant covariates)"),
            aes(label = sig_tier, 
                x = if_else(estimate > 0, estimate + (1.96 * std_error), estimate - (1.96 * std_error)),
                hjust = if_else(estimate > 0, -0.3, 1.3)), 
            vjust = 0.4, size = 3, fontface = "bold", show.legend = FALSE) +
  
  geom_text(data = subset(plot_df_final, variable == "(No significant covariates)"),
            aes(y = tidytext::reorder_within(variable, plot_sort_weight, trait), x = 0, 
                label = "No significant relationships found"),
            color = "firebrick", fontface = "italic", size = 3, show.legend = FALSE) +
  
  tidytext::scale_y_reordered() +
  facet_wrap2(~ trait, scales = "free", ncol = 2, axes = "all", remove_labels = "none") + 
  
  scale_color_manual(
    values = c(
      "Landscape & Climate"     = "#2c3e50",  
      "Water Column (Sonde)"    = "#2980b9",  
      "Zooplankton Food Supply" = "#27ae60",  
      "None"                    = "gray60"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Significant Drivers of Stickleback Morphology",
    subtitle = "Standardized Fixed-Effects Slopes (Beta \u00b1 95% CI) | Linear Regression Architecture",
    y = "Model Predictors", 
    x = "Standardized Effect Size (Continuous & Factor Indicators)",
    color = "Covariate Source",
    caption = "Significance levels: * p < 0.05, ** p < 0.01, *** p < 0.001"
  ) +
  theme_minimal(base_size = 9) + 
  theme(
    panel.border       = element_rect(color = "gray85", fill = NA, linewidth = 0.5),
    strip.background   = element_rect(fill = "gray96", color = "gray85"),
    strip.text         = element_text(face = "bold", size = 8.5), 
    panel.spacing.x    = unit(1.0, "lines"),                     
    panel.spacing.y    = unit(0.6, "lines"),                     
    axis.text.y        = element_text(size = 7, color = "gray20"), 
    axis.text.x        = element_text(size = 7.5),
    axis.title         = element_text(size = 9, face = "bold"),
    legend.position    = "bottom",
    legend.text        = element_text(size = 8)
  )

print(p_coeff)
```

![](Ecological_Drivers_of_Stickleback_Morphological_Variation_GitHUB_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggsave("figures/Figure 1. Stickleback Drivers.pdf", 
       plot = p_coeff, 
       width = 8.5, 
       height = 11, 
       units = "in", 
       dpi = 300)
```

Build summary table of results

``` r
table_list <- list()
valid_models <- data.frame()
if (nrow(plot_df_final) > 0 && "variable" %in% colnames(plot_df_final)) {
  valid_models <- subset(plot_df_final, variable != "(No significant covariates)")
}

if(nrow(valid_models) > 0) {
  for(i in 1:nrow(valid_models)) {
    t_name   <- valid_models$trait[i]
    v_name   <- valid_models$variable[i]
    cat_name <- valid_models$category[i]
    is_fact  <- valid_models$is_factor[i]
    loop_df <- data.frame(
      y_val = sticks_rf_dat[[t_name]],
      x_val = rf_dat[[v_name]],
      stringsAsFactors = FALSE
    ) %>% 
      drop_na(y_val, x_val)
    
    if(is_fact) {
      fit <- lm(y_val ~ x_val, data = loop_df)
      coefs <- as.data.frame(summary(fit)$coefficients)
      an_res <- anova(fit)
      
      b_val      <- coefs[2, "Estimate"]
      se_val     <- coefs[2, "White-level Std. Error"] # validation placeholder adjustment if needed
      se_val     <- coefs[2, "Std. Error"]
      t_stat_val <- coefs[2, "t value"]
      p_val      <- coefs[2, "Pr(>|t|)"]
      r2_val     <- summary(fit)$r.squared
      df_val     <- as.character(an_res["Residuals", "Df"])
    } else {
      loop_df$x_scaled <- scale(as.numeric(loop_df$x_val))
      fit <- lm(y_val ~ x_scaled, data = loop_df)
      coefs <- as.data.frame(summary(fit)$coefficients)
      an_res <- anova(fit)
      
      b_val      <- coefs[2, "Estimate"]
      se_val     <- coefs[2, "Std. Error"]
      t_stat_val <- coefs[2, "t value"]
      p_val      <- coefs[2, "Pr(>|t|)"]
      r2_val     <- summary(fit)$r.squared
      df_val     <- as.character(an_res["Residuals", "Df"])
    }
    
    table_list[[length(table_list) + 1]] <- data.frame(
      Trait     = t_name,
      Category  = cat_name,
      Predictor = v_name,
      Beta      = if(is.na(b_val)) "-" else as.numeric(b_val),
      SE        = if(is.na(se_val)) "-" else as.numeric(se_val),
      t_stat    = if(is.na(t_stat_val)) "-" else as.numeric(t_stat_val),
      p_val     = as.numeric(p_val),
      R2        = if(is.na(r2_val)) "-" else as.numeric(r2_val), 
      DF        = df_val,
      stringsAsFactors = FALSE
    )
  }
  pub_table_data <- do.call(rbind, table_list)
} else {
  pub_table_data <- data.frame(
    Trait = character(), Category = character(), Predictor = character(),
    Beta = character(), SE = character(), t_stat = character(), p_val = numeric(), 
    R2 = character(), DF = character(), stringsAsFactors = FALSE
  )
}

# Flags traits that lack statistically significant drivers
missing_traits <- setdiff(names(sticks_morph), unique(pub_table_data$Trait))
if(length(missing_traits) > 0) {
  empty_rows <- data.frame(
    Trait     = missing_traits, 
    Category  = "None", 
    Predictor = "(No significant predictors found)",
    Beta = NA_real_, 
    SE = NA_real_, 
    t_stat = NA_real_, 
    p_val = NA_real_, 
    R2 = NA_real_, 
    DF = "(No Data)",
    stringsAsFactors = FALSE
  )
  pub_table_data <- dplyr::bind_rows(pub_table_data, empty_rows)
}

pub_table_formatted <- pub_table_data %>%
  mutate(
    Beta = case_when(Beta == "-" ~ "-", TRUE ~ sprintf("%.3f", as.numeric(Beta))),
    SE = case_when(SE == "-" ~ "-", TRUE ~ sprintf("%.3f", as.numeric(SE))),
    t_stat = case_when(t_stat == "-" ~ "-", TRUE ~ sprintf("%.2f", as.numeric(t_stat))),
    R2 = case_when(R2 == "-" ~ "-", TRUE ~ sprintf("%.3f", as.numeric(R2))),
    p_val = case_when(
      is.na(p_val) ~ "-",
      as.numeric(p_val) < 0.001 ~ "< 0.001",
      TRUE                      ~ sprintf("%.3f", as.numeric(p_val))
    )
  )

knitr::kable(pub_table_formatted, caption = "Table 1. Linear regression models evaluating environmental drivers.")
```

| Trait | Category | Predictor | Beta | SE | t_stat | p_val | R2 | DF |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| Armor.Plate.N | Water Column (Sonde) | SD_DO_Percent | -1.039 | 0.336 | -3.09 | 0.005 | 0.268 | 26 |
| Armor.Plate.N | Landscape & Climate | mean_d2_ogdc | 1.011 | 0.340 | 2.98 | 0.006 | 0.254 | 26 |
| Armor.Plate.N | Landscape & Climate | max_depth_est | 1.422 | 0.277 | 5.13 | \< 0.001 | 0.503 | 26 |
| Armor.Plate.N | Landscape & Climate | mean_spp1_pc | -1.029 | 0.338 | -3.05 | 0.005 | 0.263 | 26 |
| Armor.Plate.N | Landscape & Climate | mean_DD_0_05 | 1.069 | 0.333 | 3.21 | 0.004 | 0.284 | 26 |
| Armor.Plate.N | Landscape & Climate | perimeter | 1.404 | 0.281 | 5.00 | \< 0.001 | 0.490 | 26 |
| Body.Depth | Landscape & Climate | mean_sp3_cvr | 0.262 | 0.098 | 2.67 | 0.013 | 0.215 | 26 |
| Body.Depth | Landscape & Climate | mean_sa_proj_ht | 0.294 | 0.095 | 3.11 | 0.004 | 0.271 | 26 |
| CP.Depth | Water Column (Sonde) | SD_DO_Percent | 0.066 | 0.015 | 4.35 | \< 0.001 | 0.421 | 26 |
| CP.Depth | Landscape & Climate | mean_RH_sp | 0.045 | 0.018 | 2.49 | 0.019 | 0.193 | 26 |
| CP.Depth | Landscape & Climate | mean_d2_paved | -0.055 | 0.017 | -3.26 | 0.003 | 0.290 | 26 |
| CP.Depth | Water Column (Sonde) | Mean_fDOM_RFU | 0.057 | 0.017 | 3.42 | 0.002 | 0.310 | 26 |
| CP.Depth | Water Column (Sonde) | Mean_pH | -0.041 | 0.018 | -2.23 | 0.035 | 0.160 | 26 |
| Eye.Diameter | Landscape & Climate | mean_d2_lake | -0.046 | 0.018 | -2.61 | 0.015 | 0.208 | 26 |
| Eye.Diameter | Landscape & Climate | mean_CMI06 | -0.045 | 0.018 | -2.54 | 0.017 | 0.199 | 26 |
| Eye.Diameter | Landscape & Climate | mean_Tmax_sp | -0.038 | 0.018 | -2.10 | 0.045 | 0.145 | 26 |
| Eye.Diameter | Landscape & Climate | mean_Tave10 | -0.055 | 0.017 | -3.30 | 0.003 | 0.296 | 26 |
| Eye.Diameter | Landscape & Climate | mean_RH02 | 0.052 | 0.017 | 3.11 | 0.005 | 0.271 | 26 |
| Eye.Diameter | Landscape & Climate | mean_Eref12 | -0.054 | 0.017 | -3.29 | 0.003 | 0.294 | 26 |
| Eye.Diameter | Landscape & Climate | mean_Tmax10 | -0.055 | 0.016 | -3.36 | 0.002 | 0.303 | 26 |
| Eye.Diameter | Landscape & Climate | mean_DD_18_10 | 0.054 | 0.017 | 3.29 | 0.003 | 0.294 | 26 |
| Eye.Diameter | Landscape & Climate | mean_RH11 | 0.058 | 0.016 | 3.59 | 0.001 | 0.331 | 26 |
| Eye.Diameter | Landscape & Climate | mean_DD5_09 | -0.050 | 0.017 | -2.90 | 0.008 | 0.244 | 26 |
| Fin.Width | Landscape & Climate | mean_RH03 | 0.286 | 0.090 | 3.17 | 0.004 | 0.279 | 26 |
| Fin.Width | Landscape & Climate | mean_RH_at | 0.297 | 0.089 | 3.34 | 0.003 | 0.300 | 26 |
| Fin.Width | Landscape & Climate | mean_Tmax09 | -0.285 | 0.090 | -3.15 | 0.004 | 0.276 | 26 |
| Fin.Width | Landscape & Climate | mean_RH | 0.280 | 0.091 | 3.07 | 0.005 | 0.266 | 26 |
| Fin.Width | Landscape & Climate | mean_Eref10 | -0.223 | 0.097 | -2.30 | 0.030 | 0.169 | 26 |
| Fin.Width | Landscape & Climate | mean_Eref09 | -0.286 | 0.090 | -3.17 | 0.004 | 0.278 | 26 |
| Fin.Width | Landscape & Climate | mean_Eref_at | -0.275 | 0.091 | -3.01 | 0.006 | 0.258 | 26 |
| Gill.raker.1.length | Landscape & Climate | mean_vert_comp | -0.056 | 0.018 | -3.09 | 0.005 | 0.269 | 26 |
| Gill.raker.1.length | Landscape & Climate | mean_basal_area | 0.065 | 0.017 | 3.84 | \< 0.001 | 0.362 | 26 |
| Gill.raker.1.length | Landscape & Climate | mean_spp1_pc | -0.045 | 0.019 | -2.34 | 0.027 | 0.174 | 26 |
| Gill.raker.1.length | Landscape & Climate | mean_d2_ancient | 0.052 | 0.018 | 2.82 | 0.009 | 0.234 | 26 |
| Gill.raker.N | Landscape & Climate | mean_d2_lake | -0.446 | 0.181 | -2.46 | 0.021 | 0.189 | 26 |
| Gill.raker.N | Water Column (Sonde) | Mean_fDOM_RFU | -0.584 | 0.165 | -3.53 | 0.002 | 0.324 | 26 |
| Gill.raker.N | Landscape & Climate | mean_DD_0_11 | 0.405 | 0.185 | 2.19 | 0.038 | 0.156 | 26 |
| Gill.raker.N | Landscape & Climate | mean_DEM_vi | 0.451 | 0.181 | 2.49 | 0.019 | 0.193 | 26 |
| Head.Length | Landscape & Climate | mean_RH02 | 0.163 | 0.057 | 2.88 | 0.008 | 0.242 | 26 |
| Head.Length | Landscape & Climate | mean_DD5_11 | -0.128 | 0.060 | -2.14 | 0.042 | 0.150 | 26 |
| Head.Length | Landscape & Climate | mean_Tmax02 | -0.141 | 0.059 | -2.40 | 0.024 | 0.182 | 26 |
| Head.Length | Landscape & Climate | mean_RH01 | 0.160 | 0.057 | 2.82 | 0.009 | 0.234 | 26 |
| Head.Length | Landscape & Climate | mean_Eref10 | -0.136 | 0.059 | -2.29 | 0.030 | 0.168 | 26 |
| Head.Length | Landscape & Climate | mean_DD5_02 | -0.131 | 0.060 | -2.20 | 0.037 | 0.157 | 26 |
| Head.Length | Landscape & Climate | mean_crwn_clsre | 0.128 | 0.060 | 2.14 | 0.042 | 0.150 | 26 |
| Opercular.4bar.KT | Landscape & Climate | mean_CMD08 | 0.160 | 0.059 | 2.71 | 0.012 | 0.221 | 26 |
| Opercular.4bar.KT | Landscape & Climate | mean_CMI08 | -0.159 | 0.059 | -2.68 | 0.013 | 0.216 | 26 |
| Opercular.4bar.KT | Landscape & Climate | mean_Tmax12 | -0.136 | 0.061 | -2.22 | 0.036 | 0.159 | 26 |
| Opercular.4bar.KT | Landscape & Climate | mean_DD18_sp | 0.157 | 0.059 | 2.64 | 0.014 | 0.211 | 26 |
| Opercular.4bar.KT | Landscape & Climate | mean_RH03 | 0.168 | 0.058 | 2.88 | 0.008 | 0.242 | 26 |
| Opercular.4bar.KT | Landscape & Climate | mean_d2_lake | -0.157 | 0.059 | -2.64 | 0.014 | 0.212 | 26 |
| Opercular.4bar.KT | Landscape & Climate | mean_CMI10 | -0.160 | 0.059 | -2.71 | 0.012 | 0.220 | 26 |
| Opercular.4bar.KT | Landscape & Climate | mean_CMI07 | -0.141 | 0.061 | -2.31 | 0.029 | 0.171 | 26 |
| Opercular.4bar.KT | Landscape & Climate | mean_Tmax_wt | -0.141 | 0.061 | -2.31 | 0.029 | 0.170 | 26 |
| Snout.Length | Landscape & Climate | mean_crwn_clsre | 0.094 | 0.043 | 2.19 | 0.038 | 0.155 | 26 |
| Snout.Length | Landscape & Climate | mean_basal_area | 0.108 | 0.042 | 2.59 | 0.015 | 0.205 | 26 |
| Standard.Length | Landscape & Climate | mean_CMD05 | 2.790 | 0.788 | 3.54 | 0.002 | 0.326 | 26 |
| Standard.Length | Landscape & Climate | mean_PPT09 | -3.118 | 0.739 | -4.22 | \< 0.001 | 0.407 | 26 |
| Standard.Length | Landscape & Climate | mean_CMI09 | -3.098 | 0.742 | -4.18 | \< 0.001 | 0.402 | 26 |
| Standard.Length | Landscape & Climate | mean_DD_18_sp | 2.203 | 0.856 | 2.57 | 0.016 | 0.203 | 26 |
| Standard.Length | Landscape & Climate | mean_PPT05 | -2.874 | 0.776 | -3.71 | 0.001 | 0.346 | 26 |
| Standard.Length | Landscape & Climate | mean_CMD_sp | 2.790 | 0.788 | 3.54 | 0.002 | 0.326 | 26 |
| Standard.Length | Landscape & Climate | mean_CMI05 | -2.952 | 0.764 | -3.86 | \< 0.001 | 0.364 | 26 |
| Standard.Length | Landscape & Climate | mean_harvest_yr | 2.875 | 0.776 | 3.71 | 0.001 | 0.346 | 26 |
| Gape.Width | None | (No significant predictors found) | NA | NA | NA | \- | NA | (No Data) |
| Gill.Raker.Density | None | (No significant predictors found) | NA | NA | NA | \- | NA | (No Data) |

Table 1. Linear regression models evaluating environmental drivers.

``` r
write.csv(pub_table_formatted, "stickleback_regression_table.csv", row.names = FALSE)
```

### Robustness validation: Distance-Based Redundancy Analysis (dbRDA)

``` r
table_list <- list()
valid_models <- data.frame()
if (nrow(plot_df_final) > 0 && "variable" %in% colnames(plot_df_final)) {
  valid_models <- subset(plot_df_final, variable != "(No significant covariates)")
}

if(nrow(valid_models) > 0) {
  for(i in 1:nrow(valid_models)) {
    t_name   <- valid_models$trait[i]
    v_name   <- valid_models$variable[i]
    cat_name <- valid_models$category[i]
    is_fact  <- valid_models$is_factor[i]
    
    loop_df <- data.frame(
      y_val = sticks_rf_dat[[t_name]],
      x_val = rf_dat[[v_name]],
      stringsAsFactors = FALSE
    ) %>% 
      drop_na(y_val, x_val)
    
    if(!is_fact) {
      loop_df$x_scaled <- scale(as.numeric(loop_df$x_val))
    }
    
    # 1. Fit initial baseline model to calculate leverage metrics
    if(is_fact) {
      initial_fit <- lm(y_val ~ x_val, data = loop_df)
    } else {
      initial_fit <- lm(y_val ~ x_scaled, data = loop_df)
    }
    
    # 2. Identify individuals exceeding standard leverage threshold (4/N)
    cooks_d <- cooks.distance(initial_fit)
    n_obs <- nrow(loop_df)
    leverage_threshold <- 4 / n_obs
    influential_rows <- which(cooks_d > leverage_threshold)
    outliers_removed <- length(influential_rows)
    
    # 3. Clean dataset and fit final robustness model if anomalies exist
    if(outliers_removed > 0) {
      robust_df <- loop_df[-influential_rows, ]
    } else {
      robust_df <- loop_df
    }
    
    if(is_fact) {
      fit <- lm(y_val ~ x_val, data = robust_df)
      coefs <- as.data.frame(summary(fit)$coefficients)
      an_res <- anova(fit)
      
      b_val      <- coefs[2, "Estimate"]
      se_val     <- coefs[2, "Std. Error"]
      t_stat_val <- coefs[2, "t value"]
      p_val      <- coefs[2, "Pr(>|t|)"]
      r2_val     <- summary(fit)$r.squared
      df_val     <- as.character(an_res["Residuals", "Df"])
    } else {
      fit <- lm(y_val ~ x_scaled, data = robust_df)
      coefs <- as.data.frame(summary(fit)$coefficients)
      an_res <- anova(fit)
      
      b_val      <- coefs[2, "Estimate"]
      se_val     <- coefs[2, "Std. Error"]
      t_stat_val <- coefs[2, "t value"]
      p_val      <- coefs[2, "Pr(>|t|)"]
      r2_val     <- summary(fit)$r.squared
      df_val     <- as.character(an_res["Residuals", "Df"])
    }

    # 4. Build table summarizing results
    table_list[[length(table_list) + 1]] <- data.frame(
      Trait     = t_name,
      Category  = cat_name,
      Predictor = v_name,
      Beta      = if(is.na(b_val)) "-" else as.numeric(b_val),
      SE        = if(is.na(se_val)) "-" else as.numeric(se_val),
      t_stat    = if(is.na(t_stat_val)) "-" else as.numeric(t_stat_val),
      p_val     = as.numeric(p_val),
      R2        = if(is.na(r2_val)) "-" else as.numeric(r2_val), 
      DF        = df_val,
      Outliers_Dropped = outliers_removed,
      stringsAsFactors = FALSE
    )
  }
  pub_table_data <- do.call(rbind, table_list)
} else {
  pub_table_data <- data.frame(
    Trait = character(), Category = character(), Predictor = character(),
    Beta = numeric(), SE = numeric(), t_stat = numeric(), p_val = numeric(), 
    R2 = numeric(), DF = character(), Outliers_Dropped = numeric(), stringsAsFactors = FALSE
  )
}

# Flag and report traits that lack statistically significant drivers for summary table
missing_traits <- setdiff(names(sticks_morph), unique(pub_table_data$Trait))
if(length(missing_traits) > 0) {
  empty_rows <- data.frame(
    Trait     = missing_traits, 
    Category  = "None", 
    Predictor = "(No significant predictors found)",
    Beta = NA_real_, SE = NA_real_, t_stat = NA_real_, p_val = NA_real_, R2 = NA_real_, 
    DF = "(No Data)", Outliers_Dropped = 0,
    stringsAsFactors = FALSE
  )
  pub_table_data <- dplyr::bind_rows(pub_table_data, empty_rows)
}

pub_table_formatted <- pub_table_data %>%
  mutate(
    Beta   = case_when(Beta == "-" | is.na(Beta) ~ "-", TRUE ~ sprintf("%.3f", as.numeric(Beta))),
    SE     = case_when(SE == "-"   | is.na(SE)   ~ "-", TRUE ~ sprintf("%.3f", as.numeric(SE))),
    t_stat = case_when(t_stat == "-" | is.na(t_stat) ~ "-", TRUE ~ sprintf("%.2f", as.numeric(t_stat))),
    R2     = case_when(R2 == "-"   | is.na(R2)   ~ "-", TRUE ~ sprintf("%.3f", as.numeric(R2))),
    p_val  = case_when(
      is.na(p_val) ~ "-",
      as.numeric(p_val) < 0.001 ~ "< 0.001",
      TRUE                      ~ sprintf("%.3f", as.numeric(p_val))
    )
  )

knitr::kable(pub_table_formatted, caption = "Table 1. Linear regression models evaluating environmental drivers.")
```

| Trait | Category | Predictor | Beta | SE | t_stat | p_val | R2 | DF | Outliers_Dropped |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|---:|
| Armor.Plate.N | Water Column (Sonde) | SD_DO_Percent | -0.345 | 0.228 | -1.51 | 0.146 | 0.094 | 22 | 4 |
| Armor.Plate.N | Landscape & Climate | mean_d2_ogdc | 0.022 | 0.310 | 0.07 | 0.944 | 0.000 | 21 | 5 |
| Armor.Plate.N | Landscape & Climate | max_depth_est | 1.214 | 0.198 | 6.13 | \< 0.001 | 0.610 | 24 | 2 |
| Armor.Plate.N | Landscape & Climate | mean_spp1_pc | -0.162 | 0.244 | -0.66 | 0.513 | 0.020 | 22 | 4 |
| Armor.Plate.N | Landscape & Climate | mean_DD_0_05 | 0.506 | 3.336 | 0.15 | 0.881 | 0.001 | 23 | 3 |
| Armor.Plate.N | Landscape & Climate | perimeter | 1.451 | 0.167 | 8.67 | \< 0.001 | 0.750 | 25 | 1 |
| Body.Depth | Landscape & Climate | mean_sp3_cvr | 0.381 | 0.143 | 2.66 | 0.014 | 0.235 | 23 | 3 |
| Body.Depth | Landscape & Climate | mean_sa_proj_ht | 0.245 | 0.112 | 2.19 | 0.038 | 0.160 | 25 | 1 |
| CP.Depth | Water Column (Sonde) | SD_DO_Percent | 0.066 | 0.015 | 4.35 | \< 0.001 | 0.421 | 26 | 0 |
| CP.Depth | Landscape & Climate | mean_RH_sp | 0.055 | 0.019 | 2.91 | 0.008 | 0.253 | 25 | 1 |
| CP.Depth | Landscape & Climate | mean_d2_paved | -0.074 | 0.017 | -4.27 | \< 0.001 | 0.431 | 24 | 2 |
| CP.Depth | Water Column (Sonde) | Mean_fDOM_RFU | 0.058 | 0.019 | 3.04 | 0.006 | 0.278 | 24 | 2 |
| CP.Depth | Water Column (Sonde) | Mean_pH | -0.072 | 0.020 | -3.65 | 0.001 | 0.366 | 23 | 3 |
| Eye.Diameter | Landscape & Climate | mean_d2_lake | -0.046 | 0.018 | -2.61 | 0.015 | 0.208 | 26 | 0 |
| Eye.Diameter | Landscape & Climate | mean_CMI06 | -0.044 | 0.016 | -2.69 | 0.013 | 0.248 | 22 | 4 |
| Eye.Diameter | Landscape & Climate | mean_Tmax_sp | -0.033 | 0.017 | -1.99 | 0.058 | 0.137 | 25 | 1 |
| Eye.Diameter | Landscape & Climate | mean_Tave10 | -0.055 | 0.015 | -3.62 | 0.001 | 0.353 | 24 | 2 |
| Eye.Diameter | Landscape & Climate | mean_RH02 | 0.049 | 0.015 | 3.31 | 0.003 | 0.323 | 23 | 3 |
| Eye.Diameter | Landscape & Climate | mean_Eref12 | -0.057 | 0.015 | -3.78 | \< 0.001 | 0.373 | 24 | 2 |
| Eye.Diameter | Landscape & Climate | mean_Tmax10 | -0.059 | 0.015 | -4.05 | \< 0.001 | 0.406 | 24 | 2 |
| Eye.Diameter | Landscape & Climate | mean_DD_18_10 | 0.055 | 0.015 | 3.60 | 0.001 | 0.350 | 24 | 2 |
| Eye.Diameter | Landscape & Climate | mean_RH11 | 0.059 | 0.014 | 4.10 | \< 0.001 | 0.412 | 24 | 2 |
| Eye.Diameter | Landscape & Climate | mean_DD5_09 | -0.053 | 0.017 | -3.07 | 0.005 | 0.282 | 24 | 2 |
| Fin.Width | Landscape & Climate | mean_RH03 | 0.251 | 0.083 | 3.01 | 0.006 | 0.275 | 24 | 2 |
| Fin.Width | Landscape & Climate | mean_RH_at | 0.258 | 0.079 | 3.25 | 0.003 | 0.306 | 24 | 2 |
| Fin.Width | Landscape & Climate | mean_Tmax09 | -0.266 | 0.085 | -3.14 | 0.004 | 0.292 | 24 | 2 |
| Fin.Width | Landscape & Climate | mean_RH | 0.251 | 0.089 | 2.82 | 0.009 | 0.249 | 24 | 2 |
| Fin.Width | Landscape & Climate | mean_Eref10 | -0.183 | 0.080 | -2.30 | 0.030 | 0.175 | 25 | 1 |
| Fin.Width | Landscape & Climate | mean_Eref09 | -0.269 | 0.087 | -3.07 | 0.005 | 0.282 | 24 | 2 |
| Fin.Width | Landscape & Climate | mean_Eref_at | -0.208 | 0.079 | -2.62 | 0.015 | 0.216 | 25 | 1 |
| Gill.raker.1.length | Landscape & Climate | mean_vert_comp | -0.074 | 0.018 | -4.03 | \< 0.001 | 0.394 | 25 | 1 |
| Gill.raker.1.length | Landscape & Climate | mean_basal_area | 0.074 | 0.015 | 5.04 | \< 0.001 | 0.524 | 23 | 3 |
| Gill.raker.1.length | Landscape & Climate | mean_spp1_pc | -0.045 | 0.019 | -2.34 | 0.027 | 0.174 | 26 | 0 |
| Gill.raker.1.length | Landscape & Climate | mean_d2_ancient | 0.052 | 0.018 | 2.82 | 0.009 | 0.234 | 26 | 0 |
| Gill.raker.N | Landscape & Climate | mean_d2_lake | -0.533 | 0.193 | -2.76 | 0.011 | 0.233 | 25 | 1 |
| Gill.raker.N | Water Column (Sonde) | Mean_fDOM_RFU | -0.584 | 0.165 | -3.53 | 0.002 | 0.324 | 26 | 0 |
| Gill.raker.N | Landscape & Climate | mean_DD_0_11 | 0.315 | 0.185 | 1.70 | 0.101 | 0.104 | 25 | 1 |
| Gill.raker.N | Landscape & Climate | mean_DEM_vi | 0.471 | 0.187 | 2.52 | 0.019 | 0.209 | 24 | 2 |
| Head.Length | Landscape & Climate | mean_RH02 | 0.141 | 0.051 | 2.76 | 0.011 | 0.234 | 25 | 1 |
| Head.Length | Landscape & Climate | mean_DD5_11 | -0.114 | 0.053 | -2.14 | 0.042 | 0.155 | 25 | 1 |
| Head.Length | Landscape & Climate | mean_Tmax02 | -0.125 | 0.053 | -2.37 | 0.026 | 0.184 | 25 | 1 |
| Head.Length | Landscape & Climate | mean_RH01 | 0.166 | 0.051 | 3.25 | 0.003 | 0.306 | 24 | 2 |
| Head.Length | Landscape & Climate | mean_Eref10 | -0.187 | 0.054 | -3.49 | 0.002 | 0.347 | 23 | 3 |
| Head.Length | Landscape & Climate | mean_DD5_02 | -0.091 | 0.052 | -1.74 | 0.095 | 0.112 | 24 | 2 |
| Head.Length | Landscape & Climate | mean_crwn_clsre | 0.109 | 0.054 | 2.02 | 0.055 | 0.140 | 25 | 1 |
| Opercular.4bar.KT | Landscape & Climate | mean_CMD08 | 0.175 | 0.053 | 3.34 | 0.003 | 0.308 | 25 | 1 |
| Opercular.4bar.KT | Landscape & Climate | mean_CMI08 | -0.174 | 0.053 | -3.29 | 0.003 | 0.303 | 25 | 1 |
| Opercular.4bar.KT | Landscape & Climate | mean_Tmax12 | -0.134 | 0.052 | -2.59 | 0.016 | 0.218 | 24 | 2 |
| Opercular.4bar.KT | Landscape & Climate | mean_DD18_sp | 0.200 | 0.070 | 2.87 | 0.008 | 0.248 | 25 | 1 |
| Opercular.4bar.KT | Landscape & Climate | mean_RH03 | 0.183 | 0.059 | 3.10 | 0.005 | 0.285 | 24 | 2 |
| Opercular.4bar.KT | Landscape & Climate | mean_d2_lake | -0.185 | 0.058 | -3.16 | 0.004 | 0.285 | 25 | 1 |
| Opercular.4bar.KT | Landscape & Climate | mean_CMI10 | -0.205 | 0.059 | -3.45 | 0.002 | 0.323 | 25 | 1 |
| Opercular.4bar.KT | Landscape & Climate | mean_CMI07 | -0.169 | 0.054 | -3.12 | 0.005 | 0.280 | 25 | 1 |
| Opercular.4bar.KT | Landscape & Climate | mean_Tmax_wt | -0.135 | 0.052 | -2.61 | 0.015 | 0.221 | 24 | 2 |
| Snout.Length | Landscape & Climate | mean_crwn_clsre | 0.081 | 0.039 | 2.06 | 0.050 | 0.145 | 25 | 1 |
| Snout.Length | Landscape & Climate | mean_basal_area | 0.090 | 0.039 | 2.33 | 0.028 | 0.178 | 25 | 1 |
| Standard.Length | Landscape & Climate | mean_CMD05 | 2.402 | 0.708 | 3.39 | 0.002 | 0.315 | 25 | 1 |
| Standard.Length | Landscape & Climate | mean_PPT09 | -2.999 | 0.722 | -4.15 | \< 0.001 | 0.418 | 24 | 2 |
| Standard.Length | Landscape & Climate | mean_CMI09 | -2.667 | 0.679 | -3.93 | \< 0.001 | 0.382 | 25 | 1 |
| Standard.Length | Landscape & Climate | mean_DD_18_sp | 2.674 | 0.756 | 3.54 | 0.002 | 0.343 | 24 | 2 |
| Standard.Length | Landscape & Climate | mean_PPT05 | -3.472 | 0.849 | -4.09 | \< 0.001 | 0.421 | 23 | 3 |
| Standard.Length | Landscape & Climate | mean_CMD_sp | 2.402 | 0.708 | 3.39 | 0.002 | 0.315 | 25 | 1 |
| Standard.Length | Landscape & Climate | mean_CMI05 | -2.897 | 0.808 | -3.59 | 0.001 | 0.349 | 24 | 2 |
| Standard.Length | Landscape & Climate | mean_harvest_yr | 2.594 | 0.670 | 3.87 | \< 0.001 | 0.375 | 25 | 1 |
| Gape.Width | None | (No significant predictors found) | \- | \- | \- | \- | \- | (No Data) | 0 |
| Gill.Raker.Density | None | (No significant predictors found) | \- | \- | \- | \- | \- | (No Data) | 0 |

Table 1. Linear regression models evaluating environmental drivers.

``` r
write.csv(pub_table_formatted, "Table 1. stickleback_regression_table.csv", row.names = FALSE)

# Run dbRDA to test zooplankton effects on global morphology
valid_rows <- which(complete.cases(sticks_morph) & complete.cases(rf_dat))

if(length(valid_rows) > 5) {
  morph_matrix_scaled <- scale(sticks_morph[valid_rows, ])
  rf_matrix_clean     <- rf_dat[valid_rows, ]
  
  morph_distance <- vegan::vegdist(morph_matrix_scaled, method = "euclidean")
  
  if("Zoop_Community_PC1" %in% colnames(rf_matrix_clean)) {
    dbrda_model <- vegan::dbrda(morph_distance ~ Zoop_Community_PC1 + Zoop_Community_PC2, data = rf_matrix_clean)
    
    r2_raw <- vegan::RsquareAdj(dbrda_model)
    r2_df <- data.frame(
      Metric = c("R-Squared (Raw)", "Adjusted R-Squared"),
      Value  = c(sprintf("%.4f", r2_raw$r.squared), sprintf("%.4f", r2_raw$adj.r.squared))
    )
    
    global_permutation_test <- vegan::anova.cca(dbrda_model, permutations = 999)
    marginal_axis_test      <- vegan::anova.cca(dbrda_model, by = "margin", permutations = 999)
    
    print(knitr::kable(r2_df, 
                       caption = "Table 2a. Morphological variance explained by the zooplankton community matrix (dbRDA R²).",
                       align = "lc"))
    cat("\n\n") 
    
    # Format and save Global Model Table
    global_df <- as.data.frame(global_permutation_test)
    global_formatted <- global_df %>%
      rownames_to_column(var = "Effect") %>%
      mutate(
        Variance = sprintf("%.4f", Variance),
        F        = ifelse(is.na(F), "-", sprintf("%.4f", F)),
        `Pr(>F)` = case_when(
          is.na(`Pr(>F)`) ~ "-",
          `Pr(>F)` < 0.001 ~ "< 0.001",
          TRUE             ~ sprintf("%.3f", `Pr(>F)`)
        )
      )
    
    print(knitr::kable(global_formatted, 
                       caption = "Table 2b. Global non-parametric permutation test for the overall dbRDA model. Results indicate that variation in the zooplankton community does not influence global variation in stickleback morphology.",
                       align = "lcccc"))
    write.csv(global_formatted, "Table 2b. dbrda_global_permutation_table.csv", row.names = FALSE)
    cat("\n\n")
    
    # Format and save Marginal Axis Table
    marginal_df <- as.data.frame(marginal_axis_test)
    marginal_formatted <- marginal_df %>%
      rownames_to_column(var = "Axis") %>%
      mutate(
        Variance = sprintf("%.4f", Variance),
        F        = ifelse(is.na(F), "-", sprintf("%.4f", F)),
        `Pr(>F)` = case_when(
          is.na(`Pr(>F)`) ~ "-",
          `Pr(>F)` < 0.001 ~ "< 0.001",
          TRUE             ~ sprintf("%.3f", `Pr(>F)`)
        )
      )
    
    print(knitr::kable(marginal_formatted, 
                       caption = "Table 2c. Marginal permutation tests evaluating individual zooplankton PCA community axes.",
                       align = "lcccc"))
    write.csv(marginal_formatted, "Table 2c. dbrda_marginal_axis_table.csv", row.names = FALSE)
    
  } else {
    cat("\n⚠️ WARNING: 'Zoop_Community_PC1' was not found in your predictor dataframe. Ensure the PCA step ran successfully.\n")
  }
} else {
  cat("\nInsufficient matching row indices available to run the multivariate validation matrix.\n")
}
```

| Metric             |  Value  |
|:-------------------|:-------:|
| R-Squared (Raw)    | 0.0665  |
| Adjusted R-Squared | -0.0082 |

Table 2a. Morphological variance explained by the zooplankton community
matrix (dbRDA R²).

| Effect   | Df  | Variance |   F    | Pr(\>F) |
|:---------|:---:|:--------:|:------:|:-------:|
| Model    |  2  |  0.8639  | 0.8898 |  0.592  |
| Residual | 25  | 12.1361  |   \-   |   \-    |

Table 2b. Global non-parametric permutation test for the overall dbRDA
model. Results indicate that variation in the zooplankton community does
not influence global variation in stickleback morphology.

| Axis               | Df  | Variance |   F    | Pr(\>F) |
|:-------------------|:---:|:--------:|:------:|:-------:|
| Zoop_Community_PC1 |  1  |  0.4353  | 0.8967 |  0.516  |
| Zoop_Community_PC2 |  1  |  0.4097  | 0.8439 |  0.578  |
| Residual           | 25  | 12.1361  |   \-   |   \-    |

Table 2c. Marginal permutation tests evaluating individual zooplankton
PCA community axes.

Ecosystem gradient plot

``` r
# Build plot evaluating baseline environmental landscape variations
sig_vars_to_plot <- plot_df_final %>%
  filter(variable != "(No significant covariates)") %>%
  select(variable, category) %>%
  distinct()

if(nrow(sig_vars_to_plot) > 0) {
  numeric_sig_vars <- rf_dat %>%
    select(all_of(intersect(sig_vars_to_plot$variable, colnames(rf_dat)))) %>%
    select(where(is.numeric)) %>%
    colnames()
    
  if(length(numeric_sig_vars) > 0) {
    variation_long_df <- rf_dat %>%
      select(all_of(numeric_sig_vars)) %>%
      pivot_longer(cols = everything(), names_to = "variable", values_to = "raw_value") %>%
      left_join(sig_vars_to_plot, by = "variable")

    p <- ggplot(variation_long_df, aes(x = "", y = raw_value, fill = category)) +
      geom_boxplot(alpha = 0.25, outlier.shape = NA, color = "gray30", width = 0.5, linewidth = 0.4) +
      geom_point(aes(color = category), 
                 position = position_jitter(width = 0.2, height = 0, seed = 42), 
                 alpha = 0.60, size = 1.2) +
      
      facet_wrap(~ variable, scales = "free_y", ncol = 7) +
      
      scale_fill_manual(values = c(
        "Landscape & Climate"     = "#2c3e50",  
        "Water Column (Sonde)"    = "#2980b9",  
        "Zooplankton Food Supply" = "#27ae60"   
      )) +
      scale_color_manual(values = c(
        "Landscape & Climate"     = "#2c3e50",  
        "Water Column (Sonde)"    = "#2980b9",  
        "Zooplankton Food Supply" = "#27ae60"
      )) +
      labs(
        title = "Ecosystem Gradient Profile: Variation Across Sample Sites",
        subtitle = "Distributions of significant phenotypic predictors; points mark individual sample observations",
        x = NULL,
        y = "Observed Environmental Metric Value (Raw Scales)",
        color = "Predictor Classification",
        fill = "Predictor Classification"
      ) +
      theme_minimal(base_size = 10) + 
      theme(
        panel.border       = element_rect(color = "gray92", fill = NA, linewidth = 0.4),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background   = element_rect(fill = "gray97", color = "gray92"),
        strip.text         = element_text(face = "bold", size = 6.5, color = "gray20"),
        axis.text.y        = element_text(size = 6.5, color = "black"), 
        axis.text.x        = element_blank(), 
        axis.ticks.x       = element_blank(),
        panel.spacing.x    = unit(0.5, "lines"),  
        panel.spacing.y    = unit(0.8, "lines"),  
        legend.position    = "bottom",
        legend.margin      = ggplot2::margin(t = 2), # FIXED HERE
        legend.background  = element_rect(fill = "white", color = "gray95"),
        plot.title         = element_text(face = "bold", size = 12),
        plot.subtitle      = element_text(color = "gray40", size = 9, margin = ggplot2::margin(b = 8)) # FIXED HERE
      ) +
      guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
    
    print(p)
    
    ggsave("figures/Figure 2. Ecosystem Gradients.pdf", 
           plot = p, 
           width = 11,      
           height = 8.5,    
           units = "in", 
           dpi = 300)
           
    } else {
    cat("\nOnly categorical significant covariates (e.g. factors) were found. Quantitative gradient plot omitted.\n")
  }
} else {
  cat("\nNo significant covariates found; ecosystem gradient plot omitted.\n")
}
```

![](Ecological_Drivers_of_Stickleback_Morphological_Variation_GitHUB_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
