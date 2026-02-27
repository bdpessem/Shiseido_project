# =============================================================================
# FIGURE 5: Sebum-stratified metabolomics analysis
# Components: Figure5a_lmm | PCA_PLSDA2 | final_figure_lmm
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(vegan)
library(mixOmics)
library(lme4)
library(lmerTest)
library(emmeans)
library(purrr)


# =============================================================================
# SHARED SETTINGS
# =============================================================================

colors_sebum <- c(
  "low"    = "#4575b4",
  "medium" = "#fee0b6",
  "high"   = "#e08214"
)

SEBUM_LOW_THRESH  <- 3.5   # µg/cm²
SEBUM_HIGH_THRESH <- 16.9  # µg/cm²

selected_classes <- c(
  "N-acyl-alpha amino acids",
  "N-acylethanolamines",
  "1-monoacylglycerols",
  "Fatty acid esters",
  "Long-chain fatty acids",
  "Fatty alcohols",
  "Alpha amino acids and derivatives",
  "Oligopeptides"
)


# =============================================================================
# PANEL A: Lipid class boxplots with LMM (Figure5a_lmm)
# =============================================================================

# --- Prepare canopus class mapping ---
canopus_mapping <- canopus %>%
  mutate(
    featureId      = as.character(featureId),
    specific_class = if_else(
      is.na(`ClassyFire#most specific class`) | `ClassyFire#most specific class` == "",
      "Unclassified",
      as.character(`ClassyFire#most specific class`)
    )
  ) %>%
  select(featureId, specific_class)

# --- Sum feature intensities per class per sample ---
# NOTE: uses data_sample_filtered (Forehead / Nose / Cheek only)
data_sample_class_sum <- data_sample_filtered %>%
  pivot_longer(-SampleID, names_to = "featureId", values_to = "intensity") %>%
  mutate(featureId = as.character(featureId)) %>%
  left_join(canopus_mapping, by = "featureId") %>%
  mutate(specific_class = if_else(is.na(specific_class), "Unclassified", specific_class)) %>%
  group_by(SampleID, specific_class) %>%
  summarise(class_intensity = sum(intensity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = specific_class, values_from = class_intensity, values_fill = 0)

# --- rCLR transform ---
data_sample_class_clr <- data_sample_class_sum %>%
  column_to_rownames("SampleID") %>%
  decostand(method = "rclr") %>%
  rownames_to_column("SampleID")

# --- Merge with metadata (sebum group per sample) ---
clr_long_subject <- data_sample_class_clr %>%
  left_join(
    metadata_with_runorder_clean %>% select(filename_2, Subject_ID, sebum_group),
    by = c("SampleID" = "filename_2")
  ) %>%
  rename(group = sebum_group) %>%
  mutate(group = factor(group, levels = c("low", "medium", "high")))

# --- Prepare long-format data for selected classes ---
clr_long_plot_df <- clr_long_subject %>%
  select(SampleID, Subject_ID, group, all_of(selected_classes)) %>%
  pivot_longer(cols = all_of(selected_classes), names_to = "Class", values_to = "rclr_intensity") %>%
  filter(!is.na(group), !is.na(rclr_intensity)) %>%
  mutate(
    group       = factor(group, levels = c("low", "medium", "high")),
    Class       = factor(Class, levels = selected_classes),
    class_index = as.numeric(Class),
    group_index = as.numeric(group),
    x           = (class_index - 1) * 4 + group_index
  )

# Class label x-positions for x-axis annotation
x_labels_df <- clr_long_plot_df %>%
  group_by(Class) %>%
  summarise(x_center = mean(x), .groups = "drop") %>%
  mutate(label = as.character(Class))

# Legend labels with group sizes
group_counts  <- clr_long_plot_df %>% distinct(SampleID, group) %>% count(group) %>% deframe()
group_labels  <- c(
  low    = paste0("Low Sebum (< ",    SEBUM_LOW_THRESH,  " µg/cm²)\n(n = ", group_counts["low"],    ")"),
  medium = paste0("Medium Sebum (",   SEBUM_LOW_THRESH,  "–", SEBUM_HIGH_THRESH, " µg/cm²)\n(n = ", group_counts["medium"], ")"),
  high   = paste0("High Sebum (> ",   SEBUM_HIGH_THRESH, " µg/cm²)\n(n = ", group_counts["high"],   ")")
)

# --- LMM pairwise comparisons (BH-adjusted) ---
lmm_pval_df <- map_dfr(selected_classes, function(class_name) {
  class_data <- clr_long_plot_df %>% filter(Class == class_name)
  model      <- lmer(rclr_intensity ~ group + (1 | Subject_ID), data = class_data)
  pairs(emmeans(model, ~ group), adjust = "BH") %>%
    as.data.frame() %>%
    mutate(
      Class    = class_name,
      Group1   = sub(" - .*", "", contrast),
      Group2   = sub(".*- ", "", contrast),
      Raw_p    = p.value
    ) %>%
    select(Class, Group1, Group2, Raw_p)
}) %>%
  mutate(
    Adj_p        = p.adjust(Raw_p, method = "BH"),
    Significance = case_when(
      Adj_p <= 0.001 ~ "***",
      Adj_p <= 0.01  ~ "**",
      Adj_p <= 0.05  ~ "*",
      TRUE           ~ "ns"
    ),
    Label       = paste0(Significance, " p = ", signif(Adj_p, 2), "\n"),
    class_index = match(Class, selected_classes),
    xmin        = (class_index - 1) * 4 + as.numeric(factor(Group1, levels = c("low", "medium", "high"))),
    xmax        = (class_index - 1) * 4 + as.numeric(factor(Group2, levels = c("low", "medium", "high"))),
    base_y      = mapply(function(cls, g1, g2) {
      max(clr_long_plot_df %>% filter(Class == cls, group %in% c(g1, g2)) %>% pull(rclr_intensity), na.rm = TRUE)
    }, Class, Group1, Group2)
  ) %>%
  group_by(Class) %>%
  arrange(base_y, .by_group = TRUE) %>%
  mutate(y_position = base_y + (row_number() - 1) * 1.2) %>%
  ungroup()

# --- Build panel A ---
Figure5a_lmm <- ggplot(clr_long_plot_df, aes(x = x, y = rclr_intensity, fill = group)) +
  geom_boxplot(aes(group = x), outlier.shape = NA, width = 0.6, color = "black") +
  geom_jitter(aes(color = group), width = 0.2, size = 0.5, alpha = 0.6) +
  geom_text(
    data        = x_labels_df,
    aes(x = x_center, y = min(clr_long_plot_df$rclr_intensity) - 0.6, label = label),
    inherit.aes = FALSE, size = 5, hjust = 0.5
  ) +
  geom_segment(
    data        = lmm_pval_df,
    aes(x = xmin, xend = xmax, y = y_position, yend = y_position),
    inherit.aes = FALSE
  ) +
  geom_text(
    data        = lmm_pval_df,
    aes(x = (xmin + xmax) / 2, y = y_position + 0.2, label = Label),
    inherit.aes = FALSE, size = 4
  ) +
  ylim(
    min(clr_long_plot_df$rclr_intensity) - 0.6,
    max(lmm_pval_df$y_position, na.rm = TRUE) + 1
  ) +
  scale_fill_manual(values  = colors_sebum, name = "Sebum Group", labels = group_labels) +
  scale_color_manual(values = colors_sebum, guide = "none") +
  labs(
    title = "Compositional abundance of lipid classes across sebum groups",
    y     = "rclr-transformed intensity",
    x     = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x        = element_blank(),
    axis.title         = element_text(size = 14),
    axis.text.y        = element_text(size = 14),
    plot.title         = element_text(size = 15, hjust = 0.5),
    panel.grid         = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_text(size = 12),
    legend.text        = element_text(size = 12),
    axis.line          = element_line(color = "black")
  )


# =============================================================================
# PANEL B: PCA of extreme sebum samples (PCA_plot_sebum_extreme)
# =============================================================================

# --- Assign extreme sebum group ---
metadata_with_runorder_clean <- metadata_with_runorder_clean %>%
  mutate(Extreme_Sebum = case_when(
    sebum < SEBUM_LOW_THRESH  ~ "Low",
    sebum > SEBUM_HIGH_THRESH ~ "High",
    TRUE                      ~ NA_character_
  ))

# --- Filter to extreme samples (Forehead / Nose / Cheek only) ---
extreme_filenames <- metadata_with_runorder_clean %>%
  filter(!is.na(Extreme_Sebum), Face_site %in% c("Forehead", "Nose", "Cheek")) %>%
  pull(filename_2)

data_sample_extreme <- data_sample %>%
  filter(SampleID %in% extreme_filenames)

# --- rCLR transform & PCA ---
data_sample_extreme_clr <- decostand(
  data_sample_extreme %>% column_to_rownames("SampleID"),
  method = "rclr"
)

PCA_extreme <- mixOmics::pca(data_sample_extreme_clr, ncomp = 2, center = TRUE, scale = TRUE)

PCA_extreme_scores <- data.frame(PCA_extreme$variates$X) %>%
  rownames_to_column("SampleID") %>%
  left_join(metadata_with_runorder_clean, by = c("SampleID" = "filename_2"))

# --- PERMANOVA ---
dist_extreme     <- vegdist(data_sample_extreme_clr, method = "euclidean")
permanova_extreme <- adonis2(
  dist_extreme ~ PCA_extreme_scores$Subject_ID + PCA_extreme_scores$Extreme_Sebum,
  PCA_extreme_scores, na.action = na.omit, by = "margin"
)

sebum_row       <- permanova_extreme[2, ]
permanova_label <- paste0(
  "PERMANOVA\nR² = ", round(sebum_row$R2 * 100, 1), "%",
  "\nF = ",           round(sebum_row$F,         2),
  "\np = ",           format(sebum_row$`Pr(>F)`, digits = 3, scientific = TRUE)
)

# --- Group sample counts for legend ---
extreme_counts <- PCA_extreme_scores %>%
  filter(!is.na(Extreme_Sebum)) %>%
  count(Extreme_Sebum) %>%
  deframe()

low_n  <- extreme_counts["Low"]
high_n <- extreme_counts["High"]

# --- Build panel B ---
PCA_plot_sebum_extreme <- PCA_extreme_scores %>%
  ggscatter(
    x = "PC1", y = "PC2", color = "Extreme_Sebum", alpha = 0.6,
    title = "PCA of RCLR-Transformed Metabolite Features by Sebum Level",
    xlab  = paste0("PC1 (", round(PCA_extreme$prop_expl_var$X[1] * 100, 1), "%)"),
    ylab  = paste0("PC2 (", round(PCA_extreme$prop_expl_var$X[2] * 100, 1), "%)"),
    ggtheme = theme_classic()
  ) +
  geom_point(
    data = PCA_extreme_scores %>%
      group_by(Extreme_Sebum) %>%
      summarise(across(matches("PC"), mean), .groups = "drop"),
    aes(PC1, PC2, color = Extreme_Sebum),
    size = 4, shape = 8
  ) +
  scale_color_manual(
    name   = "Sebum Group",
    values = c("Low" = "#4575b4", "High" = "#e08214"),
    labels = c(
      Low  = paste0("Low (< ",  SEBUM_LOW_THRESH,  " µg/cm²) (n = ", low_n,  ")"),
      High = paste0("High (> ", SEBUM_HIGH_THRESH, " µg/cm²) (n = ", high_n, ")")
    )
  ) +
  theme(
    plot.title    = element_text(size = 15),
    axis.title    = element_text(size = 12),
    axis.text     = element_text(size = 12),
    legend.title  = element_text(size = 12),
    legend.text   = element_text(size = 11)
  ) +
  coord_fixed() +
  annotate("text",
    x = Inf, y = Inf, label = permanova_label,
    hjust = 1.1, vjust = 1.2, size = 3.8
  )

# Placeholder for panel C (e.g. PLS-DA or other)
PCA_PLSDA2 <- ggarrange(
  PCA_plot_sebum_extreme, NULL,
  widths       = c(1, 1.5),
  labels       = c("b", "c"),
  font.label   = list(size = 16, face = "bold"),
  common.legend = TRUE,
  legend        = "bottom"
)


# =============================================================================
# PANELS D–K: VIP metabolite boxplots with LMM (final_figure_lmm)
# =============================================================================

# --- rCLR on all samples (used for per-feature boxplots) ---
data_sample_clr <- decostand(
  data_sample %>% column_to_rownames("SampleID"),
  method = "rclr"
)

clr_long_subject <- data_sample_clr %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(
    metadata_with_runorder_clean %>% select(filename_2, Subject_ID, sebum_group),
    by = c("SampleID" = "filename_2")
  ) %>%
  rename(group = sebum_group)

# --- Function: per-feature LMM boxplots ---
make_grouped_vip_boxplots_lmm <- function(vip_df, clr_data, canopus_df) {

  clr_data$group <- factor(clr_data$group, levels = c("low", "medium", "high"))
  vip_df         <- vip_df %>% mutate(ID = as.character(ID))

  vip_df <- vip_df %>%
    left_join(canopus_df %>% mutate(featureId = as.character(featureId)),
              by = c("ID" = "featureId"))

  # Group sizes for x-axis labels
  group_sizes  <- clr_data %>%
    filter(!is.na(group)) %>%
    distinct(SampleID, group) %>%
    count(group) %>%
    deframe()

  x_axis_labels <- c(
    low    = paste0("Low Sebum\n(< ",    SEBUM_LOW_THRESH,  " µg/cm²)\n(n = ", group_sizes["low"],    ")"),
    medium = paste0("Medium Sebum\n(",   SEBUM_LOW_THRESH,  "–", SEBUM_HIGH_THRESH, " µg/cm²)\n(n = ", group_sizes["medium"], ")"),
    high   = paste0("High Sebum\n(> ",   SEBUM_HIGH_THRESH, " µg/cm²)\n(n = ", group_sizes["high"],   ")")
  )

  plot_list  <- list()
  result_rows <- list()

  for (i in seq_len(nrow(vip_df))) {

    feature_id   <- as.character(vip_df$ID[i])
    if (!feature_id %in% colnames(clr_data)) next

    compound_name <- vip_df$Compound_Name_2[i]
    mz_val        <- round(vip_df$mz[i], 2)
    rt_val        <- round(vip_df$RT[i], 2)
    superclass    <- vip_df$`ClassyFire#superclass`[i]
    title_text    <- if (!is.na(compound_name) && compound_name != "") compound_name else paste0("Feature_", feature_id)

    temp_data <- clr_data %>%
      select(group, Subject_ID, !!sym(feature_id)) %>%
      filter(!is.na(.data[[feature_id]])) %>%
      mutate(Subject_ID = as.factor(Subject_ID))

    if (nrow(temp_data) == 0 || length(unique(temp_data$group)) < 2) next

    model <- tryCatch(
      lmer(as.formula(paste0("`", feature_id, "` ~ group + (1 | Subject_ID)")), data = temp_data),
      error = function(e) NULL
    )
    if (is.null(model)) next

    pairwise <- tryCatch({
      pairs(emmeans(model, ~ group), adjust = "BH") %>%
        as.data.frame() %>%
        mutate(
          Group1 = sub(" - .*", "", contrast),
          Group2 = sub(".*- ", "", contrast),
          signif = case_when(
            p.value < 0.001 ~ "***",
            p.value < 0.01  ~ "**",
            p.value < 0.05  ~ "*",
            TRUE            ~ "ns"
          ),
          label = paste0(signif, " p = ", formatC(p.value, format = "e", digits = 2)),
          x     = match(Group1, levels(temp_data$group)),
          xend  = match(Group2, levels(temp_data$group)),
          xmid  = (x + xend) / 2
        )
    }, error = function(e) NULL)

    y_max <- max(temp_data[[feature_id]], na.rm = TRUE)

    p <- ggplot(temp_data, aes(x = group, y = .data[[feature_id]], fill = group, color = group)) +
      geom_boxplot(outlier.shape = NA, color = "black", width = 0.5) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 0.5) +
      scale_fill_manual(values  = colors_sebum) +
      scale_color_manual(values = colors_sebum) +
      scale_x_discrete(labels = x_axis_labels) +
      labs(
        y     = paste0("rclr ", title_text, "\n(m/z: ", mz_val, ", RT: ", rt_val, ")"),
        x     = NULL,
        title = title_text
      ) +
      theme_minimal() +
      theme(
        axis.text      = element_text(size = 12),
        axis.title     = element_text(size = 12),
        plot.title     = element_text(size = 13, hjust = 0.5),
        axis.line      = element_line(color = "black"),
        panel.grid     = element_blank(),
        legend.position = "none"
      ) +
      expand_limits(y = y_max * 1.2) +
      coord_cartesian(clip = "off")

    if (!is.null(pairwise)) {
      pairwise$y <- seq(from = y_max * 1.05, length.out = nrow(pairwise), by = y_max * 0.3)
      p <- p +
        geom_segment(data = pairwise,
                     aes(x = x, xend = xend, y = y, yend = y),
                     inherit.aes = FALSE, color = "black") +
        geom_text(data = pairwise,
                  aes(x = xmid, y = y + y_max * 0.15, label = label),
                  inherit.aes = FALSE, size = 4)

      pairwise$Metabolite_ID   <- feature_id
      pairwise$Metabolite_Name <- title_text
      pairwise$mz              <- mz_val
      pairwise$RT              <- rt_val
      pairwise$Superclass      <- superclass
      result_rows[[length(result_rows) + 1]] <- pairwise
    }

    plot_key            <- if (!is.na(superclass) && superclass != "") paste0(superclass, "_", feature_id) else feature_id
    plot_list[[plot_key]] <- p
  }

  list(plots = plot_list, comparisons = do.call(rbind, result_rows))
}

# --- Selected VIP features ---
selected_ids <- c(2314, 1744, 2069, 1883, 2409, 2509, 2287, 2160)

VIPs_selected <- VIPs_extreme_Load %>%
  filter(ID %in% selected_ids, !is.na(Compound_Name_2))

# --- Generate per-feature plots ---
vip_plots_lmm <- make_grouped_vip_boxplots_lmm(VIPs_selected, clr_long_subject, canopus)

# --- Arrange panels D–K (8 plots, 4 × 2) ---
final_figure_lmm <- ggarrange(
  vip_plots_lmm$plots[[1]], vip_plots_lmm$plots[[4]],
  vip_plots_lmm$plots[[3]], vip_plots_lmm$plots[[8]],
  vip_plots_lmm$plots[[2]], vip_plots_lmm$plots[[5]],
  vip_plots_lmm$plots[[6]], vip_plots_lmm$plots[[7]],
  ncol       = 4,
  nrow       = 2,
  labels     = c("d", "e", "f", "g", "h", "i", "j", "k"),
  font.label = list(size = 16, face = "bold")
)


# =============================================================================
# ASSEMBLE FIGURE 5 (top_grid_lmm)
# =============================================================================

top_grid_lmm <- ggarrange(
  Figure5a_lmm,
  PCA_PLSDA2,
  final_figure_lmm,
  ncol       = 1,
  nrow       = 3,
  labels     = c("a", "", ""),
  heights    = c(1.25, 1.25, 3.5),
  font.label = list(size = 16, face = "bold")
)


# =============================================================================
# SAVE
# =============================================================================

ggsave("../figures/main/Figure5a-k.png",
  plot   = top_grid_lmm,
  width  = 20, height = 16, dpi = 300
)

ggsave("../figures/main/Figure5a-k.svg",
  plot   = top_grid_lmm,
  width  = 20, height = 16
)
