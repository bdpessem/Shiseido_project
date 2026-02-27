# =============================================================================
# Heatmap: Microbial taxa vs Skin biophysical properties & amino acids / NMFs
# =============================================================================


# 0. Install / load packages

required_pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "purrr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# 1. File paths 

metadata_file       <- "../data/Shiseido_metadata_2.tsv"
shiseido_df_file    <- "../data/Shiseidodfc.txt"
alpha1_file         <- "../output/alpha_diversity/217014_alpha-diversity.tsv"
alpha2_file         <- "../output/alpha_diversity/217015_alpha-diversity.tsv"
rel_abund_file      <- "../data/amplicon/16S_V4_relative_abundance_RefHit.tsv"

# 2. Load & merge metadata

Shiseido_metadata_3 <- read.csv(metadata_file, header = TRUE, sep = "\t",
                                 stringsAsFactors = FALSE)
Shiseidodf          <- read.csv(shiseido_df_file, header = TRUE, sep = "\t",
                                 stringsAsFactors = FALSE)
alpha1 <- read_tsv(alpha1_file, col_types = cols())
alpha2 <- read_tsv(alpha2_file, col_types = cols())

alpha_combined <- full_join(alpha1, alpha2, by = "#SampleID")

Shiseido_metadata <- Shiseido_metadata_3 %>%
  left_join(alpha_combined, by = c("SampleID_2" = "#SampleID")) %>%
  merge(Shiseidodf[, c("X.SampleID", "Face_site")],
        by.x = "SampleID_2", by.y = "X.SampleID", all.x = TRUE) %>%
  mutate(Subject_ID = sub("-.*", "", microbiome_ID))

# 3. Load & prepare relative abundance table

df_rel_abundance <- read_tsv(rel_abund_file)
colnames(df_rel_abundance)[1] <- "SampleID"

# Strip run-order prefix added by QIIME export (e.g. "130729.")
df_rel_abundance$SampleID <- str_remove(df_rel_abundance$SampleID, "^130729\\.")

# Helper: keep the lowest defined taxonomic level, skip "g__uncultured"
extract_lowest_defined_level <- function(tax_string) {
  parts <- str_split(tax_string, ";\\s*")[[1]]
  parts <- parts[parts != ""]
  if (length(parts) == 0) return("Unclassified")
  last <- parts[length(parts)]
  if (str_detect(last, "g__uncultured") && length(parts) > 1) {
    return(parts[length(parts) - 1])
  }
  last
}

tax_cols      <- setdiff(colnames(df_rel_abundance), "SampleID")
cleaned_names <- map_chr(tax_cols, extract_lowest_defined_level)

# Sum columns that map to the same cleaned name
df_tax_collapsed <- df_rel_abundance[, "SampleID", drop = FALSE]
for (taxon in unique(cleaned_names)) {
  matched_cols <- tax_cols[cleaned_names == taxon]
  df_tax_collapsed[[taxon]] <- rowSums(df_rel_abundance[, matched_cols, drop = FALSE])
}

# Merge with metadata (keeps shannon_entropy, faith_pd, and all skin variables)
df_merged2 <- df_tax_collapsed %>%
  left_join(Shiseido_metadata, by = c("SampleID" = "SampleID_2"))

# 4. Variable definitions 

taxa_vars <- c(
  "g__Staphylococcus", "f__Neisseriaceae", "g__Lawsonella", "g__Corynebacterium",
  "g__Streptococcus", "g__Enhydrobacter", "g__Acinetobacter", "g__Xanthomonas",
  "g__Haemophilus", "g__Paracoccus", "g__Pseudomonas", "g__Prevotella",
  "g__Neisseria", "g__Veillonella", "g__Sphingomonas", "g__Lactobacillus",
  "g__Cutibacterium", "shannon_entropy", "faith_pd"
)

x_vars <- c(
  "water_content", "TEWL_transepidermal_water_loss_", "SC_protein", "sebum",
  "a", "L", "b", "skin_temparature", "age", "cutemeter_R2", "cutemeter_R7",
  "Asp", "Glu", "hyP", "Ser", "Asn", "Gly", "Gln", "Thr", "Ala", "Tau",
  "His", "Cit", "Pro", "Arg", "Val", "Tyr", "Met", "Cys2", "Ile", "Leu",
  "Orn", "Lys", "Phe", "Trp", "Total_AAs", "Urea", "Creatinine", "UCA",
  "LA", "PCA", "D.glucosamine", "Total_org._NMF"
)

x_axis_labels <- c(
  "g__Staphylococcus"  = "Staphylococcus",
  "f__Neisseriaceae"   = "Neisseriaceae",
  "g__Lawsonella"      = "Lawsonella",
  "g__Corynebacterium" = "Corynebacterium",
  "g__Streptococcus"   = "Streptococcus",
  "g__Enhydrobacter"   = "Enhydrobacter",
  "g__Acinetobacter"   = "Acinetobacter",
  "g__Xanthomonas"     = "Xanthomonas",
  "g__Haemophilus"     = "Haemophilus",
  "g__Paracoccus"      = "Paracoccus",
  "g__Pseudomonas"     = "Pseudomonas",
  "g__Prevotella"      = "Prevotella",
  "g__Neisseria"       = "Neisseria",
  "g__Veillonella"     = "Veillonella",
  "g__Sphingomonas"    = "Sphingomonas",
  "g__Lactobacillus"   = "Lactobacillus",
  "g__Cutibacterium"   = "Cutibacterium",
  "shannon_entropy"    = "Shannon Entropy",
  "faith_pd"           = "Faith PD"
)

y_axis_labels <- c(
  "L"                              = "Skin Brightness (L*)",
  "a"                              = "Skin Redness (a*)",
  "b"                              = "Skin Yellowness (b*)",
  "cutemeter_R2"                   = "Cutemeter R2 (Elasticity Index)",
  "cutemeter_R7"                   = "Cutemeter R7 (Elasticity Recovery Index)",
  "sebum"                          = "Sebum (\u00b5g/cm\u00b2)",
  "water_content"                  = "Water Content (a.u.)",
  "TEWL_transepidermal_water_loss_"= "TEWL (g/[m\u00b2\u00b7h])",
  "SC_protein"                     = "SC Protein (\u00b5g/cm\u00b2)",
  "skin_temparature"               = "Skin Temperature (\u00b0C)",
  "age"                            = "Age (years)",
  "Asp"                            = "Aspartic Acid",
  "Glu"                            = "Glutamic Acid",
  "hyP"                            = "Hydroxyproline",
  "Ser"                            = "Serine",
  "Asn"                            = "Asparagine",
  "Gly"                            = "Glycine",
  "Gln"                            = "Glutamine",
  "Thr"                            = "Threonine",
  "Ala"                            = "Alanine",
  "Tau"                            = "Taurine",
  "His"                            = "Histidine",
  "Cit"                            = "Citrulline",
  "Pro"                            = "Proline",
  "Arg"                            = "Arginine",
  "Val"                            = "Valine",
  "Tyr"                            = "Tyrosine",
  "Met"                            = "Methionine",
  "Cys2"                           = "Cystine",
  "Ile"                            = "Isoleucine",
  "Leu"                            = "Leucine",
  "Orn"                            = "Ornithine",
  "Lys"                            = "Lysine",
  "Phe"                            = "Phenylalanine",
  "Trp"                            = "Tryptophan",
  "Total_AAs"                      = "Total Amino Acids",
  "Urea"                           = "Urea",
  "Creatinine"                     = "Creatinine",
  "UCA"                            = "Urocanic Acid",
  "LA"                             = "Lactic Acid",
  "PCA"                            = "Pyrrolidone Carboxylic Acid",
  "D.glucosamine"                  = "D-Glucosamine",
  "Total_org._NMF"                 = "Total Organic NMF"
)

# 5. Compute correlations (Pearson or Spearman based on normality) 

results <- list()

for (x in taxa_vars) {
  for (y in x_vars) {

    # Skip if either column is absent in the merged data
    if (!x %in% colnames(df_merged2) || !y %in% colnames(df_merged2)) next

    sub_df <- df_merged2 %>%
      select(all_of(c(x, y))) %>%
      filter(!is.na(.data[[x]]) & !is.na(.data[[y]]))

    if (nrow(sub_df) < 3) next

    p_x <- shapiro.test(sub_df[[x]])$p.value
    p_y <- shapiro.test(sub_df[[y]])$p.value

    method <- ifelse(p_x > 0.05 & p_y > 0.05, "pearson", "spearman")
    result <- cor.test(sub_df[[x]], sub_df[[y]], method = method)

    results[[paste(x, y, sep = "_")]] <- data.frame(
      X           = x,
      Y           = y,
      Correlation = result$estimate,
      p_value     = result$p.value,
      method      = method
    )
  }
}

cor_df <- bind_rows(results)

# 6. Add significance stars and axis labels 

cor_df <- cor_df %>%
  mutate(
    stars = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    X_label = x_axis_labels[X],
    Y_label = y_axis_labels[Y]
  )

# Preserve desired axis order
cor_df$X_label <- factor(cor_df$X_label, levels = x_axis_labels[taxa_vars])
cor_df$Y_label <- factor(cor_df$Y_label, levels = rev(y_axis_labels[x_vars]))

# 7. Plot heatmap 

heatmap_genera <- ggplot(cor_df, aes(x = X_label, y = Y_label, fill = Correlation)) +
  geom_tile(color = "white", width = 0.9, height = 0.9) +
  geom_text(aes(label = stars), color = "black", size = 5,
            hjust = 0.5, vjust = 0.8) +
  scale_fill_gradient2(
    low      = "#4575b4",
    mid      = "white",
    high     = "#e08214",
    midpoint = 0,
    limits   = c(-0.5, 0.5),
    name     = "Spearman \u03c1"
  ) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.text        = element_text(size = 12),
    axis.ticks       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_blank(),
    panel.grid       = element_blank()
  )

print(heatmap_genera)

# 8. Save outputs

ggsave("../figures/Figure6a.png", plot = heatmap_genera, width = 16, height = 9, dpi = 300)
ggsave("../figures/Figure6a.svg", plot = heatmap_genera, width = 16, height = 9)

message("Done! Outputs saved as Figure6a.png and Figure6a.svg")
