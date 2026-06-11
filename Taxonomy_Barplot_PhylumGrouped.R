
###############################################################################
# Taxonomy Composition Barplot — Phylum-grouped Hierarchical Color Palette
# Reference: StackOverflow hierarchical coloring approach
###############################################################################

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(scales)

# ── 0. Theme & setup ──
base_family <- ifelse(Sys.info()[["sysname"]] == "Windows", "Arial", "Helvetica")

Treat_colors <- c(
  "NC"  = "#1976D2", "PC"  = "#388E3C",
  "G25" = "#F57C00", "G50" = "#D32F2F"
)

theme_set(
  theme_bw(base_family = base_family) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      axis.title   = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
)

# ── 1. Load phyloseq (assume ps_rel exists) ──
# ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# ── 2. Aggregate to Genus level ──
ps_genus <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)
df <- psmelt(ps_genus)

df$Phylum <- ifelse(is.na(df$Phylum) | df$Phylum == "", "Unclassified", df$Phylum)
df$Genus  <- ifelse(is.na(df$Genus)  | df$Genus  == "", "Unclassified", df$Genus)

# ── 3. Select top N genera ──
TOP_N <- 15

top_genera <- df %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = TOP_N) %>%
  pull(Genus)

df$Genus_plot <- ifelse(df$Genus %in% top_genera, df$Genus, "Others")

# Keep phylum info for top genera
genus_phylum <- df %>%
  filter(Genus_plot != "Others") %>%
  distinct(Genus_plot, Phylum) %>%
  rename(Genus = Genus_plot)

# For "Others", assign phylum = "Others"
df$Phylum_plot <- ifelse(df$Genus_plot == "Others", "Others", df$Phylum)

# ── 4. Define phylum-based color palette ──
# Core approach: colorRampPalette within each phylum

phylum_base_colors <- c(
  "Firmicutes"        = "#3F51B5",   # Blue
  "Actinobacteriota"  = "#FF9800",   # Orange
  "Proteobacteria"    = "#4CAF50",   # Green
  "Bacteroidota"      = "#E91E63",   # Pink
  "Planctomycetota"   = "#9C27B0",   # Purple
  "Patescibacteria"   = "#00BCD4",   # Cyan
  "Chloroflexi"       = "#795548",   # Brown
  "Verrucomicrobiota" = "#607D8B",   # Grey-blue
  "Others"            = "#CCCCCC",
  "Unclassified"      = "#BDBDBD"
)

generate_shades <- function(base_hex, n) {
  if (n == 0) return(character(0))
  if (n == 1) return(base_hex)
  base_rgb <- col2rgb(base_hex) / 255
  base_hsl <- rgb2hsv(base_rgb[1], base_rgb[2], base_rgb[3])
  
  # Generate shades by varying lightness
  light_seq <- seq(0.3, 0.75, length.out = n)
  sapply(light_seq, function(l) {
    hsv(base_hsl[1], base_hsl[2], l)
  })
}

# Build color map
genera_by_phylum <- df %>%
  filter(Genus_plot != "Others") %>%
  group_by(Phylum_plot, Genus_plot) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(Phylum_plot, desc(total))

color_map <- c()
for (phylum in unique(genera_by_phylum$Phylum_plot)) {
  genera <- genera_by_phylum %>%
    filter(Phylum_plot == phylum) %>%
    pull(Genus_plot)
  
  base <- ifelse(phylum %in% names(phylum_base_colors),
                 phylum_base_colors[phylum], "#78909C")
  
  shades <- generate_shades(base, length(genera))
  names(shades) <- genera
  color_map <- c(color_map, shades)
}
color_map["Others"] <- "#CCCCCC"

# ── 5. Order factor levels ──
# Phylum order by abundance
phylum_order <- df %>%
  group_by(Phylum_plot) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(Phylum_plot)

# Move Others to end
phylum_order <- c(phylum_order[phylum_order != "Others"], "Others")

# Genus order within phylum
genus_order <- c()
for (p in phylum_order) {
  g <- genera_by_phylum %>%
    filter(Phylum_plot == p) %>%
    arrange(desc(total)) %>%
    pull(Genus_plot)
  genus_order <- c(genus_order, g)
}
genus_order <- c(genus_order, "Others")

df$Genus_plot <- factor(df$Genus_plot, levels = rev(genus_order))

# Sample order by treatment
df$Treat <- factor(df$Treat, levels = c("NC", "PC", "G25", "G50"))

# ── 6. Plot ──
p <- ggplot(df, aes(x = Sample, y = Abundance, fill = Genus_plot)) +
  geom_bar(stat = "identity", width = 0.85, color = "white", linewidth = 0.1) +
  facet_grid(~ Treat, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = color_map, name = "Genus") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Genus-level Composition (Phylum-grouped Colors)",
    x = "Sample",
    y = "Relative Abundance"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text  = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 8),
    panel.spacing = unit(0.5, "lines")
  ) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave("Fig3a_Genus_PhylumGrouped.png", p, width = 14, height = 6, dpi = 300)

# ── 7. Phylum-level plot ──
ps_phylum <- tax_glom(ps_rel, taxrank = "Phylum", NArm = FALSE)
df_phy <- psmelt(ps_phylum)
df_phy$Phylum <- ifelse(is.na(df_phy$Phylum) | df_phy$Phylum == "", "Unclassified", df_phy$Phylum)

# Top phyla
top_phyla <- df_phy %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 8) %>%
  pull(Phylum)

df_phy$Phylum_plot <- ifelse(df_phy$Phylum %in% top_phyla, df_phy$Phylum, "Others")
df_phy$Treat <- factor(df_phy$Treat, levels = c("NC", "PC", "G25", "G50"))

p_phy <- ggplot(df_phy, aes(x = Sample, y = Abundance, fill = Phylum_plot)) +
  geom_bar(stat = "identity", width = 0.85, color = "white", linewidth = 0.1) +
  facet_grid(~ Treat, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = phylum_base_colors, name = "Phylum") +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Phylum-level Composition", x = "Sample", y = "Relative Abundance") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text  = element_text(face = "bold", size = 12)
  )

ggsave("Fig3b_Phylum_Composition.png", p_phy, width = 12, height = 5, dpi = 300)

cat("[DONE] Taxonomy composition figures saved.\n")
