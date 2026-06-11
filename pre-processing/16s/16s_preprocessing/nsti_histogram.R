library(ggplot2)
library(data.table)

# read NSTI table
nsti <- fread("picrust2_out/marker_predicted_and_nsti.tsv.gz")

# rename columns for safety
colnames(nsti)[2] <- "NSTI"

# histogram
p <- ggplot(nsti, aes(x=NSTI)) +
  geom_histogram(bins=30, fill="#2E86AB", color="black") +
  theme_bw() +
  labs(
    title="NSTI distribution (PICRUSt2 prediction reliability)",
    x="NSTI score",
    y="ASV count"
  )

# save plot
ggsave("NSTI_histogram.png", p, width=6, height=4, dpi=300)

# print mean NSTI
cat("Mean NSTI:", mean(nsti$NSTI), "\n")
