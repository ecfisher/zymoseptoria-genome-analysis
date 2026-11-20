# Visualise and analyse chromosome-level statistics
# for Zymoseptoria tritici IPO323, with comparisons
# between core and accessory chromosomes.

# ---- 1. Load packages ----
library(ggplot2)
library(dplyr)

# ---- 2. Import chromosome statistics ----
chromosome_stats <- read.csv("../data/processed/chromosome_statistics.csv")

# ---- 3. Extract chromosome metadata ----
headers <- readLines("../data/raw/IPO323/ncbi_dataset/data/GCA_000219625.1/GCA_000219625.1_MYCGR_v2.0_genomic.fna")
headers <- headers[grepl("^>", headers)]
accessions <- sub("^>(\\S+).*", "\\1", headers)

# ---- 4. Assign core and accessory groups ----
chromosome_number <- sub(".*chromosome ([0-9]+),.*", "\\1", headers)
chromosome_number <- as.integer(chromosome_number)
chromosome_info <- data.frame(chromosome_id = accessions, chromosome_number = chromosome_number)
# IPO323 chromosomes 1–13 are classified as core,
# while chromosomes 14–21 are classified as accessory.
chromosome_info$chromosome_group <- ifelse(
  chromosome_info$chromosome_number <= 13,
  "Core",
  "Accessory"
)
chromosome_stats <- merge(
  chromosome_stats,
  chromosome_info,
  by = "chromosome_id"
)

# ---- 5. Exploratory analysis ----
gc_stats <- chromosome_stats %>%
  group_by(chromosome_group) %>%
  summarise(
    n = n(), 
    mean_gc = mean(gc_content), 
    median_gc = median(gc_content), 
    sd_gc = sd(gc_content)
  )

length_stats <- chromosome_stats %>%
  group_by(chromosome_group) %>%
  summarise(
    n = n(), 
    mean_length = mean(length), 
    median_length = median(length), 
    sd_length = sd(length)
  )


# ---- 6. Visualisations ----
chromosome_stats$length_mbp <- chromosome_stats$length / 1000000

#Scatter plot showing the association of chromosome length and GC content
length_vs_content <- ggplot(chromosome_stats, aes(x = length_mbp, y = gc_content, color = chromosome_group))+
  geom_point(size = 3)+
  xlab("Log of Chromosome Length (Mbp)")+
  ylab("GC Content (%)")+
  theme_classic()+
  scale_x_log10()+
  labs(title = substitute(paste("Association of chromosome length and GC content in ", italic("Z. tritici"), " IPO323")), color = "Chromosome Group")+
  scale_color_manual(values = c("cornflowerblue", "darkorange"))

#Ranked bar plot showing chromosome lengths
chrom_length_ranked <- ggplot(chromosome_stats, aes(x=reorder(chromosome_number, length), y=length_mbp, fill = chromosome_group))+
  geom_col()+
  coord_flip()+
  theme_classic()+
  labs(x = "Chromosome Number", y = "Chromosome Length (Mbp)", title = substitute(paste("Chromosome lengths in ", italic("Z. tritici"), " IPO323")))+
  scale_fill_manual(values = c("cornflowerblue", "darkorange") )+
  labs(fill = "Chromosome Group")

#Boxplot showing GC content distribution for core and accessory chromosomes
core_vs_accessory_gc <- chromosome_stats %>% 
  ggplot(aes(x = chromosome_group, y = gc_content, fill = chromosome_group)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(width = 0.1)+
  labs(x = "Chromosome Group", y = "GC Content", title = substitute(paste("GC Content in ", italic("Z. tritici"), " IPO323")))+
  scale_fill_manual(values = c("cornflowerblue", "darkorange")) +
  theme_classic()+
  theme(legend.position = "none")

#Boxplot showing chromosome length distribution for core and accessory chromosomes
core_vs_accessory_length <- ggplot(chromosome_stats, aes(x=chromosome_group, y = length_mbp, fill = chromosome_group))+
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  theme_classic()+
  labs(x = "Chromosome Group", y = "Length (Mbp)", title = substitute(paste("Chromosome Length in ", italic("Z. tritici"), " IPO323")))+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("cornflowerblue", "darkorange")) 

# ---- 7. Statistical analysis ----
wilcox.test(gc_content ~ chromosome_group, data = chromosome_stats)

chromosome_stats %>%
  group_by(chromosome_group)%>%
  summarise(
    cor(length, gc_content, method = "pearson")
  )

# ---- 8. Export figures ----
ggsave(
  filename = "../results/figures/gc_content_by_length.png",
  plot = length_vs_content,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "../results/figures/ranked_chromosome_length.png",
  plot = chrom_length_ranked,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "../results/figures/gc_content_core_accessory.png",
  plot = core_vs_accessory_gc, 
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "../results/figures/length_core_accessory.png",
  plot = core_vs_accessory_length,
  width = 7,
  height = 5,
  dpi = 300
)
