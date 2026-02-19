# 📌 1. 패키지 로드
library(tidyverse)
library(scales)
library(stringr)
library(RColorBrewer)

# 📌 2. 데이터 불러오기
setwd("C:/Users/ayoungk/Desktop")
df_long <- read.csv("Microbial_Composition__Long_Format_.csv")

# 📌 3. Taxonomy fallback
extract_taxa_fallback <- function(tax_string) {
  k <- str_extract(tax_string, "k__[^;]+") |> str_remove("k__")
  p <- str_extract(tax_string, "p__[^;]+") |> str_remove("p__")
  c <- str_extract(tax_string, "c__[^;]+") |> str_remove("c__")
  o <- str_extract(tax_string, "o__[^;]+") |> str_remove("o__")
  f <- str_extract(tax_string, "f__[^;]+") |> str_remove("f__")
  g <- str_extract(tax_string, "g__[^;]+") |> str_remove("g__")
  
  if (!is.na(g) && g != "NA" && g != "") {
    genus <- g
    family <- ifelse(!is.na(f) && f != "NA" && f != "", f, "UnknownFamily")
  } else if (!is.na(f) && f != "NA" && f != "") {
    genus <- "Unclassified"
    family <- f
  } else if (!is.na(o) && o != "NA" && o != "") {
    genus <- "Unclassified"
    family <- paste0("Order: ", o)
  } else if (!is.na(c) && c != "NA" && c != "") {
    genus <- "Unclassified"
    family <- paste0("Class: ", c)
  } else if (!is.na(p) && p != "NA" && p != "") {
    genus <- "Unclassified"
    family <- paste0("Phylum: ", p)
  } else {
    genus <- "Unclassified"
    family <- "Unknown"
  }
  
  paste0(genus, " (", family, ")")
}

df_long$Taxon_Label <- sapply(df_long$Taxonomy, extract_taxa_fallback)

# 📌 4. Top 20 Taxon 선정
top_taxa <- df_long %>%
  group_by(Taxon_Label) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 20) %>%
  pull(Taxon_Label)

# ✅ 1. Top 20 중 Unclassified와 분리
top_taxa_unclassified <- grep("^Unclassified", top_taxa, value = TRUE)
top_taxa_classified   <- setdiff(top_taxa, top_taxa_unclassified)

# ✅ 2. 최종 정렬 순서 구성
final_taxa_levels <- c(top_taxa_classified, top_taxa_unclassified, "Other")

# 📌 5. 나머지는 "Other"로 처리 + 레벨 정렬 적용
df_long$Taxon_Label_Collapsed <- ifelse(df_long$Taxon_Label %in% top_taxa,
                                        df_long$Taxon_Label,
                                        "Other")

df_long$Taxon_Label_Collapsed <- factor(df_long$Taxon_Label_Collapsed,
                                        levels = final_taxa_levels)

# 📌 6. Sample 순서
sample_order <- c(
  "6763.preop", "6764.preop", "6765.preop", "6766.preop", "6767.preop",
  "6768.preop", "6769.preop", "6770.preop", "6731.preop", "6757.preop",
  "6758.preop", "6759.preop", "6760.preop", "6737.preop", "6738.preop",
  "6740.preop", "6746.preop", "6748.preop", "6733.preop", "6734.preop",
  "6741.preop", "6742.preop", "6743.preop", "6745.preop", "6751.preop",
  "6752.preop", "6753.preop", "6754.preop", "6755.preop", "6763.MB",
  "6764.MB", "6765.MB", "6766.m", "6767.m", "6768.m", "6769.m",
  "6770.m", "6731.MB", "6757.m", "6758.m", "6759.m", "6760.m",
  "6737.MB", "6738.MB", "6740.MB", "6746.m", "6748.m", "6733.MB",
  "6734.MB", "6741.MB", "6742.MB", "6743.MB", "6745.MB", "6751.m",
  "6752.m", "6753.m", "6754.m", "6755.m"
)
df_long$Sample <- factor(df_long$Sample, levels = sample_order)

# 🎨 색상 생성
set3_colors <- brewer.pal(12, "Set3")
extra_colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                  "#E6AB02", "#A6761D", "#666666")  # 총 8개

top20_colors <- c(set3_colors, extra_colors)  # 총 20개
names(top20_colors) <- as.character(top_taxa)  # 이름 부여

# ➕ Other 색상 추가
top20_colors <- c(top20_colors, Other = "#D3D3D3")

# 📌 8. 그래프
p <- ggplot(df_long, aes(x = Sample, y = Abundance, fill = Taxon_Label_Collapsed)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_fill_manual(values = top20_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  ylab("Microbial Composition (Top 20 Taxa)") +
  xlab("Sample") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    legend.position = "right"
  ) +
  guides(fill = guide_legend(ncol = 2))

# 📌 9. 저장
ggsave("Microbial_Composition_Top20_Set3_Expanded_un.pdf", plot = p, width = 12, height = 6)
