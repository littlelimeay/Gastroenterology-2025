# 필요한 라이브러리
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(openxlsx)

# 데이터 불러오기
setwd("C:/Users/Ayoung Kim/Desktop")
data <- read_excel("Imputed_Labeled.xlsx", sheet = "Imputed")

# 그룹 설정
expr <- data[, 4:ncol(data)]
rownames(expr) <- data$Compound.ID
sample_names <- colnames(expr)
group <- gsub("\\..*", "", sample_names)
sbr_cols <- group == "SBR"
wustl_cols <- group == "SBR-WUSTL0717"

# ✅ 1.Sum normalization
expr_norm <- sweep(expr, 2, apply(expr, 2, sum), FUN = "/")

# ✅ 2. Log2 transformation (MetaboAnalyst와 동일)
expr_log2 <- log2(expr_norm + 1e-8)

# ✅ 3. Volcano 분석용 FC, p-value 계산 (autoscaling ❌)
log2_fc <- rowMeans(expr_log2[, wustl_cols]) - rowMeans(expr_log2[, sbr_cols])
p_values <- apply(expr_log2, 1, function(x) t.test(x[sbr_cols], x[wustl_cols])$p.value)
fdr_values <- p.adjust(p_values, method = "BH")
neg_log10_p <- -log10(p_values)

# Volcano 데이터프레임 생성
volcano_df <- data.frame(
  Compound_ID = rownames(expr_log2),
  Name = data$Name[match(rownames(expr_log2), data$Compound.ID)],
  log2FC = log2_fc,
  p_value = p_values,
  FDR = fdr_values,
  neg_log10_p = neg_log10_p
)

# 유의성 구분
volcano_df$Significance <- "Not significant"
volcano_df$Significance[volcano_df$FDR < 0.1 & volcano_df$log2FC > 0.585] <- "Upregulated"
volcano_df$Significance[volcano_df$FDR < 0.1 & volcano_df$log2FC < -0.585] <- "Downregulated"

# 상위 10개 라벨링
top_up <- volcano_df %>% filter(Significance == "Upregulated") %>% arrange(desc(log2FC)) %>% slice(1:10)
top_down <- volcano_df %>% filter(Significance == "Downregulated") %>% arrange(log2FC) %>% slice(1:10)
top_labeled <- bind_rows(top_up, top_down)

# 색상 설정
custom_colors <- c("Upregulated" = "#B2182B", "Downregulated" = "#2166AC", "Not significant" = "gray")

# Volcano plot 그리기
p <- ggplot(volcano_df, aes(x = log2FC, y = neg_log10_p, color = Significance)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = custom_colors) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dotted", color = "gray") +
  geom_text_repel(data = top_labeled,
                  aes(label = Name),
                  size = 2.5,
                  box.padding = 0.3,
                  max.overlaps = Inf) +
  coord_cartesian(ylim = c(0, 10)) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: SumNorm + Log2 (no autoscale)",
    x = "log2 Fold Change",
    y = "-log10(p-value)"
  ) +
  theme(legend.title = element_blank())

# 저장
ggsave("Volcano_SBR_WUSTL0717_SumNormNoAutoscale_FC1_5_no17.pdf", plot = p, width = 8, height = 5)

# 유의적인 항목만 저장
significant_metabs <- volcano_df %>%
  filter(Significance %in% c("Upregulated", "Downregulated"))

write.xlsx(significant_metabs, file = "Volcano_Significant_SumNorm_NoAutoscale_FC1_5_no17.xlsx", rowNames = FALSE)
