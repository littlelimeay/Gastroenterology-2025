# 1. 필요한 라이브러리
library(readxl)
library(ggplot2)

# 2. 데이터 불러오기
setwd("C:/Users/Ayoung Kim/Desktop")  # 경로는 필요에 맞게 수정
data <- read_excel("Imputed_Labeled.xlsx", sheet = "Imputed")

# 3. Expression matrix 설정
expr <- data[, 4:ncol(data)]
rownames(expr) <- data$Compound.ID
sample_names <- colnames(expr)
group <- gsub("\\..*", "", sample_names)

# 4. IQR filtering (상위 40%)
iqr_vals <- apply(expr, 1, IQR)
threshold <- quantile(iqr_vals, 0.6)
filtered_expr <- expr[iqr_vals > threshold, ]
cat("Number of variables after IQR filtering:", nrow(filtered_expr), "\n")

# 5. Sum normalization
#expr_norm <- sweep(filtered_expr, 2, apply(filtered_expr, 2, median), FUN = "/") #Median
expr_norm <- sweep(filtered_expr, 2, colSums(filtered_expr), FUN = "/")

# 6. log10 transformation
expr_log <- log10(expr_norm + 1e-8)

# 7. Autoscaling
expr_scaled <- t(scale(t(expr_log)))

# 8. PCA 분석
pca_result <- prcomp(t(expr_scaled), center = FALSE, scale. = FALSE)
explained_var <- round(100 * summary(pca_result)$importance[2, 1:2], 2)

# 9. 데이터프레임 생성
pca_data <- data.frame(pca_result$x[, 1:2],
                       Sample = colnames(expr_scaled),
                       Group = group)

# 10. 그룹 색상 정의
group_colors <- c("Sham" = "#008000", "SBR" = "#A85400", "SBR-WUSTL0717" = "#0000C0")

# 11. PCA plot with transparent ellipse
p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Group), size = 3) +
  stat_ellipse(aes(group = Group, color = Group, fill = Group),
               type = "norm",
               alpha = 0.3,
               linewidth = 1,
               show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "PCA (MetaboAnalyst-style IQR→Norm→Log→Auto)",
    x = paste0("PC1 (", explained_var[1], "%)"),
    y = paste0("PC2 (", explained_var[2], "%)")
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())

# 12. 저장
ggsave("PCA_Sum_Style_Ellipse.pdf", plot = p, width = 8, height = 6, units = "in")
