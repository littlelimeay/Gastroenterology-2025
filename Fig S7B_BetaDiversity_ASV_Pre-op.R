# 1. 필요한 라이브러리
library(readxl)
library(ggplot2)

# 2. 데이터 불러오기
setwd("C:/Users/ayoungk/Desktop")
data <- read_excel("Beta.xlsx", sheet = "bray_curtis_pc_Analysis", skip = 8)

# 3. 열 이름 정리
colnames(data)[1:5] <- c("Group", "Subgroup", "SampleID", "PC1", "PC2")

# 4. 필요한 열 선택 & NA 제거
pca_data <- data[, c("Group", "Subgroup", "SampleID", "PC1", "PC2")]
pca_data <- na.omit(pca_data)

# 5. 새로운 그룹 변수 생성 (6개 그룹: Pre/Post + Subgroup)
pca_data$GroupLabel <- paste(pca_data$Group, pca_data$Subgroup, sep = "-")

# 6. 색상과 도형 정의
group_colors <- c(
  "Pre-op-Sham-Vehicle" = "#008000",
  "Pre-op-SBR-Vehicle" = "#A85400",
  "Pre-op-SBR-WUSTL0717" = "#0000C0",
  "Post-op-Sham-Vehicle" = "#008000",
  "Post-op-SBR-Vehicle" = "#A85400",
  "Post-op-SBR-WUSTL0717" = "#0000C0"
)

shape_map <- c("Pre-op" = 1, "Post-op" = 16)  # open vs filled

# 7. PCA plot
p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  # 점: 색은 GroupLabel (6개 그룹), 도형은 Group (Pre/Post)
  geom_point(aes(color = GroupLabel, shape = Group), size = 3) +
  # ellipse: 6개 그룹 기준으로 그리기
  stat_ellipse(aes(group = GroupLabel, color = GroupLabel),
               type = "norm",
               alpha = 0.3,
               linewidth = 0.5,
               show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = shape_map) +
  labs(
    title = "PCA Plot with Ellipses by 6 Groups (Pre/Post + Treatment)",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")

# 8. 저장
ggsave("PCA_Ellipse_6Groups.pdf", plot = p, width = 8, height = 6, units = "in")