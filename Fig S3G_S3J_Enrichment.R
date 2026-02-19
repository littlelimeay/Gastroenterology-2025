install.packages("ggplot2")
install.packages("dplyr")

library(ggplot2)
library(dplyr)


# PDF로 저장
ggsave("GSEA_results.pdf", plot = p, width = 8, height = 6)  # 저장할 파일 이름, 크기 설정

#0에 점선
ggplot(data, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = number_genes_enriched, color = FDR.q.val)) +
  scale_color_gradient(low = "red", high = "blue", name = "FDR.q.val", limits = c(0, 0.1)) +  # FDR.q.val 범위 설정
  scale_size_continuous(range = c(3, 10), name = "Number of Genes Enriched") +
  labs(x = "Normalized Enrichment Score (NES)", y = "Pathway", title = "GSEA Results") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  # y축 글씨 크기 증가
    axis.text.x = element_text(size = 8),  # x축 글씨 크기 증가
    plot.margin = margin(10, 10, 10, 10)  # 여백 설정 (상, 우, 하, 좌)
  ) +
  coord_cartesian(xlim = c(-4, 3)) +  # NES 값 범위 설정
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")  # Add vertical dashed line at x = 0

# X(-1)
ggplot(data, aes(x = -NES, y = reorder(pathway, NES))) +  # Flip NES values here
  geom_point(aes(size = number_genes_enriched, color = FDR.q.val)) +
  scale_color_gradient(low = "red", high = "blue", name = "FDR.q.val", limits = c(0, 0.1)) +  # FDR.q.val 범위 설정
  scale_size_continuous(range = c(3, 10), name = "Number of Genes Enriched") +
  labs(x = "Normalized Enrichment Score (NES)", y = "Pathway", title = "GSEA Results") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  # y축 글씨 크기 증가
    axis.text.x = element_text(size = 8),  # x축 글씨 크기 증가
    plot.margin = margin(10, 10, 10, 10)  # 여백 설정 (상, 우, 하, 좌)
  ) +
  coord_cartesian(xlim = c(0, 4)) +  # NES 값 범위 설정
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")  # Add vertical dashed line at x = 0