# 필요한 패키지
library(tidyverse)

# 데이터 불러오기
df <- read.csv("LEfSe_filtered.csv")

# 필터링: LDA score ≥ 3.5만
df <- df %>% filter(LDA_Score >= 3.5)

# 그룹 이름 정제
df$Group <- recode(df$Group, "SBR.WUSTL0717" = "SBR-WUSTL0717")

# Taxon 축약 함수
shorten_taxon <- function(taxon) {
  parts <- unlist(strsplit(taxon, "__"))
  # 뒤에서부터 g__, f__ 등 찾기
  for (i in rev(seq_along(parts))) {
    part <- parts[i]
    if (!grepl("^k|^p|^c|^o", part) && part != "") {
      return(part)
    }
  }
  return(taxon)
}

df$Short_Taxa <- sapply(df$Taxon, shorten_taxon)

# Short_Taxa 정렬
df$Short_Taxa <- factor(df$Short_Taxa, levels = df$Short_Taxa[order(df$LDA_Score)])

# 색상 설정: 파랑 (#0000C0), 갈색 (#A45400)
custom_colors <- c("SBR-WUSTL0717" = "#0000C0", "SBR" = "#A45400")

# 그래프 그리기
p <- ggplot(df, aes(x = LDA_Score, y = Short_Taxa, fill = Group)) +
  geom_col(color = "black") +
  scale_fill_manual(values = custom_colors) +
  xlab("LDA score (log 10)") +
  ylab(NULL) +
  ggtitle("LEfSe LDA Score ≥ 3.5") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.title = element_blank(),
    panel.grid.major.y = element_blank()
  )

# PDF로 저장
ggsave("LEfSe_LDA_Barplot.pdf", plot = p, width = 8, height = 6)
