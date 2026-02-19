# 📦 라이브러리
library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(stringr)

# 1) 데이터 불러오기 (시트 하나라면 sheet 인자 생략 가능)
path <- "2_Volcano_Significant_with_Classes_FC1_5.xlsx"
dat <- read_excel(path)

# 2) 'Phospholipids' 컬럼 잡기 (이름이 정확히 있을 때 우선 사용, 없으면 D열(4번째) 사용)
col_name <- if ("Phospholipids" %in% names(dat)) "Phospholipids" else names(dat)[4]
vec <- dat[[col_name]]

# 3) 전처리: 공백 제거, 대문자화
vec_clean <- toupper(str_trim(as.character(vec)))

# 4) 대상 클래스만 필터 (정확 일치 버전)
targets <- c("PE","PC","LPC","LPE")
vec_keep <- vec_clean[vec_clean %in% targets]

# (부분 매칭이 필요하면 위 대신 이 줄을 쓰세요)
# vec_keep <- vec_clean[str_detect(vec_clean, "\\b(PE|PC|LPC|LPE)\\b")]

if (length(vec_keep) == 0) stop("Phospholipids(D) 열에서 PE/PC/LPC/LPE가 없습니다. 값/열 이름을 확인하세요.")

# 5) 빈도 및 퍼센트
pie_df <- data.frame(Class = vec_keep) |>
  count(Class, name = "n") |>
  mutate(Percentage = n / sum(n) * 100)

# 6) 레전드 순서(원하는 고정 순서) 및 시계방향(12시 시작)
legend_order <- c("PC","PE","LPC","LPE")
present <- legend_order[legend_order %in% pie_df$Class]
pie_df$Class <- factor(pie_df$Class, levels = present)

chart_df <- pie_df
chart_df$Class <- factor(chart_df$Class, levels = rev(levels(pie_df$Class)))

# 7) 색상 (항목 수에 맞춰 Set3)
cols <- RColorBrewer::brewer.pal(max(4, length(present)), "Set3")[seq_along(present)]
names(cols) <- present

# 8) 저장 및 그리기
pdf("Pie_Phospholipids.pdf", width = 7, height = 6)
ggplot(chart_df, aes(x = "", y = Percentage, fill = Class)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y", start = 0) +    # 12시 시작
  scale_fill_manual(values = cols, breaks = present, name = "Class") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%"),
                y = cumsum(Percentage) - Percentage/2),
            size = 4.2) +
  theme_void() +
  labs(title = "Phospholipids") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14))
dev.off()