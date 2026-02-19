# ==========================================================
# Heatmap (SBR vs WU: FC ≥ 1.5, FDR ≤ 0.2)
# ▶ Legend에는 실제 등장한 Class만 표시(Other 규칙 유지)
# Input : Heatmap_Prepared_Data.xlsx (log10 값 + 통계 포함)
# Output: PDF heatmap + 사용 목록 엑셀
# ==========================================================

# ---------- 0) 사용자 설정 ----------
infile  <- "Heatmap_Prepared_Data.xlsx"           # 입력 엑셀
pdf_out <- "Heatmap_SBR_vs_WU_FC1p5_FDR0p2_.pdf"   # 출력 PDF
xlsx_out<- "Heatmap_Input_List_FC1p5_FDR0p2_.xlsx" # 출력 엑셀

fdr_cut <- 0.20
fc_cut  <- log2(1.5)  # 2배면 log2(2)

# 자잘한 class는 Other로 묶기(원치 않으면 c()로 비워도 됨)
to_collapse <- c(
  "Quinolines and derivatives","Benzopyrans","Tetrapyrroles and derivatives",
  "Lactones","Benzothiazoles","Benzene and substituted derivatives"
)

# ---------- 1) 패키지 ----------
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)

dir.create(dirname(pdf_out),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(xlsx_out), recursive = TRUE, showWarnings = FALSE)

# ---------- 2) 데이터 로드 ----------
dat <- read_excel(infile)

# ---------- 3) 샘플 컬럼 자동 인식 ----------
cols <- colnames(dat)
sham_cols <- cols[grepl("^Sham($|\\.|\\.\\.\\.)", cols)]
sbr_cols  <- cols[grepl("^SBR($|\\.|\\.\\.\\.)", cols) & !grepl("^SBR-WUSTL0717", cols)]
wu_cols   <- cols[grepl("^SBR-WUSTL0717($|\\.|\\.\\.\\.)", cols)]
if (length(sham_cols)==0 | length(sbr_cols)==0 | length(wu_cols)==0) {
  stop("Sham/SBR/WU 컬럼 자동 인식 실패: colnames(dat)를 확인하세요.")
}
dat <- dat %>%
  mutate(across(all_of(c(sham_cols, sbr_cols, wu_cols)), as.numeric))

# ---------- 4) 필터: SBR vs WU (FDR ≤ 0.2 & |FC| ≥ 1.5) ----------
sig <- dat %>%
  filter(!is.na(FDR_SBR_vs_WU),
         FDR_SBR_vs_WU <= fdr_cut,
         abs(log2FC_SBR_vs_WU) >= fc_cut)
if (nrow(sig) == 0) stop("선택된 metabolite가 없습니다. 컷오프를 조정해 주세요.")

# ---------- 5) 그룹 평균(로그10 값) ----------
grp_means <- sig %>%
  mutate(
    Sham = rowMeans(across(all_of(sham_cols)), na.rm = TRUE),
    SBR  = rowMeans(across(all_of(sbr_cols)),  na.rm = TRUE),
    `SBR-WUSTL0717` = rowMeans(across(all_of(wu_cols)), na.rm = TRUE)
  )

# ---------- 6) 행 Z-score ----------
mat <- grp_means %>% select(Sham, SBR, `SBR-WUSTL0717`) %>% as.matrix()
rownames(mat) <- make.unique(paste0(grp_means$Name, " [", grp_means$Compound.ID, "]"))
keep <- apply(mat, 1, sd, na.rm = TRUE) > 0
mat  <- mat[keep, , drop = FALSE]
grp_means <- grp_means[keep, ]
z_mat <- t(scale(t(mat)))  # 행 기준 Z-score

# ---------- 7) Class 주석 (정규화 + Other 유지) ----------
# 공백/철자 정규화
class_raw <- grp_means$Class %>%
  as.character() %>%
  str_squish()
class_raw <- dplyr::recode(class_raw,
                           "Organic sulfuric acids and derivatives" = "Organic sulfonic acids and derivatives"
                           # 필요 시 여기에 더 추가
)
# Other로 접기 규칙 적용
class_collapsed <- ifelse(class_raw %in% to_collapse, "Other", class_raw)

# 사용되는 클래스만 factor 레벨로 고정 ⇒ legend에도 실제 사용된 것만 표시
row_ann <- data.frame(Class = factor(class_collapsed),
                      row.names = rownames(z_mat))

# ---------- 8) Trend 정렬 + Class 내 클러스터링 ----------
trend <- z_mat[, "SBR-WUSTL0717"] - z_mat[, "SBR"]
trend_order <- order(trend, decreasing = TRUE)

ordered_rows <- c()
ordered_classes <- unique(row_ann$Class[trend_order])
ordered_classes <- c(setdiff(ordered_classes, "Other"), "Other")  # Other는 맨 아래

for (cls in ordered_classes) {
  idx <- trend_order[row_ann$Class[trend_order] == cls]
  if (length(idx) > 1) {
    hc <- hclust(dist(z_mat[idx, , drop = FALSE]), method = "ward.D2")
    ordered_rows <- c(ordered_rows, idx[hc$order])
  } else {
    ordered_rows <- c(ordered_rows, idx)
  }
}
z_final <- z_mat[ordered_rows, , drop = FALSE]
row_ann_final <- row_ann[ordered_rows, , drop = FALSE]

# ---------- 9) 주석 색상: 실제 사용 클래스만 사용 ----------
base_colors <- c(
  "Sphingolipids" = "#00C08B",
  "Glycerophospholipids" = "#3CB44B",
  "Carboxylic acids and derivatives" = "#FDBF6F",
  "Organic sulfonic acids and derivatives" = "#B15928",
  "Organooxygen compounds" = "#B2DF8A",
  "Organonitrogen compounds" = "#1F78B4",
  "Steroids and steroid derivatives" = "#FF7F00",
  "Prenol lipids" = "#6A3D9A",
  "Fatty Acyls" = "#CAB2D6",
  "Phenanthrenes and derivatives" = "#f28e2b",
  "Phenol ethers" = "#ff9d00",
  "Other" = "#999999",
  "Unknown" = "#666666"
)
# 실제 사용된 레벨
cls_used <- levels(row_ann_final$Class)
# 팔레트에 없는 클래스는 자동 색 할당
missing <- setdiff(cls_used, names(base_colors))
if (length(missing) > 0) {
  auto_pal <- colorRampPalette(brewer.pal(12, "Paired"))(length(missing))
  names(auto_pal) <- missing
  base_colors <- c(base_colors, auto_pal)
}
# legend에는 실제 사용된 클래스만
ann_colors <- list(Class = base_colors[names(base_colors) %in% cls_used])

# ---------- 10) 라벨/테두리 정리 ----------
rownames(z_final) <- iconv(rownames(z_final), from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
rownames(row_ann_final) <- rownames(z_final)

# ---------- 11) Heatmap 저장 ----------
pdf(pdf_out, width = 9, height = 11)
pheatmap(
  z_final,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101),
  cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_row = row_ann_final,
  annotation_colors = ann_colors,
  show_rownames = FALSE, labels_row = NA,  # 행 라벨 숨김
  show_colnames = TRUE,  fontsize_col = 11,
  border_color = NA,                       # 회색 구분선 제거
  main = "Z-score Heatmap (SBR vs WU: FC≥1.5, FDR≤0.2)"
)
dev.off()

# ---------- 12) 사용된 목록 저장 ----------
out_tab <- grp_means %>%
  transmute(
    Compound.ID, Name,
    Class = class_raw,
    Class_Collapsed = class_collapsed,
    Sham, SBR, `SBR-WUSTL0717`,
    log2FC_SBR_vs_WU, FDR_SBR_vs_WU
  )
write.xlsx(out_tab, xlsx_out, asTable = TRUE)

message("Done!  ▶ ", normalizePath(pdf_out), " / ", normalizePath(xlsx_out))