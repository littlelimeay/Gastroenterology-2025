# 필요한 패키지
library(readxl)
library(dplyr)
library(stringr)
library(openxlsx)
library(rlang)

## 0) 파일 경로 (원하는 경로로 바꿔도 됩니다)
file_left  <- "25-0999-DataAnalyzed_Stat_w17.xlsx"                # Annotation 시트 있는 파일
file_right <- "Volcano_Significant_SumNorm_NoAutoscale_FC1_5_no17.xlsx" # 붙일 대상 파일
out_file   <- "Volcano_Significant_with_Classes_FC1_5.xlsx"             # 결과 저장 파일명

## 1) Annotation 시트 읽기 (왼쪽 엑셀의 4번째 시트 이름이 'Annotation'이라고 했으므로 이름으로 지정)
ann_raw <- read_excel(path = file_left, sheet = "Annotation")

## 2) 필요한 열만 선택 + 키/문자열 정리
clean_id <- function(x) x |> as.character() |> str_trim()

# 열 이름이 정확히 다음과 같다고 가정합니다: Compound.ID, Class, Sub.Class
# (혹시 대소문자/스페이스 다른 경우가 있다면 아래 select 부분에서 수정)
ann <- ann_raw %>%
  transmute(
    Compound.ID = clean_id(`Compound.ID`),
    Class       = as.character(`Class`),
    Sub.Class   = as.character(`Sub.Class`)
  ) %>%
  mutate(CompoundID_key = clean_id(Compound.ID)) %>%
  filter(!is.na(CompoundID_key) & CompoundID_key != "")

## Annotation 쪽 중복 키 체크(있으면 첫 번째만 사용)
dup_keys <- ann %>%
  count(CompoundID_key) %>%
  filter(n > 1)

if (nrow(dup_keys) > 0) {
  message("⚠️ Annotation에서 중복 키 발견: ", nrow(dup_keys), "개 (첫 행만 사용)")
  ann <- ann %>% distinct(CompoundID_key, .keep_all = TRUE)
}

## 3) 붙일 대상(오른쪽 엑셀) 읽기
df <- read_excel(path = file_right)

## 4) 키 컬럼 탐지: 'Compound_ID' 또는 'Compound.ID' 중 존재하는 것 사용
key_col <- if ("Compound_ID" %in% names(df)) {
  "Compound_ID"
} else if ("Compound.ID" %in% names(df)) {
  "Compound.ID"
} else {
  stop("오른쪽 파일에서 키 컬럼(Compound_ID 또는 Compound.ID)을 찾을 수 없습니다.")
}

## 5) 키 정리 후 left_join
df_to_join <- ann %>% select(CompoundID_key, Class, Sub.Class)

df_joined <- df %>%
  mutate(CompoundID_key = clean_id(.data[[key_col]])) %>%
  left_join(df_to_join, by = "CompoundID_key") %>%
  select(-CompoundID_key) %>%
  # Class/Sub.Class를 키 바로 뒤로 이동(원하면 위치 조정)
  relocate(Class, Sub.Class, .after = all_of(key_col))

## 6) 매칭 품질 점검
n_total <- nrow(df_joined)
n_matched <- sum(!is.na(df_joined$Class) | !is.na(df_joined$Sub.Class))
message("✅ 매칭된 행: ", n_matched, " / ", n_total)

if (n_matched < n_total) {
  not_matched_examples <- df_joined %>%
    filter(is.na(Class) & is.na(Sub.Class)) %>%
    head(10) %>%
    pull(!!sym(key_col))
  message("ℹ️ 매칭 안 된 키 예시(최대 10개): ",
          paste(not_matched_examples, collapse = ", "))
}

## 7) 엑셀로 저장
write.xlsx(df_joined, out_file)
message("💾 저장 완료: ", out_file)