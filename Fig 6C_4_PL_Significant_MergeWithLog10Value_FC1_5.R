library(readxl)
library(dplyr)
library(stringr)
library(writexl) # 결과를 다시 엑셀로 저장할 때

# 1️⃣ 파일 불러오기
# 첫 번째 파일 (Phospholipids 정보)
file1 <- read_excel("2_Volcano_Significant_with_Classes_FC1_5.xlsx")

# 두 번째 파일 (Log 값 포함)
file2 <- read_excel("25-0999-DataAnalyzed_Stat.xlsx")

# 2️⃣ 첫 번째 파일에서 Phospholipids D열에서 대상 행 필터
# 열 이름 자동 탐색
col_name <- if ("Phospholipids" %in% names(file1)) "Phospholipids" else names(file1)[4]

# 대상 클래스만 선택
targets <- c("PE", "PC", "LPC", "LPE")
file1_filtered <- file1 %>%
  filter(.data[[col_name]] %in% targets)

# 3️⃣ 두 파일의 공통 키 확인
# 열 이름 출력
names(file1_filtered)
names(file2)

# 공통으로 존재하는 키 확인 (예: Compound.ID 또는 Metabolite Name)
# 예를 들어 공통 키가 "Compound.ID" 라고 가정
common_key <- "Compound.ID"  # 실제 이름으로 바꿔야 함

# 4️⃣ 두 파일 병합 (왼쪽 파일 기준으로 Log 값 붙이기)
merged_data <- file1_filtered %>%
  left_join(file2 %>% select(all_of(common_key), Log), by = common_key)

# 5️⃣ 결과 확인
head(merged_data)

# 6️⃣ 결과 저장
write_xlsx(merged_data, "Merged_Phospholipids_with_Log.xlsx")