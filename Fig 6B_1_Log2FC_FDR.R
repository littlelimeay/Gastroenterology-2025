# 6. Log 데이터에 결과 추가 -------------------------------------------
final_data <- log_data %>%
  mutate(
    log2FC_Sh_vs_SBR = FC_Sh_vs_SBR,
    log2FC_SBR_vs_WU = FC_SBR_vs_WU
  ) %>%
  bind_cols(stats_sham_sbr, stats_sbr_wu)

# 7. Annotation 정보 매칭 --------------------------------------------
final_merged <- final_data %>%
  left_join(anno_data %>% select(Compound.ID, Class, Sub.Class),
            by = "Compound.ID")

# 8. Excel 파일로 Export --------------------------------------------
write.xlsx(final_merged, "Heatmap_Prepared_Data.xlsx")

cat("완료! 'Heatmap_Prepared_Data.xlsx' 파일이 생성되었습니다.\n")