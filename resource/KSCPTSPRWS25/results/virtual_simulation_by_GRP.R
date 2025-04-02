# library(deSolve)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)

setwd('C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/results')


# dUGEc 데이터
input_data <- read_csv("virtual_patients_for_simulation.csv")
input_data$dUGEc <- input_data$dUGEc * 7 * 0.7
# input_data <- input_data %>% filter(GRP == 3)
input_data$UID <- input_data$ID

mean_baseline_pg = input_data %>% summarise(mean_PG_ZERO = mean(PG_ZERO, na.rm = TRUE)) %>% pull(mean_PG_ZERO)
mean_baseline_HbA1c = input_data %>% summarise(mean_HbA1c = mean(HbA1c, na.rm = TRUE)) %>% pull(mean_HbA1c)
print(mean_baseline_pg)
print(mean_baseline_HbA1c)

# 공통 파라미터
params_fixed <- list(
  # FPGbaseline = mean_baseline_pg,  # 145
  # HbA1cbaseline = mean_baseline_HbA1c, # 6.72
  # FPGbaseline = 150,
  # HbA1cbaseline = 7.0,
  Pfmax = 1.9,
  Kfp = 0.34,
  DISfp = 0.0313,
  SLOPEfd = -43.3,
  Phmax = 0.051,
  Khp = 0.24,
  DIShp = 0.0031,
  Kin2 = 0.5,
  Kout = 0.2
)

# 시간 설정
times <- seq(0, 52, by = 1)

# 수정된 ODE 모델
model_fn <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    
    # 1. Placebo 모델
    FPGplacebo <- FPGbaseline + Pfmax * (1 - exp(-Kfp * t)) + DISfp * t
    FPG <- FPGplacebo + SLOPEfd * dUGEc
    HbA1cplacebo <- HbA1cbaseline - Phmax * (1 - exp(-Khp * t)) + DIShp * t
    Kin <- Kout * HbA1cbaseline - Kin2
    
    # 2. ODE: 약물 효과 (ODE 대상 변수는 HbA1cdrug)
    dHbA1cdrug <- (FPG / FPGbaseline) * Kin + Kin2 - Kout * (HbA1cplacebo + HbA1cdrug)
    
    
    # # 1. Placebo 모델
    # FPGplacebo <- PG_ZERO + Pfmax * (1 - exp(-Kfp * t)) + DISfp * t
    # FPG <- FPGplacebo + SLOPEfd * dUGEc
    # HbA1cplacebo <- HbA1c - Phmax * (1 - exp(-Khp * t)) + DIShp * t
    # Kin <- Kout * HbA1c - Kin2
    # 
    # # 2. ODE: 약물 효과 (ODE 대상 변수는 HbA1cdrug)
    # dHbA1cdrug <- (FPG / FPGbaseline) * Kin + Kin2 - Kout * (HbA1cplacebo + HbA1cdrug)
    
    list(c(dHbA1cdrug), 
         
         dHbA1c = HbA1cplacebo + HbA1cdrug - HbA1cbaseline,
         # Pfmax_term = Pfmax * (1 - exp(-Kfp * t)),
         # DISfp_term = DISfp * t,
         # dUGEc_term = SLOPEfd * dUGEc,
         # Phmax_term = Phmax * (1 - exp(-Khp * t)),
         # DIShp_term = DIShp * t,
         # HbA1cdrug = HbA1cdrug,
         HbA1ctotal = HbA1cplacebo + HbA1cdrug,
         FPGplacebo = FPGplacebo,
         FPG = FPG
    )
  })
}

# 시뮬레이션
results_list <- list()

grp_list <- unique(input_data$GRP)

for (grp in grp_list) {
  
# for (i in 1:nrow(input_data %>% filter(GRP==grp))) {
  for (i in 1:nrow(input_data)) {
  subj <- input_data[i, ]
  # parms <- c(params_fixed, dUGEc = subj$dUGEc)
  parms <- c(params_fixed, dUGEc = subj$dUGEc, FPGbaseline=subj$PG_ZERO, HbA1cbaseline=subj$HbA1c)
  state <- c(HbA1cdrug = 0)
  
  out <- ode(y = state, times = times, func = model_fn, parms = parms) %>%
    as.data.frame() %>%
    mutate(UID = subj$UID,
           GRP = subj$GRP)
  
  results_list[[i]] <- out
}

}
  
results_all <- bind_rows(results_list)
results_all



### 📊 그룹별 요약통계 (median + 5%, 95% percentile)
summary_df <- results_all %>%
  group_by(GRP, time) %>%
  summarise(
    median_dHbA1c = median(dHbA1c, na.rm = TRUE),
    lower_CI = quantile(dHbA1c, 0.25, na.rm = TRUE),
    upper_CI = quantile(dHbA1c, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

### 최종 시점에서의 그룹별 median 값 추출
last_time <- max(summary_df$time)
annotation_df <- summary_df %>%
  filter(time == last_time) %>%
  mutate(label = sprintf("%.2f", median_dHbA1c))

### 🎨 그룹별 CI 음영 + dashed line 시각화
p <- ggplot(summary_df, aes(x = time, y = median_dHbA1c, color = as.factor(GRP), fill = as.factor(GRP))) +
  geom_line(size = 1.2) +  # Median line
  # CI dashed line
  geom_line(aes(y = lower_CI), linetype = "dashed", size = 0.8, alpha = 0.8) +
  geom_line(aes(y = upper_CI), linetype = "dashed", size = 0.8, alpha = 0.8) +
  # CI 음영
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.2, color = NA) +
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "black", size = 0.8) +
  # 마지막 시점 median 값 annotation
  geom_text(
    data = annotation_df,
    aes(label = label, x = time + 3, y = median_dHbA1c),
    color = "black",
    inherit.aes = FALSE,
    size = 5
  ) +
  labs(
    title = "Simulated ΔHbA1c (%) over Time by Group",
    x = "Time (weeks)",
    y = "ΔHbA1c (%)",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal(base_size = 20)

# 저장
ggsave(
  filename = paste0("[WSCT] Simulation (dHbA1c - virtual_groups).png"),
  plot = p,
  width = 10,
  height = 10,
  dpi = 300
)

##########################################

# 
# 
# # 그룹별 반복
# for (grp in unique(summary_df$GRP)) {
# 
#   # 해당 그룹 데이터만 추출
#   grp_df <- summary_df %>% filter(GRP == grp)
# 
#   # 해당 그룹 최종 시점 mean값
#   annotation_df <- grp_df %>%
#     filter(time == max(time)) %>%
#     mutate(label = sprintf("%.2f", mean_dHbA1c))
# 
#   # 그림 그리기
#   p <- ggplot(grp_df, aes(x = time, y = mean_dHbA1c)) +
#     geom_line(color = "royalblue", size = 1.2) +
#     geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.2, fill = "royalblue", color = NA) +
#     geom_hline(yintercept = -0.5, linetype = "dashed", color = "black", size = 0.8) +
#     geom_text(
#       data = annotation_df,
#       aes(label = label, x = time + 3, y = mean_dHbA1c),
#       color = "black",
#       inherit.aes = FALSE,
#       size = 3
#     ) +
#     labs(
#       title = paste0("Simulated ΔHbA1c (%) - Group ", grp),
#       x = "Time (weeks)",
#       y = "ΔHbA1c (%)"
#     ) +
#     theme_minimal(base_size = 14)
# 
#   # 파일로 저장
#   ggsave(
#     filename = paste0("DeltaHbA1c_Group_", grp, ".png"),
#     plot = p,
#     width = 10,
#     height = 7,
#     dpi = 300
#   )
# }
# 
# 
