# Simulation in virtual patients ----------------------------

library(deSolve)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)

setwd('C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25')
input_dir_path <- file.path(getwd(), "sglt2i_dataset")
results_dir_path <- file.path(getwd(), "results_r")


# PD 데이터 (Virtual patients에서)
input_data <- read_csv(file.path(input_dir_path, "KSCPTSPRWS25_SGLT2i_Virtual_Patients.csv"))
input_data$dUGEc <- input_data$dUGEc * 7 

mean_avg_pg = input_data %>% summarise(mean_PG_avg = mean(PG_avg, na.rm = TRUE)) %>% pull(mean_PG_avg)
mean_baseline_pg = input_data %>% summarise(mean_PG_base = mean(PG_base, na.rm = TRUE)) %>% pull(mean_PG_base)
mean_baseline_HbA1c = input_data %>% summarise(mean_HbA1c = mean(HbA1c_base, na.rm = TRUE)) %>% pull(mean_HbA1c)
print(mean_avg_pg)
print(mean_baseline_pg)
print(mean_baseline_HbA1c)

# 공통 파라미터
params_fixed <- list(
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
    
    list(c(dHbA1cdrug), 
         dHbA1c = HbA1cplacebo + HbA1cdrug - HbA1cbaseline,
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
    parms <- c(params_fixed, dUGEc = subj$dUGEc, FPGbaseline=subj$PG_base, HbA1cbaseline=subj$HbA1c_base)
    state <- c(HbA1cdrug = 0)
    
    out <- ode(y = state, times = times, func = model_fn, parms = parms) %>%
      as.data.frame() %>%
      mutate(ID = subj$ID,
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
    lower_CI = quantile(dHbA1c, 0.05, na.rm = TRUE),
    upper_CI = quantile(dHbA1c, 0.95, na.rm = TRUE),
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
  filename = file.path(results_dir_path, "[WSCT] Simulation (dHbA1c - virtual_groups).png"),
  plot = p,
  width = 10,
  height = 10,
  dpi = 300,
  bg = "white"
)
