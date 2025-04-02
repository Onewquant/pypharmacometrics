# library(deSolve)
library(readr)
library(dplyr)
library(ggplot2)

setwd('C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/results')

library(deSolve)
library(tidyverse)

# dUGEc 데이터
input_data <- read_csv("KSCPTSPR25_PD_Endpoint_Sim_data.csv")
input_data$dUGEc <- input_data$dUGEc * 7 * 0.8

mean_avg_pg = input_data %>% summarise(mean_PG_AVG = mean(PG_AVG, na.rm = TRUE)) %>% pull(mean_PG_AVG)
mean_baseline_pg = input_data %>% summarise(mean_PG_ZERO = mean(PG_ZERO, na.rm = TRUE)) %>% pull(mean_PG_ZERO)
mean_baseline_HbA1c = input_data %>% summarise(mean_HbA1c = mean(HbA1c, na.rm = TRUE)) %>% pull(mean_HbA1c)
print(mean_avg_pg)
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

results_all <- bind_rows(results_list)
results_all

# results_all
#시각화
p <- ggplot(results_all, aes(x = time, y = dHbA1c, color = as.factor(GRP), group = UID)) +
  geom_line(size = 1.0, alpha = 0.8) +
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "black", size = 0.8) +
  # annotate("text", x = 50, y = -0.45, label = "Target ΔHbA1c = -0.5", color = "black", size = 5, hjust = 1) +
  labs(title = "Simulated ΔHbA1c (%) over Time by Group",
       x = "Time (weeks)", 
       y = "ΔHbA1c (%)",
       color = "Group") +
  theme_minimal(base_size = 20)

# 파일로 저장
ggsave(
  filename = paste0("[WSCT] Simulation (dHbA1c - Subjects).png"),
  plot = p,
  width = 10,
  height = 10,
  dpi = 300
)

