library(deSolve)
library(readr)
library(dplyr)
library(ggplot2)

setwd('C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/results')

library(deSolve)
library(tidyverse)

# ΔUGEc 데이터
input_data <- read_csv("PD_data.csv")

# 공통 파라미터
params_fixed <- list(
  FPGbaseline = 160,
  HbA1cbaseline = 7.92,
  Pfmax = 1.9,
  Kfp = 0.34,
  DISfp = 3.13,
  # DISfp = 3.13,
  SLOPEfd = -43.3,
  Phmax = 0.051,
  Khp = 0.24,
  DIShp = 0.31,
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
         dHbA1cdrug = dHbA1cdrug,
         # HbA1cplacebo = HbA1cplacebo, #+ HbA1cdrug,
         HbA1c = HbA1cplacebo + HbA1cdrug,
         FPGplacebo = FPGplacebo,
         FPG = FPG,
         HbA1cplacebo = HbA1cplacebo)
  })
}

# 시뮬레이션
results_list <- list()

for (i in 1:nrow(input_data)) {
  subj <- input_data[i, ]
  parms <- c(params_fixed, dUGEc = subj$dUGEc)
  state <- c(HbA1cdrug = 0)
  
  out <- ode(y = state, times = times, func = model_fn, parms = parms) %>%
    as.data.frame() %>%
    mutate(UID = subj$UID)
  
  results_list[[i]] <- out
}

results_all <- bind_rows(results_list)
results_all
# 시각화
# ggplot(results_all, aes(x = time, y = HbA1c, color = UID)) +
#   geom_line(size = 1) +
#   labs(title = "Simulated HbA1c over Time",
#        x = "Time (weeks)", y = "HbA1c (%)") +
#   theme_minimal()
