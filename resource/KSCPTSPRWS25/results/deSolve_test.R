# install.packages("deSolve")
library(deSolve)
library(readr)
library(dplyr)

setwd('C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/results')

# 대상자별 ΔUGEc 데이터 불러오기
input_data <- read_csv("PD_data.csv")  # 파일 경로는 실제 환경에 맞게 조정

# 공통 파라미터 설정 (예시값, 필요 시 조정)
common_params <- list(
  FPGbaseline = 120,
  HbA1cbaseline = 7.92,
  Pfmax = 1.90,
  Kfp = 0.34,
  DISfp = 3.13,
  SLOPEfd = -43.3,
  Phmax = 0.051,
  Khp = 0.24,
  DIShp = 0.31,
  Kin2 = 0.5,
  Kout = 0.2
)

# 시뮬레이션 시간 설정
times <- seq(0, 52, by = 1)

# ODE 모델 함수 정의
model_fn <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    FPGplacebo <- FPGbaseline + Pfmax * (1 - exp(-Kfp * t)) + DISfp * t
    FPG <- FPGplacebo + SLOPEfd * dUGEc
    HbA1cplacebo <- HbA1cbaseline - Phmax * (1 - exp(-Khp * t)) + DIShp * t
    Kin <- Kout * HbA1cbaseline - Kin2
    HbA1c <- HbA1cplacebo + HbA1cdrug
    dHbA1cdrug <- (FPG / FPGbaseline) * Kin + Kin2 - Kout * HbA1c
    
    list(c(dHbA1cdrug = dHbA1cdrug),
         FPG = FPG,
         FPGplacebo = FPGplacebo,
         HbA1cplacebo = HbA1cplacebo,
         HbA1c = HbA1cplacebo + HbA1cdrug)
  })
}

# 대상자별 시뮬레이션
results <- list()

for (i in 1:nrow(input_data)) {
  subj <- input_data[i, ]
  
  parms <- c(common_params, dUGEc = subj$dUGEc)
  state <- c(HbA1cdrug = 0)
  
  out <- ode(y = state, times = times, func = model_fn, parms = parms)
  out_df <- as.data.frame(out)
  out_df$UID <- subj$UID
  
  results[[i]] <- out_df
}

# 결과 통합
all_results <- bind_rows(results)

# 결과 확인
head(all_results)

# 결과 저장 (옵션)
# write_csv(all_results, "HbA1c_sim_results.csv")


library(ggplot2)

# HbA1c 시계열 곡선 그래프
ggplot(all_results, aes(x = time, y = HbA1c, color = UID)) +
  geom_line(size = 1) +
  labs(
    title = "Time-course of HbA1c per subject",
    x = "Time (days)",
    y = "HbA1c (%)",
    color = "Subject ID"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
