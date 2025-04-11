library(deSolve)
library(tidyverse)
library(readr)

# 1. 입력 데이터
input_data <- read_csv("PD_data.csv")  # UID, dUGEc 포함

# 2. 모델 파라미터 (TV, RSE, IIV)
param_table <- tibble::tribble(
  ~param,      ~tv,     ~rse,   ~iiv,
  "Pfmax",     15,      10,     20,
  "Kfp",       0.01,    15,     25,
  "DISfp",     0.05,    8,      30,
  "SLOPEfd",  -0.5,     12,     40,
  "Phmax",     0.5,     10,     20,
  "Khp",       0.01,    15,     25,
  "DIShp",     0.005,   10,     30,
  "Kin2",      0.01,    10,     20,
  "Kout",      0.02,    10,     20
)

# 3. 개인별 파라미터 생성 함수
generate_subject_parameters <- function(param_table) {
  param_table %>%
    rowwise() %>%
    mutate(
      est_with_rse = rnorm(1, mean = tv, sd = (rse / 100) * tv),
      value = est_with_rse * exp(rnorm(1, mean = 0, sd = sqrt(log(1 + (iiv / 100)^2))))
    ) %>%
    select(param, value) %>%
    deframe()
}

# 4. ODE 모델
model_fn <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    FPGplacebo <- FPGbaseline + Pfmax * (1 - exp(-Kfp * t)) + DISfp * t
    FPG <- FPGplacebo + SLOPEfd * dUGEc
    HbA1cplacebo <- HbA1cbaseline - Phmax * (1 - exp(-Khp * t)) + DIShp * t
    Kin <- Kout * HbA1cbaseline - Kin2
    dHbA1cdrug <- (FPG / FPGbaseline) * Kin + Kin2 - Kout * HbA1cdrug
    
    list(c(dHbA1cdrug = dHbA1cdrug),
         HbA1c = HbA1cplacebo + HbA1cdrug)
  })
}

# 5. 시뮬레이션 반복 수 설정
nrep <- 100  # 반복 횟수
times <- seq(0, 180, by = 1)

# 6. 시뮬레이션 실행
all_sims <- list()

for (i in 1:nrow(input_data)) {
  subj <- input_data[i, ]
  
  for (rep in 1:nrep) {
    base_params <- list(
      FPGbaseline = 120,
      HbA1cbaseline = 8.0,
      dUGEc = subj$dUGEc
    )
    rand_params <- generate_subject_parameters(param_table)
    parms <- c(base_params, rand_params)
    
    state <- c(HbA1cdrug = 0)
    
    out <- ode(y = state, times = times, func = model_fn, parms = parms) %>%
      as.data.frame() %>%
      mutate(UID = subj$UID, rep = rep)
    
    all_sims[[length(all_sims) + 1]] <- out
  }
}

# 7. 결과 통합
sim_data <- bind_rows(all_sims)

# 8. 신뢰구간 계산
summary_data <- sim_data %>%
  group_by(UID, time) %>%
  summarise(
    HbA1c_median = median(HbA1c),
    HbA1c_low = quantile(HbA1c, 0.025),
    HbA1c_high = quantile(HbA1c, 0.975),
    .groups = "drop"
  )

# 9. 시각화: 신뢰구간 포함 HbA1c 시뮬레이션
ggplot(summary_data, aes(x = time, y = HbA1c_median)) +
  geom_ribbon(aes(ymin = HbA1c_low, ymax = HbA1c_high), fill = "lightblue", alpha = 0.4) +
  geom_line(color = "blue", size = 1) +
  facet_wrap(~ UID) +
  labs(
    title = "HbA1c Simulation with RSE + IIV (95% CI)",
    x = "Time (days)",
    y = "HbA1c (%)"
  ) +
  theme_minimal()
