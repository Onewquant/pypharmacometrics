library(dplyr)
library(MASS)
library(readr)

set.seed(123)

### 1. Input 데이터 불러오기 & GRP를 integer로 변환
input_data <- read_csv("KSCPTSPR25_PD_Endpoint_Sim_data.csv") %>%
  mutate(GRP = as.integer(GRP))

### 2. 그룹별 통계 요약 (mean, sd)
grp_summary <- input_data %>%
  group_by(GRP) %>%
  summarise(
    TBIL_mean = mean(TBIL, na.rm = TRUE),
    TBIL_sd = sd(TBIL, na.rm = TRUE),
    ALT_mean = mean(ALT, na.rm = TRUE),
    ALT_sd = sd(ALT, na.rm = TRUE),
    HbA1c_mean = mean(HbA1c, na.rm = TRUE),
    HbA1c_sd = sd(HbA1c, na.rm = TRUE),
    PG_ZERO_mean = mean(PG_ZERO, na.rm = TRUE),
    PG_ZERO_sd = sd(PG_ZERO, na.rm = TRUE)
  )

### 3. 시뮬레이션 설정값
theta <- c(E0 = 0.025, EMAX = 0.427, EC50 = 59.9, HILL = 7.17)
omega <- c(0.157, 0.141) # ETA1, ETA2
sigma <- c(0.00066, 0.00002) # proportional, additive

n_patients_per_group <- 300
eGFR_ranges <- list(
  "1" = c(90, 120),
  "2" = c(60, 90),
  "3" = c(30, 60),
  "4" = c(15, 30)
)

### 4. 가상환자 데이터 생성 함수
generate_group_data <- function(grp_no, range, n, summary_df) {
  grp_stat <- summary_df %>% filter(GRP == as.integer(grp_no))
  
  data.frame(
    ID = seq_len(n) + (as.integer(grp_no) - 1) * n,
    GRP = as.integer(grp_no),
    eGFR = runif(n, min = range[1], max = range[2]),
    TBIL = rnorm(n, mean = grp_stat$TBIL_mean, sd = grp_stat$TBIL_sd),
    ALT = rnorm(n, mean = grp_stat$ALT_mean, sd = grp_stat$ALT_sd),
    HbA1c = rnorm(n, mean = grp_stat$HbA1c_mean, sd = grp_stat$HbA1c_sd),
    PG_ZERO = rnorm(n, mean = grp_stat$PG_ZERO_mean, sd = grp_stat$PG_ZERO_sd)
  ) %>%
    mutate(
      PG_AVG = PG_ZERO * 1.5412 - 43.3732 + rnorm(n, mean = 0, sd = (0.1 * grp_stat$PG_ZERO_sd))
    )
}

### 5. 가상환자 데이터 생성
virtual_patients <- bind_rows(
  lapply(names(eGFR_ranges), function(grp_no) {
    generate_group_data(grp_no, eGFR_ranges[[grp_no]], n_patients_per_group, grp_summary)
  })
) %>%
  mutate(eGFRxTBIL = eGFR * TBIL)

### 6. IIV 추가 (ETA)
eta <- mvrnorm(n = nrow(virtual_patients), mu = c(0, 0),
               Sigma = diag(omega))

virtual_patients$ETA1 <- eta[, 1]
virtual_patients$ETA2 <- eta[, 2]

### 7. IPRED 계산 (Sigmoid Emax model)
virtual_patients <- virtual_patients %>%
  mutate(
    E0 = theta["E0"],
    EMAX = theta["EMAX"] * exp(ETA1),
    EC50 = theta["EC50"] * exp(ETA2),
    HILL = theta["HILL"],
    IPRED = E0 + EMAX * eGFR^HILL / (EC50^HILL + eGFR^HILL)
  )

### 8. Residual variability 추가 & 최종 dUGEc 생성
virtual_patients <- virtual_patients %>%
  rowwise() %>%
  mutate(
    EPS1 = rnorm(1, mean = 0, sd = sqrt(sigma[1])),
    EPS2 = rnorm(1, mean = 0, sd = sqrt(sigma[2])),
    dUGEc = IPRED * (1 + EPS1) + EPS2
  ) %>%
  ungroup()

### 9. 결과 확인 & 저장
head(virtual_patients)

write.csv(virtual_patients, "C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25/results/virtual_patients_for_simulation.csv", row.names = FALSE)
