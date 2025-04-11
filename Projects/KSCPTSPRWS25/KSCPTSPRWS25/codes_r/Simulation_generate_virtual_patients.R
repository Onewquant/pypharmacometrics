# Simulation Data Prep---------------------------------

library(dplyr)
library(MASS)
library(readr)

setwd("C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25")
input_dir_path <- file.path(getwd(), "sglt2i_dataset")
results_dir_path <- file.path(getwd(), "results_r")

set.seed(123)

### 1. Input 데이터 불러오기 & GRP를 integer로 변환
input_data <- read_csv(file.path(input_dir_path, "KSCPTSPRWS25_SGLT2i_PD.csv")) %>%
  mutate(GRP = as.integer(GRP))

# print(mean(input_data$ALT, na.rm = TRUE))
# print(sd(input_data$ALT, na.rm = TRUE))

### 2. 그룹별 통계 요약 (mean, sd)
grp_summary <- input_data %>%
  group_by(GRP) %>%
  summarise(
    TBIL_mean = mean(TBIL, na.rm = TRUE),
    TBIL_sd = sd(TBIL, na.rm = TRUE),
    ALT_mean = mean(ALT, na.rm = TRUE),
    ALT_sd = sd(ALT, na.rm = TRUE),
    HbA1c_mean = mean(HbA1c_base, na.rm = TRUE),
    HbA1c_sd = sd(HbA1c_base, na.rm = TRUE),
    PG_base_mean = mean(PG_base, na.rm = TRUE),
    PG_base_sd = sd(PG_base, na.rm = TRUE)
  )

### 3. 시뮬레이션 설정값

theta <- c(E0 = 0.0307, EMAX = 0.512, EC50 = 58.2, HILL = 7.8)
omega <- c(0.0383, 0.0169) # ETA1, ETA2
sigma <- c(0.00004) # proportional, additive



n_patients_per_group <- 500
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
    TBIL = abs(rnorm(n, mean = grp_stat$TBIL_mean, sd = grp_stat$TBIL_sd)),
    ALT = abs(rnorm(n, mean = grp_stat$ALT_mean, sd = grp_stat$ALT_sd)),
    HbA1c_base = abs(rnorm(n, mean = grp_stat$HbA1c_mean, sd = grp_stat$HbA1c_sd)),
    PG_base = abs(rnorm(n, mean = grp_stat$PG_base_mean, sd = grp_stat$PG_base_sd))
  ) %>%
    mutate(
      PG_avg = PG_base * 1.5412 - 43.3732 + rnorm(n, mean = 0, sd = (0.1 * grp_stat$PG_base_sd))
    )
}

### 5. 가상환자 데이터 생성
virtual_patients <- bind_rows(
  lapply(names(eGFR_ranges), function(grp_no) {
    generate_group_data(grp_no, eGFR_ranges[[grp_no]], n_patients_per_group, grp_summary)
  })
)

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
    dUGEc = IPRED + EPS1
  ) %>%
  ungroup()

### 9. 결과 확인 & 저장
print(head(virtual_patients))

write.csv(virtual_patients, file.path(input_dir_path, "KSCPTSPRWS25_SGLT2i_Virtual_Patients.csv"),row.names = FALSE)
