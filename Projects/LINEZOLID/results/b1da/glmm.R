set.seed(123)

n_site <- 5              # 병원 5개
n_patient <- 200          # 전체 환자 200명
patients_per_site <- n_patient / n_site

# 병원, 환자 ID 생성
site <- rep(paste0("Site", 1:n_site), each = patients_per_site)
patient_id <- paste0(site, "_P", sprintf("%03d", 1:n_patient))

# 개체 수준 공변량
age <- rnorm(n_patient, mean = 55, sd = 12)
sex <- rbinom(n_patient, 1, 0.45)
albumin <- rnorm(n_patient, mean = 4.0, sd = 0.5)
dose <- runif(n_patient, 200, 800)

# 병원 랜덤효과
site_effect <- rnorm(n_site, 0, 0.5)
names(site_effect) <- paste0("Site", 1:n_site)
rand_intercept <- site_effect[site]

# 로짓 선형예측식: ADR 발생 확률 생성
linpred <- -5 + 
  0.003 * dose +           # 용량 ↑ → 위험 ↑
  0.5 * sex +              # 여성 ↑ → 위험 ↑
  -0.6 * albumin +         # 알부민 ↑ → 위험 ↓
  0.02 * age + 
  rand_intercept

p_ADR <- plogis(linpred)
ADR <- rbinom(n_patient, 1, p_ADR)

df <- data.frame(patient_id, site, dose, albumin, age, sex, ADR)
head(df)

# remove.packages(c("lme4", "Matrix"))
# install.packages("Matrix")   # 먼저 Matrix 설치/업데이트
# install.packages("lme4")     # 그 다음 lme4 설치


# library(lme4)
# sessionInfo()

library(lme4)

fit_glmm <- glmer(
  ADR ~ dose + albumin + age + sex + (1 | site),
  data = df,
  family = binomial(link = "logit")
)

summary(fit_glmm)
