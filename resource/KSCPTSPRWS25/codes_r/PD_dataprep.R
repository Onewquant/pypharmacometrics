library(readr)
library(dplyr)
library(ggplot2)
library(NonCompart)

setwd("C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25")
input_dir_path <- file.path(getwd(), "sglt2i_dataset")
results_dir_path <- file.path(getwd(), "results_r")

# Blood PD prep ----------------------------------------------------

# Blood PD raw 데이터 불러오기
bpd_df <- read_csv(file.path(input_dir_path, "KSCPTSPRWS25_SGLT2i_Blood_PD.csv"))

# Baseline 정보 추출
bpd_baseline_df <- bpd_df %>%
  filter(N_TIME == 0) %>%
  select(ID, N_TIME, `C_Serum GLU`, HbA1c) %>%
  rename(PG_base = `C_Serum GLU`, HbA1c_base = HbA1c)

# Glucose AUC 계산
glu_auc_df <- tblNCA(bpd_df, key = c("ID"), colTime = "N_TIME", colConc = "C_Serum GLU", dose=0.5, concUnit="ug/L", down = "Log")
glu_auc_df <- glu_auc_df %>% select(ID, AUCLST) %>% rename(ID = ID, AUC_glu = AUCLST )
glu_auc_df <- glu_auc_df %>% mutate(PG_avg = as.numeric(AUC_glu) / 24)


# baseline과 병합
glu_auc_df <- left_join(glu_auc_df, bpd_baseline_df %>% select(ID, PG_base, HbA1c_base), by = "ID")

# Urine PD prep ----------------------------------------------------

# Urine PD raw 데이터 불러오기
upd_df <- read_csv(file.path(input_dir_path, "KSCPTSPRWS25_SGLT2i_Urine_PD.csv"))

# 24시간 수집량 계산
uge24_df <- upd_df %>%
  group_by(ID, N_DAY) %>%
  summarise(
    UGE24 = sum(Amount_GLU, na.rm = TRUE),
    VOLUME24 = sum(Volume, na.rm = TRUE),
    .groups = "drop"
  )

# Baseline(-1일) 데이터
ugebase_df <- uge24_df %>%
  filter(N_DAY == -1) %>%
  rename(UGEbase = UGE24, VOLUMEbase = VOLUME24) %>%
  select(-N_DAY)

# Day 1 데이터와 병합
uge24_df <- uge24_df %>%
  filter(N_DAY == 1) %>%
  left_join(ugebase_df, by = "ID")

# COVAR prep ----------------------------------------------------

# COVAR raw 데이터 불러오기
covar_df <- read_csv(file.path(input_dir_path, "KSCPTSPRWS25_SGLT2i_COVAR.csv"))


# Generate final prep-PD df -------------------------------------

# 병합
pd_res_df <- glu_auc_df %>%
  left_join(uge24_df, by = "ID") %>%
  left_join(covar_df, by = "ID")

# Gen Col
pd_res_df = pd_res_df %>% mutate(dUGEc = (UGE24-UGEbase) / PG_avg)
pd_res_df = pd_res_df %>% mutate(EFFECT1 = dUGEc)

# 결과 확인
pd_res_df
write.csv(pd_res_df, file.path(input_dir_path, "KSCPTSPRWS25_SGLT2i_PD.csv"),row.names = FALSE)