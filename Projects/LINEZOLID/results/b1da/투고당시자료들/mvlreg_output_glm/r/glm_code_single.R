# 패키지
setwd('C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results/b1da/mvlreg_output_glm/r')

# install.packages(c("lme4", "broom.mixed", "dplyr", "readr"))  # 최초 1회
library(readr)
library(dplyr)
library(lme4)
library(broom.mixed)

# ===== 엔드포인트 리스트 =====
endpoints <- c("PLT", "WBC", "ANC", "Hb", "Lactate")

# 스케일링 대상 컬럼(존재하는 경우에만 적용)
# cols_to_scale <- c("CUM_DOSE")  # 필요 시 c("PLT", "ANC", "CUM_DOSE") 등으로 확장

# 1) 데이터 로드

# df <- read_csv("b1da_lnz_mvlreg_datasubset(PLT).csv")
# df <- read_csv("b1da_lnz_mvlreg_datasubset(WBC).csv")
# df <- read_csv("b1da_lnz_mvlreg_datasubset(ANC).csv")
# df <- read_csv("b1da_lnz_mvlreg_datasubset(Hb).csv")
df <- read_csv("b1da_lnz_mvlreg_datasubset(Lactate).csv")

# 2) 전처리
# - ENDPOINT 컬럼 제거
# - PID는 범주형(factor)으로
# - 문자형 컬럼은 factor로 변환
# - y가 0/1 형태인지 확인 (필요시 as.integer/ as.numeric 변환)
df <- df |>
  select(-ENDPOINT) |>
  mutate(
    PID = as.factor(PID),
    across(where(is.character), as.factor)
  )

# (선택) 수치형 공변량 표준화: 규모가 크게 다른 수치를 안정화하고 싶다면 주석 해제
# num_cols <- df |> select(-y, -PID) |> select(where(is.numeric)) |> names()
# df[num_cols] <- scale(df[num_cols])

# cols_to_scale <- c("PLT", "ANC", "CUM_DOSE")
# common_cols <- intersect(cols_to_scale, names(df))  # 실제 존재하는 컬럼만 선택
# df[common_cols] <- lapply(df[common_cols], scale)
df <- df %>%
  mutate(CUM_DOSE = CUM_DOSE / 600)


# 3) 예측자 목록 만들기 (y, PID 제외한 모든 컬럼을 fixed effect로)
predictors <- setdiff(names(df), c("y", "PID"))
if (length(predictors) == 0) stop("고정효과로 쓸 공변량이 없습니다.")

# 4) 결측치 제거(완전사례)
df_model <- df |>
  filter(complete.cases(y, PID, across(all_of(predictors))))

# 5) 포뮬러 생성: y ~ (고정효과들) + (1 | PID)
fml <- as.formula(
  paste0("y ~ ", paste(predictors, collapse = " + "), " + (1 | PID)")
)

# 6) GLMM 적합 (이항 로지트)
fit_glmm <- glmer(
  formula = fml,
  data    = df_model,
  family  = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl   = list(maxfun = 1e5))
)

# 7) 결과 요약
summary(fit_glmm)

# # 8) 고정효과 OR(95% CI) 테이블
# or_tbl <- tidy(fit_glmm, effects = "fixed",
#                conf.int = TRUE, conf.method = "Wald",
#                exponentiate = TRUE)
# 
# # adjusted p-value 추가 (예: Benjamini–Hochberg FDR 방식)
# or_tbl <- or_tbl %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))
# # or_tbl <- or_tbl %>% mutate(p.adj = p.adjust(p.value, method = "bonferroni"))
# 
# print(or_tbl)
# # columns: term, estimate(=OR), conf.low, conf.high, p.value
# 
# # 9) 랜덤효과(개체별 절편 분산)
# VarCorr(fit_glmm)

# (선택) 예측 확률과 분류
# df_model$pred_prob <- predict(fit_glmm, type = "response")
# df_model$pred_class <- as.integer(df_model$pred_prob >= 0.5)
# table(Observed = df_model$y, Predicted = df_model$pred_class)
