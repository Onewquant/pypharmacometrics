library(readr)
library(dplyr)
library(ggplot2)
library(NonCompart)

# setwd("C:/Users/ilma0/downloads/KSCPTSPRWS25")
input_dir_path <- file.path(getwd(), "sglt2i_dataset")
results_dir_path <- file.path(getwd(), "results_r")

# PKPD 데이터 불러오기
pk_nca_df <- read_csv(file.path(results_dir_path, "[WSCT] PK_analysis_data - NCARes.csv"))
pd_df <- read_csv(file.path(input_dir_path, "KSCPTSPRWS25_SGLT2i_PD.csv"))
pd_df <- left_join(pd_df, pk_nca_df, by = c("ID", "GRP"))

# Exposure (AUC) 와 EFFECT1 (dUGEc) 관계---------------------------------

library(broom)       # tidy model outputs
library(dplyr)
library(glue)
library(stringr)

# 반복할 변수들
x_vars <- c("AUClast", "eGFR", "TBIL", "ALT")
y_var <- "EFFECT1"
hue <- "GRP"

# 반복문
for (x in x_vars) {
  # 선형 회귀
  model <- lm(as.formula(paste(y_var, "~", x)), data = pd_df)
  model_summary <- summary(model)
  coef_vals <- coef(model)
  
  slope <- coef_vals[[x]]
  intercept <- coef_vals[["(Intercept)"]]
  r_squared <- model_summary$r.squared
  p_value <- coef(summary(model))[x, "Pr(>|t|)"]
  
  # 그래프 제목
  fig_title <- glue("[WSCT] {x} vs {y_var} by {hue}\nR-squared:{round(r_squared, 4)}, p-value: {round(p_value, 4)}\nbeta: {round(slope, 4)}, intercept: {round(intercept, 4)}")
  
  # 그래프 그리기
  p <- ggplot(pd_df, aes_string(x = x, y = y_var, color = hue)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = fig_title, x = x, y = y_var, color = hue) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right")
  
  # 파일명 만들기
  filename_base <- str_trim(str_split(fig_title, "R-squared", simplify = TRUE)[1])
  filename_safe <- str_replace(filename_base, "\\[WSCT\\]", "[WSCT] Scatter LR")
  filepath <- file.path(results_dir_path, paste0(filename_safe, ".png"))
  
  # 저장
  ggsave(filename = filepath, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
}


# Covariates heatmap ----------------------------------------------------

# install.packages("ggcorrplot")
library(ggcorrplot)

# 1. 공변량 후보 변수 지정

cov_cand <- c("AGE", "SEX", "HT", "WT", "BMI", "ALB", "eGFR", "AST", "ALT", "SODIUM", "TBIL", "AUClast")

# 2. 상관행렬 계산 (절댓값)
corr_matrix <- pd_df %>%
  dplyr::select(all_of(cov_cand)) %>%
  cor(use = "complete.obs") %>%
  abs()

# 3. 히트맵 시각화
fig_title <- "[WSCT] Correlation matrix heatmap of the covariates"

p <- ggcorrplot(
  corr_matrix,
  lab = TRUE,
  lab_size = 5,
  method = "square",
  type = "lower",
  colors = c("blue", "white", "red"),
  title = fig_title
) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15)
  )

# 4. 저장
ggsave(filename = file.path(results_dir_path, paste0(fig_title, ".png")),
       plot = p, width = 10, height = 8, dpi = 300, bg = "white")


# Covariates VIF --------------------------------------------------------

library(car)

# 1. 공변량 리스트
cov_cand <- c("AGE", "SEX", "HT", "WT", "BMI", "ALB", "eGFR", "AST", "ALT", "SODIUM", "TBIL", "AUClast")
total_cand <- cov_cand

# 2. 상관 행렬
corr_matrix <- pd_df %>%
  dplyr::select(all_of(cov_cand)) %>%
  cor(use = "complete.obs") %>%
  abs()

# 3. 상삼각 행렬 (자기 자신과 중복 제거)
lower_tri <- lower.tri(corr_matrix)
lower <- corr_matrix
lower[!lower_tri] <- NA

# 4. 상관계수 0.7 이상인 변수 제거 대상 추출
to_drop <- colnames(lower)[apply(lower, 2, function(col) any(col >= 0.7, na.rm = TRUE))]

# 5. 제거된 공변량을 제외한 X 구성
X <- pd_df %>%
  dplyr::select(setdiff(total_cand, to_drop)) %>%
  mutate(across(everything(), as.numeric))  # 모든 컬럼을 numeric으로

# 6. VIF 계산y_var <- names(X)[1]
x_vars <- names(X)[-1]

# VIF 계산
vif_values <- vif(lm(as.formula(paste('EFFECT1', "~", paste(x_vars, collapse = " + "))), data = pd_df))

df_vif <- data.frame(
  Covariates = names(vif_values),
  VIF = as.numeric(vif_values)
) %>%
  arrange(desc(VIF))

# 7. VIF 막대 그래프 그리기
fig_title <- "[WSCT] Variance Inflation Factor (VIF) of covariates"

p <- ggplot(df_vif, aes(x = reorder(Covariates, VIF), y = VIF)) +
  geom_bar(stat = "identity", fill = "royalblue") +
  geom_text(aes(label = sprintf("%.2f", VIF)), hjust = -0.1, size = 5) +
  coord_flip() +
  labs(
    title = fig_title,
    x = NULL,
    y = "VIF Value"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    panel.grid.major.x = element_line(linetype = "dashed", color = "gray70")
  )

# 8. 저장
ggsave(file.path(results_dir_path, paste0(fig_title, ".png")),
       plot = p, width = 10, height = 8, dpi = 300, bg = "white")


# Multivariable Linear Regression ---------------------------------------


effect_col = 'EFFECT1'

# 1. 종속변수 y 설정
y <- pd_df[[effect_col]]


# 3. 전체 변수 회귀 모델
formula_total <- as.formula(paste(effect_col, "~", paste(df_vif$Covariates, collapse = " + ")))
model_total <- lm(formula_total, data = pd_df)

# 5. 결과 출력
summary(model_total)


# PD fitting (최종 모델) ------------------------------------------


# 1. 기본 설정
x_col <- "eGFR"
effect_col <- "EFFECT1"

x_vals <- pd_df[[x_col]]
y_vals <- pd_df[[effect_col]]

# 2. 선형 회귀
model <- lm(y_vals ~ x_vals)
summary_model <- summary(model)

intercept <- coef(model)[1]
slope <- coef(model)[2]
r_squared <- summary_model$r.squared
p_value <- coef(summary(model))["x_vals", "Pr(>|t|)"]

cat(sprintf("Linear model: intercept = %.4f, slope = %.4f, R² = %.4f, p-value = %.4g\n",
            intercept, slope, r_squared, p_value))


# 3. Sigmoid Emax 모델 정의
sigmoid_emax <- function(conc, E0, Emax, EC50, H) {
  E0 + Emax * conc^H / (EC50^H + conc^H)
}

# 4. 초기값 설정 및 모델 피팅 (단순 fitting)

nls_fit <- nls(
  formula = EFFECT1 ~ E0 + Emax * eGFR^H / (EC50^H + eGFR^H),
  data = pd_df,
  start = list(E0 = 0.02, Emax = 0.5, EC50 = 50, H = 1),
  control = nls.control(maxiter = 500, warnOnly = TRUE)
)

# 5. 파라미터 추정값
popt <- coef(nls_fit)
E0_fit <- popt["E0"]
Emax_fit <- popt["Emax"]
EC50_fit <- popt["EC50"]
H_fit <- popt["H"]

cat(sprintf("Sigmoid Emax fit:\nE0 = %.3f\nEmax = %.3f\nEC50 = %.3f\nHill coefficient = %.2f\n",
            E0_fit, Emax_fit, EC50_fit, H_fit))

# Mixed-effect model fitting

x_col <- "eGFR"
y_col <- "EFFECT1"

# E0_fit    <- 0.0251
# Emax_fit  <- 0.435
# EC50_fit  <- 59.6
# H_fit     <- 7.36

E0_fit    <- 0.0307
Emax_fit  <- 0.512
EC50_fit  <- 58.2
H_fit     <- 7.8


# Fitting된 모델 시각화 ------------------------------------------

# 1. 관측값
observed_df <- pd_df %>%
  dplyr::select(all_of(c(x_col, y_col))) %>%
  rename(x = all_of(x_col), y = all_of(y_col)) %>%
  mutate(x = as.numeric(x), y = as.numeric(y)) %>%
  filter(!is.na(x), !is.na(y))

# 2. 모델 피팅된 예측값 생성
gfr_fit <- seq(min(observed_df$x), max(observed_df$x), length.out = 200)
effect_fit <- E0_fit + Emax_fit * gfr_fit^H_fit / (EC50_fit^H_fit + gfr_fit^H_fit)

fit_df <- data.frame(x = gfr_fit, y = effect_fit)

# 3. 제목 및 라벨 구성
fig_title <- glue::glue(
  "[WSCT] Sigmoid Emax Model Fit ({x_col} vs {y_col})\n",
  "E0: {round(E0_fit, 3)}, Emax: {round(Emax_fit, 3)}, EC50: {round(EC50_fit, 2)}, ",
  "Hill coefficient: {round(H_fit, 2)}, OFV: -112.05"
)

# 4. 그래프 그리기
p <- ggplot() +
  geom_point(data = observed_df, aes(x = x, y = y), color = "black", shape = 16, size = 2.5, alpha = 0.8) +
  geom_line(data = fit_df, aes(x = x, y = y), color = "royalblue", size = 1.2) +
  labs(
    title = fig_title,
    x = x_col,
    y = y_col
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14),
    legend.position = "none"
  )

print(p)

# 5. 파일 저장
filename_safe <- str_replace_all(str_split(fig_title, "E0:", simplify = TRUE)[1], "^[:alnum:]_ -", "")
ggsave(filename = file.path(results_dir_path, paste0(str_trim(filename_safe), ".png")),
       plot = p, width = 10, height = 8, dpi = 300, bg = "white")


