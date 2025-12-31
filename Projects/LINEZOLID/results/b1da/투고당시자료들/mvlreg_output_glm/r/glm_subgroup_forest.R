# ============================================================
# Sex-stratified risk factor analysis (R)
# - Input:  /mnt/data/b1da_lnz_mvlreg_datasubset(Lactate).csv
# - Target: y (0/1)
# - Stratify: SEX == 0, 1
# - Outputs:
#     ./sex_stratified_or_SEX0.csv
#     ./sex_stratified_or_SEX1.csv
#     ./sex_stratified_or_summary.xlsx
#     ./forest_SEX0.png, ./forest_SEX0.pdf
#     ./forest_SEX1.png, ./forest_SEX1.pdf
# - Notes:
#     * 수치형 공변량만 사용 (y, SEX, ENDPOINT, ID-like 제외)
#     * 서브그룹 내 무변동/상수 컬럼 자동 제외
#     * 중앙값 대체(median imputation)
#     * CI는 Wald(confint.default) 기반 (빠르고 안정적)
# ============================================================

# 패키지 로드 ------------------------------------------------
# install.packages(c("openxlsx", "purrr"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(broom)
  library(ggplot2)
  library(openxlsx)
  library(purrr)
})

# ===== 패키지 / 작업 폴더 =====
setwd('C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results/b1da/mvlreg_output_glm/r')

# 사용자 설정 ------------------------------------------------
in_path  <- "b1da_lnz_mvlreg_datasubset(Lactate).csv"
out_dir  <- "./subgroup_analysis"   # 원하는 경로로 변경 가능
target   <- "y"
sex_col  <- "SEX"
exclude_like_ids <- c("pid","id","uid","subject","index") # 식별자 패턴

# 유틸 함수 --------------------------------------------------
# 중앙값 대체(수치형만)
median_impute <- function(df) {
  num_cols <- names(df)[sapply(df, is.numeric)]
  for (nm in num_cols) {
    m <- stats::median(df[[nm]], na.rm = TRUE)
    if (is.finite(m)) df[[nm]][is.na(df[[nm]])] <- m
  }
  df
}

# OR 테이블 생성 (Wald CI)
or_table_from_glm <- function(fit) {
  coefs <- stats::coef(summary(fit))
  # 계수표가 없는 경우 빈 DF 반환
  if (is.null(coefs) || nrow(coefs) == 0) {
    return(tibble(Variable=character(), Beta=double(), OR=double(),
                  CI_low=double(), CI_high=double(), p_value=double(), q_value=double()))
  }
  # confint.default는 Wald 기반 CI
  ci <- suppressWarnings(confint.default(fit))
  # 계수/SE/p 추출
  df <- tibble::rownames_to_column(as.data.frame(coefs), var = "term") %>%
    rename(Estimate = Estimate, StdErr = `Std. Error`, z = `z value`, p_value = `Pr(>|z|)`) %>%
    as_tibble()
  
  # CI 결합
  ci_df <- tibble::rownames_to_column(as.data.frame(ci), var = "term") %>%
    rename(CI_low_log = `2.5 %`, CI_high_log = `97.5 %`)
  df <- df %>% left_join(ci_df, by = "term")
  
  # OR/CI 변환 및 정리
  df <- df %>%
    mutate(
      Variable = term,
      Beta     = Estimate,
      OR       = exp(Estimate),
      CI_low   = exp(CI_low_log),
      CI_high  = exp(CI_high_log)
    ) %>%
    select(Variable, Beta, OR, CI_low, CI_high, p_value)
  
  # const(Intercept) 제거
  df <- df %>%
    filter(!grepl("^(Intercept|\\(Intercept\\)|const)$", Variable))
  
  # FDR
  if (nrow(df) > 0) {
    df$q_value <- p.adjust(df$p_value, method = "BH")
  } else {
    df$q_value <- numeric(0)
  }
  df
}

# 포레스트 플롯 (단일 성별)
save_forest <- function(or_df, title, png_path, pdf_path) {
  if (nrow(or_df) == 0) return(invisible(NULL))
  # OR=1에서의 거리( |log(OR)| )로 정렬
  plot_df <- or_df %>%
    mutate(dist = abs(log(OR))) %>%
    arrange(dist) %>%
    mutate(Variable = factor(Variable, levels = Variable))
  
  p <- ggplot(plot_df, aes(x = OR, y = Variable)) +
    geom_point() +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_continuous(trans = "log10") +
    labs(x = "Odds Ratio (95% CI)", y = NULL, title = title) +
    theme_minimal(base_size = 12)
  
  ggsave(png_path, plot = p, width = 7, height = max(4, 0.35*nrow(plot_df)), dpi = 300)
  ggsave(pdf_path, plot = p, width = 7, height = max(4, 0.35*nrow(plot_df)))
}

# 메인 분석 ---------------------------------------------------
# 1) 데이터 로드
df <- readr::read_csv(in_path, show_col_types = FALSE) %>%
  mutate(across(everything(), ~ .x)) # 컬럼 이름/타입 보존

# 2) 기본 전처리
stopifnot(target %in% names(df), sex_col %in% names(df))
df <- df %>%
  mutate(
    !!target := suppressWarnings(as.numeric(.data[[target]])),
    !!sex_col := suppressWarnings(as.numeric(.data[[sex_col]]))
  )

# 제외 컬럼 규칙: 타깃, 성별, ENDPOINT, ID-like
exclude_cols <- c(target, sex_col, "ENDPOINT")
exclude_cols <- unique(c(
  exclude_cols,
  names(df)[vapply(names(df), function(nm) any(stringr::str_detect(tolower(nm), paste0(exclude_like_ids, collapse="|"))), logical(1))]
))

# 후보 예측변수: 수치형 + 제외목록 제외, 전체에서 무변동 제외
num_candidates <- names(df)[sapply(df, is.numeric)]
num_candidates <- setdiff(num_candidates, exclude_cols)
num_candidates <- num_candidates[vapply(num_candidates, function(nm) dplyr::n_distinct(df[[nm]], na.rm=TRUE) > 1, logical(1))]

# 3) 성별 서브그룹별 모델 적합
sex_values <- sort(unique(na.omit(df[[sex_col]])))
sex_values <- sex_values[sex_values %in% c(0,1)]  # 0/1만 사용

results_list <- list()
for (sx in sex_values) {
  sub <- df %>% filter(.data[[sex_col]] == sx)
  
  # 서브그룹 내 무변동 제거
  cols_ok <- num_candidates[vapply(num_candidates, function(nm) dplyr::n_distinct(sub[[nm]], na.rm=TRUE) > 1, logical(1))]
  if (length(cols_ok) == 0) {
    warning(sprintf("SEX=%s: no varying predictors after filtering.", sx))
    or_tab <- tibble(SEX = sx, N = sum(!is.na(sub[[target]])), Variable=character(),
                     Beta=double(), OR=double(), CI_low=double(), CI_high=double(),
                     p_value=double(), q_value=double(), FittedWith=character())
    results_list[[as.character(sx)]] <- or_tab
    next
  }
  
  # X/Y 구성 + 중앙값 대체
  dat <- sub %>%
    select(all_of(c(target, cols_ok))) %>%
    median_impute() %>%
    drop_na(any_of(target))
  
  # glm 적합 (binomial)
  form <- as.formula(paste(target, "~", paste(cols_ok, collapse = " + ")))
  fit <- try(glm(form, data = dat, family = binomial(link = "logit")), silent = TRUE)
  
  # 완전분리/수렴문제 등으로 실패 시, 약간의 ridge 성격을 주는 방법(작은 정규화)은 base glm에 없음.
  # 여기선 실패 시, stepwise로 차원 축소 후 재시도.
  if (inherits(fit, "try-error")) {
    # 단변량 p<0.20 변수만 선별 후 재시도
    uv <- lapply(cols_ok, function(v) {
      fm <- as.formula(paste(target, "~", v))
      mdl <- try(glm(fm, data = dat, family = binomial()), silent = TRUE)
      pv <- if (inherits(mdl, "try-error")) 1 else tryCatch(coef(summary(mdl))[2,4], error = function(e) 1)
      tibble(var=v, p=pv)
    }) %>% bind_rows()
    sel <- uv %>% filter(is.finite(p), p < 0.20) %>% pull(var)
    if (length(sel) >= 1) {
      form2 <- as.formula(paste(target, "~", paste(sel, collapse = " + ")))
      fit <- glm(form2, data = dat, family = binomial())
      cols_used <- sel
    } else {
      # 단변량에서도 아무것도 선택 안 되면 빈 결과로
      cols_used <- character(0)
    }
  } else {
    cols_used <- all.vars(form)[-1]
  }
  
  if (length(cols_used) == 0 || inherits(fit, "try-error")) {
    or_df <- tibble(Variable=character(), Beta=double(), OR=double(),
                    CI_low=double(), CI_high=double(), p_value=double(), q_value=double())
  } else {
    or_df <- or_table_from_glm(fit)
  }
  
  # 메타 컬럼 추가
  or_df <- or_df %>%
    mutate(SEX = sx, .before = 1) %>%
    mutate(N = nrow(dat), .after = "SEX") %>%
    mutate(FittedWith = "glm(binomial)", .after = "q_value")
  
  # 저장
  out_csv <- file.path(out_dir, sprintf("sex_stratified_or_SEX%d.csv", as.integer(sx)))
  readr::write_csv(or_df, out_csv)
  
  # 포레스트
  png_path <- file.path(out_dir, sprintf("forest_SEX%d.png", as.integer(sx)))
  pdf_path <- file.path(out_dir, sprintf("forest_SEX%d.pdf", as.integer(sx)))
  save_forest(or_df, title = sprintf("Sex %d - Lactic Acidosis Risk Factors", as.integer(sx)),
              png_path = png_path, pdf_path = pdf_path)
  
  results_list[[as.character(sx)]] <- or_df
}

# 4) 요약 엑셀 저장 (시트: SEX0, SEX1, ALL)
all_df <- bind_rows(results_list, .id = NULL)
xlsx_path <- file.path(out_dir, "sex_stratified_or_summary.xlsx")
wb <- createWorkbook()
for (sx in sex_values) {
  addWorksheet(wb, sprintf("SEX%d", as.integer(sx)))
  writeData(wb, sprintf("SEX%d", as.integer(sx)), results_list[[as.character(sx)]])
}
addWorksheet(wb, "ALL")
writeData(wb, "ALL", all_df)
saveWorkbook(wb, xlsx_path, overwrite = TRUE)

message("Saved files:")
message(" - ", file.path(out_dir, "sex_stratified_or_SEX0.csv"))
message(" - ", file.path(out_dir, "sex_stratified_or_SEX1.csv"))
message(" - ", xlsx_path)
message(" - ", file.path(out_dir, "forest_SEX0.png"), ", ", file.path(out_dir, "forest_SEX0.pdf"))
message(" - ", file.path(out_dir, "forest_SEX1.png"), ", ", file.path(out_dir, "forest_SEX1.pdf"))

# ------------------------------------------------------------
# (선택) 유의 시그널만 별도 저장: q<0.05 & (OR>1.5 | OR<0.67)
# ------------------------------------------------------------
# sig_list <- lapply(results_list, function(df) {
#   df %>% filter(q_value < 0.05, (OR > 1.5 | OR < 0.67))
# })
# sig_all <- bind_rows(sig_list, .id = NULL)
# write_csv(sig_all, file.path(out_dir, "sex_stratified_sig_signals.csv"))
