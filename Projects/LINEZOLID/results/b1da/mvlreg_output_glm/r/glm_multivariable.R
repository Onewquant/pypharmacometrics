# ===== 패키지 / 작업 폴더 =====
setwd('C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results/b1da/mvlreg_output_glm/r')

library(broom)
library(readr)
library(dplyr)
library(lme4)
library(broom.mixed)

# ===== 엔드포인트 리스트 =====
endpoints <- c("PLT", "WBC", "ANC", "Hb", "Lactate")

# ===== 분석 파라미터 =====
OR_upper_cutoff <- 1.5      # OR > 이 값이면 강조
OR_lower_cutoff <- 0.67     # OR < 이 값이면 강조
p_adj_cutoff    <- 0.05     # 보정된 p-value 기준
output_pdf      <- "GLM_Summary_only.pdf"
output_csv      <- "GLM_Significant_OR_Table.csv"

# ===== 시그널 모으는 DF 초기화 =====
sig_or_table <- data.frame(
  Endpoint = character(),
  Variable = character(),
  OR = numeric(),
  `2.5%` = numeric(),
  `97.5%` = numeric(),
  p.value = numeric(),
  adj.p.value = numeric(),
  sig = character(),
  check.names = FALSE
)

# ===== PDF 시작 =====
pdf(output_pdf, width = 8.5, height = 11)

for (ep in endpoints) {
  message("Processing endpoint: ", ep)
  fn <- sprintf("b1da_lnz_mvlreg_datasubset(%s).csv", ep)
  
  # 1) 데이터 로드
  df <- read_csv(fn, show_col_types = FALSE)
  if (!is.factor(df$SEX)) df$SEX <- factor(df$SEX)
  
  # 참조범주 고정: Male
  df$SEX <- factor(df$SEX, levels = c(0,1), labels = c("Male","Female"))
  df$SEX <- relevel(df$SEX, ref = "Male")
  
  # 2) 전처리 (요청: 기존 코드 유지)
  df <- df |>
    select(-ENDPOINT) |>
    mutate(
      PID = as.factor(PID),
      across(where(is.character), as.factor)
    )
  
  # CUM_DOSE만 600으로 나눈 값으로 대체 (기존 요청 반영)
  if ("CUM_DOSE" %in% names(df)) {
    df <- df %>% mutate(CUM_DOSE = CUM_DOSE / 600)
  }
  
  # 3) 예측자 목록 (y, PID 제외)
  predictors <- setdiff(names(df), c("y", "PID"))
  
  fml <- as.formula(
    paste0("y ~ ", paste(predictors, collapse = " + "), "")
  )
  
  fit_glm <- tryCatch({
    glm(
      formula = fml,
      data    = df,
      family  = binomial(link = "logit"))
  }, error = function(e) e)
  
  # ===== OR 결과 PDF 저장 =====
  plot.new()
  title_line <- paste0("Multivariable Logistic Regression (Endpoint: ", ep, ")")
  mtext(title_line, side = 3, adj = 0, line = 1, cex = 1.2, font = 2)
  mtext(paste0("File: ", fn), side = 3, adj = 0, line = 0, cex = 0.9)
  
  if (inherits(fit_glm, "glm")) {
    # 회귀계수 및 CI 계산 (Wald 방식)
    coefs <- coef(fit_glm)
    ses   <- sqrt(diag(vcov(fit_glm)))
    LCL   <- coefs - 1.96 * ses
    UCL   <- coefs + 1.96 * ses
    OR    <- exp(coefs)
    LCL_e <- exp(LCL)
    UCL_e <- exp(UCL)
    pval  <- summary(fit_glm)$coefficients[, "Pr(>|z|)"]
    padj  <- p.adjust(pval, method = "BH")
    sig   <- ifelse(padj < p_adj_cutoff, "*", "")
    
    or_tbl <- data.frame(
      Variable = names(OR),
      OR = OR,
      `2.5%` = LCL_e,
      `97.5%` = UCL_e,
      p.value = pval,
      adj.p.value = padj,
      sig = sig,
      check.names = FALSE
    )
    
    # === 유의 OR 필터링 → 누적 테이블에 저장 ===
    idx_sig <- (or_tbl$OR > OR_upper_cutoff | or_tbl$OR < OR_lower_cutoff) & 
      (or_tbl$`adj.p.value` < p_adj_cutoff)
    if (any(idx_sig, na.rm = TRUE)) {
      tmp <- or_tbl[idx_sig, , drop = FALSE]
      tmp$Endpoint <- ep
      tmp <- tmp[, c("Endpoint","Variable","OR","2.5%","97.5%","p.value","adj.p.value","sig")]
      sig_or_table <- rbind(sig_or_table, tmp)
    }
    
    # 보기 좋게 포맷팅
    fmt_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
    or_tbl_disp <- or_tbl
    or_tbl_disp$OR <- fmt_num(or_tbl_disp$OR)
    or_tbl_disp$`2.5%` <- fmt_num(or_tbl_disp$`2.5%`)
    or_tbl_disp$`97.5%` <- fmt_num(or_tbl_disp$`97.5%`)
    or_tbl_disp$p.value <- ifelse(or_tbl_disp$p.value < 0.0001, "<0.0001",
                                  fmt_num(or_tbl_disp$p.value, 4))
    or_tbl_disp$adj.p.value <- ifelse(or_tbl_disp$adj.p.value < 0.0001, "<0.0001",
                                      fmt_num(or_tbl_disp$adj.p.value, 4))
    
    # 텍스트 출력
    y <- 0.9
    mf <- paste("Model formula:", paste(deparse(fml), collapse = ""))
    wrapped <- strwrap(mf, width = 110)
    line_h <- 0.018
    for (ln in wrapped) {
      text(0.02, y, ln, adj = c(0, 1), cex = 0.8)
      y <- y - line_h
      if (y < 0.05) { plot.new(); y <- 0.95 }
    }
    y <- y - 0.012
    text(0.02, y, paste0("AIC: ", round(AIC(fit_glm),2), 
                         " | BIC: ", round(BIC(fit_glm),2),
                         " | logLik: ", round(logLik(fit_glm),2)), 
         adj = c(0,1), cex = 0.8)
    y <- y - 0.05
    text(0.02, y, "Odds Ratio (OR) with 95% Wald CI", adj = c(0,1), cex = 0.9, font = 2)
    
    # 테이블 출력
    y <- y - 0.03
    line_height <- 0.02
    header <- c("Variable", "OR", "2.5%", "97.5%", "p.value", "adj.p.value", "sig")
    text(0.02, y, paste(sprintf("%-16s", header), collapse=" "), adj=c(0,1), cex=0.8, font=2)
    y <- y - 0.02
    for (i in 1:nrow(or_tbl_disp)) {
      ln <- paste(sprintf("%-16s", or_tbl_disp[i, ]), collapse=" ")
      y <- y - line_height
      if (y < 0.05) { plot.new(); y <- 0.95 }
      text(0.02, y, ln, adj=c(0,1), cex=0.75)
    }
  }
}

# ===== 마지막 페이지: Significant OR 요약 =====
plot.new()
title_line <- "Significant OR Summary across Endpoints"
mtext(title_line, side = 3, adj = 0, line = 1, cex = 1.2, font = 2)
criteria_line <- paste0("(Criteria: OR >", OR_upper_cutoff,
                        " or OR <", OR_lower_cutoff,
                        ") AND adj.p.value <", p_adj_cutoff, ")")
mtext(criteria_line, side = 3, adj = 0, line = 0, cex = 0.9)

y <- 0.88
line_h <- 0.02

if (nrow(sig_or_table) == 0) {
  text(0.02, y, "No significant rows matched the criteria.", adj = c(0,1), cex = 0.9)
} else {
  fmt_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
  disp <- sig_or_table
  disp$OR      <- fmt_num(disp$OR)
  disp$`2.5%`  <- fmt_num(disp$`2.5%`)
  disp$`97.5%` <- fmt_num(disp$`97.5%`)
  disp$p.value <- ifelse(disp$p.value < 0.0001, "<0.0001", fmt_num(disp$p.value, 4))
  disp$`adj.p.value` <- ifelse(disp$`adj.p.value` < 0.0001, "<0.0001", fmt_num(disp$`adj.p.value`, 4))
  
  header <- c("Endpoint","Variable","OR","2.5%","97.5%","p.value","adj.p.value","sig")
  text(0.02, y, paste(sprintf("%-14s", header), collapse=" "), adj=c(0,1), cex=0.8, font=2)
  y <- y - 0.02
  
  for (i in 1:nrow(disp)) {
    ln <- paste(sprintf("%-14s", disp[i, ]), collapse=" ")
    y <- y - line_h
    if (y < 0.05) { plot.new(); y <- 0.95 }
    text(0.02, y, ln, adj=c(0,1), cex=0.75)
  }
  
  try(write.csv(sig_or_table, output_csv, row.names = FALSE), silent = TRUE)
}

dev.off()
