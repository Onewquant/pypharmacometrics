# ===== 패키지 / 작업 폴더 =====
setwd('C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results/b1da/mvlreg_output_glm/r')

# install.packages(c("lme4", "broom.mixed", "dplyr", "readr"))
library(readr)
library(dplyr)
library(lme4)
library(broom.mixed)

# ===== 엔드포인트 리스트 =====
endpoints <- c("PLT", "WBC", "ANC", "Hb", "Lactate")

# ===== (도우미) summary 텍스트를 여러 페이지로 출력하는 함수 =====
draw_text_pages <- function(lines, title = NULL, lines_per_page = 55, cex = 0.62) {
  if (!is.null(title)) {
    lines <- c(paste0("### ", title, " ###"), "", lines)
  }
  n <- length(lines)
  if (n == 0) {
    plot.new()
    title(main = title)
    text(0.02, 0.98, labels = "(no content)", adj = c(0,1), cex = cex)
    return(invisible(NULL))
  }
  idx <- seq_len(n)
  chunks <- split(lines, ceiling(idx / lines_per_page))
  opar <- par(family = "mono")
  on.exit(par(opar), add = TRUE)
  for (chunk in chunks) {
    plot.new()
    y <- 0.98
    line_height <- 0.018
    if (cex != 1) line_height <- line_height * (0.7/cex)  # 글씨 크기에 따라 라인 간격 보정
    for (ln in chunk) {
      text(0.02, y, ln, adj = c(0,1), cex = cex)
      y <- y - line_height
      if (y < 0.02) break
    }
  }
}

# ===== PDF 시작 (엔드포인트별 summary 페이지만 저장) =====
pdf_file <- "GLM_Summary_only.pdf"
pdf(pdf_file, width = 8.27, height = 11.69)  # A4 세로

for (ep in endpoints) {
  message("Processing endpoint: ", ep)
  fn <- sprintf("b1da_lnz_mvlreg_datasubset(%s).csv", ep)
  
  if (!file.exists(fn)) {
    plot.new(); title(main = paste0(ep, " - File not found"))
    next
  }
  
  # 1) 데이터 로드
  df <- read_csv(fn, show_col_types = FALSE)
  
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
  if (length(predictors) == 0) {
    plot.new(); title(main = paste0(ep, " - No predictors"))
    next
  }
  
  # 4) 완전사례
  df_model <- df |>
    filter(complete.cases(y, PID, across(all_of(predictors))))
  
  if (nrow(df_model) < 10) {
    plot.new(); title(main = paste0(ep, " - Too few complete cases (n<10)"))
    next
  }
  
  # 5) 포뮬러 생성
  # fml <- as.formula(
  #   paste0("y ~ ", paste(predictors, collapse = " + "), " + (1 | PID)")
  # )
  fml <- as.formula(
    paste0("y ~ ", paste(predictors, collapse = " + "), "")
  )
  
  # 6) GLMM 적합 (오류 발생 시에도 페이지로 남김)
  # fit_glmm <- tryCatch({
  #   glmer(
  #     formula = fml,
  #     data    = df_model,
  #     family  = binomial(link = "logit"),
  #     control = glmerControl(optimizer = "bobyqa",
  #                            optCtrl   = list(maxfun = 1e5))
  #   )
  # }, error = function(e) e)
  # 
  fit_glmm <- tryCatch({
    glm(
      formula = fml,
      data    = df_model,
      family  = binomial(link = "logit"))
  }, error = function(e) e)
  
  if (inherits(fit_glmm, "error")) {
    plot.new(); title(main = paste0(ep, " - Model failed"))
    m <- strwrap(conditionMessage(fit_glmm), width = 100)
    text(0.02, 0.95, paste(m, collapse = "\n"), adj = c(0,1), cex = 0.8)
    next
  }
  
  # 7) summary 텍스트 캡처 후 PDF에 출력
  sum_lines <- capture.output(summary(fit_glmm))
  draw_text_pages(sum_lines, title = paste0("GLMM summary - ", ep), lines_per_page = 55, cex = 0.62)
}

dev.off()
message("PDF written: ", pdf_file)
