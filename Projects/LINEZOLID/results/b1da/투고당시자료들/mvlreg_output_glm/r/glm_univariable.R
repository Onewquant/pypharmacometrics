# ===== 패키지 / 작업 폴더 =====
setwd('C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/LINEZOLID/results/b1da/mvlreg_output_glm/r')

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(broom)        # tidy() 사용
  library(lme4)
})

# ===== 엔드포인트 리스트 =====
endpoints <- c("PLT", "WBC", "ANC", "Hb", "Lactate")

# ===== 분석 파라미터 =====
OR_upper_cutoff <- 1.5
OR_lower_cutoff <- 0.67
p_adj_cutoff    <- 0.05
output_pdf_uni  <- "GLM_Univariable_Summary.pdf"
output_csv_uni  <- "GLM_Univariable_Significant_OR_Table.csv"

# ===== 결과 누적 DF 초기화 (endpoint 전역 요약) =====
sig_or_table_uni <- data.frame(
  Endpoint = character(),
  Predictor = character(),  # 원 변수명
  Term = character(),       # (factor일 때 level 포함된 term)
  OR = numeric(),
  `2.5%` = numeric(),
  `97.5%` = numeric(),
  p.value = numeric(),
  adj.p.value = numeric(),
  sig = character(),
  check.names = FALSE
)

# ===== PDF 시작 =====
pdf(output_pdf_uni, width = 8.5, height = 11)

for (ep in endpoints) {
  message("Processing (univariable) endpoint: ", ep)
  fn <- sprintf("b1da_lnz_mvlreg_datasubset(%s).csv", ep)
  
  # 1) 데이터 로드
  df <- read_csv(fn, show_col_types = FALSE)
  
  # 2) 최소 전처리 (요청 흐름 준수)
  # - 문자 → factor
  # - PID factor
  # - ENDPOINT drop
  # - SEX ref = Male (0/1인 경우 라벨링 고정)
  if (!is.factor(df$SEX)) df$SEX <- factor(df$SEX)
  df$SEX <- factor(df$SEX, levels = c(0,1), labels = c("Male","Female"))
  df$SEX <- relevel(df$SEX, ref = "Male")
  
  df <- df |>
    select(-ENDPOINT) |>
    mutate(
      PID = as.factor(PID),
      across(where(is.character), as.factor)
    )
  
  if ("CUM_DOSE" %in% names(df)) {
    df <- df %>% mutate(CUM_DOSE = CUM_DOSE / 600)
  }
  
  # 3) 예측자 목록 (y, PID 제외)
  predictors <- setdiff(names(df), c("y", "PID"))
  
  # 4) 유니버리엇: 각 predictor마다 glm(y ~ predictor)
  uni_rows <- list()
  
  for (x in predictors) {
    fml <- as.formula(paste0("y ~ ", x))
    fit <- tryCatch(
      glm(fml, data = df, family = binomial(link = "logit")),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    
    # tidy로 계수/SE/p 추출 (Intercept 제외)
    td <- broom::tidy(fit) %>%
      filter(term != "(Intercept)")
    
    if (nrow(td) == 0) next
    
    # Wald CI + OR 변환
    td <- td %>%
      mutate(
        OR    = exp(estimate),
        LCL   = exp(estimate - 1.96*std.error),
        UCL   = exp(estimate + 1.96*std.error),
        Predictor = x
      ) %>%
      transmute(
        Predictor,
        Term = term,           # factor면 "변수명레벨" 형태로 나옴 (예: SEXFemale)
        OR,
        `2.5%` = LCL,
        `97.5%` = UCL,
        p.value = p.value
      )
    
    uni_rows[[x]] <- td
  }
  
  # endpoint 내 모든 predictor 결과 결합
  uni_tbl <- dplyr::bind_rows(uni_rows)
  
  # p.adjust는 endpoint 내의 모든 term에 대해 수행
  if (nrow(uni_tbl) > 0) {
    uni_tbl <- uni_tbl %>%
      mutate(
        adj.p.value = p.adjust(p.value, method = "BH"),
        sig = ifelse(adj.p.value < p_adj_cutoff, "*", "")
      )
  } else {
    # 빈 테이블일 경우 출력만 진행
    uni_tbl <- data.frame(
      Predictor = character(), Term = character(),
      OR = numeric(), `2.5%` = numeric(), `97.5%` = numeric(),
      p.value = numeric(), adj.p.value = numeric(), sig = character(),
      check.names = FALSE
    )
  }
  
  # ===== PDF 페이지 출력 (이 endpoint용) =====
  plot.new()
  title_line <- paste0("Univariable Logistic Regression (Endpoint: ", ep, ")")
  mtext(title_line, side = 3, adj = 0, line = 1, cex = 1.2, font = 2)
  mtext(paste0("File: ", fn), side = 3, adj = 0, line = 0, cex = 0.9)
  
  y <- 0.88
  line_h <- 0.02
  
  # 설명 텍스트
  txt_lines <- c(
    "Models: y ~ predictor (one predictor per model, binomial logit)",
    paste0("N = ", nrow(df), " (rows, after NA drop inside glm)"),
    "Displayed: OR with 95% Wald CI, raw p, BH-adjusted p, significance flag (adj.p < cutoff)"
  )
  for (ln in txt_lines) {
    text(0.02, y, ln, adj = c(0,1), cex = 0.85); y <- y - line_h
  }
  crit_line <- paste0("Criteria for significance page summary: (OR > ", OR_upper_cutoff,
                      " or OR < ", OR_lower_cutoff, ") AND adj.p.value < ", p_adj_cutoff)
  text(0.02, y, crit_line, adj = c(0,1), cex = 0.85); y <- y - (line_h + 0.01)
  
  # 테이블 헤더
  header <- c("Predictor","Term","OR","2.5%","97.5%","p.value","adj.p.value","sig")
  text(0.02, y, paste(sprintf("%-16s", header), collapse=" "),
       adj=c(0,1), cex=0.8, font=2)
  y <- y - 0.02
  
  if (nrow(uni_tbl) == 0) {
    text(0.02, y, "No rows to display (no fitted terms).", adj = c(0,1), cex = 0.85)
  } else {
    # 표시용 포맷
    fmt_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
    disp <- uni_tbl
    disp$OR      <- fmt_num(disp$OR)
    disp$`2.5%`  <- fmt_num(disp$`2.5%`)
    disp$`97.5%` <- fmt_num(disp$`97.5%`)
    disp$p.value <- ifelse(disp$p.value < 0.0001, "<0.0001", fmt_num(disp$p.value, 4))
    disp$`adj.p.value` <- ifelse(disp$`adj.p.value` < 0.0001, "<0.0001", fmt_num(disp$`adj.p.value`, 4))
    
    for (i in seq_len(nrow(disp))) {
      ln <- paste(sprintf("%-16s", disp[i, ]), collapse=" ")
      y <- y - line_h
      if (y < 0.05) { plot.new(); y <- 0.95 }
      text(0.02, y, ln, adj=c(0,1), cex=0.75)
    }
  }
  
  # ===== 기준 충족(row) → 전역 요약 테이블 누적 =====
  if (nrow(uni_tbl) > 0) {
    idx_sig <- (uni_tbl$OR > OR_upper_cutoff | uni_tbl$OR < OR_lower_cutoff) &
      (uni_tbl$adj.p.value < p_adj_cutoff)
    if (any(idx_sig, na.rm = TRUE)) {
      tmp <- uni_tbl[idx_sig, , drop = FALSE]
      tmp$Endpoint <- ep
      tmp <- tmp[, c("Endpoint","Predictor","Term","OR","2.5%","97.5%","p.value","adj.p.value","sig")]
      sig_or_table_uni <- rbind(sig_or_table_uni, tmp)
    }
  }
}

# ===== 마지막 페이지: Univariable Significant OR 요약 =====
plot.new()
title_line <- "Univariable Significant OR Summary across Endpoints"
mtext(title_line, side = 3, adj = 0, line = 1, cex = 1.2, font = 2)
criteria_line <- paste0("(Criteria: OR >", OR_upper_cutoff,
                        " or OR <", OR_lower_cutoff,
                        ") AND adj.p.value <", p_adj_cutoff, ")")
mtext(criteria_line, side = 3, adj = 0, line = 0, cex = 0.9)

y <- 0.88
line_h <- 0.02

if (nrow(sig_or_table_uni) == 0) {
  text(0.02, y, "No significant rows matched the criteria.", adj = c(0,1), cex = 0.9)
} else {
  fmt_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
  disp <- sig_or_table_uni
  disp$OR      <- fmt_num(disp$OR)
  disp$`2.5%`  <- fmt_num(disp$`2.5%`)
  disp$`97.5%` <- fmt_num(disp$`97.5%`)
  disp$p.value <- ifelse(disp$p.value < 0.0001, "<0.0001", fmt_num(disp$p.value, 4))
  disp$`adj.p.value` <- ifelse(disp$`adj.p.value` < 0.0001, "<0.0001", fmt_num(disp$`adj.p.value`, 4))
  
  header <- c("Endpoint","Predictor","Term","OR","2.5%","97.5%","p.value","adj.p.value","sig")
  text(0.02, y, paste(sprintf("%-14s", header), collapse=" "), adj=c(0,1), cex=0.8, font=2)
  y <- y - 0.02
  
  for (i in seq_len(nrow(disp))) {
    ln <- paste(sprintf("%-14s", disp[i, ]), collapse=" ")
    y <- y - line_h
    if (y < 0.05) { plot.new(); y <- 0.95 }
    text(0.02, y, ln, adj=c(0,1), cex=0.75)
  }
  
  # CSV 저장 (선택)
  try(write.csv(sig_or_table_uni, output_csv_uni, row.names = FALSE), silent = TRUE)
}

# ===== PDF 종료 =====
dev.off()
