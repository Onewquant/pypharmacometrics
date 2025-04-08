# install.packages("dplyr")

library(ggplot2)
library(dplyr)


setwd("C:/Users/ilma0/PycharmProjects/pypharmacometrics/resource/KSCPTSPRWS25")
input_dir_path <- file.path(getwd(), "sglt2i_dataset")
results_dir_path <- file.path(getwd(), "results_r")


gdf <- as.data.frame(readr::read_csv(file.path(input_dir_path, "KSCPTSPRWS25_SGLT2i_PK.csv")))


# 개별 그래프 ----------------------------------------------------

# 색상 팔레트 설정
grp_levels <- sort(unique(gdf$GRP))
palette <- scales::hue_pal()(length(grp_levels))

# 개별 시간-농도 그래프
p_indiv <- ggplot(gdf, aes(x = ATIME, y = CONC, group = ID, color = factor(GRP))) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = palette, name = "GRP", labels = paste0("GRP ", grp_levels)) +
  labs(
    title = "[WSCT] Individual Time-Concentration Profiles by GRP",
    x = "Time (ATIME)",
    y = "Concentration (CONC)"
  ) +
  theme_bw(base_size = 15) +
  theme(legend.position = "right")

ggsave(filename = file.path(results_dir_path, "[WSCT] Individual Time-Concentration Profiles by GRP.png"), plot = p_indiv, width = 15, height = 12, dpi = 300)

# 그룹 평균 그래프 ----------------------------------------------------

summary_df <- gdf %>%
  group_by(GRP, NTIME) %>%
  summarise(
    mean_CONC = mean(CONC, na.rm = TRUE),
    sd_CONC = sd(CONC, na.rm = TRUE),
    .groups = "drop"
  )

p_pop <- ggplot(summary_df, aes(x = NTIME, y = mean_CONC, color = factor(GRP), group = GRP)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = mean_CONC - sd_CONC, ymax = mean_CONC + sd_CONC, fill = factor(GRP)), alpha = 0.2, color = NA) +
  scale_color_manual(values = palette, name = "GRP", labels = paste0("GRP ", grp_levels)) +
  scale_fill_manual(values = palette, guide = "none") +
  labs(
    title = "[WSCT] Population Time-Concentration Profiles by GRP",
    x = "Time (NTIME)",
    y = "Concentration (CONC)"
  ) +
  theme_bw(base_size = 15) +
  theme(legend.position = "right")

ggsave(filename = file.path(results_dir_path, "[WSCT] Population Time-Concentration Profiles by GRP.png"),
       plot = p_pop, width = 15, height = 12, dpi = 300)


# 비구획 분석(NCA) ----------------------------------------------------

library(readr)
library(NonCompart)

nca_result <- tblNCA(gdf, key = c("ID","GRP"), colTime = "ATIME", colConc = "CONC", dose=0.5, concUnit="ug/L", down = "Log")
nca_result <- nca_result %>% select(ID, GRP, TMAX, CMAX, AUCLST) %>% rename(ID = ID, GRP = GRP, Tmax = TMAX , Cmax = CMAX, AUClast = AUCLST )
write_csv(nca_result, file.path(results_dir_path, "[WSCT] NCARes_PK.csv"))

# ANOVA 및 GMR ----------------------------------------------------

nca_result <- nca_result %>%
  mutate(GRP = as.factor(GRP))

params <- c("AUClast", "Cmax")
results_list <- list()

for (param in params) {
  log_col <- paste0("log_", param)
  nca_result[[log_col]] <- log(as.numeric(nca_result[[param]]))
  
  model <- aov(as.formula(paste(log_col, "~ GRP")), data = nca_result)
  pval <- summary(model)[[1]][["Pr(>F)"]][1]
  
  coef_summary <- summary.lm(model)$coefficients
  df_resid <- df.residual(model)
  t_val <- qt(0.975, df_resid)
  
  summary_text <- c("GRP 1: reference")
  for (grp in levels(nca_result$GRP)[-1]) {
    coef_name <- paste0("GRP", grp)
    log_gmr <- coef_summary[coef_name, "Estimate"]
    se <- coef_summary[coef_name, "Std. Error"]
    ci_lower <- log_gmr - t_val * se
    ci_upper <- log_gmr + t_val * se
    
    gmr <- exp(log_gmr)
    gmr_ci_lower <- exp(ci_lower)
    gmr_ci_upper <- exp(ci_upper)
    
    summary_text <- c(summary_text,
                      sprintf("GRP %s: %.3f (%.3f-%.3f)", grp, gmr, gmr_ci_lower, gmr_ci_upper))
  }
  
  results_list[[param]] <- data.frame(
    Parameter = param,
    ANOVA_p_value = pval,
    GMR_95CI = paste(summary_text, collapse = "; ")
  )
}

anova_result_df <- bind_rows(results_list)
write.csv(anova_result_df, file.path(results_dir_path, "[WSCT] ANOVA_PK.csv"))