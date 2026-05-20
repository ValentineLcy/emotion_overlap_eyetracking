library("glmmTMB")
library("car")
library("lsmeans")
library("Rmisc")
library("ggplot2")
library("fitdistrplus")
library("dplyr")
library("interactions")
library("emmeans")
library('DHARMa')
library('interactions')
library("reshape2")

# outliers filters (±2.5 SD par participant × Age × Emotion) ----
filter_trials <- function(data) {
  data %>%
    group_by(Participant, Age, Emotion) %>%
    mutate(
      mean_dwell  = mean(DwellTimeFace, na.rm = TRUE),
      sd_dwell    = sd(DwellTimeFace,   na.rm = TRUE),
      lower_bound = mean_dwell - 2.5 * sd_dwell,
      upper_bound = mean_dwell + 2.5 * sd_dwell
    ) %>%
    filter(
      DwellTimeFace >= lower_bound &
        DwellTimeFace <= upper_bound
    ) %>%
    dplyr::select(-mean_dwell, -sd_dwell, -lower_bound, -upper_bound) %>%
    ungroup()
}


# FONCTION code.poly ----
code.poly <- function(df=NULL, predictor=NULL, poly.order=NULL,
                      orthogonal=TRUE, draw.poly=TRUE){
  require(reshape2)
  require(ggplot2)
  
  raw <- (orthogonal-1)^2
  
  if (!predictor %in% names(df)){
    warning(paste0(predictor, " is not a variable in your data frame.
                   Check spelling and try again"))
  }
  
  predictor.vector  <- df[[predictor]]
  predictor.indices <- as.numeric(as.factor(predictor.vector))
  df$temp.predictor.index <- predictor.indices
  
  predictor.polynomial <- poly(x = unique(sort(predictor.vector)),
                               degree = poly.order, raw = raw)
  
  df[, paste("poly", 1:poly.order, sep="")] <-
    predictor.polynomial[predictor.indices, 1:poly.order]
  
  if (draw.poly == TRUE){
    df.poly      <- unique(df[c(predictor, paste("poly", 1:poly.order, sep=""))])
    df.poly.melt <- melt(df.poly, id.vars = predictor)
    
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly1"] <- "Linear"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly2"] <- "Quadratic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly3"] <- "Cubic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly4"] <- "Quartic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly5"] <- "Quintic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly6"] <- "Sextic"
    
    colnames(df.poly.melt)[colnames(df.poly.melt) == "variable"] <- "Order"
    
    # poly.plot <- ggplot(df.poly.melt, aes(y = value, color = Order)) +
    #   aes_string(x = predictor) +
    #   geom_line() +
    #   xlab(paste0(predictor, " (transformed polynomials)")) +
    #   ylab("Transformed value") +
    #   scale_color_brewer(palette = "Set1") +
    #   theme_bw()
    
    #print(poly.plot)
  }
  
  colnames(df)[colnames(df) == "temp.predictor.index"] <- paste0(predictor, ".Index")
  return(df)
}


# FONCTION model fit ----
fit_mod3 <- function(predictor_time, data_input) {
  
  data_poly <- code.poly(
    df         = data_input,
    predictor  = predictor_time,
    poly.order = 2,
    orthogonal = TRUE,
    draw.poly  = TRUE
  )
  
  data_poly$y_star <- (data_poly$DwellTimeFace +1)/ 3001
  
  mod <- glmmTMB(
    y_star ~ Emotion * epds_total * (poly1 + poly2) + MeanOpticalFlow +
      (1 | Participant),
    family = beta_family(link = "logit"),
    data   = data_poly
  )
  
  return(list(model = mod, data = data_poly))
}


# FONCTION DHARMa diagnostics ----
run_diagnostics <- function(mod) {
  dev.new()
  sim <- simulateResiduals(fittedModel = mod, plot = TRUE, n = 1000)
  dev.new()
  print(testDispersion(sim))
  dev.new()
  print(testZeroInflation(sim))
}


# post-hocs ----
run_posthoc <- function(mod, data_poly, anova_results, file_suffix) {
  
  epds_levels <- c(
    mean(data_poly$epds_total, na.rm = TRUE) - sd(data_poly$epds_total, na.rm = TRUE),
    mean(data_poly$epds_total, na.rm = TRUE),
    mean(data_poly$epds_total, na.rm = TRUE) + sd(data_poly$epds_total, na.rm = TRUE)
  )
  
  epds_range <- seq(
    min(data_poly$epds_total, na.rm = TRUE),
    max(data_poly$epds_total, na.rm = TRUE),
    length.out = 100
  )
  
  if(anova_results$`Pr(>Chisq)`[13]<0.05) {
    print('3w p2')
  }
  # POST-HOC 1 - Emotion:epds_total:poly1
  cat("\n--- POST-HOC 1 : Emotion:epds_total:poly1 ---\n")
  if(anova_results$`Pr(>Chisq)`[12]<0.05) {
    emtrends_triple <- emtrends(mod, ~ Emotion | epds_total,
                                var = "poly1",
                                at  = list(epds_total = epds_levels))
    print(emtrends_triple)
    print(pairs(emtrends_triple, adjust = "tukey"))
    
    # JOHNSON-NEYMAN : pairs emotions
    cat("\n--- JOHNSON-NEYMAN : paires entre émotions ---\n")
    jn_trends <- emtrends(mod, ~ Emotion | epds_total,
                          var = "poly1",
                          at  = list(epds_total = epds_range))
    jn_df <- as.data.frame(pairs(jn_trends, adjust = "none"))
    
    jn_sig <- jn_df %>%
      mutate(contrast = as.character(contrast),
             sig      = p.value < 0.05) %>%
      dplyr::group_by(contrast) %>%
      dplyr::summarise(
        epds_min_sig = ifelse(any(sig), min(epds_total[sig]), NA),
        epds_max_sig = ifelse(any(sig), max(epds_total[sig]), NA),
        always_sig   = all(sig),
        never_sig    = !any(sig),
        .groups      = "drop"
      )
    print(jn_sig)
    
    dev.new()
    p_jn_pairs <- jn_df %>%
      mutate(contrast = as.character(contrast),
             sig      = p.value < 0.05) %>%
      ggplot(aes(x = epds_total, y = estimate)) +
      geom_ribbon(aes(ymin = estimate - 1.96 * SE,
                      ymax = estimate + 1.96 * SE,
                      fill = sig), alpha = 0.3) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey70"),
                        labels = c("TRUE" = "p < .05", "FALSE" = "p ≥ .05"),
                        name   = "") +
      facet_wrap(~ contrast, scales = "free_y") +
      labs(title    = "Johnson-Neyman : différences de pente poly1 entre émotions",
           subtitle = "Zone bleue = différence significative (p < .05)",
           x        = "Score EPDS",
           y        = "Différence de pente (poly1)") +
      theme_bw()
    print(p_jn_pairs)
    out_fname=paste0('./figures/epds_eyetracking_gca_',file_suffix,'_motion_control_emotion_x_epds_x_poly1_pairwise.pdf')
    ggsave(
      out_fname,
      p_jn_pairs,
      dpi = 300
    )
    
    # JOHNSON-NEYMAN : Emotion
    cat("\n--- JOHNSON-NEYMAN : par Emotion ---\n")
    jn_trends_byemotion <- emtrends(mod, ~ epds_total | Emotion,
                                    var = "poly1",
                                    at  = list(epds_total = epds_range))
    
    jn_df_byemotion <- as.data.frame(jn_trends_byemotion) %>%
      mutate(z     = poly1.trend / SE,
             p_val = 2 * (1 - pnorm(abs(z))),
             sig   = p_val < 0.05)
    
    jn_sig_byemotion <- jn_df_byemotion %>%
      dplyr::group_by(Emotion) %>%
      dplyr::summarise(
        epds_min_sig = ifelse(any(sig), min(epds_total[sig]), NA),
        epds_max_sig = ifelse(any(sig), max(epds_total[sig]), NA),
        always_sig   = all(sig),
        never_sig    = !any(sig),
        .groups      = "drop"
      )
    print(jn_sig_byemotion)
    
    dev.new()
    p_jn_emotion <- jn_df_byemotion %>%
      ggplot(aes(x = epds_total, y = poly1.trend)) +
      geom_ribbon(aes(ymin = poly1.trend - 1.96 * SE,
                      ymax = poly1.trend + 1.96 * SE,
                      fill = sig), alpha = 0.3) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_fill_manual(values = c("TRUE" = "deeppink2", "FALSE" = "grey70"),
                        labels = c("TRUE" = "p < .05", "FALSE" = "p ≥ .05"),
                        name   = "") +
      facet_wrap(~ Emotion) +
      labs(title    = "Johnson-Neyman : poly1 slope per Emotion by EPDS",
           subtitle = "Pink zone = slope significantly different from 0 (p < .05)",
           x        = "EPDS Score",
           y        = "poly1 Slope") +
      theme_bw()
    print(p_jn_emotion)
    out_fname=paste0('./figures/epds_eyetracking_gca_',file_suffix,'_motion_control_emotion_x_epds_x_poly1.pdf')
    ggsave(
      out_fname,
      p_jn_emotion,
      dpi = 300
    )
  }
  
  # POST-HOC 2 - epds_total:poly1
  if(anova_results$`Pr(>Chisq)`[10]<0.05){
    cat("\n--- POST-HOC 2 : epds_total:poly1 ---\n")
    emtrends_poly1 <- emtrends(mod, ~ epds_total, var = "poly1",
                               at = list(epds_total = epds_levels))
    print(emtrends_poly1)
    print(pairs(emtrends_poly1, adjust = "tukey"))
  }
  
  # POST-HOC 3 - epds_total:poly2
  if(anova_results$`Pr(>Chisq)`[11]<0.05){
    cat("\n--- POST-HOC 3 : epds_total:poly2 ---\n")
    emtrends_poly2 <- emtrends(mod, ~ epds_total, var = "poly2",
                               at = list(epds_total = epds_levels))
    print(emtrends_poly2)
    print(pairs(emtrends_poly2, adjust = "tukey"))
  }
  
  # POST-HOC 4 - Emotion
  if(anova_results$`Pr(>Chisq)`[2]<0.05){
    cat("\n--- POST-HOC 4 : Emotion ---\n")
    emmeans_emotion <- emmeans(mod, ~ Emotion, type = "response")
    print(emmeans_emotion)
    print(pairs(emmeans_emotion, adjust = "tukey"))
  }
}


# plot ----
plot_mod3 <- function(mod, data_poly, predictor_time, x_label, file_suffix) {
  
  # Fitted values
  data_poly$Fitted <- NA
  rows_used <- as.integer(rownames(model.frame(mod)))
  data_poly$Fitted[rows_used] <- predict(mod, re.form = NA, type = "response")
  
  # Moyennes par sujet
  subject_means <- data_poly %>%
    dplyr::filter(!is.na(Fitted)) %>%
    dplyr::group_by(Participant, .data[[predictor_time]], Emotion) %>%
    dplyr::summarise(
      DwellTimeFace = mean(y_star,     na.rm = TRUE),
      Fitted        = mean(Fitted,     na.rm = TRUE),
      epds_total    = mean(epds_total, na.rm = TRUE),
      .groups       = "drop"
    )
  
  # Prédictions pour 3 niveaux EPDS
  pred_3levels <- expand.grid(
    tmp        = unique(data_poly[[predictor_time]]),
    Emotion    = unique(data_poly$Emotion),
    epds_total = c(0, 12.5, 25)
  ) %>%
    dplyr::rename(!!predictor_time := tmp) %>%
    mutate(epds_label = factor(epds_total,
                               levels = c(0, 12.5, 25),
                               labels = c("EPDS = 0", "EPDS = 12.5", "EPDS = 25")))
  
  pred_3levels <- code.poly(
    df         = pred_3levels,
    predictor  = predictor_time,
    poly.order = 2,
    orthogonal = TRUE,
    draw.poly  = FALSE
  )
  
  pred_3levels$Fitted <- predict(mod,
                                 newdata = pred_3levels,
                                 re.form = NA,
                                 type    = "response")
  
  dev.new(width = 3.5, height = 2)
  p <- ggplot() +
    geom_line(data = subject_means,
              aes(x     = .data[[predictor_time]],
                  y     = DwellTimeFace,
                  group = Participant,
                  color = epds_total),
              size = 0.5, alpha = 0.25) +
    geom_point(data = subject_means,
               aes(x     = .data[[predictor_time]],
                   y     = DwellTimeFace,
                   color = epds_total),
               position = position_jitter(width = 0.1, height = 0),
               size = 2, alpha = 0.25) +
    geom_line(data = pred_3levels,
              aes(x        = .data[[predictor_time]],
                  y        = Fitted,
                  color    = epds_total,
                  group    = epds_label,
                  linetype = epds_label),
              size = 1.2) +
    facet_grid(~ Emotion) +
    scale_x_continuous(breaks = unique(data_poly[[predictor_time]])) +
    scale_color_gradient(low = "blue", high = "red", name = "Score EPDS") +
    scale_linetype_manual(values = c("EPDS = 0"    = "solid",
                                     "EPDS = 12.5" = "solid",
                                     "EPDS = 25"   = "solid"),
                          name   = "EPDS level") +
    labs(title = paste0("DwellTime ~ EPDS * (poly1 + poly2) | ", predictor_time),
         x     = x_label,
         y     = "DwellTime (proportion)") +
    theme_bw()
  
  print(p)
  out_fname=paste0('./figures/epds_eyetracking_gca_',file_suffix,'_motion_control.pdf')
  ggsave(
    out_fname,
    p,
    width = 6.5,
    height = 3.5,
    units = "in",
    dpi = 300
  )
  return(invisible(list(plot          = p,
                        subject_means = subject_means,
                        pred_3levels  = pred_3levels)))
}

time_periods<-c('first_half', 'last_half')

for(time_period in time_periods) {
  print(paste0('********************** ', time_period, ' ***************************'))
  # DATA ----
  all_data <- read.csv(paste0('./data/all_ages_face_dwell_time_', time_period, '.csv'))
  
  motion <- read.csv('./data/video_motion_brightness_contrast.csv')
  #head(all_data)
  #head(motion)
  ages     <- read.csv('./data/age_in_days.csv', sep=';')
  #epds_3   <- read.csv('./data/3_EPDS.csv',  sep = ';')
  epds_3   <- read.csv('./data/3_EPDS.csv')
  #epds_6   <- read.csv('./data/6_EPDS.csv',  sep = ';')
  epds_6   <- read.csv('./data/6_EPDS.csv')
  #epds_12  <- read.csv('./data/12_EPDS.csv', sep = ';')
  epds_12  <- read.csv('./data/12_EPDS.csv')
  
  epds_3$epds_total  <- rowSums(epds_3[,  paste0("epds", 1:10)])
  epds_6$epds_total  <- rowSums(epds_6[,  paste0("epds", 1:10)])
  epds_12$epds_total <- rowSums(epds_12[, paste0("epds", 1:10)])
  
  epds_3$Age  <- 3
  epds_6$Age  <- 6
  epds_12$Age <- 12
  
  epds_3  <- epds_3  %>% dplyr::rename(Participant = Subject)
  epds_6  <- epds_6  %>% dplyr::rename(Participant = Subject)
  epds_12 <- epds_12 %>% dplyr::rename(Participant = Subject)
  
  epds_all <- bind_rows(
    epds_3  %>% dplyr::select(Participant, Age, epds_total),
    epds_6  %>% dplyr::select(Participant, Age, epds_total),
    epds_12 %>% dplyr::select(Participant, Age, epds_total)
  )
  
  all_data <- all_data %>%
    left_join(ages %>% dplyr::select(Participant, Age_in_days_1st,
                                     Age_in_days_2nd, Age_in_days_3rd),
              by = "Participant") %>%
    mutate(Age_in_days = case_when(
      Age == 3  ~ Age_in_days_1st,
      Age == 6  ~ Age_in_days_2nd,
      Age == 12 ~ Age_in_days_3rd,
      TRUE ~ NA_real_
    )) %>%
    dplyr::select(-Age_in_days_1st, -Age_in_days_2nd, -Age_in_days_3rd) %>%
    mutate(Age = as.numeric(as.character(Age)))
  
  all_data <- merge(all_data, motion, by = c("Actor", "Emotion"), all.x = TRUE)
  
  data_epds <- all_data %>%
    left_join(epds_all, by = c("Participant", "Age"), relationship='many-to-many')
  
  
  # filtered data ----
  data_epds_raw      <- data_epds
  data_epds_filtered <- filter_trials(data_epds)
  
  cat(sprintf(
    "Trials avant filtrage : %d | après filtrage : %d | retirés : %d (%.1f%%)\n",
    nrow(data_epds_raw),
    nrow(data_epds_filtered),
    nrow(data_epds_raw) - nrow(data_epds_filtered),
    100 * (1 - nrow(data_epds_filtered) / nrow(data_epds_raw))
  ))
  
  # MODEL 1 : Age (months) — raw data ----
  cat("\n========== MODEL 1 : Age (months) — raw data ==========\n")
  file_suffix<-paste0('age_in_months_',time_period,'_unfiltered')
  res1  <- fit_mod3("Age", data_epds_raw)
  mod1  <- res1$model
  data1 <- res1$data
  
  print(summary(mod1))
  anova_results<-Anova(mod1, type = 3)
  print(anova_results)
  run_diagnostics(mod1)
  run_posthoc(mod1, data1, anova_results, file_suffix)
  # plot_mod3(mod1, data1,
  #           predictor_time = "Age",
  #           x_label        = "Age (months)",
  #           file_suffix
  # )
  
  # MODEL 2 : Age (months) — filtered data ----
  cat("\n========== MODELE 2 : Age (months) — filtered ==========\n")
  file_suffix<-paste0('age_in_months_',time_period)
  res2  <- fit_mod3("Age", data_epds_filtered)
  mod2  <- res2$model
  data2 <- res2$data
  
  print(summary(mod2))
  anova_results<-Anova(mod2, type = 3)
  print(anova_results)
  run_diagnostics(mod2)
  run_posthoc(mod2, data2, anova_results, file_suffix)
  # plot_mod3(mod2, data2,
  #           predictor_time = "Age",
  #           x_label        = "Age (months)",
  #           file_suffix
  # )
  
  # MODÈLE 3 : Age_in_days — raw data ----
  cat("\n========== MODÈLE 3 : Age_in_days — brut ==========\n")
  file_suffix<-paste0('age_in_days_',time_period,'_unfiltered')
  res3  <- fit_mod3("Age_in_days", data_epds_raw)
  mod3  <- res3$model
  data3 <- res3$data
  
  print(summary(mod3))
  anova_results<-Anova(mod3, type = 3)
  print(anova_results)
  run_diagnostics(mod3)
  run_posthoc(mod3, data3, anova_results, file_suffix)
  # plot_mod3(mod3, data3,
  #           predictor_time = "Age_in_days",
  #           x_label        = "Age (days)",
  #           file_suffix
  # )
  
  
  # MODEL 4 : Age_in_days — filtered data ----
  cat("\n========== MODÈLE 4 : Age_in_days — filtered ==========\n")
  file_suffix<-paste0('age_in_days_',time_period)
  res4  <- fit_mod3("Age_in_days", data_epds_filtered)
  mod4  <- res4$model
  data4 <- res4$data
  
  print(summary(mod4))
  anova_results<-Anova(mod4, type = 3)
  print(anova_results)
  run_diagnostics(mod4)
  run_posthoc(mod4, data4, anova_results, file_suffix)
  # plot_mod3(mod4, data4,
  #           predictor_time = "Age_in_days",
  #           x_label        = "Age (days)",
  #           file_suffix
  # )
}