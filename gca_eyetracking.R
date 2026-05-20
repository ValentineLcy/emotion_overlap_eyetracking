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

log_file <- "./results/gca_eyetracking_output.txt"

dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_con <- file(log_file, open = "wt")

sink(log_con, split = TRUE)
sink(log_con, type = "message")

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

code.poly <- function(df=NULL, predictor=NULL, poly.order=NULL, orthogonal=TRUE, draw.poly=TRUE){
  require(reshape2)
  require(ggplot2)
  
  raw <- (orthogonal-1)^2
  
  if (!predictor %in% names(df)){
    warning(paste0(predictor, " is not a variable in your data frame. Check spelling and try again"))
  }
  
  predictor.vector <- df[,which(colnames(df)==predictor)]
  predictor.vector <- df[[predictor]] 
  
  predictor.indices <- as.numeric(as.factor(predictor.vector))
  
  df$temp.predictor.index <- predictor.indices
  
  predictor.polynomial <- poly(x = unique(sort(predictor.vector)),
                               degree = poly.order, raw=raw)
  
  df[, paste("poly", 1:poly.order, sep="")] <-
    predictor.polynomial[predictor.indices, 1:poly.order]
  
  if (draw.poly == TRUE){
    df.poly <- unique(df[c(predictor, paste("poly", 1:poly.order, sep=""))])
    df.poly.melt <- melt(df.poly, id.vars=predictor)
    
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly1"] <- "Linear"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly2"] <- "Quadratic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly3"] <- "Cubic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly4"] <- "Quartic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly5"] <- "Quintic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly6"] <- "Sextic"
    
    colnames(df.poly.melt)[colnames(df.poly.melt) == "variable"] <- "Order"
    
    poly.plot <- ggplot(df.poly.melt, aes(y=value, color=Order))+
      aes_string(x=predictor)+
      geom_line()+
      xlab(paste0(predictor, " (transformed polynomials)"))+
      ylab("Transformed value")+
      scale_color_brewer(palette="Set1")+
      theme_bw()
    
    print(poly.plot)
  }
  
  colnames(df)[colnames(df) == "temp.predictor.index"] <- paste0(predictor,".Index")
  return(df)
}

# FONCTION model fit ----
fit_mod3 <- function(predictor_time, data_input, motion_control=FALSE) {
  
  data_poly <- code.poly(
    df=data_input, 
    predictor=predictor_time, 
    poly.order=2, 
    orthogonal=TRUE, 
    draw.poly=FALSE
  )
  
  # beta transformation
  data_poly$y_star <- (data_poly$DwellTimeFace +1)/ 3001
  
  form_str='y_star ~ Emotion * (poly1 + poly2)'
  if(motion_control) {
    form_str = paste0(form_str, ' + MeanOpticalFlow')
  }
  form_str = paste0(form_str, ' + (1 | Participant)')
  mod <- glmmTMB(
    as.formula(form_str),
    family = beta_family(),
    data = data_poly
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

run_posthoc <- function(mod, data_poly, anova_results, file_suffix){
  # FOLLOW-UP ANALYSES - INTERACTIONS & MAIN EFFECTS
  # ---- EMmeans for emotion principal effect ----
  #pairs comparison between emotion all ages
  emotion_main_idx=which(rownames(anova_results)=='Emotion')
  if(anova_results$`Pr(>Chisq)`[emotion_main_idx]<0.05) {
    print('Main effect of emotion')
    emm_emotion <- emmeans(mod, ~ Emotion, type = "response")
    pairs_emotion <- pairs(emm_emotion, adjust = "tukey")
    print(pairs_emotion)
    print(emm_emotion)
  }
  
  # ---- EMtrends for interaction Emotion x poly1 ----
  # slopes pairs comparison of poly1 between emotion
  emotion_poly1_idx=which(rownames(anova_results)=='Emotion:poly1')
  if(anova_results$`Pr(>Chisq)`[emotion_poly1_idx]<0.05) {
    print('Emotion x poly1')
    emt_interaction_poly1 <- emtrends(mod, ~ Emotion, var = "poly1")
    pairs_trends_linear <- pairs(emt_interaction_poly1, adjust = "tukey")
    print(pairs_trends_linear)
  }
  
  # ---- EMtrends for interaction Emotion x poly2 ----
  # slopes pairs comparison of poly2 between emotion
  emotion_poly2_idx=which(rownames(anova_results)=='Emotion:poly2')
  if(anova_results$`Pr(>Chisq)`[emotion_poly2_idx]<0.05) {
    print('Emotion x poly2')
    emt_interaction_poly2 <- emtrends(mod, ~ Emotion, var = "poly2")
    pairs_trends_quad <- pairs(emt_interaction_poly2, adjust = "tukey")
    print(pairs_trends_quad)
  }
}


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
      .groups       = "drop"
    )
  
  # Prédictions pour 3 niveaux EPDS
  pred_3levels <- expand.grid(
    tmp        = unique(data_poly[[predictor_time]]),
    Emotion    = unique(data_poly$Emotion)
  ) %>%
    dplyr::rename(!!predictor_time := tmp)
  
  pred_3levels <- code.poly(
    df         = pred_3levels,
    predictor  = predictor_time,
    poly.order = 2,
    orthogonal = TRUE,
    draw.poly  = FALSE
  )
  if ("MeanOpticalFlow" %in% names(data_poly)) {
    pred_3levels$MeanOpticalFlow <- mean(data_poly$MeanOpticalFlow, na.rm = TRUE)
  }
  
  pred_3levels$Fitted <- predict(mod,
                                 newdata = pred_3levels,
                                 re.form = NA,
                                 type    = "response")
  
  dev.new(width = 3.5, height = 2)
  g<-ggplot(data = subject_means, aes(x = as.numeric(.data[[predictor_time]]), y = DwellTimeFace, color=Emotion)) +
    geom_line(aes(group = Participant, color=Emotion), size = 0.5, alpha = 0.25) +
    geom_point(position = position_jitter(width = 0.5, height = 0), size = 2, alpha = 0.25) +
    geom_line(data = pred_3levels, aes(x = as.numeric(.data[[predictor_time]]), y = Fitted, color = Emotion), size = 1) +
    facet_grid(~Emotion) +
    labs(
      x = "Age (days)",
      y = "Dwell time (proportion)"
    ) +
    theme_bw()
  print(g)
  out_fname=paste0('./figures/eyetracking/eyetracking_gca_',file_suffix,'.pdf')
  ggsave(
    out_fname,
    g,
    width = 6.5,
    height = 3.5,
    units = "in",
    dpi = 300
  )
  
  return(invisible(list(plot          = g,
                        subject_means = subject_means,
                        pred_3levels  = pred_3levels)))
}

time_periods<-c('first_half', 'last_half')

for(time_period in time_periods) {
  print(paste0('********************** ', time_period, ' ***************************'))
  all_data <- read.csv(paste0('./data/all_ages_face_dwell_time_', time_period, '.csv'))
  ages <- read.csv('./data/age_in_days.csv', sep=';')
  motion <- read.csv('./data/video_motion_brightness_contrast.csv')
  
  all_data <- all_data %>%
    left_join(ages %>% dplyr::select(Participant, Age_in_days_1st, Age_in_days_2nd, Age_in_days_3rd),
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
  
  data_raw<-all_data
  data_filtered <- filter_trials(all_data)
  
  cat("\n========== MODEL 1 : Age (months) — raw data ==========\n")
  file_suffix<-paste0('age_in_months_',time_period,'_unfiltered')
  res1  <- fit_mod3("Age", data_raw)
  mod1  <- res1$model
  data1 <- res1$data
  
  print(summary(mod1))
  anova_results<-Anova(mod1, type = 3)
  print(anova_results)
  run_diagnostics(mod1)
  run_posthoc(mod1, data1, anova_results, file_suffix)
  plot_mod3(mod1, data1,
            predictor_time = "Age",
            x_label        = "Age (months)",
            file_suffix
  )
  
  # MODEL 2 : Age (months) — filtered data ----
  cat("\n========== MODELE 2 : Age (months) — filtered ==========\n")
  file_suffix<-paste0('age_in_months_',time_period)
  res2  <- fit_mod3("Age", data_filtered)
  mod2  <- res2$model
  data2 <- res2$data
  
  print(summary(mod2))
  anova_results<-Anova(mod2, type = 3)
  print(anova_results)
  run_diagnostics(mod2)
  run_posthoc(mod2, data2, anova_results, file_suffix)
  plot_mod3(mod2, data2,
            predictor_time = "Age",
            x_label        = "Age (months)",
            file_suffix
  )
  
  # MODÈLE 3 : Age_in_days — raw data ----
  cat("\n========== MODÈLE 3 : Age_in_days — brut ==========\n")
  file_suffix<-paste0('age_in_days_',time_period,'_unfiltered')
  res3  <- fit_mod3("Age_in_days", data_raw)
  mod3  <- res3$model
  data3 <- res3$data
  
  print(summary(mod3))
  anova_results<-Anova(mod3, type = 3)
  print(anova_results)
  run_diagnostics(mod3)
  run_posthoc(mod3, data3, anova_results, file_suffix)
  plot_mod3(mod3, data3,
            predictor_time = "Age_in_days",
            x_label        = "Age (days)",
            file_suffix
  )
  
  
  # MODEL 4 : Age_in_days — filtered data ----
  cat("\n========== MODÈLE 4 : Age_in_days — filtered ==========\n")
  file_suffix<-paste0('age_in_days_',time_period)
  res4  <- fit_mod3("Age_in_days", data_filtered)
  mod4  <- res4$model
  data4 <- res4$data
  
  print(summary(mod4))
  anova_results<-Anova(mod4, type = 3)
  print(anova_results)
  run_diagnostics(mod4)
  run_posthoc(mod4, data4, anova_results, file_suffix)
  plot_mod3(mod4, data4,
            predictor_time = "Age_in_days",
            x_label        = "Age (days)",
            file_suffix
  )
  
  
  
  cat("\n========== MODEL 5 : Age (months) — raw data with motion control==========\n")
  file_suffix<-paste0('age_in_months_',time_period,'_motion_control_unfiltered')
  res1  <- fit_mod3("Age", data_raw, motion_control=TRUE)
  mod1  <- res1$model
  data1 <- res1$data
  
  print(summary(mod1))
  anova_results<-Anova(mod1, type = 3)
  print(anova_results)
  run_diagnostics(mod1)
  run_posthoc(mod1, data1, anova_results, file_suffix)
  plot_mod3(mod1, data1,
            predictor_time = "Age",
            x_label        = "Age (months)",
            file_suffix
  )
  
  # MODEL 2 : Age (months) — filtered data ----
  cat("\n========== MODELE 6 : Age (months) — filtered with motion control==========\n")
  file_suffix<-paste0('age_in_months_',time_period,'_motion_control')
  res2  <- fit_mod3("Age", data_filtered, motion_control=TRUE)
  mod2  <- res2$model
  data2 <- res2$data
  
  print(summary(mod2))
  anova_results<-Anova(mod2, type = 3)
  print(anova_results)
  run_diagnostics(mod2)
  run_posthoc(mod2, data2, anova_results, file_suffix)
  plot_mod3(mod2, data2,
            predictor_time = "Age",
            x_label        = "Age (months)",
            file_suffix
  )
  
  # MODÈLE 3 : Age_in_days — raw data ----
  cat("\n========== MODÈLE 7 : Age_in_days — brut with motion control==========\n")
  file_suffix<-paste0('age_in_days_',time_period,'_motion_control_unfiltered')
  res3  <- fit_mod3("Age_in_days", data_raw, motion_control=TRUE)
  mod3  <- res3$model
  data3 <- res3$data
  
  print(summary(mod3))
  anova_results<-Anova(mod3, type = 3)
  print(anova_results)
  run_diagnostics(mod3)
  run_posthoc(mod3, data3, anova_results, file_suffix)
  plot_mod3(mod3, data3,
            predictor_time = "Age_in_days",
            x_label        = "Age (days)",
            file_suffix
  )
  
  
  # MODEL 4 : Age_in_days — filtered data ----
  cat("\n========== MODÈLE 8 : Age_in_days — filtered with motion control ==========\n")
  file_suffix<-paste0('age_in_days_',time_period,'_motion_control')
  res4  <- fit_mod3("Age_in_days", data_filtered, motion_control=TRUE)
  mod4  <- res4$model
  data4 <- res4$data
  
  print(summary(mod4))
  anova_results<-Anova(mod4, type = 3)
  print(anova_results)
  run_diagnostics(mod4)
  run_posthoc(mod4, data4, anova_results, file_suffix)
  plot_mod3(mod4, data4,
            predictor_time = "Age_in_days",
            x_label        = "Age (days)",
            file_suffix
  )
  
  graphics.off()
}

sink(type = "message")
sink()
close(log_con)

