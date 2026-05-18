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

# DATA ----
all_data <- read.csv('all_ages_face_dwell_time_last_half.csv')
ages     <- read.csv('age_in_days.csv', sep=';')
bai_3   <- read.csv('3_bai.csv',  sep = ';')
bai_6   <- read.csv('6_bai.csv',  sep = ';')
bai_12  <- read.csv('12_bai.csv', sep = ';')
motion   <- read.csv('video_motion_brightness_contrast.csv')

bai_3$bai_total  <- rowSums(bai_3[,  paste0("bai", 1:10)])
bai_6$bai_total  <- rowSums(bai_6[,  paste0("bai", 1:10)])
bai_12$bai_total <- rowSums(bai_12[, paste0("bai", 1:10)])

bai_3$Age  <- 3
bai_6$Age  <- 6
bai_12$Age <- 12

bai_3  <- bai_3  %>% dplyr::rename(Participant = Subject)
bai_6  <- bai_6  %>% dplyr::rename(Participant = Subject)
bai_12 <- bai_12 %>% dplyr::rename(Participant = Subject)

bai_all <- bind_rows(
  bai_3  %>% dplyr::select(Participant, Age, bai_total),
  bai_6  %>% dplyr::select(Participant, Age, bai_total),
  bai_12 %>% dplyr::select(Participant, Age, bai_total)
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

data_bai <- all_data %>%
  left_join(bai_all, by = c("Participant", "Age"))

# filtering fonction
filter_trials <- function(data) {
  data %>%
    group_by(Participant, Age, Emotion) %>%  # Group by subject, age, and condition
    mutate(
      mean_dwell = mean(DwellTimeFace, na.rm = TRUE),
      sd_dwell = sd(DwellTimeFace, na.rm = TRUE),
      lower_bound = mean_dwell - 2.5 * sd_dwell,
      upper_bound = mean_dwell + 2.5 * sd_dwell
    ) %>%
    filter(
      DwellTimeFace >= lower_bound & 
        DwellTimeFace <= upper_bound
    ) %>%
    dplyr::select(-mean_dwell, -sd_dwell, -lower_bound, -upper_bound)  # Remove intermediate columns
}

data_bai_filtered <- filter_trials(data_bai)

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
    
    poly.plot <- ggplot(df.poly.melt, aes(y = value, color = Order)) +
      aes_string(x = predictor) +
      geom_line() +
      xlab(paste0(predictor, " (transformed polynomials)")) +
      ylab("Transformed value") +
      scale_color_brewer(palette = "Set1") +
      theme_bw()
    
    print(poly.plot)
  }
  
  colnames(df)[colnames(df) == "temp.predictor.index"] <- paste0(predictor, ".Index")
  return(df)
}

# FONCTION : fit base model ----
fit_mod_base <- function(predictor_time, data_input) {
  
  data_poly <- code.poly(
    df         = data_input,
    predictor  = predictor_time,
    poly.order = 2,
    orthogonal = TRUE,
    draw.poly  = FALSE
  )
  
  data_poly$y_star <- (data_poly$DwellTimeFace + 1) /3001
  
  mod <- glmmTMB(
    y_star ~ Emotion * bai_total * (poly1 + poly2) +
      (1 | Participant),
    family = beta_family(link = "logit"),
    data   = data_poly
  )
  
  return(list(model = mod, data = data_poly))
}

# FONCTION : fit motion model-----
fit_mod_control <- function(predictor_time, data_input) {
  
  data_poly <- code.poly(
    df         = data_input,
    predictor  = predictor_time,
    poly.order = 2,
    orthogonal = TRUE,
    draw.poly  = FALSE
  )
  
  data_poly$y_star <- (data_poly$DwellTimeFace +1) / 3001
  
  mod <- glmmTMB(
    y_star ~ Emotion * bai_total * (poly1 + poly2) + MeanOpticalFlow +
      (1 | Participant),
    family = beta_family(link = "logit"),
    data   = data_poly
  )
  
  return(list(model = mod, data = data_poly))
}

# FONCTION :DHARMa diagnostic ----
run_diagnostics <- function(mod) {
  sim <- simulateResiduals(fittedModel = mod, plot = TRUE, n = 1000)
  print(testDispersion(sim))
  print(testZeroInflation(sim))
}

# FONCTION : post-hocs ----
run_posthoc <- function(mod, data_poly) {
  
  bai_levels <- c(
    mean(data_poly$bai_total, na.rm = TRUE) - sd(data_poly$bai_total, na.rm = TRUE),
    mean(data_poly$bai_total, na.rm = TRUE),
    mean(data_poly$bai_total, na.rm = TRUE) + sd(data_poly$bai_total, na.rm = TRUE)
  )
  
  bai_range <- seq(
    min(data_poly$bai_total, na.rm = TRUE),
    max(data_poly$bai_total, na.rm = TRUE),
    length.out = 100
  )
  
  # bai_total:poly1
  cat("\n--- bai_total:poly1 ---\n")
  
  jn_poly1 <- emtrends(mod, ~ bai_total, var = "poly1",
                       at = list(bai_total = bai_range))
  
  jn_poly1_df <- as.data.frame(jn_poly1) %>%
    mutate(z     = poly1.trend / SE,
           p_val = 2 * (1 - pnorm(abs(z))),
           sig   = p_val < 0.05)
  
  jn_poly1_sig <- jn_poly1_df %>%
    dplyr::summarise(
      bai_min_sig = ifelse(any(sig), min(bai_total[sig]), NA),
      bai_max_sig = ifelse(any(sig), max(bai_total[sig]), NA),
      always_sig   = all(sig),
      never_sig    = !any(sig)
    )
  print(jn_poly1_sig)
  
  p_jn_poly1 <- jn_poly1_df %>%
    ggplot(aes(x = bai_total, y = poly1.trend)) +
    geom_ribbon(aes(ymin = poly1.trend - 1.96 * SE,
                    ymax = poly1.trend + 1.96 * SE,
                    fill = sig), alpha = 0.3) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("TRUE" = "deeppink2", "FALSE" = "grey70"),
                      labels = c("TRUE" = "p < .05", "FALSE" = "p ≥ .05"),
                      name   = "") +
    labs(title    = "Johnson-Neyman : poly1 slope per Emotion by BAI",
         subtitle = "Pink zone = slope significantly different from 0 (p < .05)",
         x        = "BAI score",
         y        = "poly1 Slope") +
    theme_bw()
  print(p_jn_poly1)
  
  #bai_total:poly2
  cat("\n--- bai_total:poly2 ---\n")
  
  jn_poly2 <- emtrends(mod, ~ bai_total, var = "poly2",
                       at = list(bai_total = bai_range))
  
  jn_poly2_df <- as.data.frame(jn_poly2) %>%
    mutate(z     = poly2.trend / SE,
           p_val = 2 * (1 - pnorm(abs(z))),
           sig   = p_val < 0.05)
  
  jn_poly2_sig <- jn_poly2_df %>%
    dplyr::summarise(
      bai_min_sig = ifelse(any(sig), min(bai_total[sig]), NA),
      bai_max_sig = ifelse(any(sig), max(bai_total[sig]), NA),
      always_sig   = all(sig),
      never_sig    = !any(sig)
    )
  print(jn_poly2_sig)
  
  p_jn_poly2 <- jn_poly2_df %>%
    ggplot(aes(x = bai_total, y = poly2.trend)) +
    geom_ribbon(aes(ymin = poly2.trend - 1.96 * SE,
                    ymax = poly2.trend + 1.96 * SE,
                    fill = sig), alpha = 0.3) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("TRUE" = "darkorange", "FALSE" = "grey70"),
                      labels = c("TRUE" = "p < .05", "FALSE" = "p ≥ .05"),
                      name   = "") +
    labs(title    = "Johnson-Neyman : poly2 slope by BAI score",
         subtitle = "Orange zone = slope significantly different from 0 (p < .05)",
         x        = "Score bai",
         y        = "Pente poly2") +
    theme_bw()
  print(p_jn_poly2)
  
  # Emotion
  cat("\n--- Emotion ---\n")
  emmeans_emotion <- emmeans(mod, ~ Emotion, type = "response")
  print(emmeans_emotion)
  print(pairs(emmeans_emotion, adjust = "tukey"))
}

# FONCTION : plot
plot_mod3 <- function(mod, data_poly, predictor_time, x_label) {
  
  # Fitted values
  data_poly$Fitted <- NA
  rows_used <- as.integer(rownames(model.frame(mod)))
  data_poly$Fitted[rows_used] <- predict(mod, re.form = NA, type = "response")
  
  # mean subject
  subject_means <- data_poly %>%
    dplyr::filter(!is.na(Fitted)) %>%
    dplyr::group_by(Participant, .data[[predictor_time]], Emotion) %>%
    dplyr::summarise(
      DwellTimeFace = mean(y_star,     na.rm = TRUE),
      Fitted        = mean(Fitted,     na.rm = TRUE),
      bai_total    = mean(bai_total, na.rm = TRUE),
      .groups       = "drop"
    )
  
  # Predictions for 3 BAI levels
  pred_3levels <- expand.grid(
    tmp        = unique(data_poly[[predictor_time]]),
    Emotion    = unique(data_poly$Emotion),
    bai_total = c(0, 31.5, 63)
  ) %>%
    dplyr::rename(!!predictor_time := tmp) %>%
    mutate(bai_label = factor(bai_total,
                               levels = c(0, 31.5, 63),
                               labels = c("bai = 0", "bai = 31.5", "bai = 63")))
  
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
  
  p <- ggplot() +
    geom_line(data = subject_means,
              aes(x     = .data[[predictor_time]],
                  y     = DwellTimeFace,
                  group = Participant,
                  color = bai_total),
              size = 0.5, alpha = 0.25) +
    geom_point(data = subject_means,
               aes(x     = .data[[predictor_time]],
                   y     = DwellTimeFace,
                   color = bai_total),
               position = position_jitter(width = 0.1, height = 0),
               size = 2, alpha = 0.25) +
    geom_line(data = pred_3levels,
              aes(x        = .data[[predictor_time]],
                  y        = Fitted,
                  color    = bai_total,
                  group    = bai_label,
                  linetype = bai_label),
              size = 1.2) +
    facet_grid(~ Emotion) +
    scale_x_continuous(breaks = unique(data_poly[[predictor_time]])) +
    scale_color_gradient(low = "blue", high = "red", name = "Score bai") +
    scale_linetype_manual(values = c("bai = 0"    = "solid",
                                     "bai = 31.5" = "solid",
                                     "bai = 63"   = "solid"),
                          name   = "bai level") +
    labs(title = paste0("DwellTime ~ bai * (poly1 + poly2) | ", predictor_time),
         x     = x_label,
         y     = "DwellTime (proportion)") +
    theme_bw()
  
  print(p)
  return(invisible(list(plot          = p,
                        subject_means = subject_means,
                        pred_3levels  = pred_3levels)))
}

# HELPER 
print_header <- function(title) {
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("  ", title, "\n")
  cat(strrep("=", 60), "\n\n")
}


# 4 combinations ----
configs <- list(
  #list(label       = "BASE MODEL — raw data",
  #fit_fn      = fit_mod_base,
  #data_input  = data_bai,
  #pred_time   = "Age",
  #x_label     = "Age (months)"),
  
  #list(label       = "BASE MODEL — filtering data",
  #fit_fn      = fit_mod_base,
  #data_input  = data_bai_filtered,
  #pred_time   = "Age",
  #x_label     = "Age (months)"),
  
  #list(label       = "CONTROL MODEL — raw data",
  #fit_fn      = fit_mod_control,
  #data_input  = data_bai,
  #pred_time   = "Age",
  #x_label     = "Age (months)"),
  
  list(label       = "CONTROL MODEL — filtering data",
       fit_fn      = fit_mod_control,
       data_input  = data_bai_filtered,
       pred_time   = "Age",
       x_label     = "Age (months)")
)

# inspect model ----
inspect_model <- function(mod, data_poly, predictor_time, x_label) {
  
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("  SUMMARY\n")
  cat(strrep("=", 60), "\n")
  print(summary(mod))
  
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("  ANOVA TYPE III\n")
  cat(strrep("=", 60), "\n")
  print(Anova(mod, type = 3))
  
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("  DIAGNOSTICS DHARMa\n")
  cat(strrep("=", 60), "\n")
  run_diagnostics(mod)
  
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("  POST-HOCS\n")
  cat(strrep("=", 60), "\n")
  run_posthoc(mod, data_poly)
  
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("  PLOT\n")
  cat(strrep("=", 60), "\n")
  plot_mod3(mod, data_poly,
            predictor_time = predictor_time,
            x_label        = x_label)
  
  invisible(NULL)
}

#EXECUTION ----
results_list <- list()

for (cfg in configs) {
  
  print_header(cfg$label)
  
  res  <- cfg$fit_fn(cfg$pred_time, cfg$data_input)
  mod  <- res$model
  data <- res$data
  
  inspect_model(mod, data, cfg$pred_time, cfg$x_label)
  
  results_list[[cfg$label]] <- list(model = mod, data = data)
}