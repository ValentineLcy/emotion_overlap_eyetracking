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

filter<-TRUE

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
    
    #print(poly.plot)
  }
  
  colnames(df)[colnames(df) == "temp.predictor.index"] <- paste0(predictor,".Index")
  return(df)
}


time_periods<-c('first_half', 'last_half')

for(time_period in time_periods) {
  all_data <- read.csv(paste0('./data/all_ages_face_dwell_time_', time_period, '.csv'))
  ages <- read.csv('./data/age_in_days.csv', sep=';')
  
  all_data <- all_data %>%
    left_join(ages %>% dplyr::select(Participant, Age_in_days_1st, Age_in_days_2nd, Age_in_days_3rd),
              by = "Participant") %>%
    mutate(Age_in_days = case_when(
      Age == 3  ~ Age_in_days_1st,
      Age == 6  ~ Age_in_days_2nd,
      Age == 12 ~ Age_in_days_3rd,
      TRUE ~ NA_real_
    )) %>%
    dplyr::select(-Age_in_days_1st, -Age_in_days_2nd, -Age_in_days_3rd)
  
  all_data <- all_data %>%
    mutate(Age = as.numeric(as.character(Age)))
  
  if(filter){
    all_data <- filter_trials(all_data)
  }
  
  predictor <- "Age_in_days"
  poly.order <- 2
  orthogonal <- TRUE
  
  data.gca.beta <- code.poly(df=all_data, predictor="Age_in_days", poly.order=2, orthogonal=TRUE, draw.poly=TRUE)
  
  # beta transformation
  data.gca.beta$y_star <- (data.gca.beta$DwellTimeFace +1)/ 3001
  
  # Modèle glmmTMB avec distribution Beta
  #glmmTMB_model <- glmmTMB(y_star ~ Emotion*(poly1+poly2)+(1+Emotion|Participant),
                           #family = beta_family(link = "logit"),
                           #data = data.gca.beta)
  #print(glmmTMB_model)
  #glmmTMB_results <- Anova(glmmTMB_model, type = 3)
  #print(glmmTMB_results)
  
  mod3 <- glmmTMB(y_star ~ Emotion * (poly1 + poly2) + 
                    (1 | Participant),
                  family = beta_family(),
                  data = data.gca.beta)
  
  print(summary(mod3))
  mod3_results <- Anova(mod3, type = 3)
  print(mod3_results)
  
  subject_means <- all_data %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Participant, Age_in_days, Emotion) %>%
    dplyr::summarise(DwellTimeFace = mean(DwellTimeFace, na.rm = TRUE), .groups = "drop")
  cat("Colonnes :", names(subject_means), "\n")
  
  
  data.gca.beta$Fitted <- predict(mod3, re.form = NA, type = "response") * 3001
  
  dev.new(width = 3.5, height = 2)
  g<-ggplot(data = subject_means, aes(x = as.numeric(Age_in_days), y = DwellTimeFace, color=Emotion)) +
    geom_line(aes(group = Participant, color=Emotion), size = 0.5, alpha = 0.25) +
    geom_point(position = position_jitter(width = 0.5, height = 0), size = 2, alpha = 0.25) +
    geom_line(data = data.gca.beta, aes(x = as.numeric(Age_in_days), y = Fitted, color = Emotion), size = 1) +
    facet_grid(~Emotion) +
    labs(
      x = "Age (days)",
      y = "Mean dwell time to face"
    ) +
    theme_bw()
  print(g)
  out_fname=paste0('./figures/eyetracking_gca_age_in_days_',time_period)
  if(!filter){
    out_fname=paste0(out_fname,'_unfiltered')
  }
  ggsave(
    paste0(out_fname,'.pdf'),
    g,
    width = 6.5,
    height = 3.5,
    units = "in",
    dpi = 300
  )
  
  
  # FOLLOW-UP ANALYSES - INTERACTIONS & MAIN EFFECTS
  # ---- EMmeans for emotion principal effect ----
  #pairs comparison between emotion all ages
  if(mod3_results$`Pr(>Chisq)`[2]<0.05) {
    print('Main effect of emotion')
    emm_emotion <- emmeans(mod3, ~ Emotion, type = "response")
    pairs_emotion <- pairs(emm_emotion, adjust = "tukey")
    print(pairs_emotion)
    print(emm_emotion)
  }
  
  # ---- EMtrends for interaction Emotion x poly1 ----
  # slopes pairs comparison of poly1 between emotion
  if(mod3_results$`Pr(>Chisq)`[5]<0.05) {
    print('Emotion x poly1')
    emt_interaction_poly1 <- emtrends(mod3, ~ Emotion, var = "poly1")
    pairs_trends_linear <- pairs(emt_interaction_poly1, adjust = "tukey")
    print(pairs_trends_linear)
  }
  
  # ---- EMtrends for interaction Emotion x poly2 ----
  # slopes pairs comparison of poly2 between emotion
  if(mod3_results$`Pr(>Chisq)`[6]<0.05) {
    print('Emotion x poly2')
    emt_interaction_poly2 <- emtrends(mod3, ~ Emotion, var = "poly2")
    pairs_trends_quad <- pairs(emt_interaction_poly2, adjust = "tukey")
    print(pairs_trends_quad)
  }
}