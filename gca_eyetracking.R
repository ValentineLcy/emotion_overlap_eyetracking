library("lme4")
library("car")
library("lsmeans")
library("Rmisc")
library("ggplot2")
library("fitdistrplus")
library("dplyr")
library("interactions")

all_data<- read.csv('all_ages_face_dwell_time_last_half.csv')

# For each condition, within each subject , remove trials where
# DwellTimeDistractor is less than 2.5 SDs below the mean or
# 2.5 SDs above the mean
# Function to filter out trials based on 2.5 SDs within age, subject, and condition
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
    select(-mean_dwell, -sd_dwell, -lower_bound, -upper_bound)  # Remove intermediate columns
}

# Apply the filtering function
all_data <- filter_trials(all_data)

all_data <- all_data %>%
  mutate(Age = as.numeric(as.character(Age)))  # Ensure Age is numeric

#all_data <- all_data %>%
#mutate(EPDS = as.numeric(as.character(EPDS)))  # Ensure Age is numeric

# dummy variables to debug within the function
predictor <- "Age"
poly.order <- 2
orthogonal <- TRUE

code.poly <- function(df=NULL, predictor=NULL, poly.order=NULL, orthogonal=TRUE, draw.poly=TRUE){
  require(reshape2)
  require(ggplot2)
  
  # convert choice for orthogonal into choice for raw
  raw <- (orthogonal-1)^2
  
  # make sure that the declared predictor is actually present in the data.frame
  if (!predictor %in% names(df)){
    warning(paste0(predictor, " is not a variable in your data frame. Check spelling and try again"))
  }
  
  # Extract the vector to be used as the predictor
  predictor.vector <- df[,which(colnames(df)==predictor)]
  predictor.vector <- df[[predictor]] 
  
  # create index of predictor (e.g. numbered time bins)
  # the index of the time bin will be used later as an index to call the time sample
  predictor.indices <- as.numeric(as.factor(predictor.vector))
  
  df$temp.predictor.index <- predictor.indices
  
  #create x-order order polys (orthogonal if not raw)
  predictor.polynomial <- poly(x = unique(sort(predictor.vector)),
                               degree = poly.order, raw=raw)
  
  # use predictor index as index to align
  # polynomial-transformed predictor values with original dataset
  # (as many as called for by the polynomial order)
  df[, paste("poly", 1:poly.order, sep="")] <-
    predictor.polynomial[predictor.indices, 1:poly.order]
  
  # draw a plot of the polynomial transformations, if desired
  if (draw.poly == TRUE){
    # extract the polynomials from the df
    df.poly <- unique(df[c(predictor, paste("poly", 1:poly.order, sep=""))])
    
    # melt from wide to long format
    df.poly.melt <- melt(df.poly, id.vars=predictor)
    
    # Make level names intuitive
    # don't bother with anything above 6th order.
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly1"] <- "Linear"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly2"] <- "Quadratic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly3"] <- "Cubic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly4"] <- "Quartic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly5"] <- "Quintic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly6"] <- "Sextic"
    
    # change some column names for the output
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
  
  # restore correct column names
  colnames(df)[colnames(df) == "temp.predictor.index"] <- paste0(predictor,".Index")
  return(df)
}

data.gca <- code.poly(df=all_data, predictor="Age", poly.order=2, orthogonal=TRUE, draw.poly=TRUE)

lmer_model <- lmer(DwellTimeFace ~ Emotion*(poly1+poly2)+(1+Emotion|Participant), data = data.gca, REML = TRUE)
print(lmer_model)

lmer_results<-Anova(lmer_model, type = 3)
print(lmer_results)

subject_means <- all_data %>%
  ungroup() %>%
  group_by(Participant, Age, Emotion) %>%
  summarise(DwellTimeFace = mean(DwellTimeFace, na.rm = TRUE), .groups = "drop")

# Generate predictions from the fitted model
data.gca$Fitted <- predict(lmer_model, re.form = NA)

# Create the plot with actual fitted GCA model
ggplot(data = subject_means, aes(x = as.numeric(Age), y = DwellTimeFace, color=Emotion)) +
  geom_line(aes(group = Participant, color=Emotion), size = 0.5, alpha = 0.25) + # Faint grey lines for participants
  geom_point(position = position_jitter(width = 0.5, height = 0), size = 2, alpha = 0.25) +
  geom_line(data = data.gca, aes(x = as.numeric(Age), y = Fitted, color = Emotion), size = 1) +
  facet_grid(~Emotion) +
  scale_x_continuous(breaks = c(3, 6, 12), labels = c("3", "6", "12")) +
  labs(
    x = "Age (months)",
    y = "Mean latency to distractor"
  ) +
  theme_bw()