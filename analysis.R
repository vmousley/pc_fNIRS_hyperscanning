#########################################################################
####### Parent-child fNIRS WTC hyperscanning analysis script ############
######### Written by Victoria Mousley (v.mousley@bbk.ac.uk) #############
#########################################################################

# Libraries -----
require('R.matlab') # for reading in structs
require('ggplot2') # plotting
require('lme4') # for linear mixed effects models
require('lmerTest') # for checking model fit
require('rstatix') # t-tests
require('dplyr') # managing data
require('tidyr') # managing data
require('ggpubr') # add stats to plots
require('report') # reporting model results
require('stringr') # for cleaning some excel column names

# Paths ---------

# change to match the locations of your data 
hbo <- readMat('/set/your/path/raw_HbO.mat')
hbr <- readMat('/set/your/path/raw_HbR.mat')
ssc_hbo <- readMat('/set/your/path/ssc_HbO.mat')
ssc_hbr <- readMat('/set/your/path/ssc_HbR.mat')

pda_hbo <- readMat('/set/your/path/pda_HbO.mat')
pda_hbr <- readMat('/set/your/path/pda_HbR.mat')
  
demo <- read.csv('/set/your/path/demo.csv')
questionnaire <- read.csv('/set/your/path/questionnaire.csv')
performance <- read.csv('/set/your/path/performance.csv')

# Unpacking struct function --------- 
# Unpack structs resulting from
# fNIRS hyperscanning pipeline: https://github.com/vmousley/pc_fNIRS_hyperscanning

unpack <- function(data, level, count) {
  
  # data = aggregated struct exported from hyperscanning pipeline linked above
  # level = 'channels' for channel-wise data and 'rois' for roi-wise data
  # count = number (as integer) of channels or regions
  
  # Extract allIDs and allType from hbo
  allIDs <- unlist(data$Dyads[,])
  allType <- as.character(unlist(data$Type[1,1]))
  
  chansDF <- as.data.frame(allIDs)
  chansDF$allType <- allType
  
  list <- list()

  # row = column in $Channel.ConditionMeans that corresponds to participant 
  # next two indices to get to data we want (organised by condition type)
  # within condition type (x), indexing first by trial number within the condition, 
  # then by two more indices necessary to handle the struct format. 
  # then get to the normal rows and columns that indicate the channel pairings 
  # rows = p1, columns = p2 
  
  if (level == 'channels'){
    subdata <- data$Channel.ConditionMeans
  } else if(level == 'rois'){
    subdata <- data$ROI.ConditionMeans
  } else{
    stop('Incorrect level value. Channel-wise data = `channels`. ROI-wise data = `rois`.')}
  
  for (row in 1:nrow(chansDF)){
    condition_count <- length(grep("$", names(subdata[1,row][[1]][[1]][,,1])))
    
    for (x in 1:condition_count){
      id <- chansDF$allIDs[row]
      type <- chansDF$allType[row]
      
      for (k in 1:length(subdata[1,row][[1]][[1]][,,1][x][[1]])){
        for (i in 1:count){
          for (j in 1:count){
            cell_value <- subdata[1,row][[1]][[1]][,,1][x][[1]][[k]][[1]][i,j]
            
            if (!is.na(cell_value)){
              pair <- paste0(i, "-", j) # written as channel of p1 - channel of p2 
              condition <- x
              trial <- k
              
              df <- data.frame(id, type, meanCoherence = cell_value, pair, condition, trial)
              list[[length(list)+1]] <- df
              
            }
          }
        }
      }
    }
  }
  chans <- do.call(rbind, list)
  return(chans)
}

# unpack channel-wise data
chansHbo <- unpack(data = hbo, level = 'channels', count = 18); dim(chansHbo)
chansHbr <- unpack(data = hbr, level = 'channels', count = 18); dim(chansHbr)
chansHboSSC <- unpack(data = ssc_hbo, level = 'channels', count = 18); dim(chansHboSSC)
chansHbrSSC <- unpack(data = ssc_hbr, level = 'channels', count = 18); dim(chansHbrSSC)

# unpack ROI-wise data
roisHbo <- unpack(data = hbo, level = 'rois', count = 4); dim(roisHbo)
roisHbr <- unpack(data = hbr, level = 'rois', count = 4); dim(roisHbr)
roisHboSSC <- unpack(data = ssc_hbo, level = 'rois', count = 4); dim(roisHboSSC)
roisHbrSSC <- unpack(data = ssc_hbr, level = 'rois', count = 4); dim(roisHbrSSC)

# eight data frames - four channel-wise and four ROI-wise.
# chromophores without SSC regression should have the same dimensions,
# and both chromophores with SSC regression should the same dimensions,
# but with and without SSC regression will be different from each other

# Exclusion criteria (channel-wise data) ----

# remove dyads with fewer than 2 trials per condition
# change this if you have different criteria
removeHbo <- chansHbo %>% 
  dplyr::group_by(id, condition) %>%
  dplyr::filter(max(trial) == 1) %>%
  dplyr::distinct(id); removeHbo

chansHbo <- chansHbo[!chansHbo$id %in% removeHbo$id, ]

removeHbr <- chansHbr %>% 
  dplyr::group_by(id, condition) %>%
  dplyr::filter(max(trial) == 1) %>%
  dplyr::distinct(id); removeHbr

chansHbr <- chansHbr[!chansHbr$id %in% removeHbr$id, ]

removeHboSSC <- chansHboSSC %>% 
  dplyr::group_by(id, condition) %>%
  dplyr::filter(max(trial) == 1) %>%
  dplyr::distinct(id); removeHboSSC

chansHboSSC <- chansHboSSC[!chansHboSSC$id %in% removeHboSSC$id, ]

removeHbrSSC <- chansHbrSSC %>% 
  dplyr::group_by(id, condition) %>%
  dplyr::filter(max(trial) == 1) %>%
  dplyr::distinct(id); removeHbrSSC

chansHbrSSC <- chansHbrSSC[!chansHbrSSC$id %in% removeHbrSSC$id, ]

# remove manually two dyads (in the case of this study) 
# for failure to meet channel-wise inclusion of > 1 valid ROI per participant

chansHbo <- chansHbo %>% dplyr::filter(id != 'C06M06' & id != 'C43M43')

# check that you have the number of participants you're expecting
length(unique(chansHbo$id))

chansHbr <- chansHbr %>% dplyr::filter(id != 'C06M06' & id != 'C43M43')
length(unique(chansHbr$id))

chansHboSSC <- chansHboSSC %>% dplyr::filter(id != 'C06M06' & id != 'C43M43')
length(unique(chansHboSSC$id))

chansHbrSSC <- chansHbrSSC %>% dplyr::filter(id != 'C06M06' & id != 'C43M43')
length(unique(chansHbrSSC$id))

# Aim 1: SSC impact ------------
# conduct paired t tests for dyads who contributed both raw & ssc_reg data
chansHbo$type = "No SSC Regression"
chansHboSSC$type = "With SSC Regression"

# should change condition labels based on your conditions
chansHbo$condition <- ifelse(chansHbo$condition == 1, "collab",
                                ifelse(chansHbo$condition == 2, "collabScreen",
                                       ifelse(chansHbo$condition == 3, "ind",
                                              ifelse(chansHbo$condition == 4, "rest", chansHbo$condition))))

chansHboSSC$condition <- ifelse(chansHboSSC$condition == 1, "collab",
                                ifelse(chansHboSSC$condition == 2, "collabScreen",
                                       ifelse(chansHboSSC$condition == 3, "ind",
                                              ifelse(chansHboSSC$condition == 4, "rest", chansHboSSC$condition))))

chansHbr$type = "No SSC Regression"
chansHbrSSC$type = "With SSC Regression"

chansHbr$condition <- ifelse(chansHbr$condition == 1, "collab",
                                ifelse(chansHbr$condition == 2, "collabScreen",
                                       ifelse(chansHbr$condition == 3, "ind",
                                              ifelse(chansHbr$condition == 4, "rest", chansHbr$condition))))

chansHbrSSC$condition <- ifelse(chansHbrSSC$condition == 1, "collab",
                                ifelse(chansHbrSSC$condition == 2, "collabScreen",
                                       ifelse(chansHbrSSC$condition == 3, "ind",
                                              ifelse(chansHbrSSC$condition == 4, "rest", chansHbrSSC$condition))))

ssc_impact_hbo <- rbind(chansHbo, chansHboSSC)
common_ids_hbo <- intersect(chansHbo$id, chansHboSSC$id)
ssc_impact_hbo <- ssc_impact_hbo %>% filter(id %in% common_ids_hbo)
table(ssc_impact_hbo$id, ssc_impact_hbo$type)

ssc_impact_hbr <- rbind(chansHbr, chansHbrSSC)
common_ids_hbr <- intersect(chansHbr$id, chansHbrSSC$id)
ssc_impact_hbr <- ssc_impact_hbr %>% filter(id %in% common_ids_hbr)
table(ssc_impact_hbr$id, ssc_impact_hbr$type)

# average across trial
ssc_impact_hbo_agg <- ssc_impact_hbo %>% dplyr::group_by(id, type, pair, condition) %>%
  dplyr::mutate(allTrialAv = mean(meanCoherence)) %>%
  select(-c('meanCoherence', 'trial')) %>% distinct()

ssc_impact_hbr_agg <- ssc_impact_hbr %>% dplyr::group_by(id, type, pair, condition) %>%
  dplyr::mutate(allTrialAv = mean(meanCoherence)) %>%
  select(-c('meanCoherence', 'trial')) %>% distinct()

# paired, one-sided, fdr-corrected welch's t-tests
ssc_impact_hbo_test <- ssc_impact_hbo_agg %>% 
  filter(condition != 'rest') %>%
  dplyr::group_by(condition) %>%
  t_test(allTrialAv ~ type, alternative = 'greater', paired = TRUE, detailed = TRUE)

ssc_impact_hbr_test <- ssc_impact_hbr_agg %>% 
  filter(condition != 'rest') %>%
  dplyr::group_by(condition) %>%
  t_test(allTrialAv ~ type, alternative = 'greater', paired = TRUE, detailed = TRUE)

ssc_impact_hbo_test <- adjust_pvalue(ssc_impact_hbo_test, p.col= 'p', method = 'fdr')
ssc_impact_hbr_test <- adjust_pvalue(ssc_impact_hbr_test, p.col= 'p', method = 'fdr')

ssc_impact_hbo_test <- ssc_impact_hbo_test %>% add_xy_position(x = 'condition')
ssc_impact_hbr_test <- ssc_impact_hbr_test %>% add_xy_position(x = 'condition')

# manually add stars for plots based on your results
ssc_impact_hbo_test$stars <- c('***','***','***')
ssc_impact_hbr_test$stars <- c('***','***','')

ssc_impact_hbo_plot <- ssc_impact_hbo_agg %>%
  dplyr::filter(condition != 'rest') %>%
  ggplot(aes(x = condition, y = allTrialAv)) +
  geom_boxplot(aes(fill = type)) +
  ylim(0.2,0.8)+
  scale_fill_manual(values = c('#993333','#F33000'))+
  stat_pvalue_manual(ssc_impact_hbo_test, label = 'stars', stacked = TRUE, hide.ns = TRUE, step.increase = 0.05)+
  labs(x = '', y = 'WTC', 
       title = 'Oxyhaemoglobin', 
       fill = 'Data Type') +
  scale_x_discrete(labels = c('Full Collaboration', 'CollaborationScreen', 'Individual'))+
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = 'bold', margin = margin(t = 10)),
        #plot.subtitle = element_text(hjust = 0.5, size = 12, face = 'bold', margin = margin(t = 0, b = 10)),
        axis.title.x = element_text(hjust = 0.5, size = 14, face = 'bold'),
        axis.text.x = element_text(hjust = 0.5, size = 14),
        axis.title.y = element_text(hjust = 0.5, size = 14, face = 'bold', margin = margin(l = 10)),
        axis.text.y = element_text(hjust = 0.5, size = 14),
        legend.position = 'right',
        legend.title = element_text(hjust= 0.5, size = 14, face = 'bold'),
        legend.text = element_text(size = 14),
        plot.margin = unit(c(0, 1, 0, 0), "cm")); ssc_impact_hbo_plot

ssc_impact_hbr_plot <- ssc_impact_hbr_agg %>%
  dplyr::filter(condition != 'rest') %>%
  ggplot(aes(x = condition, y = allTrialAv)) +
  geom_boxplot(aes(fill = type)) +
  ylim(0.2,0.8)+
  scale_fill_manual(values = c('#0066CC','#99CCFF'))+
  stat_pvalue_manual(ssc_impact_hbr_test, label = 'stars', stacked = TRUE, hide.ns = TRUE, step.increase = 0.05)+
  labs(x = '', y = 'WTC', 
       title = 'Deoxyhaemoglobin', 
       fill = 'Data Type') +
  scale_x_discrete(labels = c('Full Collaboration', 'CollaborationScreen', 'Individual'))+
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = 'bold', margin = margin(t = 10)),
        #plot.subtitle = element_text(hjust = 0.5, size = 12, face = 'bold', margin = margin(t = 0, b = 10)),
        axis.title.x = element_text(hjust = 0.5, size = 14, face = 'bold'),
        axis.text.x = element_text(hjust = 0.5, size = 14),
        axis.title.y = element_text(hjust = 0.5, size = 14, face = 'bold', margin = margin(l = 10)),
        axis.text.y = element_text(hjust = 0.5, size = 14),
        legend.position = 'right',
        legend.title = element_text(hjust= 0.5, size = 14, face = 'bold'),
        legend.text = element_text(size = 14),
        plot.margin = unit(c(0, 1, 0, 0), "cm")); ssc_impact_hbr_plot

ssc_impact_hbo_hbr_plots <- ggarrange(ssc_impact_hbo_plot,
                                      ssc_impact_hbr_plot,
                                 ncol=1,
                                 nrow=2,
                                 common.legend = FALSE); ssc_impact_hbo_hbr_plots

# change to match your preferred save location
ggsave("/set/your/path/ssc_impact.png", 
       plot = ssc_impact_hbo_hbr_plots,
       width = 10, height = 10)

hbo_ssc_impact_descriptives <- ssc_impact_hbo_agg %>%
  dplyr::group_by(type, condition) %>%
  dplyr::mutate(av = mean(allTrialAv),
                sd = sd(allTrialAv),
                min = min(allTrialAv),
                max = max(allTrialAv)) %>%
  dplyr::ungroup() %>%
  select(-c('id','pair','allTrialAv')) %>% 
  distinct()

hbr_ssc_impact_descriptives <- ssc_impact_hbr_agg %>%
  dplyr::group_by(type, condition) %>%
  dplyr::mutate(av = mean(allTrialAv),
                sd = sd(allTrialAv),
                min = min(allTrialAv),
                max = max(allTrialAv)) %>%
  dplyr::ungroup() %>%
  select(-c('id','pair','allTrialAv')) %>% 
  distinct()

## Create full (mixed) datasets for rest of channel-wise analysis ----
### hbo ----
raw_data_to_add_hbo <- anti_join(chansHbo, chansHboSSC, by = "id")

# should be the number of participants who have raw but not SSC data
length(unique(raw_data_to_add_hbo$id))
hbo_df <- rbind(chansHboSSC, raw_data_to_add_hbo)

# should be the number of total participants (with & without ssc)
length(unique(hbo_df$id))

# all participants should have EITHER raw or SSC data, not both
table(hbo_df$id, hbo_df$type)

### hbr-------
raw_data_to_add_hbr <- anti_join(chansHbr, chansHbrSSC, by = "id")

# should be the number of participants who have raw but not SSC data
length(unique(raw_data_to_add_hbr$id))
hbr_df <- rbind(chansHbrSSC, raw_data_to_add_hbr)

# should be the number of total participants (with & without ssc)
length(unique(hbr_df$id))

# all participants should have EITHER raw or SSC data, not both
table(hbr_df$id, hbr_df$type)

## full analysis datasets are hbo_df and hbr_df

### Demographics -------
hbo_df$age <- NA
hbr_df$age <- NA

names(demo)[names(demo)=='Collab.ID'] <- 'id'
demo2 <- demo %>% dplyr::select('id', 'childAge') %>% na.omit()

for (unique_id in unique(demo2$id)) {
  matching_rows_hbo <- grepl(unique_id, hbo_df$id)
  matching_rows_hbr <- grepl(unique_id, hbr_df$id)
  
  if (any(matching_rows_hbo)) {
    hbo_df$age[matching_rows_hbo] <- demo2$childAge[demo2$id == unique_id]
  }
  
  if (any(matching_rows_hbr)) {
    hbr_df$age[matching_rows_hbr] <- demo2$childAge[demo2$id == unique_id]
  }
}


# number total participants (should be the same in hbo and hbr)
length(unique(hbo_df$id))

ages <- hbo_df %>% select(id, age) %>%
  distinct() %>%
  mutate(av = mean(age),
         sd = sd(age),
         min = min(age),
         max = max(age)) %>%
  select(av, sd, min, max) %>%
  distinct()

# summary statistics for tables
bytrial_summary_hbo <- hbo_df %>% 
  dplyr::filter(condition != 'rest') %>%
  dplyr::group_by(condition, trial) %>%
  dplyr::mutate(wholeBrainCoh = round(mean(meanCoherence),2), 
                sd = round(sd(meanCoherence),2),
                min = round(min(meanCoherence),2),
                max = round(max(meanCoherence),2)) %>%
  dplyr::ungroup() %>% 
  dplyr::select(-pair, -meanCoherence, -type, -id, -age) %>%
  dplyr::distinct()

agg_summary_hbo <- hbo_df %>% dplyr::filter(condition != 'rest') %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(wholeBrainCoh = round(mean(meanCoherence),2), 
                sd = round(sd(meanCoherence),2),
                min = round(min(meanCoherence),2),
                max = round(max(meanCoherence),2)) %>%
  dplyr::ungroup() %>% 
  dplyr::select(-pair, -meanCoherence, -type, -id, -age, -trial) %>%
  dplyr::distinct()

bytrial_summary_hbr <- hbr_df %>% 
  dplyr::filter(condition != 'rest') %>%
  dplyr::group_by(condition, trial) %>%
  dplyr::mutate(wholeBrainCoh = round(mean(meanCoherence),2), 
                sd = round(sd(meanCoherence),2),
                min = round(min(meanCoherence),2),
                max = round(max(meanCoherence),2)) %>%
  dplyr::ungroup() %>% 
  dplyr::select(-pair, -meanCoherence, -type, -id, -age) %>%
  dplyr::distinct()

agg_summary_hbr <- hbr_df %>% dplyr::filter(condition != 'rest') %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(wholeBrainCoh = round(mean(meanCoherence),2), 
                sd = round(sd(meanCoherence),2),
                min = round(min(meanCoherence),2),
                max = round(max(meanCoherence),2)) %>%
  dplyr::ungroup() %>% 
  dplyr::select(-pair, -meanCoherence, -type, -id, -age, -trial) %>%
  dplyr::distinct()

# Aim 2: Main condition analyses  --------

# remove rest data for analysis
# if you are analysing rest data, you'll need to change this
hbo_df <- hbo_df %>% filter(condition != 'rest')
hbr_df <- hbr_df %>% filter(condition != 'rest')

hbo_df$trialReleveled <- as.factor(hbo_df$trial)
hbo_df$conditionReleveled <- as.factor(hbo_df$condition)
hbo_df$conditionReleveled <- relevel(hbo_df$conditionReleveled, 'ind')

hbr_df$trialReleveled <- as.factor(hbr_df$trial)
hbr_df$conditionReleveled <- as.factor(hbr_df$condition)
hbr_df$conditionReleveled <- relevel(hbr_df$conditionReleveled, 'ind')

## hbo ----- 
m0 <- lmer(meanCoherence ~ (1|id), data = hbo_df); summary(m0)

m1 <- lmer(meanCoherence ~ age + (1|id), data = hbo_df); summary(m1)
anova(m0,m1)

m2 <- lmer(meanCoherence ~ trialReleveled + (1|id), data = hbo_df); summary(m2)
anova(m0,m2)

m3 <- lmer(meanCoherence ~ conditionReleveled + (1|id), data = hbo_df); summary(m3)
anova(m0,m3)

summary(m3)
report(m3)

# these lines, which appear throughout, are checking model assumptions for lmers
# standardized residuals versus fitted values by condition
plot(m3, resid(., scaled=TRUE) ~ fitted(.) | conditionReleveled, abline = 0)

# box-plots of residuals by subject
plot(m3, id ~ resid(., scaled = TRUE))

# observed versus fitted values by subject
plot(m3, meanCoherence ~ fitted(.) | id, abline = c(0.1))
res_m3 <- residuals(m3)

qqnorm(res_m3)
hist(res_m3)
ggplot(hbo_df, aes(y = res_m3, x = condition)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  xlab("Predictor Value (condition)") +
  ylab("Residual") +
  ggtitle("Residuals vs. Predictor Plot")

channelPairs_t_tests_hbo <- hbo_df %>%
  dplyr::group_by(id, pair, condition) %>%
  dplyr::mutate(PerChannelCoherence = mean(meanCoherence)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c('trial', 'trialReleveled', 'meanCoherence', 'age')) %>%
  dplyr::distinct() %>%
  group_by(pair) %>%
  t_test(PerChannelCoherence ~ condition, alternative = 'greater',
         paired = TRUE, p.adjust.method = 'fdr') %>%
  filter(p.adj.signif != 'ns')

## hbr ----- 
m0 <- lmer(meanCoherence ~ (1|id), data = hbr_df); summary(m0)

m1 <- lmer(meanCoherence ~ age + (1|id), data = hbr_df); summary(m1)
anova(m0,m1)

m2 <- lmer(meanCoherence ~ trialReleveled + (1|id), data = hbr_df); summary(m2)
anova(m0,m2)

m3 <- lmer(meanCoherence ~ trialReleveled + conditionReleveled + (1|id), data = hbr_df); summary(m3)
anova(m2,m3)

summary(m3)
report(m3)

# standardized residuals versus fitted values by condition
plot(m3, resid(., scaled=TRUE) ~ fitted(.) | conditionReleveled, abline = 0)

# box-plots of residuals by subject
plot(m3, id ~ resid(., scaled = TRUE))

# observed versus fitted values by participant
plot(m3, meanCoherence ~ fitted(.) | id, abline = c(0.1))
res_m3 <- residuals(m3)

qqnorm(res_m3)
hist(res_m3)
ggplot(hbr_df, aes(y = res_m3, x = condition)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  xlab("Predictor Value (condition)") +
  ylab("Residual") +
  ggtitle("Residuals vs. Predictor Plot")

channelPairs_t_tests_hbr <- hbr_df %>%
  dplyr::group_by(id, pair, condition) %>%
  dplyr::mutate(PerChannelCoherence = mean(meanCoherence)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c('trial', 'trialReleveled', 'meanCoherence', 'age')) %>%
  dplyr::distinct() %>%
  group_by(pair) %>%
  t_test(PerChannelCoherence ~ condition, alternative = 'greater',
         paired = TRUE, p.adjust.method = 'fdr') %>%
  filter(p.adj.signif != 'ns')

### ROI-wise analysis -----
roisHbo$type = "No SSC Regression"
roisHboSSC$type = "With SSC Regression"

# should change condition labels based on your conditions
roisHbo$condition <- ifelse(roisHbo$condition == 1, "collab",
                             ifelse(roisHbo$condition == 2, "collabScreen",
                                    ifelse(roisHbo$condition == 3, "ind",
                                           ifelse(roisHbo$condition == 4, "rest", roisHbo$condition))))

roisHboSSC$condition <- ifelse(roisHboSSC$condition == 1, "collab",
                                ifelse(roisHboSSC$condition == 2, "collabScreen",
                                       ifelse(roisHboSSC$condition == 3, "ind",
                                              ifelse(roisHboSSC$condition == 4, "rest", roisHboSSC$condition))))

roisHbr$type = "No SSC Regression"
roisHbrSSC$type = "With SSC Regression"

roisHbr$condition <- ifelse(roisHbr$condition == 1, "collab",
                             ifelse(roisHbr$condition == 2, "collabScreen",
                                    ifelse(roisHbr$condition == 3, "ind",
                                           ifelse(roisHbr$condition == 4, "rest", roisHbr$condition))))

roisHbrSSC$condition <- ifelse(roisHbrSSC$condition == 1, "collab",
                            ifelse(roisHbrSSC$condition == 2, "collabScreen",
                                   ifelse(roisHbrSSC$condition == 3, "ind",
                                          ifelse(roisHbrSSC$condition == 4, "rest", roisHbrSSC$condition))))

raw_roi_data_to_add_hbo <- anti_join(roisHbo, roisHboSSC, by = "id")
raw_roi_data_to_add_hbr <- anti_join(roisHbr, roisHbrSSC, by = "id")

# should be the number of participants who have raw but not SSC data
length(unique(raw_roi_data_to_add_hbo$id))
hbo_roi_df <- rbind(roisHboSSC, raw_roi_data_to_add_hbo)

# should be the number of total participants (with & without ssc)
length(unique(hbo_roi_df$id))

# all participants should have EITHER raw or SSC data, not both
table(hbo_roi_df$id, hbo_roi_df$type)

# should be the number of participants who have raw but not SSC data
length(unique(raw_roi_data_to_add_hbr$id))
hbr_roi_df <- rbind(roisHbrSSC, raw_roi_data_to_add_hbr)

# should be the number of total participants (with & without ssc)
length(unique(hbr_roi_df$id))

# all participants should have EITHER raw or SSC data, not both
table(hbr_roi_df$id, hbr_roi_df$type)

for (unique_id in unique(demo2$id)) {
  matching_rows_roiHbo <- grepl(unique_id, hbo_roi_df$id)
  matching_rows_roiHbr <- grepl(unique_id, hbr_roi_df$id)
  
  if (any(matching_rows_roiHbo)) {
    hbo_roi_df$age[matching_rows_roiHbo] <- demo$childAge[demo$id == unique_id]
  }
  
  if (any(matching_rows_roiHbr)) {
    hbr_roi_df$age[matching_rows_roiHbr] <- demo$childAge[demo$id == unique_id]
  }

}

hbo_roi_df <- hbo_roi_df %>% dplyr::filter(condition != 'rest')
hbr_roi_df <- hbr_roi_df %>% dplyr::filter(condition != 'rest')

#### Exclusion criteria (ROI-wise data) ----

# 2 trials per condition
removeHboROI <- hbo_roi_df %>% 
  dplyr::group_by(id, condition) %>%
  dplyr::filter(max(trial) == 1) %>%
  dplyr::distinct(id); removeHboROI

hbo_roi_df <- hbo_roi_df[!hbo_roi_df$id %in% removeHboROI$id, ]

removeHbrROI <- hbr_roi_df %>% 
  dplyr::group_by(id, condition) %>%
  dplyr::filter(max(trial) == 1) %>%
  dplyr::distinct(id); removeHbrROI

hbr_roi_df <- hbr_roi_df[!hbr_roi_df$id %in% removeHbrROI$id, ]

# you should remove any manually here that you know need to be excluded
# in this case, removed for < 1 valid ROI per participant
hbo_roi_df <- hbo_roi_df %>% dplyr::filter(id != 'C06M06' & id != 'C43M43')
length(unique(hbo_roi_df$id))

hbr_roi_df <- hbr_roi_df %>% dplyr::filter(id != 'C06M06' & id != 'C43M43')
length(unique(chansHbr$id))

hbo_roi_df$conditionReleveled <- as.factor(hbo_roi_df$condition)
hbo_roi_df$conditionReleveled <- relevel(hbo_roi_df$conditionReleveled, 'ind')

hbr_roi_df$conditionReleveled <- as.factor(hbr_roi_df$condition)
hbr_roi_df$conditionReleveled <- relevel(hbr_roi_df$conditionReleveled, 'ind')

#### hbo ------
m0 <- lmer(meanCoherence ~ (1|id), data = hbo_roi_df); summary(m0)
m1 <- lmer(meanCoherence ~ trial + (1|id), data = hbo_roi_df); summary(m1)
anova(m0, m1)

m2 <- lmer(meanCoherence ~ age + (1|id), data = hbo_roi_df); summary(m2)
anova(m0,m2)

m3 <- lmer(meanCoherence ~ age + conditionReleveled + (1|id), data = hbo_roi_df); summary(m3)
anova(m2,m3)

summary(m2)
report(m2)

#### hbr ------
m0 <- lmer(meanCoherence ~ (1|id), data = hbr_roi_df); summary(m0)
m1 <- lmer(meanCoherence ~ trial + (1|id), data = hbr_roi_df); summary(m1)
anova(m0, m1)

m2 <- lmer(meanCoherence ~ age + (1|id), data = hbr_roi_df); summary(m2)
anova(m0,m2)

m3 <- lmer(meanCoherence ~ age + conditionReleveled + (1|id), data = hbr_roi_df); summary(m3)
anova(m2,m3)

# Aim 3: With vs w/o face info ------ 
## hbo -----
hbo_df_collab <- hbo_df %>% filter(condition != 'ind')

m0 <- lmer(meanCoherence ~ (1|id), data = hbo_df_collab); summary(m0)

m1 <- lmer(meanCoherence ~ age + (1|id), data = hbo_df_collab); summary(m1)
anova(m0,m1)

m2 <- lmer(meanCoherence ~ trialReleveled + (1|id), data = hbo_df_collab); summary(m2)
anova(m0,m2)

m3 <- lmer(meanCoherence ~ conditionReleveled + (1|id), data = hbo_df_collab); summary(m3)
anova(m0,m3)

# standardized residuals versus fitted values by condition
plot(m3, resid(., scaled=TRUE) ~ fitted(.) | conditionReleveled, abline = 0)

# box-plots of residuals by subject
plot(m3, id ~ resid(., scaled = TRUE))

# observed versus fitted values by subject
plot(m3, meanCoherence ~ fitted(.) | id, abline = c(0.1))
res_m3 <- residuals(m3)

qqnorm(res_m3)
hist(res_m3)
ggplot(hbo_df_collab, aes(y = res_m3, x = condition)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  xlab("Predictor Value (condition)") +
  ylab("Residual") +
  ggtitle("Residuals vs. Predictor Plot")

## hbr ----- 
hbr_df_collab <- hbr_df %>% filter(condition != 'ind')

m0 <- lmer(meanCoherence ~ (1|id), data = hbr_df_collab); summary(m0)

m1 <- lmer(meanCoherence ~ age + (1|id), data = hbr_df_collab); summary(m1)
anova(m0,m1)

m2 <- lmer(meanCoherence ~ trialReleveled + (1|id), data = hbr_df_collab); summary(m2)
anova(m0,m2)

m3 <- lmer(meanCoherence ~ conditionReleveled + (1|id), data = hbr_df_collab); summary(m3)
anova(m0,m3)

# standardized residuals versus fitted values by condition
plot(m3, resid(., scaled=TRUE) ~ fitted(.) | conditionReleveled, abline = 0)

# box-plots of residuals by subject
plot(m3, id ~ resid(., scaled = TRUE))

# observed versus fitted values by subject
plot(m3, meanCoherence ~ fitted(.) | id, abline = c(0.1))
res_m3 <- residuals(m3)

qqnorm(res_m3)
hist(res_m3)
ggplot(hbr_df_collab, aes(y = res_m3, x = condition)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  xlab("Predictor Value (condition)") +
  ylab("Residual") +
  ggtitle("Residuals vs. Predictor Plot")

# Aim 4: Pseudodyad analysis --------
hbo_pda <- unpack(pda_hbo, level = 'channels', 18)
hbr_pda <- unpack(pda_hbr, level = 'channels', 18)

hbo_pda_roi <- unpack(pda_hbo, level = 'rois', 4)
hbr_pda_roi <- unpack(pda_hbr, level = 'rois', 4)

# should change condition labels based on your conditions
hbo_pda$condition <- ifelse(hbo_pda$condition == 1, "collab",
                             ifelse(hbo_pda$condition == 2, "collabScreen",
                                    ifelse(hbo_pda$condition == 3, "ind",
                                           ifelse(hbo_pda$condition == 4, "rest", hbo_pda$condition))))

hbo_pda <- hbo_pda %>% filter(condition != 'rest')

# should change condition labels based on your conditions
hbr_pda$condition <- ifelse(hbr_pda$condition == 1, "collab",
                            ifelse(hbr_pda$condition == 2, "collabScreen",
                                   ifelse(hbr_pda$condition == 3, "ind",
                                          ifelse(hbr_pda$condition == 4, "rest", hbr_pda$condition))))

hbr_pda <- hbr_pda %>% filter(condition != 'rest')

# should change condition labels based on your conditions
hbo_pda_roi$condition <- ifelse(hbo_pda_roi$condition == 1, "collab",
                            ifelse(hbo_pda_roi$condition == 2, "collabScreen",
                                   ifelse(hbo_pda_roi$condition == 3, "ind",
                                          ifelse(hbo_pda_roi$condition == 4, "rest", hbo_pda_roi$condition))))

hbo_pda_roi <- hbo_pda_roi %>% filter(condition != 'rest')

# should change condition labels based on your conditions
hbr_pda_roi$condition <- ifelse(hbr_pda_roi$condition == 1, "collab",
                            ifelse(hbr_pda_roi$condition == 2, "collabScreen",
                                   ifelse(hbr_pda_roi$condition == 3, "ind",
                                          ifelse(hbr_pda_roi$condition == 4, "rest", hbr_pda_roi$condition))))

hbr_pda_roi <- hbr_pda_roi %>% filter(condition != 'rest')

# 2 trials per condition
removePDAHbo <- hbo_pda %>% 
  dplyr::group_by(id, condition) %>%
  dplyr::filter(max(trial) == 1) %>%
  dplyr::distinct(id); removePDAHbo

hbo_pda <- hbo_pda[!hbo_pda$id %in% removePDAHbo$id, ]

removePDAHbr <- hbr_pda %>% 
  dplyr::group_by(id, condition) %>%
  dplyr::filter(max(trial) == 1) %>%
  dplyr::distinct(id); removePDAHbr

hbr_pda <- hbr_pda[!hbr_pda$id %in% removePDAHbr$id, ]

# 2 trials per condition
removePDAHboROI <- hbo_pda_roi %>% 
  dplyr::group_by(id, condition) %>%
  dplyr::filter(max(trial) == 1) %>%
  dplyr::distinct(id); removePDAHboROI

hbo_pda_roi <- hbo_pda_roi[!hbo_pda_roi$id %in% removePDAHboROI$id, ]

# 2 trials per condition
removePDAHbrROI <- hbr_pda_roi %>% 
  dplyr::group_by(id, condition) %>%
  dplyr::filter(max(trial) == 1) %>%
  dplyr::distinct(id); removePDAHbrROI

hbr_pda_roi <- hbr_pda_roi[!hbr_pda_roi$id %in% removePDAHbrROI$id, ]

# 06, 43 should be removed for ROI exclusion
hbo_pda <- hbo_pda %>% dplyr::filter(id != 'C06M06'& 
                                       id != 'C43M43')

hbr_pda <- hbr_pda %>% dplyr::filter(id != 'C06M06'& 
                                       id != 'C43M43')

hbo_pda_roi <- hbo_pda_roi %>% dplyr::filter(id != 'C06M06'& 
                                       id != 'C43M43')

hbr_pda_roi <- hbr_pda_roi %>% dplyr::filter(id != 'C06M06'& 
                                       id != 'C43M43')

names(hbo_pda)[names(hbo_pda)=='meanCoherence'] <- 'pda_coherence'
hbo_pda$trialReleveled <- as.factor(hbo_pda$trial)

names(hbr_pda)[names(hbr_pda)=='meanCoherence'] <- 'pda_coherence'
hbr_pda$trialReleveled <- as.factor(hbr_pda$trial)

names(hbo_pda_roi)[names(hbo_pda_roi)=='meanCoherence'] <- 'pda_coherence'
hbo_pda_roi$trialReleveled <- as.factor(hbo_pda_roi$trial)

names(hbr_pda_roi)[names(hbr_pda_roi)=='meanCoherence'] <- 'pda_coherence'
hbr_pda_roi$trialReleveled <- as.factor(hbr_pda_roi$trial)

hbo_pda_long <- hbo_pda %>% select(-c('type','trial')) %>%
  filter(condition != 'rest') %>%
  pivot_longer(cols = pda_coherence, 
               names_to = "coherence_type", 
               values_to = "coherence_value") %>%
  mutate(coherence_type = "pda") %>%
  mutate(conditionReleveled = as.factor(condition)) %>%
  select(-c(condition)) %>% 
  filter(coherence_value != 0); length(unique(hbo_pda_long$id))

hbo_df_long <- hbo_df %>% 
  pivot_longer(cols = meanCoherence,
               names_to = "coherence_type",
               values_to = "coherence_value") %>%
  mutate(coherence_type = "true") %>%
  select(-c('trial','type','age','condition')); length(unique(hbo_df_long$id))

hbr_pda_long <- hbr_pda %>% select(-c('type','trial')) %>%
  filter(condition != 'rest') %>%
  pivot_longer(cols = pda_coherence, 
               names_to = "coherence_type", 
               values_to = "coherence_value") %>%
  mutate(coherence_type = "pda") %>%
  mutate(conditionReleveled = as.factor(condition)) %>%
  select(-c(condition)) %>%
  filter(coherence_value != 0); length(unique(hbr_pda_long$id))

hbr_df_long <- hbr_df %>% 
  pivot_longer(cols = meanCoherence,
               names_to = "coherence_type",
               values_to = "coherence_value") %>%
  mutate(coherence_type = "true") %>%
  select(-c('trial','type','age','condition')); length(unique(hbr_df_long$id))

hbo_pda2 <- rbind(hbo_pda_long, hbo_df_long); length(unique(hbo_pda2$id))
hbr_pda2 <- rbind(hbr_pda_long, hbr_df_long); length(unique(hbr_pda2$id))

hbo_pda_roi_long <- hbo_pda_roi %>% select(-c('type','trial')) %>%
  pivot_longer(cols = pda_coherence, 
               names_to = "coherence_type", 
               values_to = "coherence_value") %>%
  mutate(coherence_type = "pda") %>%
  mutate(conditionReleveled = as.factor(condition)) %>%
  select(-c(condition)) %>% 
  filter(coherence_value != 0); length(unique(hbo_pda_roi_long$id))

hbo_roi_df_long <- hbo_roi_df %>% 
  pivot_longer(cols = meanCoherence,
               names_to = "coherence_type",
               values_to = "coherence_value") %>%
  mutate(coherence_type = "true") %>%
  mutate(trialReleveled = as.factor(trial)) %>% 
  select(-c('trial','type','age','condition')); length(unique(hbo_roi_df_long$id))

hbr_pda_roi_long <- hbr_pda_roi %>% select(-c('type','trial')) %>%
  filter(condition != 'rest') %>%
  pivot_longer(cols = pda_coherence, 
               names_to = "coherence_type", 
               values_to = "coherence_value") %>%
  mutate(coherence_type = "pda") %>%
  mutate(conditionReleveled = as.factor(condition)) %>%
  select(-c(condition)) %>%
  filter(coherence_value != 0); length(unique(hbr_pda_roi_long$id))

hbr_roi_df_long <- hbr_roi_df %>% 
  pivot_longer(cols = meanCoherence,
               names_to = "coherence_type",
               values_to = "coherence_value") %>%
  mutate(trialReleveled = as.factor(trial)) %>% 
  mutate(coherence_type = "true") %>%
  select(-c('trial','type','age','condition')); length(unique(hbr_roi_df_long$id))

hbo_pda_roi2 <- rbind(hbo_pda_roi_long, hbo_roi_df_long); length(unique(hbo_pda_roi2$id))
hbr_pda_roi2 <- rbind(hbr_pda_roi_long, hbr_roi_df_long); length(unique(hbr_pda_roi2$id))

# number of participants you have in PDA analysis
length(unique(hbo_pda2$id))
length(unique(hbr_pda2$id))

length(unique(hbo_pda_roi2$id))
length(unique(hbr_pda_roi2$id))

# number of coherence types per ID to triple check
table(hbo_pda2$coherence_type, hbo_pda2$id)
table(hbr_pda2$coherence_type, hbr_pda2$id)

table(hbo_pda_roi2$coherence_type, hbo_pda_roi2$id)
table(hbr_pda_roi2$coherence_type, hbr_pda_roi2$id)

hbo_pda2$coherence_type <- factor(hbo_pda2$coherence_type, level = c('true', 'pda'))
hbr_pda2$coherence_type <- factor(hbr_pda2$coherence_type, level = c('true', 'pda'))
hbo_pda_roi2$coherence_type <- factor(hbo_pda_roi2$coherence_type, level = c('true', 'pda'))
hbr_pda_roi2$coherence_type <- factor(hbr_pda_roi2$coherence_type, level = c('true', 'pda'))

# whole-brain results
hbo_pda_test <- hbo_pda2 %>%
  dplyr::group_by(id, conditionReleveled, coherence_type) %>%
  dplyr::mutate(avCoh = mean(coherence_value)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c('trialReleveled', 'coherence_value', 'pair')) %>%
  dplyr::distinct() %>%
  group_by(conditionReleveled) %>%
  t_test(avCoh ~ coherence_type, alternative = 'greater', paired = 'TRUE')

hbo_pda_test <- adjust_pvalue(hbo_pda_test, p.col= 'p', method = 'fdr')

# channel-wise 
channelPairs_t_tests_hbo_pda <- hbo_pda2 %>%
  dplyr::group_by(id, pair, coherence_type, conditionReleveled) %>%
  dplyr::mutate(avCoh = mean(coherence_value)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c('trialReleveled', 'coherence_value')) %>%
  dplyr::distinct() %>%
  group_by(conditionReleveled, pair) %>%
  t_test(avCoh ~ coherence_type, alternative = 'greater', paired = 'TRUE')

channelPairs_t_tests_hbo_pda <- adjust_pvalue(channelPairs_t_tests_hbo_pda, p.col= 'p', method = 'fdr')

channelPairs_t_tests_hbo_pda_sig <- channelPairs_t_tests_hbo_pda %>% filter(p < .05)

# ROI-wise level, C23M23, C49M49, C54M54, C45M45 are weird 
roi_t_tests_hbo_pda <- hbo_pda_roi2 %>%
  dplyr::filter(id != 'C23M23' & id != 'C49M49' & id != 'C54M54' & id != 'C45M45') %>%
  dplyr::group_by(id, pair, coherence_type, conditionReleveled) %>%
  dplyr::mutate(avCoh = mean(coherence_value)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c('trialReleveled', 'coherence_value')) %>%
  dplyr::distinct() %>%
  group_by(conditionReleveled, pair) %>%
  t_test(avCoh ~ coherence_type, alternative = 'greater', paired = 'TRUE')

roi_t_tests_hbo_pda <- adjust_pvalue(roi_t_tests_hbo_pda, p.col= 'p', method = 'fdr')

# whole-brain results
hbr_pda_test <- hbr_pda2 %>%
  dplyr::group_by(id, conditionReleveled, coherence_type) %>%
  dplyr::mutate(avCoh = mean(coherence_value)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c('trialReleveled', 'coherence_value', 'pair')) %>%
  dplyr::distinct() %>%
  group_by(conditionReleveled) %>%
  t_test(avCoh ~ coherence_type, alternative = 'greater', paired = 'TRUE')

hbr_pda_test <- adjust_pvalue(hbr_pda_test, p.col= 'p', method = 'fdr')

# channel-wise 
channelPairs_t_tests_hbr_pda <- hbr_pda2 %>%
  dplyr::group_by(id, pair, coherence_type, conditionReleveled) %>%
  dplyr::mutate(avCoh = mean(coherence_value)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c('trialReleveled', 'coherence_value')) %>%
  dplyr::distinct() %>%
  group_by(conditionReleveled, pair) %>%
  t_test(avCoh ~ coherence_type, alternative = 'greater', paired = 'TRUE')

channelPairs_t_tests_hbr_pda <- adjust_pvalue(channelPairs_t_tests_hbr_pda, p.col= 'p', method = 'fdr')

channelPairs_t_tests_hbr_pda_sig <- channelPairs_t_tests_hbr_pda %>% filter(p < .05)

# ROI-wise level, C23M23, C49M49, C54M54, C45M45 are weird 
roi_t_tests_hbr_pda <- hbr_pda_roi2 %>%
  dplyr::filter(id != 'C23M23' & id != 'C49M49' & id != 'C54M54' & id != 'C45M45') %>%
  dplyr::group_by(id, pair, coherence_type, conditionReleveled) %>%
  dplyr::mutate(avCoh = mean(coherence_value)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c('trialReleveled', 'coherence_value')) %>%
  dplyr::distinct() %>%
  group_by(conditionReleveled, pair) %>%
  t_test(avCoh ~ coherence_type, alternative = 'greater', paired = 'TRUE')

roi_t_tests_hbr_pda <- adjust_pvalue(roi_t_tests_hbr_pda, p.col= 'p', method = 'fdr')

# Aim 5a: Task performance -----------
perf <- na.omit(performance)

names(perf)[names(perf)=='ParticipantID'] <- 'id'
names(perf)[names(perf)=='Collab1....Puzzles.Given'] <- 'collab1_given'
names(perf)[names(perf)=='Collab1.....Puzzles.Finished.Correctly'] <- 'collab1_correct'
names(perf)[names(perf)=='Collab2.....Puzzles.Given'] <- 'collab2_given'
names(perf)[names(perf)=='Collab2.....Puzzles.Finished.Correctly'] <- 'collab2_correct'
names(perf)[names(perf)=='Collab3.....Puzzles.Given'] <- 'collab3_given'
names(perf)[names(perf)=='Collab3.....Puzzles.Finished.Correctly'] <- 'collab3_correct'
names(perf)[names(perf)=='CollabScreen1.....Puzzles.Given'] <- 'collabScreen1_given'
names(perf)[names(perf)=='CollabScreen1.....Puzzles.Finished.Correctly'] <- 'collabScreen1_correct'
names(perf)[names(perf)=='CollabScreen2.....Puzzles.Given'] <- 'collabScreen2_given'
names(perf)[names(perf)=='CollabScreen2.....Puzzles.Finished.Correctly'] <- 'collabScreen2_correct'
names(perf)[names(perf)=='CollabScreen3.....Puzzles.Given'] <- 'collabScreen3_given'
names(perf)[names(perf)=='CollabScreen3.....Puzzles.Finished.Correctly'] <- 'collabScreen3_correct'
names(perf)[names(perf)=='Individual1...CHILD.....Puzzles.Given'] <- 'ind1_c_given'
names(perf)[names(perf)=='Individual1...CHILD.....Puzzles.Finished.Correctly'] <- 'ind1_c_correct'
names(perf)[names(perf)=='Individual2...CHILD.....Puzzles.Given'] <- 'ind2_c_given'
names(perf)[names(perf)=='Individual2...CHILD.....Puzzles.Finished.Correctly'] <- 'ind2_c_correct'
names(perf)[names(perf)=='Individual3...CHILD.....Puzzles.Given'] <- 'ind3_c_given'
names(perf)[names(perf)=='Individual3...CHILD.....Puzzles.Finished.Correctly'] <- 'ind3_c_correct'

perf_long <- perf %>% 
  select(-matches("MUM")) %>%
  pivot_longer(!id) %>%
  dplyr::filter(!grepl('Finished', name)) %>%
  mutate(
    trial = str_extract(name, "\\d+"),
    condition = str_extract(name, "[a-zA-Z]+"),
    variable = case_when(
      grepl('given', name) ~ 'number_given',
      grepl('correct', name) ~ 'number_correct'
    )) %>% 
  select(-name)

perf_desc <- perf_long %>%
  group_by(trial, condition, variable) %>%
  summarize(
    Mean = mean(value),
    SD = sd(value),
    Min = min(value),
    Max = max(value))

perf_long2 <- perf_long %>%
  pivot_wider(names_from = variable,
              values_from = value)
  
## hbo ------
hbo_perf <- merge(hbo_df, perf_long2, by = c('id', 'condition', 'trial'))

m0 <- lmer(meanCoherence ~ conditionReleveled + (1|id), data = hbo_perf); summary(m0)
m1 <- lmer(meanCoherence ~ conditionReleveled*number_correct + (1|id), data = hbo_perf); summary(m1)
anova(m0,m1)

m2 <- lmer(meanCoherence ~ number_given + conditionReleveled*number_correct + (1|id), data = hbo_perf); summary(m2)
anova(m1, m2)

summary(m1)
report(m1)

# standardized residuals versus fitted values by number correct
plot(m1, resid(., scaled=TRUE) ~ fitted(.) | number_correct, abline = 0)

# box-plots of residuals by subject
plot(m1, id ~ resid(., scaled = TRUE))

# observed versus fitted values by subject
plot(m1, meanCoherence ~ fitted(.) | id, abline = c(0.1))
res_m1 <- residuals(m1)

qqnorm(res_m1)
hist(res_m1)
ggplot(hbo_perf, aes(y = res_m1, x = number_correct)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  xlab("Predictor Value (condition)") +
  ylab("Residual") +
  ggtitle("Residuals vs. Predictor Plot")

## hbr ------
hbr_perf <- merge(hbr_df, perf_long2, by = c('id', 'condition', 'trial'))

m0 <- lmer(meanCoherence ~ conditionReleveled + (1|id), data = hbr_perf); summary(m0)
m1 <- lmer(meanCoherence ~ conditionReleveled*number_correct + (1|id), data = hbr_perf); summary(m1)
anova(m0,m1)

m2 <- lmer(meanCoherence ~ number_given + conditionReleveled*number_correct + (1|id), data = hbr_perf); summary(m2)
anova(m1, m2)

summary(m2)
report(m2)

# standardized residuals versus fitted values by number correct
plot(m1, resid(., scaled=TRUE) ~ fitted(.) | number_correct, abline = 0)

# box-plots of residuals by subject
plot(m1, id ~ resid(., scaled = TRUE))

# observed versus fitted values by subject
plot(m1, meanCoherence ~ fitted(.) | id, abline = c(0.1))
res_m1 <- residuals(m1)

qqnorm(res_m1)
hist(res_m1)
ggplot(hbo_perf, aes(y = res_m1, x = number_correct)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  xlab("Predictor Value (condition)") +
  ylab("Residual") +
  ggtitle("Residuals vs. Predictor Plot")

# Aim 5b: Maternal Stress ----------
ques <- questionnaire %>% select('Collab.ID', 'pss_sum')

sd <- sd(ques$pss_sum)
m_pss <- mean(ques$pss_sum)

lowerlim <- m_pss-(sd*2)
upperlim <- m_pss+(sd*2)

hist(ques$pss_sum)
shapiro_test(ques$pss_sum)

outliers_pss <- ques %>%
  filter(pss_sum > upperlim | pss_sum < lowerlim); length(outliers_pss$Collab.ID)

pss_df <- ques %>%
  filter(pss_sum < upperlim & pss_sum > lowerlim)

shapiro_test(pss_df$pss_sum)
mean(pss_df$pss_sum, na.rm = T)
sd(pss_df$pss_sum, na.rm = T)
range(pss_df$pss_sum, na.rm = T)

hbo_ques <- hbo_df
hbo_ques$pss_sum <- NA

for (unique_id in unique(pss_df$Collab.ID)) {
  matching_rows <- grepl(unique_id, hbo_ques$id)
  
  if (any(matching_rows)) {
    hbo_ques$pss_sum[matching_rows] <- pss_df$pss_sum[pss_df$Collab.ID == unique_id]
  }
}

hbr_ques <- hbr_df
hbr_ques$pss_sum <- NA

for (unique_id in unique(pss_df$Collab.ID)) {
  matching_rows <- grepl(unique_id, hbr_ques$id)
  
  if (any(matching_rows)) {
    hbr_ques$pss_sum[matching_rows] <- pss_df$pss_sum[pss_df$Collab.ID == unique_id]
  }
}

hbo_ques <- hbo_ques %>% na.omit()
hbr_ques <- hbr_ques %>% na.omit()

## hbo ----------
m0 <- lmer(meanCoherence ~ conditionReleveled + (1|id), data = hbo_ques); summary(m0)
m1 <- lmer(meanCoherence ~ conditionReleveled*pss_sum + (1|id), data = hbo_ques); summary(m1)
anova(m0,m1)

report(m1)
summary(m1)

# standardized residuals versus fitted values by condition
plot(m1, resid(., scaled=TRUE) ~ fitted(.) | pss_sum, abline = 0)

# box-plots of residuals by subject
plot(m1, id ~ resid(., scaled = TRUE))

# observed versus fitted values by subject
plot(m1, meanCoherence ~ fitted(.) | id, abline = c(0.1))
res_m1 <- residuals(m3)

qqnorm(res_m1)
hist(res_m1)

## hbr ----------
m0 <- lmer(meanCoherence ~ conditionReleveled + (1|id), data = hbr_ques); summary(m0)
m1 <- lmer(meanCoherence ~ conditionReleveled*pss_sum + (1|id), data = hbr_ques); summary(m1)
anova(m0,m1)

report(m1)
summary(m1)
