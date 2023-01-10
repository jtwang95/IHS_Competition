library(tidyverse)
library(geepack)
library(mice)

fun <- function(id,dat,var_name){
  tmp = dat %>% filter(STUDY_PRTCPT_ID == id) %>% arrange(desc(Date))
  n = dim(tmp)[1]
  for (i in 1:n){
    if (!is.na(tmp[i,var_name])){
      return(tmp[i,c("STUDY_PRTCPT_ID","Date")])
    }
  }
  return(tmp[n,c("STUDY_PRTCPT_ID","Date")])
}
Indweight = read_csv("C:/Users/jtwan/Dropbox (University of Michigan)/IHS_competition/competition analysis/revision1/subject_id_weights.csv")
load("C:/Users/jtwan/Dropbox (University of Michigan)/IHS_competition/imputation data/2020_Intern_competition_datasets_imputation.Rdata")
IndQ1dat <- impute_list[[1]]$original_dat
IndBLdat <- impute_list[[1]]$baseline_original_dat
IndBLdat_complete <- impute_list[[1]]$baseline_complete_dat
IndQ1Weekdat <- IndQ1dat %>% group_by(STUDY_PRTCPT_ID,Week,CompType,TeamName,RivalName) %>% 
  dplyr::summarize(WEEK_STEP_COUNT = mean(STEP_COUNT,na.rm = TRUE),
                   WEEK_SLEEP_COUNT = mean(SLEEP_COUNT,na.rm = TRUE),
                   WEEK_MOOD_SCORE = mean(MOOD,na.rm = TRUE))

## dropout sensitivity
tmp = IndQ1dat
### step remove list
list_last_observation = do.call(rbind.data.frame,lapply(unique(tmp$STUDY_PRTCPT_ID),
                                                        FUN = fun, dat = tmp,var_name="STEP_COUNT"))
remove_list_steps = c()
for (i in 1:dim(list_last_observation)[1]){
  if (list_last_observation[i,"Date",drop=TRUE] == as.Date("2020-09-27")){
    next
  }
  remove_list_steps = rbind(remove_list_steps,
                      expand.grid(STUDY_PRTCPT_ID = list_last_observation[i,"STUDY_PRTCPT_ID",drop=TRUE], Date = seq(list_last_observation[i,"Date",drop=TRUE]+1,as.Date("2020-09-27"),1)))
}
list_last_observation = do.call(rbind.data.frame,lapply(unique(tmp$STUDY_PRTCPT_ID),
                                                        FUN = fun, dat = tmp,var_name="SLEEP_COUNT"))
### sleep remove list
remove_list_sleep = c()
for (i in 1:dim(list_last_observation)[1]){
  if (list_last_observation[i,"Date",drop=TRUE] == as.Date("2020-09-27")){
    next
  }
  remove_list_sleep = rbind(remove_list_sleep,
                            expand.grid(STUDY_PRTCPT_ID = list_last_observation[i,"STUDY_PRTCPT_ID",drop=TRUE], Date = seq(list_last_observation[i,"Date",drop=TRUE]+1,as.Date("2020-09-27"),1)))
}
### mood remove list
remove_list_mood = c()
for (i in 1:dim(list_last_observation)[1]){
  if (list_last_observation[i,"Date",drop=TRUE] == as.Date("2020-09-27")){
    next
  }
  remove_list_mood = rbind(remove_list_mood,
                            expand.grid(STUDY_PRTCPT_ID = list_last_observation[i,"STUDY_PRTCPT_ID",drop=TRUE], Date = seq(list_last_observation[i,"Date",drop=TRUE]+1,as.Date("2020-09-27"),1)))
}


## Weekly Missingness sensitivity
remove_list_steps = IndQ1dat %>%
  group_by(STUDY_PRTCPT_ID,Week) %>%
  summarise(na_count = sum(is.na(STEP_COUNT))) %>%
  filter(na_count >= 5) %>%
  select(STUDY_PRTCPT_ID,Week)
remove_list_sleep = IndQ1dat %>%
  group_by(STUDY_PRTCPT_ID,Week) %>%
  summarise(na_count = sum(is.na(SLEEP_COUNT))) %>%
  filter(na_count >= 5) %>%
  select(STUDY_PRTCPT_ID,Week)
remove_list_mood = IndQ1dat %>%
  group_by(STUDY_PRTCPT_ID,Week) %>%
  summarise(na_count = sum(is.na(MOOD))) %>%
  filter(na_count >= 5) %>%
  select(STUDY_PRTCPT_ID,Week)

### let's rerun this long stuff
remove_list = "mood"
fit.steps = list()
fit.sleeps = list()
fit.steps0 = list()
fit.sleeps0 = list()
fit.mood = list()
fit.mood0 = list()
fit.withins_steps = list()
fit.withins_sleep = list()
fit.withins_mood = list()
for (impute_iter in 1:impute_list$num_imputations){
  tmp = impute_list[[impute_iter]]$complete_dat %>%
    left_join(impute_list[[impute_iter]]$baseline_complete_dat) %>%
    inner_join(Indweight,by=c("STUDY_PRTCPT_ID"="UserID"))
  if (remove_list == "mood")
    tmp = tmp %>% anti_join(remove_list_mood)
  if (remove_list == "steps")
    tmp = tmp %>% anti_join(remove_list_steps)
  if (remove_list == "sleep")
    tmp = tmp %>% anti_join(remove_list_sleep)
  
  tmp1 <- tmp %>%
    group_by(TeamName, Week, CompType, RivalName) %>%
    summarise(
      WEEK_STEP_COUNT = mean(STEP_COUNT, na.rm = TRUE),
      WEEK_SLEEP_COUNT = mean(SLEEP_COUNT, na.rm = TRUE),
      WEEK_MOOD_SCORE = mean(MOOD, na.rm = TRUE),
      PreInternStep = mean(PreInternStep),
      PreInternSleep = mean(PreInternSleep),
      PreInternMood = mean(PreInternMood),
      PHQtot0 = mean(PHQtot0, na.rm = TRUE),
      Neu0 = mean(Neu0, na.rm = TRUE),
      depr0 = mean(depr0, na.rm = TRUE),
      EFE0 = mean(EFE0, na.rm = TRUE),
      GenderMale = mean(Gender == "Male"),
      ReportStep = sum(!is.na(STEP_COUNT)) / length(STEP_COUNT),
      ReportSleep = sum(!is.na(SLEEP_COUNT)) / length(SLEEP_COUNT),
      ReportMood = sum(!is.na(MOOD)) / length(MOOD),
      cluster_wt = sum(finalweight)
    ) %>% ungroup() %>%
    mutate(
      PreInternStep = PreInternStep - mean(PreInternStep),
      PreInternSleep = PreInternSleep - mean(PreInternSleep),
      PreInternMood = PreInternMood - mean(PreInternMood),
      PHQtot0 = PHQtot0 - mean(PHQtot0),
      Neu0 = Neu0 - mean(Neu0),
      depr0 = depr0 - mean(depr0),
      EFE0 = EFE0 - mean(EFE0)
    ) %>%
    mutate(TeamName = factor(TeamName)) %>%
    arrange(TeamName, Week) %>%
    group_by(TeamName) %>%
    mutate(
      PreStep = lag(WEEK_STEP_COUNT, 1),
      PreSleep = lag(WEEK_SLEEP_COUNT, 1),
      PreMood = lag(WEEK_MOOD_SCORE, 1)
    ) %>% ungroup() %>%
    mutate(
      PreStep = PreStep - mean(PreStep, na.rm = TRUE),
      PreSleep = PreSleep - mean(PreSleep, na.rm = TRUE),
      PreMood = PreMood - mean(PreMood, na.rm = TRUE)
    ) %>%
    ungroup %>% group_by(TeamName, Week, RivalName) %>%
    mutate(
      CompWithinInstitution = (
        str_split(TeamName, "_", simplify = TRUE)[1] == str_split(RivalName, "_", simplify = TRUE)[1]
      ),
      CompWithinProgram = (
        str_split(TeamName, "_", simplify = TRUE)[2] == str_split(RivalName, "_", simplify = TRUE)[2]
      ) &
        (str_split(TeamName, "_", simplify = TRUE)[2] != "Programs")
    ) %>%
    mutate(
      CompWithinInstitution = replace_na(as.integer(CompWithinInstitution), 0),
      CompWithinProgram = replace_na(as.integer(CompWithinProgram), 0)
    )
  ## remove dropout
  if (remove_list == "steps"){
    fit.steps0[[impute_iter]] = geeglm(
      WEEK_STEP_COUNT ~ Week + CompType + GenderMale + PreInternStep + PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
      data = tmp1 %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
      id = TeamName,
      weights = cluster_wt
    )
    fit.steps[[impute_iter]] = geeglm(
      WEEK_STEP_COUNT ~ Week * CompType + GenderMale + PreInternStep + PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
      data = tmp1 %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
      id = TeamName,
      weights = cluster_wt
    )
    fit.withins_steps[[impute_iter]] = geeglm(
      WEEK_STEP_COUNT ~ Week * CompType + CompWithinInstitution + CompWithinProgram + GenderMale + PreInternStep + PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
      data = tmp1 %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
      id = TeamName,
      weights = cluster_wt
    )
  }
  if (remove_list == "sleep"){
    fit.sleeps0[[impute_iter]] = geeglm(
      WEEK_SLEEP_COUNT ~ Week + CompType + GenderMale + PreInternSleep + PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
      data = tmp1 %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
      id = TeamName,
      weights = cluster_wt
    )
    fit.sleeps[[impute_iter]] = geeglm(
      WEEK_SLEEP_COUNT ~ Week * CompType + GenderMale + PreInternSleep + PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
      data = tmp1 %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
      id = TeamName,
      weights = cluster_wt
    )
    fit.withins_sleep[[impute_iter]] = geeglm(
      WEEK_SLEEP_COUNT ~ Week * CompType + CompWithinInstitution + CompWithinProgram + GenderMale + PreInternSleep + PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
      data = tmp1 %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
      id = TeamName,
      weights = cluster_wt
    )
  }
  if (remove_list == "mood"){
    fit.mood0[[impute_iter]] =geeglm(
      WEEK_MOOD_SCORE ~ Week + Comp + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreInternStep + PreStep + PreInternSleep + PreSleep,
      data = tmp1 %>% filter(!is.na(PreStep)) %>% filter(!is.na(PreSleep))  %>% mutate(Comp = CompType != "noncomp"),
      id = TeamName,
      weights = cluster_wt
    )
    fit.mood[[impute_iter]] = geeglm(
      WEEK_MOOD_SCORE ~ Week * Comp + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreInternStep + PreStep + PreInternSleep + PreSleep,
      data = tmp1 %>% filter(!is.na(PreStep)) %>% filter(!is.na(PreSleep))  %>% mutate(Comp = CompType != "noncomp"),
      id = TeamName,
      weights = cluster_wt
    )
  }
  #################################################################################
  
}
if (remove_list == "mood"){
  print(summary(pool(fit.mood0),conf.int=TRUE))
  print(summary(pool(fit.mood),conf.int=TRUE))
}
if (remove_list == "steps"){
  print(summary(pool(fit.steps0),conf.int=TRUE))
  print(summary(pool(fit.steps),conf.int=TRUE))
  print(summary(pool(fit.withins_steps),conf.int=TRUE))
}
if (remove_list == "sleep"){
  print(summary(pool(fit.sleeps0),conf.int=TRUE))
  print(summary(pool(fit.sleeps),conf.int=TRUE))
  print(summary(pool(fit.withins_sleep),conf.int=TRUE))
}