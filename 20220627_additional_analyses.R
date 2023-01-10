library(tidyverse)
library(ggplot2)
## setting
setwd("~/../Dropbox (University of Michigan)/IHS_competition/competition analysis/")
min_date = as.Date("2020-07-01")
max_date = as.Date("2020-12-31")
M = 20
random_seed = 20210902
set.seed(random_seed)
min_sleep = 30
max_sleep = 1440
imputation_method = "pmm"

## individual data with id, date, weekday, step, sleep, mood, name, rival_name, competition_group, notification_type

### load data
competelist <-
  read_csv('../competition/2020 access dode_Competition subjects.csv') %>%
  mutate(TeamName = paste(Institution, Program, sep = "_"))

IndQ1Q2dat <-
  read_csv("../sleep step mood/step_sleep_mood_daily_2020_Q1.csv") %>%
  bind_rows(read_csv("../sleep step mood/step_sleep_mood_daily_2020_Q2.csv")) %>%
  mutate(
    Date = as.Date(METRIC_START_DATE, "%d-%b-%y"),
    STEP_COUNT = na_if(STEP_COUNT, 0),
    SLEEP_COUNT = na_if(SLEEP_COUNT, 0),
    MOOD = na_if(MOOD, 0)
  ) %>%
  select(STUDY_PRTCPT_ID, Date, STEP_COUNT, SLEEP_COUNT, MOOD) %>%
  filter((Date >= min_date) & (Date <= max_date))
IndQ1Q2dat <-
  expand.grid(
    STUDY_PRTCPT_ID = unique(IndQ1Q2dat$STUDY_PRTCPT_ID),
    Date = seq(min_date, max_date, 1)
  ) %>%
  left_join(IndQ1Q2dat) %>%
  mutate(Week = as.integer(format(Date, "%V")) - as.integer(format(min(Date), "%V")) -
           1) %>% tibble()
IndQ1Q2dat <-
  IndQ1Q2dat %>% mutate(SLEEP_COUNT = if_else(
    (SLEEP_COUNT < min_sleep) |
      (SLEEP_COUNT > max_sleep),
    NA_real_,
    SLEEP_COUNT
  ))

IndSuvdat <-
  read_csv("../survey data/IHSdata_2020blq1q2_01042021.csv") %>%
  mutate(
    STUDY_PRTCPT_ID = UserID,
    Gender = dplyr::recode(Sex, `1` = "Male", `2` = "Female")
  ) %>%
  select(STUDY_PRTCPT_ID, Gender, PHQtot0, Neu0, depr0, EFE0) %>%
  filter(!is.na(Gender))

GrpQ1dat <-
  read_csv("../competition/Competition_long_Q1.csv") %>% bind_rows(read_csv("../competition/Competition_long_Q2.csv")) %>%
  filter((TEAM1_INSTITUTION %in% unique(competelist$Institution)) &
           (TEAM1_PROGRAM %in% unique(competelist$Program))) %>% # remove "test" data
  mutate(
    Date = as.Date(CALCULATION_DATE, "%d-%b-%y"),
    # change date format
    TeamName = paste(TEAM1_INSTITUTION, TEAM1_PROGRAM, sep = "_"),
    # combine institution and program to team name
    RivalName = paste(TEAM2_INSTITUTION, TEAM2_PROGRAM, sep = "_"),
    # the rival team's name
    CompType = NAME
  ) %>% # Competition type, either steps or sleep
  select(TeamName, Date, RivalName, CompType)
GrpQ1dat <-
  GrpQ1dat %>% bind_rows(GrpQ1dat %>% rename(RivalName = TeamName, TeamName = RivalName))

IndMsgdat <-
  read_csv("../competition/2020 Internship Notification Type - Daily Level.csv") %>%
  mutate(Date = as.Date(Date, "%d-%b-%y")) %>%
  select(STUDY_PRTCPT_ID, Date, Notification_type)

IndBLdat <-
  read_csv("../sleep step mood/step_sleep_mood_daily_2020_BL.csv") %>%
  mutate(
    Date = as.Date(METRIC_START_DATE, "%d-%b-%y"),
    STEP_COUNT = na_if(STEP_COUNT, 0),
    SLEEP_COUNT = na_if(SLEEP_COUNT, 0),
    MOOD = na_if(MOOD, 0)
  ) %>%
  group_by(STUDY_PRTCPT_ID) %>%
  summarise(
    PreInternStep = mean(STEP_COUNT, na.rm = TRUE),
    PreInternSleep = mean(SLEEP_COUNT, na.rm = TRUE),
    PreInternMood = mean(MOOD, na.rm = TRUE)
  ) %>%
  select(STUDY_PRTCPT_ID, PreInternStep, PreInternSleep, PreInternMood)

### determine which subjects will included in imputation and analysis
ParticipantIDs = Reduce(
  intersect,
  list(
    competelist$AccessCode,
    read_csv("../sleep step mood/step_sleep_mood_daily_2020_Q1.csv") %>%
      bind_rows(
        read_csv("../sleep step mood/step_sleep_mood_daily_2020_Q2.csv")
      ) %>%
      mutate(Date = as.Date(METRIC_START_DATE, "%d-%b-%y")) %>%
      filter(Date >= as.Date("2020-07-06")) %>% pull(STUDY_PRTCPT_ID),
    IndMsgdat$STUDY_PRTCPT_ID,
    IndSuvdat$STUDY_PRTCPT_ID,
    IndBLdat$STUDY_PRTCPT_ID
  )
)


IndQ1Q2dat %>% filter(Week > -1) %>% filter(STUDY_PRTCPT_ID %in% ParticipantIDs) %>% group_by(Week) %>% summarise(across(c(STEP_COUNT, SLEEP_COUNT),
                                                                                                                         .fns = ~ mean(!is.na(.x)))) %>%
  pivot_longer(cols = !starts_with("Week")) %>%
  ggplot(aes(x = Week, y = value, color = name)) + geom_line(size = 1.2) + geom_vline(xintercept = 11, linetype =
                                                                                        "dashed") +
  ylab("Percent of users with nonmissing observations") + scale_x_continuous(breaks = 0:25) +
  scale_color_discrete(name = "data type",
                       labels = c("sleep minutes", "step count"))

#########################################################################################################################################
# review 1 point 3 analyses
#########################################################################################################################################
library(geepack)
library(ggplot2)
library(tidyverse)
library(sjPlot)
library(mice)

confint.geeglm <- function(object, parm, level = 0.95, ...) {
  cc <- coef(summary(object))
  mult <- qnorm((1 + level) / 2)
  citab <- with(
    as.data.frame(cc),
    cbind(
      est = Estimate,
      lwr = Estimate - mult * Std.err,
      upr = Estimate + mult * Std.err,
      pvalue = `Pr(>|W|)`
    )
  )
  rownames(citab) <- rownames(cc)
  citab[parm,]
}

if (.Platform$OS.type == "unix")
{
  prefix = "~/Dropbox (University of Michigan)/IHS_competition/"
} else
{
  prefix = "C:/Users/jtwan/Dropbox (University of Michigan)/IHS_competition/"
}

load(
  paste(
    prefix,
    "imputation data/2020_Intern_competition_datasets_imputation.Rdata",
    sep = ""
  )
)

device_info = read_csv(paste(
  prefix,
  "survey data/Participant_WearableDevices_2020.csv",
  sep = ""
)) %>% filter(is.na(PRTCPT_DVC_END_DT)) %>% select(STUDY_PRTCPT_ID, PRTCPT_DVC_MODEL) %>%
  distinct(STUDY_PRTCPT_ID, .keep_all = TRUE)

Indweight = read_csv(paste(prefix,"competition analysis/revision1/subject_id_weights.csv",sep=""))
IndQ1dat <- impute_list[[1]]$original_dat
IndBLdat_complete <-
  impute_list[[1]]$baseline_complete_dat %>% left_join(device_info) %>%
  left_join(Indweight,by=c("STUDY_PRTCPT_ID"="UserID"))
## 125 interns not have device information
## randomly set device for these interns
set.seed(20220614)
IndBLdat_complete[is.na(IndBLdat_complete$PRTCPT_DVC_MODEL), "PRTCPT_DVC_MODEL"] = sample(c("Charge HR", "Apple Watch"),
                                                                                          size = 125,
                                                                                          replace = TRUE)
TeamQ1Weekdat <- IndQ1dat %>% left_join(IndBLdat_complete) %>%
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
    apple_pct = sum(PRTCPT_DVC_MODEL == "Apple Watch") / length(PRTCPT_DVC_MODEL),
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
    PreMood = lag(WEEK_MOOD_SCORE, 1),
    PreReportStep = lag(ReportStep, 1),
    PreReportSleep = lag(ReportSleep, 1),
    PreReportMood = lag(ReportMood, 1)
  ) %>% ungroup() %>%
  mutate(
    PreStep = PreStep - mean(PreStep, na.rm = TRUE),
    PreSleep = PreSleep - mean(PreSleep, na.rm = TRUE),
    PreMood = PreMood - mean(PreMood, na.rm = TRUE),
    PreReportStep = PreReportStep - mean(PreReportStep, na.rm = TRUE),
    PreReportSleep = PreReportSleep - mean(PreReportSleep, na.rm =
                                             TRUE),
    PreReportMood = PreReportMood - mean(PreReportMood, na.rm = TRUE)
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
  ) %>%
  ungroup()

tmp = TeamQ1Weekdat %>% select(TeamName,
                               Week,
                               CompType,
                               RivalName,
                               WEEK_STEP_COUNT,
                               WEEK_SLEEP_COUNT)
tmp1 = tmp1 = TeamQ1Weekdat %>% select(TeamName, Week, WEEK_STEP_COUNT, WEEK_SLEEP_COUNT) %>% rename(
  RivalName = TeamName,
  RIVAL_WEEK_STEP_COUNT = WEEK_STEP_COUNT,
  RIVAL_WEEK_SLEEP_COUNT = WEEK_SLEEP_COUNT
)
tmp = tmp %>% left_join(tmp1) %>% arrange(TeamName, Week) %>% mutate(
  win_loss = case_when(
    is.na(RIVAL_WEEK_STEP_COUNT) &
      is.na(RIVAL_WEEK_SLEEP_COUNT) ~ "non",
    WEEK_STEP_COUNT > RIVAL_WEEK_STEP_COUNT &
      CompType == "steps" ~ "win",
    WEEK_SLEEP_COUNT > RIVAL_WEEK_SLEEP_COUNT &
      CompType == "sleep" ~ "win",
    TRUE ~ "loss"
  )
) %>% select(TeamName, Week, win_loss) %>% mutate(win_loss = factor(win_loss, levels =
                                                                      c("non", "win", "loss"))) %>%
  mutate(win_loss_current = win_loss,
         win_loss_prev = lag(win_loss, 1))
TeamQ1Weekdat = TeamQ1Weekdat %>% left_join(tmp)

## Exploratory analysis - participation rate
tmp = TeamQ1Weekdat
fit.pr_step0 = geeglm(
  ReportStep ~ Week + win_loss_prev + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportStep,
  data = tmp %>% filter(!is.na(PreReportStep)),
  id = TeamName,
  corstr = "independence",weights = cluster_wt
)
# fit.pr_step = geeglm(
#   ReportStep ~ Week * win_loss + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportStep,
#   data = tmp %>% filter(!is.na(PreReportStep)),
#   id = TeamName,
#   corstr = "independence"
# )
fit.pr_sleep0 = geeglm(
  ReportSleep ~ Week + win_loss_prev + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportSleep,
  data = tmp %>% filter(!is.na(PreReportSleep)),
  id = TeamName,
  corstr = "independence",weights = cluster_wt
)
# fit.pr_sleep = geeglm(
#   ReportSleep ~ Week * win_loss + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportSleep,
#   data = tmp %>% filter(!is.na(PreReportSleep)),
#   id = TeamName,
#   corstr = "independence"
# )
fit.pr_mood0 = geeglm(
  ReportMood ~ Week + win_loss_prev + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportMood,
  data = tmp %>% filter(!is.na(PreReportMood)),
  id = TeamName,
  corstr = "independence",weights = cluster_wt
)
# fit.pr_mood = geeglm(
#   ReportMood ~ Week * win_loss + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportMood,
#   data = tmp %>% filter(!is.na(PreReportMood)),
#   id = TeamName,
#   corstr = "independence"
# )


## Exploratory analysis - device type
tmp = TeamQ1Weekdat %>% arrange(TeamName, Week)
### step count
fit.ind = geeglm(
  WEEK_STEP_COUNT ~ Week + CompType * apple_pct + GenderMale + PreInternStep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
  data = tmp %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
  id = TeamName,
  corstr = "independence"
)
summary(fit.ind)


### sleep minute
fit.ind = geeglm(
  WEEK_SLEEP_COUNT ~ Week + CompType * apple_pct + GenderMale + PreInternSleep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
  data = tmp %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
  id = TeamName,
  corstr = "independence"
)
summary(fit.ind)

## mood moderation on competition on step count / sleep minutes
tmp = TeamQ1Weekdat %>% arrange(TeamName, Week) %>% group_by(TeamName) %>% mutate(WEEK_MOOD_SCORE_PREV = lag(WEEK_MOOD_SCORE, 1)) %>% ungroup()
### step count
fit.ind = geeglm(
  WEEK_STEP_COUNT ~ Week + CompType * WEEK_MOOD_SCORE_PREV + GenderMale + PreInternStep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
  data = tmp %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
  id = TeamName,
  corstr = "independence",
  weights = cluster_wt
)
summary(fit.ind)

### step count
fit.ind = geeglm(
  WEEK_SLEEP_COUNT ~ Week + CompType * WEEK_MOOD_SCORE_PREV + GenderMale + PreInternSleep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
  data = tmp %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
  id = TeamName,
  corstr = "independence",
  weights = cluster_wt
)
summary(fit.ind)

# MI
fits.cluster_step0 = list()
fits.cluster_step = list()
fits.cluster_intra_step = list()
fits.cluster_sleep0 = list()
fits.cluster_sleep = list()
fits.cluster_intra_sleep = list()
fits.cluster_mood0 = list()
fits.cluster_mood = list()

## new
fits.cluster_step0_dvc = list()
fits.cluster_sleep0_dvc = list()
fits.cluster_step0_mood = list()
fits.cluster_sleep0_mood = list()
fits.cluster_step0_byspecialty = list()
fits.cluster_sleep0_byspecialty = list()

results = tibble()
for (i in 1:impute_list$num_imputations) {
  tmp_bl = impute_list[[i]]$baseline_complete_dat %>% left_join(device_info)
  ## 125 interns not have device information
  ## randomly set device for these interns
  set.seed(20220614)
  tmp_bl[is.na(tmp_bl$PRTCPT_DVC_MODEL), "PRTCPT_DVC_MODEL"] = sample(c("Charge HR", "Apple Watch"),
                                                                      size = 125,
                                                                      replace = TRUE)
  tmp = impute_list[[i]]$complete_dat %>%
    left_join(tmp_bl) %>% left_join(Indweight,by=c("STUDY_PRTCPT_ID"="UserID"))
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
      apple_pct = mean(PRTCPT_DVC_MODEL == "Apple Watch"),
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
  ####
  fit.step0 = geeglm(
    WEEK_STEP_COUNT ~ Week + CompType + GenderMale + PreInternStep + PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
    data = tmp1 %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
    id = TeamName,
    weights = cluster_wt
  )
  fit.step = geeglm(
    WEEK_STEP_COUNT ~ Week * CompType + GenderMale + PreInternStep + PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
    data = tmp1 %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
    id = TeamName,
    weights = cluster_wt
  )
  fit.within_step = geeglm(
    WEEK_STEP_COUNT ~ Week * CompType + CompWithinInstitution + CompWithinProgram + GenderMale + PreInternStep + PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
    data = tmp1 %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
    id = TeamName,
    weights = cluster_wt
  )
  
  fit.sleep0 = geeglm(
    WEEK_SLEEP_COUNT ~ Week + CompType + GenderMale + PreInternSleep + PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
    data = tmp1 %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
    id = TeamName,
    weights = cluster_wt
  )
  fit.sleep = geeglm(
    WEEK_SLEEP_COUNT ~ Week * CompType + GenderMale + PreInternSleep + PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
    data = tmp1 %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
    id = TeamName,
    weights = cluster_wt
  )
  fit.within_sleep = geeglm(
    WEEK_SLEEP_COUNT ~ Week * CompType + CompWithinInstitution + CompWithinProgram + GenderMale + PreInternSleep + PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
    data = tmp1 %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
    id = TeamName,
    weights = cluster_wt
  )
  
  tmp1 = tmp1 %>% mutate(Comp = CompType != "noncomp")
  fit.mood0 = geeglm(
    WEEK_MOOD_SCORE ~ Week + Comp + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreInternStep + PreStep + PreInternSleep + PreSleep,
    data = tmp1 %>% filter(!is.na(PreStep)) %>% filter(!is.na(PreSleep)),
    id = TeamName,
    weights = cluster_wt
  )
  fit.mood = geeglm(
    WEEK_MOOD_SCORE ~ Week * Comp + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreInternStep + PreStep + PreInternSleep + PreSleep,
    data = tmp1 %>% filter(!is.na(PreStep)) %>% filter(!is.na(PreSleep)),
    id = TeamName,
    weights = cluster_wt
  )
  fit.step0_dvc = geeglm(
    WEEK_STEP_COUNT ~ Week + CompType * apple_pct + GenderMale + PreInternStep + PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
    data = tmp1 %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
    id = TeamName,
    weights = cluster_wt
  )
  fit.sleep0_dvc = geeglm(
    WEEK_SLEEP_COUNT ~ Week + CompType * apple_pct + GenderMale + PreInternSleep + PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
    data = tmp1 %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
    id = TeamName,
    weights = cluster_wt
  )
  tmp1 = tmp1 %>% arrange(TeamName, Week) %>% group_by(TeamName) %>% mutate(WEEK_MOOD_SCORE_PREV = lag(WEEK_MOOD_SCORE, 1)) %>% ungroup()

  fit.step0_mood = geeglm(
    WEEK_STEP_COUNT ~ Week + CompType * WEEK_MOOD_SCORE_PREV + GenderMale + PreInternStep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
    data = tmp1 %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
    id = TeamName,
    corstr = "independence",
    weights = cluster_wt
  )
  fit.sleep0_mood = geeglm(
    WEEK_SLEEP_COUNT ~ Week + CompType * WEEK_MOOD_SCORE_PREV + GenderMale + PreInternSleep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
    data = tmp1 %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
    id = TeamName,
    corstr = "independence",
    weights = cluster_wt
  )
  fit.step0_byspecialty = geeglm(
    WEEK_STEP_COUNT ~ Week + CompType : Specialty + Specialty + GenderMale + PreInternStep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
    data = tmp1 %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep))  %>% 
      separate(TeamName,c("Institution","Specialty"),sep="_",remove = FALSE) %>%
      mutate(Specialty=relevel(factor(Specialty),ref="Programs")),
    id = TeamName,
    corstr = "independence",
    weights = cluster_wt
  )
  fit.sleep0_byspecialty = geeglm(
    WEEK_SLEEP_COUNT ~ Week + CompType : Specialty +Specialty + GenderMale + PreInternSleep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
    data = tmp1 %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep))%>% separate(TeamName,c("Institution","Specialty"),sep="_",remove = FALSE),
    id = TeamName,
    corstr = "independence",
    weights = cluster_wt
  )
  
  fits.cluster_step0[[i]] = fit.step0
  fits.cluster_step[[i]] = fit.step
  fits.cluster_intra_step[[i]] = fit.within_step
  fits.cluster_sleep0[[i]] = fit.sleep0
  fits.cluster_sleep[[i]] = fit.sleep
  fits.cluster_intra_sleep[[i]] = fit.within_sleep
  fits.cluster_mood0[[i]] = fit.mood0
  fits.cluster_mood[[i]] = fit.mood
  fits.cluster_step0_dvc[[i]] = fit.step0_dvc
  fits.cluster_sleep0_dvc[[i]] = fit.sleep0_dvc
  fits.cluster_step0_mood[[i]] = fit.step0_mood
  fits.cluster_sleep0_mood[[i]] = fit.sleep0_mood
  fits.cluster_step0_byspecialty[[i]] = fit.step0_byspecialty
  fits.cluster_sleep0_byspecialty[[i]] = fit.sleep0_byspecialty
}

summary(pool(fits.cluster_step0), conf.int = TRUE)
summary(pool(fits.cluster_step), conf.int = TRUE)
summary(pool(fits.cluster_intra_step), conf.int = TRUE)
summary(pool(fits.cluster_sleep0), conf.int = TRUE)
summary(pool(fits.cluster_sleep), conf.int = TRUE)
summary(pool(fits.cluster_intra_sleep), conf.int = TRUE)
summary(pool(fits.cluster_mood0), conf.int = TRUE)
summary(pool(fits.cluster_mood), conf.int = TRUE)
summary(pool(fits.cluster_step0_dvc), conf.int = TRUE)
summary(pool(fits.cluster_sleep0_dvc), conf.int = TRUE)
summary(pool(fits.cluster_step0_mood), conf.int = TRUE)
summary(pool(fits.cluster_sleep0_mood), conf.int = TRUE)
summary(pool(fits.cluster_step0_byspecialty),conf.int = TRUE)
summary(pool(fits.cluster_sleep0_byspecialty),conf.int = TRUE)


# draw specialty-specific plot
summary(pool(fits.cluster_step0_byspecialty),conf.int = TRUE) %>% 
  select(term,estimate,`2.5 %`,`97.5 %`) %>%
  filter(str_detect(term,":")) %>%
  mutate(specialty = str_replace(term,"CompTypesteps:Specialty",""),
         lb = `2.5 %`,
         ub = `97.5 %`) %>%
  select(specialty,estimate,lb,ub) %>% 
  arrange(estimate) %>% 
  mutate(specialty = factor(specialty,levels = c("Medicine-Pediatrics","Surgery","Pediatrics","Programs",
                                                 "Internal Medicine","Psychiatry","OBGYN","Family Medicine",
                                                 "Anesthesiology","Emergency Medicine","Neurology","Transitional"))) %>%
  ggplot(aes(x=specialty,y=estimate)) + geom_point(size=3.0) + 
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.3,size=1.0) +
  geom_hline(yintercept = 0,linetype="dashed",size=1.0,color="red") +
  ylab("Effect of team competition on daily step count") + 
  xlab("Specialty") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,vjust = 0.6))-> p1
  
summary(pool(fits.cluster_sleep0_byspecialty),conf.int = TRUE) %>% 
  select(term,estimate,`2.5 %`,`97.5 %`) %>%
  filter(str_detect(term,":")) %>%
  mutate(specialty = str_replace(term,"CompTypesleep:Specialty",""),
         lb = `2.5 %`,
         ub = `97.5 %`) %>%
  select(specialty,estimate,lb,ub) %>% 
  arrange(estimate) %>% 
  mutate(specialty = factor(specialty,levels = c("Medicine-Pediatrics","Anesthesiology","Neurology","Family Medicine",
                                                 "Programs","Pediatrics","Internal Medicine","Transitional",
                                                 "Surgery","Psychiatry","OBGYN","Emergency Medicine"))) %>%
  ggplot(aes(x=specialty,y=estimate)) + geom_point(size=3.0) + 
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.3,size=1.0) +
  geom_hline(yintercept = 0,linetype="dashed",size=1.0,color="red") +
  ylab("Effect of team competition on daily sleep minutes") + 
  xlab("Specialty") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.6)) -> p2
library(ggpubr)
png(file="C://Users/jtwan/Downloads/effect_specialty.png",width=1500, height=750,res=150)
ggarrange(p1,p2,ncol = 2,nrow = 1,labels = c("a","b")) + theme(plot.margin = margin(0.0,0.5,0.0,0.0, "cm")) 
dev.off()
