## Code to impute 2020 IHS computation dataset

library(mice)
library(tidyverse)
options(nwarnings = 10000000)

# set working dictionary
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("~/../Dropbox (University of Michigan)/IHS_competition/competition analysis/")
setwd("~/Dropbox (University of Michigan)/IHS_competition/competition analysis/")
# prepare all the data needed

## setting
min_date = as.Date("2020-07-01")
max_date = as.Date("2020-09-27")
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

IndQ1dat <-
  read_csv("../sleep step mood/step_sleep_mood_daily_2020_Q1.csv") %>%
  mutate(
    Date = as.Date(METRIC_START_DATE, "%d-%b-%y"),
    STEP_COUNT = na_if(STEP_COUNT, 0),
    SLEEP_COUNT = na_if(SLEEP_COUNT, 0),
    MOOD = na_if(MOOD, 0)
  ) %>%
  select(STUDY_PRTCPT_ID, Date, STEP_COUNT, SLEEP_COUNT, MOOD) %>%
  filter((Date >= min_date) & (Date <= max_date))
IndQ1dat <-
  expand.grid(
    STUDY_PRTCPT_ID = unique(IndQ1dat$STUDY_PRTCPT_ID),
    Date = seq(min_date, max_date, 1)
  ) %>%
  left_join(IndQ1dat) %>%
  mutate(Week = as.integer(format(Date, "%V")) - as.integer(format(min(Date), "%V"))-1) %>% tibble()
IndQ1dat <-
  IndQ1dat %>% mutate(SLEEP_COUNT = if_else(
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

GrpQ1dat <- read_csv("../competition/Competition_long_Q1.csv") %>%
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
      mutate(Date = as.Date(METRIC_START_DATE, "%d-%b-%y")) %>%
      filter(Date >= as.Date("2020-07-06")) %>% pull(STUDY_PRTCPT_ID),
    IndMsgdat$STUDY_PRTCPT_ID,
    IndSuvdat$STUDY_PRTCPT_ID,
    IndBLdat$STUDY_PRTCPT_ID
  )
) # find common participants in all datasets, 1779
competelist <-
  competelist %>% filter(AccessCode %in% ParticipantIDs)
IndQ1dat <- IndQ1dat %>% filter(STUDY_PRTCPT_ID %in% ParticipantIDs)
IndMsgdat <-
  IndMsgdat %>% filter(STUDY_PRTCPT_ID %in% ParticipantIDs)
IndSuvdat <-
  IndSuvdat %>% filter(STUDY_PRTCPT_ID %in% ParticipantIDs)
IndBLdat <-
  IndBLdat %>% filter(STUDY_PRTCPT_ID %in% ParticipantIDs) %>% left_join(IndSuvdat)



# impute baseline average step,sleep and mood before 07-01-2020
impute_list = list()
for (m in 1:M) {
  IndQ1dat_complete = tibble()
  ## variable in personal baseline
  ## "STUDY_PRTCPT_ID","PreInternStep","PreInternSleep","PreInternMood","Gender","PHQtot0","Neu0","depr0","EFE0"
  IndBLdat_impute <- IndBLdat %>% select(-c(STUDY_PRTCPT_ID))
  mice_output <-
    mice(IndBLdat_impute, m = 1, method = imputation_method)
  IndBLdat_complete <-
    IndBLdat %>% select(STUDY_PRTCPT_ID) %>% bind_cols(complete(mice_output, 1))
  
  # impute 07-01 to 07-05
  
  ## read data
  tmpdat <- IndQ1dat %>% right_join(IndBLdat_complete %>% distinct(STUDY_PRTCPT_ID)) %>% arrange(STUDY_PRTCPT_ID,Date)
  
  ## 07-01
  CurrentDate = as.Date("2020-07-01")
  HistoryMat = IndBLdat_complete
  CurrentMat = tmpdat %>% filter(Date == CurrentDate) %>%
    select(-c(Date)) %>%
    rename_with(
      .fn = ~ paste(.x, format(CurrentDate, "%m%d"), sep = "_"),
      .cols = c(STEP_COUNT, SLEEP_COUNT, MOOD)
    )
  ImputeMat = CurrentMat %>% left_join(HistoryMat,by = "STUDY_PRTCPT_ID") %>% select(-c(STUDY_PRTCPT_ID,Week))
  mice_output <- mice(ImputeMat, m = 1, method = imputation_method)
  CompleteMatAll <- CurrentMat %>% select(STUDY_PRTCPT_ID, Week) %>% bind_cols(complete(mice_output, 1)) %>%
     mutate(CompType = "noncomp",Notification_type="nomsg")
  IndQ1dat_complete <- CompleteMatAll %>% select(STUDY_PRTCPT_ID,
                                                 Week,
                                                 Gender,
                                                 CompType,
                                                 Notification_type,
                                                 ends_with(format(CurrentDate, "%m%d"))) %>%
    rename_with(.fn = ~ str_remove(., paste(
      "_", format(CurrentDate, "%m%d"), sep = ""
    ))) %>%
    add_column(Date = CurrentDate, .after = "STUDY_PRTCPT_ID") %>%
    bind_rows(IndQ1dat_complete)
  HistoryMat <-
    HistoryMat %>% select(STUDY_PRTCPT_ID) %>% bind_cols(complete(mice_output, 1))
  
  ## 07-02
  CurrentDate = as.Date("2020-07-02")
  CurrentMat = tmpdat %>% filter(Date == CurrentDate) %>%
    select(-c(Date)) %>%
    rename_with(
      .fn = ~ paste(.x, format(CurrentDate, "%m%d"), sep = "_"),
      .cols = c(STEP_COUNT, SLEEP_COUNT, MOOD)
    )
  ImputeMat = CurrentMat %>% left_join(HistoryMat,by = "STUDY_PRTCPT_ID") %>% select(-c(STUDY_PRTCPT_ID,Week))
  mice_output <- mice(ImputeMat, m = 1, method = imputation_method)
  CompleteMatAll <- CurrentMat %>% select(STUDY_PRTCPT_ID, Week) %>% bind_cols(complete(mice_output, 1)) %>%
    mutate(CompType = "noncomp",Notification_type="nomsg")
  IndQ1dat_complete <- CompleteMatAll %>% select(STUDY_PRTCPT_ID,
                                                 Week,
                                                 Gender,
                                                 CompType,
                                                 Notification_type,
                                                 ends_with(format(CurrentDate, "%m%d"))) %>%
    rename_with(.fn = ~ str_remove(., paste(
      "_", format(CurrentDate, "%m%d"), sep = ""
    ))) %>%
    add_column(Date = CurrentDate, .after = "STUDY_PRTCPT_ID") %>%
    bind_rows(IndQ1dat_complete)
  HistoryMat <-
    HistoryMat %>% select(STUDY_PRTCPT_ID) %>% bind_cols(complete(mice_output, 1))
  
  ## 07-03
  CurrentDate = as.Date("2020-07-03")
  CurrentMat = tmpdat %>% filter(Date == CurrentDate) %>%
    select(-c(Date)) %>%
    rename_with(
      .fn = ~ paste(.x, format(CurrentDate, "%m%d"), sep = "_"),
      .cols = c(STEP_COUNT, SLEEP_COUNT, MOOD)
    )
  ImputeMat = CurrentMat %>% left_join(HistoryMat,by = "STUDY_PRTCPT_ID") %>% select(-c(STUDY_PRTCPT_ID,Week))
  mice_output <- mice(ImputeMat, m = 1, method = imputation_method)
  CompleteMatAll <- CurrentMat %>% select(STUDY_PRTCPT_ID, Week) %>% bind_cols(complete(mice_output, 1)) %>%
    mutate(CompType = "noncomp",Notification_type="nomsg")
  IndQ1dat_complete <- CompleteMatAll %>% select(STUDY_PRTCPT_ID,
                                                 Week,
                                                 Gender,
                                                 CompType,
                                                 Notification_type,
                                                 ends_with(format(CurrentDate, "%m%d"))) %>%
    rename_with(.fn = ~ str_remove(., paste(
      "_", format(CurrentDate, "%m%d"), sep = ""
    ))) %>%
    add_column(Date = CurrentDate, .after = "STUDY_PRTCPT_ID") %>%
    bind_rows(IndQ1dat_complete)
  HistoryMat <-
    HistoryMat %>% select(STUDY_PRTCPT_ID) %>% bind_cols(complete(mice_output, 1))
  
  ## 07-04 to 07-05
  days = seq.Date(as.Date("2020-07-04"),as.Date("2020-07-05"),"day")
  for (CurrentDate in as.list(days)){
    CurrentMat = tmpdat %>% filter(Date == CurrentDate) %>%
      select(-c(Date)) %>%
      rename_with(
        .fn = ~ paste(.x, format(CurrentDate, "%m%d"), sep = "_"),
        .cols = c(STEP_COUNT, SLEEP_COUNT, MOOD)
      )
    ImputeMat = CurrentMat %>% left_join(HistoryMat,by = "STUDY_PRTCPT_ID") %>% select(-c(STUDY_PRTCPT_ID,Week))
    mice_output <- mice(ImputeMat, m = 1, method = imputation_method)
    CompleteMatAll <- CurrentMat %>% select(STUDY_PRTCPT_ID, Week) %>% bind_cols(complete(mice_output, 1)) %>%
      mutate(CompType = "noncomp",Notification_type="nomsg")
    IndQ1dat_complete <- CompleteMatAll %>% select(STUDY_PRTCPT_ID,
                                                   Week,
                                                   Gender,
                                                   CompType,
                                                   Notification_type,
                                                   ends_with(format(CurrentDate, "%m%d"))) %>%
      rename_with(.fn = ~ str_remove(., paste(
        "_", format(CurrentDate, "%m%d"), sep = ""
      ))) %>%
      add_column(Date = CurrentDate, .after = "STUDY_PRTCPT_ID") %>%
      bind_rows(IndQ1dat_complete)
    HistoryMat <-
      CompleteMatAll %>% select(-ends_with(format(CurrentDate - 3, "%m%d"))) %>%
      select(-c(CompType, Notification_type, Week))
  }
  
  
  # impute trial data
  
  ## by competition group
  IndQ1dat <- IndQ1dat %>%
    left_join(
      competelist %>% rename(STUDY_PRTCPT_ID = AccessCode) %>% select(STUDY_PRTCPT_ID, TeamName)
    ) %>%
    left_join(GrpQ1dat) %>%
    left_join(IndMsgdat) %>%
    mutate(
      CompType = replace_na(CompType, "noncomp"),
      Notification_type = replace_na(Notification_type, "nomsg")
    )
  
  
  
  TrialDates = seq(as.Date("2020-07-06"), max_date, 1)
  
  for (CurrentDate in as.list(TrialDates)) {
    CurrentMatAll = IndQ1dat %>% filter(Date == CurrentDate) %>%
      select(
        STUDY_PRTCPT_ID,
        Week,
        STEP_COUNT,
        SLEEP_COUNT,
        MOOD,
        CompType,
        Notification_type,
        TeamName
      ) %>%
      rename_with(
        .fn = ~ paste(.x, format(CurrentDate, "%m%d"), sep = "_"),
        .cols = c(STEP_COUNT, SLEEP_COUNT, MOOD)
      ) %>%
      separate(TeamName,
               into = c("Institution", "Specialty"),
               sep = "_")
    CurrentWeek = CurrentMatAll$Week[1]
    if ((weekdays(CurrentDate) == "Monday")) {
      tmp = IndQ1dat_complete %>% filter(Week == CurrentWeek - 1) %>%
        group_by(STUDY_PRTCPT_ID) %>% summarise(
          PrevWeekStep = mean(STEP_COUNT),
          PrevWeekSleep = mean(SLEEP_COUNT),
          PrevWeekMood = mean(MOOD),
          PrevWeekCompType = unique(CompType)
        )
      HistoryMat <-
        HistoryMat %>% select(-starts_with("PrevWeek")) %>%
        left_join(tmp)
    }
    CompleteMatAll = tibble()
    subgroups = unique(CurrentMatAll[, c("CompType")])
    for (i in 1:dim(subgroups)[1]) {
      CurrentMat = CurrentMatAll %>% filter(CompType == unlist(subgroups[i, "CompType"]))
      ImputeMat = CurrentMat %>% inner_join(HistoryMat,by = "STUDY_PRTCPT_ID") %>% select(-c(STUDY_PRTCPT_ID, Week))
      mice_output <-
        mice(ImputeMat, m = 1, method = imputation_method)
      CompleteMatAll <-
        CompleteMatAll %>% bind_rows(CurrentMat %>% select(STUDY_PRTCPT_ID, Week) %>% bind_cols(complete(mice_output, 1)))
    }
    # modify full dataframe
    IndQ1dat_complete <- CompleteMatAll %>% select(STUDY_PRTCPT_ID,
                                                   Week,
                                                   Gender,
                                                   CompType,
                                                   Notification_type,
                                                   ends_with(format(CurrentDate, "%m%d"))) %>%
      rename_with(.fn = ~ str_remove(., paste(
        "_", format(CurrentDate, "%m%d"), sep = ""
      ))) %>%
      add_column(Date = CurrentDate, .after = "STUDY_PRTCPT_ID") %>%
      bind_rows(IndQ1dat_complete)
    
    # update history info
    HistoryMat <-
      CompleteMatAll %>% select(-ends_with(format(CurrentDate - 3, "%m%d"))) %>%
      select(-c(CompType, Notification_type, Week))
  }
  
  IndQ1dat_complete %>% left_join(
    competelist %>% rename(STUDY_PRTCPT_ID = AccessCode) %>% select(STUDY_PRTCPT_ID, TeamName)
  ) %>%
    left_join(GrpQ1dat) %>%
    arrange(STUDY_PRTCPT_ID, Date) -> IndQ1dat_complete
  impute_list[[m]] = list(
    "complete_dat" = IndQ1dat_complete,
    "original_dat" = IndQ1dat,
    "baseline_original_dat" = IndBLdat,
    "baseline_complete_dat" = IndBLdat_complete
  )
}

# save imputation setting
impute_list["seed"] = random_seed
impute_list["num_imputations"] = M
save(impute_list, file = "../imputation data/2020_Intern_competition_datasets_imputation.Rdata")



### supp codes
tmp = read_csv("../competition/2020 access code_All subjects.csv")
length(unique(tmp$USERID)) # 2286
x = unique(tmp$USERID) # 2286
length(setdiff(x, competelist$AccessCode)) # 350
x = intersect(x, competelist$AccessCode) # 1936
length(setdiff(x, IndQ1dat$STUDY_PRTCPT_ID)) # 139
x = intersect(x, IndQ1dat$STUDY_PRTCPT_ID) # 1797
#competelist %>% filter(AccessCode %in% x) %>% distinct(Institution,Program) %>% arrange(Institution,Program)
length(setdiff(x, IndSuvdat$STUDY_PRTCPT_ID)) # 3
x = intersect(x, IndSuvdat$STUDY_PRTCPT_ID) # 1794
length(setdiff(x, IndBLdat$STUDY_PRTCPT_ID)) # 15
x = intersect(x, IndBLdat$STUDY_PRTCPT_ID) # 1779

## specialty list
c(
  "Anesthesiology",
  "Child Neurology",
  "Child Psychiatry",
  "Dermatology",
  "Emergency Medicine",
  "Family Medicine",
  "Internal Medicine",
  "Interventional Radiology",
  "Neurology",
  "OBGYN",
  "Ophthalmology",
  "Otolaryngology",
  "Pathology",
  "Pediatrics",
  "Physical Medicine & Rehabilitation",
  "Primary Care",
  "Primary Medicine",
  "Psychiatry",
  "Radiation Oncology",
  "Radiology",
  "Surgery",
  "Transitional",
  "Urology"
) #23

tmp = read_csv("../competition/2020 access code_All subjects.csv")
length(unique(tmp$USERID)) # 2286
x = unique(tmp$USERID) # 2286
tmp = read_csv("../survey data/IHSdata_2020blq1q2_01042021.csv") %>% select(UserID,Sex)
tmp %>% filter(UserID %in% x) %>% select(Sex) %>% table()
