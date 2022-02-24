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
    citab[parm, ]
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


IndQ1dat <- impute_list[[1]]$original_dat
IndBLdat_complete <- impute_list[[1]]$baseline_complete_dat
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
        ReportMood = sum(!is.na(MOOD)) / length(MOOD)
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
    )

# complete case analysis
## primary aim
tmp = TeamQ1Weekdat
### step count
fit.ind = geeglm(
    WEEK_STEP_COUNT ~ Week + CompType + GenderMale + PreInternStep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
    data = tmp %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
    id = TeamName,
    corstr = "independence"
)
summary(fit.ind)


### sleep minute
fit.ind = geeglm(
    WEEK_SLEEP_COUNT ~ Week + CompType + GenderMale + PreInternSleep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
    data = tmp %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
    id = TeamName,
    corstr = "independence"
)
summary(fit.ind)

## Moderation analysis - week-in-study
tmp = TeamQ1Weekdat
### step count
fit.ind = geeglm(
    WEEK_STEP_COUNT ~ Week * CompType + GenderMale + PreInternStep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
    data = tmp %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
    id = TeamName,
    corstr = "independence"
)
summary(fit.ind)

### sleep minute
fit.ind = geeglm(
    WEEK_SLEEP_COUNT ~ Week * CompType + GenderMale + PreInternSleep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
    data = tmp %>% filter(CompType != "steps")  %>% filter(!is.na(PreSleep)),
    id = TeamName,
    corstr = "independence"
)
summary(fit.ind)

## Moderation analysis - intra-institution, intra-specialty
tmp = TeamQ1Weekdat

fit.ind = geeglm(
    WEEK_STEP_COUNT ~ CompWithinInstitution + CompWithinProgram + Week * CompType +
        GenderMale + PreInternStep + PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
    data = tmp %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
    family = gaussian(link = "identity"),
    id = TeamName,
    corstr = "independence"
)
summary(fit.ind)


fit.ind = geeglm(
    WEEK_SLEEP_COUNT ~ CompWithinInstitution + CompWithinProgram + Week * CompType +
        GenderMale + PreInternSleep + PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
    data = tmp %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
    family = gaussian(link = "identity"),
    id = TeamName,
    corstr = "independence"
)

## Exploratory analysis - participation rate
tmp = TeamQ1Weekdat %>% mutate(Comp = CompType != "noncomp")
fit.pr_step0 = geeglm(
    ReportStep ~ Week + Comp + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportStep,
    data = tmp %>% filter(!is.na(PreReportStep)),
    id = TeamName,
    corstr = "independence"
)
fit.pr_step = geeglm(
    ReportStep ~ Week * Comp + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportStep,
    data = tmp %>% filter(!is.na(PreReportStep)),
    id = TeamName,
    corstr = "independence"
)
fit.pr_sleep0 = geeglm(
    ReportSleep ~ Week + Comp + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportSleep,
    data = tmp %>% filter(!is.na(PreReportSleep)),
    id = TeamName,
    corstr = "independence"
)
fit.pr_sleep = geeglm(
    ReportSleep ~ Week * Comp + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportSleep,
    data = tmp %>% filter(!is.na(PreReportSleep)),
    id = TeamName,
    corstr = "independence"
)
fit.pr_mood0 = geeglm(
    ReportMood ~ Week + Comp + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportMood,
    data = tmp %>% filter(!is.na(PreReportMood)),
    id = TeamName,
    corstr = "independence"
)
fit.pr_mood = geeglm(
    ReportMood ~ Week * Comp + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportMood,
    data = tmp %>% filter(!is.na(PreReportMood)),
    id = TeamName,
    corstr = "independence"
)

## Exploratory analysis - mood score
tmp = TeamQ1Weekdat %>% mutate(Comp = CompType != "noncomp")
fit.ind = geeglm(
    WEEK_MOOD_SCORE ~ Week + Comp + GenderMale + PreInternMood +  PHQtot0 + Neu0 + depr0 + EFE0 + PreMood,
    data = tmp %>% filter(!is.na(PreMood)),
    id = TeamName,
    corstr = "independence"
)
summary(fit.ind)
fit.ind = geeglm(
    WEEK_MOOD_SCORE ~ Week * Comp + GenderMale + PreInternMood +  PHQtot0 + Neu0 + depr0 + EFE0 + PreMood,
    data = tmp %>% filter(!is.na(PreMood)),
    id = TeamName,
    corstr = "independence"
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

results = tibble()
for (i in 1:impute_list$num_imputations) {
    tmp = impute_list[[i]]$complete_dat %>%
        left_join(impute_list[[i]]$baseline_complete_dat)
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
            ReportMood = sum(!is.na(MOOD)) / length(MOOD)
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
        id = TeamName
    )
    fit.step = geeglm(
        WEEK_STEP_COUNT ~ Week * CompType + GenderMale + PreInternStep + PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
        data = tmp1 %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
        id = TeamName
    )
    fit.within_step = geeglm(
        WEEK_STEP_COUNT ~ Week * CompType + CompWithinInstitution + CompWithinProgram + GenderMale + PreInternStep + PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
        data = tmp1 %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),
        id = TeamName
    )
    
    fit.sleep0 = geeglm(
        WEEK_SLEEP_COUNT ~ Week + CompType + GenderMale + PreInternSleep + PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
        data = tmp1 %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
        id = TeamName
    )
    fit.sleep = geeglm(
        WEEK_SLEEP_COUNT ~ Week * CompType + GenderMale + PreInternSleep + PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
        data = tmp1 %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
        id = TeamName
    )
    fit.within_sleep = geeglm(
        WEEK_SLEEP_COUNT ~ Week * CompType + CompWithinInstitution + CompWithinProgram + GenderMale + PreInternSleep + PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
        data = tmp1 %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),
        id = TeamName
    )
    
    tmp1 = tmp1 %>% mutate(Comp = CompType != "noncomp")
    fit.mood0 = geeglm(
        WEEK_MOOD_SCORE ~ Week + Comp + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreInternStep + PreStep + PreInternSleep + PreSleep,
        data = tmp1 %>% filter(!is.na(PreStep)) %>% filter(!is.na(PreSleep)),
        id = TeamName
    )
    fit.mood = geeglm(
        WEEK_MOOD_SCORE ~ Week * Comp + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreInternStep + PreStep + PreInternSleep + PreSleep,
        data = tmp1 %>% filter(!is.na(PreStep)) %>% filter(!is.na(PreSleep)),
        id = TeamName
    )
    fits.cluster_step0[[i]] = fit.step0
    fits.cluster_step[[i]] = fit.step
    fits.cluster_intra_step[[i]] = fit.within_step
    fits.cluster_sleep0[[i]] = fit.sleep0
    fits.cluster_sleep[[i]] = fit.sleep
    fits.cluster_intra_sleep[[i]] = fit.within_sleep
    fits.cluster_mood0[[i]] = fit.mood0
    fits.cluster_mood[[i]] = fit.mood
    
}

summary(pool(fits.cluster_step0), conf.int = TRUE)
summary(pool(fits.cluster_step), conf.int = TRUE)
summary(pool(fits.cluster_intra_step), conf.int = TRUE)
summary(pool(fits.cluster_sleep0), conf.int = TRUE)
summary(pool(fits.cluster_sleep), conf.int = TRUE)
summary(pool(fits.cluster_intra_sleep), conf.int = TRUE)
summary(pool(fits.cluster_mood0), conf.int = TRUE)
summary(pool(fits.cluster_mood), conf.int = TRUE)
#


## do you like painting
library(ggpubr)
plot_ci_bound = function(fits,var_name,title,mi=TRUE){
    # https://www.rdocumentation.org/packages/miceadds/versions/3.11-6/topics/pool_mi
    if (mi == TRUE){
        require("mitools")
        require("miceadds")
        cmod = MIextract(fits, fun = coef)
        vmod = MIextract(fits, fun = vcov)
        res <- miceadds::pool_mi(qhat = cmod, u = vmod)
    }
    else{
        res = fits
    }
    res.coef = coef(res)
    res.vcov = vcov(res)
    
    interaction_term = paste("Week:",var_name,sep="")
    p_est=unname(sapply(
        0:11,
        FUN = function(i)
            res.coef[var_name] + i * res.coef[interaction_term]
    ))
    std_est = unname(sapply(
        0:11,
        FUN = function(i)
            sqrt(res.vcov[var_name, var_name] + i ** 2 * res.vcov[interaction_term, interaction_term] + 2 *
                     i * res.vcov[var_name, interaction_term])
    ))
    gdat = data.frame(week=1:12,est=p_est,lwr=p_est-1.96*std_est,upr=p_est+1.96*std_est)
    ggplot(data = gdat,aes(x = week,y=est)) +
        theme_bw() +
        geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) + 
        geom_line(size=0.8) + geom_hline(yintercept = 0,color = "red", linetype="dashed") +
        xlab("Week") + ylab(title) + 
        scale_x_continuous(breaks = 1:12,labels = 1:12) + 
        theme(text=element_text(size=20))
}

# week-in-study analysis
p1 <- plot_ci_bound(fits.cluster_step,"CompTypesteps","Competition step effect on weekly average daily step count")
p2 <- plot_ci_bound(fits.cluster_sleep,"CompTypesleep","Competition sleep effect on weekly average daily sleep minutes")
png(file="C://Users/jtwan/Downloads/week_in_study.png",width=1500, height=750)
ggarrange(p1,p2,ncol = 2,labels = c("a","b"))
dev.off()

# participantion rate analysis
p1 <- plot_ci_bound(fit.pr_step,"CompTRUE","Competition effect on weekly average daily step count participation rate",mi=FALSE)
p2 <- plot_ci_bound(fit.pr_sleep,"CompTRUE","Competition effect on weekly average daily sleep minutes participation rate",mi=FALSE)
p3 <- plot_ci_bound(fit.pr_mood,"CompTRUE","Competition effect on weekly average daily mood score participation rate",mi=FALSE)
png(file="C://Users/jtwan/Downloads/participation_rate.png",width=1500, height=1500)
ggarrange(p1,p2,p3,ncol = 2,nrow=2,labels = c("a","b","c"))
dev.off()
