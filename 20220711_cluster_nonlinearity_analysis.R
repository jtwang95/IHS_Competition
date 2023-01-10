library(geepack)
library(ggplot2)
library(tidyverse)
library(sjPlot)
library(mice)
library(mgcv)
library(parallel)
library(splines)
library(ggpubr)
library(ClusterBootstrap)

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
my_clusbootglm <- function (model, data, clusterid, family = gaussian, B = 5000, 
                            confint.level = 0.95, n.cores = 1) 
{
  tt_cores <- detectCores()
  if (n.cores > tt_cores) {
    message(sprintf("Note: \"n.cores\" was set to %d, but only %d are available. Using all cores.", 
                    n.cores, tt_cores))
  }
  model <- as.formula(model)
  res.or <- glm(model, family = family, data = data,weights = cluster_wt)
  confint.pboundaries = c((1 - confint.level)/2, 1 - (1 - confint.level)/2)
  confint.Zboundaries = qnorm(confint.pboundaries)
  n <- nrow(data)
  p <- length(res.or$coef)
  coefs <- matrix(NA, nrow = B, ncol = p)
  arguments <- as.list(match.call())
  clusterid <- eval(arguments$clusterid, data)
  cluster <- as.character(clusterid)
  clusters <- unique(cluster)
  Obsno <- split(1:n, cluster)
  f = matrix(clusters, length(clusters), B)
  ff = matrix(f, prod(dim(f)), 1)
  fff = sample(ff)
  f = matrix(fff, length(clusters), B)
  if (is.numeric(n.cores) & n.cores > 0) {
    if (n.cores == 1) {
      for (i in 1:B) {
        j <- f[, i]
        obs <- unlist(Obsno[j])
        bootcoef <- tryCatch(coef(glm(model, family = family, 
                                      data = data[obs, ])), warning = function(x) rep(as.numeric(NA), 
                                                                                      p))
        coefs[i, which(names(res.or$coef) %in% names(bootcoef))] <- bootcoef
      }
    }
    if (n.cores > 1) {
      cl <- makeCluster(max(min(n.cores, tt_cores),2))
      previous_RNGkind <- RNGkind()[1]
      RNGkind("L'Ecuyer-CMRG")
      nextRNGStream(.Random.seed)
      clusterExport(cl, varlist = c("f", "Obsno", 
                                    "model", "family", "data", 
                                    "p", "res.or", "clusbootglm_sample_glm"), 
                    envir = environment())
      splitclusters <- 1:B
      out <- parSapplyLB(cl, splitclusters, function(x) clusbootglm_sample_glm(f, 
                                                                               x, Obsno, model, family, data, p, res.or))
      coefs <- t(out)
      stopCluster(cl)
      RNGkind(previous_RNGkind)
    }
  }
  invalid.samples <- colSums(is.na(coefs))
  names(invalid.samples) <- colnames(coefs) <- names(res.or$coef)
  samples.with.NA.coef <- which(is.na(rowSums(coefs)))
  sdcoefs <- apply(coefs, 2, sd, na.rm = TRUE)
  result <- list(call = match.call(), model = model, family = family, 
                 B = B, coefficients = coefs, data = data, bootstrap.matrix = f, 
                 subject.vector = clusterid, lm.coefs = res.or$coef, boot.coefs = colMeans(coefs, 
                                                                                           na.rm = TRUE), boot.sds = sdcoefs, ci.level = confint.level)
  class(result) <- "clusbootglm"
  return(result)
}
clusbootglm_sample_glm <-function(f, i, Obsno, model, family, data, p, res.or){
  j <- f[, i]
  obs <- unlist(Obsno[j])
  coef <- rep(NA,p) #added
  bootcoef <- tryCatch(coef(glm(model, family = family, data = data[obs,],weights = cluster_wt)),
                       warning=function(x) rep(as.numeric(NA),p))
  ifelse(length(bootcoef)==p, coef <- as.vector(bootcoef), coef[which(names(res.or$coef) %in% names(bootcoef))] <- bootcoef)
  return(coef)
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

Indweight = read_csv("C:/Users/jtwan/Dropbox (University of Michigan)/IHS_competition/competition analysis/revision1/subject_id_weights.csv")
IndQ1dat <- impute_list[[1]]$original_dat
IndBLdat_complete <- impute_list[[1]]$baseline_complete_dat %>%
  inner_join(Indweight,by=c("STUDY_PRTCPT_ID" = "UserID"))
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
  )

# complete case analysis

## Moderation analysis - week-in-study
tmp = TeamQ1Weekdat
#############################################
############### step count ##################
#############################################
fit.gam = gam(
  WEEK_STEP_COUNT ~CompType + s(Week) + s(Week,by=factor(CompType))+ GenderMale + PreInternStep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreStep,
  data = tmp %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),method = "REML",
  weights = cluster_wt
)
summary(fit.gam)

pdat <- expand.grid(Week = seq(0, 11),
                    GenderMale = c(0.5),
                    CompType = c("noncomp","steps"),PreInternStep = 0,
                    PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreStep=0)
xp <- predict(fit.gam, newdata = pdat, type = 'lpmatrix')
c1 <- grepl('Comp', colnames(xp))
r1 <- with(pdat, CompType == "steps")
r0 <- with(pdat, CompType == "noncomp")
X <- xp[r1, ] - xp[r0, ]
X[, !c1] <- 0
dif <- X %*% coef(fit.gam)
se <- sqrt(diag(X %*% vcov(fit.gam) %*% t(X)))
crit <- qt(.975, df.residual(fit.gam))
upr <- dif + (crit * se)
lwr <- dif - (crit * se)

p1 <- ggplot(data.frame(week = seq(0,11), effect = dif, effect_lb = lwr, effect_ub = upr), 
             aes(x = week, y = effect)) +
  theme_bw() + 
  geom_ribbon(aes(ymin = effect_lb, ymax = effect_ub), alpha = 0.2) +
  geom_line(size=0.8) + 
  geom_hline(yintercept = 0, color = "red", linetype="dashed")  +
  xlab("Week") + ylab("Effect of competition step on weekly average daily step count") +
  scale_x_continuous(breaks = 0:11,labels = 1:12) + theme(text=element_text(size=20)) +
  ylim(c(-125,350))


spline_model_formula = as.formula(WEEK_STEP_COUNT ~ ns(Week,knots = c(6))*CompType + GenderMale + PreInternStep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreStep)
test = my_clusbootglm(spline_model_formula,data =tmp %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),clusterid = TeamName,B=100)
calculate_effect_and_ci_bootstrap = function(i,test){
  stopifnot(class(test) == "clusbootglm")
  tmp = clusbootsample(test,i)
  fit.boots = glm(spline_model_formula,data = tmp)
  effects = sapply(0:11, function(j) predict(fit.boots,list( GenderMale =0.5,Week=j,CompType="steps",
                                                             PreInternStep = 0,
                                                             PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreStep=0)) -
                     predict(fit.boots,list(GenderMale =0.5,Week=j,CompType="noncomp",
                                            PreInternStep = 0,
                                            PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreStep=0)))
  return(effects)
}

effects = t(sapply(1:dim(test$coefficients)[1],calculate_effect_and_ci_bootstrap, test = test))
fit.boots = glm(spline_model_formula,data = tmp %>% filter(CompType != "sleep") %>% filter(!is.na(PreStep)),family = gaussian(link = "identity"))
effect = sapply(0:11,function(j) predict(fit.boots,list(GenderMale =0.5,Week=j,CompType="steps",
                                                        PreInternStep = 0,
                                                        PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreStep=0)) -
                  predict(fit.boots,list(GenderMale =0.5,Week=j,CompType="noncomp",
                                         PreInternStep = 0,
                                         PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreStep=0)))
res = data.frame(week=0:11,
                 effect = effect,
                 effect_lb = apply(effects,2,function(x) quantile(x,0.025)),
                 effect_ub = apply(effects,2,function(x) quantile(x,0.975)))

p2<- ggplot(data = res,aes(x = week,y=effect)) +
  theme_bw() +
  geom_ribbon(aes(ymin = effect_lb, ymax = effect_ub), alpha = 0.2) + 
  geom_line(size=1) + geom_hline(yintercept = 0,color = "red", linetype="dashed") +
  xlab("Week") + ylab("Effect of competition step on weekly average daily step count") + 
  scale_x_continuous(breaks = 0:11,labels = 1:12) + theme(text=element_text(size=20)) + 
  ylim(c(-125,350))

#############################################
############### sleep mins ##################
#############################################
fit.gam = gam(
  WEEK_SLEEP_COUNT ~CompType + s(Week) + s(Week,by=factor(CompType))+ GenderMale + PreInternSleep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep,
  data = tmp %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),method="REML",weights = cluster_wt)
summary(fit.gam)

pdat <- expand.grid(Week = seq(0, 11),
                    GenderMale = c(0.5),
                    CompType = c("noncomp","sleep"),PreInternSleep = 0,
                    PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreSleep=0)
xp <- predict(fit.gam, newdata = pdat, type = 'lpmatrix')
c1 <- grepl('Comp', colnames(xp))
r1 <- with(pdat, CompType == "sleep")
r0 <- with(pdat, CompType == "noncomp")
X <- xp[r1, ] - xp[r0, ]
X[, !c1] <- 0
dif <- X %*% coef(fit.gam)
se <- sqrt(diag(X %*% vcov(fit.gam) %*% t(X)))
crit <- qt(.975, df.residual(fit.gam))
upr <- dif + (crit * se)
lwr <- dif - (crit * se)

p3 <- ggplot(data.frame(week = seq(0,11), effect = dif, effect_lb = lwr, effect_ub = upr), 
             aes(x = week, y = effect)) +
  theme_bw() + 
  geom_ribbon(aes(ymin = effect_lb, ymax = effect_ub), alpha = 0.2) +
  geom_line(size=0.8) + 
  geom_hline(yintercept = 0, color = "red", linetype="dashed")  +
  xlab("Week") + ylab("Effect of competition sleep on weekly average daily sleep minutes") +
  scale_x_continuous(breaks = 0:11,labels = 1:12)+theme(text=element_text(size=20)) + 
  ylim(c(-20,15))


spline_model_formula = as.formula(WEEK_SLEEP_COUNT ~ ns(Week,knots = c(6))*CompType + GenderMale + PreInternSleep +  PHQtot0 + Neu0 + depr0 + EFE0 + PreSleep)
test = my_clusbootglm(spline_model_formula,data =tmp %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),clusterid = TeamName,B=100)
calculate_effect_and_ci_bootstrap = function(i,test){
  stopifnot(class(test) == "clusbootglm")
  tmp = clusbootsample(test,i)
  fit.boots = glm(spline_model_formula,data = tmp)
  effects = sapply(0:11, function(j) predict(fit.boots,list( GenderMale =0.5,Week=j,CompType="sleep",
                                                             PreInternSleep = 0,
                                                             PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreSleep=0)) -
                     predict(fit.boots,list(GenderMale =0.5,Week=j,CompType="noncomp",
                                            PreInternSleep = 0,
                                            PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreSleep=0)))
  return(effects)
}

effects = t(sapply(1:dim(test$coefficients)[1],calculate_effect_and_ci_bootstrap, test = test))
fit.boots = glm(spline_model_formula,data = tmp %>% filter(CompType != "steps") %>% filter(!is.na(PreSleep)),family = gaussian(link = "identity"))
effect = sapply(0:11,function(j) predict(fit.boots,list(GenderMale =0.5,Week=j,CompType="sleep",
                                                        PreInternSleep = 0,
                                                        PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreSleep=0)) -
                  predict(fit.boots,list(GenderMale =0.5,Week=j,CompType="noncomp",
                                         PreInternSleep = 0,
                                         PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreSleep=0)))
res = data.frame(week=0:11,
                 effect = effect,
                 effect_lb = apply(effects,2,function(x) quantile(x,0.025)),
                 effect_ub = apply(effects,2,function(x) quantile(x,0.975)))

p4<- ggplot(data = res,aes(x = week,y=effect)) +
  theme_bw() +
  geom_ribbon(aes(ymin = effect_lb, ymax = effect_ub), alpha = 0.2) + 
  geom_line(size=1) + geom_hline(yintercept = 0,color = "red", linetype="dashed") +
  xlab("Week") + ylab("Effect of competition sleep on weekly average daily sleep minutes") + 
  scale_x_continuous(breaks = 0:11,labels = 1:12)+ theme(text=element_text(size=20)) + 
  ylim(c(-20,15))

png(file="C://Users/jtwan/Downloads/nonlinear_week_in_study.png",width=1500, height=1500)
ggarrange(p1,p3,p2,p4,nrow=2,ncol=2,labels = c("a","b","c","d"))
dev.off()


## Exploratory analysis - participation rate
tmp = TeamQ1Weekdat %>% mutate(Comp = CompType != "noncomp")
#############################################################
################### step participation rate #################
#############################################################
fit.gam = gam(
  ReportStep ~ s(Week) + Comp + s(Week,by=factor(Comp)) + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportStep,
  data = tmp %>% filter(!is.na(PreReportStep)),method="REML",weights = cluster_wt)
summary(fit.gam)

pdat <- expand.grid(Week = seq(0, 11),
                    GenderMale = c(0.5),
                    Comp = c(TRUE,FALSE),
                    PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportStep=0)
xp <- predict(fit.gam, newdata = pdat, type = 'lpmatrix')
c1 <- grepl('Comp', colnames(xp))
r1 <- with(pdat, Comp == TRUE)
r0 <- with(pdat, Comp == FALSE)
X <- xp[r1, ] - xp[r0, ]
X[, !c1] <- 0
dif <- X %*% coef(fit.gam)
se <- sqrt(diag(X %*% vcov(fit.gam) %*% t(X)))
crit <- qt(.975, df.residual(fit.gam))
upr <- dif + (crit * se)
lwr <- dif - (crit * se)

p1 <- ggplot(data.frame(week = seq(0,11), effect = dif, effect_lb = lwr, effect_ub = upr), 
             aes(x = week, y = effect)) +
  theme_bw() + 
  geom_ribbon(aes(ymin = effect_lb, ymax = effect_ub), alpha = 0.2) +
  geom_line(size=0.8) + 
  geom_hline(yintercept = 0, color = "red", linetype="dashed")  +
  xlab("Week") + ylab("Effect of competition on step count participation rate") +
  scale_x_continuous(breaks = 0:11,labels = 1:12)+theme(text=element_text(size=18))


spline_model_formula = as.formula(ReportStep~ ns(Week,knots = c(6))*factor(Comp) +  GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportStep)
test = my_clusbootglm(spline_model_formula,data =tmp %>% filter(!is.na(PreReportStep)),clusterid = TeamName,B=100)
calculate_effect_and_ci_bootstrap = function(i,test){
  stopifnot(class(test) == "clusbootglm")
  tmp = clusbootsample(test,i)
  fit.boots = glm(spline_model_formula,data = tmp)
  effects = sapply(0:11, function(j) predict(fit.boots,list( GenderMale =0.5,Week=j,Comp=TRUE,
                                                             PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportStep=0)) -
                     predict(fit.boots,list(GenderMale =0.5,Week=j,Comp=FALSE,
                                            PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportStep=0)))
  return(effects)
}

effects = t(sapply(1:dim(test$coefficients)[1],calculate_effect_and_ci_bootstrap, test = test))
fit.boots = glm(spline_model_formula,data = tmp %>% filter(!is.na(PreReportStep)),family = gaussian(link = "identity"))
effect = sapply(0:11,function(j) predict(fit.boots,list(GenderMale =0.5,Week=j,Comp=TRUE,
                                                        PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportStep=0)) -
                  predict(fit.boots,list(GenderMale =0.5,Week=j,Comp=FALSE,
                                         PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportStep=0)))
res = data.frame(week=0:11,
                 effect = effect,
                 effect_lb = apply(effects,2,function(x) quantile(x,0.025)),
                 effect_ub = apply(effects,2,function(x) quantile(x,0.975)))

p2<- ggplot(data = res,aes(x = week,y=effect)) +
  theme_bw() +
  geom_ribbon(aes(ymin = effect_lb, ymax = effect_ub), alpha = 0.2) + 
  geom_line(size=1) + geom_hline(yintercept = 0,color = "red", linetype="dashed") +
  xlab("Week") + ylab("Effect of competition sleep on weekly average daily sleep minutes") + 
  scale_x_continuous(breaks = 0:11,labels = 1:12)+ theme(text=element_text(size=18))

#############################################################
################### sleep participation rate #################
#############################################################
fit.gam = gam(
  ReportSleep ~ s(Week) + Comp + s(Week,by=factor(Comp)) + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportSleep,
  data = tmp %>% filter(!is.na(PreReportSleep)),method="REML",weights = cluster_wt)
summary(fit.gam)

pdat <- expand.grid(Week = seq(0, 11),
                    GenderMale = c(0.5),
                    Comp = c(TRUE,FALSE),
                    PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportSleep=0)
xp <- predict(fit.gam, newdata = pdat, type = 'lpmatrix')
c1 <- grepl('Comp', colnames(xp))
r1 <- with(pdat, Comp == TRUE)
r0 <- with(pdat, Comp == FALSE)
X <- xp[r1, ] - xp[r0, ]
X[, !c1] <- 0
dif <- X %*% coef(fit.gam)
se <- sqrt(diag(X %*% vcov(fit.gam) %*% t(X)))
crit <- qt(.975, df.residual(fit.gam))
upr <- dif + (crit * se)
lwr <- dif - (crit * se)

p3 <- ggplot(data.frame(week = seq(0,11), effect = dif, effect_lb = lwr, effect_ub = upr), 
             aes(x = week, y = effect)) +
  theme_bw() + 
  geom_ribbon(aes(ymin = effect_lb, ymax = effect_ub), alpha = 0.2) +
  geom_line(size=0.8) + 
  geom_hline(yintercept = 0, color = "red", linetype="dashed")  +
  xlab("Week") + ylab("Effect of competition on sleep minutes participation rate") +
  scale_x_continuous(breaks = 0:11,labels = 1:12)+theme(text=element_text(size=18))


spline_model_formula = as.formula(ReportSleep~ ns(Week,knots = c(6))*factor(Comp) +  GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportSleep)
test = my_clusbootglm(spline_model_formula,data =tmp %>% filter(!is.na(PreReportSleep)),clusterid = TeamName,B=100)
calculate_effect_and_ci_bootstrap = function(i,test){
  stopifnot(class(test) == "clusbootglm")
  tmp = clusbootsample(test,i)
  fit.boots = glm(spline_model_formula,data = tmp)
  effects = sapply(0:11, function(j) predict(fit.boots,list( GenderMale =0.5,Week=j,Comp=TRUE,
                                                             PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportSleep=0)) -
                     predict(fit.boots,list(GenderMale =0.5,Week=j,Comp=FALSE,
                                            PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportSleep=0)))
  return(effects)
}

effects = t(sapply(1:dim(test$coefficients)[1],calculate_effect_and_ci_bootstrap, test = test))
fit.boots = glm(spline_model_formula,data = tmp %>% filter(!is.na(PreReportSleep)),family = gaussian(link = "identity"))
effect = sapply(0:11,function(j) predict(fit.boots,list(GenderMale =0.5,Week=j,Comp=TRUE,
                                                        PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportSleep=0)) -
                  predict(fit.boots,list(GenderMale =0.5,Week=j,Comp=FALSE,
                                         PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportSleep=0)))
res = data.frame(week=0:11,
                 effect = effect,
                 effect_lb = apply(effects,2,function(x) quantile(x,0.025)),
                 effect_ub = apply(effects,2,function(x) quantile(x,0.975)))

p4<- ggplot(data = res,aes(x = week,y=effect)) +
  theme_bw() +
  geom_ribbon(aes(ymin = effect_lb, ymax = effect_ub), alpha = 0.2) + 
  geom_line(size=1) + geom_hline(yintercept = 0,color = "red", linetype="dashed") +
  xlab("Week") + ylab("Effect of competition on sleep minutes participation rate") + 
  scale_x_continuous(breaks = 0:11,labels = 1:12)+ theme(text=element_text(size=18))

#############################################################
################### mood participation rate #################
#############################################################
fit.gam = gam(
  ReportMood ~ s(Week) + Comp + s(Week,by=factor(Comp)) + GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportMood,
  data = tmp %>% filter(!is.na(PreReportMood)),method="REML",weights = cluster_wt)
summary(fit.gam)

pdat <- expand.grid(Week = seq(0, 11),
                    GenderMale = c(0.5),
                    Comp = c(TRUE,FALSE),
                    PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportMood=0)
xp <- predict(fit.gam, newdata = pdat, type = 'lpmatrix')
c1 <- grepl('Comp', colnames(xp))
r1 <- with(pdat, Comp == TRUE)
r0 <- with(pdat, Comp == FALSE)
X <- xp[r1, ] - xp[r0, ]
X[, !c1] <- 0
dif <- X %*% coef(fit.gam)
se <- sqrt(diag(X %*% vcov(fit.gam) %*% t(X)))
crit <- qt(.975, df.residual(fit.gam))
upr <- dif + (crit * se)
lwr <- dif - (crit * se)

p5 <- ggplot(data.frame(week = seq(0,11), effect = dif, effect_lb = lwr, effect_ub = upr), 
             aes(x = week, y = effect)) +
  theme_bw() + 
  geom_ribbon(aes(ymin = effect_lb, ymax = effect_ub), alpha = 0.2) +
  geom_line(size=0.8) + 
  geom_hline(yintercept = 0, color = "red", linetype="dashed")  +
  xlab("Week") + ylab("Effect of competition on mood survey participation rate") +
  scale_x_continuous(breaks = 0:11,labels = 1:12)+theme(text=element_text(size=18))


spline_model_formula = as.formula(ReportMood~ ns(Week,knots = c(6))*factor(Comp) +  GenderMale + PHQtot0 + Neu0 + depr0 + EFE0 + PreReportMood)
test = my_clusbootglm(spline_model_formula,data =tmp %>% filter(!is.na(PreReportMood)),clusterid = TeamName,B=100)
calculate_effect_and_ci_bootstrap = function(i,test){
  stopifnot(class(test) == "clusbootglm")
  tmp = clusbootsample(test,i)
  fit.boots = glm(spline_model_formula,data = tmp)
  effects = sapply(0:11, function(j) predict(fit.boots,list( GenderMale =0.5,Week=j,Comp=TRUE,
                                                             PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportMood=0)) -
                     predict(fit.boots,list(GenderMale =0.5,Week=j,Comp=FALSE,
                                            PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportMood=0)))
  return(effects)
}

effects = t(sapply(1:dim(test$coefficients)[1],calculate_effect_and_ci_bootstrap, test = test))
fit.boots = glm(spline_model_formula,data = tmp %>% filter(!is.na(PreReportMood)),family = gaussian(link = "identity"))
effect = sapply(0:11,function(j) predict(fit.boots,list(GenderMale =0.5,Week=j,Comp=TRUE,
                                                        PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportMood=0)) -
                  predict(fit.boots,list(GenderMale =0.5,Week=j,Comp=FALSE,
                                         PHQtot0=0,Neu0=0,depr0=0,EFE0=0,PreReportMood=0)))
res = data.frame(week=0:11,
                 effect = effect,
                 effect_lb = apply(effects,2,function(x) quantile(x,0.025)),
                 effect_ub = apply(effects,2,function(x) quantile(x,0.975)))

p6<- ggplot(data = res,aes(x = week,y=effect)) +
  theme_bw() +
  geom_ribbon(aes(ymin = effect_lb, ymax = effect_ub), alpha = 0.2) + 
  geom_line(size=1) + geom_hline(yintercept = 0,color = "red", linetype="dashed") +
  xlab("Week") + ylab("Effect of competition on mood survey participation rate") + 
  scale_x_continuous(breaks = 0:11,labels = 1:12)+ theme(text=element_text(size=18))


png(file="C://Users/jtwan/Downloads/nonlinear_participation.png",width=1500, height=1000)
ggarrange(p1,p3,p5,p2,p4,p6,nrow=2,ncol=3,labels = c("a","b","c","d","e","f"))
dev.off()



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