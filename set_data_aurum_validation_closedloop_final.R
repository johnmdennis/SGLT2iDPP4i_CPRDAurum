# load libraries
rm(list=ls())
require(tidyverse)
library(rms)
library(tableone)
#library(plyr)
library(ggthemes)
library(forestplot)
library(broom)
library(ggplotify)
library(patchwork)
library(scales)
library(forcats)
library(marginaleffects)

output_dir <- "C:/Users/jmd237/OneDrive - University of Exeter/John/Projects/2023_SGLT2DPP4Aurum/results/" 
data_dir <- "C:/Users/jmd237/OneDrive - University of Exeter/John/CPRD/mastermind22/"

#### Data prep (run once only) ####

#Adapt from Pedros code (see Projects\2019_SGLT2vsDPP4\scripts\validation\set_data_validation.R)
dataset.type <- "diagnostics"

set_up_data_sglt2_dpp4 <- function(dataset.type) {
  ##### Input variables
  # dataset.type: a character string mentioning the type of dataset required
  
  # initial checks
  if (missing(dataset.type)) {stop("'dataset.type' needs to be supplied")}
  if (!is.character(dataset.type)) {stop("'dataset.type' must be a character string")}
  if (!(dataset.type %in% c("diagnostics", "synthetic", "full.cohort", "ps.model.train", "ps.model.test", "hba1c.train", "hba1c.test", "weight.dataset", "discontinuation.dataset", "egfr.dataset"))) {
    stop("'dataset.type' must one of: diagnostics / synthetic / full.cohort / ps.model.train / ps.model.test / hba1c.train / hba1c.test / weight.dataset / discontinuation.dataset / egfr.dataset")
  }}


# load original dataset # name - t2d_1stinstance
#load("/slade/CPRD_data/mastermind_2022/20221110_t2d_1stinstance.Rda")
load(paste0(data_dir,"20221212_t2d_1stinstance.Rda"))

cprd <- t2d_1stinstance

##### Select only SGLT-2 and DPP4

#SGLT2 vs DPP4
cprd <- cprd %>% 
  filter(drugclass == "DPP4" | drugclass == "SGLT2")

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### Select only SGLT-2 and GLP-1")
  print("################################################")
  print(nrow(cprd))
  print(table(cprd$drugclass))
  
}


##### Drop patients initiating before 1/1/2013

# Explore adjusted HbA1c repsonse by calendar year

cprd  <- cprd %>%
  mutate(yrdrugstart = format(dstartdate, format = "%Y")) %>%
  mutate(yrdrugstart = as.numeric(yrdrugstart))


cprd <- cprd %>%
  mutate(dstartdate_cutoff = ifelse(dstartdate < "2013-01-01", 1 , NA_real_))

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### Drop patients initiating before 1/1/2013")
  print("################################################")
  print(table(cprd$dstartdate_cutoff))
  print(table(cprd$dstartdate_cutoff, cprd$drugclass))
  
}

cprd <- cprd %>%
  filter(is.na(dstartdate_cutoff))


##### Drop if treated with insulin when starting new drug

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### Drop if treated with insulin when starting new drug")
  print("################################################")
  print(table(cprd$INS))
  print(table(cprd$INS, cprd$drugclass))
  
}

cprd <- cprd %>% 
  filter(INS == 0)      

##### Drop patients with ESRD

describe(cprd$preegfr)

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### Drop patients with ESRD")
  print("################################################")
  print(table(cprd$preckdstage))
  print(table(cprd$preckdstage, cprd$drugclass))
  
}

cprd <- cprd %>%
  filter(preckdstage != "stage_5" & preckdstage != "stage_4" & preckdstage != "stage_3b")

#describe(cprd$preegfr)

cprd <- cprd %>% mutate(egfrl45 = ifelse(preegfr<45,1,NA))
#describe(cprd$egfrl45)

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### Drop patients with EGFR<45")
  print("################################################")
  print(table(cprd$egfrl45))
  print(table(cprd$egfrl45, cprd$drugclass))
  
}  

cprd <- cprd %>%
  filter(is.na(egfrl45))
#describe(cprd$preegfr)

##### Drop if first-line treatment

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### Drop if first-line treatment")
  print("################################################")
  print(table(cprd$drugline_all))
  print(table(cprd$drugline_all, cprd$drugclass))
  
}

cprd <- cprd %>%
  filter(drugline_all != 1)


##### Drop if semaglutide

# cprd <- cprd %>%
#   mutate(semaglutide_drug = ifelse(str_detect(drugsubstances, "Semaglutide"), 1, NA_real_))
# 
# # printing inclusion patients
# if (dataset.type == "diagnostics") {
#   
#   print("################################################")
#   print("##### Drop if semaglutide")
#   print("################################################")
#   print(table(cprd$semaglutide_drug))
#   
# }
# 
# cprd <- cprd %>%
#   filter(is.na(semaglutide_drug))

############################# Variable Prep 

### Add variable that identifies an individual entry in the data

cprd <- cprd %>%
  mutate(pated = paste(patid, drugclass, dstartdate, sep = ".")) #%>%

##### Drug of interest

#####   - DPP4 and SGLT2: drugclass
#mutate(drugclass = factor(drugclass)) 
# 1 - DPP4; 2 - SGLT2

##### Outcome HbA1c # name: posthba1c12m (missing - 46383)
#####   - posthba1c6m but if missing
#####     - posthba1c12m

cprd <- cprd %>%
  mutate(posthba1cfinal = ifelse(is.na(posthba1c6m), posthba1c12m, posthba1c6m)) %>%
  mutate(posthba1cfinal = as.numeric(posthba1cfinal))

##### Outcome Weight # name: 
#####   - postweight6m but if missing
#####     - postweight12m

cprd <- cprd %>%
  mutate(postweightfinal = ifelse(is.na(postweight6m), postweight12m, postweight6m)) %>%
  mutate(postweightfinal = as.numeric(postweightfinal))

##### Sociodemographic variables

cprd <- cprd %>%
  #####   - Age: agetx (new var)
  mutate(agetx = as.numeric(dstartdate_age)) %>%
  #####   - Sex: sex
  mutate(sex = factor(ifelse(gender == 1, "Male", "Female"))) %>%
  
  #####   - Duration of diabetes: t2dmduration
  mutate(t2dmduration = as.numeric(dstartdate_dm_dur_all)) %>%
  #####   - Ethnicity: ethnicity
  mutate(ethnicity_5cat = ifelse(is.na(ethnicity_5cat),99,ethnicity_5cat)) %>% 
  mutate(ethnicity = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4, 99), labels = c("White", "South Asian", "Black", "Other", "Mixed", "Missing"))) %>%
  mutate(ethnicity_16cat = ifelse(is.na(ethnicity_16cat),99,ethnicity_16cat)) %>% 
  mutate(ethnicity16 = factor(ethnicity_16cat, levels = c(1:16, 99), 
                              labels = c("British","Irish","Other White","White & Black Caribbean",
                                         "White & Black African", "White & Asian", "Other Mixed", "Indian", "Pakastani",
                                         "Bangladeshi","Other Asian","Caribbean","African","Other Black","Chinese",
                                         "Other Ethnic Group","Missing"))) %>%
  #####   - Deprivation: deprivation
  mutate(deprivation = factor(imd2015_10)) %>%
  
  #####   - Smoking Status: smoke
  mutate(smoke = factor(smoking_cat)) %>%
  
  #####   - Line Therapy: drugline: turn all > 4 to 5+
  mutate(drugline = ifelse(drugline_all > 4, 5, drugline_all)) %>% 
  mutate(drugline.t = as.factor(as.numeric(drugline))) %>%
  # mutate(drugline = factor(drugline, levels = c(2, 3, 4, 5), labels = c("2", "3", "4", "5+"))) %>%
  
  #####   - Hospitalisations in previous year
  mutate(prehospitalisation = factor(hosp_admission_prev_year, levels = c(0, 1), labels = c("No", "Yes")))

##### Diabetes treatment

cprd <- cprd %>%
  #####   - Drugs taken alongside treatment
  # #####     - SU
  # #####     - MFN
  # #####     - DPP4
  # #####     - TZD
  # mutate(SU = factor(SU, levels = c(0, 1), labels = c("No", "Yes")),
  #        MFN = factor(MFN, levels = c(0, 1), labels = c("No", "Yes")),
  #        DPP4 = factor(DPP4, levels = c(0, 1), labels = c("No", "Yes")),
  #        TZD = factor(TZD, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  mutate(ncurrtx = DPP4 + SGLT2 + GLP1 + TZD + SU + MFN + Acarbose + Glinide - 1) %>%
  mutate(ncurrtx = ifelse(ncurrtx > 3, 3, ncurrtx)) %>%
  mutate(ncurrtx = factor(ncurrtx, levels = c(0, 1, 2, 3), labels = c("0", "1", "2", "3"))) %>%
  
  #####   - Outcome month: hba1cmonth
  mutate(hba1cmonth_12 = difftime(posthba1c12mdate, dstartdate, units = "days") / 30) %>%
  mutate(hba1cmonth_6 = difftime(posthba1c6mdate, dstartdate, units = "days") / 30) %>%
  mutate(hba1cmonth = ifelse(is.na(hba1cmonth_6), hba1cmonth_12, hba1cmonth_6)) %>%
  mutate(hba1cmonth = as.numeric(hba1cmonth))


##### Biomarkers

#####   - hba1c: prehba1c (Nothing to do)
#####   - BMI: prebmi (Nothing to do)
#####   - eGFR: preegfr (Nothing to do)
#####   - Albumin:Creatine ratio: preacr (Nothing to do)
#####   - Serum albumin: prealbumin_blood (remove the _ from the name)

cprd <- cprd %>%
  dplyr::rename("prealbuminblood" = "prealbumin_blood",
                "prealbuminblooddate" = "prealbumin_blooddate",
                "prealbuminblooddrugdiff" = "prealbumin_blooddrugdiff")

#####   - Alanine aminotransferase: prealt (Nothing to do)
#####   - Aspartate aminotransferase: preast (Nothing to do)
#####   - Bilirubin: prebilirubin (Nothing to do)
#####   - Fasting glucose: prefastingglucose (Nothing to do)
#####   - Fasting haematocrit: prehaematocrit (Nothing to do)
#####   - Fasting haemoglobin: prehaemoglobin (Nothing to do)
#####   - High-density lipoprotein (HDL): prehdl (Nothing to do)

#####   - Mean arterial BP: premap

cprd <- cprd %>%
  mutate(premap = predbp + ((presbp - predbp) / 3))

#####   - Total cholesterol: pretotalcholesterol (Nothing to do)
#####   - Triglycerides: pretriglyceride (Nothing to do)



##### Comorbidities

cprd <- cprd %>%
  #####   - Angina: predrug_earliest_angina
  mutate(preangina = factor(predrug_angina, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Chronic Liver Disease: predrug_earliest_cld
  mutate(precld = factor(predrug_cld, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Diabetic Nephropathy: predrug_earliest_diabeticnephropathy
  mutate(prediabeticnephropathy = factor(predrug_diabeticnephropathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Heart failure: predrug_earliest_heartfailure
  mutate(preheartfailure = factor(predrug_heartfailure, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Hypertension: predrug_earliest_hypertension
  mutate(prehypertension = factor(predrug_hypertension, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Ischaemic Heart Disease: predrug_earliest_ihd
  mutate(preihd = factor(predrug_ihd, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Myocardial Infarction: predrug_earliest_myocardialinfarction
  mutate(premyocardialinfarction = factor(predrug_myocardialinfarction, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Neuropathy: predrug_earliest_neuropathy
  mutate(preneuropathy = factor(predrug_neuropathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Peripheral Arterial Disease: predrug_earliest_pad
  mutate(prepad = factor(predrug_pad, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Retinopathy: predrug_earliest_retinopathy
  mutate(preretinopathy = factor(predrug_retinopathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Cardiac Revascularisation: predrug_earliest_revasc
  mutate(prerevasc = factor(predrug_revasc, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Stroke: predrug_earliest_stroke
  mutate(prestroke = factor(predrug_stroke, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Transient Ischaemic Attack: predrug_earliest_tia
  mutate(pretia = factor(predrug_tia, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Atrial fibrillation: predrug_earliest_af
  mutate(preaf = factor(predrug_af, levels = c(0, 1), labels = c("No", "Yes")))


#################### Final dataset - all patients 
#
# Add all variables necessary for ALL analysis in the paper.
#

final.dataset <- cprd %>%
  select(
    # information regarding patient
    patid, pated, multi_drug_start, timeprevcombo, drugsubstances,gp_record_end,death_date,
    # response hba1c
    posthba1cfinal,
    # therapies of interest
    drugclass,
    # Sociodemographic features
    agetx, sex, t2dmduration, ethnicity,  ethnicity16, deprivation, smoke, prehospitalisation,
    # Diabetes treatment 
    drugline, ncurrtx, hba1cmonth, dstartdate, dstopdate, yrdrugstart, dcstopdate,
    # Biomarkers
    prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
    prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
    # Comorbidities
    preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
    preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
    # Weight analysis
    preweight, postweight12m, postweight6m, postweight12mdate, postweight6mdate, postweightfinal,
    # Discontinuation
    stopdrug_6m_3mFU,
    # eGFR analysis
    postegfr12m, postegfr6m,
    # mv complications outcome
    postdrug_first_diabeticnephropathy, postdrug_first_neuropathy, postdrug_first_retinopathy
  ) %>%
  as.data.frame()

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### Final dataset - all patients")
  print("################################################")
  print(nrow(final.dataset))
  print(table(final.dataset$drugclass))
  
}

# if full cohort was requested
if (dataset.type == "full.cohort") {
  return(final.dataset)
}

############################### HbA1c model 

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### HbA1c model")
  print("################################################")
  
}

##### Drop duplicates (i.e. started treatment on same day)

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### Drop duplicates (i.e. started treatment on same day)")
  print("################################################")
  print(table(final.dataset$multi_drug_start))
  print(table(final.dataset$multi_drug_start, final.dataset$drugclass))
}

final.dataset <- final.dataset %>%
  filter(multi_drug_start == 0)

##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)

final.dataset <- final.dataset %>%
  mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
  print("################################################")
  print(table(final.dataset$timeprevcombo_less61))
  print(table(final.dataset$timeprevcombo_less61, final.dataset$drugclass))
  
}

final.dataset <- final.dataset %>%
  filter(is.na(timeprevcombo_less61))

##### Drop if HbA1c <53

final.dataset <- final.dataset %>%
  mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### Drop if HbA1c <53")
  print("################################################")
  print(table(final.dataset$hb_extreme_53))
  print(table(final.dataset$hb_extreme_53, final.dataset$drugclass))
  
}

final.dataset <- final.dataset %>%
  filter(is.na(hb_extreme_53))

##### Drop if HbA1c >120

final.dataset <- final.dataset %>%
  mutate(hb_extreme_120 = ifelse(prehba1c >120, 1, NA_real_))

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### Drop if HbA1c <120")
  print("################################################")
  print(table(final.dataset$hb_extreme_120))
  print(table(final.dataset$hb_extreme_53, final.dataset$drugclass))
  
}

final.dataset <- final.dataset %>%
  filter(is.na(hb_extreme_120))

##### Drop if HbA1c is missing

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### Drop if HbA1c is missing")
  print("################################################")
  print(table(is.na(final.dataset$prehba1c)))
  print(table(is.na(final.dataset$prehba1c), final.dataset$drugclass))
  
}

final.dataset <- final.dataset %>%
  filter(!is.na(prehba1c))

##### Drop if post HbA1c missing

# # printing inclusion patients
# if (dataset.type == "diagnostics") {
#   
#   print("################################################")
#   print("##### Drop if post HbA1c missing")
#   print("################################################")
#   print(table(is.na(final.dataset$posthba1cfinal)))
#   print(table(is.na(final.dataset$posthba1cfinal), final.dataset$drugclass))
#   
# }
# 
# # final.dataset <- final.dataset %>%
# #   filter(!is.na(posthba1cfinal))

##### Drop if ethnicity missing

# printing inclusion patients
if (dataset.type == "diagnostics") {
  
  print("################################################")
  print("##### Drop if ethnicity missing")
  print("################################################")
  print(table(final.dataset$ethnicity=="Missing"))
  print(table(final.dataset$ethnicity=="Missing", final.dataset$drugclass))
  
}

final.dataset <- final.dataset %>%
  filter(ethnicity!="Missing") %>% 
  mutate(ethnicity=factor(ethnicity))

print(table(final.dataset$ethnicity))
print(table(final.dataset$ethnicity, final.dataset$drugclass))
print(table(final.dataset$ethnicity16, final.dataset$drugclass))

#####  IMPORT ALL SGLT2, DPP4, GLP1 + TZD STARTS FOR LATER CENSORING (SGLT2 for DPP4SU arm only)

### Combine with main dataset to get next SGLT2/GLP1/TZD/DPP4 start date
### Also get latest SGLT2 stop date before drug start for DPP4 arm in case needed for sensitivity analysis

load(paste0(data_dir,"20221212_t2d_all_drug_periods.Rda"))

later_sglt2 <- final.dataset %>%
  select(patid, dstartdate) %>%
  inner_join((t2d_all_drug_periods %>%
                filter(drugclass=="SGLT2") %>%
                select(patid, next_sglt2=dstartdate)), by="patid") %>%
  filter(next_sglt2>dstartdate) %>%
  group_by(patid, dstartdate) %>%
  summarise(next_sglt2_start=min(next_sglt2, na.rm=TRUE)) %>%
  ungroup()

later_dpp4 <- final.dataset%>%
  select(patid, dstartdate) %>%
  inner_join((t2d_all_drug_periods %>%
                filter(drugclass=="DPP4") %>%
                select(patid, next_dpp4=dstartdate)), by="patid") %>%
  filter(next_dpp4>dstartdate) %>%
  group_by(patid, dstartdate) %>%
  summarise(next_dpp4_start=min(next_dpp4, na.rm=TRUE)) %>%
  ungroup()

later_glp1 <- final.dataset%>%
  select(patid, dstartdate) %>%
  inner_join((t2d_all_drug_periods %>%
                filter(drugclass=="GLP1") %>%
                select(patid, next_glp1=dstartdate)), by="patid") %>%
  filter(next_glp1>dstartdate) %>%
  group_by(patid, dstartdate) %>%
  summarise(next_glp1_start=min(next_glp1, na.rm=TRUE)) %>%
  ungroup()

later_tzd <- final.dataset %>%
  select(patid, dstartdate) %>%
  inner_join((t2d_all_drug_periods %>%
                filter(drugclass=="TZD") %>%
                select(patid, next_tzd=dstartdate)), by="patid") %>%
  filter(next_tzd>dstartdate) %>%
  group_by(patid, dstartdate) %>%
  summarise(next_tzd_start=min(next_tzd, na.rm=TRUE)) %>%
  ungroup()

last_sglt2_stop <- final.dataset %>%
  select(patid, dstartdate) %>%
  inner_join((t2d_all_drug_periods %>%
                filter(drugclass=="SGLT2") %>%
                select(patid, last_sglt2=dstopdate)), by="patid") %>%
  filter(last_sglt2<dstartdate) %>%
  group_by(patid, dstartdate) %>%
  summarise(last_sglt2_stop=min(last_sglt2, na.rm=TRUE)) %>%
  ungroup()

final.dataset <- final.dataset %>%
  left_join(later_sglt2, by=c("patid", "dstartdate")) %>%
  left_join(later_dpp4, by=c("patid", "dstartdate")) %>%
  left_join(later_glp1, by=c("patid", "dstartdate")) %>%
  left_join(later_tzd, by=c("patid", "dstartdate")) %>%
  left_join(last_sglt2_stop, by=c("patid", "dstartdate"))

#### Microvascular complications
#Pre-existing composite, and survival analysis setup

final.dataset <- final.dataset %>%
  mutate(predrug.mv = ifelse(prediabeticnephropathy=="Yes" | preneuropathy=="Yes" | preretinopathy=="Yes", 1, 0),
         postdrug.mv=pmin(postdrug_first_diabeticnephropathy, postdrug_first_neuropathy, postdrug_first_retinopathy, na.rm=TRUE)) 

#Check
# table(final.dataset$prediabeticnephropathy)
# table(final.dataset$preneuropathy)
# table(final.dataset$preretinopathy)
# table(final.dataset$predrug.mv)

#ITT MV complications

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Start SGLT2 (DPP4 arm), or DPP4 (SGLT2 arm)
## End of GP records
## MV complications

final.dataset <- final.dataset %>%
  mutate(five_years_post_dstart=dstartdate+(365.25*5),
         mv_itt_censdate=if_else(drugclass=="SGLT2", 
                                 pmin(five_years_post_dstart, 
                                      death_date, 
                                      next_dpp4_start, 
                                      gp_record_end, 
                                      postdrug.mv, na.rm=TRUE),
                                 if_else(drugclass=="DPP4", 
                                         pmin(five_years_post_dstart, 
                                              death_date, 
                                              next_sglt2_start, 
                                              gp_record_end, 
                                              postdrug.mv, na.rm=TRUE), 
                                         as.Date(NA))),
         
         mv_itt_censvar=ifelse((!is.na(postdrug.mv) & mv_itt_censdate==postdrug.mv), 1, 0),
         
         mv_itt_censtime_yrs=as.numeric(difftime(mv_itt_censdate, dstartdate, unit="days"))/365.25)

# describe(final.dataset$mv_itt_censdate)
# describe(final.dataset$mv_itt_censvar)
# describe(final.dataset$mv_itt_censtime_yrs)


#PP MV complications

# Find censoring dates - earliest of:
## 5 years from dstartdate
## Death
## Any change in glucose-lowering therapy
## End of GP records
## MV complications

final.dataset <- final.dataset %>%
  mutate(mv_pp_censdate=if_else(drugclass=="SGLT2", 
                                pmin(five_years_post_dstart, 
                                     death_date, 
                                     gp_record_end,
                                     postdrug.mv,
                                     dcstopdate, na.rm=TRUE),
                                if_else(drugclass=="DPP4", 
                                        pmin(five_years_post_dstart, 
                                             death_date, 
                                             gp_record_end, 
                                             postdrug.mv,
                                             dcstopdate, na.rm=TRUE), as.Date(NA))),
         
         mv_pp_censvar=ifelse(!is.na(postdrug.mv) & mv_pp_censdate==postdrug.mv, 1, 0),
         
         mv_pp_censtime_yrs=as.numeric(difftime(mv_pp_censdate, dstartdate, unit="days"))/365.25)

# describe(final.dataset$mv_pp_censdate)
# describe(final.dataset$mv_pp_censvar)
# describe(final.dataset$mv_pp_censtime_yrs)

#Save
save(final.dataset,file=paste0(data_dir,"final.dataset.sglt2.dpp4.val.Rda"))

#### Load data ####

#Dataset define in prep
load(paste0(data_dir,"final.dataset.sglt2.dpp4.val.Rda"))
#GOLD SGLT2-DPP4 model
load("C:/Users/jmd237/OneDrive - University of Exeter/John/Projects/2019_SGLT2vsDPP4/results/m1_hba1cmodel.Rdata")
#GOLD SGLT2-DPP4 model for discontinuation
load("C:/Users/jmd237/OneDrive - University of Exeter/John/Projects/2019_SGLT2vsDPP4/results/m1_discontinuationmodel.Rdata")
#GOLD SGLT2-DPP4 model for weight
load("C:/Users/jmd237/OneDrive - University of Exeter/John/Projects/2019_SGLT2vsDPP4/results/m1_weightmodel.Rdata")

#### Useful functions and global settings ####

#Global settings
  
  #model.name <- "unadj"
  #OR
  model.name <- "full"
  #OR
  #model.name <- "simple"
  
  
  pdfwidth <- 14
  pdfheight <- 10
  pngwidth <- 3200
  pngheight <- 2400
  pngres <- 200
  row_names <- matrix(c("Predicted SGLT2i glycaemic benefit",paste0(intToUtf8(8805),"5 mmol/mol"), "3-5 mmol/mol","0-3 mmol/mol",
                        "Predicted DPP4i glycaemic benefit","0-3 mmol/mol", paste0(intToUtf8(8805),"3 mmol/mol")))
  tick.hb <- c(-10,-7.5,-5,-2.5,0,2.5,5,7.5, 10)
  tick.hb.resp <- c(-18,-15,-12,-9,-6,-3,0,3)
  tick.dc <- c(0,5,10,15,20,25,30,35,40)
  tick.wt <- c(-4,-3,-2,-1,0)
  B <- 1000

#Add dummy legend
  dummy <- final.dataset %>% 
    sample_n(1000) %>% 
    mutate(drugclass = ifelse(drugclass=="SGLT2","SGLT2 inhibitor","DPP-4 inhibitor")) %>%
    ggplot(aes(x = prehba1c, y = prebmi, group=drugclass)) + 
    scale_color_manual(values=c("#4118de","#f1a340")) +
    geom_point(aes(colour=drugclass), size = 1.5) + theme_bw() +
    theme(legend.text = element_text(colour="black", size=rel(1))) + 
    theme(legend.title=element_blank())  + 
    theme(legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal"
    ) +
    guides(colour = guide_legend(override.aes = list(size=4)))
  
  # Create user-defined function, which extracts legends from ggplots #https://statisticsglobe.com/add-common-legend-to-combined-ggplot2-plots-in-r/
  extract_legend <- function(my_ggp) {
    step1 <- ggplot_gtable(ggplot_build(my_ggp))
    step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
    step3 <- step1$grobs[[step2]]
    return(step3)
  }
  # Apply user-defined function to extract legend
  shared_legend <- extract_legend(dummy)  
  
#Define HTE treatment diff plot function
hist_plot <- function(data,sx,sy,y) {
  #label for hist
  annotation <- data.frame(
    x = c(sx,sy),
    y = c(y),
    label = c("Favours SGLT2i", "Favours DPP4i")
  )
  #define data
  dat <- data %>% dplyr::select(hba1c_diff) %>% mutate(above=ifelse(hba1c_diff> 0, "Favours DPP4i", "Favours SGLT2i")) 
  c_low <- quantile(dat$hba1c_diff,.001)
  c_upp <- quantile(dat$hba1c_diff,.999)
  c_lowr  <- 2*round(c_low/2)
  c_uppr <- 2*round(c_upp/2)
  c_low <- min(dat$hba1c_diff,.001)
  c_upp <- quantile(dat$hba1c_diff,.999)
  c_lowr  <- 2*round(c_low/2)
  c_uppr <- 2*round(c_upp/2)
  
  #plot
  ggplot(data=dat, aes(x=hba1c_diff,fill=above)) +
    geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-15,10,by=1)) +
    geom_vline(aes(xintercept=0), linetype="dashed")+
    labs(title="",x="HbA1c difference (mmol/mol)", y = "Number of people") +
    #scale_x_continuous(limits=c(c_low,c_upp),breaks=c(seq(c_lowr,c_uppr,by=2))) +
    scale_fill_manual(values=c("#998ec3","#f1a340"))+
    theme_classic() + theme(legend.position = "none")
  #theme(legend.position = c(0.87, 0.97)) + theme(legend.title = element_blank())
  #geom_text(data=annotation, aes(x=x, y=y, label=label),color=c("#f1a340","#998ec3"),size=4, fontface="bold" )
}   
#Function to fit a series of models and output the coefficient(s) of interest with CIs and p-value
hte.model.coefs <- function(x,nmodels) {
  mnumber = c(1:nmodels)
  models <- as.list(1:nmodels)
  nobs <- vector()
  coef <- vector()
  lower <- vector()
  upper <- vector()
  pvalue <- vector()
  data <- x
  
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(f[[i]]),data=data)
    nobs <- append(nobs,nobs(models[[i]]))
    coef <- append(coef,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower <- append(lower,confint_all[2,1])
    upper <- append(upper,confint_all[2,2])
    pvalue <- append(pvalue,summary(models[[i]])$coefficients[2,4])
  }
  
  datasetname = c(deparse(substitute(x)),deparse(substitute(x)),deparse(substitute(x)))
  x <- data.frame(datasetname,modelname,cbind(nobs,coef,lower,upper,pvalue))
  rownames(x) <- c()
  return(x)
}  
#Function to output HTE by subgroup
hte_plot <- function(data,pred,obs,obslowerci,obsupperci,ymin.ymax) {
  
  #ymin <- min(data$lci); ymax <- max(data$uci);yminr  <- 2*round(ymin/2);  ymaxr <- 2*round(ymax/2)
  ymin  <- -14;  ymax <- 14
  
  ggplot(data=data,aes_string(x=pred,y=obs)) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey60") + geom_hline(yintercept=0, linetype="dashed", color = "grey60") +
    geom_abline(intercept=0,slope=1, color="red", lwd=0.75) + ggtitle("") +
    geom_point(alpha=1) + theme_classic() +
    geom_errorbar(aes_string(ymin=obslowerci, ymax=obsupperci), colour="black", width=.1) +
    ylab("Observed HbA1c difference (mmol/mol)*") + xlab("Predicted HbA1c difference (mmol/mol)") +
    scale_x_continuous(limits=c(ymin,ymax),breaks=c(seq(ymin,ymax,by=2))) +
    scale_y_continuous(limits=c(ymin,ymax),oob=rescale_none,breaks=c(seq(ymin,ymax,by=2))) + 
    # scale_x_continuous(limits=c(ymin,ymax),breaks=c(seq(yminr,ymaxr,by=2))) +
    # scale_y_continuous(limits=c(ymin,ymax),breaks=c(seq(yminr,ymaxr,by=2))) +
    theme_base(base_size = 8)  +
    theme(plot.background = element_blank()) + theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8)) +
    theme(text = element_text(size = 10),
          # axis.ticks.x = element_blank(),
          # axis.ticks.y = element_blank(),
          axis.ticks.x = element_line(colour =  "grey50"),
          axis.ticks.y = element_line(colour =  "grey50"),
          axis.line = element_line(colour =  "grey50" ),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}
#Model and extract HTE coefficients for subgroups
hte_model <- function(data) {
  final <- data
  #Define subsets of interest to calc HTE 
  #overall
  overall <- final
  #sglt2.best
  sglt2.best <- final %>% dplyr::filter(bestdrug=="SGLT2")
  #dpp4.best
  dpp4.best <- final %>% dplyr::filter(bestdrug=="DPP4")
  
  #match trials
  #df1 <- final %>% dplyr::filter(hba1c_diff<= -10) 
  df2 <- final %>% dplyr::filter(hba1c_diff<= -5) 
  df3 <- final %>% dplyr::filter(hba1c_diff<= -3 & hba1c_diff > -5) 
  df4 <- final %>% dplyr::filter(hba1c_diff> -3 & hba1c_diff <= 0) 
  df5 <- final %>% dplyr::filter(hba1c_diff> 0 & hba1c_diff < 3) 
  df6 <- final %>% dplyr::filter(hba1c_diff>= 3) 
  #df7 <- final %>% dplyr::filter(hba1c_diff>= 5) 
  
  #Run HTE models and extract coefs for each subset
  dflist <- list(overall,sglt2.best,dpp4.best,df2,df3,df4,df5,df6)
  res.list <- lapply(dflist, function(df) {
    hte.model.coefs(df,3)
  })
  
  res.list <- bind_rows(res.list, .id="column_label")
  res.lab <- c(rep("overall",3),rep("sglt2.best",3),rep("dpp4.best",3),rep("sglt.best5",3),rep("sglt.best3",3),rep("sglt.best0-3",3),
               rep("dpp4.best0-3",3),rep("dpp4.best3",3))
  res.list <- data.frame(res.lab,res.list) %>% dplyr::select(-column_label,-datasetname)
  res.list %>% dplyr::filter(modelname==model.name) #final adjusted list
  
  #Calibration plot HTE
  
  #Unadjusted obs vs pred 
  #Define tenths
  final <- final %>% mutate(hba1c_diff.q = ntile(hba1c_diff, 10))    
  
  #define dataset with predicted values
  t1 <- final %>% 
    group_by(hba1c_diff.q) %>%
    dplyr::summarise(N = length(hba1c_diff),
                     hba1c_diff.pred = mean(hba1c_diff))
  
  #check some patients actually prescribed both drugs in each tenth
  # ddply(final, c("hba1c_diff.q","drugclass"), dplyr::summarise,
  #       N    = length(posthba1c_final),
  #       posthba1c_final.m = mean(posthba1c_final),
  #       se = sd(posthba1c_final)/sqrt((length(posthba1c_final))))
  
  #obs vs pred, by decile of predicted treatment difference
  #For Formula 1-3
  mnumber = c(1:10)
  models  <- as.list(1:10)
  
  hba1c_diff.obs.unadj <- vector()
  lower.unadj <- vector()
  upper.unadj <- vector()
  hba1c_diff.obs.sim <- vector()
  lower.sim <- vector()
  upper.sim <- vector() 
  hba1c_diff.obs.adj <- vector()
  lower.adj <- vector()
  upper.adj <- vector() 
  
  #Unadj
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(formula1),data=final,subset=hba1c_diff.q==i)
    hba1c_diff.obs.unadj <- append(hba1c_diff.obs.unadj,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower.unadj <- append(lower.unadj,confint_all[2,1])
    upper.unadj <- append(upper.unadj,confint_all[2,2])
  }
  #Simple 
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(formula2),data=final,subset=hba1c_diff.q==i)
    hba1c_diff.obs.sim <- append(hba1c_diff.obs.sim,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower.sim <- append(lower.sim,confint_all[2,1])
    upper.sim <- append(upper.sim,confint_all[2,2])
  }
  #Full
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(formula3),data=final,subset=hba1c_diff.q==i)
    hba1c_diff.obs.adj <- append(hba1c_diff.obs.adj,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower.adj <- append(lower.adj,confint_all[2,1])
    upper.adj <- append(upper.adj,confint_all[2,2])
  }
  
  #Final data.frame  
  t1 <- data.frame(t1,cbind(hba1c_diff.obs.unadj,lower.unadj,upper.unadj,
                            hba1c_diff.obs.sim,lower.sim,upper.sim,
                            hba1c_diff.obs.adj,lower.adj,upper.adj))
  
  #unadj
  # plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.unadj,lci=lower.unadj,uci=upper.unadj)
  # hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci")
  # #simple adj
  # plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.sim,lci=lower.sim,uci=upper.sim)
  # hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci") 
  #splie adj
  plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.adj,lci=lower.adj,uci=upper.adj)
  hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci")
  
  #outputs
  hist <- hist_plot(final,-2.5,2.3,1100)
  hte <- hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci")  
  res.list <- res.list %>% filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall")  %>%
    mutate(order=c(1,5,2,3,4,6,7)) %>% 
    arrange(order) 
  return(list(hist,hte,res.list))
}
#Forest plot for calibration comparison within ethnicity
fp_plot <- function(coef,cim,cip) {
  fp <-
    forestplot(row_names,
               mean = coef,
               lower= cim,
               upper = cip,
               hrzl_lines = gpar(col="#444444"),#lineheight=unit(2,'cm'),
               is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,2)),
               #title="c) TZD - Oedema",
               xticks = tick.hb,
               zero = 0,
               #boxsize=0.1,
               # graphwidth = unit(2,"inches"),
               # lineheight = unit(0.7,"inches"),
               ci.vertices=TRUE,
               col=fpColors(box=c("#66c2a5","#fc8d62","#8da0cb"), lines=c("#66c2a5","#fc8d62","#8da0cb"), zero = "gray50"),
               #lty.ci = c(1,2,3,4),
               xlab="Observed HbA1c difference (mmol/mol)* [negative favours SGLT2i]",cex=1,
               new_page = TRUE,
               fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
               boxsize = .2, # We set the box size to better visualize the type
               #line.margin = .1, # We need to add this to avoid crowding
               txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),ticks  = gpar(cex = 1)),
               #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
               # legend_args = fpLegend(
               #   pos = list("topright"),
               #   title = "Group",
               #   r = unit(.1, "snpc"),
               #   gp = gpar(col = "#CCCCCC", lwd = 1.5))
               legend = c("Uncalibrated","Recalibrated"),#,
               legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
               #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
               #xlog = TRUE
    )
  return(fp)
}
#Forest plot for calibration comparison by ethnicity
fp_plot.eth <- function(coef,cim,cip) {
  fp <-
    forestplot(row_names,
               mean = coef,
               lower= cim,
               upper = cip,
               hrzl_lines = gpar(col="#444444"),#lineheight=unit(2,'cm'),
               is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,2)),
               #title="c) TZD - Oedema",
               xticks = tick.hb,
               zero = 0,
               #boxsize=0.1,
               # graphwidth = unit(2,"inches"),
               # lineheight = unit(0.7,"inches"),
               ci.vertices=TRUE,
               col=fpColors(box=c("#66c2a5","#fc8d62","#8da0cb","#752302"), lines=c("#66c2a5","#fc8d62","#8da0cb","#752302"), zero = "gray50"),
               #lty.ci = c(1,2,3,4),
               xlab="Observed HbA1c difference (mmol/mol)* [negative favours SGLT2i]",cex=1,
               new_page = TRUE,
               fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI, fpDrawDiamondCI, fpDrawPointCI),
               boxsize = .1, # We set the box size to better visualize the type
               #line.margin = .1, # We need to add this to avoid crowding
               txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),ticks  = gpar(cex = 1)),
               #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
               # legend_args = fpLegend(
               #   pos = list("topright"),
               #   title = "Group",
               #   r = unit(.1, "snpc"),
               #   gp = gpar(col = "#CCCCCC", lwd = 1.5))
               legend = c("White (77.8%)", "Asian (13.7%)", "Black (4.2%)", "Mixed/Other (2.4%)"),#,
               legend_args = fpLegend(pos = list(x=.850, y=0.95))#, 
               #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
               #xlog = TRUE
    )
  return(fp)
}

#### Define model formula ####
formula1 <- "posthba1c_final~drugclass"
#formula2 <- "posthba1c_final~factor(drugclass)+prehba1cmmol+ncurrtx+drugline+rcs(hba1cmonth,3)+egfr_ckdepi+prealtlog"
formula2 <- "posthba1c_final~drugclass+ncurrtx+drugline"
formula3 <- "posthba1c_final~drugclass+rcs(prehba1cmmol,3)+ncurrtx+drugline+rcs(hba1cmonth,3)+rcs(egfr_ckdepi,3)+rcs(prealtlog,3)+rcs(agetx,3)+rcs(prebmi,3)"
modelname <- c("unadj","simple","full")
f <- as.list(c(formula1,formula2,formula3))

#### Define each of the HbA1c, weight, discontinuation and mv outcome cohorts ####

#Overall cohort

#Set drugline as factor
final.dataset <- final.dataset %>% mutate(drugclass=as.factor(drugclass),
                                          drugline=as.factor(drugline))

#Collapse ethnicity
final.dataset <- final.dataset %>% mutate(ethnicity.backup = ethnicity,
                                ethnicity=fct_collapse(ethnicity,mixed.other=c("Mixed","Other")))

#HbA1c cohort
final.hb <- final.dataset %>%
  select(patid, pated, posthba1cfinal, ethnicity, ethnicity16, prehba1c, drugclass, drugline, ncurrtx, preegfr, prealt, prebmi, agetx, hba1cmonth) %>%
  dplyr::rename("prealtlog"="prealt",
                "posthba1c_final"="posthba1cfinal",
                "prehba1cmmol"="prehba1c",
                "egfr_ckdepi"="preegfr") %>%
  mutate(prealtlog = log(prealtlog)) %>% 
  drop_na() #%>% 

final.hb.backup <- final.hb

#Proportion of patients by ethnicity 
prop.tab.hb <- round(prop.table(table(final.hb$ethnicity))*100,1)

#Weight cohort
final.wt <- final.dataset %>%
  select(patid, pated, preweight, postweightfinal, ethnicity, ethnicity16, prehba1c, drugclass, drugline, ncurrtx, preegfr, prealt, prebmi, agetx) %>%
  dplyr::rename("prealtlog"="prealt",
                "prehba1cmmol"="prehba1c",
                "egfr_ckdepi"="preegfr") %>%
  mutate(prealtlog = log(prealtlog),
         wtchange=postweightfinal-preweight) %>% 
  drop_na() #%>% 

#Discontinuation cohort
final.dc <- final.dataset %>%
  select(patid, pated, stopdrug_6m_3mFU, ethnicity, ethnicity16, prehba1c, drugclass, drugline, ncurrtx, preegfr, prealt, prebmi, agetx) %>%
  dplyr::rename("prealtlog"="prealt",
                "prehba1cmmol"="prehba1c",
                "egfr_ckdepi"="preegfr") %>%
  mutate(prealtlog = log(prealtlog)) %>% 
  drop_na() #%>% 

#MV complications cohort
final.mv <- final.dataset %>%
  dplyr::rename("prealtlog"="prealt",
                "prehba1cmmol"="prehba1c",
                "egfr_ckdepi"="preegfr") %>%
  mutate(prealtlog = log(prealtlog)) %>% 
  filter(predrug.mv==0)

##### Overall calibration uncalibrated ####  

#Predict outcomes
  final.hb <- final.hb %>% mutate(drug=drugclass) %>% mutate(drugclass="DPP4")
  final.hb$DPP4.pred.lm <- predict(m1,final.hb)
  final.hb <- final.hb %>% mutate(drugclass="SGLT2")
  final.hb$SGLT2.pred.lm <- predict(m1,final.hb)
  final.hb <- final.hb %>% mutate(drugclass=drug) %>%
    mutate(hba1c_diff = SGLT2.pred.lm-DPP4.pred.lm,
           bestdrug=ifelse(hba1c_diff<=0,"SGLT2","DPP4"),
           drugclass==drug)
  head(final.hb)

#brief summary of predicted treatment difference
  describe(final.hb$hba1c_diff)
  hist(final.hb$hba1c_diff,breaks=50); abline(v = 0, col="black", lwd=3, lty=2)
  table(final.hb$bestdrug)
  
  # ddply(final.hb, "bestdrug", dplyr::summarise,
  #       N    = length(hba1c_diff),
  #       hba1c_diff.sd = sd(hba1c_diff),
  #       hba1c_diff = mean(hba1c_diff))

#Run HTE models and extract coefs
  hte.output.overall <- hte_model(final.hb)

#outputs
  hist.overall <- hte.output.overall[1]
  hte.overall <- hte.output.overall[2]  
  res.list.overall <- data.frame(hte.output.overall[3]) %>% filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall") 
  hist.overall  
  hte.overall
  res.list.overall

  val <- res.list.overall %>% filter(res.lab != "sglt2.best" & res.lab != "dpp4.best")

#Plot

  x <- rep(NA,4)
  coef = rbind(x, data.frame(cbind(val[1:5,4])))
  cim = rbind(x, data.frame(cbind(val[1:5,5])))
  cip = rbind(x, data.frame(cbind(val[1:5,6])))
  
  insertrow <- function(existingDF, newrow, r) {
    existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
    existingDF[r,] <- newrow
    existingDF
  }
  
  coef <- data.matrix(insertrow(coef, NA, 5))
  cim <- data.matrix(insertrow(cim, NA, 5))
  cip <- data.matrix(insertrow(cip, NA, 5))
  
  
  row_names.n <- matrix(c("Predicted SGLT2i benefit",paste0(intToUtf8(8805),">5 mmol/mol (n=42,057)"), "3-5 mmol/mol (n=19,808)","0-3 mmol/mol (n=21,171)",
                        "Predicted DPP4i benefit","0-3 mmol/mol (n=10,244)",paste0(intToUtf8(8805),"3 mmol/mol (n=4,207)")))
  # pdf.options(reset = TRUE, onefile = FALSE)
  # pdf("treatmenteffectval(3).pdf",width=8,height=6)
  fp <-
    forestplot(row_names.n,
               mean = coef,
               lower= cim,
               upper = cip,
               hrzl_lines = gpar(col="#444444"),lineheight=unit(2,'cm'),
               is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,2)),
               #title="c) TZD - Oedema",
               xticks = tick.hb,
               zero = 0,
               #boxsize=0.1,
               # graphwidth = unit(2,"inches"),
               # lineheight = unit(0.7,"inches"),
               ci.vertices=TRUE,
               col=fpColors(box=c("#66c2a5"), lines=c("#66c2a5"), zero = "gray50"),
               #lty.ci = c(1,2,3,4),
               xlab="Average HbA1c difference (mmol/mol; negative favours SGLT2i)",cex=1,
               new_page = TRUE,
               #fn.ci_norm = c(fpDrawNormalCI),
               boxsize = .1, # We set the box size to better visualize the type
               #line.margin = .1, # We need to add this to avoid crowding
               txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),ticks  = gpar(cex = 1)),
               #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
               # legend_args = fpLegend(
               #   pos = list("topright"),
               #   title = "Group",
               #   r = unit(.1, "snpc"),
               #   gp = gpar(col = "#CCCCCC", lwd = 1.5))
               #legend = c("White (77.8%)", "Asian (13.7%)", "Black (4.2%)", "Mixed/Other (2.4%)"),#,
               #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
               #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
               #xlog = TRUE
    )
  
  
  fp.overall <- fp
  fp.overall

  final.hb <- final.hb %>% 
    mutate(DPP4.pred.lm.uncal=DPP4.pred.lm,
           SGLT2.pred.lm.uncal=SGLT2.pred.lm,
           hba1c_diff.uncal=hba1c_diff,
           bestdrug.uncal=bestdrug) %>%
    select(-DPP4.pred.lm,-SGLT2.pred.lm,-hba1c_diff,-bestdrug)
  head(final.hb)
  
#### Closed loop HbA1c model update ####

#global settings
p.value <- 0.05
ncolx <- length(m1$coefficients)-1
set.seed(8731)
sample_frac <- 1

#Predict outcome on therapy received
final.hb$pred <- predict(m1,final.hb)

#closed testing function
closedtest <- function(cohort,dataset,observed,predicted,p.value){
  
  #Original model
  
  #Residuals
  resid <- observed-predicted
  
  #variance of residuals
  sigma2 <- var(resid)
  
  #original log-likelihood
  n <- length(resid)
  logLik.original <- -n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * sum(resid^2)
  
  #Update intercept
  
  #Model with updated intercept
  m <- lm(observed-predicted~1,data=dataset)
  
  #Extract coefficient
  m1.intercept <- cbind(m$coefficients[1],confint(m)[1],confint(m)[2])
  
  #Residuals (actual - predicted)
  resid <- residuals(m)
  
  #variance of residuals
  sigma2 <- var(resid)
  
  #log-likelihood
  n <- length(resid)
  logLik.intercept <- -n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * sum(resid^2)
  
  #Update slope & intercept
  
  #Model with updated slope &  intercept
  m <- lm(observed~predicted,data=dataset)
  
  #Extract coefficient
  m2.intercept <- cbind(m$coefficients[1],confint(m)[1,1],confint(m)[1,2])
  m2.slope <- cbind(m$coefficients[2],confint(m)[2,1],confint(m)[2,2])
  
  #Residuals (actual - predicted)
  resid <- residuals(m)
  
  #variance of residuals
  sigma2 <- var(resid)
  
  #log-likelihood
  n <- length(resid)
  logLik.recal <- -n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * sum(resid^2)
  
  #test significance
  
  #1. Test recal in the large against the original model (no extra coeffs estimated)
  #2. If 2. is significant, test full recal against the recal in the large model using p df + 1 (1 extra coef estimated)
  #3. If 3. is significant, select 3. as final model, if not select 2. If neither, select 1 
  
  #ll diff
  dev_intercept <- -2*logLik.original + 2*logLik.intercept
  dev_recal <- -2*logLik.intercept + 2*logLik.recal
  
  #Diff in ll
  ncolx <- ncolx
  test1 <- (1-pchisq(dev_intercept, ncolx)) < p.value
  test2 <- (1-pchisq(dev_recal, ncolx+1)) < p.value
  
  #p.value
  p1 <- (1-pchisq(dev_intercept, ncolx))
  p2 <- (1-pchisq(dev_recal, ncolx+1))
  
  #Which model is chosen
  test_intercept <- 1 * (!test1)
  test_recal <- 2 * ((!test1)&(!test2))
  
  index_test <- (test_intercept + test_recal)
  
  res <- data.frame(cohort=c(cohort,cohort,cohort),
                    n=c(nrow(dataset),nrow(dataset),nrow(dataset)),
                    model=c("Original","Updated intercept","Recalibrated"),
                    loglikelihood=c(logLik.original,logLik.intercept,logLik.recal),
                    intercept=c(NA,
                                m1.intercept[1],
                                m1.intercept[2]),
                    intercept.w.ci=c(NA,
                                     paste0(round(m1.intercept[1],2), " (",paste0(round(m1.intercept[2],2),", ",paste0(round(m1.intercept[3],2)),")")),
                                     paste0(round(m2.intercept[1],2), " (",paste0(round(m2.intercept[2],2),", ",paste0(round(m2.intercept[3],2)),")"))),
                    
                    slope=c(NA,NA,m2.slope[1]),
                    slope.w.ci=c(NA,
                                 NA,
                                 paste0(round(m2.slope[1],2), " (",paste0(round(m2.slope[2],2),", ",paste0(round(m2.slope[3],2)),")"))),
                    p.value=c(NA,
                              round(p1,5),
                              round(p2,5)),
                    model.selected=c(ifelse(test1==FALSE & test2==FALSE,"Yes","No"),
                                     ifelse(test1==TRUE & test2==FALSE,"Yes","No"),
                                     ifelse(test2==TRUE,"Yes","No"))
  )
  return(res)
}

#Application

#DPP4i.White

#Setup
cohort <- "DPP4i.White"

dataset <- final.hb %>% 
  filter(ethnicity=="White" & drugclass=="DPP4") %>%
  sample_frac(sample_frac)

observed <- dataset$posthba1c_final
predicted <- dataset$pred

#Test
DPP4i.White <- closedtest(cohort,dataset,observed,predicted,p.value)
DPP4i.White


#DPP4i.Asian

#Setup
cohort <- "DPP4i.Asian"

dataset <- final.hb %>% 
  filter(ethnicity=="South Asian" & drugclass=="DPP4") %>%
  sample_frac(sample_frac)

observed <- dataset$posthba1c_final
predicted <- dataset$pred

#Test
DPP4i.Asian <- closedtest(cohort,dataset,observed,predicted,p.value)
DPP4i.Asian

#DPP4i.Black

#Setup
cohort <- "DPP4i.Black"

dataset <- final.hb %>% 
  filter(ethnicity=="Black" & drugclass=="DPP4") %>%
  sample_frac(sample_frac)

observed <- dataset$posthba1c_final
predicted <- dataset$pred

#Test
DPP4i.Black <- closedtest(cohort,dataset,observed,predicted,p.value)
DPP4i.Black

#DPP4i.Mixed.other

#Setup
cohort <- "DPP4i.Mixed.Other"

dataset <- final.hb %>% 
  filter((ethnicity=="mixed.other") & drugclass=="DPP4") %>%
  sample_frac(sample_frac)

observed <- dataset$posthba1c_final
predicted <- dataset$pred

#Test
DPP4i.Mixed.Other <- closedtest(cohort,dataset,observed,predicted,p.value)
DPP4i.Mixed.Other

#SGLT2i.White

#Setup
cohort <- "SGLT2i.White"

dataset <- final.hb %>% 
  filter(ethnicity=="White" & drugclass=="SGLT2") %>%
  sample_frac(sample_frac)

observed <- dataset$posthba1c_final
predicted <- dataset$pred

#Test
SGLT2i.White <- closedtest(cohort,dataset,observed,predicted,p.value)
SGLT2i.White


#SGLT2i.Asian

#Setup
cohort <- "SGLT2i.Asian"

dataset <- final.hb %>% 
  filter(ethnicity=="South Asian" & drugclass=="SGLT2") %>%
  sample_frac(sample_frac)

observed <- dataset$posthba1c_final
predicted <- dataset$pred

#Test
SGLT2.Asian <- closedtest(cohort,dataset,observed,predicted,p.value)
SGLT2.Asian

#SGLT2.Black

#Setup
cohort <- "SGLT2i.Black"

dataset <- final.hb %>% 
  filter(ethnicity=="Black" & drugclass=="SGLT2") %>%
  sample_frac(sample_frac)

observed <- dataset$posthba1c_final
predicted <- dataset$pred

#Test
SGLT2i.Black <- closedtest(cohort,dataset,observed,predicted,p.value)
SGLT2i.Black

#SGLT2i.Mixed.other

#Setup
cohort <- "SGLT2i.Mixed.Other"

dataset <- final.hb %>% 
  filter((ethnicity=="mixed.other") & drugclass=="SGLT2") %>%
  sample_frac(sample_frac)

observed <- dataset$posthba1c_final
predicted <- dataset$pred

#Test
SGLT2i.Mixed.Other <- closedtest(cohort,dataset,observed,predicted,p.value)
SGLT2i.Mixed.Other


#final.hb results
closedtest.final.hb <- rbind(DPP4i.White,SGLT2i.White,DPP4i.Asian,SGLT2.Asian,DPP4i.Black,SGLT2i.Black,DPP4i.Mixed.Other,SGLT2i.Mixed.Other)
closedtest.final.hb 

ctfm <- closedtest.final.hb %>% filter(model.selected=="Yes")

ctfm.save <- closedtest.final.hb %>% filter(model.selected=="Yes") %>% 
  select("Therapy/Ethnicity arm"=cohort,"Number of people"=n,"Model selected" = model, "Updated intercept (if required)" = intercept.w.ci)
write.csv(ctfm,file=paste0(output_dir,"closedloopmodelupdate.csv"))


#Update predictions
# final.hb <- final.hb %>% mutate(drug=drugclass) %>% mutate(drugclass="DPP4")
# final.hb$DPP4.pred.lm <- predict(m1,final.hb)
# final.hb <- final.hb %>% mutate(drugclass="SGLT2")
# final.hb$SGLT2.pred.lm <- predict(m1,final.hb)
# final.hb <- final.hb %>% mutate(drugclass=drug,
#                           hba1c_diff = SGLT2.pred.lm-DPP4.pred.lm,
#                           bestdrug=ifelse(hba1c_diff<=0,"SGLT2","DPP4"))


final.hb <- final.hb %>% mutate(DPP4.pred.lm.recal=ifelse(ethnicity == "White",DPP4.pred.lm.uncal+ctfm$intercept[1],
                                                    ifelse(ethnicity == "South Asian",DPP4.pred.lm.uncal+ctfm$intercept[3],
                                                           ifelse(ethnicity == "Black",DPP4.pred.lm.uncal+ctfm$intercept[5],
                                                                  ifelse(ethnicity == "mixed.other",DPP4.pred.lm.uncal+ctfm$intercept[7],NA
                                                                  )))),
                          SGLT2.pred.lm.recal=ifelse(ethnicity == "White",SGLT2.pred.lm.uncal+ctfm$intercept[2],SGLT2.pred.lm.uncal),
                          # ifelse(ethnicity == "South Asian",SGLT2.pred.lm+ctfm$intercept[4],
                          #        ifelse(ethnicity == "Black",SGLT2.pred.lm+ctfm$intercept[6],
                          #               ifelse(ethnicity == "Mixed"|ethnicity=="Other",SGLT2.pred.lm+ctfm$intercept[8],NA
                          #               )))),
                          hba1c_diff.recal = SGLT2.pred.lm.recal-DPP4.pred.lm.recal,
                          bestdrug.recal=ifelse(hba1c_diff.recal<=0,"SGLT2","DPP4"))

table(final.hb$bestdrug.uncal,final.hb$bestdrug.recal,final.hb$ethnicity)                                                  



##### Model performance Overall with and without calibration by ethnicity ####  

#GOLD model 
final.hb$hba1c_diff <- final.hb$hba1c_diff.uncal #GOLD model
final.hb$bestdrug <- final.hb$bestdrug.uncal #GOLD model
overall.gold <- hte_model(final.hb)
#recalibrated
final.hb$hba1c_diff <- final.hb$hba1c_diff.recal 
final.hb$bestdrug <- final.hb$bestdrug.recal 
overall.recal <- hte_model(final.hb)
# #recalibrated slope and intercept
# final.hb$hba1c_diff <- final.hb$hba1c_diff.intlp #recalibrated slope and  intercept
# overall.intlp <- hte_model(final.hb)

final.hb$hba1c_diff <- NULL
final.hb$bestdrug <- NULL


#outputs

#Treatment diff histograms
p1 <- grid2grob(print(overall.gold[1]))
p2 <- grid2grob(print(overall.recal[1]))
#p3 <- grid2grob(print(overall.intlp[1]))

hist.overall.recal <- wrap_elements(p1) / wrap_elements(p2) #/ wrap_elements(p3)
hist.overall.recal

#Calibration plots
p1 <- grid2grob(print(overall.gold[2]))
p2 <- grid2grob(print(overall.recal[2]))
#p3 <- grid2grob(print(overall.intlp[2]))

cal.overall.recal <- wrap_elements(p1) / wrap_elements(p2) #/ wrap_elements(p3)
cal.overall.recal

#HTE subgroups
res.list.overall.gold <- data.frame(overall.gold[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Uncalibrated")
res.list.overall.recal <- data.frame(overall.recal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Recalibrated")
# res.list.overall.intlp <- data.frame(overall.intlp[3]) %>% 
#   filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
#   mutate(analysis="Recalibrated slope & intercept")
res.list.overall.comparison <- rbind(res.list.overall.gold,res.list.overall.recal)#,res.list.overall.intlp)
res.list.overall.comparison

#Plot
val <- res.list.overall.comparison

x <- rep(NA,2)
coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4]))))#,val[11:15,4]))))
cim = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5]))))#,val[11:15,5]))))
cip = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6]))))#,val[11:15,6]))))

coef <- rbind(coef[1:4,],c(NA),coef[5:6,])
cim <- rbind(cim[1:4,],c(NA),cim[5:6,])
cip <- rbind(cip[1:4,],c(NA),cip[5:6,])

fp.overall.recal <- fp_plot(coef,cim,cip)
fp.overall.recal 

##### Model performance By ethnicity ####

#Asian
final <- final.hb %>% filter(ethnicity=="South Asian")

#GOLD model (same as code above)
final$hba1c_diff <- final$hba1c_diff.uncal #GOLD model
final$bestdrug <- final$bestdrug.uncal #GOLD model
overall.uncal <- hte_model(final)
#recalibrated intercept
final$hba1c_diff <- final$hba1c_diff.recal #recalibrated
final$bestdrug <- final$bestdrug.recal 
overall.recal <- hte_model(final)
# #recalibrated slope and intercept
# final$hba1c_diff <- final$hba1c_diff.recallp #recalibrated slope and  intercept
# overall.recallp <- hte_model(final)

#outputs

#Treatment diff histograms
p1 <- grid2grob(print(overall.uncal[1]))
p2 <- grid2grob(print(overall.recal[1]))
#p3 <- grid2grob(print(overall.recallp[1]))

hist.asian.both <- wrap_elements(p1) + wrap_elements(p2) #/ wrap_elements(p3)
hist.asian.uncal <- wrap_elements(p1) #/ wrap_elements(p3)
hist.asian.recal <- wrap_elements(p2) #/ wrap_elements(p3)
hist.asian.both

#Calibration plots
p1 <- grid2grob(print(overall.uncal[2]))
p2 <- grid2grob(print(overall.recal[2]))
#p3 <- grid2grob(print(overall.recallp[2]))

cal.asian.both <- wrap_elements(p1) + wrap_elements(p2) #/ wrap_elements(p3)
cal.asian.uncal <- wrap_elements(p1) #/ wrap_elements(p3)
cal.asian.recal <- wrap_elements(p2) #/ wrap_elements(p3)
cal.asian.both

#HTE subgroups
res.list.overall.uncal <- data.frame(overall.uncal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Uncalibrated")
res.list.overall.recal <- data.frame(overall.recal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Recalibrated intercept")
# res.list.overall.recallp <- data.frame(overall.recallp[3]) %>% 
#   filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
#   mutate(analysis="Recalibrated slope & intercept")
res.list.overall.comparison <- rbind(res.list.overall.uncal,res.list.overall.recal)#,res.list.overall.recallp)
res.list.overall.comparison

val <- res.list.overall.comparison

x <- rep(NA,2)
coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4]))))#,val[11:15,4]))))
cim = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5]))))#,val[11:15,5]))))
cip = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6]))))#,val[11:15,6]))))

coef <- rbind(coef[1:4,],c(NA),coef[5:6,])
cim <- rbind(cim[1:4,],c(NA),cim[5:6,])
cip <- rbind(cip[1:4,],c(NA),cip[5:6,])

fp.asian.recal <- fp_plot(coef,cim,cip)
fp.asian.recal 
val.asian.recal <- val

#White
final <- final.hb %>% filter(ethnicity=="White")

#uncal model (same as code above)
final$hba1c_diff <- final$hba1c_diff.uncal #uncal model
final$bestdrug <- final$bestdrug.uncal #GOLD model
overall.uncal <- hte_model(final)
#recalibrated intercept
final$hba1c_diff <- final$hba1c_diff.recal #recalibrated
final$bestdrug <- final$bestdrug.recal 
overall.recal <- hte_model(final)
#recalibrated slope and intercept
# final$hba1c_diff <- final$hba1c_diff.recallp #recalibrated slope and  intercept
# overall.recallp <- hte_model(final)

#outputs

#Treatment diff histograms
p1 <- grid2grob(print(overall.uncal[1]))
p2 <- grid2grob(print(overall.recal[1]))
#p3 <- grid2grob(print(overall.recallp[1]))

hist.white.both <- wrap_elements(p1) + wrap_elements(p2) #/ wrap_elements(p3)
hist.white.uncal <- wrap_elements(p1) #/ wrap_elements(p3)
hist.white.recal <- wrap_elements(p2) #/ wrap_elements(p3)
hist.white.both

#Calibration plots
p1 <- grid2grob(print(overall.uncal[2]))
p2 <- grid2grob(print(overall.recal[2]))
#p3 <- grid2grob(print(overall.recallp[2]))

cal.white.both <- wrap_elements(p1) + wrap_elements(p2) #/ wrap_elements(p3)
cal.white.uncal <- wrap_elements(p1) #/ wrap_elements(p3)
cal.white.recal <- wrap_elements(p2) #/ wrap_elements(p3)
cal.white.both

#HTE subgroups
res.list.overall.uncal <- data.frame(overall.uncal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Uncalibrated")
res.list.overall.recal <- data.frame(overall.recal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Recalibrated intercept")
# res.list.overall.recallp <- data.frame(overall.recallp[3]) %>% 
#   filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
#   mutate(analysis="Recalibrated slope & intercept")
res.list.overall.comparison <- rbind(res.list.overall.uncal,res.list.overall.recal)#,res.list.overall.recallp)
res.list.overall.comparison

val <- res.list.overall.comparison

x <- rep(NA,2)
coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4]))))#,val[11:15,4]))))
cim = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5]))))#,val[11:15,5]))))
cip = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6]))))#,val[11:15,6]))))

coef <- rbind(coef[1:4,],c(NA),coef[5:6,])
cim <- rbind(cim[1:4,],c(NA),cim[5:6,])
cip <- rbind(cip[1:4,],c(NA),cip[5:6,])

fp.white.recal <- fp_plot(coef,cim,cip)
fp.white.recal 
val.white.recal <- val

#Black
final <- final.hb %>% filter(ethnicity=="Black")

#uncal model (same as code above)
final$hba1c_diff <- final$hba1c_diff.uncal #uncal model
final$bestdrug <- final$bestdrug.uncal #GOLD model
overall.uncal <- hte_model(final)
#recalibrated intercept
final$hba1c_diff <- final$hba1c_diff.recal #recalibrated
final$bestdrug <- final$bestdrug.recal 
overall.recal <- hte_model(final)
# #recalibrated slope and intercept
# final$hba1c_diff <- final$hba1c_diff.recallp #recalibrated slope and  intercept
# overall.recallp <- hte_model(final)

#outputs

#Treatment diff histograms
p1 <- grid2grob(print(overall.uncal[1]))
p2 <- grid2grob(print(overall.recal[1]))
#p3 <- grid2grob(print(overall.recallp[1]))

hist.black.both <- wrap_elements(p1) + wrap_elements(p2) #/ wrap_elements(p3)
hist.black.uncal <- wrap_elements(p1) #/ wrap_elements(p3)
hist.black.recal <- wrap_elements(p2) #/ wrap_elements(p3)
hist.black.both

#Calibration plots
p1 <- grid2grob(print(overall.uncal[2]))
p2 <- grid2grob(print(overall.recal[2]))
#p3 <- grid2grob(print(overall.recallp[2]))

cal.black.both <- wrap_elements(p1) + wrap_elements(p2) #/ wrap_elements(p3)
cal.black.uncal <- wrap_elements(p1) #/ wrap_elements(p3)
cal.black.recal <- wrap_elements(p2) #/ wrap_elements(p3)
cal.black.both

#HTE subgroups
res.list.overall.uncal <- data.frame(overall.uncal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Uncalibrated")
res.list.overall.recal <- data.frame(overall.recal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Recalibrated intercept")
# res.list.overall.recallp <- data.frame(overall.recallp[3]) %>% 
#   filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
#   mutate(analysis="Recalibrated slope & intercept")
res.list.overall.comparison <- rbind(res.list.overall.uncal,res.list.overall.recal)#,res.list.overall.recallp)
res.list.overall.comparison

val <- res.list.overall.comparison

x <- rep(NA,2)
coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4]))))#,val[11:15,4]))))
cim = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5]))))#,val[11:15,5]))))
cip = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6]))))#,val[11:15,6]))))

coef <- rbind(coef[1:4,],c(NA),coef[5:6,])
cim <- rbind(cim[1:4,],c(NA),cim[5:6,])
cip <- rbind(cip[1:4,],c(NA),cip[5:6,])

fp.black.recal <- fp_plot(coef,cim,cip)
fp.black.recal 
val.black.recal <- val

#Mixed.Other
final <- final.hb %>% filter(ethnicity=="mixed.other")

#uncal model (same as code above)
final$hba1c_diff <- final$hba1c_diff.uncal #uncal model
final$bestdrug <- final$bestdrug.uncal #GOLD model
overall.uncal <- hte_model(final)
#recalibrated intercept
final$hba1c_diff <- final$hba1c_diff.recal #recalibrated
final$bestdrug <- final$bestdrug.recal 
overall.recal <- hte_model(final)
# #recalibrated slope and intercept
# final$hba1c_diff <- final$hba1c_diff.recallp #recalibrated slope and  intercept
# overall.recallp <- hte_model(final)

#outputs

#Treatment diff histograms
p1 <- grid2grob(print(overall.uncal[1]))
p2 <- grid2grob(print(overall.recal[1]))
#p3 <- grid2grob(print(overall.recallp[1]))

hist.mixed.other.both <- wrap_elements(p1) + wrap_elements(p2) #/ wrap_elements(p3)
hist.mixed.other.uncal <- wrap_elements(p1) #/ wrap_elements(p3)
hist.mixed.other.recal <- wrap_elements(p2) #/ wrap_elements(p3)
hist.mixed.other.both

#Calibration plots
p1 <- grid2grob(print(overall.uncal[2]))
p2 <- grid2grob(print(overall.recal[2]))
#p3 <- grid2grob(print(overall.recallp[2]))

cal.mixed.other.both <- wrap_elements(p1) + wrap_elements(p2) #/ wrap_elements(p3)
cal.mixed.other.uncal <- wrap_elements(p1) #/ wrap_elements(p3)
cal.mixed.other.recal <- wrap_elements(p2) #/ wrap_elements(p3)
cal.mixed.other.both

#HTE subgroups
res.list.overall.uncal <- data.frame(overall.uncal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Uncalibrated")
res.list.overall.recal <- data.frame(overall.recal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Recalibrated intercept")
# res.list.overall.recallp <- data.frame(overall.recallp[3]) %>% 
#   filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
#   mutate(analysis="Recalibrated slope & intercept")
res.list.overall.comparison <- rbind(res.list.overall.uncal,res.list.overall.recal)#,res.list.overall.recallp)
res.list.overall.comparison

val <- res.list.overall.comparison

x <- rep(NA,2)
coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4]))))#,val[11:15,4]))))
cim = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5]))))#,val[11:15,5]))))
cip = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6]))))#,val[11:15,6]))))

coef <- rbind(coef[1:4,],c(NA),coef[5:6,])
cim <- rbind(cim[1:4,],c(NA),cim[5:6,])
cip <- rbind(cip[1:4,],c(NA),cip[5:6,])

fp.mixed.other.recal <- fp_plot(coef,cim,cip)
fp.mixed.other.recal 
val.mixed.other.recal <- val

##### Overall plot by ethnicity for each type of calibration #####

#Uncalibrated
val.uncal <- rbind(val.white.recal[1:5,],val.asian.recal[1:5,],val.black.recal[1:5,],val.mixed.other.recal[1:5,])

val <- val.uncal
x <- rep(NA,4)
coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4],val[11:15,4],val[16:20,4]))))
cim = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5],val[11:15,5],val[16:20,5]))))
cip = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6],val[11:15,6],val[16:20,6]))))

coef <- rbind(coef[1:4,],c(NA,NA,NA,NA),coef[5:6,])
cim <- rbind(cim[1:4,],c(NA,NA,NA,NA),cim[5:6,])
cip <- rbind(cip[1:4,],c(NA,NA,NA,NA),cip[5:6,])

fp.eth.uncal <-fp_plot.eth(coef,cim,cip)

#Intercept
val.int <- rbind(val.white.recal[6:10,],val.asian.recal[6:10,],val.black.recal[6:10,],val.mixed.other.recal[6:10,])

val <- val.int
x <- rep(NA,4)
coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4],val[11:15,4],val[16:20,4]))))
cim = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5],val[11:15,5],val[16:20,5]))))
cip = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6],val[11:15,6],val[16:20,6]))))

coef <- rbind(coef[1:4,],c(NA,NA,NA,NA),coef[5:6,])
cim <- rbind(cim[1:4,],c(NA,NA,NA,NA),cim[5:6,])
cip <- rbind(cip[1:4,],c(NA,NA,NA,NA),cip[5:6,])

fp.eth.int <-fp_plot.eth(coef,cim,cip)

# #Intercept + slope
# val.intlp <- rbind(val.white.recal[11:15,],val.asian.recal[11:15,],val.black.recal[11:15,],val.mixed.other.recal[11:15,])
# 
# val <- val.intlp
# x <- rep(NA,4)
# coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4],val[11:15,4],val[16:20,4]))))
# cim = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5],val[11:15,5],val[16:20,5]))))
# cip = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6],val[11:15,6],val[16:20,6]))))
# 
# coef <- rbind(coef[1:4,],c(NA,NA,NA,NA),coef[5:6,])
# cim <- rbind(cim[1:4,],c(NA,NA,NA,NA),cim[5:6,])
# cip <- rbind(cip[1:4,],c(NA,NA,NA,NA),cip[5:6,])
# 
# fp.eth.intlp <-fp_plot.eth(coef,cim,cip)

#Plot together
p1 <- grid2grob(print(fp.eth.uncal))
p2 <- grid2grob(print(fp.eth.int))
#p3 <- grid2grob(print(fp.eth.intlp))

fp.eth <- wrap_elements(p1) + wrap_elements(p2) + #/ wrap_elements(p3) +
  plot_annotation("Left=uncalibrated, Right=recalibrated intercept")#, Bottom=recalibrated slope & intercept")
#fp.eth 

grDevices::cairo_pdf(paste0(output_dir,"fp.eth_recal_both.pdf"),width=20,height=8)
fp.eth 
dev.off()

png(paste0(output_dir,"fp.eth_recal_both.png"),width=pngwidth,height=1250,res=pngres,restoreConsole=TRUE)
fp.eth 
dev.off()

grDevices::cairo_pdf(paste0(output_dir,"fp.eth_uncal.pdf"),width=10,height=8)
wrap_elements(p1)
dev.off()

png(paste0(output_dir,"fp.eth_uncal.png"),width=1800,height=1250,res=pngres,restoreConsole=TRUE)
wrap_elements(p1)
dev.off()

grDevices::cairo_pdf(paste0(output_dir,"fp.eth_recal.pdf"),width=10,height=8)
wrap_elements(p2)
dev.off()

png(paste0(output_dir,"fp.eth_recal.png"),width=1800,height=1250,res=pngres,restoreConsole=TRUE)
wrap_elements(p2)
dev.off()

#### save other outputs #####

#Hist

#recalibrated
thm <- theme(plot.title = element_text(face = 1, size = 12))
hist.white.recal.p <- wrap_elements(hist.white.recal+plot_annotation(title = paste0("White (",prop.tab.hb[1],"%)"), theme = thm))
hist.asian.recal.p <- wrap_elements(hist.asian.recal+plot_annotation(title = paste0("South Asian (",prop.tab.hb[2],"%)"), theme = thm))
hist.black.recal.p <- wrap_elements(hist.black.recal+plot_annotation(title = paste0("Black (",prop.tab.hb[3],"%)"), theme = thm))
hist.mixed.other.recal.p <- wrap_elements(hist.mixed.other.recal+plot_annotation(title = paste0("Mixed or Other (",prop.tab.hb[4],"%)"), theme = thm))

hist.eth <- 
  hist.white.recal.p +
  hist.asian.recal.p + 
  hist.black.recal.p + 
  hist.mixed.other.recal.p +
  plot_layout(ncol = 2) #+
  # plot_annotation(title="White,                                                      Asian,                                                      Black,                                                      Mixed/Other",
  #                 subtitle="Top=uncalibrated, Middle=recalibrated intercept, Bottom=recalibrated slope & intercept")
hist.eth

# hist.eth <- wrap_elements(hist.white.recal) + 
#   wrap_elements(hist.asian.recal) + 
#   wrap_elements(hist.black.recal) + 
#   wrap_elements(hist.mixed.other.recal) +
#   plot_layout(ncol = 4)  +
#   plot_annotation(title="White,                                                      Asian,                                                      Black,                                                      Mixed/Other",
#                   subtitle="Top=uncalibrated, Middle=recalibrated intercept, Bottom=recalibrated slope & intercept")

grDevices::cairo_pdf(paste0(output_dir,"hist.eth_recal.pdf"),width=8,height=8)
hist.eth
dev.off()

png(paste0(output_dir,"hist.eth_recal.png"),width=1800,height=1800,res=pngres,restoreConsole=TRUE)
hist.eth
dev.off()

round(prop.table(table(final.hb$ethnicity,final.hb$bestdrug.recal),1)*100,1)

#uncalibrated
hist.white.uncal.p <- wrap_elements(hist.white.uncal+plot_annotation(title = paste0("White (",prop.tab.hb[1],"%)"), theme = thm))
hist.asian.uncal.p <- wrap_elements(hist.asian.uncal+plot_annotation(title = paste0("South Asian (",prop.tab.hb[2],"%)"), theme = thm))
hist.black.uncal.p <- wrap_elements(hist.black.uncal+plot_annotation(title = paste0("Black (",prop.tab.hb[3],"%)"), theme = thm))
hist.mixed.other.uncal.p <- wrap_elements(hist.mixed.other.uncal+plot_annotation(title = paste0("Mixed or Other (",prop.tab.hb[4],"%)"), theme = thm))

hist.eth <- 
  hist.white.uncal.p +
  hist.asian.uncal.p + 
  hist.black.uncal.p + 
  hist.mixed.other.uncal.p +
  plot_layout(ncol = 2) 

grDevices::cairo_pdf(paste0(output_dir,"hist.eth_uncal.pdf"),width=8,height=8)
hist.eth
dev.off()

png(paste0(output_dir,"hist.eth_uncal.png"),width=1800,height=1800,res=pngres,restoreConsole=TRUE)
hist.eth
dev.off()

#Calibration plot

#Recalibrated
cal.white.recal.p <- wrap_elements(cal.white.recal+plot_annotation(title = paste0("White (",prop.tab.hb[1],"%)"), theme = thm))
cal.asian.recal.p <- wrap_elements(cal.asian.recal+plot_annotation(title = paste0("South Asian (",prop.tab.hb[2],"%)"), theme = thm))
cal.black.recal.p <- wrap_elements(cal.black.recal+plot_annotation(title = paste0("Black (",prop.tab.hb[3],"%)"), theme = thm))
cal.mixed.other.recal.p <- wrap_elements(cal.mixed.other.recal+plot_annotation(title = paste0("Mixed or Other (",prop.tab.hb[4],"%)"), theme = thm))

cal.eth <- 
  cal.white.recal.p +
  cal.asian.recal.p + 
  cal.black.recal.p + 
  cal.mixed.other.recal.p +
  plot_layout(ncol = 2) 

grDevices::cairo_pdf(paste0(output_dir,"cal.eth_recal.pdf"),width=8,height=8)
cal.eth
dev.off()

png(paste0(output_dir,"cal.eth_recal.png"),width=1800,height=1800,res=pngres,restoreConsole=TRUE)
cal.eth
dev.off()

#Uncalibrated
cal.white.uncal.p <- wrap_elements(cal.white.uncal+plot_annotation(title = paste0("White (",prop.tab.hb[1],"%)"), theme = thm))
cal.asian.uncal.p <- wrap_elements(cal.asian.uncal+plot_annotation(title = paste0("South Asian (",prop.tab.hb[2],"%)"), theme = thm))
cal.black.uncal.p <- wrap_elements(cal.black.uncal+plot_annotation(title = paste0("Black (",prop.tab.hb[3],"%)"), theme = thm))
cal.mixed.other.uncal.p <- wrap_elements(cal.mixed.other.uncal+plot_annotation(title = paste0("Mixed or Other (",prop.tab.hb[4],"%)"), theme = thm))

cal.eth <- 
  cal.white.uncal.p +
  cal.asian.uncal.p + 
  cal.black.uncal.p + 
  cal.mixed.other.uncal.p +
  plot_layout(ncol = 2) 

grDevices::cairo_pdf(paste0(output_dir,"cal.eth_uncal.pdf"),width=8,height=8)
cal.eth
dev.off()

png(paste0(output_dir,"cal.eth_uncal.png"),width=1800,height=1800,res=pngres,restoreConsole=TRUE)
cal.eth
dev.off()

#### HbA1c difference interaction approach ####

# 1: Subset by eth
# 2: Define deciles for that eth
# 3: Estimate predicted as the mean predicted HbA1c within each group
# 4: Estimate observed with hba1c_q*drug interaction adjusted for key variables 

white <- final.hb %>% filter(ethnicity=="White")
asian <- final.hb %>% filter(ethnicity=="South Asian")
black <- final.hb %>% filter(ethnicity=="Black")
mixed.other <- final.hb %>% filter(ethnicity=="mixed.other")

L <- list(white,asian,black,mixed.other)

names(L) <- c("White",
              "Asian",
              "Black",
              "Mixed or Other")

hb.res <- list()
hb.plot <- list()

formula3.hba1cq <- "posthba1c_final~drugclass*hba1c_diff.q+rcs(prehba1cmmol,3)+ncurrtx+drugline+rcs(hba1cmonth,3)+rcs(egfr_ckdepi,3)+rcs(prealtlog,3)+rcs(agetx,3)+rcs(prebmi,3)"

for(i in 1:4) 
{ 
  data  <- L[[i]] %>% mutate(hba1c_diff.q = factor(ntile(hba1c_diff.recal, 10))) 
  ddist <- datadist(data); options(datadist='ddist') 
  
  #define dataset with predicted values
  t1 <- final %>% 
    group_by(hba1c_diff.q) %>%
    dplyr::summarise(N = length(hba1c_diff),
                     hba1c_diff.pred = mean(hba1c_diff.recal)) %>%
    mutate(ethnicity=names(L)[[i]])
  
  #Model for adjusted observed values
  m.hb.cal <- ols(as.formula(formula3.hba1cq),data=data,x=T,y=T)
  
  #Marginal effects https://vincentarelbundock.github.io/marginaleffects/articles/comparisons.html?q=interactions#interactions-and-cross-contrasts
  comp <- avg_comparisons(m.hb.cal, variables = "drugclass", by = "hba1c_diff.q")  %>%
    arrange(as.numeric(hba1c_diff.q)) %>%
    mutate(ethnicity=names(L)[[i]]) %>%
    select(ethnicity,hba1c_diff.q,hba1c_diff=estimate,lci=conf.low,uci=conf.high) %>% 
    as.data.frame()
  
  #Merge
  hb.cal.final <- left_join(comp,t1,by=c("ethnicity","hba1c_diff.q"))
  
  #res table
  hb.res[[i]] <- hb.cal.final
  
  #Plot
  hb.plot[[i]] <-
    hte_plot(hb.cal.final,"hba1c_diff.pred","hba1c_diff","lci","uci")  
  
}   

#Plot together
p.white <- wrap_elements(grid2grob(print(hb.plot[[1]])))
p.asian <- wrap_elements(grid2grob(print(hb.plot[[2]])))
p.black <- wrap_elements(grid2grob(print(hb.plot[[3]])))
p.mixed.other <- wrap_elements(grid2grob(print(hb.plot[[4]])))

hb.eth.cal <- (
  wrap_elements(p.white+plot_annotation(title = paste0("White (",prop.tab.hb[1],"%)"), theme = thm)) + 
    wrap_elements(p.asian+plot_annotation(title = paste0("South Asian (",prop.tab.hb[2],"%)"), theme = thm))) /
  (
    wrap_elements(p.black+plot_annotation(title = paste0("Black (",prop.tab.hb[3],"%)"), theme = thm)) + 
      wrap_elements(p.mixed.other+plot_annotation(title = paste0("Mixed or Other (",prop.tab.hb[4],"%)"), theme = thm)))

hb.eth.cal

grDevices::cairo_pdf(paste0(output_dir,"cal.eth_recal_interaction.pdf"),width=8,height=8)
hb.eth.cal
dev.off()

png(paste0(output_dir,"cal.eth_recal_interaction.png"),width=1800,height=1800,res=pngres,restoreConsole=TRUE)
hb.eth.cal
dev.off()

hb.res.cal.data <- hb.res %>% bind_rows(hb.res)
write.csv(hb.res.cal.data,file=paste0(output_dir,"cal.eth_recal_interaction.csv"))  

#Global deciles - not used as won't be deciles for ethnic subgroups
# #Define tenths
# final.hb <- final.hb %>% mutate(hba1c_diff.q = factor(ntile(hba1c_diff.recal, 10)))    
# 
# #define dataset with predicted values
# t1 <- ddply(final.hb, c("ethnicity","hba1c_diff.q"), dplyr::summarise,
#             N    = length(hba1c_diff.recal),
#             hba1c_diff.pred = mean(hba1c_diff.recal))
# 
# #Define observed values #https://grantmcdermott.com/interaction-effects/
# m.hb.cal <- ols(posthba1c_final~drug*hba1c_diff.q*ethnicity+ prehba1cmmol + ncurrtx,data=final.hb,x=T,y=T)
# 
# #Marginal effects https://vincentarelbundock.github.io/marginaleffects/articles/comparisons.html?q=interactions#interactions-and-cross-contrasts
# comp <- avg_comparisons(m.hb.cal, variables = "drug", by = c("hba1c_diff.q","ethnicity"))  %>%
#   arrange(ethnicity,as.numeric(hba1c_diff.q)) %>%
#   select(ethnicity,hba1c_diff.q,hba1c_diff=estimate,lci=conf.low,uci=conf.high) %>% 
#   as.data.frame()
# 
# #Merge
# hb.cal.final <- left_join(comp,t1,by=c("ethnicity","hba1c_diff.q"))
# 
# #Plot
# plotdata <- hb.cal.final %>% filter(ethnicity=="White")
# cal.white.recal <- hte_plot(plotdata,"hba1c_diff.pred","hba1c_diff","lci","uci")
# 
# plotdata <- hb.cal.final %>% filter(ethnicity=="South Asian")
# cal.asian.recal <- hte_plot(plotdata,"hba1c_diff.pred","hba1c_diff","lci","uci")
# 
# plotdata <- hb.cal.final %>% filter(ethnicity=="Black")
# cal.black.recal <- hte_plot(plotdata,"hba1c_diff.pred","hba1c_diff","lci","uci")
# 
# plotdata <- hb.cal.final %>% filter(ethnicity=="mixed.other")
# cal.mixed.other.recal <- hte_plot(plotdata,"hba1c_diff.pred","hba1c_diff","lci","uci")
# 
# cal.white.recal.p <- wrap_elements(cal.white.recal+plot_annotation(title = "White", theme = thm))
# cal.asian.recal.p <- wrap_elements(cal.asian.recal+plot_annotation(title = "South Asian", theme = thm))
# cal.black.recal.p <- wrap_elements(cal.black.recal+plot_annotation(title = "Black", theme = thm))
# cal.mixed.other.recal.p <- wrap_elements(cal.mixed.other.recal+plot_annotation(title = "Mixed or Other", theme = thm))
# 
# cal.eth <- 
#   cal.white.recal.p +
#   cal.asian.recal.p + 
#   cal.black.recal.p + 
#   cal.mixed.other.recal.p +
#   plot_layout(ncol = 2) 
# 
# #grDevices::cairo_pdf(paste0(output_dir,"cal.eth_recal.pdf"),width=8,height=8)
# cal.eth
# #dev.off()

#### HbA1c response ####
  
  #Variable prep
  final.hb<- final.hb %>% mutate(hba1c.breaks = cut(hba1c_diff.recal, breaks=c(min(hba1c_diff.recal)-0.01,-5,-3,0,3,max(hba1c_diff.recal)+0.01)),
                                 hba1c.change=posthba1c_final-prehba1cmmol)
  
  #Unadjusted
  hb.res.unadjusted <- final.hb %>% 
    mutate(hba1c.change=posthba1c_final-prehba1cmmol) %>%
    dplyr::group_by(ethnicity,hba1c.breaks,drugclass) %>% 
    dplyr::summarise(n=length(hba1c.change),
                     hba1c.resp = mean(hba1c.change),
                     lci = mean(hba1c.change) - (1.96*(sd(hba1c.change)/sqrt(length(hba1c.change)))),
                     uci = mean(hba1c.change) + (1.96*(sd(hba1c.change)/sqrt(length(hba1c.change))))) %>%
    select(ethnicity,hba1c.breaks,drugclass,n,hba1c.resp,lci,uci)
  
  #Long to wide for plotting
  hb.res.unadjusted <- hb.res.unadjusted %>% pivot_wider(
    names_from = drugclass,
    values_from = c(n,hba1c.resp,lci,uci)
  )

  #hb.res.unadjusted$diff <- hb.res.unadjusted$hba1c.change_SGLT2- hb.res.unadjusted$hba1c.change_DPP4

  #Plot observed and 95% CI by ethnicity and HbA1c defined subgroup
  white <- hb.res.unadjusted %>% filter(ethnicity=="White") 
  asian <- hb.res.unadjusted %>% filter(ethnicity=="South Asian")
  black <- hb.res.unadjusted %>% filter(ethnicity=="Black") 
  mixed <- hb.res.unadjusted %>% filter(ethnicity=="mixed.other") 
  
  L <- list(white,asian,black,mixed)
  
  names(L) <- c("White",
                "Asian",
                "Black",
                "Mixed or Other")
  
  hb.plot <- list()
  
  
  for(i in 1:4) 
  { 
    #Subgroups by predicted treatment difference
    plotdata <- L[[i]] %>% ungroup() %>% as.data.frame()  %>%
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 1) %>% 
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 5) 
    
    #plot
    coef = data.matrix(cbind(plotdata[,6],plotdata[,5]))
    cim = data.matrix(cbind(plotdata[,8],plotdata[,7]))
    cip = data.matrix(cbind(plotdata[,10],plotdata[,9]))
    
    hb.plot[[i]] <-
      forestplot(row_names,
                 mean = coef,
                 lower= cim,
                 upper = cip,
                 hrzl_lines = gpar(col="#444444"),lineheight=unit(2,'cm'),
                 is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
                 xticks = tick.hb.resp,
                 zero = 0,
                 #boxsize=0.1,
                 # graphwidth = unit(2,"inches"),
                 # lineheight = unit(0.7,"inches"),
                 ci.vertices=TRUE,
                 col=fpColors(box=c("#f1a340","#4118de"), lines=c("#f1a340","#4118de"), zero = "gray50"),
                 lty.ci = c(1,2),           ,
                 xlab="Average glycaemic response (mmol/mol)",cex=1,
                 title = names(L[i]),
                 new_page = TRUE,
                 fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
                 boxsize = .2, # We set the box size to better visualize the type
                 #line.margin = .1, # We need to add this to avoid crowding
                 txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),summary=gpar(cex = 1), ticks  = gpar(cex = 1), title=gpar(cex=1.3))
                 #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
                 #legend = c("SGLT2i","DPP4i"),
                 #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
                 #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
                 #xlog = TRUE
      )
    
  }  
  
  
  #Plot together
  p.hb.white <- grid2grob(print(hb.plot[[1]]))
  p.hb.asian <- grid2grob(print(hb.plot[[2]]))
  p.hb.black <- grid2grob(print(hb.plot[[3]]))
  p.hb.mixed.other <- grid2grob(print(hb.plot[[4]]))
  
  hb.eth.obs <- (wrap_elements(p.hb.white) + wrap_elements(p.hb.asian)) /
    (wrap_elements(p.hb.black) + wrap_elements(p.hb.mixed.other))/
    shared_legend +
    plot_layout(height=c(10,10,1))
  
  hb.eth.obs
  
  grDevices::cairo_pdf(paste0(output_dir,"hb_eth_recal_observed_unadjusted.pdf"),width=pdfwidth,height=pdfheight)
  hb.eth.obs
  dev.off()
  
  png(paste0(output_dir,"hb_eth_recal_observed_unadjusted.png"),width=pngwidth,height=pngheight,res=pngres,restoreConsole=TRUE)
  hb.eth.obs
  dev.off()
  
  write.csv(hb.res.unadjusted,file=paste0(output_dir,"hb_eth_recal_observed_unadjusted.csv"))
  
  #rename table outputs to compare obs vs pred
  hb.res.obs <- hb.res.unadjusted

#### Adjusted

  #Mean prediction
  m.hb.obs <- ols(hba1c.change~drug*hba1c.breaks*ethnicity + prehba1cmmol + drugline + ncurrtx,data=final.hb,x=T,y=T)
  
  #Predict weight on each drug de novo
  final.hb$drug <- "DPP4"
  final.hb$DPP4.hb <- predict(m.hb.obs, newdata=final.hb, se.fit=F)
  final.hb$drug <- "SGLT2"
  final.hb$SGLT2.hb <- predict(m.hb.obs, newdata=final.hb, se.fit=F)
  final.hb$sglt2.dpp4.hb.diff <- final.hb$SGLT2.hb-final.hb$DPP4.hb
  describe(final.hb$sglt2.dpp4.hb.diff)
  hist(final.hb$sglt2.dpp4.hb.diff)
  final.hb$drug <- final.hb$drugclass
  
  hb.res.overall <- final.hb %>% 
    group_by(ethnicity,hba1c.breaks) %>% 
    dplyr::summarise(
      DPP4.hb = mean(DPP4.hb),
      SGLT2.hb = mean(SGLT2.hb),
      sglt2.dpp4.hb.diff = mean(sglt2.dpp4.hb.diff)
    )
  
  #Bootstrap predicted means
  n=nrow(final.hb)
  
  hb.res <- list()
  
  for(b in 1:B){
    i = sample(x = 1:n, size = n, replace = TRUE) ## sample indices
    temp = final.hb[i,] ## temp data set
    temp_model =   ols(hba1c.change~drug*hba1c.breaks*ethnicity + prehba1cmmol + drugline + ncurrtx,data=final.hb,x=T,y=T)
    temp$drug <- "DPP4"
    temp$DPP4.hb <-  predict(temp_model, newdata=temp, se.fit=F)
    temp$drug <- "SGLT2"
    temp$SGLT2.hb <- predict(temp_model, newdata=temp, se.fit=F)
    temp$sglt2.dpp4.hb.diff <- temp$SGLT2.hb-temp$DPP4.hb
    
    res <- temp %>% 
      group_by(ethnicity,hba1c.breaks) %>% 
      dplyr::summarise(
        DPP4.hb = mean(DPP4.hb),
        SGLT2.hb = mean(SGLT2.hb),
        sglt2.dpp4.hb.diff = mean(sglt2.dpp4.hb.diff)
      )
    
    #res table
    hb.res[[b]] <- res
  }
  
  #Derive the bootstrapped CIs for the mean predictions for each subgrup
  hb.res.ci <- bind_rows(hb.res, .id = "bs_run") %>% 
    group_by(ethnicity,hba1c.breaks) %>% 
    dplyr::summarise(
      #DPP4.hb = mean(DPP4.hb),
      DPP4.hb.lci = quantile(DPP4.hb,probs=0.025),
      DPP4.hb.uci = quantile(DPP4.hb,probs=0.975),
      #SGLT2.hb = mean(SGLT2.hb),
      SGLT2.hb.lci = quantile(SGLT2.hb,probs=0.025),
      SGLT2.hb.uci = quantile(SGLT2.hb,probs=0.975),
      #sglt2.dpp4.hb.diff = mean(sglt2.dpp4.hb.diff),
      sglt2.dpp4.hb.diff.lci = quantile(sglt2.dpp4.hb.diff,probs=0.025),
      sglt2.dpp4.hb.diff.uci = quantile(sglt2.dpp4.hb.diff,probs=0.975),
    )
  
  
  hb.res.adjusted <- left_join(hb.res.overall,hb.res.ci,by=c("ethnicity","hba1c.breaks"))
  tail(hb.res.adjusted)
  
  #Plot mean and 95% interval of prediction by ethnicity and HbA1c defined subgroup
  white <- hb.res.adjusted %>% filter(ethnicity=="White") 
  asian <- hb.res.adjusted %>% filter(ethnicity=="South Asian")
  black <- hb.res.adjusted %>% filter(ethnicity=="Black") 
  mixed <- hb.res.adjusted %>% filter(ethnicity=="mixed.other") 
  
  L <- list(white,asian,black,mixed)
  
  names(L) <- c("White",
                "Asian",
                "Black",
                "Mixed or Other")
  
  hb.plot <- list()
  
  for(i in 1:4) 
  { 
    #Subgroups by predicted treatment difference
    plotdata <- L[[i]] %>% ungroup() %>% as.data.frame()  %>%
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 1) %>% 
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 5) 
    
    #plot
    coef = data.matrix(cbind(plotdata[,4],plotdata[,3]))
    cim = data.matrix(cbind(plotdata[,8],plotdata[,6]))
    cip = data.matrix(cbind(plotdata[,9],plotdata[,7]))
    hb.plot[[i]] <-
      forestplot(row_names,
                 mean = coef,
                 lower= cim,
                 upper = cip,
                 hrzl_lines = gpar(col="#444444"),lineheight=unit(2,'cm'),
                 is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
                 xticks = tick.hb.resp,
                 zero = 0,
                 #boxsize=0.1,
                 # graphwidth = unit(2,"inches"),
                 # lineheight = unit(0.7,"inches"),
                 ci.vertices=TRUE,
                 col=fpColors(box=c("#f1a340","#4118de"), lines=c("#f1a340","#4118de"), zero = "gray50"),
                 lty.ci = c(1,2),           ,
                 xlab="Average glycaemic response (mmol/mol)",cex=1,
                 title = names(L[i]),
                 new_page = TRUE,
                 fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
                 boxsize = .2, # We set the box size to better visualize the type
                 #line.margin = .1, # We need to add this to avoid crowding
                 txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),summary=gpar(cex = 1), ticks  = gpar(cex = 1), title=gpar(cex=1.3))
                 #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
                 #legend = c("SGLT2i","DPP4i"),
                 #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
                 #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
                 #xlog = TRUE
      )
    
  }  

#Plot together
p.white <- grid2grob(print(hb.plot[[1]]))
p.asian <- grid2grob(print(hb.plot[[2]]))
p.black <- grid2grob(print(hb.plot[[3]]))
p.mixed.other <- grid2grob(print(hb.plot[[4]]))

hb.eth.obs.adj <- (wrap_elements(p.white) + wrap_elements(p.asian)) /
  (wrap_elements(p.black) + wrap_elements(p.mixed.other)) /
  shared_legend +
  plot_layout(height=c(10,10,1))

hb.eth.obs.adj

#Save
grDevices::cairo_pdf(paste0(output_dir,"hb_eth_recal_observed_adjusted.pdf"),width=pdfwidth,height=pdfheight)
hb.eth.obs.adj
dev.off()

png(paste0(output_dir,"hb_eth_recal_observed_adjusted.png"),width=pngwidth,height=pngheight,res=pngres,restoreConsole=TRUE)
hb.eth.obs.adj
dev.off()

write.csv(hb.res.adjusted,file=paste0(output_dir,"hb_eth_recal_observed_adjusted.csv"))



##### R2 and ATE recalibration (exploratory only) ######

#R2 for how well predicted HbA1c explains variation in outcomes HbA1c
w.d <- summary(lm(posthba1c_final~DPP4.pred.lm.uncal,data=final.hb,subset=drugclass=="DPP4" & ethnicity=="White"))$adj.r.squared
a.d <- summary(lm(posthba1c_final~DPP4.pred.lm.uncal,data=final.hb,subset=drugclass=="DPP4" & ethnicity=="South Asian"))$adj.r.squared
b.d <- summary(lm(posthba1c_final~DPP4.pred.lm.uncal,data=final.hb,subset=drugclass=="DPP4" & ethnicity=="Black"))$adj.r.squared
m.d <- summary(lm(posthba1c_final~DPP4.pred.lm.uncal,data=final.hb,subset=drugclass=="DPP4" & (ethnicity=="Mixed" | ethnicity=="Other")))$adj.r.squared
w.s <- summary(lm(posthba1c_final~SGLT2.pred.lm.uncal,data=final.hb,subset=drugclass=="SGLT2" & ethnicity=="White"))$adj.r.squared
a.s <- summary(lm(posthba1c_final~SGLT2.pred.lm.uncal,data=final.hb,subset=drugclass=="SGLT2" & ethnicity=="South Asian"))$adj.r.squared
b.s <- summary(lm(posthba1c_final~SGLT2.pred.lm.uncal,data=final.hb,subset=drugclass=="SGLT2" & ethnicity=="Black"))$adj.r.squared
m.s <- summary(lm(posthba1c_final~SGLT2.pred.lm.uncal,data=final.hb,subset=drugclass=="SGLT2" & (ethnicity=="Mixed" | ethnicity=="Other")))$adj.r.squared

#Miscalibration of prediction (intercept only)
w.d.i <- summary(lm(posthba1c_final-DPP4.pred.lm.uncal~1,data=final.hb,subset=drugclass=="DPP4" & ethnicity=="White")) %>% tidy() %>% select(estimate) %>% as.numeric()
a.d.i <- summary(lm(posthba1c_final-DPP4.pred.lm.uncal~1,data=final.hb,subset=drugclass=="DPP4" & ethnicity=="South Asian")) %>% tidy() %>% select(estimate) %>% as.numeric() 
b.d.i <- summary(lm(posthba1c_final-DPP4.pred.lm.uncal~1,data=final.hb,subset=drugclass=="DPP4" & ethnicity=="Black")) %>% tidy() %>% select(estimate) %>% as.numeric()
m.d.i <- summary(lm(posthba1c_final-DPP4.pred.lm.uncal~1,data=final.hb,subset=drugclass=="DPP4" & (ethnicity=="Mixed" | ethnicity=="Other"))) %>% tidy() %>% select(estimate) %>% as.numeric()
w.s.i <- summary(lm(posthba1c_final-SGLT2.pred.lm.uncal~1,data=final.hb,subset=drugclass=="SGLT2" & ethnicity=="White")) %>% tidy() %>% select(estimate) %>% as.numeric()
a.s.i <- summary(lm(posthba1c_final-SGLT2.pred.lm.uncal~1,data=final.hb,subset=drugclass=="SGLT2" & ethnicity=="South Asian")) %>% tidy() %>% select(estimate) %>% as.numeric()
b.s.i <- summary(lm(posthba1c_final-SGLT2.pred.lm.uncal~1,data=final.hb,subset=drugclass=="SGLT2" & ethnicity=="Black")) %>% tidy() %>% select(estimate) %>% as.numeric()
m.s.i <- summary(lm(posthba1c_final-SGLT2.pred.lm.uncal~1,data=final.hb,subset=drugclass=="SGLT2" & (ethnicity=="Mixed" | ethnicity=="Other"))) %>% tidy() %>% select(estimate) %>% as.numeric()

data.frame(eth=c("White","Asian","Black","Mixed/Other"),
           dpp4.pred.miscal=c(w.d.i,a.d.i,b.d.i,m.d.i),
           sglt2.pred.miscal=c(w.s.i,a.s.i,b.s.i,m.s.i),
           dpp4.pred.r2=c(w.d,a.d,b.d,m.d),
           sglt2.pred.r2=c(w.s,a.s,b.s,m.s))

#R2 and linear predictor coefficient for baseline HbA1c only (to get insight into 'added value' of predictive model)
w.d <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="DPP4" & ethnicity=="White"))$adj.r.squared
a.d <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="DPP4" & ethnicity=="South Asian"))$adj.r.squared
b.d <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="DPP4" & ethnicity=="Black"))$adj.r.squared
m.d <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="DPP4" & (ethnicity=="Mixed" | ethnicity=="Other")))$adj.r.squared
w.s <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="SGLT2" & ethnicity=="White"))$adj.r.squared
a.s <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="SGLT2" & ethnicity=="South Asian"))$adj.r.squared
b.s <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="SGLT2" & ethnicity=="Black"))$adj.r.squared
m.s <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="SGLT2" & (ethnicity=="Mixed" | ethnicity=="Other")))$adj.r.squared

w.d.lp <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="DPP4" & ethnicity=="White")) %>% tidy() %>% select(estimate) %>% slice(2) %>% as.numeric()
a.d.lp <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="DPP4" & ethnicity=="South Asian")) %>% tidy() %>% select(estimate) %>% slice(2) %>% as.numeric() 
b.d.lp <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="DPP4" & ethnicity=="Black")) %>% tidy() %>% select(estimate) %>% slice(2) %>% as.numeric()
m.d.lp <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="DPP4" & (ethnicity=="Mixed" | ethnicity=="Other"))) %>% tidy() %>% slice(2) %>% select(estimate) %>% as.numeric()
w.s.lp <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="SGLT2" & ethnicity=="White")) %>% tidy() %>% select(estimate) %>% slice(2) %>% as.numeric()
a.s.lp <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="SGLT2" & ethnicity=="South Asian")) %>% tidy() %>% select(estimate) %>% slice(2) %>% as.numeric()
b.s.lp <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="SGLT2" & ethnicity=="Black")) %>% tidy() %>% select(estimate) %>% slice(2) %>% as.numeric()
m.s.lp <- summary(lm(posthba1c_final~prehba1cmmol,data=final.hb,subset=drugclass=="SGLT2" & (ethnicity=="Mixed" | ethnicity=="Other"))) %>% tidy() %>% slice(2) %>% select(estimate) %>% as.numeric()

data.frame(eth=c("White","Asian","Black","Mixed/Other"),
           dpp4.hba1c.r2=c(w.d,a.d,b.d,m.d),
           sglt2.hba1c.r2=c(w.s,a.s,b.s,m.s),
           dpp4.hba1c.lp=c(w.d.lp,a.d.lp,b.d.lp,m.d.lp),
           sglt2.hba1c.lp=c(w.s.lp,a.s.lp,b.s.lp,m.s.lp))

##### c-for-benefit (exploratory only) ####
library(MatchIt)

describe(final.hb$hba1c_diff.recal)
describe(final.hb$drugclass)
final.hb <- final.hb %>% mutate(drugdummy = ifelse(drugclass=="SGLT2",1,0))
table(final.hb$drugdummy)

match <- matchit(drugdummy ~ hba1c_diff.recal, 
                 data = final.hb,  
                 exact = ~ ethnicity,
                 method = "nearest",ratio=1, replace=F)

match <- matchit(drugdummy ~ prehba1cmmol + ncurrtx +drugline+hba1cmonth+egfr_ckdepi+prealtlog+agetx+prebmi, 
                 data = final.hb,  
                 exact = ~ ethnicity,
                 method = "nearest",ratio=1, replace=F)

summary(match)
m.nn <- match
nn.match<-match.data(m.nn)

temp <- rownames_to_column(final.hb, var="RowID")

test1 <- tibble(pair.t= as.numeric(row.names(m.nn$match.matrix)), pair.x=as.vector(m.nn$match.matrix))

drug1 <- merge(test1,temp,by.x="pair.t",by.y="RowID")
drug2 <- merge(test1,temp,by.x="pair.x",by.y="RowID")
drug1 <- drug1 %>% select(hba1c_diff.recal.A=hba1c_diff.recal,
                          posthba1c_final.A=posthba1c_final,
                          pair.x=pair.x,
                          pair.t=pair.t)
drug2 <- drug2 %>% select(hba1c_diff.recal.B=hba1c_diff.recal,
                          posthba1c_final.B=posthba1c_final,
                          pair.x=pair.x)
cfb <- merge(drug1,drug2,by="pair.x")
head(cfb)
#describe(cfb)

cfb <- cfb %>% mutate(pred.ben.avg= (hba1c_diff.recal.A+hba1c_diff.recal.B)/2,
                      obs.ben= posthba1c_final.B-posthba1c_final.A) %>% arrange(pair.t)
plot(cfb$pred.ben.avg,cfb$obs.ben)

cfb <- cfb %>% mutate(pred.ben.avg= (hba1c_diff.recal.A+hba1c_diff.recal.B)/2,
                      obs.ben= ifelse(posthba1c_final.A-posthba1c_final.B<0,-1,
                                      ifelse(posthba1c_final.A-posthba1c_final.B==0,0,1)))
cfb <- cfb %>% mutate(pred.ben.avg= (hba1c_diff.recal.A+hba1c_diff.recal.B)/2,
                      obs.ben= ifelse(posthba1c_final.B-posthba1c_final.A<0,1,
                                      ifelse(posthba1c_final.B-posthba1c_final.A==0,0,-1)))

head(cfb)
table(cfb$obs.ben)


# Benefit c-statistic
cindex <- rcorr.cens(cfb$pred.ben.avg, cfb$obs.ben)
c.benefit <- cindex["C Index"][[1]]
c.benefit.se <- cindex["S.D."][[1]]/2	# The sd of the c-index is half the sd of Dxy

c.benefit
c.benefit - 1.96*c.benefit.se
c.benefit + 1.96*c.benefit.se

# Benefit c-statistic 5 mmol/mol SGLT2
cfb.sglt5 <- cfb %>% filter(pred.ben.avg< -5)
cindex <- rcorr.cens(cfb.sglt5$pred.ben.avg, cfb.sglt5$obs.ben)
c.benefit <- cindex["C Index"][[1]]
c.benefit.se <- cindex["S.D."][[1]]/2	# The sd of the c-index is half the sd of Dxy

c.benefit
c.benefit - 1.96*c.benefit.se
c.benefit + 1.96*c.benefit.se


library(Matching)
set.seed(1)   # For reproducibility
rr <- Match(Tr=final.hb$drugdummy, X=final.hb$hba1c_diff.recal, M=1,ties=F,replace=FALSE)
ind.A <- rr$index.control
ind.B <- rr$index.treated

### Calculation of predicted and observed benefit in matched pairs
pred.ben.A <- pred.ben[ind.A]
pred.ben.B <- pred.ben[ind.B]
pred.ben.avg <- (pred.ben.A+pred.ben.B)/2
obs.out.A <- y[ind.A]
obs.out.B <- y[ind.B]
obs.ben <- obs.out.A-obs.out.B

# Matching plot
plot(pred.ben.A, pred.ben.B)
lines(c(-1,1),c(-1,1),lwd=2)

# Benefit c-statistic
cindex <- rcorr.cens(pred.ben.avg, obs.ben)
c.benefit <- cindex["C Index"][[1]]
c.benefit.se <- cindex["S.D."][[1]]/2	# Half the sd of Dxy
c.benefit	
# [1] 0.6212		# The c-for-benefit
c.benefit - 1.96*c.benefit.se
# [1] 0.578976		# The 95% lower bound of the c-for-benefit
c.benefit + 1.96*c.benefit.se	
# [1] 0.6634239		# The 95% upper bound of the c-for-benefit



# test2 <- test1 %>% 

#   inner_join(nn.match  %>% mutate(pair.t=row_number()))%>% na.omit
# 
# tibble(pair.t= as.numeric(row.names(m.nn$match.matrix)), pair.x=as.vector(m.nn$match.matrix)) %>% 
#   inner_join(nn.match  %>% mutate(pair.t=row_number())) %>% na.omit
# 
# # Output matched individuals
# # NB easiest way to pair up with their case is to use the R rownumber which is constant across the tables
# m.out <- data.frame(match$match.matrix) %>%
#   rownames_to_column(var="CaseRowID") 
# 
# # Reshape to long
# temp <- rownames_to_column(final.hb, var="RowID")
# m.out <- merge(m.out, temp["RowID"], by.x="CaseRowID", by.y="RowID", all.x=F) %>%
#   as_tibble()
# head(m.out)
# 
# m.out_long <- gather(data=m.out,key=MatchNumber,value=RowID) %>%
#   select(-MatchNumber) 
# 
# # Matched control dataset
# matchedcontrols  <- merge(temp,m.out_long,by="RowID",all.y=TRUE)
# table(matchedcontrols$drugclass)
# 
# head(matchedcontrols)

##### Discontinuation #####

#Set up 

  #Set dummy hba1cmonth at 6 months
  final.dc <- final.dc %>% dplyr::mutate(hba1cmonth=6)
  
  #Predict HbA1c outcome uncalibrated
  final.dc <- final.dc %>% mutate(drug=drugclass,
                                  drugclass="DPP4")
  final.dc$DPP4.pred.lm.uncal <- predict(m1,final.dc)
  final.dc <- final.dc %>% mutate(drugclass="SGLT2")
  final.dc$SGLT2.pred.lm.uncal <- predict(m1,final.dc)
  final.dc <- final.dc %>% 
    mutate(hba1c_diff.uncal = SGLT2.pred.lm.uncal-DPP4.pred.lm.uncal,
           bestdrug.uncal=ifelse(hba1c_diff.uncal<=0,"SGLT2","DPP4"),
           drugclass=drug)
  head(final.dc)
 
  #Recalibrate HbA1c outcome to AURUM
  final.dc <- final.dc %>% mutate(DPP4.pred.lm.recal=ifelse(ethnicity == "White",DPP4.pred.lm.uncal+ctfm$intercept[1],
                                                            ifelse(ethnicity == "South Asian",DPP4.pred.lm.uncal+ctfm$intercept[3],
                                                                   ifelse(ethnicity == "Black",DPP4.pred.lm.uncal+ctfm$intercept[5],
                                                                          ifelse(ethnicity == "mixed.other",DPP4.pred.lm.uncal+ctfm$intercept[7],NA
                                                                          )))),
                                  SGLT2.pred.lm.recal=ifelse(ethnicity == "White",SGLT2.pred.lm.uncal+ctfm$intercept[2],SGLT2.pred.lm.uncal),
                                  hba1c_diff.recal = SGLT2.pred.lm.recal-DPP4.pred.lm.recal,
                                  bestdrug.recal=ifelse(hba1c_diff.recal<=0,"SGLT2","DPP4"))
  
  #Generate other variables for modelling
  
  #Define hba1c.breaks based on recalibrated HbA1c outcome
  final.dc <- final.dc %>% mutate(hba1c.breaks = cut(hba1c_diff.recal, breaks=c(min(hba1c_diff.recal,na.rm=T)-0.01,-5,-3,0,3,max(hba1c_diff.recal,na.rm=T)+0.01)))
  describe(final.dc$hba1c.breaks)
  
  #Define datadist for modelling
  ddist <- datadist(final.dc); options(datadist='ddist') 
  
#Predicted discontination from GOLD model (after HbA1c recalibration)  
  
  #Rename to allow fit of GOLD model
  final.dc <- final.dc %>% dplyr::rename(hba1c_diff="hba1c_diff.recal",
                                         bestdrug="bestdrug.recal")

  #Predict probabilities on each drug using GOLD model
  final.dc$drug <- "DPP4"
  final.dc$preddisonDPP4 <- plogis(predict(m.dis, newdata=final.dc, se.fit=F))
  final.dc$drug <- "SGLT2"
  final.dc$preddisonSGLT2 <- plogis(predict(m.dis, newdata=final.dc, se.fit=F))
  final.dc$drug <- final.dc$drugclass
  final.dc$preddisdiffsglt2dpp4 <- abs(final.dc$preddisonSGLT2-final.dc$preddisonDPP4)
  describe(final.dc$preddisdiffsglt2dpp4)
  hist(final.dc$preddisdiffsglt2dpp4)

  #Plot mean and 95% interval of prediction by ethnicity and HbA1c defined subgroup
  white <- final.dc %>% filter(ethnicity=="White")
  asian <- final.dc %>% filter(ethnicity=="South Asian")
  black <- final.dc %>% filter(ethnicity=="Black")
  mixed.other <- final.dc %>% filter(ethnicity=="mixed.other")
  
  L <- list(white,asian,black,mixed.other)
  
  names(L) <- c("White",
                "Asian",
                "Black",
                "Mixed or Other")
  
  dc.res <- list()
  dc.plot <- list()
  
  for(i in 1:4) 
  { 
    #Subgroups by predicted treatment difference
    data  <- L[[i]] %>% dplyr::filter(hba1c_diff<= -5) 
    # df <- expand.grid(quantile(data$preddisonSGLT2,c(0.5)),quantile(data$preddisonSGLT2,c(0.25)),quantile(data$preddisonSGLT2,c(0.75)),
    #                   quantile(data$preddisonDPP4,c(0.5)),quantile(data$preddisonDPP4,c(0.25)),quantile(data$preddisonDPP4,c(0.75)))
    df <- expand.grid(mean(data$preddisonSGLT2),quantile(data$preddisonSGLT2,c(0.025)),quantile(data$preddisonSGLT2,c(0.975)),
                      mean(data$preddisonDPP4),quantile(data$preddisonDPP4,c(0.025)),quantile(data$preddisonDPP4,c(0.975)))
    SGLT2.5 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
    
    data  <- L[[i]] %>% dplyr::filter(hba1c_diff<= -3 & hba1c_diff > -5) 
    df <- expand.grid(mean(data$preddisonSGLT2),quantile(data$preddisonSGLT2,c(0.025)),quantile(data$preddisonSGLT2,c(0.975)),
                      mean(data$preddisonDPP4),quantile(data$preddisonDPP4,c(0.025)),quantile(data$preddisonDPP4,c(0.975)))
    SGLT2.3 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
    
    data  <- L[[i]] %>% dplyr::filter(hba1c_diff> -3 & hba1c_diff <= 0) 
    df <- expand.grid(mean(data$preddisonSGLT2),quantile(data$preddisonSGLT2,c(0.025)),quantile(data$preddisonSGLT2,c(0.975)),
                      mean(data$preddisonDPP4),quantile(data$preddisonDPP4,c(0.025)),quantile(data$preddisonDPP4,c(0.975)))
    SGLT2.0 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
    
    data  <- L[[i]] %>% dplyr::filter(hba1c_diff> 0 & hba1c_diff < 3) 
    df <- expand.grid(mean(data$preddisonSGLT2),quantile(data$preddisonSGLT2,c(0.025)),quantile(data$preddisonSGLT2,c(0.975)),
                      mean(data$preddisonDPP4),quantile(data$preddisonDPP4,c(0.025)),quantile(data$preddisonDPP4,c(0.975)))
    DPP4.0 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
    
    data  <- L[[i]] %>% dplyr::filter(hba1c_diff>= 3) 
    df <- expand.grid(mean(data$preddisonSGLT2),quantile(data$preddisonSGLT2,c(0.025)),quantile(data$preddisonSGLT2,c(0.975)),
                      mean(data$preddisonDPP4),quantile(data$preddisonDPP4,c(0.025)),quantile(data$preddisonDPP4,c(0.975)))
    DPP4.3 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
    
    final.dis.pred <- rbind(SGLT2.5,SGLT2.3,SGLT2.0,DPP4.0,DPP4.3)
    df.lab <- c("SGLT2.5","SGLT2.3","SGLT2.0","DPP4.0","DPP4.3")
    final.dis.pred <- cbind(df.lab,final.dis.pred)
    final.dis.pred.p <- final.dis.pred %>% 
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 1) %>% 
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 5) 
    
    #res table
    dc.res[[i]] <- final.dis.pred.p
    
    #plot
    coef = data.matrix(cbind(final.dis.pred.p[,2]*100,final.dis.pred.p[,5]*100))
    cim = data.matrix(cbind(final.dis.pred.p[,3]*100,final.dis.pred.p[,6]*100))
    cip = data.matrix(cbind(final.dis.pred.p[,4]*100,final.dis.pred.p[,7]*100))
    dc.plot[[i]] <-
      forestplot(row_names,
                 mean = coef,
                 lower= cim,
                 upper = cip,,
                 hrzl_lines = gpar(col="#444444"),lineheight=unit(2,'cm'),
                 is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
                 xticks = tick.dc,
                 zero = 0,
                 #boxsize=0.1,
                 # graphwidth = unit(2,"inches"),
                 # lineheight = unit(0.7,"inches"),
                 ci.vertices=TRUE,
                 col=fpColors(box=c("#f1a340","#4118de"), lines=c("#f1a340","#4118de"), zero = "gray50"),
                 lty.ci = c(1,2),           ,
                 xlab="Proportion discontinuing therapy (%)",cex=1,
                 title = names(L[i]),
                 new_page = TRUE,
                 fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
                 boxsize = .2, # We set the box size to better visualize the type
                 #line.margin = .1, # We need to add this to avoid crowding
                 txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),summary=gpar(cex = 1), ticks  = gpar(cex = 1), title=gpar(cex=1.3))
                 #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
                 #legend = c("SGLT2i","DPP4i"),
                 #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
                 #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
                 #xlog = TRUE
      )
    
  }  
  
  #Plot together
  p.white <- grid2grob(print(dc.plot[[1]]))
  p.asian <- grid2grob(print(dc.plot[[2]]))
  p.black <- grid2grob(print(dc.plot[[3]]))
  p.mixed.other <- grid2grob(print(dc.plot[[4]]))
  
  dc.eth <- (wrap_elements(p.white) + wrap_elements(p.asian)) /
    (wrap_elements(p.black) + wrap_elements(p.mixed.other)) /
    shared_legend +
    plot_layout(height=c(10,10,1))
  
  dc.eth
  dc.res[1]
  dc.res[2]
  dc.res[3]
  dc.res[4]
  
  grDevices::cairo_pdf(paste0(output_dir,"dc_eth_recal_predicted.pdf"),width=pdfwidth,height=pdfheight)
  dc.eth
  dev.off()
  
  png(paste0(output_dir,"dc_eth_recal_predicted.png"),width=pngwidth,height=pngheight,res=pngres,restoreConsole=TRUE)
  dc.eth
  dev.off()
  
  dc.res[[1]]$eth = "White"
  dc.res[[2]]$eth = "Asian"
  dc.res[[3]]$eth = "Black"
  dc.res[[4]]$eth = "Mixed or Other"
  dc.res.pred <- dc.res %>% bind_rows(dc.res) %>% filter(!is.na(df.lab))
  write.csv(dc.res.pred,file=paste0(output_dir,"dc_eth_recal_predicted.csv"))
  
  #rename table outputs to compare obs vs pred
  dc.res.pred <- dc.res

#### Observed discontinuation, by ethnicity and HbA1c benefit subgroup

  # final.dc <- final.dc %>% mutate(drugline=as.factor(drugline))
  # ddist <- datadist(final.dc); options(datadist='ddist') 
  # 
  # #Model
  # m.dis.obs <- lrm(stopdrug_6m_3mFU ~ drug*rcs(hba1c_diff,3)*ethnicity + drugline + ncurrtx, data=final.dc,x=T,y=T,se.fit=TRUE)
  # m.dis.obs
  # 
  # #Predict probabilities on each drug de novo
  # final.dc$drug <- "DPP4"
  # final.dc$DPP4.dc <- plogis(predict(m.dis.obs, newdata=final.dc, se.fit=F))
  # final.dc$drug <- "SGLT2"
  # final.dc$SGLT2.dc <- plogis(predict(m.dis.obs, newdata=final.dc, se.fit=F))
  # final.dc$drug <- final.dc$drugclass
  # final.dc$sglt2.dpp4.dc.diff <- abs(final.dc$SGLT2.dc-final.dc$DPP4.dc)
  # describe(final.dc$sglt2.dpp4.dc.diff)
  # hist(final.dc$sglt2.dpp4.dc.diff)
  # 
  # #Plot mean and 95% interval of prediction by ethnicity and HbA1c defined subgroup
  # 
  # white <- final.dc %>% filter(ethnicity=="White") 
  # asian <- final.dc %>% filter(ethnicity=="South Asian")
  # black <- final.dc %>% filter(ethnicity=="Black") 
  # mixed.other <- final.dc %>% filter(ethnicity=="Mixed"|ethnicity=="Other") 
  # 
  # L <- list(white,asian,black,mixed.other)
  # 
  # names(L) <- c("White",
  #               "Asian",
  #               "Black",
  #               "Mixed or Other")
  # 
  # dc.res <- list()
  # dc.plot <- list()
  # 
  # for(i in 1:4) 
  # { 
  #   #Subgroups by predicted treatment difference
  #   data  <- L[[i]] %>% dplyr::filter(hba1c_diff<= -5) 
  #   # df <- expand.grid(quantile(data$SGLT2.dc,c(0.5)),quantile(data$SGLT2.dc,c(0.25)),quantile(data$SGLT2.dc,c(0.75)),
  #   #                   quantile(data$DPP4.dc,c(0.5)),quantile(data$DPP4.dc,c(0.25)),quantile(data$DPP4.dc,c(0.75)))
  #   df <- expand.grid(mean(data$SGLT2.dc),quantile(data$SGLT2.dc,c(0.025)),quantile(data$SGLT2.dc,c(0.975)),
  #                     mean(data$DPP4.dc),quantile(data$DPP4.dc,c(0.025)),quantile(data$DPP4.dc,c(0.975)))
  #   SGLT2.5 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
  #   
  #   data  <- L[[i]] %>% dplyr::filter(hba1c_diff<= -3 & hba1c_diff > -5) 
  #   df <- expand.grid(mean(data$SGLT2.dc),quantile(data$SGLT2.dc,c(0.025)),quantile(data$SGLT2.dc,c(0.975)),
  #                     mean(data$DPP4.dc),quantile(data$DPP4.dc,c(0.025)),quantile(data$DPP4.dc,c(0.975)))
  #   SGLT2.3 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
  #   
  #   data  <- L[[i]] %>% dplyr::filter(hba1c_diff> -3 & hba1c_diff <= 0) 
  #   df <- expand.grid(mean(data$SGLT2.dc),quantile(data$SGLT2.dc,c(0.025)),quantile(data$SGLT2.dc,c(0.975)),
  #                     mean(data$DPP4.dc),quantile(data$DPP4.dc,c(0.025)),quantile(data$DPP4.dc,c(0.975)))
  #   SGLT2.0 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
  #   
  #   data  <- L[[i]] %>% dplyr::filter(hba1c_diff> 0 & hba1c_diff < 3) 
  #   df <- expand.grid(mean(data$SGLT2.dc),quantile(data$SGLT2.dc,c(0.025)),quantile(data$SGLT2.dc,c(0.975)),
  #                     mean(data$DPP4.dc),quantile(data$DPP4.dc,c(0.025)),quantile(data$DPP4.dc,c(0.975)))
  #   DPP4.0 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
  #   
  #   data  <- L[[i]] %>% dplyr::filter(hba1c_diff>= 3) 
  #   df <- expand.grid(mean(data$SGLT2.dc),quantile(data$SGLT2.dc,c(0.025)),quantile(data$SGLT2.dc,c(0.975)),
  #                     mean(data$DPP4.dc),quantile(data$DPP4.dc,c(0.025)),quantile(data$DPP4.dc,c(0.975)))
  #   DPP4.3 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
  #   
  #   final.dis <- rbind(SGLT2.5,SGLT2.3,SGLT2.0,DPP4.0,DPP4.3)
  #   df.lab <- c("SGLT2.5","SGLT2.3","SGLT2.0","DPP4.0","DPP4.3")
  #   final.dis6m <- cbind(df.lab,final.dis)
  #   final.dis6m.p <- final.dis6m %>% 
  #     add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 1) %>% 
  #     add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 5) 
  #   
  #   #res table
  #   dc.res[[i]] <- final.dis6m.p
  #   
  #   #plot
  #   coef = data.matrix(cbind(final.dis6m.p[,2]*100,final.dis6m.p[,5]*100))
  #   cim = data.matrix(cbind(final.dis6m.p[,3]*100,final.dis6m.p[,6]*100))
  #   cip = data.matrix(cbind(final.dis6m.p[,4]*100,final.dis6m.p[,7]*100))
  #   dc.plot[[i]] <-
  #     forestplot(row_names,
  #              mean = coef,
  #              lower= cim,
  #              upper = cip,,
  #              hrzl_lines = gpar(col="#444444"),lineheight=unit(2,'cm'),
  #              is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
  #              xticks = tick,
  #              zero = 0,
  #              #boxsize=0.1,
  #              # graphwidth = unit(2,"inches"),
  #              # lineheight = unit(0.7,"inches"),
  #              ci.vertices=TRUE,
  #              col=fpColors(box=c("#f1a340","#4118de"), lines=c("#f1a340","#4118de"), zero = "gray50"),
  #              lty.ci = c(1,2),           ,
  #              xlab="Early therapy discontinuation (%)",cex=1,
  #              title = names(L[i]),
  #              new_page = TRUE,
  #              fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
  #              boxsize = .2, # We set the box size to better visualize the type
  #              #line.margin = .1, # We need to add this to avoid crowding
  #              txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),summary=gpar(cex = 1), ticks  = gpar(cex = 1), title=gpar(cex=1.3))
  #              #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
  #              #legend = c("SGLT2i","DPP4i"),
  #              #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
  #              #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
  #              #xlog = TRUE
  #   )
  #   
  # }  
  #   
  # #Plot together
  # p.white <- grid2grob(print(dc.plot[[1]]))
  # p.asian <- grid2grob(print(dc.plot[[2]]))
  # p.black <- grid2grob(print(dc.plot[[3]]))
  # p.mixed.other <- grid2grob(print(dc.plot[[4]]))
  # 
  # dc.eth <- (wrap_elements(p.white) + wrap_elements(p.asian)) /
  # (wrap_elements(p.black) + wrap_elements(p.mixed.other))
  # 
  # dc.eth
  # dc.res[1]
  # dc.res[2]
  # dc.res[3]
  # dc.res[4]
  
  #### Unadjusted
  #Note this does not correct for potential imbalance amongst patients within each interval receiving each therapy
  
  # dc.res.unadjusted <- final.dc %>% 
  #   dplyr::group_by(ethnicity,hba1c.breaks,drugclass) %>% 
  #   dplyr::summarise(n=length(stopdrug_6m_3mFU),
  #                    n.discontinued = sum(stopdrug_6m_3mFU),
  #                    prop = sum(stopdrug_6m_3mFU) / length(stopdrug_6m_3mFU),
  #                    lower = lapply(sum(stopdrug_6m_3mFU), prop.test, n = length(stopdrug_6m_3mFU)),
  #                    upper = sapply(lower, function(x) x$conf.int[2]), 
  #                    lower = sapply(lower, function(x) x$conf.int[1])) %>%
  #   select(ethnicity,hba1c.breaks,drugclass,n,n.discontinued,prop,lower,upper)
  
  #Harrell method
  dc.res.unadjusted <- final.dc %>% 
    dplyr::group_by(ethnicity,hba1c.breaks,drugclass) %>% 
    dplyr::summarise(n=length(stopdrug_6m_3mFU),
                     n.discontinued = sum(stopdrug_6m_3mFU),
                     #prop = sum(stopdrug_6m_3mFU) / length(stopdrug_6m_3mFU),
                     prop.conf = binconf(x=sum(stopdrug_6m_3mFU), n = length(stopdrug_6m_3mFU)),
                     prop=prop.conf[1],
                     lci=prop.conf[2],
                     uci=prop.conf[3]) %>%
    select(ethnicity,hba1c.breaks,drugclass,n,n.discontinued,prop,lci,uci)
  
  #Long to wide for plotting
  dc.res.unadjusted <- dc.res.unadjusted %>% pivot_wider(
    names_from = drugclass,
    values_from = c(n,n.discontinued,prop,lci,uci)
  )
  
  #Plot observed and 95% CI by ethnicity and HbA1c defined subgroup
  
  white <- dc.res.unadjusted %>% filter(ethnicity=="White") 
  asian <- dc.res.unadjusted %>% filter(ethnicity=="South Asian")
  black <- dc.res.unadjusted %>% filter(ethnicity=="Black") 
  mixed <- dc.res.unadjusted %>% filter(ethnicity=="mixed.other") 

  L <- list(white,asian,black,mixed)
  
  names(L) <- c("White",
                "Asian",
                "Black",
                "Mixed or Other")
  
  dc.plot <- list()
  

  for(i in 1:4) 
  { 
    #Subgroups by predicted treatment difference
    plotdata <- L[[i]] %>% ungroup() %>% as.data.frame()  %>%
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 1) %>% 
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 5) 
    
    #plot
    coef = data.matrix(cbind(plotdata[,8]*100,plotdata[,7]*100))
    cim = data.matrix(cbind(plotdata[,10]*100,plotdata[,9]*100))
    cip = data.matrix(cbind(plotdata[,12]*100,plotdata[,11]*100))
    dc.plot[[i]] <-
      forestplot(row_names,
                 mean = coef,
                 lower= cim,
                 upper = cip,
                 hrzl_lines = gpar(col="#444444"),lineheight=unit(2,'cm'),
                 is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
                 xticks = tick.dc,
                 zero = 0,
                 #boxsize=0.1,
                 # graphwidth = unit(2,"inches"),
                 # lineheight = unit(0.7,"inches"),
                 ci.vertices=TRUE,
                 col=fpColors(box=c("#f1a340","#4118de"), lines=c("#f1a340","#4118de"), zero = "gray50"),
                 lty.ci = c(1,2),           ,
                 xlab="Proportion discontinuing therapy (%)",cex=1,
                 title = names(L[i]),
                 new_page = TRUE,
                 fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
                 boxsize = .2, # We set the box size to better visualize the type
                 #line.margin = .1, # We need to add this to avoid crowding
                 txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),summary=gpar(cex = 1), ticks  = gpar(cex = 1), title=gpar(cex=1.3))
                 #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
                 #legend = c("SGLT2i","DPP4i"),
                 #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
                 #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
                 #xlog = TRUE
      )
    
  }  
  
  
  #Plot together
  p.dc.white <- grid2grob(print(dc.plot[[1]]))
  p.dc.asian <- grid2grob(print(dc.plot[[2]]))
  p.dc.black <- grid2grob(print(dc.plot[[3]]))
  p.dc.mixed.other <- grid2grob(print(dc.plot[[4]]))
  
  dc.eth.obs <- (wrap_elements(p.dc.white) + wrap_elements(p.dc.asian)) /
    (wrap_elements(p.dc.black) + wrap_elements(p.dc.mixed.other)) /
    shared_legend +
    plot_layout(height=c(10,10,1))
  
  dc.eth.obs
  
  grDevices::cairo_pdf(paste0(output_dir,"dc_eth_recal_observed_unadjusted.pdf"),width=pdfwidth,height=pdfheight)
  dc.eth.obs
  dev.off()
  
  png(paste0(output_dir,"dc_eth_recal_observed_unadjusted.png"),width=pngwidth,height=pngheight,res=pngres,restoreConsole=TRUE)
  dc.eth.obs
  dev.off()
  
  write.csv(dc.res.unadjusted,file=paste0(output_dir,"dc_eth_recal_observed_unadjusted.csv"))

  #rename table outputs to compare obs vs pred
  dc.res.obs <- dc.res.unadjusted
  
  #### Adjusted
  
  #Mean prediction
  #m.dis.obs <- lrm(stopdrug_6m_3mFU ~ drug*rcs(hba1c_diff,3)*ethnicity + drugline + ncurrtx, data=final.dc,x=T,y=T,se.fit=TRUE)
  #m.dis.obs <- lrm(stopdrug_6m_3mFU ~ drug*ethnicity + drugline + ncurrtx, data=final.dc,x=T,y=T,se.fit=TRUE)
  # m.dis.obs <- lrm(stopdrug_6m_3mFU ~ drug*ethnicity + drugline + ncurrtx + 
  #                    drug*agetx+
  #                    drug*prebmi+
  #                    drug*prehba1cmmol+
  #                    drug*prealtlog+
  #                    drug*egfr_ckdepi, data=final.dc,x=T,y=T,se.fit=TRUE)
  # m.dis.obs
  m.dis.obs <- lrm(stopdrug_6m_3mFU ~ drug*ethnicity*hba1c.breaks + drugline + ncurrtx, data=final.dc,x=T,y=T,se.fit=TRUE)
  
  #The above models show the predicted HbA1c difference or covariate * drug interactions are required as an input variable to stratify discontinuation
  
  #Predict probabilities on each drug de novo
  final.dc$drug <- "DPP4"
  final.dc$DPP4.dc <- plogis(predict(m.dis.obs, newdata=final.dc, se.fit=F))
  final.dc$drug <- "SGLT2"
  final.dc$SGLT2.dc <- plogis(predict(m.dis.obs, newdata=final.dc, se.fit=F))
  final.dc$drug <- final.dc$drugclass
  final.dc$sglt2.dpp4.dc.diff <- abs(final.dc$SGLT2.dc-final.dc$DPP4.dc)
  describe(final.dc$sglt2.dpp4.dc.diff)
  hist(final.dc$sglt2.dpp4.dc.diff)
  
  dc.res.overall <- final.dc %>% 
    group_by(ethnicity,hba1c.breaks) %>% 
    dplyr::summarise(
      DPP4.dc = mean(DPP4.dc),
      SGLT2.dc = mean(SGLT2.dc),
      sglt2.dpp4.dc.diff = mean(sglt2.dpp4.dc.diff)
    )
  

  n=nrow(final.dc)
  
  dc.res <- list()
  
  for(b in 1:B){
    i = sample(x = 1:n, size = n, replace = TRUE) ## sample indices
    temp = final.dc[i,] ## temp data set
    #temp_model =  lrm(stopdrug_6m_3mFU ~ drug*rcs(hba1c_diff,3)*ethnicity + drugline + ncurrtx, data=temp,x=T,y=T,se.fit=TRUE)
    temp_model =  lrm(stopdrug_6m_3mFU ~ drug*hba1c.breaks*ethnicity + drugline + ncurrtx, data=temp,x=T,y=T,se.fit=TRUE)
    
    temp$drug <- "DPP4"
    temp$DPP4.dc <- plogis(predict(temp_model, newdata=temp, se.fit=F))
    temp$drug <- "SGLT2"
    temp$SGLT2.dc <- plogis(predict(temp_model, newdata=temp, se.fit=F))
    temp$sglt2.dpp4.dc.diff <- abs(temp$SGLT2.dc-temp$DPP4.dc)
    
    res <- temp %>% 
      group_by(ethnicity,hba1c.breaks) %>% 
      dplyr::summarise(
        DPP4.dc = mean(DPP4.dc),
        SGLT2.dc = mean(SGLT2.dc),
        sglt2.dpp4.dc.diff = mean(sglt2.dpp4.dc.diff)
      )
    
    #res table
    dc.res[[b]] <- res
  }
  
  #Derive the bootstrapped CIs for the mean predictions for each subgrup
  
  dc.res.ci <- bind_rows(dc.res, .id = "bs_run") %>% 
    group_by(ethnicity,hba1c.breaks) %>% 
    dplyr::summarise(
      #DPP4.dc = mean(DPP4.dc),
      DPP4.dc.lci = quantile(DPP4.dc,probs=0.025),
      DPP4.dc.uci = quantile(DPP4.dc,probs=0.975),
      #SGLT2.dc = mean(SGLT2.dc),
      SGLT2.dc.lci = quantile(SGLT2.dc,probs=0.025),
      SGLT2.dc.uci = quantile(SGLT2.dc,probs=0.975),
      #sglt2.dpp4.dc.diff = mean(sglt2.dpp4.dc.diff),
      sglt2.dpp4.dc.diff.lci = quantile(sglt2.dpp4.dc.diff,probs=0.025),
      sglt2.dpp4.dc.diff.uci = quantile(sglt2.dpp4.dc.diff,probs=0.975),
    )
  
  
  dc.res.adjusted <- left_join(dc.res.overall,dc.res.ci,by=c("ethnicity","hba1c.breaks"))
  tail(dc.res.adjusted)
  
  #Plot mean and 95% interval of prediction by ethnicity and HbA1c defined subgroup
  
  white <- dc.res.adjusted %>% filter(ethnicity=="White") 
  asian <- dc.res.adjusted %>% filter(ethnicity=="South Asian")
  black <- dc.res.adjusted %>% filter(ethnicity=="Black") 
  mixed <- dc.res.adjusted %>% filter(ethnicity=="mixed.other") 
  #mixed <- dc.res.adjusted %>% filter(ethnicity=="Mixed") 
  #other <- dc.res.adjusted %>% filter(ethnicity=="Other") 
  
  L <- list(white,asian,black,mixed)
  
  names(L) <- c("White",
                "Asian",
                "Black",
                "Mixed or Other")
  
  # L <- list(white,asian,black,mixed,other)
  # 
  # names(L) <- c("White",
  #               "Asian",
  #               "Black",
  #               "Mixed",
  #               "Other")
  
  dc.plot <- list()
  
  for(i in 1:4) 
  { 
    #Subgroups by predicted treatment difference
    plotdata <- L[[i]] %>% ungroup() %>% as.data.frame()  %>%
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 1) %>% 
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 5) 
    
    #plot
    coef = data.matrix(cbind(plotdata[,4]*100,plotdata[,3]*100))
    cim = data.matrix(cbind(plotdata[,8]*100,plotdata[,6]*100))
    cip = data.matrix(cbind(plotdata[,9]*100,plotdata[,7]*100))
    dc.plot[[i]] <-
      forestplot(row_names,
                 mean = coef,
                 lower= cim,
                 upper = cip,
                 hrzl_lines = gpar(col="#444444"),lineheight=unit(2,'cm'),
                 is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
                 xticks = tick.dc,
                 zero = 0,
                 #boxsize=0.1,
                 # graphwidth = unit(2,"inches"),
                 # lineheight = unit(0.7,"inches"),
                 ci.vertices=TRUE,
                 col=fpColors(box=c("#f1a340","#4118de"), lines=c("#f1a340","#4118de"), zero = "gray50"),
                 lty.ci = c(1,2),           ,
                 xlab="Proportion discontinuing therapy (%)",cex=1,
                 title = names(L[i]),
                 new_page = TRUE,
                 fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
                 boxsize = .2, # We set the box size to better visualize the type
                 #line.margin = .1, # We need to add this to avoid crowding
                 txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),summary=gpar(cex = 1), ticks  = gpar(cex = 1), title=gpar(cex=1.3))
                 #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
                 #legend = c("SGLT2i","DPP4i"),
                 #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
                 #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
                 #xlog = TRUE
      )
    
  }  
  
  
  #Plot together
  p.dc.white <- grid2grob(print(dc.plot[[1]]))
  p.dc.asian <- grid2grob(print(dc.plot[[2]]))
  p.dc.black <- grid2grob(print(dc.plot[[3]]))
  p.dc.mixed.other <- grid2grob(print(dc.plot[[4]]))
  
  dc.eth.obs.adj <- (wrap_elements(p.dc.white) + wrap_elements(p.dc.asian)) /
    (wrap_elements(p.dc.black) + wrap_elements(p.dc.mixed.other)) /
    shared_legend +
    plot_layout(height=c(10,10,1))
  
  dc.eth.obs.adj
  
  grDevices::cairo_pdf(paste0(output_dir,"dc_eth_recal_observed_adjusted.pdf"),width=pdfwidth,height=pdfheight)
  dc.eth.obs.adj
  dev.off()
  
  png(paste0(output_dir,"dc_eth_recal_observed_adjusted.png"),width=pngwidth,height=pngheight,res=pngres,restoreConsole=TRUE)
  dc.eth.obs.adj
  dev.off()

  write.csv(dc.res.adjusted,file=paste0(output_dir,"dc_eth_recal_observed_adjusted.csv"))
  
 ####Explore obs vs pred difference
  # dc.res.white.obs <- dc.res.obs %>% 
  #   filter(ethnicity=="South Asian") %>% 
  #   mutate(obs=prop_DPP4-prop_SGLT2) %>%
  #   select(ethnicity,hba1c.breaks,obs)
  # #Add CIs for contrast
  # 
  # # dc.res.white.obs <- data.frame(Ethnicity="Asian",
  # #                            model="obs",
  # #                            dc.res.obs[4]) %>%
  # #   mutate(obs=DPP4.est-SGLT2.est,
  # #          obs.l = DPP4.lci-SGLT2.lci,
  # #          obs.u = DPP4.uci-SGLT2.uci) %>% 
  # #   select(Ethnicity,df.lab,model,obs,obs.l,obs.u) %>%
  # #   filter(!is.na(obs))
  # 
  # dc.res.white.pred <- data.frame(dc.res.pred[4]) %>% 
  #   mutate(pred=DPP4.est-SGLT2.est) %>% 
  #   select(pred) %>%
  #   filter(!is.na(pred))
  # 
  # dc.res.white <- cbind(dc.res.white.obs,dc.res.white.pred) 
  # 
  # ggplot(dc.res.white,aes(x=pred*100,y=obs*100)) +
  #   geom_vline(xintercept=0, linetype="dashed", color = "grey60") + 
  #   geom_hline(yintercept=0, linetype="dashed", color = "grey60") +
  #   geom_abline(intercept=0,slope=1, color="red", lwd=0.75) + ggtitle("") +
  #   geom_point(alpha=1) + theme_classic()  +
  #   #geom_errorbar(aes(ymin=obs.l, ymax=obs.u), colour="black", width=.1)  +
  #   ylab("Observed discontinuation)") + xlab("Predicted discontinuation") +
  #   theme_base(base_size = 8)  +
  #   scale_x_continuous(limits=c(-25,5)) +
  #   scale_y_continuous(limits=c(-25,5)) +
  #   theme(plot.background = element_blank()) + theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8)) +
  #   theme(text = element_text(size = 10),
  #         # axis.ticks.x = element_blank(),
  #         # axis.ticks.y = element_blank(),
  #         axis.ticks.x = element_line(colour =  "grey50"),
  #         axis.ticks.y = element_line(colour =  "grey50"),
  #         axis.line = element_line(colour =  "grey50" ),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.border = element_blank(),
  #         panel.background = element_blank())

#### Weight #####

  #Set up 
  
  #Set dummy hba1cmonth at 6 months
  final.wt <- final.wt %>% dplyr::mutate(hba1cmonth=6)
  
  #Predict HbA1c outcome uncalibrated
  final.wt <- final.wt %>% mutate(drug=drugclass,
                                  drugclass="DPP4")
  final.wt$DPP4.pred.lm.uncal <- predict(m1,final.wt)
  final.wt <- final.wt %>% mutate(drugclass="SGLT2")
  final.wt$SGLT2.pred.lm.uncal <- predict(m1,final.wt)
  final.wt <- final.wt %>% 
    mutate(hba1c_diff.uncal = SGLT2.pred.lm.uncal-DPP4.pred.lm.uncal,
           bestdrug.uncal=ifelse(hba1c_diff.uncal<=0,"SGLT2","DPP4"),
           drugclass=drug)

  #Recalibrate HbA1c outcome to AURUM
  final.wt <- final.wt %>% mutate(DPP4.pred.lm.recal=ifelse(ethnicity == "White",DPP4.pred.lm.uncal+ctfm$intercept[1],
                                                            ifelse(ethnicity == "South Asian",DPP4.pred.lm.uncal+ctfm$intercept[3],
                                                                   ifelse(ethnicity == "Black",DPP4.pred.lm.uncal+ctfm$intercept[5],
                                                                          ifelse(ethnicity == "mixed.other",DPP4.pred.lm.uncal+ctfm$intercept[7],NA
                                                                          )))),
                                  SGLT2.pred.lm.recal=ifelse(ethnicity == "White",SGLT2.pred.lm.uncal+ctfm$intercept[2],SGLT2.pred.lm.uncal),
                                  hba1c_diff.recal = SGLT2.pred.lm.recal-DPP4.pred.lm.recal,
                                  bestdrug.recal=ifelse(hba1c_diff.recal<=0,"SGLT2","DPP4"))
  
#Generate other variables for modelling
  
  #Define hba1c.breaks based on recalibrated HbA1c outcome
  final.wt <- final.wt %>% mutate(hba1c.breaks = cut(hba1c_diff.recal, breaks=c(min(hba1c_diff.recal,na.rm=T)-0.01,-5,-3,0,3,max(hba1c_diff.recal,na.rm=T)+0.01)))
  describe(final.wt$hba1c.breaks)
  
  #Set HbA1c diff to recalibrated HbA1c outcome
  final.wt <- final.wt %>% dplyr::rename(hba1c_diff="hba1c_diff.recal",
                                         bestdrug="bestdrug.recal")
  
  #Define datadist for modelling
  ddist <- datadist(final.wt); options(datadist='ddist') 
  
#### Predicted weight change from GOLD model (after HbA1c recalibration)  

#Predict weight change on each drug using GOLD model
  final.wt$drug <- "DPP4"
  final.wt$predwtonDPP4 <- predict(mwt, newdata=final.wt, se.fit=F)
  final.wt$drug <- "SGLT2"
  final.wt$predwtonSGLT2 <- predict(mwt, newdata=final.wt, se.fit=F)
  final.wt$drug <- final.wt$drugclass
  final.wt$predwtdiffsglt2dpp4 <- final.wt$predwtonSGLT2-final.wt$predwtonDPP4
  describe(final.wt$predwtdiffsglt2dpp4)
  hist(final.wt$predwtdiffsglt2dpp4)
  
#Plot mean and 95% interval of prediction by ethnicity and HbA1c defined subgroup
  white <- final.wt %>% filter(ethnicity=="White")
  asian <- final.wt %>% filter(ethnicity=="South Asian")
  black <- final.wt %>% filter(ethnicity=="Black")
  mixed.other <- final.wt %>% filter(ethnicity=="mixed.other")
  
  L <- list(white,asian,black,mixed.other)
  
  names(L) <- c("White",
                "Asian",
                "Black",
                "Mixed or Other")
  
  wt.res <- list()
  wt.plot <- list()
  
  for(i in 1:4) 
  { 
    #Subgroups by predicted treatment difference
    data  <- L[[i]] %>% dplyr::filter(hba1c_diff<= -5) 
    df <- expand.grid(mean(data$predwtonSGLT2),quantile(data$predwtonSGLT2,c(0.025)),quantile(data$predwtonSGLT2,c(0.975)),
                      mean(data$predwtonDPP4),quantile(data$predwtonDPP4,c(0.025)),quantile(data$predwtonDPP4,c(0.975)))
    SGLT2.5 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
    
    data  <- L[[i]] %>% dplyr::filter(hba1c_diff<= -3 & hba1c_diff > -5) 
    df <- expand.grid(mean(data$predwtonSGLT2),quantile(data$predwtonSGLT2,c(0.025)),quantile(data$predwtonSGLT2,c(0.975)),
                      mean(data$predwtonDPP4),quantile(data$predwtonDPP4,c(0.025)),quantile(data$predwtonDPP4,c(0.975)))
    SGLT2.3 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
    
    data  <- L[[i]] %>% dplyr::filter(hba1c_diff> -3 & hba1c_diff <= 0) 
    df <- expand.grid(mean(data$predwtonSGLT2),quantile(data$predwtonSGLT2,c(0.025)),quantile(data$predwtonSGLT2,c(0.975)),
                      mean(data$predwtonDPP4),quantile(data$predwtonDPP4,c(0.025)),quantile(data$predwtonDPP4,c(0.975)))
    SGLT2.0 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
    
    data  <- L[[i]] %>% dplyr::filter(hba1c_diff> 0 & hba1c_diff < 3) 
    df <- expand.grid(mean(data$predwtonSGLT2),quantile(data$predwtonSGLT2,c(0.025)),quantile(data$predwtonSGLT2,c(0.975)),
                      mean(data$predwtonDPP4),quantile(data$predwtonDPP4,c(0.025)),quantile(data$predwtonDPP4,c(0.975)))
    DPP4.0 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
    
    data  <- L[[i]] %>% dplyr::filter(hba1c_diff>= 3) 
    df <- expand.grid(mean(data$predwtonSGLT2),quantile(data$predwtonSGLT2,c(0.025)),quantile(data$predwtonSGLT2,c(0.975)),
                      mean(data$predwtonDPP4),quantile(data$predwtonDPP4,c(0.025)),quantile(data$predwtonDPP4,c(0.975)))
    DPP4.3 <- df %>% dplyr::rename(SGLT2.est=Var1,SGLT2.lci=Var2,SGLT2.uci=Var3,DPP4.est=Var4,DPP4.lci=Var5,DPP4.uci=Var6) 
    
    final.wt.obs <- rbind(SGLT2.5,SGLT2.3,SGLT2.0,DPP4.0,DPP4.3)
    df.lab <- c("SGLT2.5","SGLT2.3","SGLT2.0","DPP4.0","DPP4.3")
    final.wt.obs <- cbind(df.lab,final.wt.obs)
    final.wt.obs.p <- final.wt.obs %>% 
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 1) %>% 
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 5) 
    
    #res table
    wt.res[[i]] <- final.wt.obs.p
    
    #plot
    coef = data.matrix(cbind(final.wt.obs.p[,2],final.wt.obs.p[,5]))
    cim = data.matrix(cbind(final.wt.obs.p[,3],final.wt.obs.p[,6]))
    cip = data.matrix(cbind(final.wt.obs.p[,4],final.wt.obs.p[,7]))
    wt.plot[[i]] <-
      forestplot(row_names,
                 mean = coef,
                 lower= cim,
                 upper = cip,,
                 hrzl_lines = gpar(col="#444444"),lineheight=unit(2,'cm'),
                 is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
                 xticks = tick.wt,
                 zero = 0,
                 #boxsize=0.1,
                 # graphwidth = unit(2,"inches"),
                 # lineheight = unit(0.7,"inches"),
                 ci.vertices=TRUE,
                 col=fpColors(box=c("#f1a340","#4118de"), lines=c("#f1a340","#4118de"), zero = "gray50"),
                 lty.ci = c(1,2),           ,
                 xlab="Average weight change (kg)",cex=1,
                 title = names(L[i]),
                 new_page = TRUE,
                 fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
                 boxsize = .2, # We set the box size to better visualize the type
                 #line.margin = .1, # We need to add this to avoid crowding
                 txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),summary=gpar(cex = 1), ticks  = gpar(cex = 1), title=gpar(cex=1.3))
                 #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
                 #legend = c("SGLT2i","DPP4i"),
                 #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
                 #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
                 #xlog = TRUE
      )
    
  }  
  
#Plot together
  p.white <- grid2grob(print(wt.plot[[1]]))
  p.asian <- grid2grob(print(wt.plot[[2]]))
  p.black <- grid2grob(print(wt.plot[[3]]))
  p.mixed.other <- grid2grob(print(wt.plot[[4]]))
  
  wt.eth.pred <- (wrap_elements(p.white) + wrap_elements(p.asian)) /
    (wrap_elements(p.black) + wrap_elements(p.mixed.other)) /
    shared_legend +
    plot_layout(height=c(10,10,1))
  
  wt.eth.pred
  
  grDevices::cairo_pdf(paste0(output_dir,"wt_eth_recal_predicted.pdf"),width=pdfwidth,height=pdfheight)
  wt.eth.pred
  dev.off()
  
  png(paste0(output_dir,"wt_eth_recal_predicted.png"),width=pngwidth,height=pngheight,res=pngres,restoreConsole=TRUE)
  wt.eth.pred
  dev.off()
  
  wt.res[[1]]$eth = "White"
  wt.res[[2]]$eth = "Asian"
  wt.res[[3]]$eth = "Black"
  wt.res[[4]]$eth = "Mixed or Other"
  wt.res.pred <- wt.res %>% bind_rows(wt.res) %>% filter(!is.na(df.lab))
  write.csv(wt.res.pred,file=paste0(output_dir,"wt_eth_recal_predicted.csv"))  
  
#### Observed weight, by ethnicity and HbA1c benefit subgroup
  
  #Harrell method
  wt.res.unadjusted <- final.wt %>% 
    dplyr::group_by(ethnicity,hba1c.breaks,drugclass) %>% 
    dplyr::summarise(n=length(wtchange),
                     wt.change = mean(wtchange),
                     lci = mean(wtchange) - (1.96*(sd(wtchange)/sqrt(length(wtchange)))),
                     uci = mean(wtchange) + (1.96*(sd(wtchange)/sqrt(length(wtchange))))) %>%
    select(ethnicity,hba1c.breaks,drugclass,n,wt.change,lci,uci)
  
  #Long to wide for plotting
  wt.res.unadjusted <- wt.res.unadjusted %>% pivot_wider(
    names_from = drugclass,
    values_from = c(n,wt.change,lci,uci)
  )
  
  #Plot observed and 95% CI by ethnicity and HbA1c defined subgroup
  
  white <- wt.res.unadjusted %>% filter(ethnicity=="White") 
  asian <- wt.res.unadjusted %>% filter(ethnicity=="South Asian")
  black <- wt.res.unadjusted %>% filter(ethnicity=="Black") 
  mixed <- wt.res.unadjusted %>% filter(ethnicity=="mixed.other") 
  
  L <- list(white,asian,black,mixed)
  
  names(L) <- c("White",
                "Asian",
                "Black",
                "Mixed or Other")
  
  wt.plot <- list()
  
  for(i in 1:4) 
  { 
    #Subgroups by predicted treatment difference
    plotdata <- L[[i]] %>% ungroup() %>% as.data.frame()  %>%
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 1) %>% 
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 5) 
    
    #plot
    coef = data.matrix(cbind(plotdata[,6],plotdata[,5]))
    cim = data.matrix(cbind(plotdata[,8],plotdata[,7]))
    cip = data.matrix(cbind(plotdata[,10],plotdata[,9]))
    wt.plot[[i]] <-
      forestplot(row_names,
                 mean = coef,
                 lower= cim,
                 upper = cip,
                 hrzl_lines = gpar(col="#444444"),lineheight=unit(2,'cm'),
                 is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
                 xticks = tick.wt,
                 zero = 0,
                 #boxsize=0.1,
                 # graphwidth = unit(2,"inches"),
                 # lineheight = unit(0.7,"inches"),
                 ci.vertices=TRUE,
                 col=fpColors(box=c("#f1a340","#4118de"), lines=c("#f1a340","#4118de"), zero = "gray50"),
                 lty.ci = c(1,2),           ,
                 xlab="Average weight change (kg)",cex=1,
                 title = names(L[i]),
                 new_page = TRUE,
                 fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
                 boxsize = .2, # We set the box size to better visualize the type
                 #line.margin = .1, # We need to add this to avoid crowding
                 txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),summary=gpar(cex = 1), ticks  = gpar(cex = 1), title=gpar(cex=1.3))
                 #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
                 #legend = c("SGLT2i","DPP4i"),
                 #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
                 #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
                 #xlog = TRUE
      )
    
  }  
  
  
  #Plot together
  p.wt.white <- grid2grob(print(wt.plot[[1]]))
  p.wt.asian <- grid2grob(print(wt.plot[[2]]))
  p.wt.black <- grid2grob(print(wt.plot[[3]]))
  p.wt.mixed.other <- grid2grob(print(wt.plot[[4]]))
  
  wt.eth.obs <- (wrap_elements(p.wt.white) + wrap_elements(p.wt.asian)) /
    (wrap_elements(p.wt.black) + wrap_elements(p.wt.mixed.other)) /
    shared_legend +
    plot_layout(height=c(10,10,1))
  
  wt.eth.obs
  
  grDevices::cairo_pdf(paste0(output_dir,"wt_eth_recal_observed_unadjusted.pdf"),width=pdfwidth,height=pdfheight)
  wt.eth.obs
  dev.off()
  
  png(paste0(output_dir,"wt_eth_recal_observed_unadjusted.png"),width=pngwidth,height=pngheight,res=pngres,restoreConsole=TRUE)
  wt.eth.obs
  dev.off()
  
  write.csv(wt.res.unadjusted,file=paste0(output_dir,"wt_eth_recal_observed_unadjusted.csv"))
  
  #rename table outputs to compare obs vs pred
  wt.res.obs <- wt.res.unadjusted
  
#### Adjusted
  
  #overall
  summary(lm(wtchange~drug+preweight,data=final.wt))
  
  #Mean prediction
  #m.wt.obs <- ols(wtchange~drug*rcs(hba1c_diff,3)*ethnicity + preweight + drugline + ncurrtx,data=final.wt,x=T,y=T)
  m.wt.obs <- ols(wtchange~drug*hba1c.breaks*ethnicity + preweight + drugline + ncurrtx,data=final.wt,x=T,y=T)
  
  #Predict weight on each drug de novo
  final.wt$drug <- "DPP4"
  final.wt$DPP4.wt <- predict(m.wt.obs, newdata=final.wt, se.fit=F)
  final.wt$drug <- "SGLT2"
  final.wt$SGLT2.wt <- predict(m.wt.obs, newdata=final.wt, se.fit=F)
  final.wt$sglt2.dpp4.wt.diff <- final.wt$SGLT2.wt-final.wt$DPP4.wt
  describe(final.wt$sglt2.dpp4.wt.diff)
  hist(final.wt$sglt2.dpp4.wt.diff)
  final.wt$drug <- final.wt$drugclass
  
  wt.res.overall <- final.wt %>% 
    group_by(ethnicity,hba1c.breaks) %>% 
    dplyr::summarise(
      DPP4.wt = mean(DPP4.wt),
      SGLT2.wt = mean(SGLT2.wt),
      sglt2.dpp4.wt.diff = mean(sglt2.dpp4.wt.diff)
    )
  
  #Bootstrap predicted means
  n=nrow(final.wt)
  
  wt.res <- list()
  
  for(b in 1:B){
    i = sample(x = 1:n, size = n, replace = TRUE) ## sample indices
    temp = final.wt[i,] ## temp data set
    temp_model =   ols(wtchange~drug*hba1c.breaks*ethnicity + preweight + drugline + ncurrtx,data=final.wt,x=T,y=T)
    temp$drug <- "DPP4"
    temp$DPP4.wt <-  predict(temp_model, newdata=temp, se.fit=F)
    temp$drug <- "SGLT2"
    temp$SGLT2.wt <- predict(temp_model, newdata=temp, se.fit=F)
    temp$sglt2.dpp4.wt.diff <- temp$SGLT2.wt-temp$DPP4.wt
    
    res <- temp %>% 
      group_by(ethnicity,hba1c.breaks) %>% 
      dplyr::summarise(
        DPP4.wt = mean(DPP4.wt),
        SGLT2.wt = mean(SGLT2.wt),
        sglt2.dpp4.wt.diff = mean(sglt2.dpp4.wt.diff)
      )
    
    #res table
    wt.res[[b]] <- res
  }
  
  #Derive the bootstrapped CIs for the mean predictions for each subgrup
  wt.res.ci <- bind_rows(wt.res, .id = "bs_run") %>% 
    group_by(ethnicity,hba1c.breaks) %>% 
    dplyr::summarise(
      #DPP4.wt = mean(DPP4.wt),
      DPP4.wt.lci = quantile(DPP4.wt,probs=0.025),
      DPP4.wt.uci = quantile(DPP4.wt,probs=0.975),
      #SGLT2.wt = mean(SGLT2.wt),
      SGLT2.wt.lci = quantile(SGLT2.wt,probs=0.025),
      SGLT2.wt.uci = quantile(SGLT2.wt,probs=0.975),
      #sglt2.dpp4.wt.diff = mean(sglt2.dpp4.wt.diff),
      sglt2.dpp4.wt.diff.lci = quantile(sglt2.dpp4.wt.diff,probs=0.025),
      sglt2.dpp4.wt.diff.uci = quantile(sglt2.dpp4.wt.diff,probs=0.975),
    )
  
  
  wt.res.adjusted <- left_join(wt.res.overall,wt.res.ci,by=c("ethnicity","hba1c.breaks"))
  tail(wt.res.adjusted)
  
  #Plot mean and 95% interval of prediction by ethnicity and HbA1c defined subgroup
  white <- wt.res.adjusted %>% filter(ethnicity=="White") 
  asian <- wt.res.adjusted %>% filter(ethnicity=="South Asian")
  black <- wt.res.adjusted %>% filter(ethnicity=="Black") 
  mixed <- wt.res.adjusted %>% filter(ethnicity=="mixed.other") 
  #mixed <- wt.res.adjusted %>% filter(ethnicity=="Mixed") 
  #other <- wt.res.adjusted %>% filter(ethnicity=="Other") 
  
  L <- list(white,asian,black,mixed)
  
  names(L) <- c("White",
                "Asian",
                "Black",
                "Mixed or Other")
  
  # L <- list(white,asian,black,mixed,other)
  # 
  # names(L) <- c("White",
  #               "Asian",
  #               "Black",
  #               "Mixed",
  #               "Other")
  
  wt.plot <- list()
  
  for(i in 1:4) 
  { 
    #Subgroups by predicted treatment difference
    plotdata <- L[[i]] %>% ungroup() %>% as.data.frame()  %>%
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 1) %>% 
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 5) 
    
    #plot
    coef = data.matrix(cbind(plotdata[,4],plotdata[,3]))
    cim = data.matrix(cbind(plotdata[,8],plotdata[,6]))
    cip = data.matrix(cbind(plotdata[,9],plotdata[,7]))
    wt.plot[[i]] <-
      forestplot(row_names,
                 mean = coef,
                 lower= cim,
                 upper = cip,
                 #legend = c("SGLT2i","DPP4i"),
                 hrzl_lines = gpar(col="#444444"),#,lineheight=unit(2,'cm'))
                 is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
                 xticks = tick.wt,
                 zero = 0,
                 #boxsize=0.1,
                 # graphwidth = unit(2,"inches"),
                 # lineheight = unit(0.7,"inches"),
                 ci.vertices=TRUE,
                 col=fpColors(box=c("#f1a340","#4118de"), lines=c("#f1a340","#4118de"), zero = "gray50"),
                 lty.ci = c(1,2),           ,
                 xlab="Average weight change (kg)",cex=1,
                 title = names(L[i]),
                 new_page = TRUE,
                 fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
                 boxsize = .1, # We set the box size to better visualize the type
                 #line.margin = .1, # We need to add this to avoid crowding
                 txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),summary=gpar(cex = 1), ticks  = gpar(cex = 1), title=gpar(cex=1.3))
                 #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
                 #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
                 #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
                #xlog = TRUE
     #xlog = TRUE
      )
    
  }  
  
  #Plot together
  p.white <- grid2grob(print(wt.plot[[1]]))
  p.asian <- grid2grob(print(wt.plot[[2]]))
  p.black <- grid2grob(print(wt.plot[[3]]))
  p.mixed.other <- grid2grob(print(wt.plot[[4]]))
  
  wt.eth.obs.adj <- (wrap_elements(p.white) + wrap_elements(p.asian)) /
    (wrap_elements(p.black) + wrap_elements(p.mixed.other)) /
    shared_legend +
    plot_layout(height=c(10,10,1))
  
  wt.eth.obs.adj
  
  #Save
  grDevices::cairo_pdf(paste0(output_dir,"wt_eth_recal_observed_adjusted.pdf"))
  wt.eth.obs.adj
  dev.off()

  png(paste0(output_dir,"wt_eth_recal_observed_adjusted.png"),width=pngwidth,height=pngheight,res=pngres,restoreConsole=TRUE)
  wt.eth.obs.adj
  dev.off()
  
  write.csv(wt.res.adjusted,file=paste0(output_dir,"wt_eth_recal_observed_adjusted.csv"))