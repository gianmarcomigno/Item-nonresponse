#In the first part of the code, the data are imported and both the participants and the relevant variables of interest have been subset.

# Import libriaries
libraries <- c("data.table", "reshape2", "ggplot2", "plotly", "spatstat", 
               "stringr", "plyr", "dplyr", "psych", "GPArotation", "lavaan", 
               "car", "tidyverse", "tidylog", "RColorBrewer", "dichromat")
lapply(libraries, require, character.only = TRUE)

# Reading Phenotype file
d <- fread("neale_lab_parsed.tsv", header=T)
colnames(d) <- paste0("f.",gsub("x","",gsub("_",".",colnames(d))))

# Read in the sample and subsetting, excluding those who gave no permission
sample <- read.table("ukb31063.gwas_samples.both_sexes.txt", header = TRUE)
s1 <- read.csv('w31063_20181016.csv')
s2 <- setdiff(sample$s,s1$X1142430)

# Subsetting Variables - only those asked to every participant. NB only f.X.0.0 are included.
keep <- c("f.userId","f.54.0.0","f.20116.0.0","f.21000.0.0","f.670.0.0","f.6139.0.0","f.699.0.0","f.6142.0.0","f.6138.0.0","f.864.0.0",
          "f.884.0.0","f.904.0.0","f.1031.0.0","f.6160.0.0","f.1050.0.0","f.1060.0.0","f.1070.0.0","f.1080.0.0","f.1090.0.0",
          "f.1100.0.0","f.1110.0.0","f.1160.0.0","f.1170.0.0","f.1180.0.0","f.1190.0.0","f.1200.0.0","f.1210.0.0","f.1220.0.0",
          "f.1239.0.0","f.1289.0.0","f.1299.0.0","f.1309.0.0","f.1319.0.0","f.1329.0.0","f.1339.0.0","f.1349.0.0","f.1359.0.0",
          "f.1369.0.0","f.1379.0.0","f.1389.0.0","f.6144.0.0","f.1418.0.0","f.1428.0.0","f.1438.0.0","f.1458.0.0","f.1478.0.0",
          "f.1488.0.0","f.1498.0.0","f.1518.0.0","f.1528.0.0","f.1538.0.0","f.1548.0.0","f.1558.0.0","f.1647.0.0","f.1677.0.0",
          "f.1687.0.0","f.1697.0.0","f.1707.0.0","f.1717.0.0","f.1727.0.0","f.1737.0.0","f.1747.0.0","f.1757.0.0","f.1767.0.0",
          "f.1920.0.0","f.1930.0.0","f.1940.0.0","f.1950.0.0","f.1960.0.0","f.1970.0.0","f.1980.0.0","f.1990.0.0","f.2000.0.0",
          "f.2010.0.0","f.2020.0.0","f.2030.0.0","f.2040.0.0","f.2050.0.0","f.2060.0.0","f.2070.0.0","f.2080.0.0","f.2090.0.0",
          "f.2100.0.0","f.2110.0.0","f.6149.0.0","f.2129.0.0","f.2178.0.0","f.2188.0.0","f.6146.0.0","f.2207.0.0","f.2227.0.0",
          "f.2267.0.0","f.2277.0.0","f.2296.0.0","f.2306.0.0","f.2316.0.0","f.6145.0.0","f.6159.0.0","f.2335.0.0","f.2345.0.0",
          "f.6150.0.0","f.6152.0.0","f.2443.0.0","f.2453.0.0","f.2463.0.0","f.2473.0.0","f.2492.0.0","f.6154.0.0","f.6155.0.0",
          "f.2247.0.0","f.6179.0.0","f.2237.0.0", "f.53.0.0")

# Subsetting partecipants, excluding pilot-study predecessors
d1 <- d[d$f.userId %in% s2, keep, with=FALSE]
d1 <- d1[d1$f.53.0.0 > "2007-01-01",]
d1 <- d1[,!"f.53.0.0"]

# Remove partecipants with too many NAs
var_miss <- data.frame(xvar = names(colSums(is.na(d1))),
                       yvar = colSums(is.na(d1)))
var_miss_num <- var_miss[var_miss$yvar==max(colSums(is.na(d1))),]$xvar

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols, with = F])
  return(data[completeVec, ])}

alldata<- completeFun(d1, var_miss_num)

########################################################
#                    CHECKPOINT 1                      #
########################################################
saveRDS(alldata, file = 'alldata.rds')
#write.csv(alldata, file = "alldata.csv")
alldata <- readRDS('alldata.rds', refhook = NULL) 
alldata <- alldata[,-1]

# Adding labels to variables, to include the UserID = "f.1.0.0"
names(alldata)[1] <- "f.1.0.0"
infos <- rbind(infos, data.frame("FieldID" = 1, "Field" = "UserID",stringsAsFactors = TRUE))

# Read in Uk Biobank info and adding labels
info <- read.csv('Data_Dictionary_Showcase.csv')
infos <- info[c(grep("UK Biobank assessment centre",info$Field), grep("Touchscreen", info$Path)),c("FieldID","Field")]

lab_vec <- infos[paste0("f.",infos$FieldID,".0.0") %in% colnames(alldata),c("Field","FieldID")]
labvec <- lab_vec[match(names(alldata), paste0("f.",lab_vec$FieldID,".0.0")),]

colnames(alldata) <- lab_vec[match(colnames(alldata), paste0("f.",lab_vec$FieldID,".0.0")),]$Field

# Variables' information
# Type: C = Categorical, O = Ordinal, N = Numeric Integer. Multiple Instances: 1 = yes, 0 = no. Progressive.code: code in the survey
vars_info <- read.table(file = 'vars_info.tsv', sep = '\t', header = TRUE)




# ---------------------------------


# import full phenotypes dataset
d <- fread("ukb31063.phenotypes.20191008.csv", header=T)
colnames(d) <- paste0("f.",gsub("x","",gsub("_",".",colnames(d))))
d1 <- d %>% select(grep("20400", colnames(d)), grep("f.eid", colnames(d)),  # 110003 has only NAs, 110005 not found
                   grep("12188", colnames(d)),                             
                   grep("110001", colnames(d)), grep("110005", colnames(d)),
                   grep("2000", colnames(d))) 
saveRDS(d1, file = "recontact_variables.rds")

# read in all files
alldata <- readRDS("alldata.rds")
covariates <- fread("ukb31063.neale_gwas_covariates.both_sexes.csv", skip = 1)
covariates <- covariates[covariates$s %in% alldata$f.userId,] # subset
covariates <- covariates[match(alldata$f.userId, covariates$s),] # sort for userID

pna_eba <- data.frame(readRDS(file = "factscores_PNA.rds"))['g']
idk_eba <- data.frame(readRDS(file = "factscores_IDK.rds"))['g']

# COVARIATES: dataframe with all the  covariates in the model
covariate_fin <- as.data.frame(
  cbind(covariates, # female, age, age2, age*female, age2*female, PC1--20
        alldata$f.54.0.0, alldata$f.2178.0.0, alldata$f.6138.0.0, # f.54: assessment center, f.6138: education, f.2138: health
        scale(idk_eba), # scaled idk factorscores phenotypes
        scale(pna_eba)) # scaled pna factorscores phenotypes
)
colnames(covariate_fin)[c(1,27:31)] <- c("userID","assessm_center","health","eduYears", "scaled_idk", "scaled_pna")

# recode
covariate_fin$health <- car::recode(covariate_fin$health,   # overall health (1218 idk, 113 pna)
                                    "-1=NA ; -3=NA ")
covariate_fin$isFemale <- ifelse(covariate_fin$isFemale == T, 1,0) #1 is female
covariate_fin$eduYears <- car::recode(covariate_fin$eduYears,  # impute f.6138 as eduyears in paper ea gwas (james lee et al.)
                                      "1=20 ; 2=15 ; 3=13 ; 4=12 ;
                                      5=19 ; 6=17 ; -7=6 ; -3=NA ")

#--------------------------------------------------------
#--------------------------------------------------------
recontact_variables <- readRDS("recontact_variables.rds")
covariate_fin <- readRDS("covariate_fin.rds")

# Field 110001-1.0: Invitation to complete online 24-hour recall dietary questionnaire - acceptance (1st istance)
acceptance <- recontact_variables %>% select("f.eid", "f.110001.1.0") %>% 
  mutate(f.110001.1.0 = if_else(f.110001.1.0 == 0, 1, 0)) %>% #recode (0 non reponse, 2 completed)
  filter(f.eid %in% covariate_fin$userID) %>% #subset
  left_join(covariate_fin, by = c("f.eid" = "userID")) %>%
  rename(nonresponded_diet_questionn_1visit = f.110001.1.0) %>% 
  na.omit 

#logistic model
model1 <- glm(nonresponded_diet_questionn_1visit ~ isFemale + age + age_squared + age_isFemale + age_squared_isFemale + PC1 + PC2 + PC3 + PC4 +
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +  
                health +
                eduYears +
                scaled_idk +   #0.0003
                scaled_pna,     #0.0010
              data = acceptance, family = "binomial")

summary(model1)

blr_rsq_mckelvey_zavoina(model1)

# Field 110001-1.0: Invitation to complete online 24-hour recall dietary questionnaire - acceptance (1st istance) - 
# NON-RESPONDED ALL VISITS (1) vs COMPLETED ALL VISITS (0)
none_vs_all <- recontact_variables %>% select("f.eid", "f.110001.1.0", "f.110001.2.0", "f.110001.3.0", "f.110001.4.0") %>% 
  mutate(none = case_when(
    f.110001.1.0+f.110001.2.0+f.110001.3.0+f.110001.4.0 == 8 ~ 0, 
    f.110001.1.0+f.110001.2.0+f.110001.3.0+f.110001.4.0 == 0 ~ 1, 
    TRUE ~ NA_real_)) %>% #recode (1 non-reponded to all visits, 0 completed all visits, NA: completed some visits)
  filter(f.eid %in% covariate_fin$userID) %>% #subset
  left_join(covariate_fin, by = c("f.eid" = "userID")) %>%
  rename(nonresponded_all_visits = none) %>% 
  select(-c("f.110001.1.0", "f.110001.2.0", "f.110001.3.0", "f.110001.4.0")) %>%
  na.omit()

 model0 <- glm(nonresponded_all_visits ~ isFemale + age + age_squared + age_isFemale + age_squared_isFemale + PC1 + PC2 + PC3 + PC4 +
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +  
                health +
                eduYears +
               scaled_idk +   #0.0012
               scaled_pna,     #0.0027
              data = none_vs_all, family = "binomial")
 
summary(model0)

 
blr_rsq_mckelvey_zavoina(model1) #https://www.slideshare.net/MarcoDAlessandro11/pseudor-quadro for different methods

###########################
# ROC curve logistic regression
#####################"######
install.packages("pROC")
library(pROC)

mydata <- acceptance
model1 <- glm(nonresponded_diet_questionn_1visit ~ isFemale + age + age_squared + age_isFemale + age_squared_isFemale + PC1 + PC2 + PC3 + PC4 +
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + 
                eduYears + assessm_center + health + scaled_pna + scaled_idk, 
              data = acceptance, family = "binomial")
acceptance$prob <- predict(model1,type=c("response"))
plot(roc(nonresponded_diet_questionn_1visit ~ prob, data = acceptance))
