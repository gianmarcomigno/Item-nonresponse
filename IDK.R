########################################################
#                                                      #
#                 -1: "I don't know"                   #
#                                                      #
########################################################
# Analysis on the items with the "I don't know" option.

# import all data
alldata <- readRDS('alldata.rds', refhook = NULL)

# select only those with the IDK option
not1 <- c(names(which(colSums(alldata==-1)==0)))
alldata2 <- alldata[,!colnames(alldata) %in% c(not1,"f.1647.0.0"), with = F]  # remove f.1647.0.0
alldata2 <- as.data.frame(ifelse(alldata2==-1,1,0))

# plot cluster dendrogram on tetrachoric correlation matrix
alldata2_corr <- tetrachoric(alldata2)$rho
colnames(alldata2_corr) <- rownames(alldata2_corr)  <- c("...
....
....   # selected items/variables in the format of f.XXXX.0.0
....
....
....")

heatmap(alldata2_corr)
plot(eigen(alldata2_corr)$values[1:15], type = "l") # plot eigenvalues
abline(h=1, col = "red", lty = 2) 


########################################################
#                    BEFORE PRUNING                    #
########################################################
# EFA Alldata - extract off the general latent factor showing the hidden missingness substructure 
fa1 <- fa(r= alldata2,nfactors=1,fm="ols", n.obs=nrow(alldata2),rotate="biquartimin", residuals = TRUE, cor='tet')

# heatmap residuals 1 factor fa():
interactive_matrix(fa1$residual, scale= 'Residuals intensity', main = 'Heatmap residuals after 1 factor FA - IDK')

# dendrogram cluster residuals' distance - cut at height 0.775 for pruning
#rownames(fa2$residual) <- colnames(fa2$residual) <- lab_vec[match(colnames(fa2$residual), paste0("f.",lab_vec$FieldID,".0.0")),]$Field 
plot(hclust(dist(1-fa1$residual)), cex=0.62, xlab=NA, 
     main = "Cluster Dendrogram on (1-residuals) from FA 1 factor")
abline(h=.775, col= "red", lty = 4)
plot(cut(as.dendrogram(hclust(dist(1-fa2$residual))), h=.775)$upper, main="Upper tree of cut at h=1")
labels(cut(as.dendrogram(hclust(dist(1-fa2$residual)), h=7)$lower[[25]], main="branch of lower tree with cut at h=1"))


########################################################
#                     PRUNING                          #
########################################################
# pruning 1: sum # of IDK within each branch
a <- as.data.frame(alldata2)
IDK_pruned_sum <- 
            cbind(a$f.21000.0.0, a$f.6139.0.0, 
                  a$f.864.0.0 + a$f.884.0.0 + a$f.904.0.0, 
                  a$f.1050.0.0 + a$f.1060.0.0, 
                  a$f.1080.0.0 + a$f.1090.0.0,
                  a$f.1100.0.0,  
                  a$f.2277.0.0, a$f.1110.0.0, a$f.1031.0.0, a$f.1170.0.0, a$f.1220.0.0, a$f.1548.0.0, a$f.1180.0.0, 
                  a$f.1210.0.0,
                  a$f.1289.0.0 + a$f.1299.0.0 + a$f.1309.0.0 + a$f.1319.0.0,
                  a$f.1329.0.0 + a$f.1339.0.0 + a$f.1349.0.0 + a$f.1359.0.0 + a$f.1369.0.0 + a$f.1379.0.0 + a$f.1389.0.0,
                  a$f.1418.0.0 + a$f.1428.0.0, 
                  a$f.1438.0.0 + a$f.1458.0.0 + a$f.1488.0.0 + a$f.1498.0.0 + a$f.1528.0.0,
                  a$f.1070.0.0 + a$f.1160.0.0, 
                  a$f.1747.0.0, a$f.699.0.0, a$f.1677.0.0, 
                  a$f.1687.0.0 + a$f.1697.0.0,
                  a$f.1737.0.0, a$f.1757.0.0, a$f.1717.0.0, a$f.1727.0.0, 
                  a$f.1920.0.0 + a$f.1930.0.0 + a$f.1960.0.0, 
                  a$f.1950.0.0, a$f.1980.0.0, a$f.2000.0.0, a$f.2040.0.0, a$f.1940.0.0, 
                  a$f.1970.0.0 + a$f.1990.0.0 + a$f.2010.0.0, 
                  a$f.2020.0.0, a$f.2030.0.0, 
                  a$f.2050.0.0 + a$f.2060.0.0 + a$f.2070.0.0, 
                  a$f.2080.0.0, 
                  a$f.2090.0.0 + a$f.2100.0.0, 
                  a$f.2110.0.0, a$f.2178.0.0, a$f.2306.0.0, a$f.2316.0.0, a$f.2335.0.0, a$f.2345.0.0, a$f.2443.0.0, 
                  a$f.2453.0.0, a$f.2463.0.0, a$f.2492.0.0, a$f.6154.0.0, a$f.2247.0.0, a$f.2188.0.0, a$f.2473.0.0, 
                  a$f.6146.0.0, a$f.1767.0.0, a$f.2267.0.0)
IDK_pruned_sum <- as.data.frame(IDK_pruned_sum)

# colnames and labels
col_names_IDK <- c("f.21000.0.0", "f.6139.0.0", "branch_34", "branch_37", "branch_39", "f.1100.0.0", "f.2277.0.0", 
                   "f.1110.0.0","f.1031.0.0","f.1170.0.0","f.1220.0.0","f.1548.0.0","f.1180.0.0","f.1210.0.0",
                   "branch_38", "branch_33", "branch_48", "branch_41", "branch_40",  
                   "f.1747.0.0","f.699.0.0","f.1677.0.0","branch_44", 
                   "f.1737.0.0", "f.1757.0.0", "f.1717.0.0", "f.1727.0.0", "branch_18",
                   "f.1950.0.0","f.1980.0.0","f.2000.0.0","f.2040.0.0","f.1940.0.0",
                   "branch_22", "f.2020.0.0","f.2030.0.0","branch_17","f.2080.0.0","branch_12","f.2110.0.0",
                   "f.2178.0.0", 'f.2306.0.0', 'f.2316.0.0', 'f.2335.0.0','f.2345.0.0','f.2443.0.0',
                   'f.2453.0.0', 'f.2463.0.0',"f.2492.0.0","f.6154.0.0","f.2247.0.0",
                   "f.2188.0.0", "f.2473.0.0", "f.6146.0.0", "f.1767.0.0","f.2267.0.0")
labels_IDK <- c("Ethnic background", "Gas or solid-fuel cooking/heating", "Physical Activities", "Time spent outdoors", "Time spent driving/using PC", 
                "Drive faster than motorway speed limit", "Frequency of solarium/sunlamp use", "Length of mobile phone use", "Frequency of friend/family visits", "Getting up in morning",
                "Daytime dozing / sleeping (narcolepsy)", "Variation in diet", "Morning/evening person (chronotype)", "Snoring", "Fruit/Vegetable intake",
                "Fish and meat intake", "Milk and spread type", "Cereal/Bread and Coffee/Tea intake", "Time spent watching TV and sleep duration",
                "Hair colour (natural, before greying)", "Length of time at current address", "Breastfed as a baby",
                "Comparative body/height size at age 10", "Childhood sunburn occasions", "Facial ageing", "Skin colour",
                "Ease of skin tanning", "Fed-up feelings and mood swings", "Sensitivity / hurt feelings", "Worrier / anxious feelings", 
                "Worry too long after embarrassment", "Risk taking", "Irritability", "Nervous and tense feelings", "Loneliness, isolation",
                "Guilty feelings", "Frequency of unenthusiasm/tenseness feelings", "Frequency of tiredness / lethargy in last 2 weeks", "Seen a doctor or psychiatrist fro nerves, anxiety or depression",
                "Able to confide", "Overall health rating", "Weight change compared with 1 year ago", "Wheeze or whistling in the chest in last year", "Chest pain or discomfort",
                "Ever had bowel cancer screening", "Diabetes diagnosed by doctor", "Cancer diagnosed by doctor", "Fractured/broken bones in last 5 years", "Taking other prescription medications",
                "Medication for pain relief, constipation, heartburn", "Hearing difficulty/problems", "Long-standing illness, disability or infirmity", 
                "Other serious medical condition/disability diagnosed by doctor", 
                "Attendance/disability/mobility allowance", "Adopted as a child", "Use of sun/uv protection")
colnames(IDK_pruned_sum) <- col_names_IDK


########################################################
#                  FACTOR ANALYSIS                     #
########################################################
# splitting in training set and testing set (80:20) for EFA and CFA (before and after modification indeces), respectively 
set.seed(42) # Meaning of life: 42
# similar to FA in PNA, look there for more details

# order the factors for CFA
IDK_pruned_sum  <- as.data.frame(apply(IDK_pruned_sum, 2, ordered))

# CFA after modification index on all dataset (100%) FINAL
cfa_IDK <- '
      g	=~   f.21000.0.0 + f.6139.0.0 + branch_34 +  branch_37  + branch_39  + f.1100.0.0 + f.2277.0.0 +
              f.1110.0.0 + f.1031.0.0 + f.1170.0.0 + f.1220.0.0 + f.1548.0.0 + f.1180.0.0 + f.1210.0.0 +
              branch_38 + branch_33 + branch_48 + branch_41 + branch_40 + f.1747.0.0 + f.699.0.0  +
              f.1677.0.0 + branch_44 + f.1737.0.0 + f.1757.0.0 + f.1717.0.0 + f.1727.0.0 + branch_18 +  
              f.1950.0.0 + f.1980.0.0 + f.2000.0.0 + f.2040.0.0 + f.1940.0.0 + branch_22  + f.2020.0.0 +
              f.2030.0.0 + branch_17 + f.2080.0.0 + branch_12 + f.2110.0.0 + f.2178.0.0 + f.2306.0.0 +
              f.2316.0.0 + f.2335.0.0 + f.2345.0.0 + f.2443.0.0 + f.2453.0.0 + f.2463.0.0 + f.2492.0.0 +
              f.6154.0.0 + f.2247.0.0 + f.2188.0.0 + f.2473.0.0 + f.6146.0.0 + f.1767.0.0 + f.2267.0.0 
      
      # Depression/Neuroticism (fed-up/ hurt/ anxious/ guilty/ tense feelings)
      f1=~   f.1950.0.0 + f.1980.0.0 + f.2000.0.0 + f.2040.0.0 + f.1940.0.0 + branch_22  + f.2020.0.0 +
              f.2030.0.0 + branch_17 + f.2080.0.0 + branch_18

      # Lifestyle (Physical Activities, Time spent, Frequency of visits, Food intake)
      f2=~   branch_34 +  branch_37  + branch_39 + f.1110.0.0 + f.1031.0.0 + f.1170.0.0 + 
              branch_38 + branch_33 + branch_41 + branch_40 + f.1737.0.0  

      # Overall health (Chest pain, disability, medical condition, Hearing difficulties)
      f3=~  f.2316.0.0 + f.2335.0.0 + f.2247.0.0 + f.2188.0.0 + f.2473.0.0

      # General factor orthogonal to the group factors
      g  ~~  0*f1
      g  ~~  0*f2
      g  ~~  0*f3
      
      # modification index
      branch_17 ~~  f.2080.0.0
      branch_44 ~~  f.1737.0.0
'

fit_cfa_IDK <- cfa(cfa_IDK, 
                   data = IDK_pruned_sum, sample.nobs = nrow(IDK_pruned_sum),
                   estimator="WLSMV", std.lv = TRUE,
                   ordered = colnames(IDK_pruned_sum))
saveRDS(fit_cfa_IDK, file = "fit_cfa_IDK.rds")
inspect(fit_cfa_IDK,what="std")$lambda
summary(fit_cfa_IDK,fit.measures=TRUE)
fitMeasures(fit_cfa_IDK)

# extraction factor scores: EBA method
factscores_IDK <- lavPredict(fit_cfa_IDK)
#saveRDS(factscores_IDK, file = "factscores_IDK.rds")
factscores_IDK <- as.data.frame(readRDS(file = "factscores_IDK.rds"))


# gwas on idk if breastfed as a baby
idk_breastfed <- alldata %>% select("f.userId","f.1677.0.0") %>%
  mutate(f.1677.0.0 = if_else(f.1677.0.0 == -1,1,0))
write.table(idk_breastfed, file='IDK_breastfed.tsv', quote=FALSE, sep='\t', row.names=FALSE)
