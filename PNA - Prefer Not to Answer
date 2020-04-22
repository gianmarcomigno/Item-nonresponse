########################################################
#                                                      #
#             -3: "Prefer not to answer"               #
#                                                      #
########################################################
# Analysis on the items with the "Prefer not to answer" option (-3 in UK Biobank)

# import all data
alldata <- readRDS('alldata.rds', refhook = NULL)

# select only those with the PNA option
alldata1 <- alldata %>% select(-f.userId, -f.54.0.0, -f.2129.0.0) %>%
  mutate_all(funs(ifelse(.==-3,1,0)))

# plot cluster dendrogram on tetrachoric correlation matrix
alldata1_tetra_corr <- tetrachoric(alldata1)$rho
colnames(alldata1_tetra_corr) <- rownames(alldata1_tetra_corr)  <- c("...
....
....
....
....
....")
plot(hclust(dist(alldata1_tetra_corr)), cex=0.62, xlab=NA, 
     main = "Cluster Dendrogram on tetrachoric correlation matrix")

interactive_matrix(get_lower_tri(alldata1_tetra_corr), 
                   scale= 'Tetrachoric correlation', main = 'Prefer Not to Answer tetrachoric correlation')


########################################################
#                    BEFORE PRUNING                    #
########################################################
# EFA Alldata - extract off the general latent factor showing the hidden missingness substructure 
fa1 <- fa(r= alldata1_tetra_corr,nfactors=1,fm="ols", n.obs=nrow(alldata1),rotate="biquartimin", residuals = TRUE) # Oblique

# heatmap residuals 1 factor fa():
interactive_matrix(fa1$residual, scale= 'Residuals intensity', main = 'Heatmap residuals 1 factor FA - PNA')

# dendrogram cluster residuals' distance - cut at height 0.500 for pruning
#rownames(fa1$residual) <- colnames(fa1$residual) <- substr(lab_vec[match(colnames(fa1$residual), paste0("f.",lab_vec$FieldID,".0.0")),]$Field,0, 25)
plot(hclust(dist(fa1$residual)), cex=0.62, xlab=NA, main = "Cluster Dendrogram on residuals from FA with 1 factor")
abline(h=.5, col= "red", lty = 4)
plot(cut(as.dendrogram(hclust(dist(1-fa1$residual))), h=.5)$upper, main="Upper tree of cut at h=0.5")
plot(cut(as.dendrogram(hclust(dist(1-fa1$residual))), h=.5)$lower[[2]], main="Second branch of lower tree with cut at h=.5")


########################################################
#                     PRUNING                          #
########################################################
# pruning 1: sum # of PNA within each branch
a <- alldata1
PNA_pruned1 <- as.data.frame(
  cbind(a$f.20116.0.0 + a$f.1239.0.0,
        a$f.21000.0.0,
        a$f.670.0.0 + a$f.6139.0.0 + a$f.699.0.0,
        a$f.6142.0.0, a$f.6138.0.0,
        a$f.864.0.0 + a$f.884.0.0 + a$f.904.0.0,
        a$f.1031.0.0, a$f.6160.0.0,
        a$f.1050.0.0 + a$f.1060.0.0 + a$f.1070.0.0, 
        a$f.1080.0.0 + a$f.1090.0.0,
        a$f.1100.0.0, a$f.1110.0.0,
        a$f.1160.0.0 + a$f.1170.0.0 + a$f.1190.0.0 + a$f.1200.0.0 + a$f.1220.0.0,
        a$f.1180.0.0, a$f.1210.0.0,
        a$f.1289.0.0 + a$f.1299.0.0 + a$f.1309.0.0 + a$f.1319.0.0,
        a$f.1329.0.0 + a$f.1339.0.0 + a$f.1349.0.0 + a$f.1359.0.0 + a$f.1369.0.0 + a$f.1379.0.0 + a$f.1389.0.0 + a$f.1418.0.0 + 
          a$f.1428.0.0 + a$f.1438.0.0 + a$f.1458.0.0 + a$f.1478.0.0 + a$f.1488.0.0 + a$f.1498.0.0 + a$f.1518.0.0 + 
          a$f.1528.0.0 + a$f.1558.0.0 + a$f.1548.0.0,
        a$f.6144.0.0, a$f.1538.0.0,
        a$f.1647.0.0 + a$f.1707.0.0 + a$f.1747.0.0 + a$f.1767.0.0 + a$f.1677.0.0 + a$f.1687.0.0 + a$f.1697.0.0,
        a$f.1737.0.0 + a$f.1757.0.0,
        a$f.1717.0.0 + a$f.1727.0.0,
        a$f.1920.0.0 + a$f.1930.0.0 + a$f.1940.0.0 + a$f.1950.0.0 + a$f.1960.0.0 + a$f.1970.0.0 + a$f.1980.0.0 + a$f.1990.0.0 +
          a$f.2000.0.0 + a$f.2010.0.0 + a$f.2020.0.0 + a$f.2030.0.0 + a$f.2040.0.0,
        a$f.2050.0.0 + a$f.2060.0.0 + a$f.2070.0.0 + a$f.2080.0.0,
        a$f.2090.0.0 + a$f.2100.0.0,
        a$f.2110.0.0, a$f.6149.0.0,
        a$f.2178.0.0 + a$f.2267.0.0 + a$f.2306.0.0 + a$f.2316.0.0 + a$f.2335.0.0 + a$f.2345.0.0 + a$f.2443.0.0 + a$f.2453.0.0 +
          a$f.2463.0.0 + a$f.2473.0.0 + a$f.2492.0.0 + a$f.6154.0.0 + a$f.2247.0.0 + a$f.2237.0.0,
        a$f.2188.0.0, a$f.6146.0.0,
        a$f.2207.0.0 + a$f.2227.0.0,
        a$f.2277.0.0, a$f.2296.0.0, a$f.6145.0.0, a$f.6159.0.0,
        a$f.6150.0.0 + a$f.6152.0.0,
        a$f.6155.0.0 + a$f.6179.0.0))
        
# colnames and labels
column_names <- c("branch_31","f.21000.0.0","branch_3","f.6142.0.0","f.6138.0.0",
                  "branch_4", "f.1031.0.0", "f.6160.0.0", "branch_9", "branch_6",
                  "f.1100.0.0", "f.1110.0.0", "branch_33", "f.1180.0.0", "f.1210.0.0",
                  "branch_34", "branch_35", "f.6144.0.0", "f.1538.0.0", "branch_30",
                  "branch_32", "branch_29", "branch_14", "branch_13", "branch_16",
                  "f.2110.0.0", "f.6149.0.0", "branch_27", "f.2188.0.0", "f.6146.0.0",
                  "branch_22", "f.2277.0.0", "f.2296.0.0", "f.6145.0.0", "f.6159.0.0",
                  "branch_23", "branch_21")
labels <- c("Smoking status","Ethnic background","Type of accommodation and lenght at current address","Current employment status","Qualifications",
            "Physical activity", "Frequency of friend/family visits", "Leisure/social activities", "Time spent outdoors or watching TV", "Time spent using PC or driving",
            "Drive faster than motorway speed limit", "Length of mobile phone use", "Sleep quality and duration", "Morning/evening person", "Snoring",
            "Fruit and vegetable intake", "Diet, food and alcohol intake", "Never eat eggs, dairy, wheat, sugar", "Major dietary changes in the last 5 years", "Person-specific information: adopted, handedness, hair color",
            "Childhood sunburn occasions and facial ageing", "Skin color and ease of skin tanning", "Depression/neuroticism feelings", "Frequency od depressed mood/unenthusiasm/tenseness", "Seen a psychiatrist or GP for nerves, anxiety, tension or depression",
            "Able to confide", "Mouth/teeth dental problems", "Overall health status", "Long-standing illness, disability or infirmity", "Attendance/disability/mobility allowance",
            "Wears glasses and eye problems", "Frequency of solarium/sunlamp use", "Falls in the last year", "Illness, injury, bereavement, stress in last 2 years", "Pain type(s) experienced in last month",
            "Hearth and lungs problems", "Vitamin and dietary supplements")
colnames(PNA_pruned1) <- column_names

# second step pruning (PNA only): split in 4 bins: 0, 1, 2: [2,n-1], 3: n PNAs for 4+ items
for(i in c("branch_33", "branch_34", "branch_35", "branch_30", "branch_14", "branch_13", "branch_27")){
  PNA_pruned1[,i] <- ifelse(PNA_pruned1[,i] == max(as.numeric(names(table(PNA_pruned1[,i])))), 'n', PNA_pruned1[,i])}
PNA_pruned2 <- as.data.frame(apply(PNA_pruned1, 2, function(x) {car::recode(x, "'0'= 0; '1'= 1; 'n'= 3; else= 2")}))


########################################################
#                  FACTOR ANALYSIS                     #
########################################################
# Splitting in training set and testing set (80:20) for EFA and CFA (before and after modification indeces), respectively 
set.seed(42) # Meaning of life: 42
#PNA_pruned2_train_ind <- sample(seq_len(nrow(PNA_pruned2)), size = floor(0.80 * nrow(PNA_pruned2)))
#PNA_pruned2_train <- PNA_pruned2[PNA_pruned2_train_ind, ] # EFA
#PNA_pruned2_test <- PNA_pruned2[-PNA_pruned2_train_ind, ] # CFA

# Convert to numeric for EFA
PNA_pruned2  <- as.data.frame(apply(PNA_pruned2, 2, as.numeric))

# Order the factors for CFA
PNA_pruned2  <- as.data.frame(apply(PNA_pruned2, 2, as.factor))
PNA_pruned2 <- plyr::catcolwise( function(v) ordered(v, levels = c(0, 1, 2, 3)))(PNA_pruned2)

# CFA after modification index on all dataset (100%) FINAL
cfa_PNA <- '
        g	=~  branch_31 + f.21000.0.0 + branch_3 + f.6142.0.0 + f.6138.0.0 + 
              branch_4 +  f.1031.0.0 +  f.6160.0.0 +  branch_9 +  branch_6 + 
              f.1100.0.0 +  f.1110.0.0 +  branch_33 +  f.1180.0.0 +  f.1210.0.0 + 
              branch_34 +  branch_35 +  f.6144.0.0 +  f.1538.0.0 +  branch_30 + 
              branch_32 +  branch_29 +  branch_14 +  branch_13 +  branch_16 + 
              f.2110.0.0 +  f.6149.0.0 +  branch_27 +  f.2188.0.0 +  f.6146.0.0 + 
              branch_22 +  f.2277.0.0 +  f.2296.0.0 +  f.6145.0.0 +  f.6159.0.0 + 
              branch_23 +  branch_21
        
        # Social status: Type of accommodation and lenght at current address, Current employment status, Qualifications, Frequency of friend/family visits 
        f1	=~  branch_3 + f.6142.0.0 + f.6138.0.0 + f.1031.0.0
        
        # Overall health
        f2	=~  f.6149.0.0 +  branch_27 +  f.2188.0.0 +  f.6146.0.0 + branch_22 + 
                f.2296.0.0 + f.6159.0.0 + branch_23 +  branch_21
        
        # Ethnicity
        f3	=~  f.21000.0.0 + branch_30 + branch_29 
        
        # Depression/Neuroticism
        f4	=~  branch_14 +  branch_13 +  branch_16 + f.2110.0.0 
        
        # General factor orthogonal to the group factors
        g  ~~  0*f1
        g  ~~  0*f2
        g  ~~  0*f3
        g  ~~  0*f4
        
        # Modification indeces
        branch_6 ~~  f.1100.0.0
        f.1100.0.0 ~~  f.1110.0.0
'
fit_cfa_PNA <- cfa(cfa_PNA, 
                   data = PNA_pruned2, sample.nobs = nrow(PNA_pruned2),
                   estimator="WLSMV", std.lv = TRUE,
                   ordered = colnames(PNA_pruned2))
saveRDS(fit_cfa_PNA, file = "fit_cfa_PNA.rds")
inspect(fit_cfa_PNA,what="std")$lambda
summary(fit_cfa_PNA,fit.measures=TRUE)
fitMeasures(fit_cfa_PNA)

# extraction factor scores: EBA method
factscores_PNA <- lavPredict(fit_cfa_PNA)
#saveRDS(factscores_PNA, file = "factscores_PNA.rds")
factscores_PNA <- as.data.frame(readRDS(file = "factscores_PNA.rds"))
factors_name <- c("general", "Social status", "Overall health", "Ethnicity", "Depression")
colnames(factscores_PNA) <- factors_name
