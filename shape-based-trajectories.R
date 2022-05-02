library(ggplot2)
library(reshape2) # long to wide
library(tidyverse) # separate

### load long format data
ldata$F1_sex <- as.factor(ldata$F1_sex)

##### 1. Baseline centering ----------------------------------------------------

#### 1.1 ZHFA ------------------------------------------------------------------
zhfa001 <- by(ldata, ldata$id, function(x) 
  within(x, zhfa_norm001 <- zhfa - (1.002*zhfa[1])))
ldata <- do.call(rbind, zhfa001) # takes some time (like 2 min)
head(ldata)

# long to wide
# make subset
idata_zhfa001 <- cbind.data.frame(ldata %>% select("F_id", "zhfa_norm001"))
head(idata_zhfa001)
idata_zhfa001 <- idata_zhfa001 %>% separate(F_id, c("followup", "id"), sep = "_")
idata_zhfa001$followup <-  paste(idata_zhfa001$followup, "zhfa_norm001", sep = "_")
head(idata_zhfa001)
w_idata_zhfa001 <- dcast(idata_zhfa001, id ~ followup, value.var = "zhfa_norm001")
head(w_idata_zhfa001)

# add to main wide dataset
wdata <- merge(wdata, w_idata_zhfa001, by="id", all=TRUE)
colnames(wdata)

#### 1.2 ZWFA ------------------------------------------------------------------
# move everything to start at ~0.01
zwfa001 <- by(ldata, ldata$id, function(x) 
  within(x, zwfa_norm001 <- zwfa - (1.002*zwfa[1])))
ldata <- do.call(rbind, zwfa001) # takes some time (like 2 min)
head(ldata)

# long to wide
# make subset
idata_zwfa001 <- cbind.data.frame(ldata %>% select("F_id", "zwfa_norm001"))
head(idata_zwfa001)
idata_zwfa001 <- idata_zwfa001 %>% separate(F_id, c("followup", "id"), sep = "_")
idata_zwfa001$followup <-  paste(idata_zwfa001$followup, "zwfa_norm001", sep = "_")
head(idata_zwfa001)
w_idata_zwfa001 <- dcast(idata_zwfa001, id ~ followup, value.var = "zwfa_norm001")
w_idata_zwfa001 <- w_idata_zwfa001[, 1:7]
head(w_idata_zwfa001)

# add to main wide dataset
wdata <- merge(wdata, w_idata_zwfa001, by="id", all=TRUE)
colnames(wdata)

#### 1.3 ZBFA ------------------------------------------------------------------
# move everything to start at ~0.01
zbfa001 <- by(ldata, ldata$id, function(x) 
  within(x, zbfa_norm001 <- zbfa - (1.002*zbfa[1])))
ldata <- do.call(rbind, zbfa001) # takes some time (like 2 min)
head(ldata)

# long to wide
# make subset
idata_zbfa001 <- cbind.data.frame(ldata %>% select("F_id", "zbfa_norm001"))
head(idata_zbfa001)
idata_zbfa001 <- idata_zbfa001 %>% separate(F_id, c("followup", "id"), sep = "_")
idata_zbfa001$followup <-  paste(idata_zbfa001$followup, "zbfa_norm001", sep = "_")
head(idata_zbfa001)
w_idata_zbfa001 <- dcast(idata_zbfa001, id ~ followup, value.var = "zbfa_norm001")
head(w_idata_zbfa001)

# add to main wide dataset
wdata <- merge(wdata, w_idata_zbfa001, by="id", all=TRUE)
colnames(wdata)

## 2. Number of clusters -------------------------------------------------------
library(NbClust)
library(sjlabelled)

## 2.1 Create datasets for each z-score ----------------------------------------
# ** zhfa_ig_n001
my_data_zhfa_n001<-cbind.data.frame(idata %>% select("id", contains("zlen_ig_norm001"), contains("zhfa_ig_norm001")))
head(my_data_zhfa_n001)

# transpose
t_my_data_zhfa_n001 <- transpose(my_data_zhfa_n001)
# get row and colnames in order
colnames(t_my_data_zhfa_n001) <- rownames(my_data_zhfa_n001)
rownames(t_my_data_zhfa_n001) <- colnames(my_data_zhfa_n001)

# ** zwfa_ig_n001
my_data_zwfa_n001<-cbind.data.frame(idata %>% select("id", contains("zlen_ig_norm001"), contains("zwfa_ig_norm001")))
head(my_data_zwfa_n001)

# transpose (to create dissimilarity matrix) - nbClust
t_my_data_zwfa_n001 <- transpose(my_data_zwfa_n001)
# get row and colnames in order
colnames(t_my_data_zwfa_n001) <- rownames(my_data_zwfa_n001)
rownames(t_my_data_zwfa_n001) <- colnames(my_data_zwfa_n001)

# ** zbfa_ga_n001
my_data_zbfa_n001<-cbind.data.frame(idata %>% select("id", contains("zlen_ga_norm001"), contains("zbfa_ga_norm001")))
head(my_data_zbfa_n001)

# transpose (to create dissimilarity matrix) - nbClust
t_my_data_zbfa_n001 <- transpose(my_data_zbfa_n001)
# get row and colnames in order
colnames(t_my_data_zbfa_n001) <- rownames(my_data_zbfa_n001)
rownames(t_my_data_zbfa_n001) <- colnames(my_data_zbfa_n001)

## 2.2 Create dissimilarity matrix ---------------------------------------------
########## using frechet distance
my_data_zhfa_n001_frech_diss <- TSDatabaseDistances(t_my_data_zhfa_n001, distance = "frechet") #2 hr
saveRDS(my_data_zhfa_n001_frech_diss, file="~/zhfa/Frechet/NbClust/my_data_zhfa_n001_frech_diss.RData")

my_data_zwfa_n001_frech_diss <- TSDatabaseDistances(t_my_data_zwfa_n001, distance = "frechet") #2 hr
saveRDS(my_data_zwfa_n001_frech_diss, file="~/zwfa/Frechet/NbClust/my_data_zwfa_n001_frech_diss.RData")

my_data_zbfa_n001_frech_diss <- TSDatabaseDistances(t_my_data_zbfa_n001, distance = "frechet") #2 hr
saveRDS(my_data_zbfa_n001_frech_diss, file="~/zbfa/Frechet/NbClust/my_data_zbfa_n001_frech_diss.RData")

## 2.3 NbClust -----------------------------------------------------------------
## remove labels and run NbClust
my_data_zhfa_n001_frech_diss <- readRDS("~/zhfa/Frechet/NbClust/diss_matrix/my_data_zhfa_n001_frech_diss.RData")
zhfa_n001 <-  remove_all_labels(my_data_zhfa_n001) # 11 proposed 2 as the best number of clusters 20190830
zhfa_n001_frech_nbclust <- NbClust(zhfa_n001, diss = my_data_zhfa_n001_frech_diss, distance = NULL, min.nc = 2, max.nc = 15, method = "kmeans")
saveRDS(zhfa_n001_frech_nbclust, file="~/zhfa/Frechet/NbClust/zhfa_n001_frech_nbclust.RData")

my_data_zwfa_n001_frech_diss <- readRDS("~/zwfa/Frechet/NbClust/diss_matrix/my_data_zwfa_n001_frech_diss.RData")
zwfa_n001 <-  remove_all_labels(my_data_zwfa_n001) # 11 proposed 2 as the best number of clusters 20190830
zwfa_n001_frech_nbclust <- NbClust(zwfa_n001, diss = my_data_zwfa_n001_frech_diss, distance = NULL, min.nc = 2, max.nc = 15, method = "kmeans")
saveRDS(zwfa_n001_frech_nbclust, file="~/zwfa/Frechet/NbClust/zwfa_n001_frech_nbclust.RData")

my_data_zbfa_n001_frech_diss <- readRDS("~/zbfa/Frechet/NbClust/diss_matrix/my_data_zbfa_n001_frech_diss.RData")
zbfa_n001 <-  remove_all_labels(my_data_zbfa_n001) # 11 proposed 2 as the best number of clusters 20190830
zbfa_n001_frech_nbclust <- NbClust(zbfa_n001, diss = my_data_zbfa_n001_frech_diss, distance = NULL, min.nc = 2, max.nc = 15, method = "kmeans")
saveRDS(zbfa_n001_frech_nbclust, file="~/zbfa/Frechet/NbClust/zbfa_n001_frech_nbclust.RData")

## 3. Cluster trajectories  ----------------------------------------------------
library(tidyverse)
library(dplyr)
library(kmlShape)
library(reshape2) #covert wide to long. melt

## set timepoints
cldstimes <- c(0L, 3L, 12L, 24L, 49L, 80L, 130L)

## 3.1 Run clustering ----------------------------------------------------------
## create clds
## ** ZHFA
myCldsGzhfa3_ig <- cldsWide(idata_zhfa_ig, times = cldstimes)
kmlShape(myCldsGzhfa3_ig, nbClusters =  3, toPlot = "none")
saveRDS(myCldsGzhfa3_ig, file="~/zhfa/Frechet/age_corrected/myCldsGzhfa3_ig.RData")

## ** zwfa
myCldsGzwfa3_ig <- cldsWide(idata_zwfa_ig, times = cldstimes)
kmlShape(myCldsGzwfa3_ig, nbClusters =  3, toPlot = "none")
saveRDS(myCldsGzwfa3_ig, file="~/zwfa/Frechet/age_corrected/myCldsGzwfa3_ig.RData")

## ** zbfa
myCldsGzbfa3_ga <- cldsWide(idata_zbfa_ga, times = cldstimes)
kmlShape(myCldsGzbfa3_ga, nbClusters =  3, toPlot = "none")
saveRDS(myCldsGzbfa3_ga, file="~/zbfa/Frechet/age_corrected/myCldsGzbfa3_ga.RData")

## 3.2 Add results to main dataset ---------------------------------------------
# ZHFA
myCldsGzhfa3_ig <- readRDS(file="~/zhfa/Frechet/age_corrected/myCldsGzhfa3_ig.RData")
# new minidatabase for zhfa_ig no age
idata_zhfa_ig <- idata_zhfa_ig_age[, -c(2:8)]
head(idata_zhfa_ig)

## merge clusters with dataset
clus3_ig <- as.data.frame(myCldsGzhfa3_ig@id)
colnames(clus3_ig) <- c("id")
clus3_ig$zhfa_ig_3c_400 <- myCldsGzhfa3_ig@clusters
colnames(clus3_ig)
idata_zhfa_ig <- merge(idata_zhfa_ig, clus3_ig, by = "id")
head(idata_zhfa_ig)

# zwfa
myCldsGzwfa3_ig <- readRDS(file="~/zwfa/Frechet/age_corrected/myCldsGzwfa3_ig.RData")
# new minidatabase for zwfa_ig no age
idata_zwfa_ig <- idata_zwfa_ig_age[, -c(2:8)]
head(idata_zwfa_ig)

## merge clusters with dataset
clus3_ig <- as.data.frame(myCldsGzwfa3_ig@id)
colnames(clus3_ig) <- c("id")
clus3_ig$zwfa_ig_3c_400 <- myCldsGzwfa3_ig@clusters
colnames(clus3_ig)
idata_zwfa_ig <- merge(idata_zwfa_ig, clus3_ig, by = "id")
head(idata_zwfa_ig)

# ZHFA
myCldsGzbfa3_ga <- readRDS(file="~/zbfa/Frechet/age_corrected/myCldsGzbfa3_ga.RData")
# new minidatabase for zbfa_ga no age
idata_zbfa_ga <- idata_zbfa_ga_age[, -c(2:8)]
head(idata_zbfa_ga)

## merge clusters with dataset
clus3_ga <- as.data.frame(myCldsGzbfa3_ga@id)
colnames(clus3_ga) <- c("id")
clus3_ga$zbfa_ga_3c_400 <- myCldsGzbfa3_ga@clusters
colnames(clus3_ga)
idata_zbfa_ga <- merge(idata_zbfa_ga, clus3_ga, by = "id")
head(idata_zbfa_ga)

