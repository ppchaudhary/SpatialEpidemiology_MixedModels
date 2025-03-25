library(doParallel)
library(foreach)
library(lme4)
library(MuMIn)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
require(ggplot2)
require(reshape2)
library(raster)
library(sf)
library(spdep)
library(zipcodeR)
library(RVAideMemoire)

options(max.print=10000)

### Variables ####
# all or peds?
demographic<-"peds"
# Test just a subset of the data if running on regular laptop:
  test_subset<-

#### Read ####
# Read Exposure data frame
# ****** replace with 2010-2016 version when ready and ADD peds version too
# Load Exposure Data Frame
exposure_file <- ifelse(demographic == "all", "features_matrix_10_16_all.csv", "features_matrix_10_16_peds.csv")
Exposure <- read_csv(exposure_file, col_types = cols(ZipCode = col_character(), .default = col_double()), show_col_types = FALSE)
# Read the disease_df and denominator
disease_df <- read_csv("under18set4_df_5.csv")
denominator <- read_csv("under18_denominator.csv")

#Limit to diseases with >330 zip codes with data 
disease_df <- disease_df[, colSums(disease_df != 0) >= 330]


##### Clean Exposure, denominator dataframes
dim(Exposure) # 16515   544
nrow(denominator) # 16309
nrow(disease_df) # 16309
sum(is.na(disease_df)) # No missing data
sum(is.na(denominator)) # No missing data
sum(is.na(Exposure)) # Good amount of missing data. Hm what's going on.
# OK Looks like there are 67 zip codes which are uniformly missing across all chemicals
missing_index<-which(is.na(Exposure$Acetonitrile))
# Removing these yields zero problems
sum(is.na(Exposure[-missing_index,]))
# So just remove those
Exposure<-Exposure[-missing_index,]
# Same number of rows in disease_df and denominator, but more in exposure dataframe. 
# Double check that everything is numeric except zip code
colnames(Exposure)[sapply(Exposure,class)!="numeric"]
# Check duplicate zipcodes
sum(duplicated(Exposure$ZipCode)) # There are 0 duplicate zipcodes
#remove duplicates
#Exposure <- Exposure[!duplicated(Exposure), ] #now all the duplicates are gone

sum(duplicated(denominator$ZipCode)) # There are no duplicate zipcodes
# Make sure zip codes with leading zeros are 5 digits
Exposure$ZipCode<- stringr::str_pad(Exposure$ZipCode, width=5, pad = "0")
denominator$ZipCode<- stringr::str_pad(denominator$ZipCode, width=5, pad = "0")

# REMOVE ZIP CODES (79+123 = 202) where there is no overlap of data
sum(!(Exposure$ZipCode%in%denominator$ZipCode)) # 235
sum(!(denominator$ZipCode%in%Exposure$ZipCode)) # 116
# Common ZipCodes
common_zipcodes <- Reduce(intersect, list(denominator$ZipCode, disease_df$ZipCode, Exposure$ZipCode))
Exposure <- Exposure[Exposure$ZipCode %in% common_zipcodes, ]
denominator <- denominator[denominator$ZipCode %in% common_zipcodes, ]
disease_df <- disease_df[disease_df$ZipCode %in% common_zipcodes, ]
# Check to make sure it worked
sum(!(Exposure$ZipCode%in%denominator$ZipCode)) # 0
sum(!(denominator$ZipCode%in%Exposure$ZipCode)) # 0
# Merge for strictly overlapping
#This one from deepseek
Data <- merge(Exposure, denominator, by = "ZipCode") %>% arrange(ZipCode)
#Data<-merge(Exposure,denominator, by.x="ZipCode",by.y="ZipCode", all.x=FALSE, all.y=FALSE)
Data <- Data[ , -2]


ZipCodes<-Data$ZipCode
Data<-Data%>% dplyr::select(-ZipCode)
disease_df<-disease_df %>% dplyr::select(-ZipCode)

if(test_subset==TRUE) {
  Data<-Data[1:10000,(ncol(Data)-50):ncol(Data)]
  ZipCodes<-ZipCodes[1:10000]
  disease_df<-disease_df[1:10000,1:10]
}


#### Spatial Weights ####
# Matrix for Moran's I #
lonlat <- cbind(Data$Longitude, Data$Latitude)
pts <- SpatialPoints(lonlat)
crdref <- CRS('+proj=longlat +datum=NAD83')
pts <- SpatialPoints(lonlat, proj4string=crdref)
nb_dist<-dnearneigh(x=pts,d1=0,d2=50*1.609344)
lw <-nb2listwdist(neighbours=nb_dist, x=pts,type="idw", style="W", alpha = 1,zero.policy = TRUE)

#### Spatial Clusters ####
mdist  <- geosphere::distm(cbind(Data$Longitude,Data$Latitude))
hc <- hclust(as.dist(mdist ), method="complete")

Data $cluster0 <- as.factor(cutree(hc, k=round(nrow(Data )/81)))
Data $cluster1 <- as.factor(cutree(hc, k=round(nrow(Data )/27)))
Data $cluster2 <- as.factor(cutree(hc, k=round(nrow(Data )/9)))
Data $cluster3 <- as.factor(cutree(hc, k=round(nrow(Data )/3)))


# Select appropriate age range
if(demographic=="all"){
  Ages_list<-c("Ages_0_5","Ages_5_9", "Ages_10_14","Ages_15_19",
               "Ages_20_24" , "Ages_25_29", "Ages_30_34","Ages_35_39",
               "Ages_40_44" ,"Ages_45_49", "Ages_50_54" ,"Ages_55_59","Ages_60_64",
               "Ages_65_69" , "Ages_70_74" ,
               "Ages_75_79", "Ages_80_84" ,"Ages_85_plus")
}

if(demographic=="peds"){
  Ages_list<-c("Ages_0_5","Ages_5_9", "Ages_10_14","Ages_15_19")
}

## Remove last age category to prevent too much multicollinearity
if(demographic=="peds"){
  Data<-Data[,!colnames(Data)=="Ages_15_19"]
  Ages_list_short<-Ages_list[!Ages_list=="Ages_15_19"]
}
if(demographic=="all"){
  Data<-Data[,!colnames(Data)=="Ages_85_plus"]
  Ages_list_short<-Ages_list[!Ages_list=="Ages_85_plus"]
}

#### SCALE ####
# Select numeric columns to index for scaling but select out the denominator 
numeric_cols <- sapply(Data, is.numeric)
total_index<-which(colnames(Data)=="Total")
numeric_cols[total_index]<-FALSE
# Scale numeric columns
#This one from deepseek
scaled_matrix <- as.data.frame(lapply(Data[, numeric_cols], function(x) as.numeric(scale(x))))
#scaled_matrix <- as.data.frame(lapply(Data[, numeric_cols], scale))
# Set column names back to original
colnames(scaled_matrix) <- names(Data)[numeric_cols]
# Combine scaled numeric columns with non-numeric columns
Data <- cbind(Data[!numeric_cols], scaled_matrix)

# Use foreach for parallel processing
Mixed_Effect<-function(ICD){

  # This can be replaced by the code below, which preserves column names
  #Data[, sapply(Data, is.numeric)] <- lapply(Data[, sapply(Data, is.numeric)], scale)
  
numerator<-disease_df[[ICD]]
  
  MM_loop <-NULL
  MM_loop <-foreach (i = c(1:ncol(Data)),.combine="rbind",.errorhandling = 'remove')%dopar%{

    # For troubleshooting, print each iteration of for loop
    # print(paste("Iteration:", i, ":", colnames(Data)[i]))
    
    if(colnames(Data)[i] %in% c("Total","cluster0","cluster1","cluster2","cluster3")){next}
    
    tryCatch({
      
      Covariates<-c(
        Ages_list,
        "Latitude",
        # "Deprivation",
        "Population_Density",
        "cluster0","cluster1","cluster2","cluster3",
        "Total") 
      
      if(colnames(Data)[i] %in% Covariates) {
        skips<-which(colnames(Data) %in% Covariates)
      } else {
        skips<-which(colnames(Data) %in% c(colnames(Data)[i],Covariates) )}
      
      Data_sample=Data[,skips]
      Data_sample = Data_sample%>%relocate(as.name(colnames(Data)[i]))

      try ( {
        
       #  Check unique values of columns for troubleshooting
       #  sapply(lapply(Data_sample, unique), length)
        
        model  <- glmer.nb(numerator ~.
                        - Total
                        - cluster0 - cluster1 - cluster2 - cluster3 
                        + offset(log(Total)) 
                        + (1|cluster0/cluster1/cluster2/cluster3)
                        ,
                        data=Data_sample,
                        nAGQ = 0,
                        control=glmerControl(optimizer="bobyqa",
                                             optCtrl=list(maxfun=1e9)))
        
        
      }, silent=TRUE)
      
      summo<-summary(model)
      
      if(summo$coefficients[2,4]<0.1){
        I_ac<-moran.test(residuals(model,type="pearson"), lw, zero.policy = TRUE)
        MI<-I_ac$estimate[1]
        MI_p<-I_ac$p.value
      }else{
        MI<-NA
        MI_p<-NA
      }
      
      overdisp_output<-capture.output(RVAideMemoire::overdisp.glmer(model))
      ratio <- sub(".*ratio: ([0-9.]+).*", "\\1", overdisp_output)
      
    }, warning = function(w) {
      message(paste("Warning message on iteration", i, ":", conditionMessage(w)))
    })
    
    # Disease Name, Beta Coefficient, p value, AIC, Moran's I, Moran's I p value, overdispersion ratio
    c(colnames(Data)[i],summo$coefficients[2,1],summo$coefficients[2,4],summo$AICtab[1],MI,MI_p,ratio)
    
  }

  MM_loop<-as.data.frame(MM_loop)
  MM_loop <- MM_loop %>% mutate_at(vars(2:7), as.numeric)
  MM_loop<- MM_loop %>% dplyr::mutate(V2=exp(V2))
  colnames(MM_loop)<-c("Variable","Adjusted RR", "P-value","AIC","Moran's I", "Moran's I p-value","Overdispersion")
  MM_loop$Variable<-gsub("_", " ", MM_loop$Variable)
  MM_loop$Variable<-gsub("PM25", "PM2.5", MM_loop$Variable)
  MM_loop <- MM_loop[order(MM_loop$`P-value`),]
  rownames(MM_loop)<-NULL
  
  # Write results to CSV with ICD-specific filename
  write.csv(MM_loop, file = file.path("results", paste0(ICD, "_results.csv")), row.names = FALSE)
  
  return(MM_loop)
}


start_time <- Sys.time()
cl<-makeCluster(future::availableCores()-2)
registerDoParallel(cores=cl)
Loop_list<-lapply(as.list(colnames(disease_df)),FUN=Mixed_Effect)
names(Loop_list)<-colnames(disease_df)
stopCluster(cl)
# Time elapsed
end_time <- Sys.time()
end_time - start_time

saveRDS(Loop_list, "ME_results_youth.rds")
#Loop_list_df %>% write_csv("ME_results_youth.csv")

#from deep seek
# Save each entry in Loop_list as a separate CSV
#for (disease_category in names(Loop_list)) {
  # Extract results for this disease category
  #results_df <- Loop_list[[disease_category]]
  
  # Clean up filename by removing special characters
  #clean_name <- gsub("[^[:alnum:]]", "_", disease_category)
  
  # Write to CSV with demographic prefix
  #write.csv(
    #results_df,
    #file = paste0(demographic, "_", clean_name, "_results.csv"),
    #row.names = FALSE
  #)
#}
