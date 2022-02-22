# Farmland practices are driving bird populations decline across Europe.

R scripts and data for the following article: "Farmland practices are driving bird populations decline across Europe."

## R

The R scripts have been implemented on R version 3.4.4.

### Loading R packages

```{r setup, include=FALSE}

# Load packages

source("R_packages.R")

```


## Bird data

### Loading data

#### Relative abundance indices

```{r}

# Downlad Species national indices from Brlik et al. (2021)

df <- fread('https://zenodo.org/record/4590199/files/national_indices2017.csv?download=1')

# Pass from wide to long format

df <- melt(df, id.vars=c("species","euring_code","scheme","type"))
df <- dcast(df, species+euring_code+scheme+variable~type, fun.aggregate = sum)
df <- df[,c("euring_code","species","scheme","variable","index","se")]
names(df) <- c("Code","Species","CountryGroup","Year","Index","Index_SE")
df <- droplevels(na.omit(df))
df$Species <- as.factor(df$Species)
df$CountryGroup <- as.factor(df$CountryGroup)
df$Year <- as.numeric(as.character(df$Year))

```

#### Abundance data
```{r}

# Load data from the EU Bird Directive Reporting (see Supplementary material for more details)

abd <- setDT(read.table("Abundance_data_PECBMS.txt", header = T, sep="\t")) 


```

### Preparing data

#### Species and coutry names

```{r}
# Get species and coutry names

S <- levels(df$Species)
C <- levels(df$CountryGroup)

# Update species names

diff_name <- merge(data.frame(sp=levels(as.factor(abd$Species)),num=1),data.frame(sp=levels(as.factor(df$Species)),num2=2), by="sp",all=T)
diff_name

df$Species <- as.character(df$Species)
df$Species[df$Species=="Carduelis cannabina"] <- "Linaria cannabina"
df$Species[df$Species=="Carduelis chloris"] <- "Chloris chloris"
df$Species[df$Species=="Carduelis flammea"] <- "Acanthis flammea"
df$Species[df$Species=="Carduelis spinus"] <- "Spinus spinus"
df$Species[df$Species=="Corvus corone+cornix"] <- "Corvus corone"
df$Species[df$Species=="Delichon urbica"] <- "Delichon urbicum"
df$Species[df$Species=="Dendrocopos medius"] <- "Dendrocoptes medius"
df$Species[df$Species=="Dendrocopos minor"] <- "Dryobates minor"
df$Species[df$Species=="Hippolais pallida"] <- "Iduna pallida"
df$Species[df$Species=="Hirundo daurica"] <- "Cecropis daurica"
df$Species[df$Species=="Hirundo rupestris"] <- "Ptyonoprogne rupestris"
df$Species[df$Species=="Miliaria calandra"] <- "Emberiza calandra"
df$Species[df$Species=="Parus ater"] <- "Periparus ater"
df$Species[df$Species=="Parus caeruleus"] <- "Cyanistes caeruleus"
df$Species[df$Species=="Parus cristatus"] <- "Lophophanes cristatus"
df$Species[df$Species=="Parus montanus"] <- "Poecile montanus"
df$Species[df$Species=="Parus palustris"] <- "Poecile palustris"
df$Species[df$Species=="Saxicola torquata"] <- "Saxicola torquatus"
df$Species[df$Species=="Serinus citrinella"] <- "Carduelis citrinella"
df$Species[df$Species=="Tetrao tetrix"] <- "Lyrurus tetrix"
df$Species <- as.factor(df$Species)

```

#### Merge relative abundance indices and bird abundance
```{r}

# Generate abundance estimates (for a given year)

for (i in 1:nrow(abd)){
  abd[i,estimate := round((geoMean(c(Count_min,Count_max),na.rm=T))*2)]
}

# Create sub-estimates for Belgium and Germany regions
# assuming populations to be uniformely distributed across these countries
# frac corresponds to the fraction surface of the region

subd <- data.table(reg = c("Belgium-Brussels", "Belgium-Wallonia", "Germany East", "Germany West"), 
                   frac = c(161/30528, 16901/30528, 108333/357022, 248577/357022))
Country_subd <- abd[Country=="Belgium" | Country=="Germany"][rep(1:340,each=2)]
Country_subd[,Country := rep(subd[,reg],170)]

# estimating pop size from the percentage of total area of the country covered by the region

for (i in 1:4){
  Country_subd[Country==subd[,reg][i],Count_min := Count_min*subd[,frac][i]]
  Country_subd[Country==subd[,reg][i],Count_max := Count_max*subd[,frac][i]]
}
for (i in 1:680){
  Country_subd[i,estimate := round((geoMean(c(Count_min,Count_max),na.rm=T))*2)]
}

# adding Belgium and Germany regions estimates to the main abundance data.frame

abd <- rbind(abd,Country_subd)
setorder(abd,Species,Country)

# removing unrealised species-country combinations

abd <- abd[is.na(Count_min) == FALSE]

# Convertion [Relative => Absolute] abundance time series

for(i in 1:nrow(subd)){
  df$Index[which(df$CountryGroup==subd$reg[i])] <- df$Index[which(df$CountryGroup==subd$reg[i])]*subd$frac[i]
  df$Index_SE[which(df$CountryGroup==subd$reg[i])] <- df$Index_SE[which(df$CountryGroup==subd$reg[i])]*subd$frac[i]
}

df_belge <- data.table(droplevels(subset(df, CountryGroup %in% c("Belgium-Brussels","Belgium-Wallonia"))) %>%
                       group_by(Code, Species, Year) %>% summarize(Index=sum(Index),Index_SE=sum(Index_SE)))
                       
df_belge <- data.table(df_belge[,1:2],CountryGroup=rep("Belgium",nrow(df_belge)),df_belge[,3:5])

df_germany <- data.table(droplevels(subset(df, CountryGroup %in% c("Germany East","Germany West"))) %>%
                       group_by(Code, Species, Year) %>% summarize(Index=sum(Index),Index_SE=sum(Index_SE)))

df_germany <- data.table(df_germany[,1:2],CountryGroup=rep("Germany",nrow(df_germany)),df_germany[,3:5])

df <- rbind(droplevels(subset(df, !(CountryGroup %in% c("Belgium-Brussels","Belgium-Wallonia","Germany East","Germany West")))),df_belge, df_germany)

df[,start_year := min(Year), by=c("Species","CountryGroup")]
df[,end_year := max(Year), by=c("Species","CountryGroup")]

for (i in levels(df$Species)){

  for (j in levels(df$CountryGroup)){
  
    # defining the reference year for the computation of a weighing factor
    
    if (dim(df[Species==i & CountryGroup==j])[1]!=0){
    
      # if the last year of absolute abundance estimate falls within the abundance time span
      # the abs. abund. ref year matches that of the time series
      
      if (abd[Species==i & Country==j, Year_end] <= df[Species==i & CountryGroup==j, end_year][1] &
          abd[Species==i & Country==j, Year_end] >= df[Species==i & CountryGroup==j, start_year][1]){
        Yref <- abd[Species==i & Country==j, Year_end]
        popsize <- "ok"
      }
      
      # if the last year of absolute abundance estimate falls after the abundance time span
      # the abs. abund. ref year is the last of the time series
      
      if (abd[Species==i & Country==j, Year_end] > df[Species==i & CountryGroup==j, end_year][1]){
        Yref <- df[Species==i & CountryGroup==j, end_year][1]
        popsize <- "later"
      }
      
      # if the last year of absolute abundance estimate falls before the abundance time span
      # the abs. abund. ref year is the first of the time series
      
      if (abd[Species==i & Country==j, Year_end] < df[Species==i & CountryGroup==j, start_year][1]){
        Yref <- df[Species==i & CountryGroup==j, start_year][1]
        popsize <- "earlier"
      }
      
      # computation of a weighing factor for the ref year
      
      WF <- abd[Species==i & Country==j, estimate]/df[Species==i & CountryGroup==j & Year==Yref, Index]
      
      # convertion using the weighing factor
      
      df[Species==i & CountryGroup==j, Abd := round(Index*WF)]
      df[Species==i & CountryGroup==j, SE_Abd := round(Index_SE*WF)]
      df[Species==i & CountryGroup==j, ref_year := Yref]
      df[Species==i & CountryGroup==j, estimation := popsize]
    }
  }
}

# particular case of collapsed populations

df[Species=="Galerida cristata" & CountryGroup=="Czech Republic",c("Abd","SE_Abd") := 0]

# Merge index and abundance

df_pop <- droplevels(subset(df, Year %in% c(1980:2016)))
df_pop <- droplevels(subset(df_pop, Species %in% levels(droplevels(subset(df_pop, start_year<=1981))$Species)))
df_pop <- droplevels(df_pop[!which(df_pop$Species=="Passer domesticus" & df_pop$CountryGroup %in% levels(df_pop$CountryGroup)[23]),])
```

## Combining species indices and abundance at the European level: RTIM

### Importing data

See also details in `rtrim` documentation (Bogaart et al., 2016).

```{r}
# Data import & variables

S <- unique(df_pop$Species)
Code <- unique(df_pop$Code)
C <- unique(df_pop$CountryGroup)

# assign countries with a number

stratum_number <- data.frame(code_sp = C, site = 1:length(C))

# assign species with their EURING code

species_code <- df_pop[, .SD[1], by = .(Code)][, .(Code,Species)]
```

### RTRIM preprocessing
```{r}
# Set working directory to /output for intermediary results

setwd("output/")

# Create data files for each species cf. README_TRIM_SHELL.docx

for (i in S){
  subset <- df_pop[Species==i]
  subset <- droplevels(subset)
  c <- levels(subset[, CountryGroup]) 
  nyear <- max(subset$Year) - min(subset$Year)+1
  ncountry <- length(levels(subset[, CountryGroup]))
  transition <- data.frame(code_sp = rep(c,each=nyear),
                           Year = rep(c(min(subset$Year):max(subset$Year)),ncountry))
  transition <- merge(transition, stratum_number, by = "code_sp", all.x = TRUE)
  transition <- merge(transition, subset[, .(CountryGroup,Year,Abd)],
                      by.x = c("code_sp","Year"), by.y = c("CountryGroup", "Year"),
                      all = TRUE)
  dataset <- transition[, c("site","Year","Abd")]
  setnames(dataset, old = c("Year","Abd"), new = c("year","count"))
  
  write.csv2(dataset, file =paste0("BIRD_", species_code[Species==i,Code],
                                   "_0_counts.csv"), row.names = FALSE)
  arg <- data.frame(File = paste0("BIRD_",species_code[Species==i,Code],"_0"),
                    Base_year_first_year = 1, Base_year_last_year = 1,
                    Changepoints = "all", Serial_correlation = TRUE,
                    Overdispersion = TRUE, Presence_weights = FALSE,
                    Slope_from = min(subset$Year)-1 + floor(nyear/2),
                    Save_fitted_values = TRUE) # explicitation of parameters in shell_rtrim_strata.R
  write.csv2(arg, file = paste0("BIRD_",df_pop[Species==i,Code][1],
                                "_0_arg_input_stratum.csv"), row.names = FALSE)
}
```

### RTRIM functions

```{r}
# Load RTRIM functions

source("RTRIM_functions.R")

```

### Running RTRIM
```{r}
# Select the directory which contains the files with counts that need to be analysed. 
folder<-getwd()
setwd(folder)

# Select the files in the directory which contain arguments.

listSpeciesStratumCombinations <- dir(folder, pattern = "arg_input_stratum.csv")

# The number of files in the directory.                                                
numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)

# Make a table for an overview. 

overview <- makeOverview(listSpeciesStratumCombinations)

for(j in 1:numberSpeciesStratumCombinations){

# The file with arguments only contains the arguments to analyse the counts for one 
# particular combination of species and stratum. 
# The arguments are used when calling the function 'rtrim'.

  print(j)
  arguments <- read.csv2(listSpeciesStratumCombinations[j], header = TRUE, stringsAsFactors = FALSE)
  arguments$File <- gsub(arguments$File, pattern = "s:/SWAN_POP/MOO/TRIM_WERK/", replacement = "")

# The file with counts only contains the counts for one particular combination 
# of species and stratum. Weights are also in this file.

  counts <- read.csv2(paste(arguments$File, "_counts.csv", sep = ""), header = TRUE)
  overview$first_year[j] <- min(counts$year, na.rm = TRUE)
  overview$last_year[j] <- max(counts$year, na.rm = TRUE)
  
# The entire period for each combination of stratum and species is stored.
# Necessary to make one output file with indices and trends for all combinations 
# of species and strata.
  
# RUNNING TRIM
# Start with the most elaborate model, if necessary use a more simple model.                                                  
  result <- tryCatch(
    {
    
# First attempt, using the most elaborate model.

      if (arguments$Presence_weights == TRUE){
          trim(count ~ site + year, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
        } else{
          trim(count ~ site + year, data = counts,                      model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
          }
    }
    , error = warning)     
  
  if (class(result) == "trim") {
    save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
    overview$attempt_1[overview$ss_combinations == arguments$File] <- "success"
    overview$success[overview$ss_combinations == arguments$File] <- "yes"
  } else {
  
# The first attempt was unsuccessful a less elaborate model will be tried.

      overview$attempt_1[overview$ss_combinations == arguments$File] <- "error"
      overview$error_1[overview$ss_combinations == arguments$File] <- result
      result <- tryCatch(
        {
        
# Second attempt, using a less elaborate model.

        if (arguments$Presence_weights == TRUE){
          trim(count ~ site + year, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
        } else {
          trim(count ~ site + year, data = counts,                      model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
          }
      }
      , error = warning)
      if (class(result) == "trim") {
        save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
        overview$attempt_2[overview$ss_combinations == arguments$File] <- "success"
        overview$success[overview$ss_combinations == arguments$File] <- "yes"
      } else {
      
# Also the second attempt has failed, the model is further simplified. 
# No changepoints, serial correlation is taken into account.

          overview$attempt_2[overview$ss_combinations == arguments$File] <- "error"
          overview$error_2[overview$ss_combinations == arguments$File] <- result
          result <- tryCatch( 
            {
            
# Third attempt.

                if (arguments$Presence_weights == TRUE){
                  trim(count ~ site + year, data = counts, weights = "weights", model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
              } else {
                  trim(count ~ site + year, data = counts,                      model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
              }
            }
            , error = warning)
          if (class(result) == "trim") {
            save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
            overview$attempt_3[overview$ss_combinations == arguments$File] <- "success"
            overview$success[overview$ss_combinations == arguments$File] <- "yes"
          } else {
                overview$attempt_3[overview$ss_combinations == arguments$File] <- "error"
                overview$error_3[overview$ss_combinations == arguments$File] <- result
                
# The third attempt has also failed. 
# Now the most simple model is tried. No changepoints, no serial correlation.

              result <- tryCatch(
                {
# Last attempt.
                    if (arguments$Presence_weights == TRUE){
                      trim(count ~ site + year, data = counts, weights = "weights", model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
                  } else {
                      if (arguments$Presence_weights == FALSE){
                        trim(count ~ site + year, data = counts,                    model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
                      } 
                    }
                }
                , error = warning)
              if (class(result) == "trim") {
               save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
               overview$attempt_4[overview$ss_combinations == arguments$File] <- "success"
               overview$success[overview$ss_combinations == arguments$File] <- "yes"
              } else {
              
# When the analysis fails again a message is printed to the screen.
# Also in overview.csv this information is stored, together with the error message.

                  overview$attempt_4[overview$ss_combinations == arguments$File] <- "error"
                  overview$error_4[overview$ss_combinations == arguments$File] <- result
                  cat("Problem with this combination of species and stratum:", arguments$File, "\n")
        }
      }  
    } 
  }  
}

overview <- overview[order(overview$species_number, overview$stratum_number), ]
listsuccessfullAnalyses <- overview$ss_combinations[overview$success == "yes"]
listSpeciesStratumCombinations <- paste(listsuccessfullAnalyses, "_arg_input_stratum.csv", sep = "")

# Determine how many combinations of species and stratum have been analysed 
# successfully.

numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)

# File wich stores trends and indices for each combination of species and stratum

all_Indices_All_Trends <- make_All_Indices_All_Trends(overview = overview, listSpeciesStratumCombinations = listSpeciesStratumCombinations)

for(j in 1:numberSpeciesStratumCombinations){

# The file with arguments only contains the arguments to run the analysis for a 
# particular combination of species and stratum (stratum is e.g. a region). 
# These arguments are used when calling the function 'rtrim'.

  arguments <- read.csv2(listSpeciesStratumCombinations[j], header = TRUE, stringsAsFactors = FALSE)
  counts <- read.csv2(paste(arguments$File, "_counts.csv", sep = ""), header = TRUE)
  
# Also the result of the rtrim-analysis is required.
  load(paste(arguments$File, ".RData", sep = ""))

# MAKING FILES THAT WILL COME IN HANDY. 
# Three different output files are created (dataframes). Filling the files with 
# output is done later.
# To create the file containing indices and time totals (indices_TT_file). The 
# file with the results is needed as argument.
# The file arg_output contains the slopes and arguments used to call the rtrim 
# function.

  indices_TT_file <- make_Indices_TT_file(result = result)
  arg_output_file <- make_arg_output_file(arguments = arguments)

  indices_TT_file <- fill_Indices_TT_file(indices_TT_file = indices_TT_file, result = result, arguments = arguments)

  arg_output_file <- fill_arg_output_file(arg_output_file = arg_output_file, result = result, counts = counts)

  all_Indices_All_Trends <- fill_All_Indices_All_Trends(all_Indices_All_Trends = all_Indices_All_Trends, result = result, arguments = arguments, j = j, listSpeciesStratumCombinations = listSpeciesStratumCombinations)
  
  covariant_matrix <- vcov(result)
  if (arguments$Save_fitted_values){
    FI <- results(result)
  }

# STORING THE OUTPUT FILES. 
# Indices and time totals
# name & file extension

  naam_Indices_TT_file <- paste(arguments$File, "_indices_TT.csv", sep = "")
  
# storage

  write.csv2(indices_TT_file, naam_Indices_TT_file, row.names = FALSE)

# Slopes entire period and arguments
# name & file extension

  naam_arg_output_file <- paste(arguments$File, "_arg_output.csv", sep = "")
  
# storage

  write.csv2(arg_output_file, naam_arg_output_file, row.names = FALSE)

# covariant matrix
# name & file extension

  naam_covariant_matrix <- paste(arguments$File, "_ocv.csv", sep = "")
  
# storage

  write.csv2(covariant_matrix, naam_covariant_matrix, row.names = FALSE)

# File with fitted values

  if (arguments$Save_fitted_values){
  
# storage

    naam_Fitted_Values_File <- paste(arguments$File, "_fitted_values.csv", sep = "")
    write.csv2(FI, naam_Fitted_Values_File, row.names = FALSE)
  }
}

all_Indices_All_Trends <- all_Indices_All_Trends[order(all_Indices_All_Trends$Species_number, all_Indices_All_Trends$Stratum_number, all_Indices_All_Trends$Recordtype_number), ]
write.csv2(all_Indices_All_Trends, "all_Indices_All_Trends.csv", row.names = FALSE)
```

### RTRIM postprocessing
```{r}
Code <- unique(all_Indices_All_Trends$Species_number)

# gather imputed COUNTRYWIDE species' ABUNDANCE time series into one file

df.imputed <- data.table(Code=integer(), Species=factor(), CountryGroup=factor(),
                         Year=integer(), Abundance=integer(), data=character(),
                         First_year=integer())
for (k in Code){
  print(k)
  import <- setDT(read.csv2(file = paste0("BIRD_",k,"_0_fitted_values.csv")))
  import[, data := "observed"]
  import[is.na(observed), data := "imputed"]
  import <- merge(import, stratum_number, by = "site", all.x = TRUE)
  import[, Code := k]
  import[, Species := as.character(species_code[Code==k, Species])]
  import[, First_year := min(import[,time])]
  import[, imputed := round(imputed)]
  setnames(import, old = c("time","imputed","code_sp"), new = c("Year","Abundance","CountryGroup"))
  transfer <- import[,.(Code,Species,CountryGroup,Year,Abundance,data,First_year)]
  df.imputed <- rbind(df.imputed, transfer)
}

write.csv2(df.imputed, file = "Imputed-abundance.csv", row.names = FALSE)

# gather imputed SUPRANATIONAL species' INDEX time series into one file

SI.imputed <- data.table(Code=integer(), Species=factor(), Year=integer(), Abundance=integer(),
                         Abundance_SE=integer(), Index_imputed=integer(), Index_imputed_SE=integer())

for (k in Code){
  import <- setDT(read.csv2(file = paste0("BIRD_",k,"_0_indices_TT.csv")))
  import[, Code := k]
  import[, Species := as.character(species_code[Code==k, Species])]
  import[, Index_imputed := round(Index_imputed*100)][, Index_imputed_SE := round(Index_imputed_SE*100)]
  setnames(import, old = c("TT_imputed","TT_imputed_SE"), new = c("Abundance","Abundance_SE"))
  transfer <- import[,.(Code,Species,Year,Abundance,Abundance_SE,Index_imputed,Index_imputed_SE)]
  SI.imputed <- rbind(SI.imputed, transfer)
}
```

## Estiamting abundance dynamics of the avifauna

### Functions accounting for non-linearity and sampling error

See details in Rigal et al. (2020).
```{r}

# Load functions to estimate dynamics

source("Nonlinear_functions.R")

```

### Estimate dynamics for the common avifauna
```{r}
# Load data

df_pop2<-SI.imputed

# Compute dynamic estimation

niter <- 1000 # number of MCMC iterations
df_pop3b <- msi_fun3(df_pop2, ref_year = 1980, niter = niter, ref_value = "Abundance")

# Plot

ggplot(df_pop3b$msi, aes(x = c(1980:2016), y = mean_msi_final)) +
  geom_ribbon(aes(ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final), alpha=0.5, fill = "lightgrey") +
  geom_point() + theme_modern(base_size = 20) +
  labs(x ="Year", y = "Abundance") +
  stat_function(fun = function(x){df_pop3b$coef$alpha2*x^2 + df_pop3b$coef$alpha1*x + df_pop3b$coef$inter}) +
  stat_function(fun = function(x){df_pop3b$coef$alpha2*x^2 + df_pop3b$coef$alpha1*x + df_pop3b$coef$inter - 1.96*df_pop3b$coef$strd}, linetype = "dashed") +
  stat_function(fun = function(x){df_pop3b$coef$alpha2*x^2 + df_pop3b$coef$alpha1*x + df_pop3b$coef$inter + 1.96*df_pop3b$coef$strd}, linetype = "dashed")
```

### Estimate dynamics by habitat

#### Farmland birds

```{r}
# Load data

pecbms_hab <- read.csv2("Habitat_class_PECBMS.csv") # availble on https://pecbms.info/

# Select species

df_agri_pec <-  droplevels(subset(df_pop2, Species %in% pecbms_hab$Species[pecbms_hab$Habitat=="Farmland"]))

# Compute dynamic estimation

msi_agri_pec_ab<- msi_fun3(df_agri_pec, ref_year = 1980, niter = niter, ref_value = "Abundance")

# Plot

ggplot(msi_agri_pec_ab$msi, aes(x = c(1980:2016), y = mean_msi_final)) +
  geom_ribbon(aes(ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final), alpha=0.5, fill = "lightgrey") +
  geom_point() + theme_modern(base_size = 20) +
  labs(x ="Year", y = "Abundance") +
  stat_function(fun = function(x){msi_agri_pec_ab$coef$alpha2*x^2 + msi_agri_pec_ab$coef$alpha1*x + msi_agri_pec_ab$coef$inter}) +
  stat_function(fun = function(x){msi_agri_pec_ab$coef$alpha2*x^2 + msi_agri_pec_ab$coef$alpha1*x + msi_agri_pec_ab$coef$inter - 1.96*msi_agri_pec_ab$coef$strd}, linetype = "dashed") +
  stat_function(fun = function(x){msi_agri_pec_ab$coef$alpha2*x^2 + msi_agri_pec_ab$coef$alpha1*x + msi_agri_pec_ab$coef$inter + 1.96*msi_agri_pec_ab$coef$strd}, linetype = "dashed")

```

#### Forest birds

```{r}

# Select species

df_forest_pec <-  droplevels(subset(df_pop2, Species %in% pecbms_hab$Species[pecbms_hab$Habitat=="Forest"]))

# Compute dynamic estimation

msi_forest_pec_ab<- msi_fun3(df_forest_pec,ref_year = 1980, niter = niter, ref_value = "Abundance")

# Plot

ggplot(msi_forest_pec_ab$msi, aes(x = c(1980:2016), y = mean_msi_final)) +
  geom_ribbon(aes(ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final), alpha=0.5,fill = "lightgrey") +
  geom_point() + theme_modern(base_size = 20) +
  labs(x ="Year", y = "Abundance") +
  stat_function(fun = function(x){msi_forest_pec_ab$coef$alpha2*x^2 + msi_forest_pec_ab$coef$alpha1*x + msi_forest_pec_ab$coef$inter}) +
  stat_function(fun = function(x){msi_forest_pec_ab$coef$alpha2*x^2 + msi_forest_pec_ab$coef$alpha1*x + msi_forest_pec_ab$coef$inter - 1.96*msi_forest_pec_ab$coef$strd}, linetype = "dashed") +
  stat_function(fun = function(x){msi_forest_pec_ab$coef$alpha2*x^2 + msi_forest_pec_ab$coef$alpha1*x + msi_forest_pec_ab$coef$inter + 1.96*msi_forest_pec_ab$coef$strd}, linetype = "dashed")

```

#### Bird species breeding in urban area

```{r}
# Load data
# from https://www.eea.europa.eu/data-and-maps/data/linkages-of-species-and-habitat

eunis_hab<-read.csv("species_birds_maes_EU27b.csv")
eunis_hab2<-dcast(eunis_hab[eunis_hab$season=="B",], speciesname ~ codeeco)
eunis_hab3<-data.frame(Species=levels(as.factor(eunis_hab$speciesname)), is_urban=FALSE)

# Update species names

eunis_hab3$Species <- as.character(eunis_hab3$Species)
eunis_hab3$Species[eunis_hab3$Species=="Carduelis chloris"] <- "Chloris chloris"
eunis_hab3$Species[eunis_hab3$Species=="Parus caeruleus"] <- "Cyanistes caeruleus"
eunis_hab3$Species[eunis_hab3$Species=="Dendrocopos medius"] <- "Dendrocoptes medius"
eunis_hab3$Species[eunis_hab3$Species=="Dendrocopos minor"] <- "Dryobates minor"
eunis_hab3$Species[eunis_hab3$Species=="Miliaria calandra"] <- "Emberiza calandra"
eunis_hab3$Species[eunis_hab3$Species=="Carduelis cannabina"] <- "Linaria cannabina"
eunis_hab3$Species[eunis_hab3$Species=="Parus cristatus"] <- "Lophophanes cristatus"
eunis_hab3$Species[eunis_hab3$Species=="Tetrao tetrix"] <- "Lyrurus tetrix"
eunis_hab3$Species[eunis_hab3$Species=="Parus ater"] <- "Periparus ater"
eunis_hab3$Species[eunis_hab3$Species=="Parus montanus"] <- "Poecile montanus"
eunis_hab3$Species[eunis_hab3$Species=="Parus palustris"] <- "Poecile palustris"
eunis_hab3$Species[eunis_hab3$Species=="Hirundo rupestris"] <- "Ptyonoprogne rupestris"
eunis_hab3$Species[eunis_hab3$Species=="Regulus ignicapillus"] <- "Regulus ignicapilla"
eunis_hab3$Species[eunis_hab3$Species=="Carduelis spinus"] <- "Spinus spinus"
eunis_hab3$Species <- as.factor(eunis_hab3$Species)
eunis_hab3$is_urban[eunis_hab3$Species %in% levels(droplevels(eunis_hab$speciesname[eunis_hab$season=="B" & eunis_hab$codeeco=="urban"]))] <- TRUE

# Select species

df_build_eunis <- droplevels(subset(df_pop2, Species %in% levels(as.factor(eunis_hab3$Species[eunis_hab3$is_urban==T]))))

# Compute dynamic estimation

msi_build_eunis_ab<- msi_fun3(df_build_eunis, ref_year = 1980, niter = niter, ref_value = "Abundance")

# Plot

ggplot(msi_build_eunis_ab$msi, aes(x = c(1980:2016), y = mean_msi_final)) +
  geom_ribbon(aes(ymin = mean_msi_final-sd_msi_final, ymax = mean_msi_final+sd_msi_final), alpha=0.5, fill = "lightgrey") +
  geom_point() + theme_modern(base_size = 20) +
  labs(x ="Year", y = "Abundance") +
  stat_function(fun = function(x){msi_build_eunis_ab$coef$alpha2*x^2 + msi_build_eunis_ab$coef$alpha1*x + msi_build_eunis_ab$coef$inter}) +
  stat_function(fun = function(x){msi_build_eunis_ab$coef$alpha2*x^2 + msi_build_eunis_ab$coef$alpha1*x + msi_build_eunis_ab$coef$inter - 1.96*msi_build_eunis_ab$coef$strd}, linetype = "dashed") +
  stat_function(fun = function(x){msi_build_eunis_ab$coef$alpha2*x^2 + msi_build_eunis_ab$coef$alpha1*x + msi_build_eunis_ab$coef$inter + 1.96*msi_build_eunis_ab$coef$strd}, linetype = "dashed")

```

#### Hot and cold dwellers

```{r}
# Load data
# from (Devictor et al., 2012)

sti<-read.csv("STI_Devictor.csv") 
sti_eu<-droplevels(subset(sti, SPECIES %in% df_pop2$Species)) 

# Select species (hot dwellers)

df_temp <- droplevels(subset(df_pop2, Species %in% as.character(sti$SPECIES[sti_eu$STI > 12.695]))) # 20% hotter dwellers

# Compute dynamic estimation (hot dwellers)

msi_temp_ab1<- msi_fun3(df_temp,ref_year = 1980, niter = niter, ref_value = "Abundance")

# Select species (cold dwellers)

df_temp <- droplevels(subset(df_pop2, Species %in% as.character(sti$SPECIES[sti_eu$STI < 11.15]))) # 20% colder dwellers

# Compute dynamic estimation (cold dwellers)

msi_temp_ab2<- msi_fun3(df_temp,ref_year = 1980, niter = niter, ref_value = "Abundance")

# Merge results

msi_temp_ab3<-data.frame(rbind(msi_temp_ab1$msi,msi_temp_ab2$msi),
                         fact=c(rep("Hot dwellers", nrow(msi_temp_ab1$msi)),rep("Cold dwellers", nrow(msi_temp_ab2$msi))),
                         year=rep(c(1980:2016),2))

# Plot

ggplot(msi_temp_ab3, aes(x = year, y = mean_msi_final)) +
  geom_ribbon(aes(ymin = mean_msi_final - sd_msi_final, ymax = mean_msi_final + sd_msi_final, fill=fact), alpha=0.5) +
  scale_fill_grey(start=0.8,end=0.2) +
  geom_point() + theme_modern(base_size = 20)+theme(legend.title = element_blank(), legend.position = c(0.8, 0.8)) +
  labs(x ="Year", y = "Abundance")
```

# Trend analysis

### Species trends

```{r}
res_trend2<-function(dataset, niter, ref_year=NULL, mid, correction){tryCatch(res_trend(dataset, niter, ref_year=NULL, mid, correction),
                                            error=function(e) data.frame(alpha2=NA, alpha1=NA,sd_alpha1=NA, inter=NA, strd=NA, p_1=NA, 
                                                                         sd_p_1=NA, p_2=NA, sd_p_2=NA, p_3=NA, sd_p_3=NA, second_order_pvalue=NA,
                                                                         first_order_pvalue=NA,slope=NA, slope_sd=NA, ref_year=NA, max_shape=NA))}

df_trend<-droplevels(subset(df, start_year<=1997))
df_trend<-droplevels(subset(df_trend, end_year>=2014))
df_trend<-droplevels(subset(df_trend, Year %in% c(1995:2016)))
trend_sp_cty_1995<-ddply(df_trend, .(Species, CountryGroup), .fun=res_trend2, niter=1000, correction=T, mid="first", .parallel = F, .progress = "text")
trend_sp_cty_1995$linear<-substr(trend_sp_cty_1995$max_shape, 1, 8)
trend_sp_cty_1995$linear[trend_sp_cty_1995$linear=="stable_c"]<-"stable"

y<-trend_sp_cty_1995
```

### Pressure data and trends

```{r}
# county boundaries
country_name <- c("Austria","Bulgaria","Croatia","Cyprus","Germany","France","UK", "Belgium","Netherlands","Switzerland","Greece","Hungary","Iceland","Ireland","Italy","Latvia","Lithuania","Luxembourg","Malta","Norway","Poland","Portugal","Romania","Slovakia","Denmark", "Czech Republic","Finland", "Sweden", "Estonia","Slovenia","Spain")

country <- droplevels(subset(map_data("world"), region %in% country_name))
coordinates(country)=~long+lat
proj4string(country)<- CRS("+proj=longlat +datum=WGS84")
country<-spTransform(country,CRS("+init=epsg:27572"))
country2<-droplevels(subset(map_data("world"), region %in% country_name))
country2$long<-country$long
country2$lat<-country$lat
country3<-lapply(split(country2[,c(1:2)], country2$group), Polygon)
country3<-SpatialPolygons(lapply(seq_along(country3),function(i){Polygons(list(country3[[i]]),ID=row.names(country2[!duplicated(country2$group),])[i])}))
country_id<-country2 %>% group_by(region,group) %>% summarize(count=n()) %>% data.frame()
country_id<-country_id[order(country_id$group),]
country3<-unionSpatialPolygons(country3,country_id[,1])
proj4string(country3)<-CRS("+init=epsg:27572")

# artificialisation
clc_1990<-raster("U2000_CLC1990_V2020_20u1.tif") # from https://land.copernicus.eu/pan-european/corine-land-cover
clc_2000<-raster("U2006_CLC2000_V2020_20u1.tif")
clc_2006<-raster("U2012_CLC2006_V2020_20u1.tif")
clc_2012<-raster("U2018_CLC2012_V2020_20u1.tif")
clc_2018<-raster("U2018_CLC2018_V2020_20u1.tif")
country4b<-spTransform(country3, CRS(proj4string(clc_1990)))

country_data<-data.frame(t(rep(NA,length(country4b))))
names(country_data)<-levels(as.factor(country_id$region))
for(i in 1:length(country4b)){
  print(i)
  b <- extent(country4b[i])
  test<-crop(clc_1990, b)
  test[test>39]<-NA
  test[test<=11]<-1
  test[test>11]<-0
  test2<-mask(test,country4b[i])
  country_data[3,i]<-extract(test2,extent(test2), fun=mean, na.rm=T)
  
  test<-crop(clc_2000, b)
  test[test>39]<-NA
  test[test<=11]<-1
  test[test>11]<-0
  test2<-mask(test,country4b[i])
  country_data[4,i]<-extract(test2,extent(test2), fun=mean, na.rm=T)
  
  test<-crop(clc_2006, b)
  test[test>39]<-NA
  test[test<=11]<-1
  test[test>11]<-0
  test2<-mask(test,country4b[i])
  country_data[5,i]<-extract(test2,extent(test2), fun=mean, na.rm=T)
  
  test<-crop(clc_2012, b)
  test[test>39]<-NA
  test[test<=11]<-1
  test[test>11]<-0
  test2<-mask(test,country4b[i])
  country_data[6,i]<-extract(test2,extent(test2), fun=mean, na.rm=T)
  
  test<-crop(clc_2018, b)
  test[test>39]<-NA
  test[test<=11]<-1
  test[test>11]<-0
  test2<-mask(test,country4b[i])
  country_data[7,i]<-extract(test2,extent(test2), fun=mean, na.rm=T)
}
row.names(country_data)[1:5]<-c("clc_1990","clc_2000","clc_2006","clc_2012","clc_2018")

country_data[6,]<-apply(country_data[2:5,],2,function(x){mean(x, na.rm=T)})
row.names(country_data)[6]<-"clc_mean"
country_data[7,]<-apply(country_data[2:5,],2,function(x){
  summary(lm(x~c(2000,2006,2012,2018)))$coef[2,1]})
row.names(country_data)[7]<-"d_clc"
country_data[7,]<-unlist(country_data[7,])/unlist(country_data[2,])

# temperature
r_temp<-brick("tg_ens_mean_0.1deg_reg_v20.0e.nc") # from http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php
end_year<-0
for(i in 1:(2019-1950)){
  print(i)
  beg_year<-end_year+1
  end_year<-end_year+365
  year<-paste0("temp_",1949+i)
  if(i %in% c(seq(1952,2019,4)-1949)){
    end_year<-end_year+1 # account for bisextil years
  }
  assign(year, mean(r_temp[[beg_year:end_year]], na.rm=T))
}

country5<-spTransform(country3, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" ))
country6a<-droplevels(subset(map_data("world"), region %in% country_name))
coordinates(country6a)=~long+lat
proj4string(country6a)<- CRS("+proj=longlat +datum=WGS84")
country6a<-spTransform(country6a,CRS("+init=epsg:27572"))
country6<-droplevels(subset(map_data("world"), region %in% country_name))
country6$long<-country6a$long
country6$lat<-country6a$lat
country6b<-lapply(split(country6[,c(1:2)], country6$group), Polygon)
country6b<-SpatialPolygons(lapply(seq_along(country6b),function(i){Polygons(list(country6b[[i]]),ID=row.names(country6[!duplicated(country6$group),])[i])}))
country_id<-country6 %>% group_by(region,group) %>% summarize(count=n()) %>% data.frame()
country_id<-country_id[order(country_id$group),]
country6b<-unionSpatialPolygons(country6b,country_id[,1])
proj4string(country6b)<-CRS("+init=epsg:27572")

country6c<-lapply(split(country6[,c(1:2)], country6$group), Polygon)
country6c<-SpatialPolygons(lapply(seq_along(country6c),function(i){Polygons(list(country6c[[i]]),ID=row.names(country6[!duplicated(country6$group),])[i])}))
country_id<-country6 %>% group_by(region,group) %>% summarize(count=n()) %>% data.frame()
country_id<-country_id[order(country_id$group),]
country_id$ue<-as.factor(as.character(rep(1,nrow(country_id))))
country6c<-unionSpatialPolygons(country6c,country_id[,"ue"])
proj4string(country6c)<-CRS("+init=epsg:27572")
country6d<-spTransform(country6c, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" ))

temp_1950_2018<-list(temp_1950,temp_1951,temp_1952,temp_1953,temp_1954,temp_1955,temp_1956,temp_1957,temp_1958,temp_1959,
  temp_1960,temp_1961,temp_1962,temp_1963,temp_1964,temp_1965,temp_1966,temp_1967,temp_1968,temp_1969,
  temp_1970,temp_1971,temp_1972,temp_1973,temp_1974,temp_1975,temp_1976,temp_1977,temp_1978,temp_1979,
  temp_1980,temp_1981,temp_1982,temp_1983,temp_1984,temp_1985,temp_1986,temp_1987,temp_1988,temp_1989,
  temp_1990,temp_1991,temp_1992,temp_1993,temp_1994,temp_1995,temp_1996,temp_1997,temp_1998,temp_1999,
  temp_2000,temp_2001,temp_2002,temp_2003,temp_2004,temp_2005,temp_2006,temp_2007,temp_2008,temp_2009,
  temp_2010,temp_2011,temp_2012,temp_2013,temp_2014,temp_2015,temp_2016,temp_2017,temp_2018)

for(i in 1:length(temp_1950_2018)){ # mean by year by country
  for(j in 1:ncol(country_data)){
    country_data[i+7,j]<-mean(extract(temp_1950_2018[[i]],country5[j])[[1]],na.rm=T)
  }
}

row.names(country_data)[8:76]<-paste0("temp",sep="_",1950:2018)

country_data[77,]<-apply(country_data[38:76,],2,function(x){mean(x, na.rm=T)})
row.names(country_data)[77]<-"temp_mean"
country_data[78,]<-apply(country_data[38:76,],2,function(x){summary(lm(x~c(1980:2018)))$coef[2,1]})/country_data[38,]
row.names(country_data)[78]<-"d_temp"

# forest and high input cover data from https://ec.europa.eu/eurostat/fr/data/database

# high input farm cover
sau_country<-read.csv("ef_m_farmleg_1_Data.csv", header = T)
sau_country$Value<-str_replace_all(sau_country$Value, " ", "")
sau_country$Value<-as.numeric(as.character(sau_country$Value))
sau_country<-na.omit(sau_country)
sau_country<-dcast(sau_country, INDIC_AGR+TIME~GEO, fun.aggregate=sum, value.var="Value")
names(sau_country)[c(8,13,34)]<-c("Czech Republic","Germany","UK")

to_merge<-data.frame(TIME=c(2005:2016))
sau_country2<-merge(sau_country[36:40,-1], to_merge, by="TIME", all=T)
sau_country2$Iceland<-rep(sau_country2$Iceland[sau_country2$TIME==2010], nrow(sau_country2))
sau_country2$Montenegro<-rep(sau_country2$Montenegro[sau_country2$TIME==2010], nrow(sau_country2))
sau_country2[is.na(sau_country2)]<-0
rep_0<-function(x, time_vec){
  id_0<-which(x==0)
  id_1<-which(x!=0)
  if(length(id_0)>0){
    if(id_0[1]==1){
      alpha<-(x[id_1[2]]-x[id_1[1]])/(time_vec[id_1[2]]-time_vec[id_1[1]])
      beta<-x[id_1[1]]-alpha*time_vec[id_1[1]]
      x[1]<-alpha*time_vec[1]+beta
    }
    if(id_0[length(id_0)]==length(x)){
      alpha<-(x[id_1[length(id_1)]]-x[id_1[(length(id_1)-1)]])/(time_vec[id_1[length(id_1)]]-time_vec[id_1[(length(id_1)-1)]])
      beta<-x[id_1[length(id_1)]]-alpha*time_vec[id_1[length(id_1)]]
      x[length(x)]<-alpha*time_vec[length(x)]+beta
    }
    id_0<-which(x==0)
    id_1<-which(x!=0)
    for(i in 1:(length(id_0))){
      alpha<-(x[min(id_1[which(id_1>id_0[i])])]-x[max(id_1[which(id_1<id_0[i])])])/(time_vec[min(id_1[which(id_1>id_0[i])])]-time_vec[max(id_1[which(id_1<id_0[i])])])
      beta<-x[min(id_1[which(id_1>id_0[i])])]-alpha*time_vec[min(id_1[which(id_1>id_0[i])])]
      x[id_0[i]]<-alpha*time_vec[(id_0[i])]+beta
    } 
  }
  return(x)
}
for(j in 2:ncol(sau_country2)){
  sau_country2[,j]<-rep_0(sau_country2[,j], sau_country2$TIME)
}

input_country<-read.csv("aei_ps_inp_1_Data.csv", header = T)
input_country$Value<-str_replace_all(input_country$Value, " ", "")
input_country$Value<-as.numeric(as.character(input_country$Value))
input_country<-na.omit(input_country)
input_country<-dcast(input_country, INDIC_AG+TIME~GEO, fun.aggregate=sum, value.var="Value")
names(input_country)[c(8,14,31)]<-c("Czech Republic","Germany","UK")

input_country2<-input_country[1:10,-c(1,11)]
sau_country3<-sau_country2[sau_country2$TIME %in% c(2007:2016), which(names(sau_country2) %in% names(input_country2))]
for(i in 2:ncol(input_country2)){
  input_country2[,i]<-input_country2[,i]/sau_country3[,i]
  to_add<-summary(lm(input_country2[which(input_country2[,i]>0),i]~input_country2[which(input_country2[,i]>0),1]))$coef
  input_country2[11,i]<-mean(input_country2[which(input_country2[,i]>0),i], na.rm=T)
  input_country2[12,i]<-to_add[2,1]
  input_country2[13,i]<-to_add[2,4]
}
vect_trans<-input_country2[11,which(colnames(input_country2) %in% names(country6b))]
country_data[79,]<-unlist(c(vect_trans[1:13],0,vect_trans[14:20],0,vect_trans[21:27],0,vect_trans[28]))
vect_trans<-input_country2[12,which(colnames(input_country2) %in% names(country6b))]
vect_trans[which(input_country2[13,which(colnames(input_country2) %in% names(country6b))]>0.05)]<-0
country_data[80,]<-unlist(c(vect_trans[1:13],0,vect_trans[14:20],0,vect_trans[21:27],0,vect_trans[28]))
country_data[80,]<-unlist(country_data[83,])/unlist(country_data[82,])
row.names(country_data)[79:80]<-c("high_input_cover","d_hic")

# forest cover
foret_country<-read.csv("for_area_1_Data.csv", header = T)
foret_country$Value<-str_replace_all(foret_country$Value, " ", "")
foret_country$Value<-as.numeric(as.character(foret_country$Value))
foret_country<-na.omit(foret_country)
foret_country<-dcast(foret_country, INDIC_FO+TIME~GEO, fun.aggregate=sum, value.var="Value")
names(foret_country)[c(9,16,41)]<-c("Czech Republic","Germany","UK")

country_data[82,]<-t(foret_country[6,which(colnames(foret_country) %in% names(country6b))])/t(foret_country[1,which(colnames(foret_country) %in% names(country6b))])
row.names(country_data)[82]<-"forest"
country_data[83,]<-apply(foret_country[2:6,which(colnames(foret_country) %in% names(country6b))],2,function(x){summary(lm(x~c(1990,2000,2005,2010,2015)))$coef[2,1]})/foret_country[2,which(colnames(foret_country) %in% names(country6b))]
row.names(country_data)[83]<-"d_forest"

country_data2<-as.data.frame(t(country_data))
country_data2$country<-as.factor(row.names(country_data2))

for(i in c(1:83)){
  country_data2[which(country_data2[,i]==0),i]<-NA
}
```

### Preparing data for partial least square regression (PLS)
```{r}

country_data2$country2<-gsub(" ", "_", country_data2$country)
country_data2$country2[country_data2$country2=="UK"]<-"United Kingdom"

ssi_eu<-read.csv("SSI_EU.csv") # LeViol 2012
sxi<-read.csv("SXI_EU.csv")
species_name_data<-read.csv("species_name_data.csv", header=T)

global_data2<-merge(y,sxi, by.x="Species", by.y="Nom_Europe",all.x=T)
global_data2<-merge(global_data2,sti, by.x="Species", by.y="SPECIES",all.x=T)
global_data2<-merge(global_data2,pecbms_hab, by="Species",all.x=T)
global_data2<-merge(global_data2,synanthrop, by.x="Species", by.y="sp_name",all.x=T)
global_data2<-merge(global_data2,ssi_eu, by="Species",all.x=T)

trait2<-read.csv("life_history_bird2.csv",header = TRUE)
trait2$is_migrant<-rep(0, nrow(trait2))
trait2$is_migrant[which(trait2$Short.distance.migrant==1 | trait2$Long.distance.migrant==1)]<-1

global_data2b<-merge(global_data2,trait2[,c(4,68:93)], by="Species",all.x=T)
global_data2b[,c("STI")]<-scale(global_data2b[,c("STI")])
global_data2b$slope[global_data2b$first_order_pvalue>0.05]<-0 # removing non significant trends
global_data2b$slope[global_data2b$slope<0]<- -log( -global_data2b$slope[global_data2b$slope<0] + 1) # Yeo-Johnson power transformation
global_data2b$slope[global_data2b$slope>0]<- log(global_data2b$slope[global_data2b$slope>0] + 1)
global_data2b$is_farmland<-as.factor(global_data2b$Habitat=="Farmland")
global_data2b$is_forest<-as.factor(global_data2b$Habitat=="Forest")

global_data3<-merge(global_data2b,country_data2, by.x="CountryGroup", by.y="country2",all.x=T)
global_data3b<-global_data3

global_data3b[,c("temp_mean","d_temp","high_input_cover","d_hic","forest","d_forest","clc_mean","d_clc")]<-scale(global_data3b[,c("temp_mean","d_temp","high_input_cover","d_hic","forest","d_forest","clc_mean","d_clc")])
global_data3b<-merge(global_data3b, centroid_b[,c("Country","lon2","lat2")], by.x="CountryGroup",by.y="Country",all.x=T)
```

### Applying PLS
```{r}
gb_test<-global_data3b[, c("slope","clc_mean","d_clc","temp_mean","d_temp","high_input_cover","d_hic","forest","d_forest")]
gb_test$slope<-scale(gb_test$slope)

cv.modpls<-cv.plsR(gb_test$slope,gb_test[,-1],nt=10)
res.cv.modpls<-cvtable(summary(cv.modpls))
res1<-plsR(gb_test$slope,gb_test[,-1], nt=10, typeVC="adaptative", pvals.expli=TRUE) # adaptative as NA in data
colSums(res1$pvalstep)

# searching the best number of component to keep via CV
cv.modpls<-cv.plsR(slope~.,data=gb_test,nt=10,NK=100)
res.cv.modpls=cvtable(summary(cv.modpls)) 

# using CV PRESS, 1 or 2 components must be kept

# PLS with 2 components
res<-plsR(slope~.,data=gb_test,nt=2,pvals.expli=TRUE)
trend.bootYT1=bootpls(res,typeboot="fmodel_np",R=10000)
temp.ci=confints.bootpls(trend.bootYT1,indices=2:ncol(gb_test))

# PLS with 1 component
resb<-plsR(slope~.,data=gb_test,nt=1,pvals.expli=TRUE)
trend.bootYT1b=bootpls(resb,typeboot="fmodel_np",R=10000)
temp.cib=confints.bootpls(trend.bootYT1b,indices=2:ncol(gb_test))

# using the empirical distribution of the best number of component, we can obtain an empircal measure of the significance of each effect.
ind.BCa.YT1 <- (temp.ci[,7]<0&temp.ci[,8]<0)|(temp.ci[,7]>0&temp.ci[,8]>0) 
ind.BCa.YT1b <- (temp.cib[,7]<0&temp.cib[,8]<0)|(temp.cib[,7]>0&temp.cib[,8]>0)
(matind=(rbind(YT1b=ind.BCa.YT1b, YT1=ind.BCa.YT1)))
pi.e=prop.table(res.cv.modpls$CVPress)[1:2]%*%matind
signpred(t(matind),labsize=.5, plotsize = 12)

coef_plot<-data.frame(var=c("Artificialised cover","Artificialisation trend","Mean temperature","Temperature trend",
                            "High input farm cover", "High input farm cover trend",
                            "Forest cover","Forest cover trend"),val=Pine.bootYT1$t0[-1,1],
                      inf=temp.ci[,1],sup=temp.ci[,2],t(matind), sig=t(pi.e))
coef_plot$col_val<-"ns"
coef_plot$col_val[which(coef_plot$sig>=0.95 & coef_plot$val>0)]<-"pos"
coef_plot$col_val[which(coef_plot$sig>=0.95 & coef_plot$val<0)]<-"neg"
coef_plot$var<-as.character(coef_plot$var)
coef_plot$var<-factor(coef_plot$var, levels = c("Temperature trend","Mean temperature","Artificialisation trend","Artificialised cover", 
                                                "Forest cover trend","Forest cover", "High input farm cover trend","High input farm cover"))
ggplot(coef_plot, aes(y=val, x=var))+
  geom_rect(fill = "#CECEF6",xmin = -Inf,xmax = Inf,    ymin = 0,ymax = Inf, alpha = 0.1) +
  geom_rect(fill = "#F5A9A9",xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = 0, alpha = 0.1) +
  geom_bar(stat="identity", position=position_dodge(), alpha=0.7,aes(fill=var)) +
  scale_fill_manual(values=c("High input farm cover"="#D302F9","High input farm cover trend"="#D302F9",
                             "Artificialised cover"="#196DF6","Artificialisation trend"="#196DF6",
                             "Forest cover"="#1BAE20","Forest cover trend"="#1BAE20",
                             "Mean temperature"="#FA0900","Temperature trend"="#FA0900"))+
  geom_errorbar(aes(ymin=inf, ymax=sup), width=.2,position=position_dodge(.9))+
  theme_modern()+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  geom_hline(yintercept=0, linetype="dashed", size=1) + 
  coord_flip()
```

# Causal effect analysis
### Functions for CCM and Smap
```{r}
CCM_EU<-function(x, column, niter, press=c("short","mid","long")){
  
  x<-droplevels(x[which(!is.na(x$value)),])

  if(length(levels(x$CountryGroup))>1 & nrow(x)>=50){
    
    long_ts<-as.data.frame(x %>% group_by(CountryGroup) %>% summarize(count=n(),ab=sum(Index),sd_i=sd(Index)))
    x<-droplevels(subset(x, CountryGroup %in% long_ts$CountryGroup[which(long_ts$count>3 & long_ts$ab>0 & long_ts$sd_i>0)]))
    long_ts2<-as.data.frame(x %>% group_by(CountryGroup) %>% summarize(count=n()))
    
    x<-data.frame(x %>% group_by(CountryGroup) %>% mutate(Index2=scale(Index),value2=scale(value)))
    ab<-ddply(x, .(CountryGroup), .fun=add_row, .parallel = F)
    Accm<-ab$Index[-nrow(ab)]
    val<-ddply(x, .(CountryGroup), .fun=add_row, .parallel = F)
    Bccm<-val$value[-nrow(val)]
    
    if(press=="short"){maxE<-2}
    if(press=="mid"){maxE<-min(c(min(long_ts2$count),4))}
    if(press=="long"){maxE<-min(c(min(long_ts2$count),6))}
    Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("A", "B")

    for(E in 2:maxE) {
      Emat[E-1,"A"]<-SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho 
      Emat[E-1,"B"]<-SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
    }
    
    E_A<-which.max(na.omit(Emat[,1]))+1
    E_B<-which.max(na.omit(Emat[,2]))+1
    
    if(length(E_A)>0 & length(E_B)>0){
      signal_A_out<-SSR_check_signal(A=Accm, E=E_A, tau=1, predsteplist=1:10)
      signal_B_out<-SSR_check_signal(A=Bccm, E=E_B, tau=1, predsteplist=1:10)
      
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]<0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]<0){
        
        CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=niter)
        CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=niter)
        CCM_significance_test<-ccmtest(CCM_boot_A,CCM_boot_B)
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]<0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]>=0){
        
        CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=niter)
        CCM_significance_test<-ccmtest(CCM_boot_A,CCM_boot_A)
        CCM_significance_test[2]<-1
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]>=0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]<0){
        
        CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=niter)
        CCM_significance_test<-ccmtest(CCM_boot_B,CCM_boot_B)
        CCM_significance_test[1]<-1
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]>=0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]>=0){
        CCM_significance_test<-c(1,1)
      }
      
      result<-data.frame(a_cause_b=CCM_significance_test[1], # a = species, b= pressure
                         b_cause_a=CCM_significance_test[2])
      x<-as.data.frame(x %>% group_by(CountryGroup) %>% mutate(Index2=scale(Index),value2=scale(value)))
      block<-data.frame(ind=x$Index2, pressure=x$value2)
      if(result[2]>0.05 & result[1]<=0.05){ # b does not cause a but a causes b
        block<-block[,c(2,1)] # effect from a to b
        res_smap<-smap_intb(block)
        strenght1<-res_smap[4]
        strenght2<-0
        min_str2<-firstq_str2<-med_str2<-thirdq_str2<-max_str2<-NA}
      if(result[1]>0.05 & result[2]<=0.05){ # b causes a but a does not cause b
        block<-block # effcet from b to a
        res_smap<-smap_intb(block)
        strenght1<-0
        strenght2<-res_smap[4]
        min_str2<-res_smap[1];firstq_str2<-res_smap[2];med_str2<-res_smap[3];thirdq_str2<-res_smap[5];max_str2<-res_smap[6]}
      if(result[1]<=0.05 & result[2]<=0.05){ # a causes b and b causes a
        res_smap<-smap_intb(block)
        strenght2<-res_smap[4] # effect from b to a 
        min_str2<-res_smap[1];firstq_str2<-res_smap[2];med_str2<-res_smap[3];thirdq_str2<-res_smap[5];max_str2<-res_smap[6]
        block<-block[,c(2,1)]
        strenght1<-smap_intb(block)[4] # effect from a to b
      }
      if(result[1]>0.05 & result[2]>0.05){
        strenght1<-NA
        strenght2<-NA
        min_str2<-firstq_str2<-med_str2<-thirdq_str2<-max_str2<-NA
      }
    }else{
      strenght1<-strenght2<-min_str2<-firstq_str2<-med_str2<-thirdq_str2<-max_str2<-NA
      CCM_significance_test<-c(NA,NA)
    }
  }else{
    strenght1<-strenght2<-min_str2<-firstq_str2<-med_str2<-thirdq_str2<-max_str2<-NA
    CCM_significance_test<-c(NA,NA)
  }
  result2<-data.frame(a_cause_b=CCM_significance_test[1],
                      b_cause_a=CCM_significance_test[2],
                      strenght1, strenght2,min_str2,
                      firstq_str2,med_str2,thirdq_str2,max_str2,
                      E_A=ifelse(length(E_A)>0,E_A,NA),
                      E_B=ifelse(length(E_B)>0,E_B,NA))
  
  return(result2)
}

smap_intb <- function(block){
   # Determine the best theta
    theta.examined <- seq(0, 10, by = 0.1)
    th.test <-
      pforeach(
        i      = theta.examined,
        .c     = rbind,
        .cores = 7
      )({
        th.test0 <- block_lnlp(block, method = "s-map", tp = 1,columns = 2,
                               theta  = i, silent = T, num_neighbors = 0)
      })
    
    best.th <- th.test[th.test$mae == min(th.test$mae), 'theta']
    
    # Perform multivariate S-map to quantify interaction strengths
    smapc.res <- block_lnlp(block, method = "s-map", tp = 1, columns=2,
                            theta  = best.th, num_neighbors = 0, silent = T,
                            save_smap_coefficients = T)
    smapc.tmp <- smapc.res[[1]]$smap_coefficients
    smapc.tmp <- as.data.frame(smapc.tmp)
    
    # Column names
    colnames(smapc.tmp) <- c("effect_of_pressure","constant") 
    smapc.tmp<-smapc.tmp[,1]
    smapc.tmp<-c(summary(na.omit(smapc.tmp)))
  return(smapc.tmp)
}
```
### Pressure time-series
```{r}
country_data_temp<-data.frame(year=1950:2018,country_data[8:76,])
country_data_temp2<-melt(country_data_temp, id.vars="year")

country_data_clc<-data.frame(year=c(1990,2000,2006,2012,2016),country_data[1:5,])
country_data_clc2<-melt(country_data_clc, id.vars="year")
country_data_clc2[c(11,12,121,122,146,151),3]<-NA # remove absurd values

country_data_hic<-data.frame(year=2007:2016,input_country2[1:10,-1])
country_data_hic2<-melt(country_data_hic, id.vars="year")

country_data_forest<-data.frame(year=c(1990,2000,2005,2010,2015),foret_country[2:6,-c(1:2)])
x<-as.data.frame(t(country_data_forest[,-1]))
x$area<-t(foret_country[1,-c(1:2)])
x2<-apply(x[,1:5], 2, function(z)(z/x$area))
country_data_forest[,-1]<-t(x2)
country_data_forest2<-melt(country_data_forest, id.vars="year")
```

### Appying CCM and Smap on species/pressure time-series
```{r}
df_press<-as.data.frame(df)

press<-country_data_temp2
press$country<-as.character(press$variable)
press$country[press$country=="Czech.Republic"]<-"Czech_Republic"
press$country[press$country=="UK"]<-"United_Kingdom"
press$country[press$country=="Ireland"]<-"Republic_of_Ireland"
df_press2<-merge(df_press, press, by.x=c("CountryGroup","Year"), by.y=c("country","year"), all.x=T)
df_press2<-droplevels(df_press2[which(df_press2$Index>=df_press2$Index_SE),])
ccm_temp<-ddply(df_press2, .(Species), .fun = CCM_EU, column="Index", niter=300, press="long", .parallel = F, .progress = "text")
names(ccm_temp)[5]<-"temp"

press<-country_data_clc2
press$country<-as.character(press$variable)
press$country[press$country=="Czech.Republic"]<-"Czech_Republic"
press$country[press$country=="UK"]<-"United_Kingdom"
press$country[press$country=="Ireland"]<-"Republic_of_Ireland"
df_press2<-merge(df_press, press, by.x=c("CountryGroup","Year"), by.y=c("country","year"), all.x=T)
df_press2<-droplevels(df_press2[which(df_press2$Index>=df_press2$Index_SE),])
ccm_clc<-ddply(df_press2, .(Species), .fun = CCM_EU, column="Index", niter=1000, .parallel = F, .progress = "text")
names(ccm_clc)[5]<-"clc"

press<-country_data_hic2
press$country<-as.character(press$variable)
press$country[press$country=="Czech.Republic"]<-"Czech_Republic"
press$country[press$country=="UK"]<-"United_Kingdom"
press$country[press$country=="Ireland"]<-"Republic_of_Ireland"
df_press2<-merge(df_press, press, by.x=c("CountryGroup","Year"), by.y=c("country","year"), all.x=T)
df_press2<-droplevels(df_press2[which(df_press2$Index>=df_press2$Index_SE),])
ccm_hic<-ddply(df_press2, .(Species), .fun = CCM_EU, column="Index", niter=1000, press="mid", .parallel = F, .progress = "text")
names(ccm_hic)[5]<-"hic"

press<-country_data_forest2
press$country<-as.character(press$variable)
press$country[press$country=="Czech.Republic"]<-"Czech_Republic"
press$country[press$country=="UK"]<-"United_Kingdom"
press$country[press$country=="Ireland"]<-"Republic_of_Ireland"
df_press2<-merge(df_press, press, by.x=c("CountryGroup","Year"), by.y=c("country","year"), all.x=T)
df_press2<-droplevels(df_press2[which(df_press2$Index>=df_press2$Index_SE),])
ccm_for<-ddply(df_press2, .(Species), .fun = CCM_EU, column="Index", niter=1000, .parallel = F, .progress = "text")
names(ccm_for)[5]<-"forest"

sp_data<-merge(ccm_temp[,c(1,5)],ccm_clc[,c(1,5)], by="Species", all.x=T)
sp_data<-merge(sp_data,ccm_hic[,c(1,5)], by="Species", all.x=T)
sp_data<-merge(sp_data,ccm_for[,c(1,5)], by="Species", all.x=T)

sp_data<-merge(sp_data,sxi, by.x="Species", by.y="Nom_Europe",all.x=T)
sp_data<-merge(sp_data,sti, by.x="Species", by.y="SPECIES",all.x=T)
sp_data<-merge(sp_data,trait2[,c(4,51,52,68:93)], by="Species",all.x=T)
sp_data<-merge(sp_data,pecbms_hab, by="Species",all.x=T)
sp_data<-merge(sp_data,synanthrop, by.x="Species", by.y="sp_name",all.x=T)
sp_data<-merge(sp_data,ssi_eu, by="Species",all.x=T)

sp_data$temp[is.na(sp_data$temp)]<-0
sp_data$clc[is.na(sp_data$clc)]<-0
sp_data$hic[is.na(sp_data$hic)]<-0
sp_data$forest[is.na(sp_data$forest)]<-0
sp_data$is_forest<-as.factor(sp_data$Habitat=="Forest")
sp_data$is_farmland<-as.factor(sp_data$Habitat=="Farmland")
```

### Applying PLS on pressure influence vs. species traits
```{r}
# Temperature vs traits
trait_inter_data_temp<-sp_data[which(sp_data$temp!=0 | !is.na(sp_data$STI)),c("temp","is_farmland","is_forest","STI","SSI","is_migrant","Granivore_B",
                                                 "Arthropods_B","Other.invertebrates_B","s1")]
trait_inter_data_temp$STI<-scale(trait_inter_data_temp$STI)
trait_inter_data_temp$SSI<-scale(trait_inter_data_temp$SSI)
trait_inter_data_temp$s1<-scale(trait_inter_data_temp$s1)
trait_inter_data_temp$is_farmland<-as.numeric(trait_inter_data_temp$is_farmland)-1
trait_inter_data_temp$is_forest<-as.numeric(trait_inter_data_temp$is_forest)-1

cv.modpls_temp<-cv.plsR(trait_inter_data_temp$temp,trait_inter_data_temp[,-1],nt=10)
res.cv.modpls_temp<-cvtable(summary(cv.modpls_temp))
res1<-plsR(trait_inter_data_temp$temp,trait_inter_data_temp[,-1], nt=10, typeVC="adaptative", pvals.expli=TRUE)
colSums(res1$pvalstep)
cv.modpls_temp<-cv.plsR(temp~.,data=trait_inter_data_temp,nt=10,NK=100)
res.cv.modpls_temp<-cvtable(summary(cv.modpls_temp))
res1<-plsR(temp~.,data=trait_inter_data_temp,nt=1,pvals.expli=TRUE)

temp.bootYT1=bootpls(res1,typeboot="fmodel_np",R=10000)
boxplots.bootpls(temp.bootYT1,indices=2:ncol(trait_inter_data_temp))
temper.ci=confints.bootpls(temp.bootYT1,indices=2:ncol(trait_inter_data_temp))
plots.confints.bootpls(temper.ci,typeIC="BCa",colIC=c("blue","blue","blue","blue"),
                       legendpos ="topright")

res1b<-plsR(temp~.,data=trait_inter_data_temp,nt=2,pvals.expli=TRUE)
temp.bootYT1b=bootpls(res1b,typeboot="fmodel_np",R=10000)
temper.cib=confints.bootpls(temp.bootYT1b,indices=2:ncol(trait_inter_data_temp))

res1c<-plsR(temp~.,data=trait_inter_data_temp,nt=3,pvals.expli=TRUE)
temp.bootYT1c=bootpls(res1c,typeboot="fmodel_np",R=10000)
temper.cic=confints.bootpls(temp.bootYT1c,indices=2:ncol(trait_inter_data_temp))

res1d<-plsR(temp~.,data=trait_inter_data_temp,nt=4,pvals.expli=TRUE)
temp.bootYT1d=bootpls(res1d,typeboot="fmodel_np",R=10000)
temper.cid=confints.bootpls(temp.bootYT1d,indices=2:ncol(trait_inter_data_temp))

res1e<-plsR(temp~.,data=trait_inter_data_temp,nt=5,pvals.expli=TRUE)
temp.bootYT1e=bootpls(res1e,typeboot="fmodel_np",R=10000)
temper.cie=confints.bootpls(temp.bootYT1e,indices=2:ncol(trait_inter_data_temp))

ind.BCa.YT1 <- (temper.ci[,7]<0&temper.ci[,8]<0)|(temper.ci[,7]>0&temper.ci[,8]>0)
ind.BCa.YT1b <- (temper.cib[,7]<0&temper.cib[,8]<0)|(temper.cib[,7]>0&temper.cib[,8]>0)
ind.BCa.YT1c <- (temper.cic[,7]<0&temper.cic[,8]<0)|(temper.cic[,7]>0&temper.cic[,8]>0)
ind.BCa.YT1d <- (temper.cid[,7]<0&temper.cid[,8]<0)|(temper.cid[,7]>0&temper.cid[,8]>0)
ind.BCa.YT1e <- (temper.cie[,7]<0&temper.cie[,8]<0)|(temper.cie[,7]>0&temper.cie[,8]>0)

(matind=(rbind(YT1=ind.BCa.YT1, YT1b=ind.BCa.YT1b, YT1c=ind.BCa.YT1c, YT1d=ind.BCa.YT1d, YT1e=ind.BCa.YT1e)))
pi.e=(prop.table(res.cv.modpls_temp$CVPress)[c(1:5)]/sum(prop.table(res.cv.modpls_temp$CVPress)[c(1:5)]))%*%matind
signpred(t(matind),labsize=.5, plotsize = 12)

coef_plot_temp<-data.frame(var=c("Farmland","Forest","STI","SSI","Migrant",
                                 "Granivorous diet","Arthropod diet","Other invertebrate diet","Synanthropy"),val=temp.bootYT1c$t0[-1,1],
                           inf=temper.cic[,7],sup=temper.cic[,8],t(matind),sig=t(pi.e))
coef_plot_temp$col_val<-"ns"
coef_plot_temp$col_val[which(coef_plot_temp$sig>=0.95 & coef_plot_temp$val>0)]<-"pos"
coef_plot_temp$col_val[which(coef_plot_temp$sig>=0.95 & coef_plot_temp$val<0)]<-"neg"

# Artificialisation vs traits

trait_inter_data_clc<-sp_data[which(sp_data$clc!=0 | !is.na(sp_data$STI)),c("clc","is_farmland","is_forest","STI","SSI","is_migrant","Granivore_B",
                                                 "Arthropods_B","Other.invertebrates_B","s1")] 
trait_inter_data_clc$STI<-scale(trait_inter_data_clc$STI)
trait_inter_data_clc$SSI<-scale(trait_inter_data_clc$SSI)
trait_inter_data_clc$s1<-scale(trait_inter_data_clc$s1)
trait_inter_data_clc$is_farmland<-as.numeric(trait_inter_data_clc$is_farmland)-1
trait_inter_data_clc$is_forest<-as.numeric(trait_inter_data_clc$is_forest)-1

cv.modpls_clc<-cv.plsR(trait_inter_data_clc$clc,trait_inter_data_clc[,-1],nt=10)
res.cv.modpls_clc<-cvtable(summary(cv.modpls_clc))
res2<-plsR(trait_inter_data_clc$clc,trait_inter_data_clc[,-1], nt=10, typeVC="adaptative", pvals.expli=TRUE) 
colSums(res2$pvalstep)
cv.modpls_clc<-cv.plsR(clc~.,data=trait_inter_data_clc,nt=10,NK=100)
res.cv.modpls_clc<-cvtable(summary(cv.modpls_clc))
res2<-plsR(clc~.,data=trait_inter_data_clc,nt=1,pvals.expli=TRUE)

clc.bootYT1=bootpls(res2,typeboot="fmodel_np",R=10000)
boxplots.bootpls(clc.bootYT1,indices=2:ncol(trait_inter_data_clc))
clcer.ci=confints.bootpls(clc.bootYT1,indices=2:ncol(trait_inter_data_clc))
plots.confints.bootpls(clcer.ci,typeIC="BCa",colIC=c("blue","blue","blue","blue"),
                       legendpos ="topright")

res2b<-plsR(clc~.,data=trait_inter_data_clc,nt=3,pvals.expli=TRUE)
clc.bootYT1b=bootpls(res2b,typeboot="fmodel_np",R=10000)
clcer.cib=confints.bootpls(clc.bootYT1b,indices=2:ncol(trait_inter_data_clc))

res2c<-plsR(clc~.,data=trait_inter_data_clc,nt=4,pvals.expli=TRUE)
clc.bootYT1c=bootpls(res2c,typeboot="fmodel_np",R=10000)
clcer.cic=confints.bootpls(clc.bootYT1c,indices=2:ncol(trait_inter_data_clc))

res2d<-plsR(clc~.,data=trait_inter_data_clc,nt=5,pvals.expli=TRUE)
clc.bootYT1d=bootpls(res2d,typeboot="fmodel_np",R=10000)
clcer.cid=confints.bootpls(clc.bootYT1d,indices=2:ncol(trait_inter_data_clc))

ind.BCa.YT1 <- (clcer.ci[,7]<0&clcer.ci[,8]<0)|(clcer.ci[,7]>0&clcer.ci[,8]>0)
ind.BCa.YT1b <- (clcer.cib[,7]<0&clcer.cib[,8]<0)|(clcer.cib[,7]>0&clcer.cib[,8]>0)
ind.BCa.YT1c <- (clcer.cic[,7]<0&clcer.cic[,8]<0)|(clcer.cic[,7]>0&clcer.cic[,8]>0)
ind.BCa.YT1d <- (clcer.cid[,7]<0&clcer.cid[,8]<0)|(clcer.cid[,7]>0&clcer.cid[,8]>0)
(matind=(rbind(YT1=ind.BCa.YT1,YT1b=ind.BCa.YT1b,YT1c=ind.BCa.YT1c,YT1d=ind.BCa.YT1d)))
pi.e=(prop.table(res.cv.modpls_clc$CVPress)[c(1,3:5)]/sum(prop.table(res.cv.modpls_clc$CVPress)[c(1,3:5)]))%*%matind
signpred(t(matind),labsize=.5, plotsize = 12)

coef_plot_clc<-data.frame(var=c("Farmland","Forest","STI","SSI","Migrant",
                                 "Granivorous diet","Arthropod diet","Other invertebrate diet","Synanthropy"),val=clc.bootYT1c$t0[-1,1],
                           inf=clcer.cic[,7],sup=clcer.cic[,8],t(matind),sig=t(pi.e))
coef_plot_clc$col_val<-"ns"
coef_plot_clc$col_val[which(coef_plot_clc$sig>=0.95 & coef_plot_clc$val>0)]<-"pos"
coef_plot_clc$col_val[which(coef_plot_clc$sig>=0.95 & coef_plot_clc$val<0)]<-"neg"

# High input farm cover vs traits

trait_inter_data_hic<-sp_data[which(sp_data$hic!=0| !is.na(sp_data$STI)), c("hic","is_farmland","is_forest","STI","SSI","is_migrant","Granivore_B",
                                                 "Arthropods_B","Other.invertebrates_B","s1")]
trait_inter_data_hic$STI<-scale(trait_inter_data_hic$STI)
trait_inter_data_hic$SSI<-scale(trait_inter_data_hic$SSI)
trait_inter_data_hic$s1<-scale(trait_inter_data_hic$s1)
trait_inter_data_hic$is_farmland<-as.numeric(trait_inter_data_hic$is_farmland)-1
trait_inter_data_hic$is_forest<-as.numeric(trait_inter_data_hic$is_forest)-1

cv.modpls_hic<-cv.plsR(trait_inter_data_hic$hic,trait_inter_data_hic[,-1],nt=10)
res.cv.modpls_hic<-cvtable(summary(cv.modpls_hic))
res3<-plsR(trait_inter_data_hic$hic,trait_inter_data_hic[,-1], nt=10, typeVC="adaptative", pvals.expli=TRUE) 
colSums(res3$pvalstep)
cv.modpls_hic<-cv.plsR(hic~.,data=trait_inter_data_hic,nt=10,NK=100)
res.cv.modpls_hic<-cvtable(summary(cv.modpls_hic))
plot(res.cv.modpls_hic)
res3<-plsR(hic~.,data=trait_inter_data_hic,nt=1,pvals.expli=TRUE)
biplot(res3$tt,res3$pp)

hic.bootYT1=bootpls(res3,typeboot="fmodel_np",R=10000)
boxplots.bootpls(hic.bootYT1,indices=2:ncol(trait_inter_data_hic))
hicer.ci=confints.bootpls(hic.bootYT1,indices=2:ncol(trait_inter_data_hic))
plots.confints.bootpls(hicer.ci,typeIC="BCa",colIC=c("blue","blue","blue","blue"),
                       legendpos ="topright")

res3b<-plsR(hic~.,data=trait_inter_data_hic,nt=2,pvals.expli=TRUE)
hic.bootYT1b=bootpls(res3b,typeboot="fmodel_np",R=10000)
hicer.cib=confints.bootpls(hic.bootYT1b,indices=2:ncol(trait_inter_data_hic))

ind.BCa.YT1 <- (hicer.ci[,7]<0&hicer.ci[,8]<0)|(hicer.ci[,7]>0&hicer.ci[,8]>0)
ind.BCa.YT1b <- (hicer.cib[,7]<0&hicer.cib[,8]<0)|(hicer.cib[,7]>0&hicer.cib[,8]>0)

(matind=(rbind(YT1=ind.BCa.YT1,YT1b=ind.BCa.YT1b)))
pi.e=(prop.table(res.cv.modpls_hic$CVPress)[1:2]/sum(prop.table(res.cv.modpls_hic$CVPress)[1:2]))%*%matind
pi.e
signpred(t(matind),labsize=.5, plotsize = 12)

coef_plot_hic<-data.frame(var=c("Farmland","Forest","STI","SSI","Migrant",
                                 "Granivorous diet","Arthropod diet","Other invertebrate diet","Synanthropy"),val=hic.bootYT1$t0[-1,1],
                           inf=hicer.ci[,7],sup=hicer.ci[,8],t(matind),sig=t(pi.e))
coef_plot_hic$col_val<-"ns"
coef_plot_hic$col_val[which(coef_plot_hic$sig>=0.95 & coef_plot_hic$val>0)]<-"pos"
coef_plot_hic$col_val[which(coef_plot_hic$sig>=0.95 & coef_plot_hic$val<0)]<-"neg"

# Forest vs traits

trait_inter_data_forest<-sp_data[which(sp_data$forest!=0 | !is.na(sp_data$STI)),c("forest","is_farmland","is_forest","STI","SSI","is_migrant","Granivore_B",
                                                 "Arthropods_B","Other.invertebrates_B","s1")]
trait_inter_data_forest$STI<-scale(trait_inter_data_forest$STI)
trait_inter_data_forest$SSI<-scale(trait_inter_data_forest$SSI)
trait_inter_data_forest$s1<-scale(trait_inter_data_forest$s1)
trait_inter_data_forest$is_farmland<-as.numeric(trait_inter_data_forest$is_farmland)-1
trait_inter_data_forest$is_forest<-as.numeric(trait_inter_data_forest$is_forest)-1

cv.modpls_forest<-cv.plsR(trait_inter_data_forest$forest,trait_inter_data_forest[,-1],nt=10)
res.cv.modpls_forest<-cvtable(summary(cv.modpls_forest))
res4<-plsR(trait_inter_data_forest$forest,trait_inter_data_forest[,-1], nt=10, typeVC="adaptative", pvals.expli=TRUE) 
colSums(res4$pvalstep)
cv.modpls_forest<-cv.plsR(forest~.,data=trait_inter_data_forest,nt=10,NK=100)
res.cv.modpls_forest<-cvtable(summary(cv.modpls_forest))
res4<-plsR(forest~.,data=trait_inter_data_forest,nt=1,pvals.expli=TRUE)

forest.bootYT1=bootpls(res4,typeboot="fmodel_np",R=10000)
boxplots.bootpls(forest.bootYT1,indices=2:ncol(trait_inter_data_forest))
forester.ci=confints.bootpls(forest.bootYT1,indices=2:ncol(trait_inter_data_forest))
plots.confints.bootpls(forester.ci,typeIC="BCa",colIC=c("blue","blue","blue","blue"),
                       legendpos ="topright")

res4b<-plsR(forest~.,data=trait_inter_data_forest,nt=3,pvals.expli=TRUE)
forest.bootYT1b=bootpls(res4b,typeboot="fmodel_np",R=10000)
forester.cib=confints.bootpls(forest.bootYT1b,indices=2:ncol(trait_inter_data_forest))

ind.BCa.YT1 <- (forester.ci[,7]<0&forester.ci[,8]<0)|(forester.ci[,7]>0&forester.ci[,8]>0)
ind.BCa.YT1b <- (forester.cib[,7]<0&forester.cib[,8]<0)|(forester.cib[,7]>0&forester.cib[,8]>0)

(matind=(rbind(YT1=ind.BCa.YT1, YT1b=ind.BCa.YT1b)))
pi.e=(prop.table(res.cv.modpls_forest$CVPress)[c(1,3)]/sum(prop.table(res.cv.modpls_forest$CVPress)[c(1,3)]))%*%matind
signpred(t(matind),labsize=.5, plotsize = 12)

coef_plot_forest<-data.frame(var=c("Farmland","Forest","STI","SSI","Migrant",
                                 "Granivorous diet","Arthropod diet","Other invertebrate diet","Synanthropy"),val=forest.bootYT1$t0[-1,1],
                           inf=forester.ci[,7],sup=forester.ci[,8],t(matind),sig=t(pi.e))
coef_plot_forest$col_val<-"ns"
coef_plot_forest$col_val[which(coef_plot_forest$sig>=0.95 & coef_plot_forest$val>0)]<-"pos"
coef_plot_forest$col_val[which(coef_plot_forest$sig>=0.95 & coef_plot_forest$val<0)]<-"neg"

ggplot(coef_plot_temp, aes(y=val, x=var))+
  geom_bar(stat="identity", position=position_dodge(), aes(alpha=abs(val), fill=col_val)) +
  scale_fill_manual(values=c("ns"="#959393","neg"="#FF0000","pos"="#0000FF"))+
  geom_errorbar(aes(ymin=inf, ymax=sup), width=.2, position=position_dodge(.9))+
  theme_modern()+labs(y="Temperature influence on species dynamic")+
  theme(legend.position = "none", axis.title.y = element_blank())+
  geom_hline(yintercept=0, linetype="dashed", size=1) + 
  coord_flip()

ggplot(coef_plot_clc, aes(y=val, x=var))+
  geom_bar(stat="identity", position=position_dodge(), aes(alpha=abs(val), fill=col_val)) +
  scale_fill_manual(values=c("ns"="#959393","neg"="#FF0000","pos"="#0000FF"))+
  geom_errorbar(aes(ymin=inf, ymax=sup), width=.2, position=position_dodge(.9))+
  theme_modern()+labs(y="Artificialisation influence on species dynamic")+
  theme(legend.position = "none", axis.title.y = element_blank())+
  geom_hline(yintercept=0, linetype="dashed", size=1) + 
  coord_flip()

ggplot(coef_plot_hic, aes(y=val, x=var))+
  geom_bar(stat="identity", position=position_dodge(), aes(alpha=abs(val), fill=col_val)) +
  scale_fill_manual(values=c("ns"="#959393","neg"="#FF0000","pos"="#0000FF"))+
  geom_errorbar(aes(ymin=inf, ymax=sup), width=.2, position=position_dodge(.9))+
  theme_modern()+labs(y="High input farm influence on species dynamic")+
  theme(legend.position = "none", axis.title.y = element_blank())+
  geom_hline(yintercept=0, linetype="dashed", size=1) + 
  coord_flip()

ggplot(coef_plot_forest, aes(y=val, x=var))+
  geom_bar(stat="identity", position=position_dodge(), aes(alpha=abs(val), fill=col_val)) +
  scale_fill_manual(values=c("ns"="#959393","neg"="#FF0000","pos"="#0000FF"))+
  geom_errorbar(aes(ymin=inf, ymax=sup), width=.2, position=position_dodge(.9))+
  theme_modern()+labs(y="Forest influence on species dynamic")+
  theme(legend.position = "none", axis.title.y = element_blank())+
  geom_hline(yintercept=0, linetype="dashed", size=1) + 
  coord_flip()
```

# References
Bogaart, P., van der Loo, M., Pannekoek, J., & Bogaart, M. P. (2016). Package ‘rtrim’.
Devictor, V., Van Swaay, C., Brereton, T., Brotons, L., Chamberlain, D., Heliölä, J., ... & Jiguet, F. (2012). Differences in the climatic debts of birds and butterflies at a continental scale. Nature climate change, 2(2), 121-124.
Brlík, V., Šilarová, E., Škorpilová, J. et al. Long-term and large-scale multispecies dataset tracking population changes of common European breeding birds. Sci Data 8, 21 (2021). https://doi.org/10.1038/s41597-021-00804-2
Guetté, A., Gaüzère, P., Devictor, V., Jiguet, F., & Godet, L. (2017). Measuring the synanthropy of species and communities to monitor the effects of urbanization on biodiversity. Ecological Indicators, 79, 139-154.
Le Viol, I., Jiguet, F., Brotons, L., Herrando, S., Lindström, Å., Pearce-Higgins, J. W., ... & Devictor, V. (2012). More and more generalists: two decades of changes in the European avifauna. Biology letters, 8(5), 780-782.
Rigal, S., Devictor, V., & Dakos, V. (2020). A method for classifying and comparing non-linear trajectories of ecological variables. Ecological Indicators, 112, 106113.