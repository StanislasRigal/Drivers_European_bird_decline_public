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

abd <- setDT(read.table("raw_data/Abundance_data_PECBMS.txt", header = T, sep="\t")) 

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

# Assign countries with a number

stratum_number <- data.frame(code_sp = C, site = 1:length(C))

# Assign species with their EURING code

species_code <- df_pop[, .SD[1], by = .(Code)][, .(Code,Species)]

```

### RTRIM functions

```{r}

# Load RTRIM functions

source("RTRIM_functions.R")

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

# Determine how many combinations of species and stratum have been analysed successfully.

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
  
# Storage

  write.csv2(indices_TT_file, naam_Indices_TT_file, row.names = FALSE)

# Slopes entire period and arguments
# name & file extension

  naam_arg_output_file <- paste(arguments$File, "_arg_output.csv", sep = "")
  
# Storage

  write.csv2(arg_output_file, naam_arg_output_file, row.names = FALSE)

# Covariant matrix
# name & file extension

  naam_covariant_matrix <- paste(arguments$File, "_ocv.csv", sep = "")
  
# Storage

  write.csv2(covariant_matrix, naam_covariant_matrix, row.names = FALSE)

# File with fitted values

  if (arguments$Save_fitted_values){
  
# Storage

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

# Gather imputed COUNTRYWIDE species' ABUNDANCE time series into one file

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

# Gather imputed SUPRANATIONAL species' INDEX time series into one file

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

See Rigal et al. (2020) for more details about the following functions.
```{r}

# Load functions to estimate dynamics

setwd("..")
source("Nonlinear_functions.R")

```

### Estimate dynamics for the common avifauna
```{r}

# Load data

df_pop2 <- SI.imputed

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
  
# Select species

df_agri_pec <-  droplevels(subset(df_pop2, Species %in% pecbms_hab$Species[pecbms_hab$Habitat=="Farmland"]))

# Compute dynamic estimation

msi_agri_pec_ab <- msi_fun3(df_agri_pec, ref_year = 1980, niter = niter, ref_value = "Abundance")

# Plot

ggplot(msi_agri_pec_ab$msi[17:37,], aes(x = c(1996:2016), y = mean_msi_final)) +
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

msi_forest_pec_ab <- msi_fun3(df_forest_pec,ref_year = 1980, niter = niter, ref_value = "Abundance")

# Plot

ggplot(msi_forest_pec_ab$msi[17:37,], aes(x = c(1996:2016), y = mean_msi_final)) +
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

eunis_hab <- read.csv("raw_data/species_birds_maes_EU27b.csv")
eunis_hab2 <- dcast(eunis_hab[eunis_hab$season=="B",], speciesname ~ codeeco)
eunis_hab3 <- data.frame(Species=levels(as.factor(eunis_hab$speciesname)), is_urban=FALSE)

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
eunis_hab3$is_forest2 <- FALSE
eunis_hab3$is_forest2[eunis_hab3$Species %in% levels(droplevels(eunis_hab2$speciesname[eunis_hab2$woodlandForest!=0 & eunis_hab2$heathlandShrub!=0]))] <- TRUE
eunis_hab3$is_farmland2 <- FALSE
eunis_hab3$is_farmland2[eunis_hab3$Species %in% levels(droplevels(eunis_hab2$speciesname[eunis_hab2$cropland!=0 & eunis_hab2$grassland!=0]))] <- TRUE

# Select species

df_build_eunis <- droplevels(subset(df_pop2, Species %in% levels(droplevels(eunis_hab3$Species[eunis_hab3$is_urban==T]))))

# Compute dynamic estimation

msi_build_eunis_ab <- msi_fun3(df_build_eunis, ref_year = 1980, niter = niter, ref_value = "Abundance")

# Plot

ggplot(msi_build_eunis_ab$msi[17:37,], aes(x = c(1996:2016), y = mean_msi_final)) +
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

sti <- read.csv("raw_data/STI_Devictor.csv") 
sti_eu <- droplevels(subset(sti, SPECIES %in% df_pop2$Species)) 

# Select species (hot dwellers)

df_temp <- droplevels(subset(df_pop2, Species %in% as.character(sti_eu$SPECIES[sti_eu$STI > 12.695]))) # 20% hotter dwellers

# Compute dynamic estimation (hot dwellers)

msi_temp_ab1 <- msi_fun3(df_temp,ref_year = 1980, niter = niter, ref_value = "Abundance")

# Select species (cold dwellers)

df_temp <- droplevels(subset(df_pop2, Species %in% as.character(sti_eu$SPECIES[sti_eu$STI < 11.15]))) # 20% colder dwellers

# Compute dynamic estimation (cold dwellers)

msi_temp_ab2 <- msi_fun3(df_temp,ref_year = 1980, niter = niter, ref_value = "Abundance")

# Merge results

msi_temp_ab3 <- data.frame(rbind(msi_temp_ab1$msi,msi_temp_ab2$msi),
                         fact=c(rep("Hot dwellers", nrow(msi_temp_ab1$msi)),rep("Cold dwellers", nrow(msi_temp_ab2$msi))),
                         year=rep(c(1980:2016),2))

# Plot

ggplot(msi_temp_ab3[msi_temp_ab3$year %in% 1996:2016,], aes(x = year, y = mean_msi_final)) +
  geom_ribbon(aes(ymin = mean_msi_final - sd_msi_final, ymax = mean_msi_final + sd_msi_final, fill=fact), alpha=0.5) +
  scale_fill_grey(start=0.8,end=0.2) +
  geom_point() + theme_modern(base_size = 20)+theme(legend.title = element_blank(), legend.position = c(0.8, 0.8)) +
  labs(x ="Year", y = "Abundance")
  
```

## Pressures

### Geographical data

```{r}

# Select countries involved in the PECBMS

country_name <- c("Austria","Bulgaria","Croatia","Cyprus","Germany","France","UK", "Belgium","Netherlands","Switzerland","Greece","Hungary","Iceland","Ireland","Italy","Latvia","Lithuania","Luxembourg","Malta","Norway","Poland","Portugal","Romania","Slovakia","Denmark", "Czech Republic","Finland", "Sweden", "Estonia","Slovenia","Spain")

# County boundaries in WGS84 and Lambert II

country <- droplevels(subset(map_data("world"), region %in% country_name))
coordinates(country) = ~long+lat
proj4string(country) <- CRS("+proj=longlat +datum=WGS84")
country <- spTransform(country,CRS("+init=epsg:27572"))

country2 <- droplevels(subset(map_data("world"), region %in% country_name))
country2$long <- country$long
country2$lat <- country$lat

country3 <- lapply(split(country2[,c(1:2)], country2$group), Polygon)
country3 <- SpatialPolygons(lapply(seq_along(country3),function(i){Polygons(list(country3[[i]]),ID=row.names(country2[!duplicated(country2$group),])[i])}))

country_id <- country2 %>% group_by(region,group) %>% summarize(count=n()) %>% data.frame()
country_id <- country_id[order(country_id$group),]

country3 <- unionSpatialPolygons(country3,country_id[,1])
proj4string(country3) <- CRS("+init=epsg:27572")

```

### Urban cover

```{r}

# Load data
# from https://www.fao.org/faostat/en/#data/LC

fao_data_landcover <- read.csv("raw_data/FAOSTAT_data_LC.csv", header = T)
fao_data_landcover$Area <- as.character(fao_data_landcover$Area)
fao_data_landcover[fao_data_landcover=="United Kingdom of Great Britain and Northern Ireland"] <- "UK"
fao_data_landcover[fao_data_landcover=="Czechia"] <- "Czech Republic"
fao_data_landcover <- droplevels(fao_data_landcover[fao_data_landcover$Area %in% country_name,])
write.csv(fao_data_landcover,"output/fao_data_landcover.csv", row.names = F)

# from https://www.fao.org/faostat/en/#data/RL

fao_data_landuse <- read.csv("raw_data/FAOSTAT_data_RL.csv", header = T)
fao_data_landuse$Area <- as.character(fao_data_landuse$Area)
fao_data_landuse[fao_data_landuse=="United Kingdom of Great Britain and Northern Ireland"] <- "UK"
fao_data_landuse[fao_data_landuse=="Czechia"] <- "Czech Republic"
fao_data_landuse <- droplevels(fao_data_landuse[fao_data_landuse$Area %in% country_name,])
write.csv(fao_data_landuse,"output/fao_data_landuse.csv", row.names = F)

area_country <- read.csv("output/fao_data_landuse.csv", header = T)
area_country <- area_country[area_country$Item=="Country area" & area_country$Year==2016, c("Area","Value")]
area_country$Value <- 10*area_country$Value
area_country <- area_country[order(area_country$Area),]
area_country$Value[area_country$Area=="Norway"] <- 385178 # there is an abrubt change for Norway in data series of surface area from 385,000 to 625,000 sq km after 2013

# Load data
# from https://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=lan_lcv_art&lang=en

urban_country <- read.csv("raw_data/lan_lcv_art_1_Data.csv", header = T)
urban_country$Value <- as.numeric(as.character(urban_country$Value))
urban_country <- dcast(urban_country, TIME~GEO, fun.aggregate=sum, value.var="Value")
names(urban_country)[c(7,12,29)] <- c("Czech Republic", "Germany","UK")
urban_country$TIME <- NULL

urban_country <- as.data.frame(apply(urban_country,2,function(x){
                                      a1 <- (x[2]-x[1])/(2012-2009)
                                      b1 <- x[2]-a1*2012
                                      a2 <- (x[3]-x[2])/(2015-2012)
                                      b2 <- x[3]-a2*2015
                                      a3 <- (x[4]-x[3])/(2018-2015)
                                      b3 <- x[4]-a3*2018
                                      y <- c(x[1], a1*2010+b1, a1*2011+b1, x[2], a2*2013+b2, a2*2014+b2, x[3], a3*2016+b3)
                                      return(y)
}))

row.names(urban_country)[1:8] <- paste0("urb_", 2009:2016)

# Mean value over the period

urban_country[9,] <- apply(urban_country,2,function(x){mean(x, na.rm=T)})
row.names(urban_country)[9] <- "urb_mean"

# Trend over the period

urban_country[10,] <- apply(urban_country[1:8,],2,function(x){
  summary(lm(x~c(2009:2016)))$coef[2,1]})/urban_country[1,]
row.names(urban_country)[10] <- "d_urb"

urban_country <- data.frame(urban_country[,c(1:13)], Iceland=0, urban_country[,c(14:20)], Norway=0, urban_country[,c(21:27)], Switzerland=0, UK=urban_country$UK)

names(urban_country)[6] <- "Czech Republic"

# Dataset to merge with other pressures

country_data <- urban_country

```

### Temperature

```{r}

# Load data
# from http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php

r_temp <- brick("raw_data/tg_ens_mean_0.1deg_reg_v20.0e.nc") 

# Average daily data by year

end_year <- 0
for(i in 1:(2019-1950)){
  print(i)
  beg_year <- end_year+1
  end_year <- end_year+365
  year <- paste0("temp_",1949+i)
  if(i %in% c(seq(1952,2019,4)-1949)){
    end_year <- end_year+1 # account for bisextil years
  }
  map_year <- assign(year, mean(r_temp[[beg_year:end_year]], na.rm=T))
  path <- paste0("output/temp_",1949+i,".tif")
  writeRaster(x=map_year, filename=path, overwrite=T)
}

for(i in 1:(2019-1950)){
  print(i)
  year <- paste0("temp_",1949+i)
  path <- paste0("output/temp_",1949+i,".tif")
  assign(year, raster(path))
}

# Reproject county boundaries

country5 <- spTransform(country3, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" ))

country6a <- droplevels(subset(map_data("world"), region %in% country_name))
coordinates(country6a) <- ~long+lat
proj4string(country6a) <- CRS("+proj=longlat +datum=WGS84")
country6a <- spTransform(country6a,CRS("+init=epsg:27572"))

country6 <- droplevels(subset(map_data("world"), region %in% country_name))
country6$long <- country6a$long
country6$lat <- country6a$lat

country6b <- lapply(split(country6[,c(1:2)], country6$group), Polygon)
country6b <- SpatialPolygons(lapply(seq_along(country6b),function(i){Polygons(list(country6b[[i]]),ID=row.names(country6[!duplicated(country6$group),])[i])}))

country_id <- country6 %>% group_by(region,group) %>% summarize(count=n()) %>% data.frame()
country_id <- country_id[order(country_id$group),]

country6b <- unionSpatialPolygons(country6b,country_id[,1])
proj4string(country6b) <- CRS("+init=epsg:27572")

# List dataset of yearly temperature

temp_1950_2018 <- list(temp_1950,temp_1951,temp_1952,temp_1953,temp_1954,temp_1955,temp_1956,temp_1957,temp_1958,temp_1959,
  temp_1960,temp_1961,temp_1962,temp_1963,temp_1964,temp_1965,temp_1966,temp_1967,temp_1968,temp_1969,
  temp_1970,temp_1971,temp_1972,temp_1973,temp_1974,temp_1975,temp_1976,temp_1977,temp_1978,temp_1979,
  temp_1980,temp_1981,temp_1982,temp_1983,temp_1984,temp_1985,temp_1986,temp_1987,temp_1988,temp_1989,
  temp_1990,temp_1991,temp_1992,temp_1993,temp_1994,temp_1995,temp_1996,temp_1997,temp_1998,temp_1999,
  temp_2000,temp_2001,temp_2002,temp_2003,temp_2004,temp_2005,temp_2006,temp_2007,temp_2008,temp_2009,
  temp_2010,temp_2011,temp_2012,temp_2013,temp_2014,temp_2015,temp_2016,temp_2017,temp_2018)
  
# Extract temperature values by year by country

for(i in 1:length(temp_1950_2018)){ # mean by year by country
  for(j in 1:ncol(country_data)){
    country_data[i+10,j] <- mean(extract(temp_1950_2018[[i]],country5[j])[[1]],na.rm=T)
  }
}
row.names(country_data)[11:79] <- paste0("temp",sep="_",1950:2018)

# Mean value over the period

country_data[80,] <- apply(country_data[57:77,], 2, function(x){mean(x, na.rm=T)})
row.names(country_data)[80] <- "temp_mean"

# Trend over the period

country_data[81,] <- apply(country_data[57:77,], 2, function(x){summary(lm(x~c(1996:2016)))$coef[2,1]})/country_data[57,]
row.names(country_data)[81] <- "d_temp"

```

### High input cover data

```{r}

# Load data
# from https://www.fao.org/faostat/en/#data/LC

# Utilised Agricultural Area

uaa_country <- read.csv("raw_data/FAOSTAT_data_UAA.csv", header = T)
uaa_country$Area <- as.character(uaa_country$Area)
uaa_country[uaa_country=="United Kingdom of Great Britain and Northern Ireland"] <- "UK"
uaa_country[uaa_country=="Czechia"] <- "Czech Republic"
uaa_country <- droplevels(uaa_country[uaa_country$Area %in% country_name,c("Area","Year","Value")])
uaa_country$Value <- 1000*uaa_country$Value
uaa_country <- dcast(uaa_country[uaa_country$Year>1999,], Year~Area, fun.aggregate=sum, value.var="Value")
row.names(uaa_country) <- paste0("uaa",sep="_", uaa_country$Year)
uaa_country$Year <- NULL
uaa_country[21,] <- apply(uaa_country[1:17,], 2, function(x){mean(x, na.rm=T)})
row.names(uaa_country)[21] <- "uaa_mean"
country_data <- rbind(country_data, uaa_country)
 

# High input farm cover

input_country <- read.csv("raw_data/aei_ps_inp_1_Data.csv", header = T)
# from https://appsso.eurostat.ec.europa.eu/nui/submitViewTableAction.do 
# Notice that this dataset has been downloaded in 2020, the current version (2022) contains some slight changes, which lead to opposite trends for Belgium, Greece, Latvia and Spain. Based on previous releases about agriculture intensification in the EU (https://www.eea.europa.eu/publications/eea_report_2005_6 and https://link.springer.com/article/10.1007/s11356-021-17655-4), we focus on the dataset from 2020 ("aei_ps_inp_1_Data.csv") which is more consistant with these releases.
 
input_country <- na.omit(droplevels(input_country[input_country$INDIC_AG=="High-input farms",]))
input_country$Value <- as.character(input_country$Value)
input_country$Value <- gsub(" ","",input_country$Value)
input_country$Value <- as.numeric(input_country$Value)
input_country <- dcast(input_country, TIME~GEO, fun.aggregate=sum, value.var="Value")

names(input_country) <- c("Year","Austria","Belgium","Bulgaria","Croatia","Cyprus","Czech Republic","Denmark","Estonia","EU","Finland","France","Germany","Greece","Hungary","Ireland","Italy","Latvia","Lithuania","Luxembourg","Malta","Netherlands","Poland","Portugal","Romania","Slovakia","Slovenia","Spain","Sweden","UK")

input_country[input_country==0] <- NA
 
row.names(input_country) <- paste0("hico",sep="_",input_country$Year)
 
input_country$Year <- input_country$EU <- NULL
 
input_country[11,] <- apply(input_country[1:10,], 2, function(x){mean(x, na.rm=T)})
 
row.names(input_country)[11] <- "hico_mean"

input_country[12,] <- apply(input_country[1:10,], 2, function(x){summary(lm(x~c(2007:2016)))$coef[2,1]})/apply(input_country[1:3,], 2,function(x){mean(x, na.rm=T)})
 
row.names(input_country)[12] <- "d_hico"


input_country <- input_country[,sort(names(input_country))]
 
input_country <- data.frame(input_country[,c(1:13)], Iceland=0, input_country[,c(14:20)], Norway=0, input_country[,c(21:27)], Switzerland=0, UK=input_country$UK)
 
names(input_country)[6] <- "Czech Republic"

country_data <- rbind(country_data, input_country) 

country_data[115,] <- country_data[113,]/country_data[102,]
row.names(country_data)[115] <- "hico_mean_perc"

```

### Forest cover data

```{r}

# Load data

forest_country <- read.csv("output/fao_data_landuse.csv", header = T)
forest_country <- forest_country[forest_country$Item!="Country area", c("Area","Item","Year","Value")]
forest_country$Value <- 10*forest_country$Value
forest_country <- dcast(forest_country, Item+Year~Area, fun.aggregate=sum, value.var="Value")

# Mean value over the period

forest_country$Item <- as.character(forest_country$Item)
forest_country[forest_country=="Forest land"] <- "for"
forest_country[forest_country=="Naturally regenerating forest"] <- "nat"
forest_country[forest_country=="Planted Forest"] <- "pla"
forest_country[forest_country==0] <- NA
row.names(forest_country) <- paste0(forest_country$Item,sep="_",forest_country$Year)
forest_country$Year <- forest_country$Item <- NULL
forest_country <- as.data.frame(t(apply(forest_country, 1, function(x){x/area_country$Value})))
forest_country[91,] <- apply(forest_country[1:27,], 2, function(x){mean(x, na.rm=T)})
row.names(forest_country)[91] <- "for_mean"
forest_country[92,] <- apply(forest_country[31:57,], 2, function(x){mean(x, na.rm=T)})
row.names(forest_country)[92] <- "nat_mean"
forest_country[93,] <- apply(forest_country[61:87,], 2, function(x){mean(x, na.rm=T)})
row.names(forest_country)[93] <- "pla_mean"

# Trend over the period

forest_country[94,] <- apply(forest_country[1:27,], 2, function(x){summary(lm(x~c(1990:2016)))$coef[2,1]})/forest_country[11,]
row.names(forest_country)[94] <- "d_for"
for_sig_trend <- apply(forest_country[1:27,], 2, function(x){summary(lm(x~c(1990:2016)))$coef[2,4]})
forest_country[94,which(for_sig_trend>0.05)] <- 0

# Merge data

country_data <- rbind(country_data, forest_country)

# Save data (you can directly use this file to skip the previous steps)

write.csv(country_data,"output/country_data.csv", row.names = F)
# country_data <- read.csv("output/country_data.csv", header = T, row.names=1)

# Clean up the final dataset

country_data2 <- as.data.frame(t(country_data))
country_data2$country <- row.names(country_data2)
country_data2$country[country_data2$country=="UK"]<-"United Kingdom"
country_data2[country_data2==0] <- NA
country_data2$country <- as.factor(country_data2$country)

```

## Plot pressures and species trends

### Map pressures

```{r}

library(sf)
library(rnaturalearth)

worldmap <- ne_countries(scale = 'medium', type = 'countries',returnclass = 'sf')

europe_cropped <- st_crop(worldmap[worldmap$sovereignt %in% c("Austria","Belgium","Bulgaria","Cyprus","Czech Republic","Denmark",
                                                           "Estonia","Finland","France","Germany","Greece","Hungary","Ireland",
                                                           "Italy","Latvia","Lithuania","Netherlands","Norway","Portugal","Poland","Romania",
                                                           "Slovakia","Spain","Sweden","Switzerland","United Kingdom",
                                                           "Croatia","Republic of Serbia","Albania","Slovenia","Bosnia and Herzegovina","Kosovo",
                                                           "Montenegro", "Belarus","Ukraine","Russia","Moldova","Macedonia","Luxembourg"),],
                          xmin = -12, xmax = 35,ymin = 30, ymax = 73)
                          
europe_cropped$fill_param <- rep(NA,nrow(europe_cropped))

europe_cropped$fill_param[europe_cropped$sovereignt %in% c("Austria","Belgium","Bulgaria","Cyprus","Czech Republic","Denmark",
                                                          "Estonia","Finland","France","Germany","Greece","Hungary","Ireland",
                                                          "Italy","Latvia","Lithuania","Netherlands","Norway","Portugal","Poland","Romania",
                                                          "Slovakia","Spain","Sweden","Switzerland","United Kingdom","Slovenia","Luxembourg")] <- "PECBMS member (in 2016)"
                                                          

europe_cropped$Urbanisation <- country_data2$urb_mean[match(europe_cropped$sovereignt,country_data2$country)]

europe_cropped$High_input_farm_cover <- country_data2$hico_mean_perc[match(europe_cropped$sovereignt,country_data2$country)]

europe_cropped$Temperature <- country_data2$temp_mean[match(europe_cropped$sovereignt,country_data2$country)]

europe_cropped$Forest_cover <- country_data2$for_mean[match(europe_cropped$sovereignt,country_data2$country)]

europe_cropped[europe_cropped$sovereignt=="Croatia",c("Urbanisation","High_input_farm_cover","Temperature","Forest_cover")] <- NA
                                                        
```

### Map species trends

```{r}

# Slope by country for farmland species
                                                        
country_bird_trend_agri <- data.frame(year=NA,variable=NA,value=NA, slope=NA)

for(ii in levels(df_pop$CountryGroup)){

  df_agri_1 <- droplevels(df_pop[df_pop$Species %in% pecbms_hab$Species[pecbms_hab$Habitat=="Farmland"] & df_pop$CountryGroup==ii,])
  
  df_agri_1b <- as.data.frame(df_agri_1 %>% group_by(Species) %>% summarize(sum_ab=sum(Abd)))
  
  df_agri_1 <- as.data.frame(df_agri_1[!(df_agri_1$Species %in% levels(droplevels(df_agri_1b$Species[df_agri_1b$sum_ab==0]))),c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")])
  
  if(ii=="Netherlands"){
    df_agri_1 <- as.data.frame(df_agri_1[!(df_agri_1$Species %in% c("Galerida cristata","Emberiza calandra",levels(droplevels(df_agri_1b$Species[df_agri_1b$sum_ab==0])))),c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")]) 
  }
  
  if(ii=="Czech Republic"){
    df_agri_1 <- as.data.frame(df_agri_1[!(df_agri_1$Species %in% c("Motacilla flava",levels(droplevels(df_agri_1b$Species[df_agri_1b$sum_ab==0])))),c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")]) 
  }
  
  if(ii=="Latvia"){df_agri_1 <- droplevels(df_agri_1[df_agri_1$Year>2004,])}
  if(ii=="Lithuania"){df_agri_1 <- droplevels(df_agri_1[df_agri_1$Year>2010,])}
  if(ii=="Sweden"){df_agri_1 <- droplevels(df_agri_1[df_agri_1$Year>1990,])}
  if(ii=="Switzerland"){df_agri_1 <- droplevels(df_agri_1[df_agri_1$Year>1999,])}
  if(ii=="Germany"){df_agri_1 <- droplevels(df_agri_1[df_agri_1$Year>1990,])}
  
  names(df_agri_1) <- c("Species","Code","Year","Abundance","Abundance_SE","Index_imputed","Index_imputed_SE")
  msi_agri_1 <- msi_fun3b(df_agri_1,ref_year = min(df_agri_1$Year), niter = niter, ref_value = "Abundance")
  country_bird_trend_agri2 <- data.frame(year=c(min(df_agri_1$Year):max(df_agri_1$Year)),variable=ii,
                                    value=msi_agri_1$msi$mean_msi_final,slope=msi_agri_1$coef$slope/msi_agri_1$msi$mean_msi_final[1])
                                    
  country_bird_trend_agri2$slope <- summary(lm(value~year,country_bird_trend_agri2))$coef[2,1]/msi_agri_1$msi$mean_msi_final[1]
  
  if(ii=="Czech Republic"){country_bird_trend_agri2$variable <- "Czech.Republic"}
  if(ii=="Republic of Ireland"){country_bird_trend_agri2$variable <- "Ireland"}
  if(ii=="United Kingdom"){country_bird_trend_agri2$variable <- "UK"}
  country_bird_trend_agri <- rbind(country_bird_trend_agri,country_bird_trend_agri2)
}


# Slope by country for woodland species

country_bird_trend_forest <- data.frame(year=NA,variable=NA,value=NA, slope=NA)

for(ii in levels(df_pop$CountryGroup)){

  df_forest_1 <-  droplevels(df_pop[df_pop$Species %in% pecbms_hab$Species[pecbms_hab$Habitat=="Forest"] & df_pop$CountryGroup==ii,])
  
  df_forest_1b <- as.data.frame(df_forest_1 %>% group_by(Species) %>% summarize(sum_ab=sum(Abd), count=n()))
  
  df_forest_1 <- as.data.frame(df_forest_1[!(df_forest_1$Species %in% levels(droplevels(df_forest_1b$Species[df_forest_1b$sum_ab==0]))) & df_forest_1$Year %in% as.numeric(levels(as.factor(as.character(df_forest_1$Year[df_forest_1$Species %in% droplevels(df_forest_1b$Species[which(df_forest_1b$count==max(df_forest_1b$count))])])))),c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")])
  
  if(ii=="Latvia"){df_forest_1 <- droplevels(df_forest_1[df_forest_1$Year>2004,])}
  if(ii=="Lithuania"){df_forest_1 <- droplevels(df_forest_1[df_forest_1$Year>2010,])}
  if(ii=="Sweden"){df_forest_1 <- droplevels(df_forest_1[df_forest_1$Year>1990,])}
  if(ii=="Switzerland"){df_forest_1 <- droplevels(df_forest_1[df_forest_1$Year>1999,])}
  if(ii=="Germany"){df_forest_1 <- droplevels(df_forest_1[df_forest_1$Year>1990,])}
  
  names(df_forest_1) <- c("Species","Code","Year","Abundance","Abundance_SE","Index_imputed","Index_imputed_SE")
  msi_forest_1 <- msi_fun3b(df_forest_1,ref_year = min(df_forest_1$Year), niter = niter, ref_value = "Abundance")
  
  country_bird_trend_forest2 <- data.frame(year=c(min(df_forest_1$Year):max(df_forest_1$Year)),variable=ii,value=msi_forest_1$msi$mean_msi_final,slope=msi_forest_1$coef$slope/msi_forest_1$msi$mean_msi_final[1])
  
  country_bird_trend_forest2$slope <- summary(lm(value~year,country_bird_trend_forest2))$coef[2,1]/msi_forest_1$msi$mean_msi_final[1]
  
  if(ii=="Czech Republic"){country_bird_trend_forest2$variable <- "Czech.Republic"}
  if(ii=="Republic of Ireland"){country_bird_trend_forest2$variable <- "Ireland"}
  if(ii=="United Kingdom"){country_bird_trend_forest2$variable <- "UK"}
  
  country_bird_trend_forest <- rbind(country_bird_trend_forest,country_bird_trend_forest2)
}


# Slope by country for urban species

country_bird_trend_build <- data.frame(year=NA,variable=NA,value=NA, slope=NA)

for(ii in levels(df_pop$CountryGroup)){

  df_build_1 <- droplevels(df_pop[df_pop$Species %in% levels(droplevels(eunis_hab$speciesname[eunis_hab$season=="B" & eunis_hab$codeeco=="urban"])) & df_pop$CountryGroup==ii,])
  
  df_build_1b <- as.data.frame(df_build_1 %>% group_by(Species) %>% summarize(sum_ab=sum(Abd), count=n()))
  
  df_build_1 <- as.data.frame(df_build_1[!(df_build_1$Species %in% levels(droplevels(df_build_1b$Species[df_build_1b$sum_ab==0]))) & df_build_1$Year %in% as.numeric(levels(as.factor(as.character(df_build_1$Year[df_build_1$Species %in% droplevels(df_build_1b$Species[which(df_build_1b$count==max(df_build_1b$count))])])))),c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")])
  
    if(ii=="Netherlands"){
    df_build_1 <- as.data.frame(df_build_1[!(df_build_1$Species %in% c("Galerida cristata","Emberiza calandra",levels(droplevels(df_build_1b$Species[df_build_1b$sum_ab==0])))),c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")]) 
    }
    
    if(ii=="Czech Republic"){
    df_build_1 <- as.data.frame(df_build_1[!(df_build_1$Species %in% c("Motacilla flava",levels(droplevels(df_build_1b$Species[df_build_1b$sum_ab==0])))),c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")]) 
    }
    
  if(ii=="Latvia"){df_build_1 <- droplevels(df_build_1[df_build_1$Year>2004,])}
  if(ii=="Lithuania"){df_build_1 <- droplevels(df_build_1[df_build_1$Year>2010,])}
  if(ii=="Sweden"){df_build_1 <- droplevels(df_build_1[df_build_1$Year>1990,])}
  if(ii=="Switzerland"){df_build_1 <- droplevels(df_build_1[df_build_1$Year>1999,])}
  if(ii=="Germany"){df_build_1 <- droplevels(df_build_1[df_build_1$Year>1990,])}
  
  names(df_build_1) <- c("Species","Code","Year","Abundance","Abundance_SE","Index_imputed","Index_imputed_SE")
  msi_build_1 <- msi_fun3b(df_build_1,ref_year = min(df_build_1$Year), niter = niter, ref_value = "Abundance")
  
  country_bird_trend_build2 <- data.frame(year=c(min(df_build_1$Year):max(df_build_1$Year)),variable=ii,value=msi_build_1$msi$mean_msi_final,slope=msi_build_1$coef$slope/msi_build_1$msi$mean_msi_final[1])
  country_bird_trend_build2$slope <- summary(lm(value~year,country_bird_trend_build2))$coef[2,1]/msi_build_1$msi$mean_msi_final[1]
  
  if(ii=="Czech Republic"){country_bird_trend_build2$variable <- "Czech.Republic"}
  if(ii=="Republic of Ireland"){country_bird_trend_build2$variable <- "Ireland"}
  if(ii=="United Kingdom"){country_bird_trend_build2$variable <- "UK"}
  country_bird_trend_build <- rbind(country_bird_trend_build,country_bird_trend_build2)
}


# Slope by country for hot dwellers

country_bird_trend_temp_warm <- data.frame(year=NA,variable=NA,value=NA,slope=NA)

for(ii in levels(df_pop$CountryGroup)){

  df_temp_warm_1 <- droplevels(df_pop[df_pop$Species %in% as.character(sti$SPECIES[sti$STI>quantile(sti$STI[sti$SPECIES %in% levels(droplevels(df_pop$Species))], 0.7)]) & df_pop$CountryGroup==ii,])
  
  df_temp_warm_1b <- as.data.frame(df_temp_warm_1 %>% group_by(Species) %>% summarize(sum_ab=sum(Abd)))
  
  df_temp_warm_1 <- as.data.frame(df_temp_warm_1[!(df_temp_warm_1$Species %in% levels(droplevels(df_temp_warm_1b$Species[df_temp_warm_1b$sum_ab==0]))),c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")])
  
  df_temp_warm_1b <- as.data.frame(df_temp_warm_1 %>% group_by(Species) %>% summarize(sum_ab=sum(Abd), count=n()))
  
  df_temp_warm_1 <- as.data.frame(df_temp_warm_1[!(df_temp_warm_1$Species %in% levels(droplevels(df_temp_warm_1b$Species[df_temp_warm_1b$sum_ab==0]))) & df_temp_warm_1$Year %in% as.numeric(levels(as.factor(as.character(df_temp_warm_1$Year[df_temp_warm_1$Species %in% droplevels(df_temp_warm_1b$Species[which(df_temp_warm_1b$count==max(df_temp_warm_1b$count))])])))),c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")])
  
  if(ii=="Netherlands"){
    df_temp_warm_1 <- as.data.frame(df_temp_warm_1[!(df_temp_warm_1$Species %in% c("Galerida cristata","Emberiza calandra",levels(droplevels(df_temp_warm_1b$Species[df_temp_warm_1b$sum_ab==0])))),c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")]) 
  }
  
  if(ii=="Latvia"){df_temp_warm_1 <- droplevels(df_temp_warm_1[df_temp_warm_1$Year>2004,])}
  if(ii=="Lithuania"){df_temp_warm_1 <- droplevels(df_temp_warm_1[df_temp_warm_1$Year>2010,])}
  if(ii=="Switzerland"){
    df_temp_warm_1 <- as.data.frame(df_temp_warm_1[df_temp_warm_1$Year>1999,c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")]) 
  }
  if(ii=="Sweden"){
    df_temp_warm_1 <- as.data.frame(df_temp_warm_1[df_temp_warm_1$Year>1990,c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")]) 
  }
  if(ii=="Germany"){
    df_temp_warm_1 <- as.data.frame(df_temp_warm_1[df_temp_warm_1$Year>1990,c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")]) 
  }
  names(df_temp_warm_1) <- c("Species","Code","Year","Abundance","Abundance_SE","Index_imputed","Index_imputed_SE")
  msi_temp_warm_1 <- msi_fun3b(df_temp_warm_1,ref_year = min(df_temp_warm_1$Year), niter = niter, ref_value = "Abundance")
  
  country_bird_trend_temp_warm2 <- data.frame(year=c(min(df_temp_warm_1$Year):max(df_temp_warm_1$Year)),variable=ii,value=msi_temp_warm_1$msi$mean_msi_final,slope=msi_temp_warm_1$coef$slope/msi_temp_warm_1$msi$mean_msi_final[1])
  country_bird_trend_temp_warm2$slope <- summary(lm(value~year,country_bird_trend_temp_warm2))$coef[2,1]/msi_temp_warm_1$msi$mean_msi_final[1]
  
  if(ii=="Czech Republic"){country_bird_trend_temp_warm2$variable <- "Czech.Republic"}
  if(ii=="Republic of Ireland"){country_bird_trend_temp_warm2$variable <- "Ireland"}
  if(ii=="United Kingdom"){country_bird_trend_temp_warm2$variable <- "UK"}
  
  country_bird_trend_temp_warm <- rbind(country_bird_trend_temp_warm,country_bird_trend_temp_warm2)
}


# Slope by country for cold dwellers

country_bird_trend_temp_cold <- data.frame(year=NA,variable=NA,value=NA, slope=NA)

for(ii in levels(df_pop$CountryGroup)[-c(2,3,9)]){

  df_temp_cold_1 <-  droplevels(df_pop[df_pop$Species %in% as.character(sti$SPECIES[sti$STI<quantile(sti$STI[sti$SPECIES %in% levels(droplevels(df_pop$Species))], 0.3)]) & df_pop$CountryGroup==ii,])
  
  df_temp_cold_1b <- as.data.frame(df_temp_cold_1 %>% group_by(Species) %>% summarize(sum_ab=sum(Abd)))
  
  df_temp_cold_1 <- as.data.frame(df_temp_cold_1[!(df_temp_cold_1$Species %in% levels(droplevels(df_temp_cold_1b$Species[df_temp_cold_1b$sum_ab==0]))),c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")])
  
  df_temp_cold_1b <- as.data.frame(df_temp_cold_1 %>% group_by(Species) %>% summarize(sum_ab=sum(Abd), count=n()))
  
  df_temp_cold_1 <- as.data.frame(df_temp_cold_1[!(df_temp_cold_1$Species %in% levels(droplevels(df_temp_cold_1b$Species[df_temp_cold_1b$sum_ab==0]))) & df_temp_cold_1$Year %in% as.numeric(levels(as.factor(as.character(df_temp_cold_1$Year[df_temp_cold_1$Species %in% droplevels(df_temp_cold_1b$Species[which(df_temp_cold_1b$count==max(df_temp_cold_1b$count))])])))),c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")])
  
  if(ii=="Netherlands"){
    df_temp_cold_1 <- as.data.frame(df_temp_cold_1[!(df_temp_cold_1$Species %in% c("Galerida cristata","Emberiza calandra",levels(droplevels(df_temp_cold_1b$Species[df_temp_cold_1b$sum_ab==0])))),c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")]) 
  }
  if(ii=="Latvia"){df_temp_cold_1 <- droplevels(df_temp_cold_1[df_temp_cold_1$Year>2004,])}
  if(ii=="Lithuania"){df_temp_cold_1 <- droplevels(df_temp_cold_1[df_temp_cold_1$Year>2010,])}
  if(ii=="Sweden"){
    df_temp_cold_1 <- as.data.frame(df_temp_cold_1[df_temp_cold_1$Year>1990,c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")]) 
  }
  if(ii=="Switzerland"){
    df_temp_cold_1 <- as.data.frame(df_temp_cold_1[df_temp_cold_1$Year>1999,c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")]) 
  }
  if(ii=="Germany"){
    df_temp_cold_1 <- as.data.frame(df_temp_cold_1[df_temp_cold_1$Year>1990,c("Species","Code","Year","Abd","SE_Abd","Index","Index_SE")]) 
  }
  names(df_temp_cold_1) <- c("Species","Code","Year","Abundance","Abundance_SE","Index_imputed","Index_imputed_SE")
  if(!(ii %in% c(2,3,9))){
    msi_temp_cold_1 <- msi_fun3b(df_temp_cold_1,ref_year = min(df_temp_cold_1$Year), niter = niter, ref_value = "Abundance")
    
    country_bird_trend_temp_cold2 <- data.frame(year=c(min(df_temp_cold_1$Year):max(df_temp_cold_1$Year)),variable=ii,value=msi_temp_cold_1$msi$mean_msi_final,slope=msi_temp_cold_1$coef$slope/msi_temp_cold_1$msi$mean_msi_final[1])
    
    country_bird_trend_temp_cold2$slope <- summary(lm(value~year,country_bird_trend_temp_cold2))$coef[2,1]/msi_temp_cold_1$msi$mean_msi_final[1]
    
  }else{country_bird_trend_temp_cold2 <- data.frame(year=NA,variable=NA,value=NA, slope=NA)
  }
 
  if(ii=="Czech Republic"){country_bird_trend_temp_cold2$variable <- "Czech.Republic"}
  if(ii=="Republic of Ireland"){country_bird_trend_temp_cold2$variable <- "Ireland"}
  if(ii=="United Kingdom"){country_bird_trend_temp_cold2$variable <- "UK"}
  
  country_bird_trend_temp_cold <- rbind(country_bird_trend_temp_cold,country_bird_trend_temp_cold2)
}


# Merge hot and cold dwellers

country_bird_trend_temp <- merge(country_bird_trend_temp_warm,country_bird_trend_temp_cold, by=c("year","variable"), all=T)
names(country_bird_trend_temp)[3:6]<-c("value_warm","slope_warm","value_cold","slope_cold")


# Summarise data

country_bird_trend_agrib <- data.frame(droplevels(country_bird_trend_agri) %>% group_by(variable) %>% summarize(slope_agri2=mean(slope)))
country_bird_trend_buildb <- data.frame(droplevels(country_bird_trend_build) %>% group_by(variable) %>% summarize(slope_build2=mean(slope)))
country_bird_trend_forestb <- data.frame(droplevels(country_bird_trend_forest) %>% group_by(variable) %>% summarize(slope_forest2=mean(slope)))
country_bird_trend_tempb <- data.frame(droplevels(country_bird_trend_temp) %>% group_by(variable) %>% summarize(slope_warm2=mean(slope_warm),slope_cold2=mean(slope_cold)))


# Update country names

country_bird_trend_agrib$variable[country_bird_trend_agrib$variable=="UK"] <- country_bird_trend_buildb$variable[country_bird_trend_buildb$variable=="UK"] <-
  country_bird_trend_forestb$variable[country_bird_trend_forestb$variable=="UK"] <- country_bird_trend_tempb$variable[country_bird_trend_tempb$variable=="UK"] <- "United Kingdom"
country_bird_trend_agrib$variable[country_bird_trend_agrib$variable=="Czech.Republic"] <- country_bird_trend_buildb$variable[country_bird_trend_buildb$variable=="Czech.Republic"] <-
  country_bird_trend_forestb$variable[country_bird_trend_forestb$variable=="Czech.Republic"] <- country_bird_trend_tempb$variable[country_bird_trend_tempb$variable=="Czech.Republic"] <- "Czech Republic"


# Add slopes to geographical data

europe_cropped$slope_agri <- country_bird_trend_agrib$slope_agri2[match(europe_cropped$sovereignt,country_bird_trend_agrib$variable)]

europe_cropped$slope_build <- country_bird_trend_buildb$slope_build2[match(europe_cropped$sovereignt,country_bird_trend_buildb$variable)]
europe_cropped$slope_build[europe_cropped$slope_build>0.05] <- 0.05

europe_cropped$slope_forest <- country_bird_trend_forestb$slope_forest2[match(europe_cropped$sovereignt,country_bird_trend_forestb$variable)]
europe_cropped$slope_forest[europe_cropped$slope_forest>0.05] <- 0.05

europe_cropped$slope_warm <- country_bird_trend_tempb$slope_warm2[match(europe_cropped$sovereignt,country_bird_trend_tempb$variable)]
europe_cropped$slope_warm[europe_cropped$slope_warm>0.07] <- 0.07

europe_cropped$slope_cold <- country_bird_trend_tempb$slope_cold2[match(europe_cropped$sovereignt,country_bird_trend_tempb$variable)]
europe_cropped$slope_cold[europe_cropped$slope_cold>0.07] <- 0.07

# Reproject data

europe_map <- sf::st_transform(
  europe_cropped,
  "+init=epsg:27572"
)

# Simplify map

library(rmapshaper)
europe_map_simpl <- ms_simplify(europe_map, keep = 0.07,
                                keep_shapes = FALSE)
                                
                                
# Map of high input farm cover

pres1 <- ggplot() + geom_sf(data = europe_map_simpl, aes(fill = High_input_farm_cover)) + scale_fill_gradient(low="white",high="#D302F9",na.value="lightgrey")+
  theme_void()+ coord_sf(datum = NA)+
  annotation_custom(grob = gtable_filter(ggplot_gtable(ggplot_build(ggplot() + geom_sf(data = europe_map_simpl, aes(fill = High_input_farm_cover)) + scale_fill_gradient(low="white",high="#D302F9",na.value="lightgrey")+
                                                                      theme_void()+  theme(legend.title=element_blank())+ coord_sf(datum = NA))), "guide-box") , xmin = 0, xmax = 432734.7, ymin = 4000000, ymax = 5000000) +theme(legend.position = "none")


# Map of urbanisation
                                                                      
pres2 <- ggplot() + geom_sf(data = europe_map_simpl, aes(fill = Urbanisation)) + scale_fill_gradient(low="white",high="#196DF6",na.value="lightgrey")+
  theme_void()+ coord_sf(datum = NA)+
  annotation_custom(grob = gtable_filter(ggplot_gtable(ggplot_build(ggplot() + geom_sf(data = europe_map_simpl, aes(fill = Urbanisation)) +scale_fill_gradient(low="white",high="#196DF6",na.value="lightgrey")+
                                                                      theme_void()+  theme(legend.title=element_blank())+ coord_sf(datum = NA))), "guide-box") , xmin = 0, xmax = 432734.7, ymin = 4000000, ymax = 5000000) +theme(legend.position = "none")


# Map of temperature
                                                                      
pres3 <- ggplot() + geom_sf(data = europe_map_simpl, aes(fill = Temperature)) + scale_fill_gradient(low="white",high="#FA0900",na.value="lightgrey")+
  theme_void()+ coord_sf(datum = NA)+
  annotation_custom(grob = gtable_filter(ggplot_gtable(ggplot_build(ggplot() + geom_sf(data = europe_map_simpl, aes(fill = Temperature)) + scale_fill_gradient(low="white",high="#FA0900",na.value="lightgrey")+
                                                                      theme_void()+  theme(legend.title=element_blank())+ coord_sf(datum = NA))), "guide-box") , xmin = 0, xmax = 432734.7, ymin = 4000000, ymax = 5000000) +theme(legend.position = "none")


# Map of forest cover

pres4 <- ggplot() + geom_sf(data = europe_map_simpl, aes(fill = Forest_cover)) + scale_fill_gradient(low="white",high="#1BAE20",na.value="lightgrey")+
  theme_void()+ coord_sf(datum = NA)+
  annotation_custom(grob = gtable_filter(ggplot_gtable(ggplot_build(ggplot() + geom_sf(data = europe_map_simpl, aes(fill = Forest_cover)) + scale_fill_gradient(low="white",high="#1BAE20",na.value="lightgrey")+
                                                                      theme_void()+  theme(legend.title=element_blank())+ coord_sf(datum = NA))), "guide-box") , xmin = 0, xmax = 432734.7, ymin = 4000000, ymax = 5000000) +theme(legend.position = "none")


# Map of farmland bird trends

pres_tend_agri <- ggplot() + geom_sf(data = europe_map_simpl, aes(fill = slope_agri)) + scale_fill_gradient2(low="#FA0900",high="#1BAE20", mid="white",midpoint = 0,na.value="lightgrey")+theme_void()+ coord_sf(datum = NA)+
  annotation_custom(grob = gtable_filter(ggplot_gtable(ggplot_build(ggplot() + geom_sf(data = europe_map_simpl, aes(fill = slope_agri)) + scale_fill_gradient2(low="#FA0900",high="#1BAE20",mid="white",midpoint = 0,na.value="white")+#scale_fill_viridis(option = "C",na.value="lightgrey")+
                                                                      theme_void()+  theme(legend.title=element_blank())+ coord_sf(datum = NA))), "guide-box") , xmin = 0, xmax = 432734.7, ymin = 4000000, ymax = 5000000) +theme(legend.position = "none")


# Map of woodland bird trends

pres_tend_forest <- ggplot() + geom_sf(data = europe_map_simpl, aes(fill = slope_forest)) + scale_fill_gradient2(low="#FA0900",high="#1BAE20", mid="white",midpoint = 0,na.value="lightgrey")+
  theme_void()+ coord_sf(datum = NA)+
  annotation_custom(grob = gtable_filter(ggplot_gtable(ggplot_build(ggplot() + geom_sf(data = europe_map_simpl, aes(fill = slope_forest)) + scale_fill_gradient2(low="#FA0900",high="#1BAE20",mid="white",midpoint = 0,na.value="white")+                                                                      theme_void()+  theme(legend.title=element_blank())+ coord_sf(datum = NA))), "guide-box") , xmin = 0, xmax = 432734.7, ymin = 4000000, ymax = 5000000) +theme(legend.position = "none")


# Map of urban bird trends

pres_tend_build <- ggplot() + geom_sf(data = europe_map_simpl, aes(fill = slope_build)) + scale_fill_gradient2(low="#FA0900",high="#1BAE20", mid="white",midpoint = 0,na.value="lightgrey")+
  theme_void()+ coord_sf(datum = NA)+
  annotation_custom(grob = gtable_filter(ggplot_gtable(ggplot_build(ggplot() + geom_sf(data = europe_map_simpl, aes(fill = slope_build)) + scale_fill_gradient2(low="#FA0900",high="#1BAE20",mid="white",midpoint = 0,na.value="lightgrey")+                                                                      theme_void()+  theme(legend.title=element_blank())+ coord_sf(datum = NA))), "guide-box") , xmin = 0, xmax = 432734.7, ymin = 4000000, ymax = 5000000) +theme(legend.position = "none")


# Map of hot dweller trends

pres_tend_warm <- ggplot() + geom_sf(data = europe_map_simpl, aes(fill = slope_warm)) + scale_fill_gradient2(low="#FA0900",high="#1BAE20", mid="white",midpoint = 0,na.value="lightgrey")+
  theme_void()+ coord_sf(datum = NA)+
  annotation_custom(grob = gtable_filter(ggplot_gtable(ggplot_build(ggplot() + geom_sf(data = europe_map_simpl, aes(fill = slope_warm)) + scale_fill_gradient2(low="#FA0900",high="#1BAE20",mid="white",midpoint = 0,na.value="lightgrey")+                                                                      theme_void()+  theme(legend.title=element_blank())+ coord_sf(datum = NA))), "guide-box") , xmin = 0, xmax = 432734.7, ymin = 4000000, ymax = 5000000) +theme(legend.position = "none")

# Map of cold dweller trends

pres_tend_cold <- ggplot() + geom_sf(data = europe_map_simpl, aes(fill = slope_cold)) + scale_fill_gradient2(low="#FA0900",high="#1BAE20", mid="white",midpoint = 0,na.value="lightgrey")+
  theme_void()+ coord_sf(datum = NA)+
  annotation_custom(grob = gtable_filter(ggplot_gtable(ggplot_build(ggplot() + geom_sf(data = europe_map_simpl, aes(fill = slope_cold)) + scale_fill_gradient2(low="#FA0900",high="#1BAE20",mid="white",midpoint = 0,na.value="lightgrey")+                                                                      theme_void()+  theme(legend.title=element_blank())+ coord_sf(datum = NA))), "guide-box") , xmin = 0, xmax = 432734.7, ymin = 4000000, ymax = 5000000) +theme(legend.position = "none")
  
pres_tend_temp <- ggplot() + geom_sf(data = europe_map_simpl, fill="transparent")+theme_void()+ coord_sf(datum = NA)

```

### Plot trends
```{r}

# Compute coordinate centroids of each country

require("pacman") 
p_load(ggplot2, ggtree, dplyr, tidyr, sp, maps, pipeR, grid, XML, gtable)

getLabelPoint <- function(county) {Polygon(county[c('long', 'lat')])@labpt}
map_base_data <- droplevels(subset(map_data("world"), region %in% c("Austria","Belgium","Bulgaria","Cyprus","Czech Republic","Denmark", #2000
                                                                  "Estonia","Finland","France","Germany","Greece","Hungary","Italy",
                                                                  "Latvia","Lithuania","Luxembourg","Netherlands","Norway",
                                                                  "Poland","Portugal","Ireland","Romania","Slovakia","Slovenia","Spain","Sweden","Switzerland","UK")))
                                                                  
map_base_data2 <- droplevels(map_base_data[which(is.na(map_base_data$subregion) | map_base_data$region=="UK"),])

centroids <- by(map_base_data2, map_base_data2$region, getLabelPoint)     # Returns list

centroids <- do.call("rbind.data.frame", centroids)  # Convert to Data Frame

names(centroids) <- c('long', 'lat')                 # Appropriate Header

rownames(centroids) <- c("Austria","Belgium","Bulgaria","Cyprus","Czech.Republic","Denmark","Estonia","Finland","France","Germany",
                       "Greece","Hungary","Ireland","Italy","Latvia","Lithuania","Luxembourg","Netherlands","Norway",
                       "Poland","Portugal","Romania","Slovakia","Slovenia","Spain","Sweden","Switzerland","UK")

# Get trend by country for high input farm cover

data_pres_wide <- melt(data.frame(year=sub(".*_","",row.names(input_country[2:10,])),apply(rbind(input_country[1:10,],country_data[102,]),2,function(x){y <- x[2:10]/x[11]})))

names(data_pres_wide)[2] <- "Country"
data_pres_long <- add_rownames(centroids, "Country") %>% left_join(data_pres_wide) %>% split(., .$Country)

data_pres_long$Switzerland <- data_pres_long$Norway <- NULL

graph_hico <- setNames(lapply(1:length(data_pres_long), function(i){
  test <- data_pres_long[[i]]
  
  bound <- max(test$value)-min(test$value)
  med_val <- bound/2+min(test$value)
  if(bound<0.10){ylimit <- c((med_val-0.05),(med_val+0.05))
  }else{ylimit <- c(min(test$value),max(test$value))}
  
  ggplot(na.omit(test), aes(x=year, y=value, group=Country)) + # ici
      geom_line(col="black" ,size=1.5, alpha=0.5)+
      xlab(NULL) + 
      ylab(NULL) + 
      theme_modern() + theme_transparent()+
      scale_y_continuous(limits = ylimit) +
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)}), names(data_pres_long))
            
            
# Get trend by country for urbanisation

data_pres_wide <- melt(data.frame(year=sub(".*_","",row.names(urban_country[1:8,])),urban_country[1:8,]))

names(data_pres_wide)[2] <- "Country"
data_pres_long <- add_rownames(centroids, "Country") %>% left_join(data_pres_wide) %>% split(., .$Country)

data_pres_long$Switzerland <- data_pres_long$Norway <- NULL
data_pres_long$Lithuania$value[8] <- data_pres_long$Luxembourg$value[8] <- NA

graph_urb <- setNames(lapply(1:length(data_pres_long), function(i){
  test <- data_pres_long[[i]]
  
  bound <- max(test$value, na.rm=T)-min(test$value,na.rm=T)
  med_val <- bound/2+min(test$value, na.rm=T)
  if(bound<1){ylimit <- c((med_val-0.5),(med_val+0.5))
  }else{ylimit <- c(min(test$value, na.rm=T),max(test$value, na.rm=T))}
  
  ggplot(na.omit(test), aes(x=year, y=value, group=Country)) + # ici
      geom_line(col="black" ,size=1.5, alpha=0.5)+
      xlab(NULL) + 
      ylab(NULL) + 
      theme_modern() + theme_transparent()+
      scale_y_continuous(limits = ylimit) +
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)}), names(data_pres_long))
            
            
# Get trend by country for forest

data_pres_wide <- melt(data.frame(year=sub(".*_","",row.names(forest_country[7:27,])),forest_country[7:27,]))

names(data_pres_wide)[2] <- "Country"
data_pres_long <- add_rownames(centroids, "Country") %>% left_join(data_pres_wide) %>% split(., .$Country)

graph_forest <- setNames(lapply(1:length(data_pres_long), function(i){
  test <- data_pres_long[[i]]
  
  bound <- max(test$value, na.rm=T)-min(test$value,na.rm=T)
  med_val <- bound/2+min(test$value, na.rm=T)
  if(bound<0.04){ylimit <- c((med_val-0.02),(med_val+0.02))
  }else{ylimit <- c(min(test$value, na.rm=T),max(test$value, na.rm=T))}
  
  ggplot(na.omit(test), aes(x=year, y=value, group=Country)) + # ici
      geom_line(col="black" ,size=1.5, alpha=0.5)+
      xlab(NULL) + 
      ylab(NULL) + 
      theme_modern() + theme_transparent()+
      scale_y_continuous(limits = ylimit) +
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)}), names(data_pres_long))
                        
                        
# Get trend by country for temperature

data_pres_wide <- melt(data.frame(year=sub(".*_","",row.names(country_data[57:77,])),country_data[57:77,]))

names(data_pres_wide)[2] <- "Country"
data_pres_long <- add_rownames(centroids, "Country") %>% left_join(data_pres_wide) %>% split(., .$Country)

graph_temp <- setNames(lapply(1:length(data_pres_long), function(i){
  test <- data_pres_long[[i]]
  
  bound <- max(test$value, na.rm=T)-min(test$value,na.rm=T)
  med_val <- bound/2+min(test$value, na.rm=T)
  if(bound<4){ylimit <- c((med_val-2),(med_val+2))
  }else{ylimit <- c(min(test$value, na.rm=T),max(test$value, na.rm=T))}
  
  ggplot(na.omit(test), aes(x=year, y=value, group=Country)) + # ici
      geom_line(col="black" ,size=1.5, alpha=0.5)+
      xlab(NULL) + 
      ylab(NULL) + 
      theme_modern() + theme_transparent()+
      scale_y_continuous(limits = ylimit) +
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)}), names(data_pres_long))
                                    

# Get trend by country for farmland species

data_pres_wide <- na.omit(country_bird_trend_agri[,c("year","variable","value")])

names(data_pres_wide)[2] <- "Country"
data_pres_long <- add_rownames(centroids, "Country") %>% left_join(data_pres_wide) %>% split(., .$Country)

graph_agri <- setNames(lapply(1:length(data_pres_long), function(i){
  test <- data_pres_long[[i]]
  
  test <- test[test$year %in% c(1996:2016),]
  
  test$value <- test$value/test$value[1]
  
  bound <- max(test$value)-min(test$value)
  med_val <- bound/2+min(test$value)
  if(bound<0.30){ylimit <- c((med_val-0.15),(med_val+0.15))
  }else{ylimit <- c(min(test$value),max(test$value))}
  
  ggplot(na.omit(test), aes(x=year, y=value, group=Country)) + # ici
      geom_line(col="black" ,size=1.5, alpha=0.5)+
      xlab(NULL) + 
      ylab(NULL) + 
      theme_modern() + theme_transparent()+
      scale_y_continuous(limits = ylimit) +
      scale_x_continuous(limits = c(1996,2016)) +
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)}), names(data_pres_long))


# Get trend by country for forest species

data_pres_wide <- na.omit(country_bird_trend_forest[,c("year","variable","value")])

names(data_pres_wide)[2] <- "Country"
data_pres_long <- add_rownames(centroids, "Country") %>% left_join(data_pres_wide) %>% split(., .$Country)

graph_wood <- setNames(lapply(1:length(data_pres_long), function(i){
  test <- data_pres_long[[i]]
  
  test <- test[test$year %in% c(1996:2016),]
  
  test$value <- test$value/test$value[1]
  
  bound <- max(test$value)-min(test$value)
  med_val <- bound/2+min(test$value)
  if(bound<0.30){ylimit <- c((med_val-0.15),(med_val+0.15))
  }else{ylimit <- c(min(test$value),max(test$value))}
  
  ggplot(na.omit(test), aes(x=year, y=value, group=Country)) + # ici
      geom_line(col="black" ,size=1.5, alpha=0.5)+
      xlab(NULL) + 
      ylab(NULL) + 
      theme_modern() + theme_transparent()+
      scale_y_continuous(limits = ylimit) +
      scale_x_continuous(limits = c(1996,2016)) +
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)}), names(data_pres_long))


# Get trend by country for urban species

data_pres_wide <- na.omit(country_bird_trend_build[,c("year","variable","value")])

names(data_pres_wide)[2] <- "Country"
data_pres_long <- add_rownames(centroids, "Country") %>% left_join(data_pres_wide) %>% split(., .$Country)

graph_build <- setNames(lapply(1:length(data_pres_long), function(i){
  test <- data_pres_long[[i]]
  
  test <- test[test$year %in% c(1996:2016),]
  
  test$value <- test$value/test$value[1]
  
  bound <- max(test$value)-min(test$value)
  med_val <- bound/2+min(test$value)
  if(bound<0.30){ylimit <- c((med_val-0.15),(med_val+0.15))
  }else{ylimit <- c(min(test$value),max(test$value))}
  
  ggplot(na.omit(test), aes(x=year, y=value, group=Country)) + # ici
      geom_line(col="black" ,size=1.5, alpha=0.5)+
      xlab(NULL) + 
      ylab(NULL) + 
      theme_modern() + theme_transparent()+
      scale_y_continuous(limits = ylimit) +
      scale_x_continuous(limits = c(1996,2016)) +
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)}), names(data_pres_long))


# Get trend by country for hot and cold species

data_pres_wide <- country_bird_trend_temp[,c("year","variable","value_warm","value_cold")]

names(data_pres_wide)[2] <- "Country"
data_pres_long <- add_rownames(centroids, "Country") %>% left_join(data_pres_wide) %>% split(., .$Country)

graph_warm_cold <- setNames(lapply(1:length(data_pres_long), function(i){
  test <- data_pres_long[[i]]
  
  test <- test[test$year %in% c(1996:2016),]
  
  test$value_warm <- test$value_warm/test$value_warm[1]
  test$value_cold <- test$value_cold/test$value_cold[1]
  
  if(test$Country[1]=="Sweden"){test$value_warm[test$value_warm>1.5] <- 1.5}
  
  bound <- max(c(test$value_warm,test$value_cold), na.rm=T)-min(c(test$value_warm,test$value_cold), na.rm=T)
  med_val <- bound/2+min(c(test$value_warm,test$value_cold), na.rm=T)
  if(bound<0.30){ylimit <- c((med_val-0.15),(med_val+0.15))
  }else{ylimit <- c(min(c(test$value_warm,test$value_cold), na.rm=T),max(c(test$value_warm,test$value_cold), na.rm=T))}
  
  ggplot(test, aes(x=year, group=Country)) + # ici
      geom_line(aes(y=value_warm),col="black" ,size=1.5, alpha=0.7)+
      geom_line(aes(y=value_cold),col="black" ,size=1.5, alpha=0.3)+
      xlab(NULL) + 
      ylab(NULL) + 
      theme_modern() + theme_transparent()+
      scale_y_continuous(limits = ylimit) +
      scale_x_continuous(limits = c(1996,2016)) +
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)}), names(data_pres_long))
            
# Reproject centroids of countries
            
centroid_b <- data.frame(add_rownames(centroids))
names(centroid_b)[1] <- "Country"
centroid_coord <- centroid_b[,c("long","lat")]
centroid_coord <- SpatialPoints(centroid_coord,proj4string = CRS("+proj=longlat +datum=WGS84"))
centroid_coord <- spTransform(centroid_coord, CRS("+init=epsg:27572"))
centroid_b$lon2 <- as.data.frame(centroid_coord)[,1]
centroid_b$lat2 <- as.data.frame(centroid_coord)[,2]

# Avoid overlap

centroid_b$lat2[centroid_b$Country=="Norway"] <- centroid_b$lat2[centroid_b$Country=="Norway"] - 500000
centroid_b$lon2[centroid_b$Country=="Norway"] <- centroid_b$lon2[centroid_b$Country=="Norway"] - 200000

centroid_b$lat2[centroid_b$Country=="Slovenia"] <- centroid_b$lat2[centroid_b$Country=="Slovenia"] - 100000

centroid_b$lat2[centroid_b$Country=="Hungary"] <- centroid_b$lat2[centroid_b$Country=="Hungary"] - 100000

centroid_b$lat2[centroid_b$Country=="Netherlands"] <- centroid_b$lat2[centroid_b$Country=="Netherlands"] + 100000

centroid_b$lon2[centroid_b$Country=="Belgium"] <- centroid_b$lon2[centroid_b$Country=="Belgium"] - 100000

centroid_b$lat2[centroid_b$Country=="Luxembourg"] <- centroid_b$lat2[centroid_b$Country=="Luxembourg"] - 100000
centroid_b$lon2[centroid_b$Country=="Luxembourg"] <- centroid_b$lon2[centroid_b$Country=="Luxembourg"] + 100000

centroid_b$lat2[centroid_b$Country=="Estonia"] <- centroid_b$lat2[centroid_b$Country=="Estonia"] + 100000

centroid_b$lat2[centroid_b$Country=="Lithuania"] <- centroid_b$lat2[centroid_b$Country=="Lithuania"] - 100000

centroid_b$lat2[centroid_b$Country=="Italy"] <- centroid_b$lat2[centroid_b$Country=="Italy"] - 100000
centroid_b$lon2[centroid_b$Country=="Italy"] <- centroid_b$lon2[centroid_b$Country=="Italy"] + 100000

# Remove countries without data 

centroid_b2 <- centroid_b[-c(19,27),]

# Pressures
# Plot high input farm cover

centroid_hico <- tibble(x=centroid_b2$lon2,
                   y=centroid_b2$lat2,
                   width=350000,
                   pie = graph_hico)
library(ggimage)
pres1 + geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=centroid_hico)

# Plot urbanisation

centroid_urb <- tibble(x=centroid_b2$lon2,
                   y=centroid_b2$lat2,
                   width=350000,
                   pie = graph_urb)
                   
pres2 + geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=centroid_urb)

# Plot temperature

centroid_temp <- tibble(x=centroid_b$lon2,
                   y=centroid_b$lat2,
                   width=350000,
                   pie = graph_temp)

pres3 + geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=centroid_temp)

# Plot forest

centroid_forest <- tibble(x=centroid_b$lon2,
                   y=centroid_b$lat2,
                   width=350000,
                   pie = graph_forest)

pres4 + geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=centroid_forest)

# species
# Plot farmland species

centroid_agri <- tibble(x=centroid_b$lon2,
                   y=centroid_b$lat2,
                   width=350000,
                   pie = graph_agri)

pres_tend_agri + geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=centroid_agri)

# Plot woodland species

centroid_wood <- tibble(x=centroid_b$lon2,
                   y=centroid_b$lat2,
                   width=350000,
                   pie = graph_wood)

pres_tend_forest + geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=centroid_wood)

# Plot urban species

centroid_build <- tibble(x=centroid_b$lon2,
                   y=centroid_b$lat2,
                   width=350000,
                   pie = graph_build)

pres_tend_build + geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=centroid_build)

# Plot hot and cold species

centroid_warm_cold <- tibble(x=centroid_b$lon2,
                   y=centroid_b$lat2,
                   width=350000,
                   pie = graph_warm_cold)
                   
pres_tend_temp + geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=centroid_warm_cold)

# Get trend for all countries

# High input farm cover

all_hico <- read.csv("raw_data/aei_ps_inp_1_Data_eu2.csv", header = T)
all_hico <- na.omit(droplevels(all_hico[all_hico$INDIC_AG=="High-input farms" & all_hico$GEO=="EU28" & all_hico$UNIT=="Percentage of area" & all_hico$TIME < 2017,]))

ggplot(all_hico, aes(x=TIME, y=Value)) + # ici
      geom_line(col="black" ,size=1.5, alpha=0.7)+
      xlab(NULL) + ylab(NULL) + 
      theme_modern() +
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)

# Urbanisation

all_urb <- urban_country[1:8,]
all_urb[9,] <- area_country$Value
all_urb2 <- apply(all_urb,2,FUN=function(y){y2 <- y[1:8]*y[9]})
all_urb3 <- data.frame(year=2009:2016,value=apply(all_urb2,1,FUN=function(y){y2 <- sum(y[-c(3:5,20,25)])})/(100*sum(area_country$Value[-c(3:5,14,20,22,25,30)])))

ggplot(all_urb3, aes(x=year, y=value)) + # ici
      geom_line(col="black" ,size=1.5, alpha=0.7)+
      xlab(NULL) + ylab(NULL) + 
      theme_modern() +
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)

# Forest

all_forest <- forest_country[7:27,]
all_forest[22,] <- area_country$Value
all_forest2 <- apply(all_forest,2,FUN=function(y){y2 <- y[1:21]*y[22]})
all_forest3 <- data.frame(year=1996:2016,value=apply(all_forest2,1,FUN=function(y){y2 <- sum(y[-c(2,14,19)])})/sum(area_country$Value[-c(2,14,19)]))

ggplot(all_forest3, aes(x=year, y=value)) + # ici
      geom_line(col="black" ,size=1.5, alpha=0.7)+
      xlab(NULL) + ylab(NULL) + 
      theme_modern() +
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)
            
# Temperature

temp_1996_2016 <- list(temp_1996,temp_1997,temp_1998,temp_1999,
  temp_2000,temp_2001,temp_2002,temp_2003,temp_2004,temp_2005,temp_2006,temp_2007,temp_2008,temp_2009,
  temp_2010,temp_2011,temp_2012,temp_2013,temp_2014,temp_2015,temp_2016)

all_temp <- data.frame(year=1996:2016,value=NA)  

for(i in 1:length(temp_1996_2016)){
  print(i)
  all_temp[i,2] <- mean(unlist(extract(temp_1996_2016[[i]],country5[-c(4,14,20)])),na.rm=T)
}

ggplot(all_temp, aes(x=year, y=value)) + # ici
      geom_line(col="black" ,size=1.5, alpha=0.7)+
      xlab(NULL) + ylab(NULL) + 
      theme_modern() +
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)

```

## Trend analysis

### Estimate trends for each species

```{r}

# Select species with data bewteen 1995 and 2016 (+- tw0 years)

df_trend <- droplevels(subset(df, start_year<=1997))
df_trend <- droplevels(subset(df_trend, end_year>=2014))
df_trend <- droplevels(subset(df_trend, Year %in% c(1995:2016)))

# Estimate trends

# Accounting for uncertainty

trend_species <- ddply(df_trend, .(Species, CountryGroup), .fun=res_trend2, niter=1000, correction=T, mid="first", .parallel = F, .progress = "text")

# Linear slope (similar results and faster)

trend_species <- ddply(df_trend, .(Species, CountryGroup), .fun=function(x){re=summary(lm(Index~Year, x))$coef;return(data.frame(slope=re[2,1],pval=re[2,4],slope_pe=re[2,1]/x$Index[1]))}, .parallel = F, .progress = "text")

```

### Preparing data for partial least square regression (PLS)

#### Species data

```{r}

# Load data

# Species specialisation index
ssi_eu <- read.csv("raw_data/SSI_EU.csv") # from LeViol et al. (2012)
sxi <- read.csv("raw_data/SXI_EU.csv") # from Godet et al. (2015)

# Species complete name
species_name_data <- read.csv("raw_data/species_name_data.csv", header=T)

# Species temperature index
sti<-read.csv("raw_data/STI_Devictor.csv") 

# Species habitat
pecbms_hab <- read.csv2("raw_data/Habitat_class_PECBMS.csv") # available on https://pecbms.info/

# Species life history traits
trait <- read.csv("raw_data/life_history_bird_2018.csv",header = TRUE) # from Storchov?? et al. (2018)

trait$is_migrant <- rep(0, nrow(trait))
trait$is_migrant[which(trait$Short.distance.migrant==1 | trait$Long.distance.migrant==1)] <- 1

trait$is_insectivore <- rep(0, nrow(trait))
trait$is_insectivore[which(trait$Arthropods_B==1 & trait$Other.invertebrates_B==1 & trait$Granivore_B==0 & trait$Folivore_B==0 & trait$Frugivore_B==0 & trait$Fish_B==0 & trait$Other.vertebrates_B==0 & trait$Carrion_B==0 & trait$Omnivore_B==0)] <- 1

trait$is_granivore <- rep(0, nrow(trait))
trait$is_granivore[which(trait$Granivore_B==1 & trait$Folivore_B==0 & trait$Frugivore_B==0 & trait$Arthropods_B==0 & trait$Other.invertebrates_B==0 & trait$Fish_B==0 & trait$Other.vertebrates_B==0 & trait$Carrion_B==0 & trait$Omnivore_B==0)] <- 1

```

#### Merge species trends, traits and pressure data

```{r}

# Species trends and traits

global_data <- merge(trend_species, sxi, by.x="Species", by.y="Name",all.x=T)
global_data <- merge(global_data, sti, by.x="Species", by.y="SPECIES",all.x=T)
global_data <- merge(global_data, pecbms_hab, by="Species",all.x=T)
global_data <- merge(global_data, eunis_hab3, by="Species",all.x=T)
global_data <- merge(global_data,ssi_eu, by="Species",all.x=T)
global_data <- merge(global_data,trait[,c(3,67:90)], by="Species",all.x=T)
global_data$is_farmland <- as.factor(global_data$Habitat=="Farmland")
global_data$is_forest <- as.factor(global_data$Habitat=="Forest")

# Remove outlier trends

global_data <- global_data[abs(global_data$slope) < 50,]

# Merge with pressures

global_data <- merge(global_data,country_data2, by.x="CountryGroup", by.y="country",all.x=T)

# Scale function

Zscore <- function(x){
  return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
}

```

### Applying PLS
```{r}

# Selecting data

data_pls <- global_data[, c("slope","hico_mean","d_hico","for_mean","d_for","urb_mean","d_urb","temp_mean","d_temp")]
data_pls <- data.frame(apply(data_pls,2,Zscore))

# Initiate PLS

cv.modpls <- cv.plsR(data_pls$slope,data_pls[,-1],K=10,nt=10, grouplist = createFolds(data_pls[,1], k = 10, list = F, returnTrain = FALSE))
res.cv.modpls <- cvtable(summary(cv.modpls))
res1 <- plsR(data_pls$slope,data_pls[,-1], nt=10, typeVC="adaptative", pvals.expli=TRUE) # adaptative as NA in data
colSums(res1$pvalstep)

# Searching the best number of component to keep via CV

cv.modpls <- cv.plsR(slope~.,data=data_pls,K=10,nt=10, grouplist = createFolds(data_pls[,1], k = 10, list = F, returnTrain = FALSE),NK=100)
res.cv.modpls <- cvtable(summary(cv.modpls))

# Using CV PRESS, 2 or 3 components must be kept

# Run PLS with 2 components
res <- plsR(slope~.,data=data_pls,nt=2,pvals.expli=TRUE)
trend.bootYT1 <- bootpls(res,typeboot="fmodel_np",R=2000)
pls.ci <- confints.bootpls(trend.bootYT1,indices=2:ncol(data_pls))
plots.confints.bootpls(pls.ci,typeIC="BCa",colIC=c("blue","blue","blue","blue"),legendpos ="topright")

# Run PLS with 3 component
resb <- plsR(slope~.,data=data_pls,nt=3,pvals.expli=TRUE)
trend.bootYT1b=bootpls(resb,typeboot="fmodel_np",R=2000)
pls.cib=confints.bootpls(trend.bootYT1b,indices=2:ncol(data_pls))
plots.confints.bootpls(pls.cib,typeIC="BCa",colIC=c("blue","blue","blue","blue"),legendpos ="topright")

# Using the empirical distribution of the best number of component, we can obtain an empircal measure of the significance of each effect
ind.BCa.YT1 <- (pls.ci[,7] < 0 & pls.ci[,8] < 0)|(pls.ci[,7] > 0 & pls.ci[,8] > 0) 
ind.BCa.YT1b <- (pls.cib[,7] < 0 & pls.cib[,8] < 0)|(pls.cib[,7] > 0 & pls.cib[,8] > 0)
matind <- rbind(YT1b=ind.BCa.YT1b, YT1=ind.BCa.YT1)
pi.e <- prop.table(res.cv.modpls$CVPress)[2:3] %*% matind

# Plot PLS results
coef_plot <- data.frame(var=c("High input farm cover", "High input farm cover trend","Forest cover","Forest cover trend","Urbanisation cover","Urbanisation trend","Mean temperature","Temperature trend"),
val=(cbind(trend.bootYT1$t0[-1,1],trend.bootYT1b$t0[-1,1]) %*% prop.table(res.cv.modpls$CVPress)[2:3]),
                      inf=apply(cbind(pls.ci[,1],pls.cib[,1]),1,min),
                      sup=apply(cbind(pls.ci[,2],pls.cib[,2]),1,max),t(matind), sig=t(pi.e))
coef_plot$col_val <- "ns"
coef_plot$col_val[which(coef_plot$sig>=0.95 & coef_plot$val>0)] <- "pos"
coef_plot$col_val[which(coef_plot$sig>=0.95 & coef_plot$val<0)] <- "neg"
coef_plot$var <- as.character(coef_plot$var)
coef_plot$var <- factor(coef_plot$var, levels = c("Temperature trend","Mean temperature","Urbanisation trend","Urbanisation cover", 
                                                "Forest cover trend","Forest cover", "High input farm cover trend","High input farm cover"))
                                                
ggplot(coef_plot, aes(y=val, x=var)) +
  geom_rect(fill = "#DDF9DD", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, alpha = 0.1) +
  geom_rect(fill = "#F5A9A9", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, alpha = 0.1) +
  geom_bar(stat="identity", position=position_dodge(), alpha=0.7,aes(fill=var)) +
  scale_fill_manual(values=c("High input farm cover"="#D302F9","High input farm cover trend"="#D302F9",
                             "Urbanisation cover"="#196DF6","Urbanisation trend"="#196DF6",
                             "Forest cover"="#1BAE20","Forest cover trend"="#1BAE20",
                             "Mean temperature"="#FA0900","Temperature trend"="#FA0900"))+
  geom_errorbar(aes(ymin=inf, ymax=sup), width=.2,position=position_dodge(.9)) +
  theme_modern() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_hline(yintercept=0, linetype="dashed", size=1) + 
  coord_flip()
  
```

# Causal effect analysis
### Functions for CCM and Smap
```{r}

source("CCM_Smap_function.R")

```

### Pressure time-series

```{r}

# Select time-series

country_data_urb <- data.frame(year=2009:2016,country_data[1:8,])
country_data_urb[country_data_urb==0]<-NA
country_data_urb2 <- melt(country_data_urb, id.vars="year")
names(country_data_urb2)[3] <- "urb"

country_data_temp <- data.frame(year=2007:2016,country_data[68:77,])
country_data_temp2 <- melt(country_data_temp, id.vars="year")
names(country_data_temp2)[3] <- "temp"

country_data_hico <- data.frame(year=2007:2016,country_data[103:112,])
country_data_hico[country_data_hico==0]<-NA
country_data_hico2 <- melt(country_data_hico, id.vars="year")
names(country_data_hico2)[3] <- "hico"

country_data_forest <- data.frame(year=c(2007:2016),country_data[c(133:142),])
country_data_forest2 <- melt(country_data_forest, id.vars="year")
names(country_data_forest2)[3] <- "forest"

# Group

country_data_press <- merge(country_data_temp2, country_data_urb2, by=c("variable","year"),all.x=T)
country_data_press <- merge(country_data_press, country_data_hico2, by=c("variable","year"))
country_data_press <- merge(country_data_press, country_data_forest2, by=c("variable","year"))

df_press <- as.data.frame(df)

press <- country_data_press
press$country <- as.character(press$variable)
press$variable <- NULL
press$country[press$country=="Czech.Republic"] <- "Czech Republic"
press$country[press$country=="UK"] <- "United Kingdom"
press$country[press$country=="Ireland"] <- "Republic of Ireland"
df_press2 <- merge(press, df_press[,c("Species","CountryGroup","Year","Index","Index_SE","Abd")], by.x=c("country","year"), by.y=c("CountryGroup","Year"), all.x=T)
df_press2 <- droplevels(df_press2[which(df_press2$Index >= df_press2$Index_SE),])

# Specify country or species with not enough data (at least 5 year with data, an non null index over the period and a change in pressure through time)

to_remove <- data.frame(na.omit(df_press2) %>% group_by(Species, country) %>% summarize(count=n()))
df_press3 <- merge(df_press2,to_remove, by=c("Species","country"))
to_remove2 <- data.frame(na.omit(df_press2) %>% group_by(Species, country) %>% summarize(sum_ab=sum(Index)))
df_press3 <- merge(df_press3,to_remove2, by=c("Species","country"))
df_press3 <- df_press3[order(df_press3$Species, df_press3$country, df_press3$year),]

# Detrend when needed

df_press4 <- data.frame(droplevels(df_press3[df_press3$count > 4 & df_press3$sum_ab > 20 & df_press3$country!="Luxembourg",]) %>% group_by(Species, country) %>% mutate(temp_std=detrend_data(temp), urb_std=detrend_data(urb), hico_std=detrend_data(hico), forest_std=detrend_data(forest), Index_std=detrend_data(Index), Abd_std=detrend_data(Abd)))

```

### Apply multispatial CCM

It can take several hours.

```{r}

# MultispatialCCM

ccm_sp <- ddply(df_press4, .(Species), .fun = multisp_CCM, niter=1000, .parallel = F, .progress = "text")
ccm_sp2 <- ccm_sp
ccm_sp2[is.na(ccm_sp2)] <- 1

```

### Apply multispatial S-map

```{r}

# S-map

smap_sp <- dlply(droplevels(df_press4), .(Species, country), .fun = smap_fun_signif, ccm_sp2, .parallel = F, .progress = "text")

smap_sp_res <- data.frame(.id=NA,temp_smap=NA,urb_smap=NA,hico_smap=NA,forest_smap=NA)

for(i in levels(df_press4$Species)){
 sub_smap_sp <- smap_sp[grepl(i,names(smap_sp))]
 
 sub_smap_sp_res <- ldply(sub_smap_sp, .fun=function(x){
 
 if(is.na(x$res_ccm)){
  temp_smap <- urb_smap <- hico_smap <- forest_smap <- NA
 }else{
 
  if(x$pvalue<0.05 & x$res_ccm$temp_cause_species<0.05 & !is.null(x$coefficients$temp)){
   temp_smap <- na.omit(x$coefficients$temp)
  }else{temp_smap <- NA}
  
  if(x$pvalue<0.05 & x$res_ccm$urb_cause_species<0.05 & !is.null(x$coefficients$urb)){
   urb_smap <- na.omit(x$coefficients$urb)
  }else{urb_smap <- NA}
  
  if(x$pvalue<0.05 & x$res_ccm$hico_cause_species<0.05 & !is.null(x$coefficients$hico)){
   hico_smap <- na.omit(x$coefficients$hico)
  }else{hico_smap <- NA}
  
  if(x$pvalue<0.05 & x$res_ccm$forest_cause_species<0.05 & !is.null(x$coefficients$forest)){
   forest_smap <- na.omit(x$coefficients$forest)
  }else{forest_smap <- NA}
 }
 
 return(data.frame(temp_smap, urb_smap, hico_smap, forest_smap))
 
 })
 
 smap_sp_res <- rbind(smap_sp_res,sub_smap_sp_res)

}

smap_sp_res <- smap_sp_res[-1,]
smap_sp_res$Species <- sub("\\..*","",smap_sp_res$.id)

# Plot by species

smap_sp_res_long <- melt(smap_sp_res[,-1], id.vars="Species")
for(i in levels(as.factor(smap_sp_res_long$Species))){
 print(ggplot(smap_sp_res_long[smap_sp_res_long$Species==i,], aes(variable, value)) + 
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) + theme_modern() +
    theme(
        legend.position="none",
        plot.title = element_text(size=11)
    )+xlab(i))
}

# Summarise data by species

smap_sp_mean <- data.frame(smap_sp_res[,-1] %>% group_by(Species) %>% summarize(temp=mean(temp_smap, na.rm=T),temp_sd=sd(temp_smap, na.rm=T),
urb=mean(urb_smap, na.rm=T),urb_sd=sd(urb_smap, na.rm=T),
hico=mean(hico_smap, na.rm=T),hico_sd=sd(hico_smap, na.rm=T),
forest=mean(forest_smap, na.rm=T),forest_sd=sd(forest_smap, na.rm=T)))

# Plot sumarised s-map results

data_density <- data.frame(pressure=c(rep("hico",length(na.omit(smap_sp_mean$hico[smap_sp_mean$hico!=0]))),
rep("urb",length(na.omit(smap_sp_mean$urb[smap_sp_mean$urb!=0]))),
rep("forest",length(na.omit(smap_sp_mean$forest[smap_sp_mean$forest!=0]))),
rep("temp",length(na.omit(smap_sp_mean$temp[smap_sp_mean$temp!=0])))),
value=c(scale(na.omit(smap_sp_mean$hico[smap_sp_mean$hico!=0]), center=F),
scale(na.omit(smap_sp_mean$urb[smap_sp_mean$urb!=0]), center=F),
scale(na.omit(smap_sp_mean$forest[smap_sp_mean$forest!=0]), center=F),
scale(na.omit(smap_sp_mean$temp[smap_sp_mean$temp!=0]), center=F)))
                                 
ggplot(data_density, aes(x=value, fill=pressure)) +
    geom_rect(fill = "#E8FAE8", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_rect(fill = "#FACCCC", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_histogram(bins=50,data=data_density[data_density$pressure=="hico",],alpha=0.7, aes(col=pressure), position = position_nudge(y=15)) +
    geom_histogram(bins=50,data=data_density[data_density$pressure=="forest",],alpha=0.4, aes(col=pressure), position = position_nudge(y=12)) +  geom_histogram(bins=50,data=data_density[data_density$pressure=="urb",],alpha=0.4, aes(col=pressure), position = position_nudge(y=4.5)) +  geom_histogram(bins=50,data=data_density[data_density$pressure=="temp",],alpha=0.4, aes(col=pressure)) +
    scale_fill_manual(values=c("hico"="#D302F9","urb"="#196DF6","forest"="#1BAE20","temp"="#FA0900")) +
    scale_color_manual(values=c("hico"="#D302F9","urb"="#196DF6","forest"="#1BAE20","temp"="#FA0900")) +
    theme_modern() +
    labs(x ="Correlation",y="Density") +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
    
```

### Test ccm and smap

#### Influence of countries

```{r}
# Testing how one-by-one country removal affects ccm and smap output

test_country_remove <- function(df_press4, niter=100){
  
  test_ccm_sp_tot <- data.frame(Species=levels(df_press4$Species))
  n <- levels(as.factor(df_press4$country))
    
  for(j in n){
   
    df_press4j <- droplevels(subset(df_press4, country!=j))

    test_ccm_sp <- ddply(df_press4j, .(Species), .fun = multisp_CCM, niter=niter, .parallel = F, .progress = "text")
    test_ccm_sp2 <- test_ccm_sp
    test_ccm_sp2[is.na(test_ccm_sp2)] <- 1

    test_smap_sp <- dlply(droplevels(df_press4j), .(Species, country), .fun = smap_fun_signif, test_ccm_sp2, .parallel = F, .progress = "text")

    test_smap_sp_res <- data.frame(.id=NA,temp_smap=NA,urb_smap=NA,hico_smap=NA,forest_smap=NA)

    for(i in levels(df_press4j$Species)){
      sub_smap_sp <- test_smap_sp[grepl(i,names(test_smap_sp))]
 
      sub_smap_sp_res <- ldply(sub_smap_sp, .fun=function(x){
 
        if(is.na(x$res_ccm)){
          temp_smap <- urb_smap <- hico_smap <- forest_smap <- NA
        }else{
 
        if(x$pvalue<0.05 & x$res_ccm$temp_cause_species<0.05 & !is.null(x$coefficients$temp)){
          temp_smap <- na.omit(x$coefficients$temp)
        }else{temp_smap <- NA}
  
        if(x$pvalue<0.05 & x$res_ccm$urb_cause_species<0.05 & !is.null(x$coefficients$urb)){
          urb_smap <- na.omit(x$coefficients$urb)
        }else{urb_smap <- NA}
  
        if(x$pvalue<0.05 & x$res_ccm$hico_cause_species<0.05 & !is.null(x$coefficients$hico)){
          hico_smap <- na.omit(x$coefficients$hico)
        }else{hico_smap <- NA}
  
        if(x$pvalue<0.05 & x$res_ccm$forest_cause_species<0.05 & !is.null(x$coefficients$forest)){
          forest_smap <- na.omit(x$coefficients$forest)
        }else{forest_smap <- NA}
      }
 
    return(data.frame(temp_smap, urb_smap, hico_smap, forest_smap))
 
  })
 
 test_smap_sp_res <- rbind(test_smap_sp_res,sub_smap_sp_res)

  }

  test_smap_sp_res <- test_smap_sp_res[-1,]
  test_smap_sp_res$Species <- sub("\\..*","",test_smap_sp_res$.id)

  test_smap_sp_mean <- data.frame(test_smap_sp_res[,-1] %>% group_by(Species) %>% summarize(temp=mean(temp_smap, na.rm=T),temp_sd=sd(temp_smap, na.rm=T),
urb=mean(urb_smap, na.rm=T),urb_sd=sd(urb_smap, na.rm=T),
hico=mean(hico_smap, na.rm=T),hico_sd=sd(hico_smap, na.rm=T),
forest=mean(forest_smap, na.rm=T),forest_sd=sd(forest_smap, na.rm=T)))

    names(test_smap_sp_mean)[2:9] <- paste0(j,sep="_",names(test_smap_sp_mean)[2:9])
    
    test_ccm_sp_tot <- merge(test_ccm_sp_tot, test_smap_sp_mean, by="Species", all.x=T)
    
    print(j)
  }
  return(test_ccm_sp_tot)
}

influence_country <- test_country_remove(df_press4, niter=100)
ncol_test <- ncol(influence_country)
influence_country <- merge(influence_country, smap_sp_mean, by="Species", all.x=T)

influence_country$temp_sign_m <- apply(sign(influence_country[,seq(from=2,to=ncol_test,by=8)]), 1,function(x){sum(x,na.rm=T)})
influence_country$temp_nb_neg <- apply(sign(influence_country[,seq(from=2,to=ncol_test,by=8)]), 1,function(x){length(na.omit(x[x<0]))})
influence_country$temp_nb_pos <- apply(sign(influence_country[,seq(from=2,to=ncol_test,by=8)]), 1,function(x){length(na.omit(x[x>0]))})
influence_country$temp_nb_nna <- apply(sign(influence_country[,seq(from=2,to=ncol_test,by=8)]), 1,function(x){sum(abs(x),na.rm=T)})
influence_country$temp_m <- apply(influence_country[,seq(from=2,to=ncol_test,by=8)], 1,function(x){mean(x,na.rm=T)})

influence_country$urb_sign_m <- apply(sign(influence_country[,seq(from=4,to=ncol_test,by=8)]), 1,function(x){sum(x,na.rm=T)})
influence_country$urb_nb_neg <- apply(sign(influence_country[,seq(from=4,to=ncol_test,by=8)]), 1,function(x){length(na.omit(x[x<0]))})
influence_country$urb_nb_pos <- apply(sign(influence_country[,seq(from=4,to=ncol_test,by=8)]), 1,function(x){length(na.omit(x[x>0]))})
influence_country$urb_nb_nna <- apply(sign(influence_country[,seq(from=4,to=ncol_test,by=8)]), 1,function(x){sum(abs(x),na.rm=T)})
influence_country$urb_m <- apply(influence_country[,seq(from=4,to=ncol_test,by=8)], 1,function(x){mean(x,na.rm=T)})

influence_country$hico_sign_m <- apply(sign(influence_country[,seq(from=6,to=ncol_test,by=8)]), 1,function(x){sum(x,na.rm=T)})
influence_country$hico_nb_neg <- apply(sign(influence_country[,seq(from=6,to=ncol_test,by=8)]), 1,function(x){length(na.omit(x[x<0]))})
influence_country$hico_nb_pos <- apply(sign(influence_country[,seq(from=6,to=ncol_test,by=8)]), 1,function(x){length(na.omit(x[x>0]))})
influence_country$hico_nb_nna <- apply(sign(influence_country[,seq(from=6,to=ncol_test,by=8)]), 1,function(x){sum(abs(x),na.rm=T)})
influence_country$hico_m <- apply(influence_country[,seq(from=6,to=ncol_test,by=8)], 1,function(x){mean(x,na.rm=T)})

influence_country$forest_sign_m <- apply(sign(influence_country[,seq(from=8,to=ncol_test,by=8)]), 1,function(x){sum(x,na.rm=T)})
influence_country$forest_nb_neg <- apply(sign(influence_country[,seq(from=8,to=ncol_test,by=8)]), 1,function(x){length(na.omit(x[x<0]))})
influence_country$forest_nb_pos <- apply(sign(influence_country[,seq(from=8,to=ncol_test,by=8)]), 1,function(x){length(na.omit(x[x>0]))})
influence_country$forest_nb_nna <- apply(sign(influence_country[,seq(from=8,to=ncol_test,by=8)]), 1,function(x){sum(abs(x),na.rm=T)})
influence_country$forest_m <- apply(influence_country[,seq(from=8,to=ncol_test,by=8)], 1,function(x){mean(x,na.rm=T)})

influence_country_temp <- apply(influence_country[,seq(from=2,to=(ncol_test+8),by=8)], 2, function(x){length(na.omit(x[x!=0]))})
influence_country_urb <- apply(influence_country[,seq(from=4,to=(ncol_test+8),by=8)], 2, function(x){length(na.omit(x[x!=0]))})
influence_country_hico <- apply(influence_country[,seq(from=6,to=(ncol_test+8),by=8)], 2, function(x){length(na.omit(x[x!=0]))})
influence_country_forest <- apply(influence_country[,seq(from=8,to=(ncol_test+8),by=8)], 2, function(x){length(na.omit(x[x!=0]))})


influence_plot <- melt(influence_country[,c("Species","temp_nb_nna","temp","urb_nb_nna","urb","forest_nb_nna","forest","hico_nb_nna","hico")], id.vars="Species")

influence_plot$group <- "A"
influence_plot$group[influence_plot$Species %in% levels(influence_plot$Species)[85:168]] <- "B"

influence_plot$value2 <- 0
influence_plot$value2[influence_plot$variable=="hico" & !is.na(influence_plot$value)] <- 3
influence_plot$value2[influence_plot$variable=="forest" & !is.na(influence_plot$value)] <- 1
influence_plot$value2[influence_plot$variable=="urb" & !is.na(influence_plot$value)] <- 2
influence_plot$value2[influence_plot$variable=="temp" & !is.na(influence_plot$value)] <- 4

influence_plot$value2 <- as.character(influence_plot$value2)

diag_sp1a <- ggplot(data = droplevels(influence_plot[influence_plot$variable %in% c("hico","forest","urb","temp") & influence_plot$group=="A",]), aes(x=Species, y=variable, fill=value2)) + 
  geom_tile(col="white") + theme_minimal()+ coord_fixed()+
  scale_y_discrete(labels=c("Temperature","Urbanisation","Forest cover","High input farm cover"))+
 scale_fill_manual(values=c("0"="white","1"="#1BAE20","2"="#196DF6","3"="#D302F9","4"="#FA0900"))+
  theme(axis.title = element_blank(), legend.position = "none",plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        
diag_sp1b <- ggplot(data = droplevels(influence_plot[influence_plot$variable %in% c("hico","forest","urb","temp") & influence_plot$group=="B",]), aes(x=Species, y=variable, fill=value2)) + 
  geom_tile(col="white") + theme_minimal()+ coord_fixed()+
  scale_y_discrete(labels=c("Temperature","Urbanisation","Forest cover","High input farm cover"))+
 scale_fill_manual(values=c("0"="white","1"="#1BAE20","2"="#196DF6","3"="#D302F9","4"="#FA0900"))+
  theme(axis.title = element_blank(), legend.position = "none",plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        
library(egg) 
diag_sp1 <- ggarrange(diag_sp1a,diag_sp1b,ncol = 1)

influence_plot$value[influence_plot$variable %in% c("hico_nb_nna","forest_nb_nna","urb_nb_nna","temp_nb_nna")] <- influence_plot$value[influence_plot$variable %in% c("hico_nb_nna","forest_nb_nna","urb_nb_nna","temp_nb_nna")]/24

influence_plot$value[influence_plot$variable %in% c("hico","forest","urb","temp")] <- sign(abs(influence_plot$value[influence_plot$variable %in% c("hico","forest","urb","temp")]))

influence_plot[is.na(influence_plot)] <- 0

diag_sp2a <- ggplot(data = droplevels(influence_plot[influence_plot$group=="A",]), aes(x=Species, y=variable, fill=value)) + 
  geom_tile(col="white")+ theme_minimal()+ coord_fixed()+
  scale_y_discrete(labels=c("Temperature (removing 1 country)","Temperature (all countries)",
  "Urbanisation (removing 1 country)","Urbanisation (all countries)",
  "Forest cover (removing 1 country)","Forest cover (all countries)",
  "High input farm cover (removing 1 country)","High input farm cover (all countries)"))+
  scale_fill_gradient2(low = "blue", high = "black", mid = "white", midpoint = 0, limit = c(-1,1))+
  theme(axis.title = element_blank(), legend.position = "none",plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        
diag_sp2b <- ggplot(data = droplevels(influence_plot[influence_plot$group=="B",]), aes(x=Species, y=variable, fill=value)) + 
  geom_tile(col="white")+ theme_minimal()+ coord_fixed()+
  scale_y_discrete(labels=c("Temperature (removing 1 country)","Temperature (all countries)",
  "Urbanisation (removing 1 country)","Urbanisation (all countries)",
  "Forest cover (removing 1 country)","Forest cover (all countries)",
  "High input farm cover (removing 1 country)","High input farm cover (all countries)"))+
  scale_fill_gradient2(low = "blue", high = "black", mid = "white", midpoint = 0, limit = c(-1,1))+
  theme(axis.title = element_blank(), legend.position = "none",plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

diag_sp2 <- ggarrange(diag_sp2a,diag_sp2b,ncol = 1)

```

#### Changing pressure strength

```{r}
# test avec force changeante pour chaque pays

test_template <- function(df_press4, nb_iter2=1000, force_int=1.25){

  nb_iter <- 100

  lg_ts <- 10
  
  nb_country <- length(levels(as.factor(df_press4$country)))
  
  country_name <- levels(as.factor(df_press4$country))
  
  result_tempo_tot <- data.frame(species_cause_pressure=NA, pressure_cause_species=NA)
  
  for(row_num in 1:nb_iter2){
    
      rancor <- data.frame(country=NA,Year=NA,Index=NA,value=NA)
  
      for(u in country_name){
        
          ccm_data_out <- make_ccm_data(number_of_chains = 1, time=lg_ts, Sstr = force_int, seednum = round(runif(1,min=1, max=10000)))
          
          rancor0 <- data.frame(country=u,Year=na.omit(ccm_data_out$time_ccm),Index=na.omit(ccm_data_out$Accm),value=na.omit(ccm_data_out$Bccm))
          
          rancor <- as.data.frame(rbind(rancor,rancor0))
        }
        
      rancor <- na.omit(rancor)
      rancor$country <- as.factor(rancor$country)
      
      result_tempo <- multisp_CCM_test(rancor,niter=nb_iter)
      
      result_tempo_tot <- rbind(result_tempo_tot,result_tempo)
        
    }
    
    return(result_tempo_tot)
}

# with varying proportion of countries where strength = 0

test_ccm_force1 <- function(df_press4, nb_iter2=1000){

  nb_iter <- 100

  lg_ts <- 10
  
  nb_country <- length(levels(as.factor(df_press4$country)))
  
  country_name <- levels(as.factor(df_press4$country))
  
  result_tempo_tot <- data.frame(species_cause_pressure=NA, pressure_cause_species=NA, proba=NA)
  
  for(row_num in 1:nb_iter2){
      print(row_num)
  
      if(row_num <= round(nb_iter2/6)){proba <- 0}
      if(row_num > round(nb_iter2/6) & row_num <= round(nb_iter2/3)){proba <- 0.2}
      if(row_num > round(nb_iter2/3) & row_num <= round(nb_iter2/2)){proba <- 0.4}
      if(row_num > round(nb_iter2/2) & row_num <= round(2*nb_iter2/3)){proba <- 0.6}
      if(row_num > round(2*nb_iter2/3) & row_num <= round(5*nb_iter2/6)){proba <- 0.8}
      if(row_num > round(5*nb_iter2/6)){proba <- 1}
      
      rancor <- data.frame(country=NA,Year=NA,Index=NA,value=NA)
  
      for(u in country_name){
      
          force_int <- 1.25*rbinom(1,1,proba)
        
          ccm_data_out <- make_ccm_data(number_of_chains = 1, time=lg_ts, Sstr = force_int, seednum = round(runif(1,min=1, max=10000)))
          
          rancor0 <- data.frame(country=u,Year=na.omit(ccm_data_out$time_ccm),Index=na.omit(ccm_data_out$Accm),value=na.omit(ccm_data_out$Bccm))
          
          rancor <- as.data.frame(rbind(rancor,rancor0))
        }
        
      rancor <- na.omit(rancor)
      rancor$country <- as.factor(rancor$country)
      
      result_tempo <- multisp_CCM_test(rancor,niter=nb_iter)
      
      result_tempo_tot <- rbind(result_tempo_tot,data.frame(result_tempo, proba=proba))
        
    }
    
    return(result_tempo_tot)
}

# with varying values of strength

test_ccm_force2 <- function(df_press4, nb_iter2=1000){

  nb_iter <- 100

  lg_ts <- 10
  
  nb_country <- length(levels(as.factor(df_press4$country)))
  
  country_name <- levels(as.factor(df_press4$country))
  
  result_tempo_tot <- data.frame(species_cause_pressure=NA, pressure_cause_species=NA,low_bound=NA)
  
  for(row_num in 1:nb_iter2){
      print(row_num)
  
      if(row_num <= round(nb_iter2/5)){low_bound <- 1.2; high_bound <- 1.3}
      if(row_num > round(nb_iter2/5) & row_num <= round(2*nb_iter2/5)){low_bound <- 0.9; high_bound <- 1.6}
      if(row_num > round(2*nb_iter2/5) & row_num <= round(3*nb_iter2/5)){low_bound <- 0.6; high_bound <- 1.9}
      if(row_num > round(3*nb_iter2/5) & row_num <= round(4*nb_iter2/5)){low_bound <- 0.3; high_bound <- 2.2}
      if(row_num > round(4*nb_iter2/5)){low_bound <- 0; high_bound <- 2.5}
      
      rancor <- data.frame(country=NA,Year=NA,Index=NA,value=NA)
  
      for(u in country_name){
      
          force_int <- runif(1,min=low_bound, max=high_bound)
        
          ccm_data_out <- make_ccm_data(number_of_chains = 1, time=lg_ts, Sstr = force_int, seednum = round(runif(1,min=1, max=10000)))
          
          rancor0 <- data.frame(country=u,Year=na.omit(ccm_data_out$time_ccm),Index=na.omit(ccm_data_out$Accm),value=na.omit(ccm_data_out$Bccm))
          
          rancor <- as.data.frame(rbind(rancor,rancor0))
        }
        
      rancor <- na.omit(rancor)
      rancor$country <- as.factor(rancor$country)
      
      result_tempo <- multisp_CCM_test(rancor,niter=nb_iter)
      
      result_tempo_tot <- rbind(result_tempo_tot,data.frame(result_tempo,low_bound=low_bound))
        
    }
    
    return(result_tempo_tot)
}

# with one time serie with strength != 0

test_ccm_force3 <- function(df_press4, nb_iter2=1000){

  nb_iter <- 100

  lg_ts <- 10
  
  nb_country <- length(levels(as.factor(df_press4$country)))
  
  country_name <- levels(as.factor(df_press4$country))
  
  result_tempo_tot <- data.frame(species_cause_pressure=NA, pressure_cause_species=NA, force_int2=NA)
  
  for(row_num in 1:nb_iter2){
      print(row_num)
  
      if(row_num <= round(nb_iter2/5)){force_int2 <- 0.3}
      if(row_num > round(nb_iter2/5) & row_num <= round(2*nb_iter2/5)){force_int2 <- 0.85}
      if(row_num > round(2*nb_iter2/5) & row_num <= round(3*nb_iter2/5)){force_int2 <- 1.4}
      if(row_num > round(3*nb_iter2/5) & row_num <= round(4*nb_iter2/5)){force_int2 <- 1.95}
      if(row_num > round(4*nb_iter2/5)){force_int2 <- 2.5}
      
      rancor <- data.frame(country=NA,Year=NA,Index=NA,value=NA)
      
      which_diff_0 <- round(runif(1,-0.49, (nb_country+0.49)))
  
      for(u in 1:nb_country){
      
          if(u==which_diff_0){
            force_int <- force_int2
          }else{
            force_int <- 0
          }
        
          ccm_data_out <- make_ccm_data(number_of_chains = 1, time=lg_ts, Sstr = force_int, seednum = round(runif(1,min=1, max=10000)))
          
          rancor0 <- data.frame(country=country_name[u],Year=na.omit(ccm_data_out$time_ccm),Index=na.omit(ccm_data_out$Accm),value=na.omit(ccm_data_out$Bccm))
          
          rancor <- as.data.frame(rbind(rancor,rancor0))
        }
        
      rancor <- na.omit(rancor)
      rancor$country <- as.factor(rancor$country)
      
      result_tempo <- multisp_CCM_test(rancor,niter=nb_iter)
      
      result_tempo_tot <- rbind(result_tempo_tot,data.frame(result_tempo,force_int2=force_int2))
        
    }
    
    return(result_tempo_tot)
}

res_test1 <- test_ccm_force1(df_press4, nb_iter2=1000)
res_test2 <- test_ccm_force2(df_press4, nb_iter2=1000)
res_test3 <- test_ccm_force3(df_press4, nb_iter2=1000)

```

### Applying PLS on pressure influence vs. species traits

#### Prepare data

```{r}

sp_data <- smap_sp_mean

sp_data <- merge(sp_data,sxi, by.x="Species", by.y="Name",all.x=T)
sp_data <- merge(sp_data,sti, by.x="Species", by.y="SPECIES",all.x=T)
sp_data <- merge(sp_data,trait[,c("Species","is_granivore","is_migrant","is_insectivore")], by="Species",all.x=T)
sp_data <- merge(sp_data,pecbms_hab, by="Species",all.x=T)
sp_data <- merge(sp_data,ssi_eu, by="Species",all.x=T)
sp_data <- merge(sp_data,eunis_hab3, by="Species",all.x=T)

sp_data$temp[is.na(sp_data$temp)] <- 0
sp_data$urb[is.na(sp_data$urb)] <- 0
sp_data$hico[is.na(sp_data$hico)] <- 0
sp_data$forest[is.na(sp_data$forest)] <- 0
sp_data$is_urban[is.na(sp_data$is_urban)] <- 0
sp_data$is_forest <- as.factor(sp_data$Habitat=="Forest")
sp_data$is_farmland <- as.factor(sp_data$Habitat=="Farmland")

```

#### Temperature vs traits

```{r}

# Select data for PLS

trait_inter_data_temp <- sp_data[, c("temp","is_farmland","is_forest","STI","SSI","is_migrant","is_granivore","is_insectivore","is_urban")]

# Scale data

trait_inter_data_temp$STI <- scale(trait_inter_data_temp$STI)
trait_inter_data_temp$SSI <- scale(trait_inter_data_temp$SSI)
trait_inter_data_temp$is_farmland <- as.numeric(trait_inter_data_temp$is_farmland)-1
trait_inter_data_temp$is_forest <- as.numeric(trait_inter_data_temp$is_forest)-1
trait_inter_data_temp$is_urban <- as.numeric(trait_inter_data_temp$is_urban)

# Find the number of latent value

cv.modpls_temp <- cv.plsR(trait_inter_data_temp$temp,trait_inter_data_temp[,-1],nt=10,grouplist = createFolds(trait_inter_data_temp[,1], k = 10, list = F, returnTrain = FALSE))
res.cv.modpls_temp <- cvtable(summary(cv.modpls_temp))
res_pls_temp <- plsR(trait_inter_data_temp$temp,trait_inter_data_temp[,-1], nt=10, typeVC="adaptative", pvals.expli=TRUE) # adaptative because NA in some columns
colSums(res1$pvalstep)
cv.modpls_temp <- cv.plsR(temp~.,data=trait_inter_data_temp,nt=10,grouplist = createFolds(trait_inter_data_temp[,1], k = 10, list = F, returnTrain = FALSE),NK=100)
res.cv.modpls_temp <- cvtable(summary(cv.modpls_temp))

# Using CV PRESS, 1 component must be kept
# Run PLS

res_pls_temp <- plsR(temp~.,data=trait_inter_data_temp,nt=1,pvals.expli=TRUE)

# Plot PLS

temp.bootYT1 <- bootpls(res_pls_temp,typeboot="fmodel_np",R=10000)
boxplots.bootpls(temp.bootYT1,indices=2:ncol(trait_inter_data_temp))
temp.ci <- confints.bootpls(temp.bootYT1,indices=2:ncol(trait_inter_data_temp))
plots.confints.bootpls(temp.ci,typeIC="BCa",colIC=c("blue","blue","blue","blue"),
                       legendpos ="topright")

ind.BCa.tempYT1 <- (temp.ci[,7] < 0 & temp.ci[,8] < 0) | (temp.ci[,7] > 0 & temp.ci[,8] > 0)

# Save results

matind <- rbind(YT1=ind.BCa.tempYT1)
pi.e <- (prop.table(res.cv.modpls_temp$CVPress)[c(1)]/sum(prop.table(res.cv.modpls_temp$CVPress)[c(1)])) %*% matind

coef_plot_temp <- data.frame(var=c("Farmland","Forest","STI","SSI","Migrant",
                                 "Granivorous diet","Invertebrate diet","Synanthropy"),
                                 val=temp.bootYT1$t0[-1,1],
                                 inf=temp.ci[,7],sup=temp.ci[,8],t(matind),
                                 sig=t(pi.e))
coef_plot_temp$col_val <- "ns"
coef_plot_temp$col_val[which(coef_plot_temp$sig >= 0.95 & coef_plot_temp$val > 0)] <- "pos"
coef_plot_temp$col_val[which(coef_plot_temp$sig >= 0.95 & coef_plot_temp$val < 0)] <- "neg"

```

#### Urbanisation vs traits

```{r}

# Select data for PLS
set.seed(1)
trait_inter_data_urb <- sp_data[, c("urb","is_farmland","is_forest","STI","SSI","is_migrant","is_granivore","is_insectivore","is_urban")]

# Scale data

trait_inter_data_urb$STI <- scale(trait_inter_data_urb$STI)
trait_inter_data_urb$SSI <- scale(trait_inter_data_urb$SSI)
trait_inter_data_urb$is_farmland <- as.numeric(trait_inter_data_urb$is_farmland)-1
trait_inter_data_urb$is_forest <- as.numeric(trait_inter_data_urb$is_forest)-1
trait_inter_data_urb$is_urban <- as.numeric(trait_inter_data_urb$is_urban)

# Find the number of latent value

cv.modpls_urb <- cv.plsR(trait_inter_data_urb$urb,trait_inter_data_urb[,-1],nt=10,grouplist = createFolds(trait_inter_data_urb[,1], k = 10, list = F, returnTrain = FALSE))
res.cv.modpls_urb <- cvtable(summary(cv.modpls_urb))
res_pls_urb <- plsR(trait_inter_data_urb$urb,trait_inter_data_urb[,-1], nt=10, typeVC="adaptative", pvals.expli=TRUE) # adaptative because NA in some columns
colSums(res1$pvalstep)
cv.modpls_urb <- cv.plsR(urb~.,data=trait_inter_data_urb,nt=10,grouplist = createFolds(trait_inter_data_urb[,1], k = 10, list = F, returnTrain = FALSE),NK=100)
res.cv.modpls_urb <- cvtable(summary(cv.modpls_urb))

# Using CV PRESS, 1 component must be kept
# Run PLS
res_pls_urb <- plsR(urb~., data=trait_inter_data_urb, nt=1, pvals.expli=TRUE)

# Plot PLS

urb.bootYT1 <- bootpls(res_pls_urb,typeboot="fmodel_np",R=10000)
boxplots.bootpls(urb.bootYT1,indices=2:ncol(trait_inter_data_urb))
urb.ci <- confints.bootpls(urb.bootYT1,indices=2:ncol(trait_inter_data_urb))
plots.confints.bootpls(urb.ci,typeIC="BCa",colIC=c("blue","blue","blue","blue"),
                       legendpos ="topright")

ind.BCa.urbYT1 <- (urb.ci[,7] < 0 & urb.ci[,8] < 0) | (urb.ci[,7] > 0 & urb.ci[,8] > 0)

# Save results

matind <- rbind(YT1=ind.BCa.urbYT1)
pi.e <- (prop.table(res.cv.modpls_urb$CVPress)[c(1)]/sum(prop.table(res.cv.modpls_urb$CVPress)[c(1)])) %*% matind

coef_plot_urb <- data.frame(var=c("Farmland","Forest","STI","SSI","Migrant",
                                 "Granivorous diet","Invertebrate diet","Synanthropy"),
                                 val=urb.bootYT1$t0[-1,1],
                                 inf=urb.ci[,7],sup=urb.ci[,8],t(matind),
                                 sig=t(pi.e))
coef_plot_urb$col_val <- "ns"
coef_plot_urb$col_val[which(coef_plot_urb$sig >= 0.95 & coef_plot_urb$val > 0)] <- "pos"
coef_plot_urb$col_val[which(coef_plot_urb$sig >= 0.95 & coef_plot_urb$val < 0)] <- "neg"

```

#### High input farm cover vs traits

```{r}

# Select data for PLS

trait_inter_data_hico <- sp_data[, c("hico","is_farmland","is_forest","STI","SSI","is_migrant","is_granivore","is_insectivore","is_urban")]

# Scale data

trait_inter_data_hico$STI <- scale(trait_inter_data_hico$STI)
trait_inter_data_hico$SSI <- scale(trait_inter_data_hico$SSI)
trait_inter_data_hico$is_farmland <- as.numeric(trait_inter_data_hico$is_farmland)-1
trait_inter_data_hico$is_forest <- as.numeric(trait_inter_data_hico$is_forest)-1
trait_inter_data_hico$is_urban <- as.numeric(trait_inter_data_hico$is_urban)

# Find the number of latent value

cv.modpls_hico <- cv.plsR(trait_inter_data_hico$hico,trait_inter_data_hico[,-1],nt=10,grouplist = createFolds(trait_inter_data_hico[,1], k = 10, list = F, returnTrain = FALSE))
res.cv.modpls_hico <- cvtable(summary(cv.modpls_hico))
res_pls_hico <- plsR(trait_inter_data_hico$hico,trait_inter_data_hico[,-1], nt=10, typeVC="adaptative", pvals.expli=TRUE) # adaptative because NA in some columns
colSums(res1$pvalstep)
cv.modpls_hico <- cv.plsR(hico~.,data=trait_inter_data_hico,nt=10,grouplist = createFolds(trait_inter_data_hico[,1], k = 10, list = F, returnTrain = FALSE),NK=100)
res.cv.modpls_hico <- cvtable(summary(cv.modpls_hico))

# Using CV PRESS, 1 component must be kept
# Run PLS

res_pls_hico <- plsR(hico~.,data=trait_inter_data_hico,nt=1,pvals.expli=TRUE)

# Plot PLS

hico.bootYT1 <- bootpls(res_pls_hico,typeboot="fmodel_np",R=10000)
boxplots.bootpls(hico.bootYT1,indices=2:ncol(trait_inter_data_hico))
hico.ci <- confints.bootpls(hico.bootYT1,indices=2:ncol(trait_inter_data_hico))
plots.confints.bootpls(hico.ci,typeIC="BCa",colIC=c("blue","blue","blue","blue"),
                       legendpos ="topright")

ind.BCa.hicoYT1 <- (hico.ci[,7] < 0 & hico.ci[,8] < 0) | (hico.ci[,7] > 0 & hico.ci[,8] > 0)

# Save results

matind <- rbind(YT1=ind.BCa.hicoYT1)
pi.e <- (prop.table(res.cv.modpls_hico$CVPress)[c(1)]/sum(prop.table(res.cv.modpls_hico$CVPress)[c(1)])) %*% matind

coef_plot_hico <- data.frame(var=c("Farmland","Forest","STI","SSI","Migrant",
                                 "Granivorous diet","Invertebrate diet","Synanthropy"),
                                 val=hico.bootYT1$t0[-1,1],
                                 inf=hico.ci[,7],sup=hico.ci[,8],t(matind),
                                 sig=t(pi.e))
coef_plot_hico$col_val <- "ns"
coef_plot_hico$col_val[which(coef_plot_hico$sig >= 0.95 & coef_plot_hico$val > 0)] <- "pos"
coef_plot_hico$col_val[which(coef_plot_hico$sig >= 0.95 & coef_plot_hico$val < 0)] <- "neg"

```

#### Forest vs traits

```{r}

# Select data for PLS
set.seed(1)
trait_inter_data_forest <- sp_data[, c("forest","is_farmland","is_forest","STI","SSI","is_migrant","is_granivore","is_insectivore","is_urban")]

# Scale data

trait_inter_data_forest$STI <- scale(trait_inter_data_forest$STI)
trait_inter_data_forest$SSI <- scale(trait_inter_data_forest$SSI)
trait_inter_data_forest$is_farmland <- as.numeric(trait_inter_data_forest$is_farmland)-1
trait_inter_data_forest$is_forest <- as.numeric(trait_inter_data_forest$is_forest)-1
trait_inter_data_forest$is_urban <- as.numeric(trait_inter_data_forest$is_urban)

# Find the number of latent value

cv.modpls_forest <- cv.plsR(trait_inter_data_forest$forest,trait_inter_data_forest[,-1],nt=10,grouplist = createFolds(trait_inter_data_forest[,1], k = 10, list = F, returnTrain = FALSE))
res.cv.modpls_forest <- cvtable(summary(cv.modpls_forest))
res_pls_forest <- plsR(trait_inter_data_forest$forest,trait_inter_data_forest[,-1], nt=10, typeVC="adaptative", pvals.expli=TRUE) # adaptative because NA in some columns
colSums(res1$pvalstep)
cv.modpls_forest <- cv.plsR(forest~.,data=trait_inter_data_forest,nt=10,grouplist = createFolds(trait_inter_data_forest[,1], k = 10, list = F, returnTrain = FALSE),NK=100)
res.cv.modpls_forest <- cvtable(summary(cv.modpls_forest))

# Using CV PRESS, 1 or 2 components must be kept
# Run PLS with 1 component
res_pls_forest <- plsR(forest~.,data=trait_inter_data_forest,nt=1,pvals.expli=TRUE)

# Plot PLS

forest.bootYT1 <- bootpls(res_pls_forest,typeboot="fmodel_np",R=10000)
boxplots.bootpls(forest.bootYT1,indices=2:ncol(trait_inter_data_forest))
forest.ci <- confints.bootpls(forest.bootYT1,indices=2:ncol(trait_inter_data_forest))
plots.confints.bootpls(forest.ci,typeIC="BCa",colIC=c("blue","blue","blue","blue"),
                       legendpos ="topright")
                       
# Run PLS with 2 components
res_pls_forest2 <- plsR(forest~.,data=trait_inter_data_forest,nt=2,pvals.expli=TRUE)

# Plot PLS

forest.bootYT2 <- bootpls(res_pls_forest2,typeboot="fmodel_np",R=10000)
boxplots.bootpls(forest.bootYT2,indices=2:ncol(trait_inter_data_forest))
forest.ci2 <- confints.bootpls(forest.bootYT2,indices=2:ncol(trait_inter_data_forest))
plots.confints.bootpls(forest.ci2,typeIC="BCa",colIC=c("blue","blue","blue","blue"),
                       legendpos ="topright")                      

ind.BCa.forestYT1 <- (forest.ci[,7] < 0 & forest.ci[,8] < 0) | (forest.ci[,7] > 0 & forest.ci[,8] > 0)
ind.BCa.forestYT2 <- (forest.ci2[,7] < 0 & forest.ci2[,8] < 0) | (forest.ci2[,7] > 0 & forest.ci2[,8] > 0)

# Save results

matind <- rbind(YT1=ind.BCa.forestYT1, YT2=ind.BCa.forestYT2)
pi.e <- (prop.table(res.cv.modpls_forest$CVPress)[c(1:2)]/sum(prop.table(res.cv.modpls_forest$CVPress)[c(1:2)])) %*% matind

coef_plot_forest <- data.frame(var=c("Farmland","Forest","STI","SSI","Migrant",
                                 "Granivorous diet","Invertebrate diet","Synanthropy"),
                                 val=forest.bootYT2$t0[-1,1],
                                 inf=forest.ci2[,7],sup=forest.ci2[,8],t(matind),
                                 sig=t(pi.e))
coef_plot_forest$col_val <- "ns"
coef_plot_forest$col_val[which(coef_plot_forest$sig >= 0.95 & coef_plot_forest$val > 0)] <- "pos"
coef_plot_forest$col_val[which(coef_plot_forest$sig >= 0.95 & coef_plot_forest$val < 0)] <- "neg"

```

#### Plot pressure vs traits

```{r}

library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(networkD3)

# Group all data to get flows

data_trait_pression <- rbind(data.frame(coef_plot_hico[,c("var","val","inf","sup","col_val")],pressure="High input farm cover",value2=coef_plot_hico$val[coef_plot_hico$col_val!="ns"]/sum(abs(coef_plot_hico$val[coef_plot_hico$col_val!="ns"]))),
                           data.frame(coef_plot_forest[,c("var","val","inf","sup","col_val")],pressure="Forest cover",value2=coef_plot_forest$val/sum(abs(coef_plot_forest$val[coef_plot_forest$col_val!="ns"]))),
                           
data.frame(coef_plot_urb[,c("var","val","inf","sup","col_val")],pressure="Urbanisation",value2=coef_plot_urb$val[coef_plot_urb$col_val!="ns"]/sum(abs(coef_plot_urb$val[coef_plot_urb$col_val!="ns"]))),

data.frame(coef_plot_temp[,c("var","val","inf","sup","col_val")],pressure="Temperature",value2=coef_plot_temp$val[coef_plot_temp$col_val!="ns"]/sum(abs(coef_plot_temp$val[coef_plot_temp$col_val!="ns"]))))
              
data_trait_pression$var <- as.character(data_trait_pression$var) 
data_trait_pression[data_trait_pression=="STI"] <- "Species Temperature Index"
data_trait_pression[data_trait_pression=="SSI"] <- "Species Specialisation Index"

data_trait_pression2 <- data.frame(source=data_trait_pression$pressure,target=data_trait_pression$var,
                                 value=abs(data_trait_pression$value2), col_link=data_trait_pression$col_val)

# From these flows we need to create a node data frame: it lists every entities involved in the flow

nodes <- data.frame(name=c(as.character(data_trait_pression2$source), as.character(data_trait_pression2$target)) %>% unique())
nodes$group <- as.factor(c("Input","Forestc","Urbanisation","Temperature","trait","trait","trait","trait","trait","trait","trait","trait"))

# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.

data_trait_pression2$IDsource <- match(data_trait_pression2$source, nodes$name)-1 
data_trait_pression2$IDtarget <- match(data_trait_pression2$target, nodes$name)-1

# Prepare colour scale

ColourScal <-  'd3.scaleOrdinal() .domain(["neg","ns","pos","Temperature","Urbanisation","Input","Forestc","trait","trait","trait","trait","trait","trait","trait","trait"]) .range(["#F6CECE","#E6E6E6", "#CEF6CE", "#FA0900","#196DF6","#D302F9","#1BAE20", "black", "black", "black", "black", "black", "black", "black", "black"])'

# Make the Network

sankeyNetwork(Links = data_trait_pression2, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name",  LinkGroup = "col_link", NodeGroup="group",
              sinksRight=FALSE, colourScale=ColourScal, nodeWidth=60, fontSize=13, nodePadding=20, iteration=0)

```

## Supplementary material

### Data correlation

```{r}

library(GGally)

data_sp_trait <- merge(data.frame(Species=levels(df$Species)), sxi, by.x="Species", by.y="Name",all.x=T)
data_sp_trait <- merge(data_sp_trait, sti, by.x="Species", by.y="SPECIES",all.x=T)
data_sp_trait <- merge(data_sp_trait, pecbms_hab, by="Species",all.x=T)
data_sp_trait <- merge(data_sp_trait, eunis_hab3, by="Species",all.x=T)
data_sp_trait <- merge(data_sp_trait,ssi_eu, by="Species",all.x=T)
data_sp_trait <- merge(data_sp_trait,trait[,c(3,67:90)], by="Species",all.x=T)
data_sp_trait$is_farmland <- as.factor(data_sp_trait$Habitat=="Farmland")
data_sp_trait$is_forest <- as.factor(data_sp_trait$Habitat=="Forest")

pair_to_plot <- data_sp_trait[!duplicated(data_sp_trait$Species),]
pair_to_plot <- pair_to_plot[,c("is_farmland","is_forest","is_urban","STI","SSI","is_migrant","is_granivore","is_insectivore")]

pair_to_plot$is_migrant <- as.character(pair_to_plot$is_migrant)
pair_to_plot$is_migrant[pair_to_plot$is_migrant=="0"] <- "FALSE"
pair_to_plot$is_migrant[pair_to_plot$is_migrant=="1"] <- "TRUE"
pair_to_plot$is_migrant <- as.factor(pair_to_plot$is_migrant)

pair_to_plot$is_granivore <- as.character(pair_to_plot$is_granivore)
pair_to_plot$is_granivore[pair_to_plot$is_granivore=="0"] <- "FALSE"
pair_to_plot$is_granivore[pair_to_plot$is_granivore=="1"] <- "TRUE"
pair_to_plot$is_granivore <- as.factor(pair_to_plot$is_granivore)

pair_to_plot$is_insectivore <- as.character(pair_to_plot$is_insectivore)
pair_to_plot$is_insectivore[pair_to_plot$is_insectivore=="0"] <- "FALSE"
pair_to_plot$is_insectivore[pair_to_plot$is_insectivore=="1"] <- "TRUE"
pair_to_plot$is_insectivore <- as.factor(pair_to_plot$is_insectivore)

pair_to_plot$is_urban[is.na(pair_to_plot$is_urban)] <- "FALSE"
pair_to_plot$SSI[is.na(pair_to_plot$SSI)] <- median(na.omit(pair_to_plot$SSI))

names(pair_to_plot) <- c("Farmland","Woodland","Urban","Species Temperature Index","Species Specialisation Index","Long distance migrant","Granivore diet",
                        "Invertebrate-based diet")
                        
ggpairs(droplevels(na.omit(pair_to_plot)))

ggsave("colin1.png",
       width = 15,
       height = 15,
       dpi = 400)
pair_to_plot2<-global_data3b[,c("high_input_cover2","d_hic","forest","d_forest","clc_mean","d_clc","temp_mean","d_temp")]
names(pair_to_plot2)<-c("High input farm cover","High input cover trend","Forest cover","Forest cover trend","Artificialised cover",
                       "Artificialisation trend","Temperature","Temperature trend")
ggpairs(pair_to_plot2)


```

### PLSR with farm size as a pressure

```{r}

# Load data
# https://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=ef_m_farmleg&lang=en

size_farm_country <- read.csv("raw_data/ef_m_farmleg_1_Data.csv", header = T)
 
size_farm_country <- na.omit(droplevels(size_farm_country[size_farm_country$INDIC_AGR %in% c("Utilised agricultural area - hectare","Farm - number"),]))
size_farm_country$Value <- as.character(size_farm_country$Value)
size_farm_country$Value <- gsub(",","",size_farm_country$Value)
size_farm_country$Value <- as.numeric(size_farm_country$Value)
size_farm_country <- dcast(size_farm_country, INDIC_AGR+TIME~GEO, fun.aggregate=sum, value.var="Value")

size_farm_country$INDIC_AGR <- as.character(size_farm_country$INDIC_AGR)
size_farm_country$INDIC_AGR[size_farm_country$INDIC_AGR=="Farm - number"] <- "Number"
size_farm_country$INDIC_AGR[size_farm_country$INDIC_AGR=="Utilised agricultural area - hectare"] <- "UAA"

for(i in 1:5){
size_farm_country[(i+10),1] <- "Farm size"
size_farm_country[(i+10),2] <- size_farm_country$TIME[i]
size_farm_country[(i+10),3:34] <- size_farm_country[(i+5),3:34]/size_farm_country[i,3:34]
}

size_farm_country <- droplevels(size_farm_country[11:15,2:34])

names(size_farm_country) <- c("Year","Austria","Belgium","Bulgaria","Croatia","Cyprus","Czech Republic","Denmark","Estonia","Finland","France","Germany","Greece","Hungary","Iceland","Ireland","Italy","Latvia","Lithuania","Luxembourg","Malta","Montenegro","Netherlands","Norway","Poland","Portugal","Romania","Slovakia","Slovenia","Spain","Sweden","Switzerland","UK")

row.names(size_farm_country) <- paste0("size_farm",sep="_",size_farm_country$Year)
 
size_farm_country$Year <- size_farm_country$Iceland <- size_farm_country$Montenegro <- NULL
 
size_farm_country[6,] <- apply(size_farm_country[1:5,], 2, function(x){mean(x, na.rm=T)})
 
row.names(size_farm_country)[6] <- "size_farm_mean"

size_farm_country[7,] <- apply(size_farm_country[1:5,], 2, function(x){summary(lm(x~c(2005,2007,2010,2013,2016)))$coef[2,1]})/size_farm_country[1,]
 
row.names(size_farm_country)[7] <- "d_size_farm"

size_farm_country <- size_farm_country[,sort(names(size_farm_country))]
 
size_farm_country <- data.frame(size_farm_country[,c(1:13)], Iceland=NA, size_farm_country[,c(14:30)])
 
names(size_farm_country)[6] <- "Czech Republic"

# Merge data

country_data3 <- rbind(country_data, size_farm_country)

# Clean up the final dataset

country_data3 <- as.data.frame(t(country_data3))
country_data3$country <- row.names(country_data3)
country_data3$country[country_data3$country=="UK"]<-"United Kingdom"
country_data3[country_data3==0] <- NA
country_data3$country <- as.factor(country_data3$country)

# Merge species trends, traits and pressures

global_data2 <- merge(trend_species, sxi, by.x="Species", by.y="Name",all.x=T)
global_data2 <- merge(global_data2, sti, by.x="Species", by.y="SPECIES",all.x=T)
global_data2 <- merge(global_data2, pecbms_hab, by="Species",all.x=T)
global_data2 <- merge(global_data2, eunis_hab3, by="Species",all.x=T)
global_data2 <- merge(global_data2,ssi_eu, by="Species",all.x=T)
global_data2 <- merge(global_data2,trait[,c(3,67:88)], by="Species",all.x=T)
global_data2$is_farmland <- as.factor(global_data2$Habitat=="Farmland")
global_data2$is_forest <- as.factor(global_data2$Habitat=="Forest")
global_data2 <- global_data2[abs(global_data2$slope) < 50,]
global_data2 <- merge(global_data2,country_data3, by.x="CountryGroup", by.y="country",all.x=T)

# Selecting data

data_pls2 <- global_data2[, c("slope","hico_2007","d_hico","for_2000","d_for","urb_2009","d_urb","temp_2000","d_temp","size_farm_2005","d_size_farm")]

data_pls2 <- data.frame(apply(data_pls2,2,Zscore))

# Initiate PLS

cv.modpls <- cv.plsR(data_pls2$slope,data_pls2[,-1],K=10,nt=10)
res.cv.modpls <- cvtable(summary(cv.modpls))
res1 <- plsR(data_pls2$slope,data_pls2[,-1], nt=10, typeVC="adaptative", pvals.expli=TRUE) # adaptative as NA in data
colSums(res1$pvalstep)

# Searching the best number of component to keep via CV

cv.modpls <- cv.plsR(slope~.,data=data_pls2,K=10,nt=10,NK=100)
res.cv.modpls <- cvtable(summary(cv.modpls))

# Using CV PRESS, 4 components must be kept

# Run PLS with 4 components
res <- plsR(slope~.,data=data_pls2,nt=4,pvals.expli=TRUE)
trend.bootYT1 <- bootpls(res,typeboot="fmodel_np",R=2000)
pls.ci <- confints.bootpls(trend.bootYT1,indices=2:ncol(data_pls2))
plots.confints.bootpls(pls.ci,typeIC="BCa",colIC=c("blue","blue","blue","blue"),legendpos ="topright")

# Using the empirical distribution of the best number of component, we can obtain an empircal measure of the significance of each effect
ind.BCa.YT1 <- (pls.ci[,7] < 0 & pls.ci[,8] < 0)|(pls.ci[,7] > 0 & pls.ci[,8] > 0) 
matind <- rbind(YT1=ind.BCa.YT1)
pi.e <- (prop.table(res.cv.modpls$CVPress)[4] %*% matind)/prop.table(res.cv.modpls$CVPress)[4]

# Plot PLS results
coef_plot <- data.frame(var=c("High input farm cover", "High input farm cover trend","Forest cover","Forest cover trend","Urbanisation cover","Urbanisation trend","Mean temperature","Temperature trend","Farm size","Farm size trend"),val=trend.bootYT1$t0[-1,1],
                      inf=pls.ci[,1],sup=pls.ci[,2],t(matind), sig=t(pi.e))
coef_plot$col_val <- "ns"
coef_plot$col_val[which(coef_plot$sig>=0.95 & coef_plot$val>0)] <- "pos"
coef_plot$col_val[which(coef_plot$sig>=0.95 & coef_plot$val<0)] <- "neg"
coef_plot$var <- as.character(coef_plot$var)
coef_plot$var <- factor(coef_plot$var, levels = c("Temperature trend","Mean temperature","Urbanisation trend","Urbanisation cover", 
                                                "Forest cover trend","Forest cover", "High input farm cover trend","High input farm cover", "Farm size", "Farm size trend"))
                                                
ggplot(coef_plot, aes(y=val, x=var)) +
  geom_rect(fill = "#DDF9DD", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, alpha = 0.1) +
  geom_rect(fill = "#F5A9A9", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, alpha = 0.1) +
  geom_bar(stat="identity", position=position_dodge(), alpha=0.7,aes(fill=var)) +
  scale_fill_manual(values=c("High input farm cover"="#D302F9","High input farm cover trend"="#D302F9",
                             "Urbanisation cover"="#196DF6","Urbanisation trend"="#196DF6",
                             "Forest cover"="#1BAE20","Forest cover trend"="#1BAE20",
                             "Mean temperature"="#FA0900","Temperature trend"="#FA0900", "Farm size" = "grey", "Farm size trend"="grey"))+
  geom_errorbar(aes(ymin=inf, ymax=sup), width=.2,position=position_dodge(.9)) +
  theme_modern() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_hline(yintercept=0, linetype="dashed", size=1) + 
  coord_flip()

```

### Trait syndrome

```{r}
# Merge species trends, traits and pressures

global_data3 <- merge(trend_species, sxi, by.x="Species", by.y="Name",all.x=T)
global_data3 <- merge(global_data3, sti, by.x="Species", by.y="SPECIES",all.x=T)
global_data3 <- merge(global_data3, pecbms_hab, by="Species",all.x=T)
global_data3 <- merge(global_data3, eunis_hab3, by="Species",all.x=T)
global_data3 <- merge(global_data3,ssi_eu, by="Species",all.x=T)
global_data3 <- merge(global_data3,trait[,c(3,67:90)], by="Species",all.x=T)
global_data3$is_farmland <- as.factor(global_data3$Habitat=="Farmland")
global_data3$is_forest <- as.factor(global_data3$Habitat=="Forest")
global_data3 <- global_data3[abs(global_data3$slope) < 50,]

# Selecting data

data_pls3 <- global_data3[, c("slope","is_farmland","is_forest","is_urban","STI","SSI","is_migrant","is_granivore","is_insectivore")]
data_pls3$is_farmland <- as.numeric(data_pls3$is_farmland)-1
data_pls3$is_forest <- as.numeric(data_pls3$is_forest)-1
data_pls3$is_urban <- as.numeric(data_pls3$is_urban)
data_pls3$is_urban[is.na(data_pls3$is_urban)] <- 0

data_pls3[,c(1,5,6)] <- data.frame(apply(data_pls3[,c(1,5,6)],2,Zscore))

# Initiate PLS

cv.modpls <- cv.plsR(data_pls3$slope,data_pls3[,-1],K=10,nt=10)
res.cv.modpls <- cvtable(summary(cv.modpls))
res1 <- plsR(data_pls3$slope,data_pls3[,-1], nt=10, typeVC="adaptative", pvals.expli=TRUE) # adaptative as NA in data
colSums(res1$pvalstep)

# Searching the best number of component to keep via CV

cv.modpls <- cv.plsR(slope~.,data=data_pls3,K=10,nt=10,NK=100)
res.cv.modpls <- cvtable(summary(cv.modpls))

# Using CV PRESS, 3 components must be kept

# Run PLS with 3 components
res <- plsR(slope~.,data=data_pls3,nt=3,pvals.expli=TRUE)
trend.bootYT1 <- bootpls(res,typeboot="fmodel_np",R=2000)
pls.ci <- confints.bootpls(trend.bootYT1,indices=2:ncol(data_pls3))
plots.confints.bootpls(pls.ci,typeIC="BCa",colIC=c("blue","blue","blue","blue"),legendpos ="topright")

# Using the empirical distribution of the best number of component, we can obtain an empircal measure of the significance of each effect
ind.BCa.YT1 <- (pls.ci[,7] < 0 & pls.ci[,8] < 0)|(pls.ci[,7] > 0 & pls.ci[,8] > 0) 
matind <- rbind(YT1=ind.BCa.YT1)
pi.e <- (prop.table(res.cv.modpls$CVPress)[3] %*% matind)/prop.table(res.cv.modpls$CVPress)[3]

# Plot PLS results
coef_plot <- data.frame(var=c("Farmland","Woodland","Urban","Species Temperature Index","Species Specialisation Index","Long distance migrant","Granivorous","Invertebrate-based diet"),val=trend.bootYT1$t0[-1,1],
                      inf=pls.ci[,1],sup=pls.ci[,2],t(matind), sig=t(pi.e))
coef_plot$col_val <- "ns"
coef_plot$col_val[which(coef_plot$sig>=0.95 & coef_plot$val>0)] <- "pos"
coef_plot$col_val[which(coef_plot$sig>=0.95 & coef_plot$val<0)] <- "neg"
coef_plot$var <- as.character(coef_plot$var)
coef_plot$var <- factor(coef_plot$var, levels = c("Farmland","Invertebrate-based diet","Urban","Long distance migrant","Granivorous","Woodland","Species Specialisation Index","Species Temperature Index"))
                                                
ggplot(coef_plot, aes(y=val, x=var)) +
  geom_rect(fill = "#DDF9DD", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, alpha = 0.1) +
  geom_rect(fill = "#F5A9A9", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, alpha = 0.1) +
  geom_bar(stat="identity", position=position_dodge(), alpha=0.7,fill="grey") +
  geom_errorbar(aes(ymin=inf, ymax=sup), width=.2,position=position_dodge(.9)) +
  theme_modern() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_hline(yintercept=0, linetype="dashed", size=1) + 
  coord_flip()
```

# References
Bogaart, P., van der Loo, M., Pannekoek, J., & Bogaart, M. P. (2016). Package ???rtrim???.

Brl??k, V., ??ilarov??, E., ??korpilov??, J. et al. Long-term and large-scale multispecies dataset tracking population changes of common European breeding birds. Sci Data 8, 21 (2021). https://doi.org/10.1038/s41597-021-00804-2

Devictor, V., Van Swaay, C., Brereton, T., Brotons, L., Chamberlain, D., Heli??l??, J., ... & Jiguet, F. (2012). Differences in the climatic debts of birds and butterflies at a continental scale. Nature climate change, 2(2), 121-124.

Godet, L., Ga??zere, P., Jiguet, F., & Devictor, V. (2015). Dissociating several forms of commonness in birds sheds new light on biotic homogenization. Global Ecology and Biogeography, 24(4), 416-426.

Le Viol, I., Jiguet, F., Brotons, L., Herrando, S., Lindstr??m, ??., Pearce-Higgins, J. W., ... & Devictor, V. (2012). More and more generalists: two decades of changes in the European avifauna. Biology letters, 8(5), 780-782.

Rigal, S., Devictor, V., & Dakos, V. (2020). A method for classifying and comparing non-linear trajectories of ecological variables. Ecological Indicators, 112, 106113.

Storchov??, L., & Ho????k, D. (2018). Life???history characteristics of European birds. Global Ecology and Biogeography, 27(4), 400-406.
