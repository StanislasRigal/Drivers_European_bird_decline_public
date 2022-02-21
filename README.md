# Farmland practices are driving bird populations decline across Europe.

R scripts and data for the following article: "Farmland practices are driving bird populations decline across Europe."

## R

The R scripts have been implemented on R version 3.4.4.

### Loading R packages

```{r setup, include=FALSE}

# Load packages

source("R_packages")

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

# estimating pop size from the percentage of total area of the country covered by the region #
for (i in 1:4){
  Country_subd[Country==subd[,reg][i],Count_min := Count_min*subd[,frac][i]]
  Country_subd[Country==subd[,reg][i],Count_max := Count_max*subd[,frac][i]]
}
for (i in 1:680){
  Country_subd[i,estimate := round((geoMean(c(Count_min,Count_max),na.rm=T))*2)]
}
# adding Belgium and Germany regions estimates to the main abundance data.frame #
abd <- rbind(abd,Country_subd)
setorder(abd,Species,Country)
# removing unrealised species-country combinations #
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

# index and abundance
df_pop<-droplevels(subset(df, Year %in% c(1980:2016)))
df_pop<-droplevels(subset(df_pop, Species %in% levels(droplevels(subset(df_pop, start_year<=1981))$Species)))
df_pop<-droplevels(df_pop[!which(df_pop$Species=="Passer domesticus" & df_pop$CountryGroup %in% levels(df_pop$CountryGroup)[23]),])
```

# Combining species indices and abundance at the European level: RTIM

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
# Create data files for each species cf. README_TRIM_SHELL.docx
for (i in S){
  subset <- df_pop[Species==i]
  subset <- droplevels(subset)
  c <- levels(subset[, CountryGroup]) # possible need to exclude aberrant time series? cf. list in "10_Linear-regressions.R"
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
make_Indices_TT_file <- function(result) {
# This function creates a dataframe which contains indices and time totals.
  
 indices_TT_file <- data.frame(Year = integer(result$nyear), 
             # Years for which indices and time totals have been calculated.
             TT_model = numeric(result$nyear), 
             # Time totals based on the estimated model.
             TT_model_SE = numeric (result$nyear), 
             # Standard errors of model time totals.
             TT_imputed = numeric(result$nyear), 
             # Imputed time totals.
             TT_imputed_SE = numeric(result$nyear), 
             # Standard errors of imputed time totals.
             Index_model = numeric(result$nyear), 
             # Indices based on the estimated model.
             Index_model_SE = numeric(result$nyear), 
             # Standard errors of model indices.
             Index_imputed = numeric(result$nyear), 
             # Indices based on the imputed time totals.
             Index_imputed_SE = numeric(result$nyear))
             # Standard errors of imputed indices.
  return(indices_TT_file)
}

make_arg_output_file <- function(arguments) {
# This function creates a dataframe with slopes for the entire period and a subperiod.
# Both additive and multiplicative slopes are included.
  arg_output_file <- data.frame(N_sites = integer(1), 
                 # Number of unique sites.
                 N_time_values = integer(1), 
                 # Number of unique years.
                 N_observed_zero_counts = integer (1), 
                 # Number of zero counts.
                 N_observed_positive_counts = integer(1), 
                 # Number of positive counts.
                 N_missing_Counts = integer(1), 
                 # Number of missing counts.
                 N_counts_Total = integer(1), 
                 # Total number of counts.
                 Base_year_first_year = numeric(1),
                 # Calendar year used as base year for indices. 
                 # If base_year_first_year equals Base_year_last_year a single year is used as base year
                 # If Base_year_first_year < Base_year_last_year, a period is used as base time
                 # In the latter case, Base_year_first_year is the first year of the period.
                 Base_year_last_year = numeric(1),
                 # Calendar year used as base year for indices. 
                 # If base_year_first_year equals Base_year_last_year a single year is used as base year
                 # If Base_year_first_year < Base_year_last_year, a period is used as base time
                 # In the latter case, Base_year_last_year is the last year of the period.
                 Changepoints = numeric(1),
                 # Changepoints used.
                 Overdispersion = numeric(1),
                 # Estimated overdispersion.
                 Serial_correlation = numeric(1),
                 # Estimated serial correlation.
                 Slope_model_mul = numeric(1),
                 # Multiplicative model slope, calculated for the entire period.
                 Slope_model_mul_SE = numeric(1),
                 # Standard error of multiplicative model slope for the entire period.
                 Slope_model_classification = character(1),
                 # Trend classification of model slope for the entire period.
                 Slope_imputed_mul = numeric(1),
                 # Multiplicative imputed slope for the entire period.
                 Slope_imputed_mul_SE = numeric(1),
                 # Standard error of multiplicative imputed slope for the entire period.
                 Slope_imputed_classification = character(1),
                 # Trend classification of multiplicative imputed slope for the entire period.
                 Slope_model_add = numeric(1),
                 # Additive model slope for the entire period.
                 Slope_model_add_SE = numeric(1),
                 # Standard error of additive model slope for the entire period.
                 Slope_imputed_add = numeric(1),
                 # Additive imputed slope for the entire period.
                 Slope_imputed_add_SE = numeric(1),
                 # Standard error of additive imputed slope for the entire period.
                 Slope_from = numeric(1),
                 # First year of the subperiod over which a slope has been calculated. 
                 Slope_to = numeric(1),
                 # Last year of the subperiod over which a slope has been calculated.
                 # Usually this year is equal to the last year of the entire period.
                 Slope_from_model_mul = numeric(1),
                 # Multiplicative model slope for the subperiod.
                 Slope_from_model_mul_SE = numeric(1),
                 # Standard error of multiplicative model slope for the subperiod.
                 Slope_from_model_classification = character(1),
                 # Trend classification of model slope for the subperiod.
                 Slope_from_imputed_mul = numeric(1),
                 # Multiplicative imputed slope for the subperiod.
                 Slope_from_imputed_mul_SE = numeric(1),
                 # Standard error of multiplicative imputed slope for the subperiod.
                 Slope_from_imputed_classification = character(1),
                 # Trend classification of multiplicative imputed slope for the subperiod.
                 Slope_from_model_add = numeric(1),
                 # Additive model slope for the subperiod.
                 Slope_from_model_add_SE = numeric(1),
                 # Standard error of the additive model slope for the subperiod.
                 Slope_from_imputed_add = numeric(1),
                 # Additive imputed slope for the subperiod.
                 Slope_from_imputed_add_SE = numeric(1)
                 # Standard error of additive imputed slope for the subperiod.
             )
    arg_output_file$Slope_from <- arguments$Slope_from
    # A value is assigned to Slope_from. 
    # This to enable to make the function that fills the arg_output file, as simple as possible.
    # The variabele 'arguments' is not needed due to this, when calling the function.
    # Attention: 
    # It is assumed that the first year of the subperiod is closer to the present than the first year of the entire period.
    # It is also assumed that the last year of the subperiod is equal to the last year of the entire period. 
    arg_output_file$Base_year_first_year <- arguments$Base_year_first_year
    arg_output_file$Base_year_last_year <- arguments$Base_year_last_year
  return(arg_output_file)
}

fill_Indices_TT_file <- function(indices_TT_file, result, arguments) {
    
    indices <- index(result, which = "both", covars = FALSE, base = arguments$Base_year_first_year:arguments$Base_year_last_year)
    time_totals <- totals(result, which = "both")
    indices_TT_file$Year <- indices$time
    # The years over which the indices and time totals have been calculated.
    indices_TT_file$TT_model <- time_totals$fitted 
    # Time totals based on the estimated model. 
    indices_TT_file$TT_model_SE <- time_totals$se_fit
    # Standard errors of model time totals.
    indices_TT_file$TT_imputed <- time_totals$imputed
    # Imputed time totals.
    indices_TT_file$TT_imputed_SE <- time_totals$se_imp
    # Standard errors of imputed time totals.
    indices_TT_file$Index_model <- indices$fitted
    # Indices based on the estimated model.
    indices_TT_file$Index_model_SE <- indices$se_fit
    # Standard errors of model indices.
    indices_TT_file$Index_imputed <- indices$imputed
    # Indices based on the imputed time totals.
    indices_TT_file$Index_imputed_SE <- indices$se_imp
    # Standard errors of imputed indices.
    return(indices_TT_file)
  }
  
fill_arg_output_file <- function(arg_output_file, result, counts) {
    overviewCounts <- count_summary(counts, count_col = "count", site_col = "site", year_col = "time")
    
    slopes_imputed <- overall(result, which = "imputed")
    slopes_fitted <- overall(result, which = "fitted")
    # The slopes over the entire period. Both variables are lists (specific class of R-objects).
    
    slopes_deelperiode_imputed <- overall(result, which = "imputed", changepoints = c(arg_output_file$Slope_from))
    slopes_deelperiode_fitted <- overall(result, which = "fitted", changepoints = c(arg_output_file$Slope_from))
    # The slopes over a subperiod. Both variables are lists.
    # Attention: it is assumed that the first year of the subperiod is closer to the present year than the first year of the entire period.
    # The last year of a subperiod is equal to the last year of the entire period.
    arg_output_file$N_sites <- result$nsite
    # Number of unique sites in the counts file used to call rtrim.
    arg_output_file$N_time_values <- result$ntime
    # Number of unique years in the counts file used to call rtrim.
    arg_output_file$N_observed_zero_counts <- overviewCounts$zero_counts
    # Number of zero counts, in the entire periode, in the counts file used to call rtrim.
    arg_output_file$N_observed_positive_counts <- overviewCounts$positive_counts
    # Number of positive counts, in the entire periode, in the counts file used to call rtrim.
    arg_output_file$N_missing_counts <- overviewCounts$missing_counts
    # Number of missing counts, in the entire period, in the counts file used to call rtrim.
    arg_output_file$N_counts_total <- overviewCounts$total_counts
    # Total number of counts, in the entire period, in the counts file used to call rtrim.
    arg_output_file$Changepoints <- paste(result$changepoints, collapse = ", ")
    # Changepoints used to calculate the indices (years that were not important have been left out).
    arg_output_file$Overdispersion <- overdispersion(result)
    # Calculated overdispersion, over the entire period.
    arg_output_file$Serial_correlation <- serial_correlation(result)
    # Calculated serial correlation, over the entire period.
    arg_output_file$Slope_model_mul <- slopes_fitted$slope$mul
    # Multiplicative model slope, calculated for the entire period.
    arg_output_file$Slope_model_mul_SE <- slopes_fitted$slope$se_mul
    # Standard error of multiplicative model slope for the entire period.
    arg_output_file$Slope_model_classification <- slopes_fitted$slope$meaning
    # # Trend classification of model slope for the entire period.
    arg_output_file$Slope_imputed_mul <- slopes_imputed$slope$mul
    # Multiplicative imputed slope for the entire period.
    arg_output_file$Slope_imputed_mul_SE <- slopes_imputed$slope$se_mul
    # Standard error of multiplicative imputed slope for the entire period.
    arg_output_file$Slope_imputed_classification <- slopes_imputed$slope$meaning
    # Trend classification of multiplicative imputed slope for the entire period.
    arg_output_file$Slope_model_add <- slopes_fitted$slope$add
    # Additive model slope for the entire period.
    arg_output_file$Slope_model_add_SE <- slopes_fitted$slope$se_add
    # Standard error of additive model slope for the entire period.
    arg_output_file$Slope_imputed_add <-slopes_imputed$slope$add
    # Additive imputed slope for the entire period.
    arg_output_file$Slope_imputed_add_SE <- slopes_imputed$slope$add
    # Standard error of additive imputed slope for the entire period.
    arg_output_file$Slope_from_model_mul <- slopes_deelperiode_fitted$slope$mul[slopes_deelperiode_fitted$slope$from == arg_output_file$Slope_from]
    # Multiplicative model slope for the subperiod.
    # Attention: it is assumed that a subperiod has a different first year than the entire period.
    # The first year of the subperiod is closer to the present than the first year of the entire period.
    # The last year of a subperiod is the same year as the last year of the entire period. 
    # Attention: the last year with counts cannot be a changepoint!
    arg_output_file$Slope_from_model_mul_SE <- slopes_deelperiode_fitted$slope$se_mul[slopes_deelperiode_fitted$slope$from == arg_output_file$Slope_from]
    # Standard error of multiplicative model slope for the subperiod.
    arg_output_file$Slope_from_model_classification <- slopes_deelperiode_fitted$slope$meaning[slopes_deelperiode_fitted$slope$from == arg_output_file$Slope_from]
    # Trend classification of model slope for the subperiod.
    arg_output_file$Slope_from_imputed_mul <- slopes_deelperiode_imputed$slope$mul[slopes_deelperiode_imputed$slope$from == arg_output_file$Slope_from]
    # Multiplicative imputed slope for the subperiod.
    arg_output_file$Slope_from_imputed_mul_SE <- slopes_deelperiode_imputed$slope$se_mul[slopes_deelperiode_imputed$slope$from == arg_output_file$Slope_from]
    # Standard error of multiplicative imputed slope for the subperiod.
    arg_output_file$Slope_from_imputed_classification <- slopes_deelperiode_imputed$slope$meaning[slopes_deelperiode_imputed$slope$from == arg_output_file$Slope_from]
    # Trend classification of imputed slope for the subperiod.
    arg_output_file$Slope_from_model_add <- slopes_deelperiode_fitted$slope$add[slopes_deelperiode_fitted$slope$from == arg_output_file$Slope_from]
    # Additive model slope for the subperiod.
    arg_output_file$Slope_from_model_add_SE <- slopes_deelperiode_fitted$slope$se_add[slopes_deelperiode_fitted$slope$from == arg_output_file$Slope_from]
    # Standard error of additive model slope for the subperiod.
    arg_output_file$Slope_from_imputed_add <- slopes_deelperiode_imputed$slope$add[slopes_deelperiode_imputed$slope$from == arg_output_file$Slope_from]
    # Additive imputed slope for the subperiod.
    arg_output_file$Slope_from_imputed_add_SE <- slopes_deelperiode_imputed$slope$se_add[slopes_deelperiode_imputed$slope$from == arg_output_file$Slope_from]
    # Standard error of additive imputed slope for the subperiod.
    return(arg_output_file)
}
  
makeOverview <- function(listSpeciesStratumCombinations) {
  
  numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)
  overview <- data.frame(ss_combinations = character(numberSpeciesStratumCombinations), 
                            # The combinations of species and stratum in the working directory.
                            species_group = character(numberSpeciesStratumCombinations),
                            # Short for specific species group.
                            species_number = integer(numberSpeciesStratumCombinations),
                            # Unique number for a species.
                            first_year = integer(numberSpeciesStratumCombinations),
                            # The first year of which counts are available for this specific combination of species and stratum.
                            last_year = integer(numberSpeciesStratumCombinations),
                            # The last year of which counts are available for this specific combination of species and stratum.
                            stratum_number = integer(numberSpeciesStratumCombinations),
                            # Unique number for a stratum.
                            success = character(numberSpeciesStratumCombinations), 
                            # Shows whether the analysis was analysed successful, or not. 
                            # yes for a successful analysis, no for an unsuccessful analysis. 
                            attempt_1 = character(numberSpeciesStratumCombinations), 
                            # Outcome of the first attempt. The attempt was successful or not.
                            attempt_2 = character (numberSpeciesStratumCombinations), 
                            # Outcome of the second attempt. 
                            # The attempt has not been executed (because the first attempt was already successful), successful or not successful. 
                            # When the first attempt was successful, the value that has allready been assigned (n.a., not applicable) will not be altered.  
                            attempt_3 = character(numberSpeciesStratumCombinations), 
                            # Outcome of the third attempt. 
                            # The attempt has not been executed (because the second attempt was successful), successful or not successful.
                            attempt_4 = character(numberSpeciesStratumCombinations),
                            # Outcome of the last attempt. 
                            # The attempt has not been executed (because the third attempt was successful), successful or not successful.
                            error_1 = character(numberSpeciesStratumCombinations),
                            # Error message of the first attempt.
                            error_2 = character(numberSpeciesStratumCombinations),
                            # Error message of the second attempt.
                            error_3 = character(numberSpeciesStratumCombinations),
                            # Error message of the third attempt.
                            error_4 = character(numberSpeciesStratumCombinations))
                            # Error message of the last attempt.
  overview$ss_combinations <- gsub(listSpeciesStratumCombinations, pattern = "_arg_input_stratum.csv", replacement = "")
  overview$success <- "no"
  overview$attempt_1 <- "n.a."                                                                                       
  overview$attempt_2 <- "n.a."
  overview$attempt_3 <- "n.a."
  overview$attempt_4 <- "n.a."
  overview$error_1 <- "n.a."
  overview$error_2 <- "n.a."
  overview$error_3 <- "n.a."
  overview$error_4 <- "n.a."
  overview$species_group <- gsub(listSpeciesStratumCombinations, pattern = "_[0-9]+_[0-9]+_arg_input_stratum.csv", replacement = "")
  
  zonderSoortgroep <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_", replacement = "")
  
  overview$species_number <- as.integer(gsub(zonderSoortgroep, pattern = "_[0-9]_arg_input_stratum.csv", replacement = ""))
  
  zonderSoortgroepZonderSoortnummer <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_[0-9]+_", replacement = "")
  
  overview$stratum_number <- gsub(zonderSoortgroepZonderSoortnummer, pattern = "_arg_input_stratum.csv", replacement = "")
  
  overview <- overview[order(overview$species_number, overview$stratum_number), ]
  return(overview)
}


make_All_Indices_All_Trends <- function(overview, listSpeciesStratumCombinations){
  
  numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)
  first_year_over_all_species <- min(overview$first_year)
  last_year_over_all_species <- max(overview$last_year)
  
  complete_period <- first_year_over_all_species:last_year_over_all_species
  columnnames <- as.character(complete_period)
  number_of_columns <- length(complete_period)
  number_of_rows <- numberSpeciesStratumCombinations * 4
  # For rows for each combination of species and stratum (one row for each record type)
  temporary_matrix <- matrix(nrow = number_of_rows, ncol = number_of_columns)
  temporary_dataframe <- as.data.frame(temporary_matrix)
  colnames(temporary_dataframe) <- columnnames
  
  all_Indices_All_Trends <- data.frame(year_of_analysis = integer(numberSpeciesStratumCombinations * 4),
    # Dataframe, names all_Indices_All_Trends, is created.
    # First column (year_of_analysis): Most recent year of which counts are available.
    Species_number = integer(numberSpeciesStratumCombinations),
    # Unique number for a species.
    Stratum_number = integer(numberSpeciesStratumCombinations),
    # Unique number for a stratum.
    Recordtype_number = integer(numberSpeciesStratumCombinations),
    # Unique number for record types. There are four record types. 
    # 1: indices 
    # 2: standard errors of indices
    # 3: time totals
    # 4: standard errors of time totals
    Recordtype_name = character(numberSpeciesStratumCombinations),	
    # Name of the record type.
    N_sites = integer(numberSpeciesStratumCombinations),	
    # Number of unique sites in the counts file used to call rtrim.
    Slope_imputed_mul =	numeric(numberSpeciesStratumCombinations),
    # Multiplicative imputed slope for the entire period.
    Slope_imputed_mul_SE = numeric(numberSpeciesStratumCombinations),	
    # Standard error of the multiplicative imputed slope for the entire period.
    Slope_imputed_classification = character(numberSpeciesStratumCombinations),	
    # Trend classification of the imputed slope for the entire period.
    Slope_from_imputed_mul =	numeric(numberSpeciesStratumCombinations),
    # Multiplicative imputed slope for the subperiod.
    Slope_from_imputed_mul_SE = numeric(numberSpeciesStratumCombinations),	
    # Standard error of the multiplicative imputed slope for the subperiod.
    Slope_from_imputed_classification = character(numberSpeciesStratumCombinations),	
    # Trend classification of the imputed slope for the subperiod.
    Slope_from = integer(numberSpeciesStratumCombinations),	
    # First year of the subperiod over which a slope has been calculated. 
    Date_analysis = character(numberSpeciesStratumCombinations))
    # Date at which the analysis was performed.
  without_Species_Group <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_", replacement = "")
  
  species_Code <- as.integer(gsub(without_Species_Group, pattern = "_[0-9]_arg_input_stratum.csv", replacement = ""))
  
  without_Species_Group_Without_Species_Number <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_[0-9]+_", replacement = "")
  
  stratum_Number <- gsub(without_Species_Group_Without_Species_Number, pattern = "_arg_input_stratum.csv", replacement = "")
  
  all_Indices_All_Trends$Recordtype_number <- rep(1:4, each = numberSpeciesStratumCombinations)
  # Unique number for record types. There are four record types. 
  # 1: indices 
  # 2: standard errors of indices
  # 3: time totals
  # 4: standard errors of time totals
  number_of_record_types <- length(unique(all_Indices_All_Trends$Recordtype_number))
  
  all_Indices_All_Trends$Species_number <- rep(species_Code, number_of_record_types)
  all_Indices_All_Trends$Stratum_number <- rep(stratum_Number, number_of_record_types)
  
  names_record_types <- c("indices", "se_indices", "time_totals", "se_time_totals")
  all_Indices_All_Trends$Recordtype_name <- rep(names_record_types, each = numberSpeciesStratumCombinations)	
  # Name of the record type.
  all_Indices_All_Trends$Slope_imputed_classification <- ""
  all_Indices_All_Trends$Slope_from_imputed_classification <- ""
  all_Indices_All_Trends$Date_analysis <- format(Sys.Date(), "%d-%m-%Y")
  # Date at which the analysis was done.
  all_Indices_All_Trends <- cbind(all_Indices_All_Trends, temporary_dataframe)
  
return(all_Indices_All_Trends)  
}

fill_All_Indices_All_Trends <- function(result, arguments, j, listSpeciesStratumCombinations, all_Indices_All_Trends){
  
  numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)
  
  indices <- index(result, which = "both", covars = FALSE, base = arguments$Base_year_first_year:arguments$Base_year_last_year)
  time_totals <- totals(result, which = "both")
  slopes_imputed <- overall(result, which = "imputed")
  slopes_deelperiode_imputed <- overall(result, which = "imputed", changepoints = c(arguments$Slope_from))
  
  names_columns <- as.character(indices$time)
  position_colums <- colnames(all_Indices_All_Trends) %in% names_columns
  all_Indices_All_Trends$year_of_analysis <- max(indices$time)
  # Most recent year of which counts are available.
  all_Indices_All_Trends$N_sites[j] <- result$nsite
  all_Indices_All_Trends$N_sites[j + numberSpeciesStratumCombinations] <- result$nsite
  all_Indices_All_Trends$N_sites[j + numberSpeciesStratumCombinations * 2] <- result$nsite
  all_Indices_All_Trends$N_sites[j + numberSpeciesStratumCombinations * 3] <- result$nsite
  # Number of unique sites in the counts file used to call rtrim.
  all_Indices_All_Trends$Slope_imputed_mul[j] <- slopes_imputed$slope$mul
  # Multiplicative imputed slope for the entire period.
  all_Indices_All_Trends$Slope_imputed_mul_SE[j] <-	slopes_imputed$slope$se_mul
  # Standard error of the multiplicative imputed slope for the entire period.
  all_Indices_All_Trends$Slope_imputed_classification[j] <-	slopes_imputed$slope$meaning
  # Trend classification of the imputed slope for the entire period.
  all_Indices_All_Trends$Slope_from_imputed_mul[j] <- slopes_deelperiode_imputed$slope$mul[slopes_deelperiode_imputed$slope$from == arguments$Slope_from]
  # Multiplicative imputed slope for the subperiod.
  all_Indices_All_Trends$Slope_from_imputed_mul_SE[j] <- slopes_deelperiode_imputed$slope$se_mul[slopes_deelperiode_imputed$slope$from == arguments$Slope_from]	
  # Standard error of the multiplicative imputed slope for the subperiod.
  all_Indices_All_Trends$Slope_from_imputed_classification[j] <- slopes_deelperiode_imputed$slope$meaning[slopes_deelperiode_imputed$slope$from == arguments$Slope_from]	
  # Trend classification of the imputed slope for the subperiod.
  all_Indices_All_Trends$Slope_from[j] <-	arguments$Slope_from
  # First year of the subperiod over which a slope has been calculated. 
  all_Indices_All_Trends[j, position_colums] <- round(100 * indices$imputed, 1)
  all_Indices_All_Trends[j + numberSpeciesStratumCombinations, position_colums] <- round(indices$se_imp, 2)
  all_Indices_All_Trends[j + numberSpeciesStratumCombinations * 2, position_colums] <- round(time_totals$imputed, 2)
  all_Indices_All_Trends[j + numberSpeciesStratumCombinations * 3, position_colums] <- round(time_totals$se_imp, 2)
  return(all_Indices_All_Trends)  
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

# gather imputed COUNTRYWIDE species' ABUNDANCE time series into one file #
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

# Estiamting abundance dynamic of the avifauna

### Functions accounting for non-linearity and sampling error

See details in Rigal et al. (2020).
```{r}

class.trajectory <- function (Y = NULL, X = NULL, dataset = NULL, interval_size = 0.5)
{
  if (is.null(Y) == TRUE & is.null(Y) == TRUE & is.null(dataset) == TRUE){
    stop("either 'dataset' or at least 'Y' and 'X' must be specified")
  }
  if (is.null(Y) == TRUE & is.null(Y) == TRUE) {
    Y <- dataset[,1]
    X <- dataset[,2]
  }else{
    if (class(Y) == "character" & class(X) == "character") {
      if (is.null(dataset) == TRUE) {
        stop("if 'Y' and 'X' are character, 'dataset' must exist")
      }else{
        Y <- dataset[, Y]
        X <- dataset[, X]
      }
    }else{
      if (!(class(Y) %in% c("numeric","integer")) == TRUE & !(class(X) %in% c("numeric","integer")) == TRUE) {stop("'Y' and 'X' must be either characters or vector but 'class' must be similar")}
    }
  }
  
  data <- data.frame(cbind(Y, X))
  data <- data[order(data$X),]                                                                      # ordering the X values
  
  if (length(X)<4){
    stop("time series length must be at least 4")
  }
  
  Y <- data$Y
  X <- data$X
  
  linear.model <- lm(Y~X)
  
  orthogonal_polynomial <- lm(Y~poly(X,2, raw=F))                                                   # After getting Y = gamma*chi + delta*X' + epsilon with orthogonal polynomial
  # we have to perform a variable change to obtain relevant values in the X interval 
  # for first_order_coefficient, second_order_coefficient and intercept,
  # knowing that X'= alpha*X + beta 
  # and chi = eta*X'^2 + theta
  
  gammab  <-  orthogonal_polynomial$coefficients[3]
  delta  <-  orthogonal_polynomial$coefficients[2]
  epsilon  <-  orthogonal_polynomial$coefficients[1]
  
  alpha  <-  lm(orthogonal_polynomial$model[, 2][, 1]~X)$coef[2]
  beta  <-  lm(orthogonal_polynomial$model[, 2][, 1]~X)$coef[1]
  
  eta  <-  1/lm((orthogonal_polynomial$model[, 2][, 1])^2~orthogonal_polynomial$model[, 2][, 2])$coef[2]
  theta  <-  (-lm((orthogonal_polynomial$model[, 2][, 1])^2~orthogonal_polynomial$model[, 2][, 2])$coef[1])*eta
  
  Y2<-Y*(max(X)-min(X))/(max(Y)-min(Y))                                                             # p2 and p3 are relevant when Y and X amplitudes are equivalent,
  # in particular when studying scaled-to-1 indices, Y and X amplitudes
  # may be very different, so we scaled the amplitudes to calculate p2 and p3 
  polynomial_orthonormal_basis<-lm(Y2~poly(X,2, raw=T))$coefficients
  
  if(summary(orthogonal_polynomial)$coefficients[3, 4] <= 0.05){                                     # non linear case
    classification <- data.frame(first_order_coefficient = (delta+2*beta*gammab*eta)*alpha,
                                 first_order_pvalue = summary(orthogonal_polynomial)$coefficients[2, 4],
                                 second_order_coefficient = (alpha^2)*gammab*eta,
                                 second_order_pvalue = summary(orthogonal_polynomial)$coefficients[3, 4],
                                 strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
                                 intercept = epsilon+beta*delta+(beta^2)*gammab*eta+gammab*theta,
                                 x_m = (X[length(X)]-X[1])/2+X[1],
                                 p1 = -(delta+2*beta*gammab*eta)/(2*alpha*gammab*eta),                    # points of interest
                                 p2 = (-polynomial_orthonormal_basis[2]+1)/(2*polynomial_orthonormal_basis[3]),
                                 p3 = (-polynomial_orthonormal_basis[2]-1)/(2*polynomial_orthonormal_basis[3]))
  }else{                                                                                            # linear case
    classification <- data.frame(first_order_coefficient = delta*alpha,
                                 first_order_pvalue = summary(orthogonal_polynomial)$coefficients[2, 4],
                                 second_order_coefficient = 0,
                                 second_order_pvalue = summary(orthogonal_polynomial)$coefficients[3, 4],
                                 strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
                                 intercept = epsilon+delta*beta,
                                 x_m = (X[length(X)]-X[1])/2+X[1],
                                 p1 = NA,
                                 p2 = NA,
                                 p3 = NA)
  }
  
  classification$r.sq <- summary(orthogonal_polynomial)$adj.r.squared                                # retrieve the adjusted coefficient of determination
  
  # compute the derivaive at xm-delta and at xm + delta with delta being half of the input interval size
  derivative  <-  2*(classification$x_m-(X[length(X)]-X[1])*(interval_size/2))*classification$second_order_coefficient+classification$first_order_coefficient
  derivative2  <-  2*(classification$x_m+(X[length(X)]-X[1])*(interval_size/2))*classification$second_order_coefficient+classification$first_order_coefficient
  
  
  if(sign(derivative) != sign(derivative2)){                                                        # non consistent direction around x_m
    classification$derivative  <-  NA
    classification$intercept_derivative  <-  NA
  }else{                                                                                            # consistent direction around x_m
    classification$derivative  <-  mean(c(derivative, derivative2))
    classification$intercept_derivative  <-  (classification$second_order_coefficient*classification$x_m^2+classification$first_order_coefficient*classification$x_m+classification$intercept)-classification$x_m*classification$derivative
  }
  
  # compute the derivative of the curvature function
  classification$derivated_curvature  <-  -12*(classification$second_order_coefficient^2)*(2*classification$second_order_coefficient*classification$x_m+classification$first_order_coefficient)*(classification$second_order_coefficient/abs(classification$second_order_coefficient))/
    ((1+(2*classification$second_order_coefficient*classification$x_m+classification$first_order_coefficient)^2)^(2.5))
  
  if(classification$second_order_pvalue>0.05){classification$derivated_curvature <- NA}
  
  classification$direction <- NA                                                                    # classify the direction
  classification$direction[which(classification$derivative > 0)] <- "increase"
  classification$direction[which(classification$derivative < 0)] <- "decrease"
  classification$direction[which(is.na(classification$derivative))] <- "stable"
  classification$direction[which(as.numeric(classification$first_order_pvalue)>0.05 & as.numeric(classification$second_order_pvalue)>0.05)] <- "stable"
  
  classification$acceleration <- NA                                                                 # classify the acceleration
  classification$acceleration[which(classification$derivated_curvature < 0)] <- "accelerated"
  classification$acceleration[which(classification$derivated_curvature > 0)] <- "decelerated"
  classification$acceleration[which(classification$direction == "stable" &
                                      classification$second_order_coefficient < 0)] <- "concave"
  classification$acceleration[which(classification$direction == "stable" &
                                      classification$second_order_coefficient > 0)] <- "convex"
  classification$acceleration[which(is.na(classification$derivated_curvature))] <- "constant"
  
  classification$shape_class <- paste(classification$direction,                                       # give the final classification combining direction and acceleration
                                      classification$acceleration,
                                      sep="_")
  
  linear.model.summary <- summary(linear.model)                                                       # provide the linear approach results for comparison
  
  classification$linear_slope <- linear.model.summary$coefficients[2, 1]
  classification$linear_slope_pvalue <- linear.model.summary$coefficients[2, 4]
  classification$linear_intercept <- linear.model.summary$coefficients[1, 1]
  
  classification$first_X_value <- X[1]
  classification$last_X_value <- X[length(X)]
  
  row.names(classification) <- "Y"
  
  return(classification)
}

mc_trend <- function(dataset,         # data
                     niter,           # number of MC simulations
                     ref_year=NULL,   # reference year
                     mid="mid",       # by defaultif ref_year is missing, equal to the mid year of the interval or to first year
                     correction=TRUE) # set the reference value to 100 and correct values below 0 before logtransformation
{
  
  b <- data.frame(t(rep(NA, 13)))
  attributes(b)$names <- c("second_order_coefficient",
                           "first_order_coefficient",
                           "strd_error",
                           "shape_class",
                           "intercept",
                           "p_1",
                           "p_2",
                           "p_3",
                           "second_order_pvalue",
                           "slope_p_value",
                           "slope",
                           "r.sq",
                           "ref_year")
  
  if(is.null(ref_year)){
    if(mid=="mid"){
      ref_year <- dataset$Year[round(nrow(dataset)/2)+1]
    }
    if(mid=="first"){
      ref_year <- min(dataset$Year)
    }
  }
  
  if(correction == TRUE){
    ref_value <- dataset$Index[dataset$Year == ref_year]
    if(ref_value == 1){
      dataset$Index <- 100*dataset$Index             # set reference year value to 100
      dataset$Index[which(dataset$Index <= 1)] <- 1    # set values < 1 to 1
      dataset$Index_SE[which(dataset$Index <= 1)] <- 0 # and their SE to 0
      dataset$Index_SE <- 100*dataset$Index_SE
    }
    if(ref_value!=1 & ref_value!=100){
      if(ref_value>1){
        dataset$Index <- dataset$Index/ref_value
        dataset$Index_SE <- dataset$Index_SE/ref_value
      }
      if(ref_value<1){
        stop("use 'correction = FALSE' when value of the reference year is strictly below 1")
      }
      dataset$Index <- 100*dataset$Index
      if(length(which(dataset$Index <= 1))>0){print("caution, low values corrected, if strongly decreasing or increasing trajectory, use respectively  first or last year as referential")}
      dataset$Index[which(dataset$Index <= 1)] <- 1
      dataset$Index_SE[which(dataset$Index <= 1)] <- 0
      dataset$Index_SE <- 100*dataset$Index_SE
    }
    if(ref_value == 100){
      dataset$Index[which(dataset$Index <= 1)] <- 1
      dataset$Index_SE[which(dataset$Index <= 1)] <- 0
    }
    
    dataset$sd <- dataset$Index_SE/dataset$Index                               # SE of log transformed data
    dataset$log <- log(dataset$Index)                                          # log transforme Y
    for(j in 1:nrow(dataset)){
      if(dataset$sd[j]>(dataset$log[j])){dataset$sd[j] <- (dataset$log[j])}    # set SE to value amplitude if SE > value (if not, it leads to huge values when resampling in the next loop)
    }
  }
  
  if(correction == FALSE){
    if(min(dataset$Index)<0){
      min_value <- abs(min(dataset$Index))
      dataset$Index <- dataset$Index+min_value
    }else{min_value <- 0}
    dataset$sd <- dataset$Index_SE
    dataset$log <- dataset$Index
  }
  
  for(i in 1:niter){
    
    a <- rnorm(nrow(dataset), mean=dataset$log, sd=dataset$sd)                 # simulate Y values from normal distribution (mean= original Y value, sd = original SE)
    if(correction == TRUE){
      a <- exp(a)/exp(a[which(dataset$Year == ref_year)])*100                    # set reference year value to 100 and retransform values if logtranformed
    }else{a <- a-min_value}
    
    a <- class.trajectory(a, dataset$Year)
    b[i, 1] <- a$second_order_coefficient
    b[i, 2] <- a$first_order_coefficient
    b[i, 3] <- a$strd_error
    b[i, 4] <- a$shape_class
    b[i, 5] <- a$intercept
    if(a$second_order_coefficient!=0){
      if(findInterval(a$p1,  c(min(dataset$Year), max(dataset$Year))) == 1){ # record changing point inside time series
        b[i, 6] <- a$p1}else{b[i, 6] <- NA}
      if(findInterval(a$p2,  c(min(dataset$Year), max(dataset$Year))) == 1){
        b[i, 7] <- a$p2}else{b[i, 7] <- NA}
      if(findInterval(a$p3,  c(min(dataset$Year), max(dataset$Year))) == 1){
        b[i, 8] <- a$p3}else{b[i, 8] <- NA}
    }else{
      b[i, 6] <- NA
      b[i, 7] <- NA
      b[i, 8] <- NA
    }
    b[i, 9] <- a$second_order_pvalue
    b[i, 10] <- a$linear_slope_pvalue
    b[i, 11] <- a$linear_slope
    b[i, 12] <- max(0,a$r.sq)
  }
  b[, 4] <- as.factor(b[, 4])
  b[, 13] <- rep(ref_year, nrow(b))
  return(b)
}

res_trend<-function(dataset,
                    niter,
                    ref_year=NULL,
                    mid="mid",
                    correction=TRUE){
  
  if(nrow(dataset)>3 & anyNA(dataset$Index_SE) == FALSE){
    
    simulated <- mc_trend(dataset, niter, ref_year, mid, correction)
    max_max<-NULL
    
    if(length(levels(simulated$shape_class))>1){
      test<-multinomial.theo.multcomp(simulated$shape_class, p = rep(1/length(levels(simulated$shape_class)),
                                                             length(levels(simulated$shape_class))), prop=TRUE)
      if(min(test$p.value2[test$observed>test$expected])<0.05){
        max_shape <- row.names(test$p.value)[which(test$observed == max(test$observed[test$observed>test$expected]))]
        if(max_shape %in% c("increase_decelerated","stable_concave","decrease_accelerated")){
          max_max<-max_shape
          interval_mid<-(max(dataset$Year)-min(dataset$Year))/2+min(dataset$Year)
          interval_quart1<-interval_mid-(max(dataset$Year)-min(dataset$Year))*0.25
          interval_quart3<-interval_mid+(max(dataset$Year)-min(dataset$Year))*0.25
          if(in.interval.lo(interval_quart1, lo=(mean(simulated$p_1, na.rm=T)-sd(simulated$p_1, na.rm=T)), hi=(mean(simulated$p_1, na.rm=T)+sd(simulated$p_1, na.rm=T)))){
            max_max<-"decrease_accelerated"
          }
          if(in.interval.lo(interval_quart3, lo=(mean(simulated$p_1, na.rm=T)-sd(simulated$p_1, na.rm=T)), hi=(mean(simulated$p_1, na.rm=T)+sd(simulated$p_1, na.rm=T)))){
            max_max<-"increase_decelerated"
          }
          max_shape <- c("increase_decelerated","stable_concave","decrease_accelerated")
        }
        if(max_shape %in% c("decrease_decelerated","stable_convex","increase_accelerated")){
          max_max<-max_shape
          interval_mid<-(max(dataset$Year)-min(dataset$Year))/2+min(dataset$Year)
          interval_quart1<-interval_mid-(max(dataset$Year)-min(dataset$Year))*0.25
          interval_quart3<-interval_mid+(max(dataset$Year)-min(dataset$Year))*0.25
          if(in.interval.lo(interval_quart1, lo=(mean(simulated$p_1, na.rm=T)-sd(simulated$p_1, na.rm=T)), hi=(mean(simulated$p_1, na.rm=T)+sd(simulated$p_1, na.rm=T)))){
            max_max<-"increase_accelerated"
          }
          if(in.interval.lo(interval_quart3, lo=(mean(simulated$p_1, na.rm=T)-sd(simulated$p_1, na.rm=T)), hi=(mean(simulated$p_1, na.rm=T)+sd(simulated$p_1, na.rm=T)))){
            max_max<-"decrease_decelerated"
          }
          max_shape <- c("decrease_decelerated","stable_convex","increase_accelerated")
        }
      }else{
        max_shape <- c("increase_constant","decrease_constant","stable_constant")[which.max(c(length(grep("increase",simulated$shape_class)),
                                                                                              length(grep("decrease",simulated$shape_class)),
                                                                                              length(grep("stable",simulated$shape_class))))]
      }
    }
    if(length(levels(simulated$shape_class)) == 1){max_shape <- levels(simulated$shape_class)}
    if(is.null(max_max)){max_max<-max_shape}
    
    alpha2 <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 1]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]))
    sd_alpha2 <- sd(as.numeric(simulated[simulated$shape_class %in% max_shape, 1]))
    alpha1 <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 2]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]))
    sd_alpha1 <- sd(as.numeric(simulated[simulated$shape_class %in% max_shape, 2]))
    inter <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 5]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]))
    strd <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 3]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]))
    p_1 <- weighted.mean(as.numeric(simulated[simulated$shape_class==max_max, 6]),as.numeric(simulated[simulated$shape_class==max_max, 12]), na.rm=T)
    sd_p_1 <- sd(as.numeric(simulated[simulated$shape_class==max_max, 6]), na.rm=T)
    if( !is.na(p_1) && findInterval(p_1, c(min(dataset$Year), max(dataset$Year))) != 1){p_1 <- sd_p_1 <- as.numeric(NA)}
    p_2 <- weighted.mean(as.numeric(simulated[simulated$shape_class==max_max, 7]),as.numeric(simulated[simulated$shape_class==max_max, 12]), na.rm=T)
    sd_p_2 <- sd(as.numeric(simulated[simulated$shape_class==max_max, 7]), na.rm=T)
    if( !is.na(p_2) && findInterval(p_2, c(min(dataset$Year), max(dataset$Year))) != 1){p_2 <- sd_p_2 <- as.numeric(NA)}
    p_3 <- weighted.mean(as.numeric(simulated[simulated$shape_class==max_max, 8]),as.numeric(simulated[simulated$shape_class==max_max, 12]), na.rm=T)
    sd_p_3 <- sd(as.numeric(simulated[simulated$shape_class==max_max, 8]), na.rm=T)
    if( !is.na(p_3) && findInterval(p_3, c(min(dataset$Year), max(dataset$Year))) != 1){p_3 <- sd_p_3 <- as.numeric(NA)}
    second_order_pvalue <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 9]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]), na.rm=T)
    first_order_pvalue <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 10]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]), na.rm=T)
    slope <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 11]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]), na.rm=T)
    slope_sd <- sd(as.numeric(simulated[simulated$shape_class %in% max_shape, 11]), na.rm=T)
    
    if(!is.null(max_max)){max_shape<-max_max}
    
  }else{alpha2 <- alpha1 <- sd_alpha1 <- inter <- strd <- p_1 <- sd_p_1 <- p_2 <- sd_p_2 <- p_3 <- sd_p_3 <- max_shape <- second_order_pvalue <- first_order_pvalue <- slope <- slope_sd <- NA}
  
  ref_year <- simulated[1,13]
  
  return(data.frame(alpha2, alpha1,sd_alpha1, inter, strd, p_1, sd_p_1, p_2, sd_p_2, p_3, sd_p_3,
                    second_order_pvalue, first_order_pvalue, slope, slope_sd, ref_year,max_shape = as.factor(max_shape)))
}


random_index_bis <- function(dataset, ref_year, dataset2, ref_value="Index"){
  dataset<-droplevels(dataset)

  a <- data.frame(matrix(NA, ncol=1, nrow=nrow(dataset)))
  base_year <- as.numeric(names(table(dataset$Year))[1])
  last_year <- as.numeric(names(table(dataset$Year))[length(table(dataset$Year))])
  
  if(anyNA(dataset$Index)){
    years <- dataset$Year[is.na(dataset$Index)]
    for(z in 1:length(years)){
      if(years[z]==base_year){
        k <- z
        while(k<=length(years) & is.na(dataset$Index[which(dataset$Year == (base_year+k))])){
          k <- k+1
        }
        year_t <- mean(dataset2$Index[which(dataset2$Year == years[z])])
        year_t_1 <- mean(dataset2$Index[which(dataset2$Year == years[k])])
        rat <- year_t/year_t_1
        dataset$Index[which(dataset$Year == years[z])] <- rat*dataset$Index[which(dataset$Year == (base_year+k))]
        dataset$Index_SE[which(dataset$Year == years[z])] <- 0
      }else{
        year_t <- mean(dataset2$Index[which(dataset2$Year == (years[z]-1))])
        year_t_1 <- mean(dataset2$Index[which(dataset2$Year == years[z])])
        rat <- year_t_1/year_t
        dataset$Index[which(dataset$Year == years[z])] <- rat*dataset$Index[which(dataset$Year == (years[z]-1))]
        dataset$Index_SE[which(dataset$Year == years[z])] <- 0
      }
    }
  }
  
  correction <- NULL
  
  if(length(dataset$Index[dataset$Year == ref_year])==0){
    ref_year<-base_year
  }
  
  if(dataset$Index[dataset$Year == last_year]>(5*dataset$Index[dataset$Year == ref_year])){
    correction <- dataset$Index[dataset$Year == last_year]
    dataset$Index <- dataset$Index/correction*100
    dataset$Index_SE <- dataset$Index_SE/correction*100
    ref_year2 <- last_year
  }
  if(dataset$Index[dataset$Year == base_year]>(5*dataset$Index[dataset$Year == ref_year])){
    correction <- dataset$Index[dataset$Year == base_year]
    dataset$Index <- dataset$Index/correction*100
    dataset$Index_SE <- dataset$Index_SE/correction*100
    ref_year2 <- base_year
  }
  
  dataset$sd <- dataset$Index_SE/dataset$Index
  dataset$log <- log(dataset$Index)
  for(j in 1:nrow(dataset)){
    if(dataset$sd[j]>(dataset$log[j])){dataset$sd[j] <- (dataset$log[j])}
  }
  a[,1] <- rnorm(nrow(dataset), mean = dataset$log, sd = dataset$sd)
  if(is.null(correction)){
    a_1 <- a[which(dataset$Year == ref_year), 1]
    a[,1] <- a[,1]-a_1+log(100)
  }else{
    a_1 <- a[which(dataset$Year == ref_year2), 1]
    a[,1] <- a[,1]-a_1+log(100)
  }
  
  if(ref_value %in% c("Abundance","Biomass")){
    a<-log(exp(a)/100*dataset$Index[dataset$Year==ref_year])
  }
  res<-data.frame(Year=dataset$Year,Index=a[,1], Abundance=dataset$Abundance)
  
  return(res)
}

mc_trend3 <-  function(dataset, ref_year,ref_value="Index"){
  dataset <-  droplevels(dataset)
  
  dataset$Index<-dataset$Index_imputed
  dataset$Index_SE<-dataset$Index_imputed_SE
  
  if(ref_value=="Abundance"){
    dataset$Index<-dataset$Abundance
    dataset$Index_SE<-dataset$Abundance_SE
  }
  
  if(ref_value=="Biomass"){
    dataset$Index<-dataset$biomass
    dataset$Index_SE<-dataset$biomass_SE
  }

  result <-  ddply(dataset, .(Species), .fun=random_index_bis, ref_year,subset(dataset, !(Species %in% levels(droplevels(as.factor(as.character(dataset$Species))[is.na(dataset$Index)])))),ref_value)
  result2 <- as.data.frame(result %>% group_by(Year) %>% summarize(Index_SE=sd(Index, na.rm=T), Index=mean(Index, na.rm=T)))
  
  if(ref_value %in% c("Abundance","Biomass")){
    result2 <- as.data.frame(result %>% group_by(Year) %>% summarize(Index_SE=sd(Index, na.rm=T), Index=log(sum(exp(Index), na.rm=T))))
  }
  return(result2)
}

msi_fun3 <-  function(dataset, ref_year, niter, ref_value="Index"){
  aaa <-  mc_trend3(dataset, ref_year, ref_value)
  for(i in 2:niter){
    aaa <-  rbind(aaa,mc_trend3(dataset, ref_year, ref_value))
  }
  
  aaa_finish<- as.data.frame(aaa %>% group_by(Year) %>%summarize(Index_SE=sd(Index, na.rm=T), Index=mean(Index, na.rm=T)))
  mean_msi <-  aaa_finish$Index
  sd_msi <-  aaa_finish$Index_SE
  
  aaa2 <-  data.frame(matrix(NA, ncol=length(table(dataset$Year)), nrow=niter))
  
  b <- data.frame(t(rep(NA, 12)))
  attributes(b)$names <- c("second_order_coef", "first_order_coef", "strd_error", "shape_class", "intercept", "p_1", "p_2", "p_3", "second_order_pvalue", "slope_p_value","linear_slope","r.sq")
  
  
  if(ref_value=="Abundance"){
    ab_init<-sum(as.numeric(dataset$Abundance[dataset$Year==ref_year]))
  }
  
  if(ref_value=="Biomass"){
    bio_init<-sum(as.numeric(dataset$biomass[dataset$Year==ref_year]))
  }
  
  for(i in 1:niter){
    aaa2[i,] <- rnorm(length(table(dataset$Year)), mean=mean_msi, sd=sd_msi)
    
    aaa2_mod <- exp(aaa2[i,]-aaa2[i, which(as.numeric(names(table(dataset$Year))) == ref_year)] + log(100))
    
    if(ref_value=="Abundance"){
      aaa2_mod<-aaa2_mod*ab_init/100
    }
    
    if(ref_value=="Biomass"){
      aaa2_mod<-aaa2_mod*bio_init/100
    }
    
    a <- class.trajectory(unlist(aaa2_mod), as.numeric(names(table(dataset$Year))))
    b[i, 1] <- a$second_order_coef
    b[i, 2] <- a$first_order_coef
    b[i, 3] <- a$strd_error
    b[i, 4] <- a$shape_class
    b[i, 5] <- a$intercept
    if(a$second_order_coef!=0){
      if(findInterval(a$p1, c(min(dataset$Year), max(dataset$Year))) == 1){
        b[i, 6] <- a$p1}else{b[i, 6] <- NA}
      if(findInterval(a$p2, c(min(dataset$Year), max(dataset$Year))) == 1){
        b[i, 7] <- a$p2}else{b[i, 7] <- NA}
      if(findInterval(a$p3, c(min(dataset$Year), max(dataset$Year))) == 1){
        b[i, 8] <- a$p3}else{b[i, 8] <- NA}
    }else{
      b[i, 6] <- NA
      b[i, 7] <- NA
      b[i, 8] <- NA
    }
    b[i, 9] <- a$second_order_pvalue
    b[i, 10] <- a$first_order_pvalue
    b[i, 11] <- a$linear_slope
    b[i, 12] <- max(0,a$r.sq)
  }
  b[, 4] <- as.factor(b[, 4])
  
  max_max<-NULL
  
  if(length(levels(b$shape_class))>1){
    test<-multinomial.theo.multcomp(b$shape_class, p = rep(1/length(levels(b$shape_class)),
                                                           length(levels(b$shape_class))), prop=TRUE)
    if(min(test$p.value2[test$observed>test$expected])<0.05){
      max_shape <- row.names(test$p.value)[which(test$observed == max(test$observed[test$observed>test$expected]))]
      if(max_shape %in% c("increase_decelerated","stable_concave","decrease_accelerated")){
        max_max<-max_shape
        interval_mid<-(max(dataset$Year)-min(dataset$Year))/2+min(dataset$Year)
        interval_quart1<-interval_mid-(max(dataset$Year)-min(dataset$Year))*0.25
        interval_quart3<-interval_mid+(max(dataset$Year)-min(dataset$Year))*0.25
        if(in.interval.lo(interval_quart1, lo=(mean(b$p_1, na.rm=T)-sd(b$p_1, na.rm=T)), hi=(mean(b$p_1, na.rm=T)+sd(b$p_1, na.rm=T)))){
          max_max<-"decrease_accelerated"
        }
        if(in.interval.lo(interval_quart3, lo=(mean(b$p_1, na.rm=T)-sd(b$p_1, na.rm=T)), hi=(mean(b$p_1, na.rm=T)+sd(b$p_1, na.rm=T)))){
          max_max<-"increase_decelerated"
        }
        max_shape <- c("increase_decelerated","stable_concave","decrease_accelerated")
      }
      if(max_shape %in% c("decrease_decelerated","stable_convex","increase_accelerated")){
        max_max<-max_shape
        interval_mid<-(max(dataset$Year)-min(dataset$Year))/2+min(dataset$Year)
        interval_quart1<-interval_mid-(max(dataset$Year)-min(dataset$Year))*0.25
        interval_quart3<-interval_mid+(max(dataset$Year)-min(dataset$Year))*0.25
        if(in.interval.lo(interval_quart1, lo=(mean(b$p_1, na.rm=T)-sd(b$p_1, na.rm=T)), hi=(mean(b$p_1, na.rm=T)+sd(b$p_1, na.rm=T)))){
          max_max<-"increase_accelerated"
        }
        if(in.interval.lo(interval_quart3, lo=(mean(b$p_1, na.rm=T)-sd(b$p_1, na.rm=T)), hi=(mean(b$p_1, na.rm=T)+sd(b$p_1, na.rm=T)))){
          max_max<-"decrease_decelerated"
        }
        max_shape <- c("decrease_decelerated","stable_convex","increase_accelerated")
      }
    }else{
      max_shape <- c("increase_constant","decrease_constant","stable_constant")[which.max(c(length(grep("increase",b$shape_class)),
                                                                                            length(grep("decrease",b$shape_class)),
                                                                                            length(grep("stable",b$shape_class))))]
    }
  }
  if(length(levels(b$shape_class)) == 1){max_shape <- levels(b$shape_class)}
  if(is.null(max_max)){max_max<-max_shape}

  alpha2 <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 1]),as.numeric(b[b$shape_class %in% max_shape, 12]))
  sd_alpha2 <- sd(as.numeric(b[b$shape_class %in% max_shape, 1]))
  alpha1 <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 2]),as.numeric(b[b$shape_class %in% max_shape, 12]))
  sd_alpha1 <- sd(as.numeric(b[b$shape_class %in% max_shape, 2]))
  inter <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 5]),as.numeric(b[b$shape_class %in% max_shape, 12]))
  strd <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 3]),as.numeric(b[b$shape_class %in% max_shape, 12]))
  p_1 <- weighted.mean(as.numeric(b[b$shape_class==max_max, 6]),as.numeric(b[b$shape_class==max_max, 12]), na.rm=T)
  sd_p_1 <- sd(as.numeric(b[b$shape_class==max_max, 6]), na.rm=T)
  if( !is.na(p_1) && !in.interval.lo(p_1, lo=min(dataset$Year), hi=max(dataset$Year))){p_1 <- sd_p_1 <- as.numeric(NA)}
  p_2 <- weighted.mean(as.numeric(b[b$shape_class==max_max, 7]),as.numeric(b[b$shape_class==max_max, 12]), na.rm=T)
  sd_p_2 <- sd(as.numeric(b[b$shape_class==max_max, 7]), na.rm=T)
  if( !is.na(p_2) && !in.interval.lo(p_2, lo=min(dataset$Year), hi=max(dataset$Year))){p_2 <- sd_p_2 <- as.numeric(NA)}
  p_3 <- weighted.mean(as.numeric(b[b$shape_class==max_max, 8]),as.numeric(b[b$shape_class==max_max, 12]), na.rm=T)
  sd_p_3 <- sd(as.numeric(b[b$shape_class==max_max, 8]), na.rm=T)
  if( !is.na(p_3) && !in.interval.lo(p_3, lo=min(dataset$Year), hi=max(dataset$Year))){p_3 <- sd_p_3 <- as.numeric(NA)}
  second_order_pvalue <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 9]),as.numeric(b[b$shape_class %in% max_shape, 12]), na.rm=T)
  first_order_pvalue <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 10]),as.numeric(b[b$shape_class %in% max_shape, 12]), na.rm=T)
  slope <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 11]),as.numeric(b[b$shape_class %in% max_shape, 12]), na.rm=T)
  slope_sd <- sd(as.numeric(b[b$shape_class %in% max_shape, 11]), na.rm=T)
  
  mean_msi_final<-apply(aaa2, 2, function(x){exp(mean(x,na.rm = T))})
  sd_msi_final<-apply(aaa2, 2, function(x){sd(x,na.rm = T)})
  sd_msi_final<-sd_msi_final*mean_msi_final
  sd_msi_final<-sd_msi_final/mean_msi_final[which(levels(as.factor(df_pop2$Year))==ref_year)]*100
  mean_msi_final<-mean_msi_final/mean_msi_final[which(levels(as.factor(df_pop2$Year))==ref_year)]*100

  if(ref_value=="Abundance"){
    mean_msi_final<-mean_msi_final*ab_init/100
    sd_msi_final<-sd_msi_final*ab_init/100
  }
  
  if(ref_value=="Biomass"){
    mean_msi_final<-mean_msi_final*bio_init/100
    sd_msi_final<-sd_msi_final*bio_init/100
  }
  
  if(!is.null(max_max)){max_shape<-max_max}
  return(list(msi = data.frame(mean_msi_final, sd_msi_final),
              coef = data.frame(alpha2, alpha1,sd_alpha1, inter, strd, p_1, sd_p_1, p_2, sd_p_2, p_3, sd_p_3,
                                second_order_pvalue, first_order_pvalue, slope, slope_sd, max_shape = as.factor(max_shape))))
}
```

### Overall avifauna dynamic
```{r}
df_pop2<-SI.imputed 
niter <- 1000
df_pop3 <- msi_fun3(df_pop2,ref_year = 1980, niter = niter, ref_value = "Index")
df_pop3b <- msi_fun3(df_pop2,ref_year = 1980, niter = niter, ref_value = "Abundance")

ggplot(df_pop3$msi, aes(x = c(1980:2016), y = mean_msi_final))+
  geom_ribbon(aes(ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final), alpha=0.5,fill = "lightgrey")+
  geom_point() + theme_modern(base_size = 20)+
  labs(x ="Year", y = "Relative abundance")+
  stat_function(fun = function(x){df_pop3$coef$alpha2*x^2+df_pop3$coef$alpha1*x+df_pop3$coef$inter})+
  stat_function(fun = function(x){df_pop3$coef$alpha2*x^2+df_pop3$coef$alpha1*x+df_pop3$coef$inter-1.96*df_pop3$coef$strd}, linetype = "dashed")+
  stat_function(fun = function(x){df_pop3$coef$alpha2*x^2+df_pop3$coef$alpha1*x+df_pop3$coef$inter+1.96*df_pop3$coef$strd}, linetype = "dashed")

ggplot(df_pop3b$msi, aes(x = c(1980:2016), y = mean_msi_final))+
  geom_ribbon(aes(ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final), alpha=0.5, fill = "lightgrey")+
  geom_point() + theme_modern(base_size = 20)+
  labs(x ="Year", y = "Abundance")+
  stat_function(fun = function(x){df_pop3b$coef$alpha2*x^2+df_pop3b$coef$alpha1*x+df_pop3b$coef$inter})+
  stat_function(fun = function(x){df_pop3b$coef$alpha2*x^2+df_pop3b$coef$alpha1*x+df_pop3b$coef$inter-1.96*df_pop3b$coef$strd}, linetype = "dashed")+
  stat_function(fun = function(x){df_pop3b$coef$alpha2*x^2+df_pop3b$coef$alpha1*x+df_pop3b$coef$inter+1.96*df_pop3b$coef$strd}, linetype = "dashed")
```

### Dynamics by habitat
```{r}
# Farmland species
pecbms_hab<-read.csv2("Habitat_class_PECBMS.csv") # availble on https://pecbms.info/
df_agri_pec <-  droplevels(subset(df_pop2, Species %in% pecbms_hab$Species[pecbms_hab$Habitat=="Farmland"]))
msi_agri_pec_ab<- msi_fun3(df_agri_pec,ref_year = 1980, niter = niter, ref_value = "Abundance")

ggplot(msi_agri_pec_ab$msi, aes(x = c(1980:2016), y = mean_msi_final))+
  geom_ribbon(aes(ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final), alpha=0.5, fill = "lightgrey")+
  geom_point() + theme_modern(base_size = 20)+
  labs(x ="Year", y = "Abundance")+
  stat_function(fun = function(x){msi_agri_pec_ab$coef$alpha2*x^2+msi_agri_pec_ab$coef$alpha1*x+msi_agri_pec_ab$coef$inter})+
  stat_function(fun = function(x){msi_agri_pec_ab$coef$alpha2*x^2+msi_agri_pec_ab$coef$alpha1*x+msi_agri_pec_ab$coef$inter-1.96*msi_agri_pec_ab$coef$strd}, linetype = "dashed")+
  stat_function(fun = function(x){msi_agri_pec_ab$coef$alpha2*x^2+msi_agri_pec_ab$coef$alpha1*x+msi_agri_pec_ab$coef$inter+1.96*msi_agri_pec_ab$coef$strd}, linetype = "dashed")

# Forest species
df_forest_pec <-  droplevels(subset(df_pop2, Species %in% pecbms_hab$Species[pecbms_hab$Habitat=="Forest"]))
msi_forest_pec_ab<- msi_fun3(df_forest_pec,ref_year = 1980, niter = niter, ref_value = "Abundance")

ggplot(msi_forest_pec_ab$msi, aes(x = c(1980:2016), y = mean_msi_final))+
  geom_ribbon(aes(ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final), alpha=0.5,fill = "lightgrey")+
  geom_point() + theme_modern(base_size = 20)+
  labs(x ="Year", y = "Abundance")+
  stat_function(fun = function(x){msi_forest_pec_ab$coef$alpha2*x^2+msi_forest_pec_ab$coef$alpha1*x+msi_forest_pec_ab$coef$inter})+
  stat_function(fun = function(x){msi_forest_pec_ab$coef$alpha2*x^2+msi_forest_pec_ab$coef$alpha1*x+msi_forest_pec_ab$coef$inter-1.96*msi_forest_pec_ab$coef$strd}, linetype = "dashed")+
  stat_function(fun = function(x){msi_forest_pec_ab$coef$alpha2*x^2+msi_forest_pec_ab$coef$alpha1*x+msi_forest_pec_ab$coef$inter+1.96*msi_forest_pec_ab$coef$strd}, linetype = "dashed")

# Urban species
# from https://www.eea.europa.eu/data-and-maps/data/linkages-of-species-and-habitat

eunis_hab<-read.csv("species_birds_maes_EU27b.csv")
eunis_hab2<-dcast(eunis_hab[eunis_hab$season=="B",],speciesname~codeeco)
eunis_hab3<-data.frame(Species=levels(as.factor(eunis_hab$speciesname)), is_urban=FALSE)
eunis_hab3$Species<-as.character(eunis_hab3$Species)
eunis_hab3$Species[eunis_hab3$Species=="Carduelis chloris"]<-"Chloris chloris"
eunis_hab3$Species[eunis_hab3$Species=="Parus caeruleus"]<-"Cyanistes caeruleus"
eunis_hab3$Species[eunis_hab3$Species=="Dendrocopos medius"]<-"Dendrocoptes medius"
eunis_hab3$Species[eunis_hab3$Species=="Dendrocopos minor"]<-"Dryobates minor"
eunis_hab3$Species[eunis_hab3$Species=="Miliaria calandra"]<-"Emberiza calandra"
eunis_hab3$Species[eunis_hab3$Species=="Carduelis cannabina"]<-"Linaria cannabina"
eunis_hab3$Species[eunis_hab3$Species=="Parus cristatus"]<-"Lophophanes cristatus"
eunis_hab3$Species[eunis_hab3$Species=="Tetrao tetrix"]<-"Lyrurus tetrix"
eunis_hab3$Species[eunis_hab3$Species=="Parus ater"]<-"Periparus ater"
eunis_hab3$Species[eunis_hab3$Species=="Parus montanus"]<-"Poecile montanus"
eunis_hab3$Species[eunis_hab3$Species=="Parus palustris"]<-"Poecile palustris"
eunis_hab3$Species[eunis_hab3$Species=="Hirundo rupestris"]<-"Ptyonoprogne rupestris"
eunis_hab3$Species[eunis_hab3$Species=="Regulus ignicapillus"]<-"Regulus ignicapilla"
eunis_hab3$Species[eunis_hab3$Species=="Carduelis spinus"]<-"Spinus spinus"
eunis_hab3$Species<-as.factor(eunis_hab3$Species)
eunis_hab3$is_urban[eunis_hab3$Species %in% levels(droplevels(eunis_hab$speciesname[eunis_hab$season=="B" & eunis_hab$codeeco=="urban"]))]<-TRUE


df_build_eunis <- droplevels(subset(df_pop2, Species %in% levels(as.factor(eunis_hab3$Species[eunis_hab3$is_urban==T]))))
msi_build_eunis_ab<- msi_fun3(df_build_eunis,ref_year = 1980, niter = niter, ref_value = "Abundance")

ggplot(msi_build_eunis_ab$msi, aes(x = c(1980:2016), y = mean_msi_final))+
  geom_ribbon(aes(ymin = mean_msi_final-sd_msi_final, ymax = mean_msi_final+sd_msi_final), alpha=0.5, fill = "lightgrey")+
  geom_point() + theme_modern(base_size = 20)+
  labs(x ="Year", y = "Abundance")+
  stat_function(fun = function(x){msi_build_eunis_ab$coef$alpha2*x^2+msi_build_eunis_ab$coef$alpha1*x+msi_build_eunis_ab$coef$inter})+
  stat_function(fun = function(x){msi_build_eunis_ab$coef$alpha2*x^2+msi_build_eunis_ab$coef$alpha1*x+msi_build_eunis_ab$coef$inter-1.96*msi_build_eunis_ab$coef$strd}, linetype = "dashed")+
  stat_function(fun = function(x){msi_build_eunis_ab$coef$alpha2*x^2+msi_build_eunis_ab$coef$alpha1*x+msi_build_eunis_ab$coef$inter+1.96*msi_build_eunis_ab$coef$strd}, linetype = "dashed")


# Hot and cold dwellers
sti<-read.csv("STI_Devictor.csv") # from (Devictor et al., 2012)
sti_france<-droplevels(subset(sti, SPECIES %in% df_pop2$Species)) 

df_temp <- droplevels(subset(df_pop2, Species %in% as.character(sti$SPECIES[sti_france$STI>12.695]))) # 20% hotter dwellers
msi_temp_ab1<- msi_fun3(df_temp,ref_year = 1980, niter = niter, ref_value = "Abundance")
df_temp <- droplevels(subset(df_pop2, Species %in% as.character(sti$SPECIES[sti_france$STI<11.15]))) # 20% colder dwellers
msi_temp_ab2<- msi_fun3(df_temp,ref_year = 1980, niter = niter, ref_value = "Abundance")

msi_temp_ab3<-data.frame(rbind(msi_temp_ab1$msi,msi_temp_ab2$msi),
                         fact=c(rep("Hot dwellers", nrow(msi_temp_ab1$msi)),rep("Cold dwellers", nrow(msi_temp_ab2$msi))),
                         year=rep(c(1980:2016),2))

ggplot(msi_temp_ab3, aes(x = year, y = mean_msi_final))+
  geom_ribbon(aes(ymin = mean_msi_final-sd_msi_final, ymax = mean_msi_final+sd_msi_final, fill=fact), alpha=0.5)+
  scale_fill_grey(start=0.8,end=0.2)+
  geom_point() + theme_modern(base_size = 20)+theme(legend.title = element_blank(), legend.position = c(0.8,0.8))+
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