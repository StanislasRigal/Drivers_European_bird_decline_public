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