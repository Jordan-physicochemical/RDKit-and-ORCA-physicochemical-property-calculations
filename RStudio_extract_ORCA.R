# ORCA output file parser extracts quantum chemistry properties from multiple ORCA .out file. 
# Date: 2025-10-16
# Author: Jordan Campbell

library(tidyverse)
library(stringr)
setwd("C:/")
getwd()


# Function to extract dipole moment (SCF only for consistency)
extract_dipole_moment <- function(text) {
  # Look for SCF dipole moment section specifically
  dipole_start <- grep("^DIPOLE MOMENT$", text)
  
  if (length(dipole_start) > 0) {
    # Search through each dipole moment section
    for (start in dipole_start) {
      # Check if this section is SCF method (within next 10 lines)
      search_range <- start:(min(start + 15, length(text)))
      
      is_scf <- FALSE
      for (i in search_range) {
        if (grepl("Method.*:.*SCF", text[i])) {
          is_scf <- TRUE
          break
        }
      }
      
      # If SCF section found, extract magnitude
      if (is_scf) {
        # Search for "Magnitude (Debye)" within this section
        mag_search <- start:(min(start + 25, length(text)))
        for (j in mag_search) {
          if (grepl("Magnitude.*\\(Debye\\)", text[j])) {
            nums <- str_extract_all(text[j], "\\d+\\.\\d+")[[1]]
            if (length(nums) > 0) {
              return(as.numeric(nums[1]))
            }
          }
        }
      }
    }
  }
  
  return(NA)
}

# Function to extract a single value after a pattern
extract_value <- function(text, pattern, offset = 0, which_occurrence = "last") {
  line_idx <- grep(pattern, text, fixed = TRUE)
  if (length(line_idx) == 0) return(NA)
  
  if (which_occurrence == "last") {
    line_idx <- line_idx[length(line_idx)] + offset
  } else if (which_occurrence == "first") {
    line_idx <- line_idx[1] + offset
  }
  
  if (line_idx > length(text) || line_idx < 1) return(NA)
  
  # Extract numeric value from line
  nums <- str_extract_all(text[line_idx], "-?\\d+\\.\\d+|-?\\d+")[[1]]
  if (length(nums) == 0) return(NA)
  return(as.numeric(nums[1]))  # Changed from nums[length(nums)] to nums[1]
}

# Function to extract HOMO and LUMO energies
extract_homo_lumo <- function(text) {
  homo <- NA
  lumo <- NA
  
  # Look for orbital energies section
  orb_start <- grep("ORBITAL ENERGIES", text)
  if (length(orb_start) == 0) return(list(homo = NA, lumo = NA))
  
  orb_start <- orb_start[length(orb_start)]
  
  # Parse the orbital energies section
  # Format: NO   OCC          E(Eh)            E(eV)
  search_range <- (orb_start + 3):(min(orb_start + 300, length(text)))
  
  last_occupied_energy <- NA
  first_virtual_energy <- NA
  found_transition <- FALSE
  
  for (i in search_range) {
    line <- text[i]
    
    # Stop at empty line or next section
    if (grepl("^\\s*$", line) || grepl("^[A-Z]", line)) {
      if (found_transition) break
    }
    
    # Parse orbital line: orbital_number  occupancy  energy_Eh  energy_eV
    matches <- str_match(line, "^\\s*(\\d+)\\s+(\\d\\.\\d+)\\s+(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)")
    
    if (!is.na(matches[1])) {
      occupancy <- as.numeric(matches[3])
      energy_eV <- as.numeric(matches[5])  # Extract eV (4th column)
      
      # Occupied orbital (2.0000 for closed shell, 1.0000 for open shell)
      if (occupancy > 0.5) {
        last_occupied_energy <- energy_eV
      } 
      # Virtual orbital (0.0000)
      else if (occupancy < 0.5 && is.na(first_virtual_energy)) {
        first_virtual_energy <- energy_eV
        found_transition <- TRUE
      }
    }
  }
  
  homo <- last_occupied_energy
  lumo <- first_virtual_energy
  
  return(list(homo = homo, lumo = lumo))
}

# Function to extract thermochemistry data
extract_thermochem <- function(text) {
  thermo <- list(
    zpve = NA,
    total_enthalpy = NA,
    total_entropy = NA,
    gibbs_free_energy = NA,
    heat_capacity = NA
  )
  
  # Find thermochemistry section
  thermo_patterns <- c("THERMOCHEMISTRY", "INNER ENERGY", "THERMAL")
  thermo_start <- NULL
  
  for (pattern in thermo_patterns) {
    matches <- grep(pattern, text)
    if (length(matches) > 0) {
      thermo_start <- matches[length(matches)]
      break
    }
  }
  
  if (is.null(thermo_start)) return(thermo)
  
  search_end <- min(thermo_start + 150, length(text))
  
  for (i in thermo_start:search_end) {
    line <- text[i]
    
    # Zero point energy patterns
    if (grepl("Zero point energy", line, ignore.case = TRUE) || 
        grepl("Zero point vibrational energy", line, ignore.case = TRUE)) {
      nums <- as.numeric(str_extract_all(line, "-?\\d+\\.\\d+")[[1]])
      if (length(nums) > 0) thermo$zpve <- nums[length(nums)]
    }
    
    # Total Enthalpy
    if (grepl("Total Enthalpy.*Eh", line, ignore.case = TRUE)) {
      nums <- as.numeric(str_extract_all(line, "-?\\d+\\.\\d+")[[1]])
      if (length(nums) > 0) thermo$total_enthalpy <- nums[length(nums)]
    }
    
    # Total Entropy
    if (grepl("Total entropy", line, ignore.case = TRUE) && grepl("cal", line, ignore.case = TRUE)) {
      nums <- as.numeric(str_extract_all(line, "-?\\d+\\.\\d+")[[1]])
      if (length(nums) > 0) thermo$total_entropy <- nums[length(nums)]
    }
    
    # Gibbs free energy
    if (grepl("Final Gibbs free energy", line, ignore.case = TRUE) || 
        grepl("G-E\\(el\\)", line)) {
      nums <- as.numeric(str_extract_all(line, "-?\\d+\\.\\d+")[[1]])
      if (length(nums) > 0) thermo$gibbs_free_energy <- nums[length(nums)]
    }
    
    # Heat capacity (Cv)
    if (grepl("Total heat capacity", line, ignore.case = TRUE) || 
        grepl("Cv.*cal", line, ignore.case = TRUE)) {
      nums <- as.numeric(str_extract_all(line, "-?\\d+\\.\\d+")[[1]])
      if (length(nums) > 0) thermo$heat_capacity <- nums[length(nums)]
    }
  }
  
  return(thermo)
}

# Function to extract Mulliken charges statistics
extract_mulliken_stats <- function(text) {
  charges <- list(
    min_charge = NA,
    max_charge = NA,
    mean_abs_charge = NA
  )
  
  mulliken_start <- grep("MULLIKEN ATOMIC CHARGES", text)
  if (length(mulliken_start) == 0) return(charges)
  
  mulliken_start <- mulliken_start[length(mulliken_start)]
  
  # Read charges (typically start 2 lines after header)
  charge_values <- c()
  for (i in (mulliken_start + 2):(min(mulliken_start + 1000, length(text)))) {
    if (grepl("^\\s*$", text[i]) || grepl("Sum of atomic charges", text[i])) break
    
    # Extract charge value (typically last number on line)
    nums <- str_extract_all(text[i], "-?\\d+\\.\\d+")[[1]]
    if (length(nums) > 0) {
      charge_values <- c(charge_values, as.numeric(nums[length(nums)]))
    }
  }
  
  if (length(charge_values) > 0) {
    charges$min_charge <- min(charge_values)
    charges$max_charge <- max(charge_values)
    charges$mean_abs_charge <- mean(abs(charge_values))
  }
  
  return(charges)
}

# Function to extract dispersion correction
extract_dispersion <- function(text) {
  disp_lines <- grep("Dispersion correction", text)
  if (length(disp_lines) == 0) return(NA)
  
  line <- text[disp_lines[length(disp_lines)]]
  nums <- as.numeric(str_extract_all(line, "-?\\d+\\.\\d+")[[1]])
  if (length(nums) > 0) return(nums[length(nums)])
  return(NA)
}

# Function to count vibrational frequencies
count_frequencies <- function(text) {
  # Look for frequency section
  freq_patterns <- c("VIBRATIONAL FREQUENCIES", "NORMAL MODES", "FREQUENCIES")
  freq_start <- NULL
  
  for (pattern in freq_patterns) {
    matches <- grep(pattern, text)
    if (length(matches) > 0) {
      freq_start <- matches[length(matches)]
      break
    }
  }
  
  if (is.null(freq_start)) return(0)
  
  count <- 0
  search_range <- (freq_start + 1):(min(freq_start + 2000, length(text)))
  
  for (i in search_range) {
    line <- text[i]
    
    # Stop at next major section
    if (grepl("NORMAL MODES", line) || grepl("THERMOCHEMISTRY", line)) {
      if (i > freq_start + 10) break
    }
    
    # Count lines with frequency format (number followed by cm^-1 value)
    if (grepl("^\\s*\\d+:", line) && grepl("cm\\*\\*-1", line)) {
      count <- count + 1
    }
  }
  
  return(count)
}

# Function to extract correlation energy
extract_correlation <- function(text) {
  # Try multiple patterns for correlation energy
  patterns <- c(
    "E\\(CORR\\)\\(RI\\)",
    "E\\(CORR\\)",
    "RI-MP2 CORRELATION ENERGY",
    "Correlation energy"
  )
  
  for (pattern in patterns) {
    lines <- grep(pattern, text)
    if (length(lines) > 0) {
      line <- text[lines[length(lines)]]
      nums <- as.numeric(str_extract_all(line, "-?\\d+\\.\\d+")[[1]])
      if (length(nums) > 0) return(nums[length(nums)])
    }
  }
  
  return(NA)
}

# Main function to parse a single ORCA output file
parse_orca_output <- function(filepath) {
  # Read file
  text <- tryCatch({
    readLines(filepath, warn = FALSE)
  }, error = function(e) {
    return(NULL)
  })
  
  # If reading failed, return NULL
  if (is.null(text)) return(NULL)
  
  # Extract filename
  filename <- basename(filepath)
  fungicide_name <- str_replace(filename, "\\.out$", "")
  
  # Determine method based on filename
  method <- ifelse(grepl("_EnCalc\\.out$", filename), "Method 2", "Method 1")
  
  # Initialize results WITH UNITS IN COLUMN NAMES
  results <- list(
    filename = filename,
    fungicide_name = fungicide_name,
    method = method
  )
  
  # Extract energetic properties
  results$`total_energy_Eh` <- extract_value(text, "FINAL SINGLE POINT ENERGY")
  # Extract SCF energy (first value is Eh, second is eV)
  scf_line <- grep("Total Energy       :", text, fixed = TRUE)
  if (length(scf_line) > 0) {
    nums <- str_extract_all(text[scf_line[length(scf_line)]], "-?\\d+\\.\\d+")[[1]]
    results$`scf_energy_Eh` <- if(length(nums) >= 1) as.numeric(nums[1]) else NA
  } else {
    results$`scf_energy_Eh` <- NA
  }
  results$`dispersion_correction_Eh` <- extract_dispersion(text)
  results$`correlation_energy_Eh` <- extract_correlation(text)
  
  # Extract dipole moment
  results$`dipole_moment_Debye` <- extract_dipole_moment(text)
  
  # Extract HOMO/LUMO
  homo_lumo <- extract_homo_lumo(text)
  results$`homo_energy_eV` <- homo_lumo$homo
  results$`lumo_energy_eV` <- homo_lumo$lumo
  if (!is.na(results$`homo_energy_eV`) && !is.na(results$`lumo_energy_eV`)) {
    results$`homo_lumo_gap_eV` <- results$`lumo_energy_eV` - results$`homo_energy_eV`
  } else {
    results$`homo_lumo_gap_eV` <- NA
  }
  
  # Extract thermochemistry
  thermo <- extract_thermochem(text)
  results$zpve <- thermo$zpve
  results$total_enthalpy <- thermo$total_enthalpy
  results$total_entropy <- thermo$total_entropy
  results$gibbs_free_energy <- thermo$gibbs_free_energy
  results$heat_capacity <- thermo$heat_capacity
  
  # Extract Mulliken charge statistics
  mulliken <- extract_mulliken_stats(text)
  results$`min_mulliken_charge_e` <- mulliken$min_charge
  results$`max_mulliken_charge_e` <- mulliken$max_charge
  results$`mean_abs_mulliken_charge_e` <- mulliken$mean_abs_charge
  
  # Extract geometry optimization info
  results$optimization_converged <- any(grepl("THE OPTIMIZATION HAS CONVERGED", text))
  results$num_optimization_cycles <- length(grep("GEOMETRY OPTIMIZATION CYCLE", text))
  
  # Count vibrational frequencies
  results$num_vibrational_modes <- count_frequencies(text)
  
  # SCF convergence
  scf_lines <- grep("SCF CONVERGED AFTER", text)
  if (length(scf_lines) > 0) {
    scf_text <- text[scf_lines[length(scf_lines)]]
    scf_cycles <- as.numeric(str_extract(scf_text, "\\d+"))
    results$scf_cycles <- ifelse(length(scf_cycles) > 0, scf_cycles, NA)
  } else {
    results$scf_cycles <- NA
  }
  
  return(as.data.frame(results, stringsAsFactors = FALSE))
}

# Main execution
main <- function() {
  # Set the directory containing ORCA output files
  orca_dir <- "C:/20251014 217 Outputs with removed double ups"
  
  # Check if directory exists
  if (!dir.exists(orca_dir)) {
    stop("Directory not found: ", orca_dir, "\nCurrent working directory: ", getwd())
  }
  
  # Get all .out files
  out_files <- list.files(orca_dir, pattern = "\\.out$", full.names = TRUE)
  
  cat("Found", length(out_files), "ORCA output files\n")
  
  if (length(out_files) == 0) {
    stop("No .out files found in the specified directory!")
  }
  
  # Parse all files
  cat("Parsing files...\n")
  results_list <- list()
  failed_files <- c()
  
  for (i in seq_along(out_files)) {
    if (i %% 10 == 0) cat("Processed", i, "files...\n")
    
    result <- tryCatch({
      parse_orca_output(out_files[i])
    }, error = function(e) {
      failed_files <<- c(failed_files, basename(out_files[i]))
      cat("Error processing file", i, ":", basename(out_files[i]), "\n")
      cat("Error message:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(result)) {
      results_list[[i]] <- result
    } else {
      failed_files <- c(failed_files, basename(out_files[i]))
    }
  }
  
  # Report failed files
  if (length(failed_files) > 0) {
    cat("\nWarning: Failed to process", length(failed_files), "files:\n")
    print(failed_files)
  }
  
  # Combine all results
  cat("Combining results...\n")
  all_results <- bind_rows(results_list)
  
  # Save to CSV
  output_file <- "C:/Users/289840H/Documents/fungicide_properties_compiled_v5.csv"
  write_csv(all_results, output_file)
  
  cat("\nExtraction complete!\n")
  cat("Results saved to:", output_file, "\n")
  cat("Total fungicides processed:", nrow(all_results), "\n")
  cat("Method 1 (direct):", sum(all_results$method == "Method 1"), "\n")
  cat("Method 2 (two-step):", sum(all_results$method == "Method 2"), "\n")
  
  # Print column names
  cat("\nExtracted properties:\n")
  print(names(all_results))
  
  return(all_results)
}

# Run the main function
results <- main()


# Optional: View summary statistics
cat("\n=== Summary Statistics ===\n")
summary(results[, sapply(results, is.numeric)])

# Optional: Check for missing data
cat("\n=== Missing Data Summary ===\n")
missing_summary <- data.frame(
  Property = names(results),
  Missing_Count = sapply(results, function(x) sum(is.na(x))),
  Percentage = sapply(results, function(x) round(100 * sum(is.na(x)) / length(x), 2))
)
print(missing_summary)
