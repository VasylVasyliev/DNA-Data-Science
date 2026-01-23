# R script to count nucleotides (A, C, G, T)
# Using standard R vector operations

input_path <- "data/ini.txt"

# Reading the sequence
if (file.exists(input_path)) {
    dna <- readLines(input_path, warn = FALSE)
    
    # Split string into individual characters
    nucleotides <- strsplit(dna, "")[[1]]
    
    # Create a frequency table
    counts <- table(nucleotides)
    
    # Print results in the specific Rosalind order: A C G T
    # We use paste to format them with spaces
    cat(paste(counts["A"], counts["C"], counts["G"], counts["T"]), "\n")
} else {
    cat("Error: File not found\n")
}