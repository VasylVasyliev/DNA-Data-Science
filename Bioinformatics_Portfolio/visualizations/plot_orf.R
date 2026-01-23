# Load necessary library
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)

# Set path to the data
# Assuming you run this from the 'visualizations' folder
file_path <- "../data/protein_lengths.csv"

# Check if data exists
if (!file.exists(file_path)) {
  stop("Error: protein_lengths.csv not found! Please run orf.py first.")
}

# Read the protein lengths
data <- read.csv(file_path)

# Create the visualization
plot <- ggplot(data, aes(x = length)) +
  # Create histogram with the nice teal color you liked
  geom_histogram(binwidth = 10, fill = "#69b3a2", color = "white", alpha = 0.9) +
  
  # Add the RED DASHED LINE at 100 amino acids
  geom_vline(xintercept = 100, linetype = "dashed", color = "red", linewidth = 0.8) +
  
  # Add a text label near the line
  annotate("text", x = 110, y = 35, label = "Biological Threshold (100 aa)", 
           color = "red", hjust = 0, fontface = "italic") +
  
  # Scientific styling
  theme_minimal() +
  labs(
    title = "Protein Length Distribution in PhiX174 Genome",
    subtitle = paste("Total potential ORFs identified:", nrow(data)),
    x = "Protein Length (Amino Acids)",
    y = "Frequency (Count)",
    caption = "Data Source: DNA-Data-Science Analysis"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 12)
  )

# Save the plot
ggsave("orf_distribution.png", plot, width = 10, height = 6, dpi = 300)

print("Success! The plot 'orf_distribution.png' has been updated with the red line.")