library(ggplot2)

# 1. Автоматически определяем, где лежит скрипт
# Это позволит запускать скрипт из любой папки
script_path <- getwd()
cat("Current working directory: ", script_path, "\n")

# 2. Настраиваем пути относительно корня проекта
# Если ты запускаешь из папки Bioinformatics_Portfolio, пути будут такими:
input_file  <- "data/nucleotide_counts.csv"
output_plot <- "results/dna_distribution.png"

# Проверяем файл
if (file.exists(input_file)) {
    df <- read.csv(input_file)
    
    dna_plot <- ggplot(df, aes(x = Nucleotide, y = Count, fill = Nucleotide)) +
        geom_bar(stat = "identity", color = "black", width = 0.7) +
        scale_fill_manual(values = c("A"="#FF9999", "C"="#99CCFF", "G"="#99FF99", "T"="#FFCC99")) +
        theme_minimal() +
        labs(
            title = "Nucleotide Distribution Analysis",
            subtitle = "Based on Rosalind DNA task results",
            x = "Nucleotide Type",
            y = "Frequency"
        ) +
        theme(legend.position = "none")

    # Создаем папку results, если её нет
    if (!dir.exists("results")) dir.create("results")
    
    ggsave(output_plot, plot = dna_plot, width = 8, height = 5, dpi = 300)
    cat("[SUCCESS] Graph saved to: ", output_plot, "\n")

} else {
    cat("[ERROR] File '", input_file, "' not found!\n")
    cat("Try running: python scripts/dna.py first.\n")
}