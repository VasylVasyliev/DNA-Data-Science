# Загружаем библиотеку для графиков
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)

# Читаем наш свежий отчет
data <- read.csv("results/virus_report.csv")

# Строим график
plot <- ggplot(data, aes(x = Length_bp, y = Proteins_Found, label = Accession)) +
  geom_point(aes(size = GC_Content_.), color = "steelblue") +
  geom_text(vjust = -1, size = 3) +
  theme_minimal() +
  labs(title = "Virus Genome Complexity",
       x = "Genome Length (bp)",
       y = "Number of Proteins Found",
       size = "GC Content %")

# Сохраняем график
ggsave("results/summary_plot.png", plot, width = 8, height = 6)
print("✅ Scientific plot saved to results/summary_plot.png")