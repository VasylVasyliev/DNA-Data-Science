from fpdf import FPDF
import pandas as pd
import os

class BioReport(FPDF):
    def header(self):
        self.set_font('Arial', 'B', 15)
        self.cell(0, 10, 'Viral Genome Analysis Report', 0, 1, 'C')
        self.ln(5)

    def footer(self):
        self.set_y(-15)
        self.set_font('Arial', 'I', 8)
        self.cell(0, 10, f'Page {self.page_no()}', 0, 0, 'C')

def generate_pdf():
    # 1. Загружаем данные
    csv_path = "results/virus_report.csv"
    img_path = "results/summary_plot.png"
    
    if not os.path.exists(csv_path) or not os.path.exists(img_path):
        print("❌ Error: Run mass_analyzer.py and plot_summary.R first!")
        return

    df = pd.read_csv(csv_path)
    pdf = BioReport()
    pdf.add_page()
    
    # 2. Добавляем таблицу
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(0, 10, 'Summary Table:', 0, 1)
    
    pdf.set_font('Arial', '', 10)
    # Заголовки таблицы
    pdf.cell(30, 10, 'Accession', 1)
    pdf.cell(100, 10, 'Organism', 1)
    pdf.cell(30, 10, 'Proteins', 1)
    pdf.ln()
    
    # Данные таблицы
    for _, row in df.iterrows():
        pdf.cell(30, 10, str(row['Accession']), 1)
        pdf.cell(100, 10, str(row['Organism'])[:50], 1)
        pdf.cell(30, 10, str(row['Proteins_Found']), 1)
        pdf.ln()

    # 3. Добавляем график
    pdf.ln(10)
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(0, 10, 'Visual Complexity Analysis:', 0, 1)
    pdf.image(img_path, x=10, y=None, w=180)

    # 4. Сохраняем
    output_pdf = "results/final_report.pdf"
    pdf.output(output_pdf)
    print(f"✅ PDF Report generated: {output_pdf}")

if __name__ == "__main__":
    generate_pdf()