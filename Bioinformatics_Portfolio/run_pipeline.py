import subprocess
import os

def run():
    print("🚀 Starting Bioinformatics Pipeline...")
    
    # 1. Запуск анализатора (Python)
    print("\n[1/3] Downloading and analyzing sequences...")
    subprocess.run(["python3", "scripts/mass_analyzer.py"], check=True)
    
    # 2. Запуск визуализации (R)
    print("\n[2/3] Generating scientific plots in R...")
    subprocess.run(["Rscript", "scripts/plot_summary.R"], check=True)
    
    # 3. Запуск генератора отчета (Python)
    print("\n[3/3] Packing everything into a final PDF report...")
    subprocess.run(["python3", "scripts/report_generator.py"], check=True)
    
    print("\n✨ PIPELINE COMPLETE!")
    print(f"📁 Find your results in: {os.path.abspath('results/final_report.pdf')}")

if __name__ == "__main__":
    run()