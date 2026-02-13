import pandas as pd

def calculate_tes():
    # Данные из наших предыдущих этапов
    # Affinity: Энергия (чем меньше число, тем сильнее связь, берем модуль)
    # Survival: Вероятность дойти до цели сквозь GSH (из теста конкуренции)
    
    data = {
        "Strategy": ["Gold (Au)", "Platinum (Pt) Wild", "Platinum (Pt) Smart"],
        "Affinity_Score": [0.493, 0.838, 0.838], # Сила удара
        "Survival_Rate":  [0.65, 0.22, 0.88],  # Процент дошедших частиц
        "Stability":      [0.6, 0.9, 0.9]      # Насколько долго держит блок
    }

    df = pd.DataFrame(data)
    
    # Финальная формула эффективности (TES)
    # Мы перемножаем факторы, так как если один равен 0, то и всё лекарство бесполезно
    df['Final_TES'] = (df['Affinity_Score'] * df['Survival_Rate'] * df['Stability']) * 100
    
    # Сортировка по успеху
    df = df.sort_values(by='Final_TES', ascending=False)
    
    print("====================================================")
    print("   ИТОГОВЫЙ РЕЙТИНГ ЭФФЕКТИВНОСТИ ДОКИНГА (2026)    ")
    print("====================================================")
    print(df.to_string(index=False))
    print("====================================================")
    
    return df

results = calculate_tes()