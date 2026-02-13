import pandas as pd

# Коэффициенты связывания с Альбумином (1AO6)
# Низкое значение = частица легко отцепляется и идет к мишени
# Высокое значение = частица "застревает" в крови
albumin_affinity = {
    "Au": 0.35, 
    "Pt": 0.72  
}

def analyze_adme():
    results = []
    for metal, affinity in albumin_affinity.items():
        # Индекс распределения (Distribution Index)
        # Чем выше, тем лучше частица доходит до тканей
        dist_index = 1.0 - (affinity * 0.8) 
        
        results.append({
            "Metal": metal,
            "Albu_Binding": affinity,
            "Distribution_Eff": round(dist_index, 3),
            "Status": "Excellent" if dist_index > 0.6 else "Slow Transit"
        })
    
    df = pd.DataFrame(results)
    print("\n=== ADME АНАЛИЗ: ТРАНСПОРТ ЧЕРЕЗ 1AO6 ===")
    print(df.to_string(index=False))
    return df

analyze_adme()