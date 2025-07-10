import pandas as pd

# Cargar el archivo CSV
df = pd.read_csv("DatasetGLP1R.csv", sep=";")

# Verificar que estén las columnas necesarias
columnas_necesarias = {"Molecule ChEMBL ID", "Standard Type", "Standard Value"}
if not columnas_necesarias.issubset(df.columns):
    raise ValueError("Faltan columnas necesarias en el archivo CSV")

# Limpiar espacios y saltos en 'Standard Value'
df["Standard Value"] = df["Standard Value"].astype(str).str.strip()

# Convertir a numérico con coerción (valores no convertibles serán NaN)
df["Standard Value"] = pd.to_numeric(df["Standard Value"], errors='coerce')

# Definir tipos válidos
tipos_validos = ["Ki", "IC50", "EC50"]

# Filtrar por tipo válido y valor menor a 1000 (sin valores nulos)
df_filtrado = df[(df["Standard Type"].isin(tipos_validos)) & 
                 (df["Standard Value"] < 1000) & 
                 (df["Standard Value"].notna())]

# Mostrar resultados filtrados
print("Moléculas filtradas (Standard Type = Ki, IC50, EC50 y Standard Value < 1000):\n")
print(df_filtrado[["Molecule ChEMBL ID", "Standard Type", "Standard Value"]].head(20))
print(f"\nTotal moléculas filtradas: {len(df_filtrado)}")

# Guardar resultados
df_filtrado.to_csv("filtrado.csv", sep=";", index=False)
print("\nResultado guardado como filtrado.csv")

