import pandas as pd

# Cargar el archivo CSV
df = pd.read_csv("DatasetGLP1R.csv")

# Verificar que estén las columnas necesarias
columnas_necesarias = {"mol_id", "type", "value"}
if not columnas_necesarias.issubset(df.columns):
    raise ValueError("Faltan columnas necesarias en el archivo CSV")

# Filtrar por tipo y valor menor a 1000
tipos_validos = ["Ki", "IC50", "EC50"]
df_filtrado = df[(df["type"].isin(tipos_validos)) & (df["value"] < 1000)]

# Mostrar resultados
print("Moléculas filtradas (type = Ki, IC50, EC50 y value < 1000):\n")
print(df_filtrado)


