import requests
import csv

# --- Paso 1: Extraer nombres de proteína desde TF_DB.csv ---
protein_names = []
with open("TF_DB.csv", "r", encoding="utf-8") as tf_db:
    reader = csv.DictReader(tf_db)
    for row in reader:
        protein_name = row["TF_name"].strip().upper()  # columna con el nombre de proteína
        protein_names.append(protein_name)

# --- Paso 2: Configurar la API de STRING ---
string_api_url = "https://version-12-0.string-db.org/api"
output_format = "tsv-no-header"
method = "interaction_partners"
request_url = "/".join([string_api_url, output_format, method])

params = {
    "identifiers": "%0d".join(protein_names),
    "species": 9606,  # humano
    "caller_identity": "www.awesome_app.org"
}

# --- Paso 3: Solicitud POST a la API ---
response = requests.post(request_url, data=params)

# --- Paso 4: Procesar la respuesta y guardar CSV ---
with open("string_TF_interactions_name_ENSP.csv", "w", newline='', encoding="utf-8") as outfile:
    writer = csv.writer(outfile)
    writer.writerow([
        "Query_Name", "Query_ENSP", "Query_STRINGID",
        "Partner_Name", "Partner_ENSP", "Partner_STRINGID",
        "Raw_Combined_Score"
    ])
    
    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")
        if len(l) >= 6:
            # STRING IDs
            query_stringid = l[0]
            partner_stringid = l[1]

            # Extraer ENSP de los STRING IDs (formato: 9606.ENSPxxxx)
            query_ensp = query_stringid.split(".")[1]
            partner_ensp = partner_stringid.split(".")[1]

            # Nombres de proteína
            query_name = l[2].strip().upper()
            partner_name = l[3].strip().upper()

            # Puntaje combinado
            combined_score = l[5]

            writer.writerow([
                query_name, query_ensp, query_stringid,
                partner_name, partner_ensp, partner_stringid,
                combined_score
            ])