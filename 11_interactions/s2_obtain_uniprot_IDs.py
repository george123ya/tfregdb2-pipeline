
import requests
import csv
import os
from tqdm import tqdm
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep

# --- Función para obtener UniProt ID desde Ensembl ---
def obtain_uniprot_id(ensembl_id):
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/id/{ensembl_id}?"
    headers = {"Content-Type": "application/json"}
    try:
        r = requests.get(server + ext, headers=headers, timeout=10)
        if not r.ok:
            return None
        data = r.json()
        for xref in data:
            if xref.get('db_display_name') == "UniProtKB/Swiss-Prot":
                display_id = xref.get('display_id')
                if display_id and "." in display_id:
                    display_id = display_id.split(".")[0]
                return display_id
        return None
    except requests.RequestException:
        return None

# --- Archivos ---
input_csv = "/home/edgardo/Documentos/string_TF_NEW.csv"
output_csv = "TF_Uniprot_NEW.csv"
cache_file = "ensembl_to_uniprot.pkl"

# --- Leer CSV completo ---
with open(input_csv, newline='', encoding='utf-8') as infile:
    reader = csv.DictReader(infile)
    rows = list(reader)

# --- Extraer ENSPs únicos ---
unique_ensp = sorted({row['Partner_EnsemblID'].strip() for row in rows if row['Partner_EnsemblID'].strip()})
print(f"🔹 Total de Partner_EnsemblID únicos: {len(unique_ensp)}")

# --- Cargar cache si existe ---
if os.path.exists(cache_file):
    with open(cache_file, 'rb') as f:
        ensembl_to_uniprot = pickle.load(f)
else:
    ensembl_to_uniprot = {}

# --- Configuración ---
BATCH_SIZE = 50
MAX_THREADS = 10

# --- Función para procesar un ENSP individual con reintentos ---
def process_ensp(ensp):
    if ensp in ensembl_to_uniprot:
        return ensp, ensembl_to_uniprot[ensp]
    for attempt in range(3):
        result = obtain_uniprot_id(ensp)
        if result:
            return ensp, result
        sleep(1 + attempt)
    return ensp, None

# --- Procesar por batches ---
for i in range(0, len(unique_ensp), BATCH_SIZE):
    batch = unique_ensp[i:i+BATCH_SIZE]
    print(f"\n📦 Procesando batch {i//BATCH_SIZE + 1} ({len(batch)} ENSPs)")

    with ThreadPoolExecutor(max_workers=MAX_THREADS) as executor:
        futures = {executor.submit(process_ensp, ensp): ensp for ensp in batch}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Resolviendo ENSPs", ncols=100):
            ensp, uniprot_id = future.result()
            ensembl_to_uniprot[ensp] = uniprot_id

    # Guardar cache después de cada batch para resumir
    with open(cache_file, 'wb') as f:
        pickle.dump(ensembl_to_uniprot, f, protocol=pickle.HIGHEST_PROTOCOL)

# --- Escribir CSV final con UniProtID distribuido ---
with open(output_csv, "w", newline='', encoding='utf-8') as outfile:
    fieldnames = rows[0].keys() | {"UniProtID"}
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()

    for row in rows:
        ensp = row['Partner_EnsemblID'].strip()
        row['UniProtID'] = ensembl_to_uniprot.get(ensp, "")
        writer.writerow(row)

print(f"\n✅ CSV final guardado en {output_csv}")
print(f"💾 Cache de ENSP -> UniProtID guardado en {cache_file}")