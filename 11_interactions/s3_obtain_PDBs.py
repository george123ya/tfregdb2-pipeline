import os
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep
from tqdm import tqdm
import pickle  # ✅ Para cache resumible

# ============================================================
# ⚙️ CONFIGURACIÓN
# ============================================================
INPUT_CSV = "TF_Uniprot_NEW.csv"   # Tu CSV con UniProtIDs
OUTPUT_CSV = "TF_PDB_NEW.csv"
PDB_FOLDER = "./pdb_NEW_files"
MAX_THREADS = 10
CACHE_FILE = "uniprot_to_pdb_cache.pkl"
RETRY_LIMIT = 3


# ============================================================
# 🔹 Funciones para obtener PDB experimental o AlphaFold
# ============================================================
def get_first_pdb_from_uniprot(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}?fields=xref_pdb"
    for attempt in range(RETRY_LIMIT):
        try:
            r = requests.get(url, timeout=10)
            if r.ok:
                pdb_ids = [xref["id"] for xref in r.json().get("uniProtKBCrossReferences", []) if xref.get("database") == "PDB"]
                return pdb_ids[0] if pdb_ids else None
        except Exception:
            sleep(1 + attempt)
    return None


def get_pdburl_alphafold_api(uniprot_id):
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    for attempt in range(RETRY_LIMIT):
        try:
            r = requests.get(url, timeout=10)
            if r.ok:
                for model in r.json():
                    pdb_url = model.get("pdbUrl")
                    if pdb_url:
                        return pdb_url
        except Exception:
            sleep(1 + attempt)
    return None


def download_pdb_from_url(url, filepath):
    try:
        if os.path.exists(filepath):
            return filepath
        r = requests.get(url, timeout=20)
        if r.ok:
            os.makedirs(os.path.dirname(filepath), exist_ok=True)
            with open(filepath, "w") as f:
                f.write(r.text)
            return filepath
    except Exception:
        pass
    return None


# ============================================================
# 🔹 Procesar UniProtID (experimental o AlphaFold)
# ============================================================
def process_uniprot(uniprot_id):
    if not uniprot_id:
        return (None, None, None)

    # ✅ Revisar si ya existe PDB descargado
    pdb_file = os.path.join(PDB_FOLDER, f"Interactor_{uniprot_id}.pdb")
    if os.path.exists(pdb_file):
        return (uniprot_id, "Exists", "cached")

    # 1️⃣ PDB experimental
    pdb_id = get_first_pdb_from_uniprot(uniprot_id)
    if pdb_id:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        if download_pdb_from_url(url, pdb_file):
            return (uniprot_id, pdb_id, "experimental")

    # 2️⃣ AlphaFold
    pdb_url = get_pdburl_alphafold_api(uniprot_id)
    if pdb_url:
        if download_pdb_from_url(pdb_url, pdb_file):
            return (uniprot_id, f"AlphaFold_{uniprot_id}", "predicted")

    return (uniprot_id, None, None)


# ============================================================
# 🚀 MAIN
# ============================================================
if __name__ == "__main__":
    # --- Leer CSV ---
    df = pd.read_csv(INPUT_CSV)
    uniprot_ids = df["UniProtID"].dropna().astype(str).unique()

    # --- Cargar cache resumible ---
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, "rb") as f:
            cache = pickle.load(f)
    else:
        cache = {}

    os.makedirs(PDB_FOLDER, exist_ok=True)

    # --- Ejecutar descargas con ThreadPool ---
    with ThreadPoolExecutor(max_workers=MAX_THREADS) as executor:
        futures = {executor.submit(process_uniprot, uid): uid for uid in uniprot_ids}

        for future in tqdm(as_completed(futures), total=len(futures), desc="Descargando PDBs", ncols=100):
            uid, pdb_id, method = future.result()
            if uid:
                cache[uid] = (pdb_id, method)
                # Guardar cache cada vez para resumible
                with open(CACHE_FILE, "wb") as f:
                    pickle.dump(cache, f)

    # --- Mapear resultados al DataFrame ---
    df["PDB_ID"] = df["UniProtID"].map(lambda uid: cache.get(uid, (None, None))[0])
    df["Method"] = df["UniProtID"].map(lambda uid: cache.get(uid, (None, None))[1])

    # --- Guardar CSV final ---
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\n✅ Proceso completado. Resultados guardados en {OUTPUT_CSV}")
    print(f"💾 Archivos PDB descargados en {PDB_FOLDER}/")