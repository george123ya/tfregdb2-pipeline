import os

# Carpeta de entrada con archivos PDB originales
input_folder = os.environ.get("PDB_DIR", "pdb_files")
# Carpeta de salida para guardar archivos PDB filtrados
output_folder = os.environ.get("PDB_FILTERED_DIR", "pdb_files_filtered")

# Crear carpeta de salida si no existe
os.makedirs(output_folder, exist_ok=True)

# Códigos de aminoácidos estándar
AA_CODES = {
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY',
    'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
    'THR','TRP','TYR','VAL'
}

# Función para filtrar cadenas proteicas en un archivo PDB
def filter_protein_chains(pdb_in, pdb_out):
    previous_line_kept = False
    has_kept_lines = False  # Variable para verificar si se guarda alguna línea
    
    with open(pdb_in, "r") as fin, open(pdb_out, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                res_name = line[17:20].strip().upper()
                if res_name in AA_CODES:
                    fout.write(line)
                    previous_line_kept = True
                    has_kept_lines = True
                else:
                    previous_line_kept = False
            elif line.startswith("TER"):
                if previous_line_kept:
                    fout.write(line)
                previous_line_kept = False
            elif line.startswith(("MODEL", "ENDMDL", "END")):
                fout.write(line)
                previous_line_kept = False
            else:
                previous_line_kept = False

    # Verificar si el archivo quedó vacío
    if not has_kept_lines:
        os.remove(pdb_out)
        return False  # Archivo vacío
    return True  # Archivo con proteínas filtradas

# Listas para registrar resultados
filtered_files = []
empty_files = []

# Iterar sobre archivos PDB en la carpeta de entrada
for fname in os.listdir(input_folder):
    if not fname.lower().endswith(".pdb"):
        continue  # Saltar archivos que no sean PDB

    pdb_in = os.path.join(input_folder, fname)
    pdb_out = os.path.join(output_folder, fname)

    try:
        success = filter_protein_chains(pdb_in, pdb_out)
        if success:
            filtered_files.append(fname)
            print(f"{fname} filtered → {pdb_out}")
        else:
            empty_files.append(fname)
            print(f"{fname} resulted in empty file. Skipping output.")
    except Exception as e:
        print(f"Error filtering {fname}: {e}")

# Resumen final
print("\n=== Resumen final ===")
print(f"Archivos filtrados correctamente: {len(filtered_files)}")
print(f"Archivos vacíos (omitidos): {len(empty_files)}")
if empty_files:
    print("Lista de archivos vacíos:", empty_files)
