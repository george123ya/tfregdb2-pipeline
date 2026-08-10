import pandas as pd
import pickle
import os
import glob
import re
import multiprocessing as mp
from tqdm import tqdm
from finches.utils import folded_domain_utils
from finches import CALVADOS_frontend

# ==========================
# WINDOW REGIME FUNCTION
# ==========================

def get_regime_window_size(seq_len):

    def make_odd(n):
        return n if n % 2 == 1 else n - 1

    if seq_len <= 10:
        return make_odd(seq_len)
    elif 11 <= seq_len <= 20:
        return 9
    elif 21 <= seq_len <= 40:
        return 15
    else:
        return 31

# ==========================
# WORKER FUNCTION
# ==========================

def process_row(row):

    results = {}
    errors = []
    no_interactors = []

    ensembl_id = row['EnsemblID']
    interactor_uniprot = row['UniProtID']

    matching_idrs = idrs_df[idrs_df['EnsemblID'] == ensembl_id]

    if len(matching_idrs) == 0:
        no_interactors.append(ensembl_id)
        return results, errors, no_interactors

    # Cada worker inicializa su propio CALVADOS
    cf = CALVADOS_frontend()

    for _, idr_row in matching_idrs.iterrows():

        genename = idr_row.get('TF_name', idr_row.get('Gene', None))
        coords = idr_row['Coordinates']
        idr_name = f"{genename}_{coords}"
        idr_seq = idr_row['Sequence']
        window_size = get_regime_window_size(len(idr_seq))

        pdb_path = os.path.join(
            pdb_folder,
            f"Interactor_{interactor_uniprot}.pdb"
        )

        try:
            FD = folded_domain_utils.FoldedDomain(pdb_path)

            X = FD.calculate_idr_surface_patch_interactions(
                idr_seq,
                cf.IMC_object,
                window_size
            )

            matrix = X[1]

            if idr_name not in results:
                results[idr_name] = {}

            results[idr_name][interactor_uniprot] = matrix

        except Exception as e:
            errors.append({
                'IDR': idr_name,
                'Interactor': interactor_uniprot,
                'Error': str(e)
            })

    return results, errors, no_interactors

# ==========================
# MAIN
# ==========================

if __name__ == "__main__":

    pdb_folder = "pdb_NEW_files_filtered"

    print("Loading data...")
    idrs_df = pd.read_csv('TF_DB.csv')
    pairs_df = pd.read_csv('TF_PDB_NEW.csv')

    num_cores = mp.cpu_count() - 4
    print(f"Using {num_cores} cores")

    chunk_size = 50

    # Detectar batches existentes y reanudar
    existing_batches = glob.glob("IDR_interactions_batch_*.pkl")
    if existing_batches:
        batch_numbers = [
            int(re.search(r'batch_(\d+)', f).group(1))
            for f in existing_batches
        ]
        last_completed_batch = max(batch_numbers)
        already_processed = (last_completed_batch + 1) * chunk_size
        batch_counter = last_completed_batch + 1

        # Intentar recuperar logs previos
        error_log_file = f'intermap_error_log_batch_{last_completed_batch}.txt'
        no_interactor_log_file = f'no_interactor_log_batch_{last_completed_batch}.txt'

        error_log = []
        if os.path.exists(error_log_file):
            with open(error_log_file, 'r') as f:
                for line in f:
                    parts = line.strip().split(',')
                    if len(parts) == 3:
                        error_log.append({'IDR': parts[0], 'Interactor': parts[1], 'Error': parts[2]})

        no_interactor_log = []
        if os.path.exists(no_interactor_log_file):
            with open(no_interactor_log_file, 'r') as f:
                no_interactor_log = [line.strip() for line in f]

        print(f"Resuming from interaction {already_processed}")

    else:
        already_processed = 0
        batch_counter = 0
        error_log = []
        no_interactor_log = []

    processed = already_processed
    total_tasks = len(pairs_df)

    interactions = {}

    with mp.Pool(num_cores) as pool:

        results_iterator = pool.imap_unordered(
            process_row,
            pairs_df.iloc[already_processed:].to_dict("records")
        )

        for res, err, no_int in tqdm(
            results_iterator,
            total=total_tasks - already_processed,
            desc="Processing interactions",
            unit="interaction"
        ):

            processed += 1

            error_log.extend(err)
            no_interactor_log.extend(no_int)

            for idr_name in res:
                if idr_name not in interactions:
                    interactions[idr_name] = {}
                interactions[idr_name].update(res[idr_name])

            # Guardado por chunks
            if processed % chunk_size == 0:
                print(f"\nSaving batch {batch_counter} (processed {processed})")
                with open(f'IDR_interactions_batch_{batch_counter}.pkl', 'wb') as f:
                    pickle.dump(interactions, f, protocol=pickle.HIGHEST_PROTOCOL)

                with open(f'intermap_error_log_batch_{batch_counter}.txt', 'w') as f:
                    for e in error_log:
                        f.write(f"{e['IDR']},{e['Interactor']},{e['Error']}\n")

                with open(f'no_interactor_log_batch_{batch_counter}.txt', 'w') as f:
                    for ni in no_interactor_log:
                        f.write(f"{ni}\n")

                interactions.clear()
                error_log.clear()
                no_interactor_log.clear()

                batch_counter += 1

    # Guardar remanente final
    if interactions:
        print(f"\nSaving final batch {batch_counter}")
        with open(f'IDR_interactions_batch_{batch_counter}.pkl', 'wb') as f:
            pickle.dump(interactions, f, protocol=pickle.HIGHEST_PROTOCOL)

        with open(f'intermap_error_log_batch_{batch_counter}.txt', 'w') as f:
            for e in error_log:
                f.write(f"{e['IDR']},{e['Interactor']},{e['Error']}\n")

        with open(f'no_interactor_log_batch_{batch_counter}.txt', 'w') as f:
            for ni in no_interactor_log:
                f.write(f"{ni}\n")

    print("Finished successfully.")