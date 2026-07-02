# Author: Gergely Zahoranszky-Kohalmi, PhD, Nathan Miller
#
# Email: gergely.zahoranszky-kohalmi@nih.gov, millernar@nih.gov
#
# Organization: National Center for Advancing Translational Sciences
#

import json
import os
import gzip
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchiKey, MolToInchi
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

SCRIPT_DIR = Path(__file__).parent
ASKCOS_DIR_NAME = SCRIPT_DIR / '../data/input/asckos_buyables_raw'
ASKCOS_FNAME_OUT = SCRIPT_DIR / '../data/input/askcos_inv.tsv'

rows = []

for filename in os.listdir(ASKCOS_DIR_NAME):
    if filename.endswith(".json.gz"):
        filepath = os.path.join(ASKCOS_DIR_NAME, filename)
        with gzip.open(filepath, "rt") as f:
            data = json.load(f)

        print(f"Parsed {len(data)} compounds from {filename}")
        error_count = 0

        for compound in data:
            smiles = compound.get("smiles")
            source = compound.get("source", "asckos_buyables")
            properties = compound.get("properties", {})

            url = None
            if properties:
                for prop in properties:
                    try:
                        if "link" in prop:
                            url = prop["link"]
                            break
                    except (TypeError, AttributeError):
                        continue

            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    try:
                        inchikey = MolToInchiKey(mol)
                        inchi = MolToInchi(mol)
                    except Exception:
                        error_count += 1
                        inchikey = None
                        inchi = None

                    rows.append({
                        "source_file": filename,
                        "source": source,
                        "smiles": smiles,
                        "inchikey": inchikey,
                        "inchi": inchi,
                        "URL": url,
                    })

        if error_count > 0:
            print(f"Encountered {error_count} errors while processing {filename}")

df = pd.DataFrame(rows)
df.to_csv(ASKCOS_FNAME_OUT, index=False, sep="\t")
print(f"Wrote {len(df)} compounds to {ASKCOS_FNAME_OUT}")
