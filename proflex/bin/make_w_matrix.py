import pandas as pd
import glob
import os
from pathlib import Path

# import skbio

# from Bio import AlignIO
from Bio.PDB.PDBList import PDBList

# from biopandas.pdb import PandasPdb

# from skbio.alignment import TabularMSA


def main():
    """
    必要なpdb_idのデータを全て取得。
    実行完了まで数時間かかる。
    """
    prosite_msa_df = pd.read_csv("aligned_nmr_protein.csv")
    print(prosite_msa_df["pdbID"].nunique())
    os.makedirs("nmr_pdb_dir", exist_ok=True)
    PDBList().download_pdb_files(
        pdb_codes=prosite_msa_df["pdbID"].unique(),
        file_format="pdb",
        pdir="nmr_pdb_dir",
    )
    for file_path in glob.glob("nmr_pdb_dir/*"):
        file_name = Path(file_path)
        file_name.rename(file_name.with_suffix(".pdb"))


if __name__ == "__main__":
    main()
