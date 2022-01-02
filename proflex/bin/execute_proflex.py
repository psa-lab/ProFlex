import subprocess
import glob


def main():
    for pdb_file_path in glob.glob("nmr_pdb_dir/*"):
        subprocess.run(["proflex", "-h", "-nonf", f"{pdb_file_path}"])


if __name__ == "__main__":
    main()
