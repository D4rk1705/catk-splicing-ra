import argparse
import csv


def _clean_peptide(pep):
    return pep.replace("(+57.02)", "").replace("(+15.99)", "")


def get_peptides(fname, protein_ids, count_cols):
    """Peptides matching any of protein_ids (col 39), detected in every replicate (count_cols)."""
    peptide_list = []
    with open(fname, mode="r") as f:
        csv_file = csv.reader(f)
        next(csv_file)
        for row in csv_file:
            if any(pid in row[39] for pid in protein_ids):
                if all(row[c] and float(row[c]) > 0 for c in count_cols):
                    pep = _clean_peptide(row[0])
                    if pep not in peptide_list:
                        peptide_list.append(pep)
    return peptide_list


def get_cleavage_sites(peptide_list, database):
    """8-residue windows whose N- or C-terminal half matches a detected peptide's end."""
    protein_segments = [p[i:i + 8] for p in database for i in range(len(p) - 7)]
    cleavage_sites = []
    for peptide in peptide_list:
        for segment in protein_segments:
            if peptide[0:4] == segment[4:8] and segment not in cleavage_sites:
                cleavage_sites.append(segment)
            if peptide[-4:] == segment[0:4] and segment not in cleavage_sites:
                cleavage_sites.append(segment)
    return cleavage_sites


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", required=True, help="Combined PEAKS de novo candidate CSV for the substrate")
    parser.add_argument("--protein-id", nargs="+", required=True, help="Protein accession substring(s), column 39")
    parser.add_argument("--tryp-cols", nargs="+", type=int, default=[15, 16, 17], help="Trypsin, no hCatK")
    parser.add_argument("--gluc-cols", nargs="+", type=int, default=[12, 13, 14], help="Glu-C, no hCatK")
    parser.add_argument("--cat-tryp-cols", nargs="+", type=int, default=[21, 22, 23], help="Trypsin + hCatK")
    parser.add_argument("--cat-gluc-cols", nargs="+", type=int, default=[18, 19, 20], help="Glu-C + hCatK")
    parser.add_argument("--out", default="hCatK_cleavage_sites.csv", help="Output CSV path")
    args = parser.parse_args()

    protein_number = input("How many proteins do you want to add?: ")
    proteins_list = [input(f"Protein sequence {i + 1}: ") for i in range(int(protein_number))]

    tryp_peptides = get_peptides(args.db, args.protein_id, args.tryp_cols)
    gluc_peptides = get_peptides(args.db, args.protein_id, args.gluc_cols)
    cat_tryp_peptides = get_peptides(args.db, args.protein_id, args.cat_tryp_cols)
    cat_gluc_peptides = get_peptides(args.db, args.protein_id, args.cat_gluc_cols)

    tryp_sites = get_cleavage_sites(tryp_peptides, proteins_list)
    gluc_sites = get_cleavage_sites(gluc_peptides, proteins_list)
    cat_tryp_sites = get_cleavage_sites(cat_tryp_peptides, proteins_list)
    cat_gluc_sites = get_cleavage_sites(cat_gluc_peptides, proteins_list)

    tryp_only = [s for s in cat_tryp_sites if s not in tryp_sites]
    gluc_only = [s for s in cat_gluc_sites if s not in gluc_sites]
    cleavage_sites = list(dict.fromkeys(tryp_only + gluc_only))

    with open(args.out, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["hcatK cleavage sites"])
        for site in cleavage_sites:
            writer.writerow([site])

    print(f"Wrote {len(cleavage_sites)} cleavage sites to {args.out}")


if __name__ == "__main__":
    main()
