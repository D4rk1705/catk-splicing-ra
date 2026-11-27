import argparse
import csv


def _read_peptides(fname):
    peptide_list = []
    with open(fname, mode="r") as f:
        csv_file = csv.reader(f)
        next(csv_file)
        for row in csv_file:
            pep = row[2].replace("(+57.02)", "").replace("(+15.99)", "")
            if pep not in peptide_list:
                peptide_list.append(pep)
    return peptide_list


def get_split_peptides(fname, database):
    """Peptides splittable into two non-contiguous fragments (min 4 aa each) across
    the supplied protein sequences."""
    peptide_list = _read_peptides(fname)
    products = []

    for peptide in peptide_list:
        for i in range(4, len(peptide)):
            for protein1 in database:
                for protein2 in database:
                    if (
                        peptide[0:i] in protein1
                        and peptide[i:] in protein2
                        and len(peptide[i:]) >= i
                        and peptide not in protein1
                        and peptide not in protein2
                        and peptide not in products
                        and peptide[: len(peptide) - 3] not in protein1
                        and peptide[3:] not in protein1
                        and peptide[: len(peptide) - 3] not in protein2
                        and peptide[3:] not in protein2
                    ):
                        products.append(peptide)

    print(f"{fname} analyzed")
    return products


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--candidates", nargs="+", required=True, help="PEAKS de novo candidate CSVs, active hCatK co-incubation")
    parser.add_argument("--controls", nargs="+", required=True, help="PEAKS de novo candidate CSVs, E64-inactivated hCatK co-incubation")
    parser.add_argument("--exclude-cis", nargs="+", required=True, help="PEAKS de novo candidate CSVs for each substrate alone")
    parser.add_argument("--label", required=True, help="Substrate-pair label for the output header, e.g. DTIIC+SPIKE")
    parser.add_argument("--out", default="trans_fusion_sites.csv", help="Output CSV path")
    args = parser.parse_args()

    protein_number = input("How many proteins do you want to add? (both substrates' subunits/chains): ")
    proteins_list = [input(f"Protein sequence {i + 1}: ") for i in range(int(protein_number))]

    candidate_hits = [get_split_peptides(f, proteins_list) for f in args.candidates]
    control_hits = [get_split_peptides(f, proteins_list) for f in args.controls]
    cis_hits = [get_split_peptides(f, proteins_list) for f in args.exclude_cis]

    combined = [pep for hits in candidate_hits for pep in hits]
    false_positives = {pep for hits in control_hits for pep in hits}
    cis_explained = {pep for hits in cis_hits for pep in hits}

    fusion_products = {}
    for peptide in combined:
        if peptide not in false_positives and peptide not in cis_explained:
            fusion_products[peptide] = combined.count(peptide)

    with open(args.out, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([f"Trans fusion sites for {args.label}", "Replicates"])
        for product, rep in fusion_products.items():
            writer.writerow([product, rep])

    print(f"Wrote {len(fusion_products)} trans fusion candidates to {args.out}")


if __name__ == "__main__":
    main()
