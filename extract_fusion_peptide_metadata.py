import argparse
import csv


def get_infos(flist, fdatabase):
    with open(flist, mode="r") as fname1, open(fdatabase, mode="r") as fname2:
        csv_file1 = csv.reader(fname1)
        csv_file2 = csv.reader(fname2)
        next(csv_file1)
        header = next(csv_file2)
        flist_peptides = {row[0] for row in csv_file1}
        peptide_dict = {}
        for peptide in csv_file2:
            pep = peptide[2].replace("(+57.02)", "").replace("(+15.99)", "")
            if pep in flist_peptides and pep not in peptide_dict:
                peptide_dict[pep] = peptide
    return peptide_dict, header


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fusion-sites", required=True, help="CSV from identify_cis_fusion_sites.py / identify_trans_fusion_sites.py")
    parser.add_argument("--candidates", nargs="+", required=True, help="PEAKS de novo candidate CSVs")
    parser.add_argument("--labels", nargs="+", required=True, help="Column-prefix label per --candidates file, same order")
    parser.add_argument("--out", required=True, help="Output CSV path")
    args = parser.parse_args()

    if len(args.candidates) != len(args.labels):
        parser.error("--candidates and --labels must have the same number of entries")

    results, header = [], None
    for fname in args.candidates:
        data, header = get_infos(args.fusion_sites, fname)
        results.append(data)

    all_peptides = set()
    for d in results:
        all_peptides.update(d.keys())

    n_cols = len(header)
    empty = [""] * n_cols

    out_header = ["Peptide"]
    for col in header:
        out_header += [f"{label}_{col}" for label in args.labels]

    with open(args.out, mode="w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(out_header)
        for pep in sorted(all_peptides):
            rows = [d.get(pep, empty) for d in results]
            row = [pep]
            for i in range(n_cols):
                row += [r[i] for r in rows]
            writer.writerow(row)

    print(f"Done! {len(all_peptides)} unique peptides written to {args.out}")


if __name__ == "__main__":
    main()
