import argparse
import csv


def fusion_products_to_windows(proteins, fusion_products, left, right,
                               min_left, min_right, pad_char="X", allow_overlap=True):
    def clean(s):
        return s.strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")

    def find_all(seq, sub):
        hits, i = [], 0
        step = 1 if allow_overlap else max(1, len(sub))
        while True:
            j = seq.find(sub, i)
            if j == -1:
                return hits
            hits.append(j)
            i = j + step

    def extract(seq, start, end):
        left_pad, right_pad = max(0, -start), max(0, end - len(seq))
        s, e = max(0, start), min(len(seq), end)
        return (pad_char * left_pad) + seq[s:e] + (pad_char * right_pad)

    seqs = [clean(p) for p in proteins]
    out, seen = [], set()

    for pep in fusion_products:
        pep = clean(pep)
        n = len(pep)
        if n < min_left + min_right:
            continue

        for k in range(min_left, n - min_right + 1):
            acc, don = pep[:k], pep[k:]

            for a_idx, a_seq in enumerate(seqs, start=1):
                for acc_start in find_all(a_seq, acc):
                    acc_end = acc_start + len(acc)
                    acc_flank = extract(a_seq, acc_end - left, acc_end)
                    acc_source = extract(a_seq, acc_end - left, acc_end + right)

                    for d_idx, d_seq in enumerate(seqs, start=1):
                        for don_start in find_all(d_seq, don):
                            don_flank = extract(d_seq, don_start, don_start + right)
                            don_source = extract(d_seq, don_start - left, don_start + right)
                            mer = acc_flank + don_flank

                            key = (pep, k, a_idx, d_idx, mer, acc_source, don_source, acc_start, don_start)
                            if key in seen:
                                continue
                            seen.add(key)
                            out.append((pep, k, a_idx, d_idx, mer, acc_source, don_source, acc_start, don_start))

    return out


def read_fusion_products(fname):
    with open(fname, mode="r", newline="") as f:
        rows = list(csv.reader(f))
    if not rows:
        return []
    first = rows[0][0].strip() if rows[0] else ""
    start = 1 if not (first.isalpha() and first.isupper()) else 0
    return [row[0].strip() for row in rows[start:] if row and row[0].strip()]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fusion-products", required=True, help="Peptide list: plain text or CSV (col 1)")
    parser.add_argument("--left", type=int, default=7, help="Residues taken from the acceptor side of the junction")
    parser.add_argument("--right", type=int, default=8, help="Residues taken from the donor side of the junction")
    parser.add_argument("--min-left", type=int, default=3, help="Shortest acceptor fragment to consider")
    parser.add_argument("--min-right", type=int, default=3, help="Shortest donor fragment to consider")
    parser.add_argument("--pad-char", default="X", help="Padding for windows running off a protein end")
    parser.add_argument("--no-overlap", action="store_true", help="Disallow overlapping acceptor/donor matches")
    parser.add_argument("--out", default="epitope_15mers.tsv", help="Output TSV path")
    args = parser.parse_args()

    n = input("How many proteins do you want to add? (all subunits/chains): ")
    proteins = [input(f"Protein sequence {i + 1}: ") for i in range(int(n))]

    fusion_products = read_fusion_products(args.fusion_products)
    rows = fusion_products_to_windows(proteins, fusion_products, args.left, args.right,
                                      args.min_left, args.min_right,
                                      pad_char=args.pad_char, allow_overlap=not args.no_overlap)

    width = args.left + args.right
    with open(args.out, "w", newline="") as f:
        f.write("fusion_product\tsplit_k\tacc_subunit\tdon_subunit\t"
                f"fused_{width}mer\tacc_source_{width}mer\tdon_source_{width}mer\tacc_start\tdon_start\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    print(f"Wrote {len(rows)} {width}-mers to {args.out}")


if __name__ == "__main__":
    main()
