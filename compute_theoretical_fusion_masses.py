import argparse
import sys
import csv
import pandas as pd

aa_list = {
    'A': [71.0371, {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'S': 0}],
    'R': [156.1011, {'C': 6, 'H': 12, 'N': 4, 'O': 1, 'S': 0}],
    'N': [114.0429, {'C': 4, 'H': 6, 'N': 2, 'O': 2, 'S': 0}],
    'D': [115.0269, {'C': 4, 'H': 5, 'N': 1, 'O': 3, 'S': 0}],
    'C': [103.0092, {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'S': 1}],
    'E': [129.0426, {'C': 5, 'H': 7, 'N': 1, 'O': 3, 'S': 0}],
    'Q': [128.0586, {'C': 5, 'H': 8, 'N': 2, 'O': 2, 'S': 0}],
    'G': [57.0215, {'C': 2, 'H': 3, 'N': 1, 'O': 1, 'S': 0}],
    'H': [137.0589, {'C': 6, 'H': 7, 'N': 3, 'O': 1, 'S': 0}],
    'I': [113.0841, {'C': 6, 'H': 11, 'N': 1, 'O': 1, 'S': 0}],
    'L': [113.0841, {'C': 6, 'H': 11, 'N': 1, 'O': 1, 'S': 0}],
    'K': [128.0950, {'C': 6, 'H': 12, 'N': 2, 'O': 1, 'S': 0}],
    'M': [131.0405, {'C': 5, 'H': 9, 'N': 1, 'O': 1, 'S': 1}],
    'F': [147.0684, {'C': 9, 'H': 9, 'N': 1, 'O': 1, 'S': 0}],
    'P': [97.0528, {'C': 5, 'H': 7, 'N': 1, 'O': 1, 'S': 0}],
    'S': [87.0320, {'C': 3, 'H': 5, 'N': 1, 'O': 2, 'S': 0}],
    'T': [101.0477, {'C': 4, 'H': 7, 'N': 1, 'O': 2, 'S': 0}],
    'W': [186.0793, {'C': 11, 'H': 10, 'N': 2, 'O': 1, 'S': 0}],
    'Y': [163.0633, {'C': 9, 'H': 9, 'N': 1, 'O': 2, 'S': 0}],
    'V': [99.0684, {'C': 5, 'H': 9, 'N': 1, 'O': 1, 'S': 0}],
    'X': [157.0852, {'C': 6, 'H': 11, 'N': 3, 'O': 2, 'S': 0}],
    '@': [794.3058, {'C': 37, 'H': 46, 'N': 8, 'O': 10, 'S': 1}],
    '$': [708.3344, {'C': 33, 'H': 44, 'N': 10, 'O': 8, 'S': 0}],
    '&': [683.3062, {'C': 28, 'H': 45, 'N': 9, 'O': 9, 'S': 1}],
}


def compute_fusion_products(input1, input2, round_number):
    if len(input1) == 0 or len(input2) == 0:
        print('\n At least one of the two peptides is empty, please try again')
        sys.exit()

    try:
        rounds = int(round_number)
    except (TypeError, ValueError):
        print('\n The number of rounds must be a whole number, please try again')
        sys.exit()

    if rounds < 1:
        print('\n The number of rounds must be at least 1, please try again')
        sys.exit()

    mass_input1 = sum(aa_list[a][0] for a in input1) + 18.0105
    mass_input2 = sum(aa_list[a][0] for a in input2) + 18.0105

    inputs = {}
    fusion_products = {}
    i = 0
    inputs[input1] = [mass_input1, {'C': 78, 'H': 106, 'N': 20, 'O': 23, 'S': 1}]
    inputs[input2] = [mass_input2, {'C': 36, 'H': 61, 'N': 11, 'O': 14, 'S': 1}]

    while i != rounds:
        pep = {}
        prefix, suffix = {}, {}

        for p in inputs:
            pep[p] = inputs[p]
            n = len(p)

            while n != 0:
                if p[0:n] in pep.keys() or p[0:n] in list(aa_list.keys())[0:20]:
                    pass
                else:
                    trunc = p[0:n]
                    initial = p[0:n + 1]
                    pep[trunc] = [0, {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}]
                    pep[trunc][0] = pep[initial][0] - aa_list[p[n]][0]
                    for atom in pep[trunc][1]:
                        pep[trunc][1][atom] = pep[initial][1][atom] - aa_list[p[n]][1][atom]

                if p[len(p) - n:len(p)] in pep.keys() or p[len(p) - n:len(p)] in list(aa_list.keys())[0:20]:
                    pass
                else:
                    trunc = p[len(p) - n:len(p)]
                    initial = p[len(p) - n - 1:len(p)]
                    pep[trunc] = [0, {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}]
                    pep[trunc][0] = pep[initial][0] - aa_list[p[len(p) - n - 1]][0]
                    for atom in pep[trunc][1]:
                        pep[trunc][1][atom] = pep[initial][1][atom] - aa_list[p[len(p) - n - 1]][1][atom]

                n -= 1

        for peptide in pep:
            inputs[peptide] = pep[peptide]

            if peptide[0] == '@' and peptide[-1] != '$' and peptide[-1] != '&':
                prefix[peptide] = pep[peptide]
            elif peptide[0] != '@' and (peptide[-1] == '$' or peptide[-1] == '&'):
                suffix[peptide] = pep[peptide]
            elif peptide[0] == '@' and (peptide[-1] == '$' or peptide[-1] == '&'):
                pass
            else:
                prefix[peptide] = pep[peptide]
                suffix[peptide] = pep[peptide]

        for pep_pre in prefix.keys():
            for pep_suf in suffix.keys():
                fusion_prod = pep_pre + pep_suf
                if fusion_prod in fusion_products.keys():
                    continue

                mass = prefix[pep_pre][0] + suffix[pep_suf][0] - 18.0105
                inputs[fusion_prod] = [mass, {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}]

                for atom in inputs[fusion_prod][1]:
                    if atom == 'H':
                        inputs[fusion_prod][1][atom] = prefix[pep_pre][1][atom] + suffix[pep_suf][1][atom] - 2
                    elif atom == 'O':
                        inputs[fusion_prod][1][atom] = prefix[pep_pre][1][atom] + suffix[pep_suf][1][atom] - 1
                    else:
                        inputs[fusion_prod][1][atom] = prefix[pep_pre][1][atom] + suffix[pep_suf][1][atom]

                fusion_products[fusion_prod] = [0, 0, {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}]
                fusion_products[fusion_prod][0] = mass
                fusion_products[fusion_prod][1] = i + 1

                for atom in fusion_products[fusion_prod][2]:
                    if atom == 'H':
                        fusion_products[fusion_prod][2][atom] = prefix[pep_pre][1][atom] + suffix[pep_suf][1][atom] - 2
                    elif atom == 'O':
                        fusion_products[fusion_prod][2][atom] = prefix[pep_pre][1][atom] + suffix[pep_suf][1][atom] - 1
                    else:
                        fusion_products[fusion_prod][2][atom] = prefix[pep_pre][1][atom] + suffix[pep_suf][1][atom]

        i += 1

    return fusion_products


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="fusion_products.csv", help="Output CSV path (an .xlsx copy is also written)")
    args = parser.parse_args()

    input1 = input('Enter your first sequence: ')
    input2 = input('Enter your second sequence: ')
    round_number = input('Enter the number of rounds (1 round = cleavage followed by fusion): ')

    fusion_products = compute_fusion_products(input1, input2, round_number)

    with open(args.out, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['Product', 'MW', 'Formula', 'Round'])
        for fusion_prod, (mass, rnd, formula) in fusion_products.items():
            f = formula
            formula_str = f"C{f.get('C')}H{f.get('H')}N{f.get('N')}O{f.get('O')}S{f.get('S')}"
            writer.writerow([fusion_prod, str(round(mass, 4)), formula_str, str(rnd)])

    xlsx_out = args.out.rsplit('.', 1)[0] + '.xlsx'
    pd.read_csv(args.out).to_excel(xlsx_out, index=None, header=True)
    print(f"Wrote {len(fusion_products)} theoretical fusion products to {args.out} and {xlsx_out}")


if __name__ == "__main__":
    main()
