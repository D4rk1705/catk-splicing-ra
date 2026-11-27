# catk-splicing-ra

Analysis code accompanying *"Cathepsin K Mediates the Formation of Potential
Rheumatoid Arthritis-Relevant Cis- and Trans-Spliced Peptides Compatible With
HLA-DR4 Presentation."*

Scripts run on PEAKS Studio exports. Steps 2 to 4 take de novo candidate exports
(`*.denovoAllCandidates.csv`, 17 columns); step 1 takes the DB-search peptide
export (`db.peptides_*.csv`, 41 columns); steps 6 and 7 take a peptide list, the
output of step 2 or 3. Protein sequences are entered at runtime, so each script
works for any substrate. File paths are set via CLI flags.

## Pipeline

1. `identify_hcatk_cleavage_sites.py` - hCatK cleavage sites in a substrate.
2. `identify_cis_fusion_sites.py` - cis-spliced (intra-protein) candidates.
3. `identify_trans_fusion_sites.py` - trans-spliced (inter-protein) candidates.
4. `extract_fusion_peptide_metadata.py` - full PEAKS metadata for peptides from step 2 or 3.
5. `compute_theoretical_fusion_masses.py` - theoretical fusion-product masses (small peptides/FRET, synthetic standards).
6. `generate_fusion_20mers.py` - 20-mer windows (10 + 10) centred on each fusion product's splice junction, for HLA-binding prediction.
7. `generate_fusion_15mers.py` - the same, with a 15-mer window (7 + 8).

Run any script with `-h` for its full argument list.

## Requirements

```
pip install -r requirements.txt
```

Only `compute_theoretical_fusion_masses.py` needs pandas/openpyxl; the rest
use the standard library.

## Notes

- Modification tags stripped from PEAKS peptides: `(+57.02)` (carbamidomethylation), `(+15.99)` (oxidation/hydroxylation).

- Steps 6 and 7 report every occurrence of a fragment in its parent protein. A short
  acceptor or donor in a repetitive region therefore yields many windows that differ
  only in their flanking residues: the acceptor `GPPG` occurs 40 times in collagen, so
  a product splitting as `GPPG` + donor produces 40 windows.

  The published analysis used fewer, because the genomic position of an ambiguous short
  fragment was assigned from the longer co-detected peptides covering the same site
  rather than from the short fragment alone. If a sample contains both `GPSGGPPG|FRTSF`
  and `GPPG|FRTSF`, the site is taken from `GPSGGPPG`, and the other `GPPG` positions
  are not counted as separate sites.

- Run steps 6 and 7 with the sequences the peptides actually derive from, which for the
  collagen products is bovine type II collagen (UniProt P02459).

  For the published HLA-DR4 predictions the collagen-derived flanks were then substituted
  to the human orthologue (UniProt P02458) after the windows were generated, so that the
  predicted ligands correspond to sequences an HLA-DR4 carrier could present. That
  substitution is a separate post-processing step.