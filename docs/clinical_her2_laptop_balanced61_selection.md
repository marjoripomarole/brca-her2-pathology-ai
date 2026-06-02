# Clinical HER2 Cohort Selection

This file summarizes the balanced clinical HER2-positive / HER2-low / HER2-zero pilot cohort selected for the next GigaTIME run.

Status after download:

- 183 selected cases/slides total: 61 HER2-positive, 61 HER2-low, 61 HER2-zero.
- All 183 selected diagnostic H&E slide files are present locally.
- The cohort occupies about 31 GB under `data/tcga_brca/slides/`.
- This is a prepared input cohort, not an analyzed result yet. GigaTIME still needs to be run on these 183 slides.

## Counts

| Cohort group | Candidate cases | Selected cases | Already-downloaded slides |
|---|---:|---:|---:|
| HER2-positive | 174 | 61 | 61 |
| HER2-low | 407 | 61 | 61 |
| HER2-zero | 61 | 61 | 61 |

## Selection Priority

- Clinical HER2 group in requested groups.
- Direct clinical label before inferred label.
- Already-downloaded slide before not-yet-downloaded slide.
- Smaller slide file before larger slide file.
- Case submitter ID for deterministic tie-breaking.

## Selected Cases

| Group | Rank | Case | HER2 rule | IHC score | ISH status | ER | PR | ERBB2 TPM | Slide downloaded |
|---|---:|---|---|---|---|---|---|---:|---|
| HER2-low | 1 | TCGA-A7-A26J | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | 70.96 | yes |
| HER2-low | 2 | TCGA-WT-AB44 | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | 65.5 | yes |
| HER2-low | 3 | TCGA-B6-A3ZX | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Negative | Negative | 45.48 | yes |
| HER2-low | 4 | TCGA-A7-A13G | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | 54.14 | yes |
| HER2-low | 5 | TCGA-A7-A26I | IHC score 2+ and ISH negative | 2+ | Negative | Negative | Negative | 52.73 | yes |
| HER2-low | 6 | TCGA-B6-A40B | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | 193.4 | yes |
| HER2-low | 7 | TCGA-B6-A409 | IHC score 1+ with no positive ISH | 1+ | Negative | Negative | Negative | 53.4 | yes |
| HER2-low | 8 | TCGA-A7-A26F | IHC score 2+ and ISH negative | 2+ | Negative | Negative | Negative | 32.49 | yes |
| HER2-low | 9 | TCGA-B6-A40C | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | 140.7 | yes |
| HER2-low | 10 | TCGA-E2-A108 | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | 48.31 | yes |
| HER2-low | 11 | TCGA-A7-A0DB | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 12 | TCGA-AC-A3OD | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 13 | TCGA-B6-A402 | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-low | 14 | TCGA-B6-A400 | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-low | 15 | TCGA-AR-A2LO | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 16 | TCGA-E2-A10F | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 17 | TCGA-AN-A0XW | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 18 | TCGA-S3-A6ZF | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 19 | TCGA-AN-A03X | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 20 | TCGA-GM-A5PX | IHC score 1+ with no positive ISH | 1+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 21 | TCGA-A7-A0DC | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Negative | nan | yes |
| HER2-low | 22 | TCGA-A7-A13D | IHC score 2+ and ISH negative | 2+ | Negative | Negative | Positive | nan | yes |
| HER2-low | 23 | TCGA-E2-A106 | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 24 | TCGA-D8-A27P | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 25 | TCGA-A7-A26E | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 26 | TCGA-E2-A10C | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 27 | TCGA-D8-A73U | IHC score 1+ with no positive ISH | 1+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 28 | TCGA-E2-A1LK | IHC score 2+ and ISH negative | 2+ | Negative | Negative | Negative | nan | yes |
| HER2-low | 29 | TCGA-AC-A8OS | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 30 | TCGA-AC-A3QP | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 31 | TCGA-A2-A0SV | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | 110.2 | yes |
| HER2-low | 32 | TCGA-A2-A0CT | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Negative | 257.9 | yes |
| HER2-low | 33 | TCGA-D8-A27G | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 34 | TCGA-A7-A6VX | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 35 | TCGA-AR-A2LK | IHC score 1+ with no positive ISH | 1+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 36 | TCGA-AO-A0JG | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 37 | TCGA-AO-A0J7 | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 38 | TCGA-D8-A1XR | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 39 | TCGA-E2-A14T | IHC score 1+ with no positive ISH | 1+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 40 | TCGA-LL-A442 | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 41 | TCGA-LL-A5YO | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-low | 42 | TCGA-AR-A24N | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 43 | TCGA-A8-A07F | IHC score 1+ with no positive ISH | 1+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 44 | TCGA-E9-A295 | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 45 | TCGA-AR-A2LH | IHC score 2+ and ISH negative | 2+ | Negative | Negative | Negative | nan | yes |
| HER2-low | 46 | TCGA-E2-A1LA | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 47 | TCGA-LL-A8F5 | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Negative | nan | yes |
| HER2-low | 48 | TCGA-C8-A1HM | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 49 | TCGA-AR-A0U2 | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 50 | TCGA-A7-A13E | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Negative | nan | yes |
| HER2-low | 51 | TCGA-EW-A1PA | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 52 | TCGA-D8-A27K | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-low | 53 | TCGA-A8-A08O | IHC score 1+ with no positive ISH | 1+ | Negative | Positive | Positive | nan | yes |
| HER2-low | 54 | TCGA-A1-A0SJ | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | 157.4 | yes |
| HER2-low | 55 | TCGA-A2-A04T | IHC score 2+ and ISH negative | 2+ | Negative | Negative | Negative | 74.48 | yes |
| HER2-low | 56 | TCGA-5L-AAT0 | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | 156.8 | yes |
| HER2-low | 57 | TCGA-A2-A0ES | IHC score 1+ with no positive ISH | 1+ | Negative | Positive | Positive | 151.7 | yes |
| HER2-low | 58 | TCGA-A2-A04Q | IHC score 2+ and ISH negative | 2+ | Negative | Negative | Negative | 35.86 | yes |
| HER2-low | 59 | TCGA-5T-A9QA | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Negative | 216.2 | yes |
| HER2-low | 60 | TCGA-A2-A0EN | IHC score 2+ and ISH negative | 2+ | Negative | Positive | Positive | 92.41 | yes |
| HER2-low | 61 | TCGA-A2-A0T6 | IHC score 1+ with no positive ISH | 1+ | [Not Evaluated] | Positive | Positive | 252.8 | yes |
| HER2-positive | 1 | TCGA-OL-A5S0 | ISH positive | [Not Available] | Positive | Positive | Negative | 100 | yes |
| HER2-positive | 2 | TCGA-AC-A3QQ | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Positive | 264.2 | yes |
| HER2-positive | 3 | TCGA-A2-A04X | IHC score 3+ | 3+ | Positive | Positive | Positive | 877.1 | yes |
| HER2-positive | 4 | TCGA-OL-A5RZ | ISH positive | [Not Available] | Positive | Positive | Negative | 3533 | yes |
| HER2-positive | 5 | TCGA-AC-A23G | ISH positive | [Not Available] | Positive | Positive | Positive | 118.6 | yes |
| HER2-positive | 6 | TCGA-AC-A23C | ISH positive | [Not Available] | Positive | Positive | Positive | 680.2 | yes |
| HER2-positive | 7 | TCGA-AQ-A04L | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Negative | 100.3 | yes |
| HER2-positive | 8 | TCGA-OL-A5RY | ISH positive | [Not Available] | Positive | Positive | Negative | 1876 | yes |
| HER2-positive | 9 | TCGA-EW-A1PD | ISH positive | 2+ | Positive | Positive | Positive | 344.3 | yes |
| HER2-positive | 10 | TCGA-A8-A08C | ISH positive | 2+ | Positive | Positive | Positive | 157.3 | yes |
| HER2-positive | 11 | TCGA-AQ-A1H2 | ISH positive | 2+ | Positive | Positive | Positive | 191.2 | yes |
| HER2-positive | 12 | TCGA-AC-A3W5 | ISH positive | [Not Available] | Positive | Positive | Positive | nan | yes |
| HER2-positive | 13 | TCGA-A2-A1G1 | ISH positive | 2+ | Positive | Negative | Negative | nan | yes |
| HER2-positive | 14 | TCGA-AC-A23H | ISH positive | [Not Available] | Positive | Positive | Negative | nan | yes |
| HER2-positive | 15 | TCGA-AN-A04C | IHC score 3+ | 3+ | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-positive | 16 | TCGA-A7-A425 | IHC score 3+ | 3+ | Negative | Positive | Positive | nan | yes |
| HER2-positive | 17 | TCGA-A8-A09G | ISH positive | 2+ | Positive | Positive | Negative | nan | yes |
| HER2-positive | 18 | TCGA-AR-A2LJ | IHC score 3+ | 3+ | Positive | Positive | Positive | nan | yes |
| HER2-positive | 19 | TCGA-LL-A5YL | IHC score 3+ | 3+ | Negative | Positive | Negative | nan | yes |
| HER2-positive | 20 | TCGA-A7-A4SF | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Negative | nan | yes |
| HER2-positive | 21 | TCGA-LL-A7T0 | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-positive | 22 | TCGA-C8-A3M8 | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-positive | 23 | TCGA-D8-A1J9 | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Negative | nan | yes |
| HER2-positive | 24 | TCGA-A2-A3XV | ISH positive | 2+ | Positive | Positive | Negative | nan | yes |
| HER2-positive | 25 | TCGA-A8-A06U | ISH positive | 2+ | Positive | Positive | Positive | nan | yes |
| HER2-positive | 26 | TCGA-A8-A08T | ISH positive | 2+ | Positive | Positive | Positive | nan | yes |
| HER2-positive | 27 | TCGA-AR-A24U | IHC score 3+ | 3+ | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-positive | 28 | TCGA-C8-A12P | IHC score 3+ | 3+ | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-positive | 29 | TCGA-EW-A2FR | ISH positive | 2+ | Positive | Negative | Negative | nan | yes |
| HER2-positive | 30 | TCGA-AR-A0TX | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-positive | 31 | TCGA-A8-A0AB | ISH positive | 2+ | Positive | Positive | Positive | nan | yes |
| HER2-positive | 32 | TCGA-A2-A0EY | IHC score 3+ | 3+ | Positive | Positive | Negative | 1681 | yes |
| HER2-positive | 33 | TCGA-A8-A08B | IHC score 3+ | 3+ | Positive | Positive | Negative | nan | yes |
| HER2-positive | 34 | TCGA-A7-A4SC | IHC score 3+ | 3+ | Negative | Positive | Negative | nan | yes |
| HER2-positive | 35 | TCGA-A8-A09N | IHC score 3+ | 3+ | Positive | Positive | Positive | nan | yes |
| HER2-positive | 36 | TCGA-A8-A08H | IHC score 3+ | 3+ | Positive | Positive | Positive | nan | yes |
| HER2-positive | 37 | TCGA-EW-A1OZ | ISH positive | 2+ | Positive | Positive | Negative | nan | yes |
| HER2-positive | 38 | TCGA-A8-A07P | ISH positive | 2+ | Positive | Positive | Positive | nan | yes |
| HER2-positive | 39 | TCGA-A8-A075 | ISH positive | 2+ | Positive | Positive | Positive | nan | yes |
| HER2-positive | 40 | TCGA-A2-A0T1 | IHC score 3+ | 3+ | [Not Evaluated] | Negative | Negative | 1699 | yes |
| HER2-positive | 41 | TCGA-E2-A1LE | IHC score 3+ | 3+ | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-positive | 42 | TCGA-A8-A08G | IHC score 3+ | 3+ | Positive | Positive | Positive | nan | yes |
| HER2-positive | 43 | TCGA-D8-A27N | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-positive | 44 | TCGA-AN-A0AJ | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-positive | 45 | TCGA-AC-A3W6 | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-positive | 46 | TCGA-D8-A27W | ISH positive | 2+ | Positive | Positive | Positive | nan | yes |
| HER2-positive | 47 | TCGA-S3-AA14 | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-positive | 48 | TCGA-A8-A0A7 | IHC score 3+ | 3+ | Positive | Negative | Negative | nan | yes |
| HER2-positive | 49 | TCGA-A1-A0SM | IHC score 3+ | 3+ | Positive | Positive | Negative | 3101 | yes |
| HER2-positive | 50 | TCGA-E2-A14W | ISH positive | [Not Available] | Positive | Positive | Positive | nan | yes |
| HER2-positive | 51 | TCGA-AQ-A0Y5 | ISH positive | 2+ | Positive | Positive | Positive | nan | yes |
| HER2-positive | 52 | TCGA-AO-A0JM | IHC score 3+ | 3+ | Positive | Positive | Positive | nan | yes |
| HER2-positive | 53 | TCGA-A7-A2KD | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-positive | 54 | TCGA-AR-A255 | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-positive | 55 | TCGA-A2-A04U | ISH positive | 1+ | Positive | Negative | Negative | 28.21 | yes |
| HER2-positive | 56 | TCGA-A2-A0D1 | IHC score 3+ | 3+ | Positive | Negative | Negative | 1123 | yes |
| HER2-positive | 57 | TCGA-AC-A2FB | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-positive | 58 | TCGA-A2-A04W | ISH positive | [Not Available] | Positive | Negative | Negative | 1133 | yes |
| HER2-positive | 59 | TCGA-A2-A0SY | ISH positive | [Not Available] | Positive | Positive | Positive | 420 | yes |
| HER2-positive | 60 | TCGA-A2-A0EQ | IHC score 3+ | 3+ | Negative | Negative | Negative | 78.19 | yes |
| HER2-positive | 61 | TCGA-A2-A0CX | IHC score 3+ | 3+ | [Not Evaluated] | Positive | Negative | 1466 | yes |
| HER2-zero | 1 | TCGA-AO-A03R | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | 14.58 | yes |
| HER2-zero | 2 | TCGA-AN-A04D | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | 63.86 | yes |
| HER2-zero | 3 | TCGA-AN-A049 | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | 279.5 | yes |
| HER2-zero | 4 | TCGA-AO-A03V | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | 109.2 | yes |
| HER2-zero | 5 | TCGA-A1-A0SP | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | 51.02 | yes |
| HER2-zero | 6 | TCGA-A8-A0A2 | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | 126.2 | yes |
| HER2-zero | 7 | TCGA-AO-A03T | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | 29.25 | yes |
| HER2-zero | 8 | TCGA-A8-A09C | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | 80.34 | yes |
| HER2-zero | 9 | TCGA-AO-A03P | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | 101.1 | yes |
| HER2-zero | 10 | TCGA-AO-A03U | IHC score 0 with no positive ISH | 0 | Negative | Negative | Negative | 61.1 | yes |
| HER2-zero | 11 | TCGA-A2-A0EU | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | 117.1 | yes |
| HER2-zero | 12 | TCGA-A8-A09V | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | 148.2 | yes |
| HER2-zero | 13 | TCGA-AN-A0AS | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Negative | nan | yes |
| HER2-zero | 14 | TCGA-A2-A0T2 | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | 40.53 | yes |
| HER2-zero | 15 | TCGA-AN-A03Y | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 16 | TCGA-BH-A0DQ | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 17 | TCGA-AO-A12E | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 18 | TCGA-AO-A03N | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 19 | TCGA-AN-A0AL | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-zero | 20 | TCGA-A2-A0CM | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | 42.18 | yes |
| HER2-zero | 21 | TCGA-AO-A12H | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 22 | TCGA-A2-A0YH | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 23 | TCGA-AN-A04A | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 24 | TCGA-BH-A0BM | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Negative | nan | yes |
| HER2-zero | 25 | TCGA-A8-A086 | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | nan | yes |
| HER2-zero | 26 | TCGA-AO-A0JJ | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 27 | TCGA-BH-A0HK | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Negative | nan | yes |
| HER2-zero | 28 | TCGA-BH-A0DK | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 29 | TCGA-BH-A0DH | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 30 | TCGA-AO-A124 | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-zero | 31 | TCGA-AO-A12B | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 32 | TCGA-AO-A128 | IHC score 0 with no positive ISH | 0 | Negative | Negative | Negative | nan | yes |
| HER2-zero | 33 | TCGA-A2-A0EV | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | 226.6 | yes |
| HER2-zero | 34 | TCGA-A2-A0D2 | IHC score 0 with no positive ISH | 0 | Negative | Negative | Negative | 61.59 | yes |
| HER2-zero | 35 | TCGA-AQ-A04J | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-zero | 36 | TCGA-AN-A0FY | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 37 | TCGA-AN-A0AR | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-zero | 38 | TCGA-AO-A0JC | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 39 | TCGA-A2-A0YF | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Negative | nan | yes |
| HER2-zero | 40 | TCGA-A2-A0EW | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | 141.1 | yes |
| HER2-zero | 41 | TCGA-BH-A0AY | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 42 | TCGA-A2-A0D0 | IHC score 0 with no positive ISH | 0 | Negative | Negative | Negative | 38.39 | yes |
| HER2-zero | 43 | TCGA-AO-A125 | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | nan | yes |
| HER2-zero | 44 | TCGA-A2-A0YE | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-zero | 45 | TCGA-AO-A129 | IHC score 0 with no positive ISH | 0 | Negative | Negative | Negative | nan | yes |
| HER2-zero | 46 | TCGA-A8-A09W | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | nan | yes |
| HER2-zero | 47 | TCGA-BH-A0HF | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 48 | TCGA-A8-A09K | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | nan | yes |
| HER2-zero | 49 | TCGA-A2-A0YJ | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Negative | nan | yes |
| HER2-zero | 50 | TCGA-A8-A07L | IHC score 0 with no positive ISH | 0 | Negative | Positive | Positive | nan | yes |
| HER2-zero | 51 | TCGA-AN-A0FF | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 52 | TCGA-BH-A0BC | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 53 | TCGA-A1-A0SK | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | 6.718 | yes |
| HER2-zero | 54 | TCGA-AO-A0JB | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 55 | TCGA-A8-A06N | IHC score 0 with no positive ISH | 0 | Negative | Positive | Negative | nan | yes |
| HER2-zero | 56 | TCGA-A2-A0T0 | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | 54.48 | yes |
| HER2-zero | 57 | TCGA-BH-A0DP | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 58 | TCGA-AO-A0J6 | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-zero | 59 | TCGA-BH-A0BV | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |
| HER2-zero | 60 | TCGA-BH-A0B9 | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Negative | Negative | nan | yes |
| HER2-zero | 61 | TCGA-BH-A0H0 | IHC score 0 with no positive ISH | 0 | [Not Evaluated] | Positive | Positive | nan | yes |

## Local Outputs

- Cases CSV: `data/tcga_brca/clinical_her2_laptop_balanced61_cases.csv`
- Slide table: `data/tcga_brca/clinical_her2_laptop_balanced61_slides_files.csv`
- Slide manifest: `data/tcga_brca/clinical_her2_laptop_balanced61_slide_manifest.tsv`

These CSV/TSV files are under `data/`, so they are local reproducible outputs rather than tracked Git files.
