# HER2 Isoform Validation Feasibility Audit

Status: current feasibility audit after correcting the target paper to Guardia et al., Genome Research 2025, PMID 40664477.

## Bottom Line

The current local project can support gene-level ERBB2 and bulk RNA program checks, but it cannot directly validate HER2 isoforms. The local RNA files are GDC STAR augmented gene-count TSVs. They contain gene-level counts/TPM, not transcript-level isoform quantification, raw reads, BAMs, or junction-count outputs.

Therefore, we can say GigaTIME features associate with HER2-low/HER2-zero state hypotheses, but we cannot say GigaTIME detects HER2 isoforms unless we obtain or generate isoform labels from appropriate RNA-seq data.

## Correct Paper Requirement

Guardia et al. used unprocessed TCGA short-read RNA-seq data, quantified transcripts with kallisto against an expanded HER2 transcriptome, computed isoform PSI with SUPPA2, and confirmed splicing events with rMATS junction counts. That is a different data level from our current local GDC STAR gene-count TSV files.

## Local Data Audit

| Item | Count / value |
| --- | --- |
| Strict high-trust slides | 171 |
| Strict high-trust cases | 171 |
| ERBB2 expression table cases | 80 |
| STAR manifest cases | 80 |
| Local STAR gene-count files | 110 |
| Local STAR gene-count cases | 110 |
| High-trust cases with local STAR gene counts | 56 |
| Low/zero high-trust cases with local STAR gene counts | 40 |
| Local BAM files under data/tcga_brca | 0 |
| Local FASTQ files under data/tcga_brca | 0 |
| Local junction files under data/tcga_brca | 0 |
| Local isoform files under data/tcga_brca | 0 |

## Example Local Expression File Schema

| Field | Value |
| --- | --- |
| Example file | data/tcga_brca/expression_files/TCGA-3C-AAAU/253aa5dc-9853-462a-9bcd-c2e44817833b.rna_seq.augmented_star_gene_counts.tsv |
| Gene model | GENCODE v36 |
| Has gene_id | True |
| Has gene_name | True |
| Has tpm_unstranded | True |
| Has transcript_id | False |
| Has junction columns | False |

## Feasibility Table

| Analysis | Status | What it can support | Next action |
| --- | --- | --- | --- |
| Gene-level ERBB2 expression check | available_now | Sanity check that clinical HER2 labels broadly track ERBB2 RNA expression. | Use only as context, not as isoform validation. |
| Bulk RNA immune/stromal program validation | available_now_limited | Indirect validation of immune, epithelial, stromal, endothelial, and proliferation programs. | Continue treating as indirect support only; adjust for confounders. |
| HER2 isoform quantification like Guardia et al. | not_available_from_current_local_files | Nothing directly; gene-level STAR counts cannot estimate ERBB2 isoform proportions. | Requires short-read RNA-seq reads or an existing transcript-level isoform matrix from the paper/authors. |
| rMATS junction confirmation like Guardia et al. | not_available_from_current_local_files | Nothing directly right now. | Would require controlled-access BAM/FASTQ or author-provided junction/count matrices. |
| GigaTIME association with HER2 isoform/state labels | blocked_until_isoform_labels | Once labels exist, can test whether image features associate with isoform/state groups. | Ask Galante/Guardia group for TCGA sample-level HER2 isoform group/PSI table or reproduce their RNA-seq pipeline on appropriate reads. |

## Practical Recommendation

For the paper proposal, the strongest honest next move is not to claim isoform detection from H&E. Instead:

1. Present the GigaTIME result as a HER2-low/HER2-zero tissue-context association with strong TCGA confounding caveats.
2. Use gene-level ERBB2 and RNA programs only as indirect context.
3. Request or reproduce sample-level HER2 isoform labels from the Guardia/Galante workflow if we want to test isoform biology.
4. Once isoform labels exist, test whether GigaTIME slide features, tumor-rich tile features, or embeddings predict isoform/state groups under source-site/slide-size controls.
5. Keep the language: image AI predicts or associates with HER2 isoform/state hypotheses; it does not detect HER2 isoforms from the current data.

## Output Files

- `docs/her2_isoform_validation_feasibility.md`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/her2_isoform_validation_feasibility/her2_isoform_validation_feasibility_summary.json`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/her2_isoform_validation_feasibility/her2_isoform_validation_feasibility_table.csv`
