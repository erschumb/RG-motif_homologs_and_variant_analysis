#!/usr/bin/env python3

#### HERE WE WILL CLASSIFY VARIANTS OF THE GNOMAD DATA TO LATER ANALYSE AND VISUALISE THEM
## INPUT: tsv data from gnomad (per group (e.g. pos and neg))
## OUTPUT: dataframe with annotated variants (metrics for RG counts, amino acid metrics etc..) --> can always be expanded what metrics we will use


import json
import pandas as pd
import pyranges as pr
import re
from localcider.sequenceParameters import SequenceParameters as SeqParams
from typing import Literal

VariantType = Literal[
    "silent",
    "missense",
    "nonsense",
    "inframe_insertion",
    "inframe_deletion",
    "complex",
    "frameshift"
]

def compute_variant_region_overlap_full(json_file: str, variant_dfs: dict):
    """
    End-to-end workflow:
        1. Load region JSON into a PyRanges object
        2. Merge and label multiple variant datasets
        3. Convert variants to PyRanges
        4. Compute overlaps
        5. Return matched and unmatched variants

    Parameters
    ----------
    json_file : str
        Path to JSON file with region definitions.
    
    variant_dfs : dict
        Mapping of label -> pandas DataFrame.
        Each DF must contain CHROM and POS.
        Example:
            {
                "pos": df_pos,
                "neg": df_neg,
                "random": df_random
            }

    Returns
    -------
    pr_regions : pr.PyRanges
        PyRanges object of loaded regions.

    pr_overlap : pr.PyRanges
        PyRanges object of overlapping variants.

    df_matched : pd.DataFrame
        Variants that overlap at least one region.

    df_unmatched : pd.DataFrame
        Variants that do NOT overlap any region.
    """

    ###########################################################################
    # Step 1. Load JSON regions into PyRanges
    ###########################################################################

    with open(json_file, "r") as f:
        data = json.load(f)

    rows = []
    for entry in data:
        protein = entry.get("protein")
        region_id = entry.get("region_id")
        prot_seq = entry.get("prot_seq")
        dna = entry.get("dna")
        group = entry.get("group")
        prot_region = entry.get("prot_region", [None, None])

        for iv in entry["intervals"]:
            rows.append({
                "Chromosome": str(iv["chrom"]),
                "Start": int(iv["start"]),
                "End": int(iv["end"]),
                "Strand": iv.get("strand", "."),
                "protein": protein,
                "region_id": region_id,
                "prot_region_start": prot_region[0],
                "prot_region_end": prot_region[1],
                "prot_seq": prot_seq,
                "dna": dna,
                "group": group
            })

    df_regions = pd.DataFrame(rows)
    pr_regions = pr.PyRanges(df_regions)

    ###########################################################################
    # Step 2. Merge variant DFs, add labels
    ###########################################################################
    dfs = []
    for label, df in variant_dfs.items():
        tmp = df.copy()
        tmp["group"] = label
        dfs.append(tmp)

    variants = pd.concat(dfs, ignore_index=True)

    ###########################################################################
    # Step 3. Convert variants to PyRanges intervals
    ###########################################################################
    variants_interval = variants.rename(columns={"CHROM": "Chromosome", "POS": "Start"})
    variants_interval["Start"] = variants_interval["Start"]
    variants_interval["End"]   = variants_interval["Start"] + 1

    pr_variants = pr.PyRanges(variants_interval)
    # pr_overlap_relaxed = pr_variants.join(pr_regions.extend(1))

    # print(len(pr_overlap_relaxed))
    # print("____")

    # print(len(pr_variants))
    # print(len(pr_regions))
    ###########################################################################
    # Step 4. Overlap
    ###########################################################################
    # pr_overlap = pr_variants.join(pr_regions.extend(1))
    pr_overlap = pr_variants.join(pr_regions)
    df_matched = pr_overlap.as_df()

    ###########################################################################
    # Step 5. Identify unmatched variants
    ###########################################################################
    key_cols = ["Chromosome", "Start", "End", "REF", "ALT", "group"]
    # print(len(df_matched))
    df_matched_keys = df_matched[key_cols].drop_duplicates()
    # print(len(df_matched_keys))


    df_unmatched = (
        variants_interval
        .merge(df_matched_keys, on=key_cols, how="left", indicator=True)
        .query('_merge == "left_only"')
        .drop(columns=["_merge"])
    )

    return pr_regions, pr_overlap, df_matched, df_unmatched


# def classify_variant(before_dna: str, after_dna: str,
#                      before_aa: str, after_aa: str) -> VariantType:
#     """
#     Classify a protein/DNA change into mutation types:
#     silent, missense, nonsense, inframe_insertion, 
#     inframe_deletion, complex, frameshift.

#     Assumes coding DNA and correct translation.
#     """
#     if after_dna is None:
#         return None
#     if after_aa is None:
#         return None
#     if before_aa == after_aa:
#         return "silent"

#     # 3. Same AA length → could be silent, missense, nonsense, complex
#     # Compare residue by residue
#     diffs = [(i, a, b) for i, (a, b) in enumerate(zip(before_aa, after_aa), 1) if a != b]

#     # 4. If any change introduces a stop codon
#     if any(new == "*" for _, _, new in diffs):
#         return "nonsense"

#     # 5. Single AA change → missense
#     if len(diffs) == 1 and len(before_aa)==len(after_aa):
#         return "missense"
#     # 1. Frameshift check (DNA length change not multiple of 3)
#     dna_len_change = len(after_dna) - len(before_dna)
#     if dna_len_change % 3 != 0:
#         return "frameshift"

#     # 2. In-frame insertion or deletion (no frameshift)
#     aa_len_change = len(after_aa) - len(before_aa)

#     if aa_len_change > 0:
#         return "inframe_insertion"
#     elif aa_len_change < 0:
#         return "inframe_deletion"
#     # 6. Multiple AA changes but no frameshift → complex substitution
#     return "complex"


def classify_variant(before_dna: str, after_dna: str,
                     before_aa: str, after_aa: str) -> VariantType:

    if after_dna is None or after_aa is None:
        return None

    # --------------------------------------------------
    # 1. DNA-level classification FIRST (critical)
    # --------------------------------------------------
    dna_len_change = len(after_dna) - len(before_dna)

    # Frameshift
    if dna_len_change % 3 != 0:
        return "frameshift"

    # In-frame indels
    if dna_len_change > 0:
        return "inframe_insertion"

    if dna_len_change < 0:
        return "inframe_deletion"

    # --------------------------------------------------
    # 2. Now substitutions (same DNA length)
    # --------------------------------------------------

    # Silent
    if before_aa == after_aa:
        return "silent"

    diffs = [
        (i, a, b)
        for i, (a, b) in enumerate(zip(before_aa, after_aa), 1)
        if a != b
    ]

    # Nonsense
    if any(new == "*" for _, _, new in diffs):
        return "nonsense"

    # Missense
    if len(diffs) == 1 and len(before_aa) == len(after_aa):
        return "missense"

    # Complex substitution
    return "complex"


def get_physchem_metrics(before_aa: str, after_aa: str, category: str):
    

    # ---------- per-residue property changes ----------
    # charge_change = 0
    # polarity_change = 0
    # aromatic_change = 0
    # print(before_aa)
    # print(after_aa)
    # print(category)
    if category == "nonsense" or before_aa is None or after_aa is None or '*' in after_aa or after_aa == "":
        return {
        "NCPR_change": None,
        "FCR_change": None,
        "hydropathy_change": None,
        "kappa_change": None,
        "pos_count_change": None,
        "neg_count_change": None,
        "aromaticity_change": None
        }
    # align by min length (simple assumption; this is intentional)
    # L = min(len(before_aa), len(after_aa))

    NCPR_change = SeqParams(after_aa).get_NCPR() - SeqParams(before_aa).get_NCPR()
    FCR_change = SeqParams(after_aa).get_FCR() - SeqParams(before_aa).get_FCR()
    hydropathy_change = SeqParams(after_aa).get_mean_hydropathy() - SeqParams(before_aa).get_mean_hydropathy()
    kappa_change = SeqParams(after_aa).get_kappa() - SeqParams(before_aa).get_kappa()
    pos_count_change = SeqParams(after_aa).get_countPos() - SeqParams(before_aa).get_countPos()
    neg_count_change = SeqParams(after_aa).get_countNeg() - SeqParams(before_aa).get_countNeg()
    aromaticity_change =    ((SeqParams(after_aa).get_amino_acid_fractions()["Y"] + 
                            SeqParams(after_aa).get_amino_acid_fractions()["F"] + 
                            SeqParams(after_aa).get_amino_acid_fractions()["W"]) - 
                            (SeqParams(before_aa).get_amino_acid_fractions()["Y"] + 
                            SeqParams(before_aa).get_amino_acid_fractions()["F"] + 
                            SeqParams(before_aa).get_amino_acid_fractions()["W"] ))


    return {

        "NCPR_change": NCPR_change,
        "FCR_change": FCR_change,
        "hydropathy_change": hydropathy_change,
        "kappa_change": kappa_change,
        "pos_count_change": pos_count_change,
        "neg_count_change": neg_count_change,
        "aromaticity_change": aromaticity_change
    }


def count_RG_positions(seq):
    """Return all start positions of 'RG' motifs (0-based)."""
    return [m.start() for m in re.finditer("RG", seq)]


def rg_change_from_category(category, before_aa, after_aa,
                            mut_pos_dna, ref_dna, alt_dna):
    """
    category: one of
        'silent', 'missense', 'inframe_indel', 'frameshift',
        'nonsense', 'noncoding'
    before_aa: AA sequence before mutation
    after_aa: AA sequence after mutation (None if not applicable)
    mut_pos_dna: 1-based genomic DNA start position of the REF allele
    ref_dna, alt_dna: provided but only needed for frameshift logic

    Returns:
        {
            'rg_before': [...],
            'rg_after': [... or None],
            'gained': int,
            'lost': int,
            'unchanged': int
        }
    """
    # if category == None:
    #     return {
    #         'category': category,
    #         'rg_before': rg_before,
    #         'rg_after': rg_after,
    #         'gained': len(new),
    #         'lost': len(lost),
    #         'unchanged': len(unchanged)
    #     }


    # Count RG before
    rg_before = count_RG_positions(before_aa)

    # # ==============================================================
    # # CASE 1 — Noncoding variants
    # # ==============================================================
    # if category == "noncoding":
    #     return {
    #         'category': category,
    #         'rg_before': rg_before,
    #         'rg_after': rg_before,
    #         'gained': 0,
    #         'lost': 0,
    #         'unchanged': len(rg_before)
    #     }

    # ==============================================================
    # CASE 2 — Silent variants
    # (AA sequences identical)
    # ==============================================================
    if category in ("silent", None):
        return {
            # 'category': category,
            'rg_before': rg_before,
            'rg_after': rg_before,
            'gained': 0,
            'lost': 0,
            'unchanged': len(rg_before)
        }

    # ==============================================================
    # CASE 3 — Missense / In-frame Indels / Nonsense
    # (Direct AA comparison)
    # ==============================================================
    if category in ("missense", "nonsense",  "inframe_insertion", "inframe_deletion"):
        rg_after = count_RG_positions(after_aa)

        lost = len([pos for pos in rg_before if pos not in rg_after])
        gained = len([pos for pos in rg_after if pos not in rg_before])
        unchanged = len(rg_before) - lost

        return {
            # 'category': category,
            'rg_before': rg_before,
            'rg_after': rg_after,
            'gained': gained,
            'lost': lost,
            'unchanged': unchanged
        }

    # ==============================================================
    # CASE 4 — Frameshift
    # Only RG motifs upstream of the mutation position remain valid.
    # ==============================================================
    if category == "frameshift":
        aa_mut_pos = mut_pos_dna // 3
        # print(mut_pos_dna)
        # print(aa_mut_pos)
        rg_after = count_RG_positions(after_aa)

        # RGs before the mutation are unchanged
        unchanged = [pos for pos in rg_before if pos < aa_mut_pos]

        # RGs in the original sequence that overlap or are after the mutation are lost
        lost = [pos for pos in rg_before if pos >= aa_mut_pos]

        # RGs in the mutated sequence that are at or after the mutation are new
        new = [pos for pos in rg_after if pos >= aa_mut_pos]

        return {
            'rg_before': rg_before,
            'rg_after': rg_after,
            'gained': len(new),
            'lost': len(lost),
            'unchanged': len(unchanged)
        }

    # ==============================================================
    # Unknown category
    # ==============================================================
    # print(category)
    raise ValueError(f"Unknown category: {category}")



