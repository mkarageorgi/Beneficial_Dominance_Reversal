#!/usr/bin/env python

import os
import numpy as np
import numpy.ma as ma
import multiprocess as mp
import pandas as pd
from tqdm import tqdm

THRESH = 0.03


def assign_linked_status(d, r2_thresh):
    r2_thresh_lo = THRESH
    ace_r2s = d.copy()

    ace_r2s["S_linked"] = "drop"
    ace_r2s["R1_linked"] = "drop"
    ace_r2s["R2_linked"] = "drop"
    ace_r2s["R3_linked"] = "drop"

    # mark sites linked exclusively to each ace allele
    ace_r2s.loc[
        (ace_r2s["R1_r2"] > r2_thresh)
        & (ace_r2s["R2_r2"] <= r2_thresh_lo)
        & (ace_r2s["R3_r2"] <= r2_thresh_lo)
        & (ace_r2s["S_r2"] <= r2_thresh_lo),
        "R1_linked",
    ] = "linked"
    ace_r2s.loc[
        (ace_r2s["R2_r2"] > r2_thresh)
        & (ace_r2s["R1_r2"] <= r2_thresh_lo)
        & (ace_r2s["R3_r2"] <= r2_thresh_lo)
        & (ace_r2s["S_r2"] <= r2_thresh_lo),
        "R2_linked",
    ] = "linked"
    ace_r2s.loc[
        (ace_r2s["R3_r2"] > r2_thresh)
        & (ace_r2s["R1_r2"] <= r2_thresh_lo)
        & (ace_r2s["R2_r2"] <= r2_thresh_lo)
        & (ace_r2s["S_r2"] <= r2_thresh_lo),
        "R3_linked",
    ] = "linked"
    ace_r2s.loc[
        (ace_r2s["S_r2"] > r2_thresh)
        & (ace_r2s["R1_r2"] <= r2_thresh_lo)
        & (ace_r2s["R2_r2"] <= r2_thresh_lo)
        & (ace_r2s["R3_r2"] <= r2_thresh_lo),
        "S_linked",
    ] = "linked"

    # mark sites unlinked to all alleles
    ace_r2s.loc[
        (ace_r2s["R1_r2"] <= r2_thresh_lo)
        & (ace_r2s["R2_r2"] <= r2_thresh_lo)
        & (ace_r2s["R3_r2"] <= r2_thresh_lo)
        & (ace_r2s["S_r2"] <= r2_thresh_lo),
        "R1_linked",
    ] = "unlinked"
    ace_r2s.loc[
        (ace_r2s["R2_r2"] <= r2_thresh_lo)
        & (ace_r2s["R1_r2"] <= r2_thresh_lo)
        & (ace_r2s["R3_r2"] <= r2_thresh_lo)
        & (ace_r2s["S_r2"] <= r2_thresh_lo),
        "R2_linked",
    ] = "unlinked"
    ace_r2s.loc[
        (ace_r2s["R3_r2"] <= r2_thresh_lo)
        & (ace_r2s["R1_r2"] <= r2_thresh_lo)
        & (ace_r2s["R2_r2"] <= r2_thresh_lo)
        & (ace_r2s["S_r2"] <= r2_thresh_lo),
        "R3_linked",
    ] = "unlinked"
    ace_r2s.loc[
        (ace_r2s["S_r2"] <= r2_thresh_lo)
        & (ace_r2s["R1_r2"] <= r2_thresh_lo)
        & (ace_r2s["R2_r2"] <= r2_thresh_lo)
        & (ace_r2s["R3_r2"] <= r2_thresh_lo),
        "S_linked",
    ] = "unlinked"

    ace_r2s_melt = ace_r2s.melt(
        id_vars=["chrom", "pos"],
        value_vars=[
            "S_r",
            "R1_r",
            "R2_r",
            "R3_r",
            "S_linked",
            "R1_linked",
            "R2_linked",
            "R3_linked",
        ],
    )

    ace_r2s_melt["ace_allele"] = ace_r2s_melt["variable"].apply(
        lambda x: x.replace("_r", "").replace("_linked", "")
    )
    ace_r2s_melt["quantity"] = ace_r2s_melt["variable"].apply(
        lambda x: "r" if "_r" in x else "ace_linked"
    )

    ace_r2s_long = (
        ace_r2s_melt.drop(columns=["variable"])
        .pivot(index=["chrom", "pos", "ace_allele"], columns="quantity", values="value")
        .reset_index()
    )
    ace_r2s_long.loc[ace_r2s_long["r"] < 0, "flip"] = "yes"
    ace_r2s_long["flip"] = ace_r2s_long["flip"].fillna("no")
    ace_r2s_long["threshold"] = r2_thresh

    return ace_r2s_long


if __name__ == "__main__":
    snptables = [
        f"data/snptables/Orchard2021/{f}"
        for f in os.listdir("data/snptables/Orchard2021/")
        if f.endswith(".snpTable.numeric")
    ]

    # read csv files
    snptables = [pd.read_csv(f) for f in snptables]
    # add the chrom column
    snptables = [df.assign(chrom=df.columns[0]) for df in snptables]
    # rename the first column to pos
    snptables = [df.rename(columns={df.columns[0]: "pos"}) for df in snptables]

    snptable = pd.concat(snptables)

    # replace -1 with na
    snptable = snptable.replace(-1, np.nan)
    # snptable = snptable[snptable["chrom"] != "X"]

    snptable.to_csv("data/processed/snptable.csv", index=False)

    ace_table = pd.read_csv("data/raw/ace_haplotypes.csv")

    ace_table = (
        ace_table.assign(
            l1=ace_table["Haplotype.aa"].apply(lambda x: x.split(".")[0]),
            l2=ace_table["Haplotype.aa"].apply(lambda x: x.split(".")[1]),
        )
        .drop(columns=["Haplotype.aa"])
        .rename(columns={"Inbred_Line_ID": "line_id"})
        .melt(
            id_vars=["line_id"],
            value_vars=["l1", "l2"],
            var_name="copy",
            value_name="allele",
        )
        .assign(present=1)
        .groupby(["line_id", "allele"])
        .agg({"present": "sum"})
        .reset_index()
        .pivot(index="line_id", columns="allele", values="present")
        .loc[:, ["IGFG", "VGFG", "VGFA", "VAYG"]]
        .rename(
            columns={
                "IGFG": "Ace_S",
                "VGFG": "Ace_R1",
                "VGFA": "Ace_R2",
                "VAYG": "Ace_R3",
            }
        )
        .fillna(0)
        .T
        / 2
    )

    ace_table.to_csv("data/processed/ace_snptable.csv")

    # make sure columns are ordered the right way in both tables
    snptable = snptable[["chrom", "pos"] + list(ace_table.columns)]

    # compute all r^2 values for each ace allele
    ace_r2s = snptable[["chrom", "pos"]].copy()

    s_row = ma.masked_invalid(ace_table.loc["Ace_S"].values)
    r1_row = ma.masked_invalid(ace_table.loc["Ace_R1"].values)
    r2_row = ma.masked_invalid(ace_table.loc["Ace_R2"].values)
    r3_row = ma.masked_invalid(ace_table.loc["Ace_R3"].values)

    snp_rows = []

    for _, row in tqdm(
        snptable.iterrows(), total=snptable.shape[0], desc="Preparing data..."
    ):
        snp_row = ma.masked_invalid(row[2:].astype(float).values)
        mask = ~snp_row.mask

        snp_rows.append((snp_row, mask))

    def compute_r2(row):
        snp_row, mask = row

        s = ma.corrcoef(snp_row[mask], s_row[mask])[0, 1]
        r1 = ma.corrcoef(snp_row[mask], r1_row[mask])[0, 1]
        r2 = ma.corrcoef(snp_row[mask], r2_row[mask])[0, 1]
        r3 = ma.corrcoef(snp_row[mask], r3_row[mask])[0, 1]

        return np.array([s, r1, r2, r3])

    print("Computing r^2 values in parallel...")
    with mp.Pool(10) as pool:
        results = pool.map(compute_r2, snp_rows)

    results = np.array(results)

    ace_r2s.loc[:, "S_r"] = results[:, 0]
    ace_r2s.loc[:, "R1_r"] = results[:, 1]
    ace_r2s.loc[:, "R2_r"] = results[:, 2]
    ace_r2s.loc[:, "R3_r"] = results[:, 3]

    # also save r^2
    ace_r2s["S_r2"] = ace_r2s["S_r"] ** 2
    ace_r2s["R1_r2"] = ace_r2s["R1_r"] ** 2
    ace_r2s["R2_r2"] = ace_r2s["R2_r"] ** 2
    ace_r2s["R3_r2"] = ace_r2s["R3_r"] ** 2

    ace_r2s_03 = assign_linked_status(ace_r2s, 0.03)
    ace_r2s_1 = assign_linked_status(ace_r2s, 0.1)
    ace_r2s_2 = assign_linked_status(ace_r2s, 0.2)

    d = pd.concat([ace_r2s_03, ace_r2s_1, ace_r2s_2])

    d.to_parquet("data/processed/ace_r2s.parquet")
