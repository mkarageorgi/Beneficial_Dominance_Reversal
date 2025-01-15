#!/usr/bin/env python

import pandas as pd
import polars as pl
import numpy as np
import numpy.ma as ma
from tqdm import tqdm
import multiprocess as mp

if __name__ == "__main__":
    # SNPTABLES and R2S
    # compute linkage with R3 vs rest
    ace_table = pd.read_csv("data/processed/ace_snptable.csv", index_col=0)
    snptable = pd.read_csv("data/processed/snptable.csv")

    snptable = snptable[["chrom", "pos"] + list(ace_table.columns)]
    sweep_r2s = snptable[["chrom", "pos"]].copy()

    s_row = ma.masked_invalid(ace_table.loc["Ace_S"].values) + ma.masked_invalid(
        ace_table.loc["Ace_R1"].values
    )
    r_row = ma.masked_invalid(ace_table.loc["Ace_R2"].values) + ma.masked_invalid(
        ace_table.loc["Ace_R3"].values
    )

    snp_rows = []
    for _, row in tqdm(
        snptable.iterrows(), total=snptable.shape[0], desc="Preparing data..."
    ):
        snp_row = ma.masked_invalid(row[2:].astype(float).values)
        mask = ~snp_row.mask

        snp_rows.append((snp_row, mask))

    def compute_r2(row):
        snp_row, mask = row
        r = ma.corrcoef(snp_row[mask], r_row[mask])[0, 1]
        return r

    print("Computing r^2 values in parallel...")
    with mp.Pool(10) as pool:
        results = pool.map(compute_r2, snp_rows)

    results = np.array(results)

    sweep_r2s.loc[:, "r"] = results
    sweep_r2s.loc[:, "r2"] = sweep_r2s["r"] ** 2
    sweep_r2s.to_csv("data/processed/sweep_r2s.csv")

    snptable = pd.read_csv("data/processed/snptable.csv")
    snptable = snptable[["chrom", "pos"] + list(snptable.columns[1:-1])]

    snp_positions = snptable[["chrom", "pos"]]
    snptable = snptable.drop(columns=["chrom", "pos"]).T

    ace_snptable = pd.read_csv("data/processed/ace_snptable.csv")

    # COUNT CONDITION

    def count_thresh(snpcounts, thresh):
        # first, we want many ref in R and many derived in S
        a = (snpcounts["Ri"] > thresh) & (snpcounts["S"] > thresh)
        # or, we want many derived in R and many ref in S
        b = (snpcounts["R"] > thresh) & (snpcounts["Si"] > thresh)
        return a | b

    snptable = pd.read_csv("data/processed/snptable.csv")
    snptable = snptable[["chrom", "pos"] + list(snptable.columns[1:-1])]

    snp_positions = snptable[["chrom", "pos"]]
    snptable = snptable.drop(columns=["chrom", "pos"]).T

    ace_snptable = pd.read_csv("data/processed/ace_snptable.csv")

    snptable["resistance"] = ace_snptable.iloc[2:, 1:].sum()
    snptable["resistance"] = snptable["resistance"].apply(
        lambda x: "R" if x > 0.5 else "S"
    )

    snpcounts = pd.concat(
        [snp_positions, snptable.groupby("resistance").sum().T], axis=1
    )

    print(ace_snptable.iloc[2:, 1:].sum().value_counts())

    snpcounts["Si"] = 59 - snpcounts["S"]
    snpcounts["Ri"] = 17 - snpcounts["R"]

    snpcounts["S"] = snpcounts["S"] / 59
    snpcounts["Si"] = snpcounts["Si"] / 59
    snpcounts["R"] = snpcounts["R"] / 17
    snpcounts["Ri"] = snpcounts["Ri"] / 17

    # JOIN ALL

    # sites = (
    #     pl.read_parquet("data/sites_main.parquet")
    #     .join(
    #         pl.from_pandas(sweep_r2s[["chrom", "pos", "r", "r2"]]),
    #         on=["chrom", "pos"],
    #         how="left",
    #     )
    #     .with_columns(pl.col("r").gt(0).alias("r2_flip"))
    #     .filter(pl.col("sweep") == "post_trt")
    # )


    # OVERWRITE freq0 for post_trt with freq0 for trt
    sites = pl.read_parquet("data/sites_main.parquet")

    sites = (
        sites.filter(pl.col("sweep") == "post_trt")
        .drop("freq0")
        .join(
            sites.filter(pl.col("sweep") == "trt").select(
                ["chrom", "pos", "treatment", "cage", "freq0"]
            ),
            on=["chrom", "pos", "treatment", "cage"],
            how="left",
        )
        .drop_nulls()
        .join(
            pl.from_pandas(sweep_r2s[["chrom", "pos", "r", "r2"]]),
            on=["chrom", "pos"],
            how="left",
        )
        .with_columns(pl.col("r").gt(0).alias("r2_flip"))
    )

    thresh = 0.60
    snps = snpcounts[count_thresh(snpcounts, thresh)][["chrom", "pos"]]

    positions_union = pl.from_pandas(snps[["chrom", "pos"]])
    sites_linked = (
        positions_union.join(sites, on=["chrom", "pos"], how="left")
        .to_pandas()
        .dropna()
    )
    sites_linked = sites_linked.query("r2 > 0.03")
    # sites_linked['pos'] = sites_linked['pos'] / 1e6

    # select matched snps
    sites_unlinked = sites.filter(pl.col("r2") < 0.01).to_pandas()

    # flip sites based on the sign of r, i.e. the column r2_flip
    # we flip the lm_slope, lm_effect, total_delta columns, and set freq0 to 1-freq 0 when r2_flip is True
    sites_linked.loc[sites_linked["r2_flip"], "lm_slope"] = -sites_linked.loc[
        sites_linked["r2_flip"], "lm_slope"
    ]
    sites_linked.loc[sites_linked["r2_flip"], "lm_effect"] = -sites_linked.loc[
        sites_linked["r2_flip"], "lm_effect"
    ]
    sites_linked.loc[sites_linked["r2_flip"], "total_delta"] = -sites_linked.loc[
        sites_linked["r2_flip"], "total_delta"
    ]
    sites_linked.loc[sites_linked["r2_flip"], "freq0"] = (
        1 - sites_linked.loc[sites_linked["r2_flip"], "freq0"]
    )

    # same with unlinked
    sites_unlinked.loc[sites_unlinked["r2_flip"], "lm_slope"] = -sites_unlinked.loc[
        sites_unlinked["r2_flip"], "lm_slope"
    ]
    sites_unlinked.loc[sites_unlinked["r2_flip"], "lm_effect"] = -sites_unlinked.loc[
        sites_unlinked["r2_flip"], "lm_effect"
    ]
    sites_unlinked.loc[sites_unlinked["r2_flip"], "total_delta"] = -sites_unlinked.loc[
        sites_unlinked["r2_flip"], "total_delta"
    ]
    sites_unlinked.loc[sites_unlinked["r2_flip"], "freq0"] = (
        1 - sites_unlinked.loc[sites_unlinked["r2_flip"], "freq0"]
    )

    linked_initials = (
        pl.from_pandas(sites_linked)
        .sort("treatment", "cage")
        .group_by("chrom", "pos", "treatment")
        .agg([pl.col("freq0")])
    )

    unlinked_initials = (
        pl.from_pandas(sites_unlinked)
        .sort("treatment", "cage")
        .group_by("chrom", "pos", "treatment")
        .agg([pl.col("freq0")])
        # .with_columns(pl.col('pos').truediv(1e6).alias('pos'))
    )

    print(f"Linked: \n{linked_initials.head()}")

    positions_E = []

    for chrom in ["2L", "2R", "3L", "3R", "X"]:
        chrom_ids = []
        l = linked_initials.filter(
            (pl.col("chrom") == chrom) & (pl.col("treatment") == "E")
        ).to_pandas()
        u = (
            unlinked_initials.filter(
                (pl.col("chrom") == chrom) & (pl.col("treatment") == "E")
            )
            .to_pandas()
            .reset_index()
        )
        for _, lr in tqdm(l.iterrows(), total=l.shape[0], desc=f"Chromosome {chrom}"):
            pos, lfreqs = lr["pos"], lr["freq0"]
            uclose = u.query("pos > @pos - 500000 and pos < @pos + 500000")
            close_positions = []
            for _, ur in uclose.iterrows():
                upos, ufreqs = ur["pos"], ur["freq0"]
                if np.allclose(lfreqs, ufreqs, atol=0.05):
                    close_positions.append(ur["pos"])
            if len(close_positions) > 0:
                chrom_ids.append(np.random.choice(close_positions))
        positions_E.append(list(set(chrom_ids)))

    positions_P = []

    for chrom in ["2L", "2R", "3L", "3R", "X"]:
        chrom_ids = []
        l = linked_initials.filter(
            (pl.col("chrom") == chrom) & (pl.col("treatment") == "P")
        ).to_pandas()
        u = (
            unlinked_initials.filter(
                (pl.col("chrom") == chrom) & (pl.col("treatment") == "P")
            )
            .to_pandas()
            .reset_index()
        )
        for _, lr in tqdm(l.iterrows(), total=l.shape[0], desc=f"Chromosome {chrom}"):
            pos, lfreqs = lr["pos"], lr["freq0"]
            uclose = u.query("pos > @pos - 500000 and pos < @pos + 500000")
            close_positions = []
            for _, ur in uclose.iterrows():
                upos, ufreqs = ur["pos"], ur["freq0"]
                if np.allclose(lfreqs, ufreqs, atol=0.05):
                    close_positions.append(ur["pos"])
            if len(close_positions) > 0:
                chrom_ids.append(np.random.choice(close_positions))
        positions_P.append(list(set(chrom_ids)))

    sites_unlinked_sample = []
    for i, chrom in enumerate(["2L", "2R", "3L", "3R", "X"]):
        pos_P = positions_P[i]
        pos_E = positions_E[i]
        sites_unlinked_sample.append(
            sites.to_pandas().query(
                'sweep == "post_trt" and treatment == "P" and chrom == @chrom and pos in @pos_P'
            )
        )
        sites_unlinked_sample.append(
            sites.to_pandas().query(
                'sweep == "post_trt" and treatment == "E" and chrom == @chrom and pos in @pos_E'
            )
        )

    sites_unlinked_sample = pd.concat(sites_unlinked_sample)

    sites_linked_plot = (
        sites_linked[
            [
                "chrom",
                "pos",
                "treatment",
                "r2",
                "significance_level",
                "lm_slope",
                "lm_effect",
            ]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    sites_unlinked_sample_plot = (
        sites_unlinked_sample[
            [
                "chrom",
                "pos",
                "treatment",
                "r2",
                "significance_level",
                "lm_slope",
                "lm_effect",
            ]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    pd.concat(
        [
            sites_linked_plot.assign(link="linked"),
            sites_unlinked_sample_plot.assign(link="unlinked"),
        ]
    ).to_csv("plot_data/sites_post.csv")
