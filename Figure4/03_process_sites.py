#!/usr/bin/env python

import pandas as pd
import polars as pl
import numpy as np
import numpy.ma as ma
from tqdm import tqdm


if __name__ == "__main__":
    samps_initial = (
        pl.read_csv("data/processed/samps.csv")
        .with_columns(pl.col("freq_idx").sub(1))
        .filter(pl.col("tpt").is_in([2, 6]))
        .with_columns(
            pl.concat_str(
                [
                    pl.col("treatment"),
                    pl.lit("_"),
                    pl.col("cage"),
                    pl.lit("-"),
                    pl.col("tpt"),
                ]
            ).alias("cage_label-tpt")
        )
        .select(["cage_label-tpt", "freq_idx"])
    )

    afmat_indices = [
        ":".join(map(str, x))
        for x in zip(samps_initial["cage_label-tpt"], samps_initial["freq_idx"])
    ]

    afmat = np.load("data/processed/afmat.npy")

    sites = (
        pl.scan_csv(
            "data/processed/sites.csv",
        )
        # convert from 1-indexed to 0-indexed
        .with_columns(pl.col("site_idx").sub(1))
        # add r2 values and tag sites by linked status
        .drop_nulls()
        # add the other index needed to use the afmat
        .with_columns(pl.lit(afmat_indices).alias("afmat_indices"))
        .explode("afmat_indices")
        #   split the afmat_indices into cage_label-tpt and freq_idx on :
        .with_columns(
            pl.col("afmat_indices")
            .str.split_exact(":", 1)
            .struct.rename_fields(["cage_label-tpt", "freq_idx"])
        )
        .unnest("afmat_indices")
        .cast({"freq_idx": pl.Int64})
        #   split on - to get cage_label and tpt
        .with_columns(
            pl.col("cage_label-tpt")
            .str.split_exact("-", 1)
            .struct.rename_fields(["cage_label", "tpt"])
        )
        .unnest("cage_label-tpt")
        .cast({"tpt": pl.Int64})
        # split on _ to get treatment and cage
        .with_columns(
            pl.col("cage_label")
            .str.split_exact("_", 1)
            .struct.rename_fields(["treatment", "cage"])
        )
        .unnest("cage_label")
        .cast({"cage": pl.Int64})
        # .filter([pl.col("tpt") == 2])  # , pl.col('treatment') == 'P'])
        # add the afmat values
        .with_columns(
            pl.struct([pl.col("site_idx"), pl.col("freq_idx")]).alias("afmat_indices")
        )
        .with_columns(
            pl.col("afmat_indices")
            .map_elements(lambda s: afmat[s["site_idx"], s["freq_idx"]])
            .alias("initial_freq")
        )
        .drop(["freq_idx", "afmat_indices"])
    ).collect()

    samps = (
        pl.scan_csv("data/processed/samps.csv")
        .sort("tpt")
        .with_columns(pl.col("freq_idx").sub(1))
    )

    samps_trt = (
        samps.filter(pl.col("tpt").is_in([2, 3, 4, 5, 6]))
        .group_by(["cage", "treatment"])
        .agg([pl.col("freq_idx")])
        .with_columns(pl.lit("trt").alias("sweep"))
    )
    samps_post_trt = (
        samps.filter(pl.col("tpt").is_in([6, 7, 8]))
        .group_by(["cage", "treatment"])
        .agg([pl.col("freq_idx")])
        .with_columns(pl.lit("post_trt").alias("sweep"))
    )
    samps = pl.concat([samps_trt, samps_post_trt])

    lm_sites = (
        pl.scan_csv("data/raw/sigsite_malathion.csv")
        .drop(["treatment"])
        .rename(
            {
                "comparison": "sweep",
                "cage": "treatment",
                "estimate": "lm_slope",
                "sigLevel": "significance_level",
                "effect_size": "lm_effect",
            }
        )
        .filter([pl.col("sweep") != "1_2"])
        .with_columns(
            [
                pl.col("sweep").replace({"2_6": "trt", "6_8": "post_trt"}),
            ]
        )
        .select(
            [
                "chrom",
                "pos",
                "treatment",
                "sweep",
                "lm_slope",
                "significance_level",
                "lm_effect",
            ]
        )
    )

    sites = (
        sites.lazy()
        .drop(["initial_freq"])
        # relabel tpt to mean one of the two sweeps: trt or post_trt
        .rename({"tpt": "sweep"})
        .with_columns(
            pl.col("sweep").cast(pl.String).replace({2: "trt", 6: "post_trt"})
        )
        # the E2 cage has no data, so drop data for it
        .filter(~((pl.col("treatment") == "E") & (pl.col("cage") == 2)))
        .join(lm_sites, on=["chrom", "pos", "treatment", "sweep"], how="left")
        .drop_nulls()
        .join(samps, on=["cage", "treatment", "sweep"], how="left")
        .with_columns(
            pl.struct([pl.col("site_idx"), pl.col("freq_idx")]).alias("afmat_indices")
        )
        .with_columns(
            pl.col("afmat_indices")
            .map_elements(
                lambda s: list(afmat[s["site_idx"], s["freq_idx"]]),
                return_dtype=pl.List(pl.Float64),
            )
            .alias("freq")
        )
        .drop(["site_idx", "freq_idx", "afmat_indices"])
        # compute deltas: total and diff
        .with_columns(
            (pl.col("freq").list.last() - pl.col("freq").list.first()).alias(
                "total_delta"
            ),
            pl.col("freq").list.diff(1, null_behavior="drop").alias("delta"),
        )
        # .join(pl.from_pandas(sweep_r2s[['chrom','pos', 'r2']]).lazy(), on=["chrom", "pos"], how="left")
        .with_columns(pl.col("freq").list.first().alias("freq0"))
        .drop("delta", "freq")
    ).collect().write_parquet("data/sites_main.parquet")
