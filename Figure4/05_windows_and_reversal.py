#!/usr/bin/env python

import pandas as pd
import polars as pl
import numpy as np
import numpy.ma as ma
from tqdm import tqdm
import multiprocess as mp


def compute_windows(d, w=1e6):
    max_pos = d["pos"].max()
    breakpoints = np.arange(0, max_pos, w)
    breakpoints = list(breakpoints) + [breakpoints[-1] + w]
    windows = zip(breakpoints[:-1], breakpoints[1:])

    rows = []
    for start, end in tqdm(windows, total=len(breakpoints) - 1):
        dw = d.query("pos >= @start and pos <= @end")
        # print(dw.shape)
        lm_e = dw["lm_effect"]  # .abs()

        boots_lm = np.zeros(100_000)
        # bootstrap 1000 times, get quantiles and median
        for i in range(100_000):
            b_lm = np.random.choice(lm_e, size=len(lm_e), replace=True)
            boots_lm[i] = np.median(b_lm)

        rows.append(
            {
                "mid": (start + end) / 2,
                "lm_median": np.median(boots_lm),
                "lm_lower": np.quantile(boots_lm, 0.025),
                "lm_upper": np.quantile(boots_lm, 0.975),
                "nsnp": dw.shape[0],
            }
        )

    return pd.DataFrame(rows)


if __name__ == "__main__":
    d = pd.read_csv("plot_data/sites.csv").drop(columns=["Unnamed: 0"])
    dwin = (
        d.groupby(["treatment", "chrom", "link"])
        .apply(compute_windows)
        .reset_index()
        .drop(columns=["level_3"])
    )
    dwin.to_csv("plot_data/windows.csv")

    d_post = pd.read_csv("plot_data/sites_post.csv").drop(columns=["Unnamed: 0"])
    dwin_post = (
        d_post.groupby(["treatment", "chrom", "link"])
        .apply(compute_windows)
        .reset_index()
        .drop(columns=["level_3"])
    )
    dwin_post.to_csv("plot_data/windows_post.csv")

    # also process lm sites here to get reversal data
    lmd = lm_sites = (
        (
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
                    "lm_effect",
                ]
            )
        )
        .collect()
        .to_pandas()
        .drop_duplicates()
    )

    d = (
        pd.read_csv("plot_data/sites.csv")[["chrom", "pos", "link", "treatment"]]
        .drop_duplicates()
        .merge(lmd, on=["chrom", "pos", "treatment"], how="left")
        .dropna()
        .pivot(
            index=["chrom", "pos", "link", "treatment"],
            columns="sweep",
            values="lm_effect",
        )
        .reset_index()
        .sort_values(["treatment", "link", "chrom", "pos"])
        .drop_duplicates()
    )

    d.to_csv("plot_data/reversal.csv")


