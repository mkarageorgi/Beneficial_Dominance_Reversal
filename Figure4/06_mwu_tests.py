#!/usr/bin/env python
import pandas as pd
from scipy import stats

d = (
    pd.read_csv("plot_data/sites.csv")
    .drop(columns=["Unnamed: 0"])
    .query('treatment == "P" and chrom == "3R"')
)
d["lm_effect"] = -d["lm_effect"]

binsize = 1e6

# use the pos column to assign bin id to each site with 500000 bp bins
d["bin"] = (d["pos"] / binsize).astype(int)
d["binmid"] = d["bin"] * binsize + binsize / 2


def mwu(df):
    linked = df.query('link == "linked"')["lm_effect"]
    unlinked = df.query('link == "unlinked"')["lm_effect"]
    return pd.Series(stats.mannwhitneyu(linked, unlinked, alternative="greater").pvalue)


dmwuR = (
    d.groupby(["chrom", "bin", "binmid"])
    .apply(mwu)
    .reset_index()
    .rename(columns={0: "pval"})
)

d = (
    pd.read_csv("plot_data/sites.csv")
    .drop(columns=["Unnamed: 0"])
    .query('treatment == "P" and chrom == "3L"')
)
d["lm_effect"] = -d["lm_effect"]

# use the pos column to assign bin id to each site with 500000 bp bins
d["bin"] = (d["pos"] / binsize).astype(int)
d["binmid"] = d["bin"] * binsize + binsize / 2


def mwu(df):
    linked = df.query('link == "linked"')["lm_effect"]
    unlinked = df.query('link == "unlinked"')["lm_effect"]
    return pd.Series(stats.mannwhitneyu(linked, unlinked, alternative="greater").pvalue)


dmwuL = (
    d.groupby(["chrom", "bin", "binmid"])
    .apply(mwu)
    .reset_index()
    .rename(columns={0: "pval"})
)

d = (
    pd.read_csv("plot_data/sites_post.csv")
    .drop(columns=["Unnamed: 0"])
    .query('treatment == "P" and chrom == "3R"')
)

d["lm_effect"] = -d["lm_effect"]

# use the pos column to assign bin id to each site with 500000 bp bins
d["bin"] = (d["pos"] / binsize).astype(int)
d["binmid"] = d["bin"] * binsize + binsize / 2


def mwu(df):
    linked = df.query('link == "linked"')["lm_effect"]
    unlinked = df.query('link == "unlinked"')["lm_effect"]
    return pd.Series(stats.mannwhitneyu(linked, unlinked, alternative="less").pvalue)


dmwuRev = (
    d.groupby(["chrom", "bin", "binmid"])
    .apply(mwu)
    .reset_index()
    .rename(columns={0: "pval"})
)

pd.concat(
    [
        dmwuL.assign(sweep="trt", chrom="3L"),
        dmwuR.assign(sweep="trt", chrom="3R"),
        dmwuRev.assign(sweep="rev", chrom="3R"),
    ]
).to_csv("plot_data/mwu.csv", index=False)
