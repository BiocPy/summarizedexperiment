from random import random

import genomicranges
import numpy as np
import pandas as pd
import pytest
from genomicranges import GenomicRanges, IRanges

from summarizedexperiment.RangedSummarizedExperiment import RangedSummarizedExperiment
from summarizedexperiment.SummarizedExperiment import SummarizedExperiment

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


nrows = 200
ncols = 6
counts = np.random.rand(nrows, ncols)
df_gr = pd.DataFrame({
    "seqnames": [
        "chr1",
        "chr2",
        "chr2",
        "chr2",
        "chr1",
        "chr1",
        "chr3",
        "chr3",
        "chr3",
        "chr3",
    ]
    * 20,
    "starts": range(100, 300),
    "ends": range(110, 310),
    "strand": ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"] * 20,
    "score": range(0, 200),
    "GC": [random() for _ in range(10)] * 20,
})

gr = genomicranges.GenomicRanges.from_pandas(df_gr)

col_data = pd.DataFrame({
    "treatment": ["ChIP", "Input"] * 3,
})


def test_RSE_props():
    tse = RangedSummarizedExperiment(assays={"counts": counts}, row_ranges=gr, column_data=col_data)

    assert tse is not None
    assert isinstance(tse, RangedSummarizedExperiment)

    assert tse.assay_names is not None
    assert len(tse.assay_names) == 1

    assert tse.col_data is not None
    assert tse.row_ranges is not None

    assert tse.dims == tse.shape

    assert tse.metadata is not None

    assert tse.start is not None
    assert tse.end is not None

    assert tse.seq_info is not None
    assert tse.seqnames is not None
    assert tse.strand is not None
    assert tse.width is not None

    assert tse.rownames is None
    assert tse.colnames is not None


def test_RSE_subset():
    tse = RangedSummarizedExperiment(assays={"counts": counts}, row_ranges=gr, column_data=col_data)

    assert tse is not None
    assert isinstance(tse, RangedSummarizedExperiment)

    subset_tse = tse[0:10, 2:5]
    assert subset_tse is not None
    assert isinstance(subset_tse, RangedSummarizedExperiment)

    assert len(subset_tse.row_ranges) == 10
    assert len(subset_tse.col_data) == 3

    assert subset_tse.assay("counts").shape == (10, 3)


def test_RSE_subset_assays():
    tse = RangedSummarizedExperiment(assays={"counts": counts}, row_ranges=gr, column_data=col_data)

    assert tse is not None
    assert isinstance(tse, RangedSummarizedExperiment)

    subset_asys = tse.subset_assays(rows=slice(1, 10), columns=[0, 1, 2])
    assert subset_asys is not None
    assert isinstance(subset_asys, type(tse.assays))

    assert len(subset_asys.keys()) == 1
    assert subset_asys["counts"].shape == (9, 3)


def test_RSE_coverage():
    tse = RangedSummarizedExperiment(assays={"counts": counts}, row_ranges=gr, column_data=col_data)

    assert tse is not None
    assert isinstance(tse, RangedSummarizedExperiment)

    cov = tse.coverage()
    assert cov is not None
    assert len(cov.keys()) == 3


def test_RSE_distance_methods():
    tse = RangedSummarizedExperiment(assays={"counts": counts}, row_ranges=gr, column_data=col_data)

    assert tse is not None
    assert isinstance(tse, RangedSummarizedExperiment)

    nearest = tse.nearest(tse)
    assert nearest is not None

    precede = tse.precede(tse)
    assert precede is not None

    follow = tse.follow(tse)
    assert follow is not None


def test_RSE_range_methods():
    tse = RangedSummarizedExperiment(assays={"counts": counts}, row_ranges=gr, column_data=col_data)

    assert tse is not None
    assert isinstance(tse, RangedSummarizedExperiment)

    res = tse.flank(width=10)
    assert res is not None

    res = tse.resize(width=15)
    assert res is not None

    res = tse.restrict(start=400)
    assert res is not None

    res = tse.shift(shift=25)
    assert res is not None

    res = tse.promoters()
    assert res is not None

    res = tse.narrow(start=40)
    assert res is not None


def test_RSE_search():
    tse = RangedSummarizedExperiment(assays={"counts": counts}, row_ranges=gr, column_data=col_data)

    assert tse is not None
    assert isinstance(tse, RangedSummarizedExperiment)

    res = tse.find_overlaps(tse)
    assert res is not None

    res = tse.subset_by_overlaps(tse)
    assert res is not None


def test_RSE_sort_order():
    tse = RangedSummarizedExperiment(assays={"counts": counts}, row_ranges=gr, column_data=col_data)

    assert tse is not None
    assert isinstance(tse, RangedSummarizedExperiment)

    res = tse.order()
    assert res is not None
    assert len(res) == tse.shape[0]

    res = tse.sort(decreasing=True)
    assert res is not None
    assert len(res.row_ranges) == len(tse.row_ranges)


def test_RSE_assay_getters_and_setters():
    tse = RangedSummarizedExperiment(assays={"counts": counts}, row_ranges=gr, column_data=col_data)

    assert tse is not None
    assert isinstance(tse, RangedSummarizedExperiment)

    assert tse.assay(0) is not None

    new_tse = tse.set_assay("new_counts", assay=np.random.rand(nrows, ncols), in_place=False)
    assert new_tse.get_assay_names() != tse.get_assay_names()
    with pytest.raises(Exception):
        tse.get_assay("new_counts")
    assert new_tse.get_assay("new_counts") is not None

    tse.set_assay("new_counts", assay=np.random.rand(nrows, ncols), in_place=True)
    assert new_tse.get_assay_names() == tse.get_assay_names()
    assert tse.get_assay("new_counts") is not None
    assert new_tse.get_assay("new_counts") is not None


def test_RSE_to_se():
    rse = RangedSummarizedExperiment(assays={"counts": counts}, row_ranges=gr, column_data=col_data)

    se = rse.to_summarizedexperiment()

    assert se is not None
    assert isinstance(se, SummarizedExperiment)
    assert rse.shape == se.shape


def test_RSE_to_se_with_ranges():
    tse = RangedSummarizedExperiment(assays={"counts": counts}, row_ranges=gr, column_data=col_data)

    se = tse.to_summarizedexperiment()

    assert se is not None
    assert isinstance(se, SummarizedExperiment)
    assert se.shape == tse.shape

    assert se.row_data is not None
    assert "seqnames" in se.row_data.column_names
    assert "starts" in se.row_data.column_names
    assert "ends" in se.row_data.column_names


def test_RSE_combine():
    gr1 = GenomicRanges(seqnames=["chr1", "chr2"], ranges=IRanges([1, 100], [10, 50]), strand=["+", "-"])
    rse1 = RangedSummarizedExperiment(
        assays={"counts": np.array([[1, 2], [3, 4]])}, row_ranges=gr1, column_names=["sample1", "sample2"]
    )

    gr2 = GenomicRanges(seqnames=["chr3"], ranges=IRanges([200], [20]), strand=["*"])
    rse2 = RangedSummarizedExperiment(
        assays={"counts": np.array([[5, 6]])}, row_ranges=gr2, column_names=["sample1", "sample2"]
    )

    combined_rows = rse1.combine_rows(rse2)
    assert combined_rows.shape == (3, 2)
    assert list(combined_rows.row_ranges.get_seqnames()) == ["chr1", "chr2", "chr3"]

    combined_rows_rel = rse1.relaxed_combine_rows(rse2)
    assert combined_rows_rel.shape == (3, 2)

    rse3 = RangedSummarizedExperiment(assays={"counts": np.array([[7], [8]])}, row_ranges=gr1, column_names=["sample3"])
    combined_cols = rse1.combine_columns(rse3)
    assert combined_cols.shape == (2, 3)
    assert list(combined_cols.column_names) == ["sample1", "sample2", "sample3"]

    combined_cols_rel = rse1.relaxed_combine_columns(rse3)
    assert combined_cols_rel.shape == (2, 3)


def test_RSE_genomic_slice():
    df_gr = pd.DataFrame({
        "seqnames": ["chr1", "chr2", "chr3"],
        "starts": [100, 200, 300],
        "ends": [150, 250, 350],
        "strand": ["+", "-", "*"],
    })
    gr = genomicranges.GenomicRanges.from_pandas(df_gr)

    rse = RangedSummarizedExperiment(
        assays={"counts": np.random.rand(3, 2)},
        row_ranges=gr,
    )

    query_df = pd.DataFrame({
        "seqnames": ["chr2"],
        "starts": [220],
        "ends": [230],
        "strand": ["-"],
    })
    query_gr = genomicranges.GenomicRanges.from_pandas(query_df)

    sliced = rse[query_gr, :]
    assert sliced is not None
    assert sliced.shape == (1, 2)
    assert sliced.row_ranges.get_seqnames() == ["chr2"]

    sliced_1d = rse[query_gr]
    assert sliced_1d.shape == (1, 2)
    assert sliced_1d.row_ranges.get_seqnames() == ["chr2"]
