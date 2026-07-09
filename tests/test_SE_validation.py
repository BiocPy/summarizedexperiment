from copy import deepcopy

import numpy as np
import pytest
from biocframe import BiocFrame
from genomicranges import GenomicRanges
from iranges import IRanges

from summarizedexperiment import RangedSummarizedExperiment, SummarizedExperiment

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def test_assays_validation():
    with pytest.raises(Exception):
        SummarizedExperiment(assays=123)

    with pytest.raises(TypeError):
        SummarizedExperiment(assays={"counts": "invalid"})

    with pytest.raises(ValueError):
        SummarizedExperiment(assays={"counts": np.random.rand(2, 2, 2)})

    with pytest.raises(ValueError):
        SummarizedExperiment(assays={"counts": np.random.rand(2, 3), "log": np.random.rand(2, 4)})


def test_rows_cols_validation():
    with pytest.raises(TypeError):
        SummarizedExperiment(assays={"counts": np.random.rand(2, 2)}, row_data="invalid")

    with pytest.raises(ValueError):
        SummarizedExperiment(assays={"counts": np.random.rand(2, 2)}, row_data=BiocFrame({}, number_of_rows=3))

    with pytest.raises(ValueError):
        SummarizedExperiment(assays={"counts": np.random.rand(2, 2)}, row_names=["a", "b", "c"])

    with pytest.raises(TypeError):
        SummarizedExperiment(assays={"counts": np.random.rand(2, 2)}, column_data="invalid")

    with pytest.raises(ValueError):
        SummarizedExperiment(assays={"counts": np.random.rand(2, 2)}, column_data=BiocFrame({}, number_of_rows=3))

    with pytest.raises(ValueError):
        SummarizedExperiment(assays={"counts": np.random.rand(2, 2)}, column_names=["a", "b", "c"])


def test_deepcopy():
    se = SummarizedExperiment(
        assays={"counts": np.random.rand(2, 2)},
        row_names=["gene1", "gene2"],
        column_names=["sample1", "sample2"],
        metadata={"study": "cancer"},
    )

    se_copy = deepcopy(se)
    assert se_copy is not se
    assert se_copy.shape == se.shape
    assert list(se_copy.row_names) == list(se.row_names)
    assert list(se_copy.column_names) == list(se.column_names)
    assert se_copy.metadata["study"] == "cancer"

    gr = GenomicRanges(seqnames=["chr1", "chr2"], ranges=IRanges([1, 100], [10, 50]), strand=["+", "-"])
    rse = RangedSummarizedExperiment(
        assays={"counts": np.random.rand(2, 2)},
        row_ranges=gr,
        column_names=["sample1", "sample2"],
        metadata={"study": "cancer"},
    )
    rse_copy = deepcopy(rse)
    assert rse_copy is not rse
    assert rse_copy.shape == rse.shape
    assert len(rse_copy.row_ranges) == len(rse.row_ranges)


def test_properties_get_set():
    se = SummarizedExperiment(
        assays={"counts": np.random.rand(2, 2)}, row_names=["gene1", "gene2"], column_names=["sample1", "sample2"]
    )

    se.assays = {"counts2": np.random.rand(2, 2)}
    assert "counts2" in se.assay_names

    row_df = BiocFrame({"feat": [1, 2]})
    se.row_data = row_df
    assert "feat" in se.row_data.colnames
    se.rowdata = row_df
    assert "feat" in se.rowdata.colnames

    col_df = BiocFrame({"sample_id": [1, 2]})
    se.col_data = col_df
    assert "sample_id" in se.col_data.colnames
    se.coldata = col_df
    assert "sample_id" in se.coldata.colnames
    se.columndata = col_df
    assert "sample_id" in se.columndata.colnames

    rdata = se.get_row_data(replace_row_names=False)
    cdata = se.get_column_data(replace_row_names=False)
    assert rdata.row_names is None
    assert cdata.row_names is None

    new_r = BiocFrame({"feat": [3, 4]}, row_names=["g1", "g2"])
    se2 = se.set_row_data(new_r, replace_row_names=True)
    assert list(se2.row_names) == ["g1", "g2"]

    new_c = BiocFrame({"sample_id": [3, 4]}, row_names=["s1", "s2"])
    se3 = se.set_column_data(new_c, replace_column_names=True)
    assert list(se3.column_names) == ["s1", "s2"]


def test_hasattr_overlaps_slice():
    import unittest.mock as mock

    gr = GenomicRanges(seqnames=["chr1", "chr2"], ranges=IRanges([1, 100], [10, 50]), strand=["+", "-"])
    rse = RangedSummarizedExperiment(
        assays={"counts": np.random.rand(2, 2)}, row_ranges=gr, column_names=["sample1", "sample2"]
    )

    class OverlapMock:
        def find_overlaps(self, query):
            pass

    with mock.patch.object(rse.row_ranges, "find_overlaps", return_value=BiocFrame({"self_hits": [0]})):
        res = rse[OverlapMock(), :]
        assert res.shape == (1, 2)
