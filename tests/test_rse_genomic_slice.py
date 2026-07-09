import pytest
import numpy as np
import pandas as pd
import genomicranges
from summarizedexperiment.RangedSummarizedExperiment import RangedSummarizedExperiment

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def test_RSE_genomic_slice():
    df_gr = pd.DataFrame(
        {
            "seqnames": ["chr1", "chr2", "chr3"],
            "starts": [100, 200, 300],
            "ends": [150, 250, 350],
            "strand": ["+", "-", "*"],
        }
    )
    gr = genomicranges.GenomicRanges.from_pandas(df_gr)

    rse = RangedSummarizedExperiment(
        assays={"counts": np.random.rand(3, 2)},
        row_ranges=gr,
    )

    query_df = pd.DataFrame(
        {
            "seqnames": ["chr2"],
            "starts": [220],
            "ends": [230],
            "strand": ["-"],
        }
    )
    query_gr = genomicranges.GenomicRanges.from_pandas(query_df)

    sliced = rse[query_gr, :]
    assert sliced is not None
    assert sliced.shape == (1, 2)
    assert sliced.row_ranges.get_seqnames() == ["chr2"]

    sliced_1d = rse[query_gr]
    assert sliced_1d.shape == (1, 2)
    assert sliced_1d.row_ranges.get_seqnames() == ["chr2"]
