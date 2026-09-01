from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Literal, Union
from warnings import warn

import biocframe
import biocutils as ut
import numpy as np
from genomicranges import CompressedGenomicRangesList, GenomicRanges, SeqInfo

from ._combineutils import (
    check_assays_are_equal,
    merge_assays,
    merge_se_colnames,
    merge_se_rownames,
    relaxed_merge_assays,
)
from ._frameutils import is_pandas
from .SummarizedExperiment import SummarizedExperiment

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

GRangesOrGRangesList = Union[GenomicRanges, CompressedGenomicRangesList]
GRangesOrRangeSE = Union[GRangesOrGRangesList, "RangedSummarizedExperiment"]


def _check_gr_or_rse(x: GRangesOrRangeSE):
    """Check if ``x`` is either a `RangedSummarizedExperiment` or `GenomicRanges`.

    Args:
        x:
            Object to check.

    Raises:
        TypeError:
            Object is not a `RangedSummarizedExperiment` or `GenomicRanges`.
    """
    if not isinstance(x, (RangedSummarizedExperiment, GenomicRanges, CompressedGenomicRangesList)):
        raise TypeError("'x' is not a `RangedSummarizedExperiment`, `GenomicRanges` or `CompressedGenomicRangesList`.")

    if isinstance(x, CompressedGenomicRangesList):
        raise NotImplementedError("'range' operations are not implemented for `CompressedGenomicRangesList`.")


def _access_granges(x: GRangesOrRangeSE) -> GenomicRanges:
    """Access ranges from the object.

    Args:
        x:
            Input object.

    Returns:
        `GenomicRanges` object.
    """
    qranges = x
    if isinstance(x, RangedSummarizedExperiment):
        qranges = x.row_ranges

    return qranges


def _validate_rowranges(row_ranges, shape):
    if not (isinstance(row_ranges, (GenomicRanges, CompressedGenomicRangesList))):
        raise TypeError(
            f"`row_ranges` must be a `GenomicRanges` or `CompressedGenomicRangesList` , provided {type(row_ranges)}."
        )

    if len(row_ranges) != shape[0]:
        raise ValueError("Number of features in `row_ranges` and number of rows in assays do not match.")


def _sanitize_ranges_frame(frame, num_rows: int):
    frame = frame if frame is not None else CompressedGenomicRangesList.empty(n=num_rows)

    if is_pandas(frame):
        frame = GenomicRanges.from_pandas(frame)

    return frame


class RangedSummarizedExperiment(SummarizedExperiment):
    """RangedSummarizedExperiment class to represent genomic experiment data, genomic features as
    :py:class:`~genomicranges.GenomicRanges.GenomicRanges` or
    :py:class:`~genomicranges.CompressedGenomicRangesList.CompressedGenomicRangesList` sample data and any additional experimental metadata.

    Note: If ``row_ranges`` is empty, None or not a
    :py:class:`genomicranges.GenomicRanges.GenomicRanges` object, use a
    :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment` instead.
    """

    def __init__(
        self,
        assays: dict[str, Any] = None,
        row_ranges: GRangesOrGRangesList | None = None,
        row_data: biocframe.BiocFrame | None = None,
        column_data: biocframe.BiocFrame | None = None,
        row_names: list[str] | None = None,
        column_names: list[str] | None = None,
        metadata: dict[str, Any] | ut.NamedList | None = None,
        _validate: bool = True,
    ) -> None:
        """Initialize a `RangedSummarizedExperiment` (RSE) object.

        Args:
            assays:
                A dictionary containing matrices, with assay names as keys
                and 2-dimensional matrices represented as either
                :py:class:`~numpy.ndarray` or :py:class:`~scipy.sparse.spmatrix`.

                Alternatively, you may use any 2-dimensional matrix that has
                the ``shape`` property and implements the slice operation
                using the ``__getitem__`` dunder method.

                All matrices in assays must be 2-dimensional and have the
                same shape (number of rows, number of columns).

            row_ranges:
                Genomic features, must be the same length as the number of rows of
                the matrices in assays.

            row_data:
                Features, must be the same length as the number of rows of
                the matrices in assays.

                Feature information is coerced to a
                :py:class:`~biocframe.BiocFrame.BiocFrame`. Defaults to None.

            column_data:
                Sample data, must be the same length as the number of
                columns of the matrices in assays.

                Sample information is coerced to a
                :py:class:`~biocframe.BiocFrame.BiocFrame`. Defaults to None.

            row_names:
                A list of strings, same as the number of rows.

                If ``row_names`` are not provided, these are inferred from
                ``row_data``.

                Defaults to None.

            column_names:
                A list of string, same as the number of columns.

                if ``column_names`` are not provided, these are inferred from
                ``column_data``.

                Defaults to None.

            metadata:
                Additional experimental metadata describing the methods.
                Defaults to None.

            _validate:
                Internal use only.
        """
        super().__init__(
            assays,
            row_data=row_data,
            column_data=column_data,
            row_names=row_names,
            column_names=column_names,
            metadata=metadata,
            _validate=_validate,
        )

        self._row_ranges = _sanitize_ranges_frame(row_ranges, self._shape[0])
        if _validate:
            _validate_rowranges(self._row_ranges, self._shape)

    #########################
    ######>> Copying <<######
    #########################

    def __deepcopy__(self, memo=None, _nil=[]):
        """
        Returns:
            A deep copy of the current ``RangedSummarizedExperiment``.
        """
        from copy import deepcopy

        _assays_copy = deepcopy(self._assays)
        _rows_copy = deepcopy(self._rows)
        _rowranges_copy = deepcopy(self._row_ranges)
        _cols_copy = deepcopy(self._cols)
        _row_names_copy = deepcopy(self._row_names)
        _col_names_copy = deepcopy(self._column_names)
        _metadata_copy = deepcopy(self.metadata)

        current_class_const = type(self)
        return current_class_const(
            assays=_assays_copy,
            row_ranges=_rowranges_copy,
            row_data=_rows_copy,
            column_data=_cols_copy,
            row_names=_row_names_copy,
            column_names=_col_names_copy,
            metadata=_metadata_copy,
            _validate=False,
        )

    def __copy__(self):
        """
        Returns:
            A shallow copy of the current ``RangedSummarizedExperiment``.
        """
        current_class_const = type(self)
        return current_class_const(
            assays=self._assays,
            row_ranges=self._row_ranges,
            row_data=self._rows,
            column_data=self._cols,
            row_names=self._row_names,
            column_names=self._column_names,
            metadata=self._metadata,
            _validate=False,
        )

    def copy(self):
        """Alias for :py:meth:`~__copy__`."""
        return self.__copy__()

    ##########################
    ######>> Printing <<######
    ##########################

    def __repr__(self) -> str:
        """
        Returns:
            A string representation.
        """
        output = f"{type(self).__name__}(number_of_rows={self.shape[0]}"
        output += f", number_of_columns={self.shape[1]}"
        output += ", assays=" + ut.print_truncated_list(self.assay_names)

        output += ", row_data=" + self._rows.__repr__()
        if self._row_names is not None:
            output += ", row_names=" + ut.print_truncated_list(self._row_names)

        output += ", column_data=" + self._cols.__repr__()
        if self._column_names is not None:
            output += ", column_names=" + ut.print_truncated_list(self._column_names)

        if self._row_ranges is not None:
            output += "row_ranges=" + self._row_ranges.__repr__()

        if len(self._metadata) > 0:
            output += ", metadata=" + ut.print_truncated_dict(self._metadata)

        output += ")"
        return output

    def __str__(self) -> str:
        """
        Returns:
            A pretty-printed string containing the contents of this object.
        """
        output = f"class: {type(self).__name__}\n"

        output += f"dimensions: ({self.shape[0]}, {self.shape[1]})\n"

        output += f"assays({len(self.assay_names)}): {ut.print_truncated_list(self.assay_names)}\n"

        output += (
            f"row_data columns({len(self._rows.column_names)}): {ut.print_truncated_list(self._rows.column_names)}\n"
        )
        output += f"row_names({0 if self._row_names is None else len(self._row_names)}): {' ' if self._row_names is None else ut.print_truncated_list(self._row_names)}\n"

        output += (
            f"column_data columns({len(self._cols.column_names)}): {ut.print_truncated_list(self._cols.column_names)}\n"
        )
        output += f"column_names({0 if self._column_names is None else len(self._column_names)}): {' ' if self._column_names is None else ut.print_truncated_list(self._column_names)}\n"

        output += f"metadata({len(self.metadata)!s}): {ut.print_truncated_list(list(self.metadata.keys()), sep=' ', include_brackets=False, transform=lambda y: y)}"

        return output

    ############################
    ######>> row_ranges <<######
    ############################

    def get_row_ranges(self) -> GRangesOrGRangesList:
        """Get genomic feature information.

        Returns:
            Genomic feature information.
        """
        return self._row_ranges

    def set_row_ranges(
        self, row_ranges: GRangesOrGRangesList | None, in_place: bool = False
    ) -> RangedSummarizedExperiment:
        """Set new genomic features.

        Args:
            row_ranges:
                Genomic features, must be the same length as the number of rows of
                the matrices in assays.

            in_place:
                Whether to modify the ``RangeSummarizedExperiment`` in place.

        Returns:
            A modified ``RangeSummarizedExperiment`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        rows = _sanitize_ranges_frame(row_ranges, self._shape[0])
        _validate_rowranges(rows, self._shape)

        output = self._define_output(in_place)
        output._row_ranges = rows
        return output

    @property
    def row_ranges(self) -> GRangesOrGRangesList:
        """Alias for :py:meth:`~get_rowranges`."""
        return self.get_row_ranges()

    @row_ranges.setter
    def row_ranges(self, row_ranges: GRangesOrGRangesList) -> None:
        """Alias for :py:meth:`~set_rowranges`."""
        warn(
            "Setting property 'row_ranges' is an in-place operation, use 'set_rowranges' instead",
            UserWarning,
        )
        return self.set_row_ranges(row_ranges=row_ranges, in_place=True)

    #################################
    ######>> range accessors <<######
    #################################

    @property
    def end(self) -> np.ndarray:
        """Get genomic end positions for each feature or row in experimental data.

        Returns:
            A :py:class:`numpy.ndarray` of end positions.
        """
        return self.row_ranges.end

    @property
    def start(self) -> np.ndarray:
        """Get genomic start positions for each feature or row in experimental data.

        Returns:
            A :py:class:`numpy.ndarray` of start positions.
        """
        return self.row_ranges.start

    @property
    def seqnames(self) -> list[str]:
        """Get sequence or chromosome names.

        Returns:
            List of all chromosome names.
        """
        return self.row_ranges.seqnames

    @property
    def strand(self) -> np.ndarray:
        """Get strand information.

        Returns:
            A :py:class:`numpy.ndarray` of strand information.
        """
        return self.row_ranges.strand

    @property
    def width(self) -> np.ndarray:
        """Get widths of ``row_ranges``.

        Returns:
            A :py:class:`numpy.ndarray` of widths for each interval.
        """
        return self.row_ranges.width

    @property
    def seq_info(self) -> SeqInfo:
        """Get sequence information object (if available).

        Returns:
            Sequence information.
        """
        return self.row_ranges.seqinfo

    ##########################
    ######>> slicers <<#######
    ##########################

    # rest of them are inherited from BaseSE.

    def _normalize_row_slice(self, rows: str | int | bool | Sequence):

        if isinstance(rows, (GenomicRanges, CompressedGenomicRangesList)) or hasattr(rows, "find_overlaps"):
            hits = self.row_ranges.find_overlaps(query=rows)
            rows = hits.get_column("self_hits")

        return super()._normalize_row_slice(rows)

    def get_slice(
        self,
        rows: str | int | bool | Sequence | None,
        columns: str | int | bool | Sequence | None,
    ) -> RangedSummarizedExperiment:
        """Alias for :py:attr:`~__getitem__`, for back-compatibility."""

        slicer = self._generic_slice(rows=rows, columns=columns)

        new_row_ranges = None
        if not (isinstance(slicer.row_indices, slice) and slicer.row_indices == slice(None)):
            new_row_ranges = ut.subset_sequence(self.row_ranges, slicer.row_indices)

        current_class_const = type(self)
        return current_class_const(
            assays=slicer.assays,
            row_ranges=new_row_ranges,
            row_data=slicer.rows,
            column_data=slicer.columns,
            row_names=slicer.row_names,
            column_names=slicer.column_names,
            metadata=self._metadata,
        )

    ############################
    ######>> range ops <<#######
    ############################

    def coverage(self, shift: int = 0, width: int | None = None, weight: int = 1) -> dict[str, np.ndarray]:
        """Calculate coverage for each chromosome.

        Args:
            shift:
                Shift all genomic positions. Defaults to 0.

            width:
                Restrict the width of all
                chromosomes. Defaults to None.

            weight:
                Weight to use. Defaults to 1.

        Returns:
            A dictionary with chromosome names as keys and the
            coverage vector as value.
        """
        return self.row_ranges.coverage(shift=shift, width=width, weight=weight)

    def nearest(
        self,
        query: GRangesOrRangeSE,
        select: Literal["all", "arbitrary"] = "all",
        ignore_strand: bool = False,
    ) -> list[int | None] | None:
        """Search nearest positions both upstream and downstream that overlap with each range in ``query``.

        Args:
            query:
                Query intervals to find nearest positions.

                ``query`` may be a :py:class:`~genomicranges.GenomicRanges.GenomicRanges` or a
                :py:class:`~summarizedexperiment.RangedSummarizedExperiment.RangedSummarizedExperiment`
                object.

            select:
                Determine what hit to choose when there are
                multiple hits for an interval in ``query``.

            ignore_strand:
                Whether to ignore strand. Defaults to False.

        Raises:
            TypeError:
                If ``query`` is not a ``RangedSummarizedExperiment``
                or ``GenomicRanges``.

        Returns:
            A list with the same length as ``query``,
            containing hits to nearest indices.
        """

        _check_gr_or_rse(query)

        qranges = _access_granges(query)

        res = self.row_ranges.nearest(query=qranges, select=select, ignore_strand=ignore_strand)
        return res

    def precede(
        self,
        query: GRangesOrRangeSE,
        select: Literal["all", "arbitrary"] = "all",
        ignore_strand: bool = False,
    ) -> list[int | None] | None:
        """Search nearest positions only downstream that overlap with each range in ``query``.

        Args:
            query:
                Query intervals to find nearest positions.

                ``query`` may be a :py:class:`~genomicranges.GenomicRanges.GenomicRanges` or a
                :py:class:`~summarizedexperiment.RangedSummarizedExperiment.RangedSummarizedExperiment`
                object.

            select:
                Determine what hit to choose when there are
                multiple hits for an interval in ``query``.

            ignore_strand:
                Whether to ignore strand. Defaults to False.

        Raises:
            TypeError:
                If ``query`` is not a ``RangedSummarizedExperiment`` or
                ``GenomicRanges``.

        Returns:
            A List with the same length as ``query``,
            containing hits to nearest indices.
        """

        _check_gr_or_rse(query)

        qranges = _access_granges(query)

        res = self.row_ranges.precede(query=qranges, select=select, ignore_strand=ignore_strand)
        return res

    def follow(
        self,
        query: GRangesOrRangeSE,
        select: Literal["all", "arbitrary"] = "all",
        ignore_strand: bool = False,
    ) -> list[int | None] | None:
        """Search nearest positions only upstream that overlap with each range in ``query``.

        Args:
            query:
                Query intervals to find nearest positions.

                ``query`` may be a :py:class:`~genomicranges.GenomicRanges.GenomicRanges` or a
                :py:class:`~summarizedexperiment.RangedSummarizedExperiment.RangedSummarizedExperiment`
                object.

            select:
                Determine what hit to choose when there are
                multiple hits for an interval in ``query``.

            ignore_strand:
                Whether to ignore strand. Defaults to False.

        Raises:
            If ``query`` is not a ``RangedSummarizedExperiment`` or
            ``GenomicRanges``.

        Returns:
            A List with the same length as ``query``,
            containing hits to nearest indices.
        """
        _check_gr_or_rse(query)

        qranges = _access_granges(query)

        res = self.row_ranges.follow(query=qranges, select=select, ignore_strand=ignore_strand)
        return res

    def flank(
        self,
        width: int,
        start: bool = True,
        both: bool = False,
        ignore_strand: bool = False,
        in_place: bool = False,
    ) -> RangedSummarizedExperiment:
        """Compute flanking ranges for each range.

        Refer to either :py:meth:`~genomicranges.GenomicRanges.GenomicRanges.flank` or the
        Bioconductor documentation for more details.

        Args:
            width:
                Width to flank by. May be negative.

            start:
                Whether to only flank starts. Defaults to True.

            both:
                Whether to flank both starts and ends. Defaults to False.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A new `RangedSummarizedExperiment` object with the flanked ranges,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """
        new_ranges = self.row_ranges.flank(
            width=width,
            start=start,
            both=both,
            ignore_strand=ignore_strand,
            in_place=False,
        )

        output = self._define_output(in_place=in_place)
        output._row_ranges = new_ranges
        return output

    def resize(
        self,
        width: int | list[int] | np.ndarray,
        fix: Literal["start", "end", "center"] = "start",
        ignore_strand: bool = False,
        in_place: bool = False,
    ) -> RangedSummarizedExperiment:
        """Resize ranges to the specified ``width`` where either the ``start``, ``end``, or ``center`` is used as an
        anchor.

        Args:
            width:
                Width to resize, cannot be negative!

            fix:
                Fix positions by "start", "end", or "center".
                Defaults to "start".

            ignore_strand:
                Whether to ignore strands. Defaults to False.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Raises:
            ValueError:
                If ``fix`` is neither ``start``, ``center``, or ``end``.

        Returns:
            A new `RangedSummarizedExperiment` object with the resized ranges,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """
        new_ranges = self.row_ranges.resize(width=width, fix=fix, ignore_strand=ignore_strand, in_place=False)

        output = self._define_output(in_place=in_place)
        output._row_ranges = new_ranges
        return output

    def shift(self, shift: int | list[int] | np.ndarray = 0, in_place: bool = False) -> RangedSummarizedExperiment:
        """Shift all intervals.

        ``shift`` may be be negative.

        Args:
            shift:
                Shift interval. If shift is 0, the current
                object is returned. Defaults to 0.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A new `RangedSummarizedExperiment` object with the shifted ranges,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """
        new_ranges = self.row_ranges.shift(shift=shift, in_place=False)

        output = self._define_output(in_place=in_place)
        output._row_ranges = new_ranges
        return output

    def promoters(
        self, upstream: int = 2000, downstream: int = 200, in_place: bool = False
    ) -> RangedSummarizedExperiment:
        """Extend intervals to promoter regions.

        Args:
            upstream:
                Number of positions to extend in the 5' direction.
                Defaults to 2000.

            downstream:
                Number of positions to extend in the 3' direction.
                Defaults to 200.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A new `RangedSummarizedExperiment` object with the extended ranges for
            promoter regions, either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """
        new_ranges = self.row_ranges.promoters(upstream=upstream, downstream=downstream, in_place=False)

        output = self._define_output(in_place=in_place)
        output._row_ranges = new_ranges
        return output

    def restrict(
        self,
        start: int | list[int] | np.ndarray | None = None,
        end: int | list[int] | np.ndarray | None = None,
        keep_all_ranges: bool = False,
        in_place: bool = False,
    ) -> RangedSummarizedExperiment:
        """Restrict ranges to a given start and end positions.

        Args:
            start:
                Start position. Defaults to None.

            end:
                End position. Defaults to None.

            keep_all_ranges:
                Whether to keep intervals that do not overlap with start and end.
                Defaults to False.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A new `RangedSummarizedExperiment` object with restricted intervals,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """
        new_ranges = self.row_ranges.restrict(start=start, end=end, keep_all_ranges=keep_all_ranges)

        output = self._define_output(in_place=in_place)
        output._row_ranges = new_ranges
        return output

    def narrow(
        self,
        start: int | list[int] | np.ndarray | None = None,
        width: int | list[int] | np.ndarray | None = None,
        end: int | list[int] | np.ndarray | None = None,
        in_place: bool = False,
    ) -> RangedSummarizedExperiment:
        """Narrow genomic positions by provided ``start``, ``width`` and ``end`` parameters.

        Important: these parameters are relative shift in positions for each range.

        Args:
            start:
                Relative start position. Defaults to None.

            width:
                Relative end position. Defaults to None.

            end:
                Relative width of the interval. Defaults to None.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A new `RangedSummarizedExperiment` object with narrow positions,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """
        new_ranges = self.row_ranges.narrow(start=start, width=width, end=end, in_place=False)

        output = self._define_output(in_place=in_place)
        output._row_ranges = new_ranges
        return output

    def find_overlaps(
        self,
        query: GRangesOrRangeSE,
        query_type: Literal["any", "start", "end", "within"] = "any",
        select: Literal["all", "first", "last", "arbitrary"] = "all",
        max_gap: int = -1,
        min_overlap: int = 1,
        ignore_strand: bool = False,
    ) -> biocframe.BiocFrame:
        """Find overlaps between subject (self) and query ranges.

        Args:
            query:
                Query intervals to find nearest positions.

                ``query`` may be a :py:class:`~genomicranges.GenomicRanges.GenomicRanges` or a
                :py:class:`~summarizedexperiment.RangedSummarizedExperiment.RangedSummarizedExperiment`
                object.

            query_type:
                Overlap query type, must be one of

                - "any": Any overlap is good
                - "start": Overlap at the beginning of the intervals
                - "end": Must overlap at the end of the intervals
                - "within": Fully contain the query interval

                Defaults to "any".

            select:
                Determine what hit to choose when
                there are multiple hits for an interval in ``subject``.

            max_gap:
                Maximum gap allowed in the overlap.
                Defaults to -1 (no gap allowed).

            min_overlap:
                Minimum overlap with query. Defaults to 1.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Raises:
            TypeError:
                If `query` is not a `RangedSummarizedExperiment` or `GenomicRanges`.

        Returns:
            A list with the same length as ``query``,
            containing hits to overlapping indices.
        """
        _check_gr_or_rse(query)

        qranges = _access_granges(query)

        return self.row_ranges.find_overlaps(
            query=qranges,
            query_type=query_type,
            select=select,
            max_gap=max_gap,
            min_overlap=min_overlap,
            ignore_strand=ignore_strand,
        )

    def subset_by_overlaps(
        self,
        query: GRangesOrRangeSE,
        query_type: Literal["any", "start", "end", "within"] = "any",
        max_gap: int = -1,
        min_overlap: int = 1,
        ignore_strand: bool = False,
    ) -> RangedSummarizedExperiment:
        """Subset a `RangedSummarizedExperiment` by feature overlaps.

        Args:
            query:
                Query `GenomicRanges`.

                ``query`` may be a :py:class:`~genomicranges.GenomicRanges.GenomicRanges` or a
                :py:class:`~summarizedexperiment.RangedSummarizedExperiment.RangedSummarizedExperiment`
                object.

            query_type:
                Overlap query type, must be one of

                - "any": Any overlap is good
                - "start": Overlap at the beginning of the intervals
                - "end": Must overlap at the end of the intervals
                - "within": Fully contain the query interval

                Defaults to "any".

            max_gap:
                Maximum gap allowed in the overlap.
                Defaults to -1 (no gap allowed).

            min_overlap:
                Minimum overlap with query. Defaults to 1.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Raises:
            TypeError:
                If query is not a `RangedSummarizedExperiment` or `GenomicRanges`.

        Returns:
            A new `RangedSummarizedExperiment`
            object. None if there are no indices to slice.
        """
        _check_gr_or_rse(query)

        qranges = _access_granges(query)

        result = self.row_ranges.find_overlaps(
            query=qranges,
            query_type=query_type,
            max_gap=max_gap,
            min_overlap=min_overlap,
            ignore_strand=ignore_strand,
        )

        return self[result.get_column("self_hits"), :]

    def order(self, decreasing: bool = False) -> np.ndarray:
        """Get the order of indices to sort.

        Args:
            decreasing:
                Whether to sort in descending order. Defaults to False.

        Returns:
            NumPy vector containing index positions in the sorted order.
        """
        return self.row_ranges.order(decreasing=decreasing)

    def sort(self, decreasing: bool = False, in_place: bool = False) -> RangedSummarizedExperiment:
        """Sort by ranges.

        Args:
            decreasing:
                Whether to sort in descending order. Defaults to False.

            in_place:
                Whether to modify the object in place. Defaults to False.

        Returns:
            A new sorted `RangedSummarizedExperiment` object.
        """
        _order = self.row_ranges.order(decreasing=decreasing)

        output = self._define_output(in_place=in_place)
        return output[list(_order),]

    ############################
    ######>> combine ops <<#####
    ############################

    def relaxed_combine_rows(self, *other) -> RangedSummarizedExperiment:
        """Wrapper around :py:func:`~relaxed_combine_rows`."""
        return relaxed_combine_rows(self, *other)

    def relaxed_combine_columns(self, *other) -> RangedSummarizedExperiment:
        """Wrapper around :py:func:`~relaxed_combine_columns`."""
        return relaxed_combine_columns(self, *other)

    def combine_rows(self, *other) -> RangedSummarizedExperiment:
        """Wrapper around :py:func:`~combine_rows`."""
        return combine_rows(self, *other)

    def combine_columns(self, *other) -> RangedSummarizedExperiment:
        """Wrapper around :py:func:`~combine_columns`."""
        return combine_columns(self, *other)

    ######################
    ######>> to se <<#####
    ######################

    @classmethod
    def from_summarizedexperiment(cls, se: SummarizedExperiment) -> RangedSummarizedExperiment:
        return cls(
            assays=se.get_assays(),
            row_data=se.get_row_data(),
            column_data=se.get_column_data(),
            row_names=se.get_row_names(),
            column_names=se.get_column_names(),
            metadata=se.get_metadata(),
        )

    @classmethod
    def from_se(cls, se: SummarizedExperiment) -> RangedSummarizedExperiment:
        """Alias for :py:meth:`~to_summarizedexperiment`."""
        return RangedSummarizedExperiment.from_summarizedexperiment(se)

    def to_summarizedexperiment(self):
        """Coerce to :py:class:`summarizedexperiment.SummarizedExperiment.SummarizedExperiment`"""

        warn("SE does not containg a slot for ranges, merging ranges with row_data.", UserWarning)

        new_row_data = self._rows
        if self._row_ranges is not None:
            try:
                import copy

                if isinstance(self._row_ranges, CompressedGenomicRangesList):
                    gr = self._row_ranges.get_unlist_data()
                else:
                    gr = self._row_ranges

                range_dict = {
                    "seqnames": gr.get_seqnames(),
                    "starts": gr.get_start(),
                    "ends": gr.get_end(),
                    "strand": gr.get_strand(as_type="list"),
                }

                if gr.mcols is not None:
                    for col in gr.mcols.column_names:
                        range_dict[col] = gr.mcols.column(col)

                range_bf = biocframe.BiocFrame(range_dict, number_of_rows=len(gr), row_names=gr.get_names())

                if new_row_data is None or new_row_data.shape[0] == 0 or len(new_row_data.column_names) == 0:
                    new_row_data = range_bf
                else:
                    new_row_data = copy.copy(new_row_data)
                    for col in range_bf.column_names:
                        new_row_data.set_column(col, range_bf.column(col), in_place=True)
            except Exception as e:
                warn(f"Failed to transfer row_ranges to row_data: {e}", UserWarning)

        return SummarizedExperiment(
            assays=self._assays,
            row_data=new_row_data,
            column_data=self._cols,
            row_names=self._row_names,
            column_names=self._column_names,
            metadata=self._metadata,
            _validate=False,
        )

    def to_se(self):
        """Alias for :py:meth:`~to_summarizedexperiment`."""
        return self.to_summarizedexperiment()


############################
######>> combine ops <<#####
############################


@ut.combine_rows.register(RangedSummarizedExperiment)
def combine_rows(*x: RangedSummarizedExperiment) -> RangedSummarizedExperiment:
    """Combine multiple ``RangedSummarizedExperiment`` objects by row.

    All assays must contain the same assay names. If you need a
    flexible combine operation, checkout :py:func:`~relaxed_combine_rows`.

    Returns:
        A combined ``RangedSummarizedExperiment``.
    """
    first = x[0]
    _all_assays = [y.assays for y in x]
    check_assays_are_equal(_all_assays)
    _new_assays = merge_assays(_all_assays, by="row")

    _all_rows = [y._rows for y in x]
    _new_rows = ut.combine_rows(*_all_rows)
    _new_row_names = merge_se_rownames(x)

    _all_row_ranges = [y._row_ranges for y in x]
    _new_row_ranges = ut.combine_sequences(*_all_row_ranges)

    current_class_const = type(first)
    return current_class_const(
        assays=_new_assays,
        row_ranges=_new_row_ranges,
        row_data=_new_rows,
        column_data=first._cols,
        row_names=_new_row_names,
        column_names=first._column_names,
        metadata=first._metadata,
    )


@ut.combine_columns.register(RangedSummarizedExperiment)
def combine_columns(*x: RangedSummarizedExperiment) -> RangedSummarizedExperiment:
    """Combine multiple ``RangedSummarizedExperiment`` objects by column.

    All assays must contain the same assay names. If you need a
    flexible combine operation, checkout :py:func:`~relaxed_combine_columns`.

    Returns:
        A combined ``RangedSummarizedExperiment``.
    """
    first = x[0]
    _all_assays = [y.assays for y in x]
    check_assays_are_equal(_all_assays)
    _new_assays = merge_assays(_all_assays, by="column")

    _all_cols = [y._cols for y in x]
    _new_cols = ut.combine_rows(*_all_cols)
    _new_col_names = merge_se_colnames(x)

    current_class_const = type(first)
    return current_class_const(
        assays=_new_assays,
        row_ranges=first._row_ranges,
        row_data=first._rows,
        column_data=_new_cols,
        row_names=first._row_names,
        column_names=_new_col_names,
        metadata=first._metadata,
    )


@ut.relaxed_combine_rows.register(RangedSummarizedExperiment)
def relaxed_combine_rows(*x: RangedSummarizedExperiment) -> RangedSummarizedExperiment:
    """A relaxed version of the :py:func:`~biocutils.combine_rows.combine_rows` method for
    :py:class:`~RangedSummarizedExperiment` objects.  Whereas ``combine_rows`` expects that all objects have the same
    columns, ``relaxed_combine_rows`` allows for different columns. Absent columns in any object are filled in with
    appropriate placeholder values before combining.

    Args:
        x:
            One or more ``RangedSummarizedExperiment`` objects, possibly with differences in the
            number and identity of their columns.

    Returns:
        A ``RangedSummarizedExperiment`` that combines all ``experiments`` along their rows and contains
        the union of all columns. Columns absent in any ``x`` are filled in
        with placeholders consisting of Nones or masked NumPy values.
    """
    first = x[0]
    _new_assays = relaxed_merge_assays(x, by="row")

    _all_rows = [y._rows for y in x]
    _new_rows = biocframe.relaxed_combine_rows(*_all_rows)
    _new_row_names = merge_se_rownames(x)

    _all_row_ranges = [y._row_ranges for y in x]
    _new_row_ranges = ut.combine_sequences(*_all_row_ranges)

    current_class_const = type(first)
    return current_class_const(
        assays=_new_assays,
        row_ranges=_new_row_ranges,
        row_data=_new_rows,
        column_data=first._cols,
        row_names=_new_row_names,
        column_names=first._column_names,
        metadata=first._metadata,
    )


@ut.relaxed_combine_columns.register(RangedSummarizedExperiment)
def relaxed_combine_columns(
    *x: RangedSummarizedExperiment,
) -> RangedSummarizedExperiment:
    """A relaxed version of the :py:func:`~biocutils.combine_rows.combine_columns` method for
    :py:class:`~RangedSummarizedExperiment` objects.  Whereas ``combine_columns`` expects that all objects have the same
    rows, ``relaxed_combine_columns`` allows for different rows. Absent columns in any object are filled in with
    appropriate placeholder values before combining.

    Args:
        x:
            One or more ``RangedSummarizedExperiment`` objects, possibly with differences in the
            number and identity of their rows.

    Returns:
        A ``RangedSummarizedExperiment`` that combines all ``experiments`` along their columns and contains
        the union of all rows. Rows absent in any ``x`` are filled in
        with placeholders consisting of Nones or masked NumPy values.
    """
    first = x[0]
    _new_assays = relaxed_merge_assays(x, by="column")

    _all_cols = [y._cols for y in x]
    _new_cols = biocframe.relaxed_combine_rows(*_all_cols)
    _new_col_names = merge_se_colnames(x)

    current_class_const = type(first)
    return current_class_const(
        assays=_new_assays,
        row_ranges=first._row_ranges,
        row_data=first._rows,
        column_data=_new_cols,
        row_names=first._row_names,
        column_names=_new_col_names,
        metadata=first._metadata,
    )


@ut.extract_row_names.register(RangedSummarizedExperiment)
def _rownames_rse(x: RangedSummarizedExperiment):
    return x.get_row_names()


@ut.extract_column_names.register(RangedSummarizedExperiment)
def _colnames_rse(x: RangedSummarizedExperiment):
    return x.get_column_names()
