# Extending summarized experiments

`SummarizedExperiment` and `RangedSummarizedExperiment` serve as the foundational container classes for representing experimental data and its metadata. Developers can extend these classes to create specialized data structures tailored for specific scientific applications. In fact, our implementation of `SingleCellExperiment` extends the `RangedSummarizedExperiment` class.

## Define the new class

As a simple example, let's create a new class called `BioSampleSE` that stores biosample information on which the experiment was conducted. This may contain anonymized information about the patient(s) or sample(s). First, we extend the `SummarizedExperiment` class:


```python
from summarizedexperiment import SummarizedExperiment


class BioSampleSE(SummarizedExperiment):
    pass
```

## Add a new slot or attribute

To add a new slot to this class, we accept a new parameter `bio_sample_information` when the class is initialized through the __init__ method. We also provide type hints to set expectations on the accepted types for these arguments. Type hints are helpful for users and are automatically annotated in the documentation. More information on type hints can be found in our [developer guide](https://github.com/BiocPy/developer_guide).

We forward the default parameters to the base `SummarizedExperiment` class (using the `super` method) and store the new attribute.

```python
from summarizedexperiment import SummarizedExperiment
import biocframe
from typing import Dict, Any, List, Optional


class BioSampleSE(SummarizedExperiment):
    def __init__(
        self,
        assays: Dict[str, Any] = None,
        row_data: Optional[biocframe.BiocFrame] = None,
        column_data: Optional[biocframe.BiocFrame] = None,
        row_names: Optional[List[str]] = None,
        column_names: Optional[List[str]] = None,
        metadata: Optional[dict] = None,
        bio_sample_information: Optional[biocframe.BiocFrame] = None,  # NEW SLOT
        validate: bool = True,
    ) -> None:
        super().__init__(
            assays,
            row_data=row_data,
            column_data=column_data,
            row_names=row_names,
            column_names=column_names,
            metadata=metadata,
            validate=validate,
        )

        self._bio_sample_information = bio_sample_information
```

The new slot can be validated using a dedicated validator:

```python
def _validate_bio_sample_information(bio_sample_info):
    if not isinstance(bio_sample_info, biocframe.BiocFrame):
        raise Exception("Biosample information must be a BiocFrame object.")

    # any other validations. for example, if you have expectations on the columns
    # of this frame or the number of rows.
```

Our class now validates the new slot:

```python
class BioSampleSE(SummarizedExperiment):
    def __init__(
        self,
        assays: Dict[str, Any] = None,
        row_data: Optional[biocframe.BiocFrame] = None,
        column_data: Optional[biocframe.BiocFrame] = None,
        row_names: Optional[List[str]] = None,
        column_names: Optional[List[str]] = None,
        metadata: Optional[dict] = None,
        bio_sample_information: Optional[biocframe.BiocFrame] = None,  # NEW SLOT
        validate: bool = True,
    ) -> None:
        super().__init__(
            assays,
            row_data=row_data,
            column_data=column_data,
            row_names=row_names,
            column_names=column_names,
            metadata=metadata,
            validate=validate,
        )

        self._bio_sample_information = bio_sample_information

        if validate:
            _validate_bio_sample_information(self._bio_sample_information)
```


## Define getters/setters

We need accessors and setters so users can interact with this new property. We provide these accessors in both functional and property-based approaches. More details can be found in our [class design](https://biocpy.github.io/tutorial/chapters/philosophy.html#functional-discipline) document.

```python
def get_bio_sample_information(self) -> Optional[biocframe.BiocFrame]:
    """Get biosample information.

    Returns:
        biosample information or None if not availabl.
    """
    return self._bio_sample_information


def set_bio_sample_information(
    self, bio_sample_information: Optional[biocframe.BiocFrame], in_place: bool = False
) -> "BioSampleSE":
    """Set new biosample information.

    Args:
        bio_sample_information:
            A new `BiocFrame` object containing biosample information.

        in_place:
            Whether to modify the ``BioSampleSE`` in place.

    Returns:
        A modified ``BioSampleSE`` object, either as a copy of the original
        or as a reference to the (in-place-modified) original.
    """
    _validate_bio_sample_information(bio_sample_info)

    output = self._define_output(in_place)  # MAKES A SHALLOW COPY
    output._bio_sample_information = bio_sample_info
    return output
```

Additionally, let's provide property-based accessors to the new attribute:

```python
@property
def bio_sample_information(self) -> biocframe.BiocFrame:
    """Alias for :py:meth:`~get_bio_sample_info`."""
    return self.get_bio_sample_info()


@bio_sample_information.setter
def bio_sample_information(self, bio_sample_info: biocframe.BiocFrame) -> None:
    """Alias for :py:meth:`~set_bio_sample_info`."""
    warn(
        "Setting property 'bio_sample_information' is an in-place operation, use 'set_bio_sample_info' instead",
        UserWarning,
    )
    return self.set_bio_sample_information(row_ranges=row_ranges, in_place=True)
```

This allows users to easily access the new property using the **dot** notation on an instance, for example, `obj.bio_sample_info` provides access to the attribute.

## Define shallow and deep copy methods

To avoid mutating objects in-place, methods for making shallow and deep copies of the class attributes are implemented.

```python
def __deepcopy__(self, memo=None, _nil=[]):
    """
    Returns:
        A deep copy of the current ``BioSampleSE``.
    """
    from copy import deepcopy

    _assays_copy = deepcopy(self._assays)
    _rows_copy = deepcopy(self._rows)
    _cols_copy = deepcopy(self._cols)
    _row_names_copy = deepcopy(self._row_names)
    _col_names_copy = deepcopy(self._column_names)
    _bio_sample_information_copy = deepcopy(self._bio_sample_information)
    _metadata_copy = deepcopy(self.metadata)

    current_class_const = type(self)
    return current_class_const(
        assays=_assays_copy,
        row_data=_rows_copy,
        column_data=_cols_copy,
        row_names=_row_names_copy,
        column_names=_col_names_copy,
        bio_sample_information=_bio_sample_information_copy,
        metadata=_metadata_copy,
    )


def __copy__(self):
    """
    Returns:
        A shallow copy of the current ``BioSampleSE``.
    """
    current_class_const = type(self)
    return current_class_const(
        assays=self._assays,
        row_data=self._rows,
        column_data=self._cols,
        row_names=self._row_names,
        column_names=self._column_names,
        bio_sample_information=self._bio_sample_information,
        metadata=self._metadata,
    )


def copy(self):
    """Alias for :py:meth:`~__copy__`."""
    return self.__copy__()
```

## Subset operation

When the experiment is subsetted using the `[]` operator, we may also have to subset the new slots added to the extended class. Developers only need to extend the `get_slice` function from the base class that accepts row and column subsets. The `_generic_slice` function returns a `slicer` object containing normalized row and column indices, and properties from the base class.

```python
from typing import Sequence

def get_slice(
    self,
    rows: Optional[Union[str, int, bool, Sequence]],
    columns: Optional[Union[str, int, bool, Sequence]],
) -> "BioSampleSE":
    """Alias for :py:attr:`~__getitem__`, for back-compatibility."""

    slicer = self._generic_slice(rows=rows, columns=columns)

    # An illustrative example
    new_bio_sample_info = None
    if slicer.row_indices != slice(None):
        new_bio_sample_info = self.bio_sample_info[slicer.row_indices]

    current_class_const = type(self)
    return current_class_const(
        assays=slicer.assays,
        row_data=slicer.rows,
        column_data=slicer.columns,
        row_names=slicer.row_names,
        column_names=slicer.column_names,
        bio_sample_info=self.bio_sample_info.
        metadata=self._metadata,
    )
```

## Putting it all together


```python
def _validate_bio_sample_information(bio_sample_info):
    if not isinstance(bio_sample_info, biocframe.BiocFrame):
        raise Exception("Biosample information must be a BiocFrame object.")


class BioSampleSE(SummarizedExperiment):
    def __init__(
        self,
        assays: Dict[str, Any] = None,
        row_data: Optional[biocframe.BiocFrame] = None,
        column_data: Optional[biocframe.BiocFrame] = None,
        row_names: Optional[List[str]] = None,
        column_names: Optional[List[str]] = None,
        metadata: Optional[dict] = None,
        bio_sample_information: Optional[biocframe.BiocFrame] = None,  # NEW SLOT
        validate: bool = True,
    ) -> None:
        super().__init__(
            assays,
            row_data=row_data,
            column_data=column_data,
            row_names=row_names,
            column_names=column_names,
            metadata=metadata,
            validate=validate,
        )

        self._bio_sample_information = bio_sample_information

        if validate:
            _validate_bio_sample_information(self._bio_sample_information)

    def __deepcopy__(self, memo=None, _nil=[]):
        """
        Returns:
            A deep copy of the current ``BioSampleSE``.
        """
        from copy import deepcopy

        _assays_copy = deepcopy(self._assays)
        _rows_copy = deepcopy(self._rows)
        _cols_copy = deepcopy(self._cols)
        _row_names_copy = deepcopy(self._row_names)
        _col_names_copy = deepcopy(self._column_names)
        _bio_sample_information_copy = deepcopy(self._bio_sample_information)
        _metadata_copy = deepcopy(self.metadata)

        current_class_const = type(self)
        return current_class_const(
            assays=_assays_copy,
            row_data=_rows_copy,
            column_data=_cols_copy,
            row_names=_row_names_copy,
            column_names=_col_names_copy,
            bio_sample_information=_bio_sample_information_copy,
            metadata=_metadata_copy,
        )

    def __copy__(self):
        """
        Returns:
            A shallow copy of the current ``BioSampleSE``.
        """
        current_class_const = type(self)
        return current_class_const(
            assays=self._assays,
            row_data=self._rows,
            column_data=self._cols,
            row_names=self._row_names,
            column_names=self._column_names,
            bio_sample_information=self._bio_sample_information,
            metadata=self._metadata,
        )

    def copy(self):
        """Alias for :py:meth:`~__copy__`."""
        return self.__copy__()

    def get_bio_sample_information(self) -> Optional[biocframe.BiocFrame]:
        """Get biosample information.

        Returns:
            biosample information or None if not availabl.
        """
        return self._bio_sample_information

    def set_bio_sample_information(
        self, bio_sample_information: Optional[biocframe.BiocFrame], in_place: bool = False
    ) -> "BioSampleSE":
        """Set new biosample information.

        Args:
            bio_sample_information:
                A new `BiocFrame` object containing biosample information.

            in_place:
                Whether to modify the ``BioSampleSE`` in place.

        Returns:
            A modified ``BioSampleSE`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        _validate_bio_sample_information(bio_sample_info)

        output = self._define_output(in_place)  # MAKES A SHALLOW COPY
        output._bio_sample_information = bio_sample_info
        return output

    @property
    def bio_sample_information(self) -> biocframe.BiocFrame:
        """Alias for :py:meth:`~get_bio_sample_info`."""
        return self.get_bio_sample_info()

    @bio_sample_information.setter
    def bio_sample_information(self, bio_sample_info: biocframe.BiocFrame) -> None:
        """Alias for :py:meth:`~set_bio_sample_info`."""
        warn(
            "Setting property 'bio_sample_information' is an in-place operation, use 'set_bio_sample_info' instead",
            UserWarning,
        )
        return self.set_bio_sample_information(row_ranges=row_ranges, in_place=True)
```

That's the minimum required to extend a `SummarizedExperiment` and adapt it to new use cases. Please follow the [developer guide](https://github.com/BiocPy/developer_guide), which provides information on class design, package setup, and documentation to ensure consistency in how BiocPy-related packages are developed.
