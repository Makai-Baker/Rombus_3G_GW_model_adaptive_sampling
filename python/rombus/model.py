import importlib
import sys
import os
import shutil
import timeit
from abc import ABCMeta, abstractmethod
from collections import Counter
from typing import Any, Dict, Self, Optional, Tuple, TYPE_CHECKING

import numpy as np

import rombus._core.hdf5 as hdf5
import rombus.exceptions as exceptions
from rombus.params import Params
from rombus._core.log import log
from typing import NamedTuple

import warnings

warnings.simplefilter("ignore", np.ComplexWarning)

# Need to put Samples in quotes below and check TYPE_CHECKING here to
# manage circular imports with models.py
if TYPE_CHECKING:
    from rombus.samples import Samples


class _Ordinate(object):
    def __init__(self):
        self.name = None
        self.dtype = float
        self.label = None

    def set(
        self, name: Optional[str], dtype: type | np.dtype = float, label: str = ""
    ) -> None:
        """Set the details of the model's ordinate.

        Parameters
        ----------
        name : str
            Name of ordinate
        dtype : type|np.dtype
            Datatype used to represent ordinate values
        label : str
            A n optional label to use for the ordinate in plots, etc.
        """
        self.name = name
        self.dtype = dtype
        if label == "":
            self.label = self.name
        else:
            self.label = label


class _Coordinate(object):
    def __init__(self):
        self.name = None
        self.dtype = float
        self.min = None
        self.max = None
        self.label = None
        self._values = None
        self._n_values = 0

    def set(
        self,
        name: str,
        min: Any,
        max: Any,
        n_values: int,
        dtype: type | np.dtype = float,
        label: str = "",
        values: arr = []
    ) -> None:
        self.name = name
        self.dtype = dtype
        self.min = min
        self.max = max
        if label == "":
            self.label = self.name
        else:
            self.label = label
        if self.dtype != type(min):
            raise exceptions.RombusModelCoordinateError(
                f"Coordinate datatype ({self.dtype}) and min-value datatype ({type(self.min)}) don't match."
            )
        if self.dtype != type(max):
            raise exceptions.RombusModelCoordinateError(
                f"Coordinate datatype ({self.dtype}) and max-value datatype ({type(self.max)}) don't match."
            )
        self._n_values = n_values
        if values == []:
            self._values = np.linspace(self.min, self.max, self._n_values, self.dtype)
        else:
            self._values = values
            
    def get(self):
        return self._values

    def __len__(self):
        return self._n_values


class _RombusModelMeta(type):
    def __prepare__(name, *args, **kwargs):
        """Initialise the dictionary that gets passed to __new___.

        This is needed here because we don't want the user to have to
        initialise the member(s) that we are adding.  This is the only
        method that gets sourced before the class code is executed, so
        it needs to be done here, not in __new__.
        """

        result = dict()

        # Initialise the following members↵
        result["coordinate"] = _Coordinate()
        result["ordinate"] = _Ordinate()
        result["params"] = Params()

        return result

    def __new__(mcs, name, bases, dct):

        # Perform super-metaclass construction↵
        return super(_RombusModelMeta, mcs).__new__(mcs, name, bases, dct)


class _RombusModelABCMeta(_RombusModelMeta, ABCMeta):
    """Turn _RombusModelMeta into an abstract base class"""

    pass


class RombusModel(metaclass=_RombusModelABCMeta):
    """Baseclass from which all RombusModels must inherit."""

    # These members are instantiated by the metaclass
    coordinate: _Coordinate
    """The domain on which this model is defined."""
    ordinate: _Ordinate
    """The ordinate that this model maps its domain to."""
    params: Params
    """The parameters defined for this model"""

    def __init__(self, model: str):

        # Keep track of the model string so we can reinstantiate from a saved state
        self.model_str = model
        self.basename = self.model_str.split(":")[0].split(".")[-1]

        # Initialise the domain
        self.domain = self.coordinate.get()
        self.n_domain = len(self.coordinate)

        # Check that the domain has been suitably set
        assert self.n_domain > 0

        # Check that at least one parameter has beed defined
        if self.params.count <= 0:
            raise exceptions.RombusModelParamsError(
                f"Invalid number of parameters ({self.params.count}) specified for Rombus model ({self})."
            )

    def __str__(self):
        return f"<RombusModel from {self.model_str}>"

    @abstractmethod  # make sure this is the inner-most decorator
    def compute(self, params: NamedTuple, domain: np.ndarray) -> np.ndarray:
        """Abstract method which computes the user's model.

        This method does all the work of computing the user's model.  It takes a parameter set as a named tuple
        with N elements given by the names given to the N calls made to params.add() as well as the array set
        by coordinate.set() and returns a numpy array.

        Parameters
        ----------
        params : NamedTuple
            The parameters to be used when computing the model
        domain : np.ndarray
            The domain on which the model is to be computed

        Returns
        -------
        np.ndarray
            The user's model, computed for the given parameter set and domain
        """
        pass

    @classmethod
    @log.callable("Instantiating model from file")
    def from_file(cls, file_in: hdf5.FileOrFilename) -> Self:
        """Generate a RombusModel instance from a Rombus HDF5 file.

        Parameters
        ----------
        file_in : hdf5.FileOrFilename
            The Rombus HDF5 file to read from.

        Returns
        -------
        Self
            The generated RombusModel instance.
        """

        try:
            h5file, close_file = hdf5.ensure_open(file_in)
            model_str = h5file["model/model_str"].asstr()[()]
            if close_file:
                h5file.close()
        except IOError as e:
            log.handle_exception(e)
        return cls.load(model_str)

    @classmethod
    @log.callable("Loading model from file")
    def load(cls, model: str | Self) -> Self:

        """Ensure that a model has been imported for use by Rombus.

        Parameters
        ----------
        model : str | Self
            A string of format 'sub,module.name:ClassName' or a RombusModel instance
            (trivially returned in the later case).

        Returns
        -------
        Self
            A RombusModel instance
        """
        if isinstance(model, str):
            try:
                model_class = _import_from_string(model)
            except exceptions.RombusException as e:
                log.handle_exception(e)
            else:
                return model_class(model)
        elif not isinstance(model, RombusModel):
            raise exceptions.RombusModelInitError(
                "Invalid type ({type(model)}) specified when loading model {model}."
            )
        return model  # type: ignore

    @log.callable("Writing model to file")
    def write(self, h5file: hdf5.File) -> None:
        """Write a RombusModel to a Rombus HDF5 file.

        Parameters
        ----------
        h5file : hdf5.File
            An open HDF5 file
        """

        try:
            h5_group = h5file.create_group("model")
            h5_group.create_dataset("model_str", data=self.model_str)
        except IOError as e:
            log.handle_exception(e)

    # Need to put Samples in quotes and check TYPE_CHECKING above to manage circular import with models.py
    def generate_model_set(self, samples: "Samples") -> np.ndarray:
        """Generate a set of models for a given set of parameter samples.

        Parameters
        ----------
        samples : "Samples"
            A set of Samples

        Returns
        -------
        np.ndarray
            An array of model results: 1 for each given sample.
        """

        my_ts: np.ndarray = np.zeros(
            shape=(samples.n_samples, self.n_domain), dtype=self.ordinate.dtype
        )
        with log.progress("Generating training set", samples.n_samples) as progress:
            for i, params_numpy in enumerate(samples.samples):
                model_i = self.compute(self.params.np2param(params_numpy), self.domain)
                my_ts[i] = model_i / np.sqrt(np.vdot(model_i, model_i).real)
                progress.update(i)

        return my_ts

    def parse_cli_params(self, args: Tuple[str, ...]) -> Dict[str, Any]:
        """Parse parameters given as a tuple of form 'param0=val0', 'param1=val1', ... etc to a dictionary of
        form {'param0':val0, 'param1':val1, ... }.

        Generally used to parse the optional arguments recieved from Click
        into a format that can be converted into a Params or Numpy object

        Parameters
        ----------
        args : Tuple[str, ...]
            Given tuple of parameters

        Returns
        -------
        Dict[str, Any]
            Resulting dict of parameters
        """

        model_params = dict()
        for param_i in args:
            if not param_i.startswith("-"):
                res = param_i.split("=")
                if len(res) == 2:
                    # NOTE: for now, all parameters are assumed to be floats
                    model_params[res[0]] = float(res[1])
                else:
                    raise Exception(f"Don't know what to do with argument '{param_i}'")
            else:
                raise Exception(f"Don't know what to do with option '{param_i}'")

        # Check that all parameters are specified and that they match what is
        # defined in the model
        assert Counter(model_params.keys()) == Counter(self.params.names)

        return model_params

    def sample(self, kwargs: Dict[str, Any]) -> NamedTuple:
        """Create a Sample from a dictionary of the form {"param0:val0,"param1":val1,...}

        Parameters
        ----------
        kwargs : Dict[str, Any]
            Dict specifying parameter values

        Returns
        -------
        NamedTuple
            Named tuple specifying paramter values
        """
        return self.params.params_dtype(**kwargs)  # type: ignore

    def timing(self, samples: "Samples") -> float:
        """Generate timing information for the original source model.  Particularly useful when compared to
        similar timing information computed for ROMs derived from it.

        Parameters
        ----------
        samples : "Samples"
            A set of parameters to generate timing information for.  Should be the same as those used when
            timiing a ROM, if comparisons are to be made.

        Returns
        -------
        float
            Seconds elapsed
        """

        with log.context(
            f"Computing timing information for model using {samples.n_samples} samples",
            time_elapsed=False,
        ):
            start_time = timeit.default_timer()
            for i, sample in enumerate(samples.samples):
                params_numpy = self.params.np2param(sample)
                _ = self.compute(params_numpy, self.domain)
        return timeit.default_timer() - start_time

    @classmethod
    @log.callable("Writing project template")
    def write_project_template(cls, project_name: str) -> None:
        """Write a project model to the current working directory to start a new project from.

        Two files are written to the current working direcory: a Python file and a set of samples.  These can
        then be modified to suit the needs of the user.

        Parameters
        ----------
        project_name : str
            Base name to use for the project.
        """

        # Set the model we will template from
        model_name = "sinc"

        # Set source file paths
        pkgdir = sys.modules["rombus"].__path__[0]
        model_file_source = os.path.join(pkgdir, "models", f"{model_name}.py")
        samples_file_source = os.path.join(
            pkgdir, "models", f"{model_name}_samples.csv"
        )

        # Set output file paths
        model_file_out = os.path.join(os.getcwd(), f"{project_name}.py")
        samples_file_out = os.path.join(os.getcwd(), f"{project_name}_samples.csv")

        # Copy files
        with log.context("Writing files"):
            shutil.copy(model_file_source, model_file_out)
            log.comment(f"Written: {os.path.split(model_file_out)[1]}")
            shutil.copy(samples_file_source, samples_file_out)
            log.comment(f"Written: {os.path.split(samples_file_out)[1]}")


RombusModelType = RombusModel | str

# The code that follows is modified from code copied from the Uvicorn codebase:
#     https://github.com/encode/uvicorn (commit: d613cbea388bafafb6f642077c035ed137deea61)
# Copyright © 2017-present, [Encode OSS Ltd](https://www.encode.io/).
# All rights reserved.
@log.callable("Importing model")
def _import_from_string(import_str: str) -> Any:

    """Import a RombusModel class from a given string of the form 'python.module.name:ClassName'.

    Generally, the user model will be defined in a file in the current working directory with filename (for
    example) of 'model_name.py', with a class inheriting from RombusModel with name 'ClassName'.  It should
    then be referred to in this context as 'my_model:ClassName'.  More generally, the model can be anywhere
    in the user's PYTHONPATH.

    Parameters
    ----------
    import_str : str
        Given string of the form 'python.module.name:ClassName'

    Returns
    -------
    Any
        An instance of the user-defined model class
    """

    log.append(f"({import_str})...")

    if not isinstance(import_str, str):
        raise exceptions.RombusModelImportFromStringError(
            f'Import string must be a string with format "<module>:<attribute>".  It is actually of type {type(import_str)}.'
        )

    # Make sure the CWD is in the import path
    sys.path.append(os.getcwd())
    sys.path = list(dict.fromkeys(sys.path))

    # Split the string
    module_str, _, attrs_str = import_str.partition(":")
    if not module_str or not attrs_str:
        raise exceptions.RombusModelImportFromStringError(
            f'Import string "{import_str}" must be in format "<module>:<attribute>".'
        )

    # Try to import the module
    try:
        module = importlib.import_module(module_str)
    except ImportError as exc:
        if exc.name != module_str:
            raise exc from None
        raise exceptions.RombusModelImportFromStringError(
            f'Could not import module "{module_str}".\n'
        )
    instance = module

    # Try to grab the specified class
    try:
        for attr_str in attrs_str.split("."):
            instance = getattr(instance, attr_str)
    except AttributeError:
        raise exceptions.RombusModelImportFromStringError(
            f'Attribute "{attrs_str}" not found in module "{module_str}".'
        )

    return instance
