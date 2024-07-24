from rombus._core.log import log


class RombusException(Exception):
    """Base class for exceptions raised from Rombus"""

    def __init__(self, message: str):
        self._message = message
        self._is_rombus_exception = True

    def __str__(self) -> str:
        return f"{self._message}"

    def handle_exception(self) -> None:
        log.handle_exception(self)


# rombus.ei exceptions


class EmpiricalInterpolantNotComputedError(RombusException):
    """Raised when an EmpiricalInterpolant operation is attempted on a ROM whose EI has not yet been computed."""

    pass


# rombus.model exceptions


class RombusModelOrdinateError(RombusException):
    """Raised when an error is encountered when initialising a RombusModel ordinate."""

    pass


class RombusModelCoordinateError(RombusException):
    """Raised when an error is encountered when initialising a RombusModel coordinate."""

    pass


class RombusModelParamsError(RombusException):
    """Raised when a Rombus model is instantiated with no parameters specified."""

    pass


class RombusModelImportFromStringError(RombusException):
    """Raised when the import of a RombusModel from a string fails."""

    pass


class RombusModelInitError(RombusException):
    """Raised when the instantiation of a RombusModel fails."""

    pass


# rombus.plots exceptions


class RombusPlotError(RombusException):
    """Raised when there is an error generating a plot."""

    pass


# rombus.samples exceptions


# class ReducedSamplesError(RombusException):
#    """Raised when ."""
#
#    pass


# rombus.reduced_basis exceptions


class ReducedBasisInitError(RombusException):
    """Raised when a ReducedBasis object can not be initialised."""

    pass


class ReducedBasisNotComputedError(RombusException):
    """Raised when an EmpiricalInterpolant operation is attempted on a ROM whose EI has not yet been computed."""

    pass


class ReducedBasisComputeError(RombusException):
    """Raised when a ReducedBasis object can not be computed."""

    pass


class RombusModelLoadError(RombusException):
    """Raised when a Rombus model fails to load."""

    pass


# rombus.rom exceptions


class RomNotInitialised(RombusException):
    """Raised when an uninitialised ROM is accessed."""

    pass


# rombus._core.hdf5 exceptions


class RombusHDF5Error(RombusException):
    """Raised when an HDF5 error is raised when accessing a ROM file."""

    pass
