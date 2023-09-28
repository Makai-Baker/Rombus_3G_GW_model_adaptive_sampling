from typing import Callable, List, NamedTuple, Optional

import numpy as np
import rombus.exceptions as exceptions
from rombus._core.log import log

MAX_N_TRIES = 1000

speed_of_light = 299792458.0  # m/s
solar_mass = 1.988409870698050731911960804878414216e30  # Kg
gravitational_constant = 6.6743e-11  # m^3 kg^-1 s^-2

def freq_to_time(chirp_mass, frequency):
    return (5/256)*(chirp_mass*solar_mass)**(-5/3) * \ 
    (np.pi*frequency)**(-8/3)*speed_of_light**5 * gravitational_constant**(-5/3)

class Params(object):
    def __init__(self) -> None:
        self.params: List[NamedTuple] = []
        self.names: List[str] = []
        self.count: int = 0
        self._idx: int = -1
        self.params_dtype: Optional[NamedTuple] = None
        self._validate: Callable = lambda param: True

    def add(
        self, name: str, range_min: float, range_max: float, dtype: type = float
    ) -> None:
        """Add a parameter to this parameter set.

        Parameters
        ----------
        name : str
            Name of the parameter.  This will become the reference name when accessing NamedTuple & dict versions of the parameters
        range_min : float
            The minimum-allowed value for this parameter
        range_max : float
            The maximum-allowed value for this parameter
        dtype : type
            The datatype used to represent this parameter
        """

        with log.context(f"Adding parameter ({name})"):
            # N.B.: mypy struggles with NamedTuples, so typing is turned off for the following
            Param = NamedTuple(f"Param{len(self.params)}", [("name", str), ("min", dtype), ("max", dtype), ("dtype", dtype)])  # type: ignore
            self.params.append(Param(name, dtype(range_min), dtype(range_max), dtype))  # type: ignore
            self.count = self.count + 1
            self.names.append(name)

            # Update the named tuple that will be used to convert numpy arrays to Param objects
            # N.B.: mypy struggles with NamedTuples, so typing is turned off for the following
            self.params_dtype = NamedTuple("params_dtype", [(name, dtype) for name in self.names])  # type: ignore

    def set_validation(self, func: Callable[[NamedTuple], bool]) -> None:
        """Specify a callable which accepts a parameter set (in the form of a NamedTuple) and
        returns a bool indicating if the parameter set is valid.

        By default, no validation is performed and all passed samples are assumed to be valid.

        Parameters
        ----------
        func : Callable
            Callable which will perform the validation
        """
        self._validate = func

    def generate_random_sample(
        self, random_generator: np.random._generator.Generator
    ) -> np.ndarray:

        """Generate a random parameter set, with values uniformly distributed between their specified min and max values.

        Parameters
        ----------
        random_generator : np.random._generator.Generator
            An initialised numpy random number generator

        Returns
        -------
        np.ndarray
            A parameter set in numpy form
        """
        new_sample: np.ndarray = np.ndarray(self.count, dtype=np.float64)
        n_tries = 0
        try:
            while True:
                for i, param in enumerate(self.params):
                    # N.B.: mypy struggles with NamedTuples, so typing is turned off for this next line
                    new_sample[i] = random_generator.uniform(low=param.min, high=param.max)  # type: ignore
                param = self.np2param(new_sample)
                if self._validate(param):
                    break
                else:
                    n_tries = n_tries + 1
                    if n_tries >= MAX_N_TRIES:
                        raise exceptions.RombusModelParamsError(
                            f"Max number of tries ({MAX_N_TRIES}) reached when trying to generate a valid random parameter set"
                        )
        except exceptions.RombusException as e:
            log.handle_exception(e)
        return new_sample

    def np2param(self, params_np: np.ndarray) -> NamedTuple:
        """Convert a numpy version of a parameter set to a NamedTuple format.

        Parameters
        ----------
        params_np : np.ndarray
            Numpy version of parameter set

        Returns
        -------
        NamedTuple
            NamedTuple version of parameter set
        """
        if not self.params_dtype:
            return NamedTuple("params_empty", [])
        else:
            # N.B.: mypy struggles with NamedTuples, so typing is turned off for this next line
            return self.params_dtype(**dict(zip(self.names, np.atleast_1d(params_np))))  # type: ignore

    def __iter__(self):
        self._idx = -1
        return self

    def __next__(self):
        self._idx = self._idx + 1
        if self._idx >= self.count:
            raise StopIteration
        return self.params[self._idx]

    def __len__(self):
        return self.count
