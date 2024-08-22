import os
import timeit


import h5py  # type: ignore
import mpi4py
import numpy as np

from typing import Optional, Self, NamedTuple

import rombus._core.mpi as mpi
import rombus._core.hdf5 as hdf5
import rombus.exceptions as exceptions
from rombus._core.log import log
from rombus.model import RombusModel, RombusModelType
from rombus.samples import Samples
from rombus.ei import EmpiricalInterpolant
from rombus.reduced_basis import ReducedBasis

DEFAULT_TOLERANCE: float = 1e-14
DEFAULT_REFINE_N_RANDOM: int = 100


def weighted_product(weights, a, b):
    if b.ndim > 1 or a.ndim > 1: #sometimes want dot(h, B matrix), not just dot(h^*, g).
        return np.dot(a / weights, b)
    else:
        return np.vdot(a / weights, b)

class ReducedOrderModel(object):
    """Class for managing the creation, updating and subsequent use of a Reduced Order Model (ROM)."""

    def __init__(
        self,
        model: RombusModelType,
        samples: Samples,
        reduced_basis: Optional[ReducedBasis] = None,
        empirical_interpolant: Optional[EmpiricalInterpolant] = None,
        basename: Optional[str] = None,
        tol: float = DEFAULT_TOLERANCE,
        psd = 1.,
        product = None
    ):

        self.model: RombusModel = RombusModel.load(model)
        """Model used to generate the ROM"""

        self.samples = samples
        """Samples fed to the greedy algorithm to generate the ROM"""

        self.reduced_basis = reduced_basis
        """ReducedBasis generated for the ROM"""

        self.empirical_interpolant = empirical_interpolant
        """EmpiricalInterpolant generated for the ROM"""

        if basename is None:
            basename = self.model.basename
        self.basename = basename
        """Set when reading from files and provides a base name for writing plots to file, etc."""

        self.tol = tol 
        """User defined error tolerance for the RB greedy algorithm"""

        self.psd = psd 
        """Power spectral density used for noise-weighted inner products. Default is 1., i.e. unweighted inner product."""

        if product is None:
            product = weighted_product
        self.product = product
        """Defines the inner product used to build the Reduced Basis, Training set, and Empirical Interpolant."""

    @classmethod
    @log.callable("Instantiating ROM from file")
    def from_file(cls, file_in: hdf5.FileOrFilename) -> Self:
        """Instantiate a ROM from a Rombus HDF5 file.

        Parameters
        ----------
        file_in : hdf5.FileOrFilename
            Rombus file (filename or opened file) to read from

        Returns
        -------
        Self
            Return a reference to self so that methods can be chained.
        """

        h5file, close_file = hdf5.ensure_open(file_in)

        model: RombusModel = RombusModel.from_file(h5file)
        samples: Samples = Samples.from_file(h5file)
        reduced_basis: ReducedBasis = ReducedBasis.from_file(h5file)
        empirical_interpolant: EmpiricalInterpolant = EmpiricalInterpolant.from_file(
            h5file
        )

        basename = os.path.splitext(os.path.basename(h5file.filename))[0]

        if close_file:
            h5file.close()

        return cls(
            model,
            samples,
            reduced_basis=reduced_basis,
            empirical_interpolant=empirical_interpolant,
            basename=basename,
        )

    @log.callable("Building ROM")
    def build(
        self, do_step: Optional[str] = None, tol: float = DEFAULT_TOLERANCE
    ) -> Self:
        """(Re)build a ReducedOrderModel.

        Parameters
        ----------
        do_step : str|None
            Specify whether to just compute the ReducedBasis ('RB') or the EmpiricalInterpolant ('EI') or both (None)
        tol : float
            Absolute error tolerance when building the reduced basis

        Returns
        -------
        Self
            Returns a reference to self, so that method calls can be chained
        """

        if do_step is None or do_step == "RB":
            try:
                self.reduced_basis = ReducedBasis().compute(
                    self.model, self.samples, tol=tol, psd=self.psd, product=self.product
                )
            except exceptions.RombusException as e:
                e.handle_exception()
        
        if do_step is None or do_step == "EI":
            if self.reduced_basis is None:
                raise exceptions.ReducedBasisNotComputedError(
                    "A ROM whose ReducedBasis has not been computed has been asked to comput its EmpiricalInterpolant.  Compute the ReducedBasis first and try again."
                )
            self.empirical_interpolant = EmpiricalInterpolant().compute(
                self.reduced_basis
            )

        return self

    def evaluate(self, params: NamedTuple) -> np.ndarray:
        """Evaluate the ROM for a given set of parameters.

        Parameters
        ----------
        params : NamedTuple
            The parameters to evaluate the model for

        Returns
        -------
        np.ndarray
            The ROM evaluation of the model
        """
        if self.empirical_interpolant is None:
            raise exceptions.EmpiricalInterpolantNotComputedError(
                "An attempt has been made to evaluate a ROM whose EmpiricalInterpolant has not been computed.  Compute the EmpiricalInterpolant and try again."
            )
        _signal_at_nodes = self.model.compute(params, self.empirical_interpolant.nodes)
        mask = np.in1d(self.model.domain, self.empirical_interpolant.nodes)
        #self.product(self.psd[mask], _signal_at_nodes, self.empirical_interpolant.B_matrix)

        #This may need to be noise-weighted. Will have to check after passing product to EI.
        return np.dot(_signal_at_nodes, self.empirical_interpolant.B_matrix) #Used to be np.real(b_matrix), but this is wrong.

    @log.callable("Refining ROM")
    def refine(
        self,
        n_random: int = DEFAULT_REFINE_N_RANDOM,
        tol: float = DEFAULT_TOLERANCE,
        iterate: bool = True,
    ) -> Self:
        """Refine the model by attempting to add new samples to it.

        Parameters
        ----------
        n_random : int
            Number of random samples to generate per iteration
        tol : float
            The absolute tolerance to use when evaluating the errors of each sample
        iterate : bool
            Flag that sets whether to iteratively refine until no new samples are added.

        Returns
        -------
        Self
            Returns self so that methods can be chained
        """

        if self.reduced_basis is None:
            self.reduced_basis = ReducedBasis().compute(
                self.model, Samples(self.model, n_random=n_random), tol=tol, psd=self.psd, product=self.product
            )
        self._validate_and_refine_basis(n_random, tol=tol, iterate=iterate)

        self.empirical_interpolant = EmpiricalInterpolant().compute(self.reduced_basis)

        return self

    def write(self, filename: str) -> None:
        """Save the ROM to a Rombus HDF5 file.

        Parameters
        ----------
        filename : str
            Filename of the output file
        """
        with log.context(f"Writing ROM to file ({filename})"), h5py.File(
            filename, "w"
        ) as h5file:
            self.model.write(h5file)
            self.samples.write(h5file)
            if self.reduced_basis is not None:
                self.reduced_basis.write(h5file)
            if self.empirical_interpolant is not None:
                self.empirical_interpolant.write(h5file)

    def timing(self, samples: Samples) -> float:
        """Generate timing information for the ROM.  Particularly useful when compared to
        similar timing information computed for the source model it is derived from.

        Parameters
        ----------
        samples : "Samples"
            A set of parameters to generate timing information for.  Should be the same as those used when
            timiing the source model, if comparisons are to be made.

        Returns
        -------
        float
            Seconds elapsed
        """

        with log.context(
            f"Computing timing information for ROM using {samples.n_samples} samples",
            time_elapsed=False,
        ):
            start_time = timeit.default_timer()
            for i, sample in enumerate(samples.samples):
                params_numpy = self.model.params.np2param(sample)
                _ = self.evaluate(params_numpy)
        return timeit.default_timer() - start_time

    def _validate_and_refine_basis(
        self, n_random: int, tol: float = DEFAULT_TOLERANCE, iterate: bool = True
    ) -> None:

        """Perform ROM refinement.

        Parameters
        ----------
        n_random : int
            Number of random samples to generate per iteration
        tol : float
            Absolute tolerance to use when assessing errors
        iterate : bool
            Flag that sets whether to iteratively refine until no new samples are added.
        """
        if not self.reduced_basis:
            self.reduced_basis = ReducedBasis().compute(
                self.model, self.samples, tol=tol, psd=self.psd, product=self.product
            )

        if self.reduced_basis is None:
            raise exceptions.ReducedBasisNotComputedError(
                "A ROM's reduced basis uncomputed when trying to refine basis"
            )
        else:
            n_selected_greedy_points_global = np.iinfo(np.int32).max
            n_greedy_last = len(self.reduced_basis.greedypoints)
            
            n_greedy_last_global = mpi.COMM.allreduce(n_greedy_last, op=mpi4py.MPI.SUM)
            while True:
                # generate validation set by randomly sampling the parameter space
                new_samples = Samples(self.model, n_random=n_random)
                my_vs = self.model.generate_model_set(new_samples, self.product, self.psd)

                # test validation set
                RB_transpose = np.transpose(self.reduced_basis.matrix)
                selected_greedy_points = []
                for i, validation_sample in enumerate(new_samples.samples):
                    if self.model.ordinate.dtype == complex:
                        proj_error = 1 - np.sum(
                            [
                                np.real(np.conjugate(d_i) * d_i)
                                for d_i in self.product(self.psd, my_vs[i], RB_transpose)
                            ]
                        )
                    else:
                        proj_error = 1 - np.sum(self.product(self.psd, my_vs[i], RB_transpose) ** 2)
                    if proj_error > tol:
                        selected_greedy_points.append(validation_sample)
                n_selected_greedy_points_global = mpi.COMM.allreduce(
                    len(selected_greedy_points), op=mpi4py.MPI.SUM
                )
                log.comment(
                    f"Number of samples added: {n_selected_greedy_points_global}"
                )

                # add the inaccurate points to the original selected greedy
                # points and remake the basis
                self.samples.extend(selected_greedy_points)
                self.reduced_basis = ReducedBasis().compute(
                    self.model, self.samples, tol=tol, psd=self.psd, product=self.product
                )
                n_greedy_new = len(self.reduced_basis.greedypoints)
                n_greedy_new_global = mpi.COMM.allreduce(
                    n_greedy_new, op=mpi4py.MPI.SUM
                )

                if not iterate or n_greedy_new_global == n_greedy_last_global:
                    break
                else:
                    log.comment(
                        f"Current number of accepted greedy points: {n_greedy_new_global}"
                    )
                    n_greedy_last_global = n_greedy_new_global
