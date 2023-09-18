import numpy as np

import scipy as sp

from typing import Self

from scipy.linalg.lapack import get_lapack_funcs  # type: ignore

from rombus.params import time_to_freq

import rombus._core.hdf5 as hdf5

from rombus._core.log import log

from rombus.reduced_basis import ReducedBasis

DEFAULT_TOLERANCE = 1e-14
DEFAULT_REFINE_N_RANDOM = 100


def _malloc(dtype, *nums):
    """Allocate some memory with given dtype"""
    return np.zeros(tuple(nums), dtype=dtype)


class _LinAlg:
    """Linear algebra functions needed for empirical interpolation class"""

    def __init__(self):
        pass

    # This is scipy.linalg's source code. For some reason my scipy doesn't recognize
    # the check_finite option, which may help with speeding up. Because this may be an
    # issue related to the latest scipy.linalg version, I'm reproducing that code here
    # to guarantee future compatibility.
    def solve_triangular(
        self,
        a,
        b,
        trans=0,
        lower=False,
        unit_diagonal=False,
        overwrite_b=False,
        debug=False,
        check_finite=True,
    ):
        """
        Solve the equation `a x = b` for `x`, assuming a is a triangular matrix.

        Parameters
        ----------
        a : (M, M) array_like
            A triangular matrix
        b : (M,) or (M, N) array_like
            Right-hand side matrix in `a x = b`
        lower : boolean
            Use only data contained in the lower triangle of `a`.
            Default is to use upper triangle.
        trans : {0, 1, 2, 'N', 'T', 'C'}, optional
            Type of system to solve:

            ========  =========
            trans     system
            ========  =========
            0 or 'N'  a x  = b
            1 or 'T'  a^T x = b
            2 or 'C'  a^H x = b
            ========  =========
        unit_diagonal : bool, optional
            If True, diagonal elements of `a` are assumed to be 1 and
            will not be referenced.
        overwrite_b : bool, optional
            Allow overwriting data in `b` (may enhance performance)
        check_finite : bool, optional
            Whether to check that the input matrices contain only finite numbers.
            Disabling may give a performance gain, but may result in problems
            (crashes, non-termination) if the inputs do contain infinities or NaNs.

        Returns
        -------
        x : (M,) or (M, N) ndarray
            Solution to the system `a x = b`.  Shape of return matches `b`.

        Raises
        ------
        Exception
            If `a` is singular

        Notes
        -----
        .. versionadded:: 0.9.0

        """

        if check_finite:
            a1, b1 = map(np.asarray_chkfinite, (a, b))
        else:
            a1, b1 = map(np.asarray, (a, b))
        if len(a1.shape) != 2 or a1.shape[0] != a1.shape[1]:
            raise ValueError("expected square matrix")
        if a1.shape[0] != b1.shape[0]:
            raise ValueError("incompatible dimensions")
        overwrite_b = False  # overwrite_b or _datacopied(b1, b)
        if debug:
            print("solve:overwrite_b=", overwrite_b)
        trans = {"N": 0, "T": 1, "C": 2}.get(trans, trans)
        (trtrs,) = get_lapack_funcs(("trtrs",), (a1, b1))
        x, info = trtrs(
            a1,
            b1,
            overwrite_b=overwrite_b,
            lower=lower,
            trans=trans,
            unitdiag=unit_diagonal,
        )

        if info == 0:
            return x
        if info > 0:
            raise Exception(
                "singular matrix: resolution failed at diagonal %s" % (info - 1)
            )
        raise ValueError("illegal value in %d-th argument of internal trtrs" % -info)

    def transpose(self, a):
        dim = a.shape
        if len(dim) != 2:
            raise ValueError("Expected a matrix")
        aT = np.zeros(dim[::-1], dtype=a.dtype)
        for ii, aa in enumerate(a):
            aT[:, ii] = a[ii]
        return aT


class _EmpiricalInterpolation(_LinAlg):
    """
    Class for building an empirical interpolant
    """

    def __init__(self):
        _LinAlg.__init__(self)

    def _malloc_ei(self, Nbasis, Nquads, Nmodes=1, dtype="complex"):
        self.indices = _malloc("int", Nbasis)
        self.invV = _malloc(dtype, Nbasis, Nbasis)
        self.R = _malloc(dtype, Nbasis, Nquads)
        self.B = _malloc(dtype, Nbasis, Nquads)
        pass

    def coefficient(self, invV, e, indices):
        return np.dot(invV.T, e[indices])

    def residual(self, e, c, R):
        """Difference between a basis function 'e' and its empirical interpolation"""
        return e - np.dot(c, R)

    def next_invV_col(self, R, indices, check_finite=False):
        b = np.zeros(len(indices), dtype=R.dtype)
        b[-1] = 1.0
        return self.solve_triangular(
            R[:, indices], b, lower=False, check_finite=check_finite
        )

    def eim_interpolant(self, invV, R):
        """The empirical interpolation matrix 'B'"""
        return np.dot(invV, R)

    def eim_interpolate(self, h, indices, B):
        """Empirically interpolate a function"""
        dim = np.shape(h)
        if len(dim) == 1:
            return np.dot(h[indices], B)
        elif len(dim) > 1:
            return np.array([np.dot(h[ii][indices], B) for ii in range(dim[0])])


class _StandardEIM(_EmpiricalInterpolation):
    def __init__(self, Nbasis, Nquads, Nmodes=1, dtype="complex"):
        _EmpiricalInterpolation.__init__(self)

        # Allocate memory for numpy arrays
        self._malloc_ei(Nbasis, Nquads, Nmodes=Nmodes, dtype=dtype)
        self.modes = Nmodes
        self.quads = Nquads

    def seed(self, e):
        """Seed the algorithm"""
        self.indices[0] = np.argmax(np.abs(e))
        self.R[0] = e
        self.invV[:1, :1] = self.next_invV_col(self.R[:1], self.indices[:1])
        pass

    def iter(self, step, e):
        """One iteration in the empirical interpolation greedy algorithm"""

        ctr = step + 1

        # Compute interpolant residual
        c = self.coefficient(self.invV[:ctr, :ctr], e, self.indices[:ctr])
        r = self.residual(e, c, self.R[:ctr])

        # Update
        self.indices[ctr] = np.argmax(np.abs(r))
        self.R[ctr] = r
        self.invV[:, ctr][: ctr + 1] = self.next_invV_col(
            self.R[: ctr + 1], self.indices[: ctr + 1]
        )
        pass

    def make(self, basis, verbose=True):
        """Make an empirical interpolant using the standard greedy algorithm"""

        # Seed the greedy algorithm
        self.seed(basis[0])

        # EIM algorithm with reduced complexity for inverting the van der Monde matrix
        dim = len(basis)
        with log.progress("Iterating", dim - 1) as progress:
            for ctr in range(dim - 1):

                # Single iteration
                self.iter(ctr, basis[ctr + 1])
                progress.update(ctr)

        self.trim(ctr + 2)
        self.make_interpolant()
        self.size = len(self.indices)
        pass

    def make_interpolant(self):
        self.B = self.eim_interpolant(self.invV, self.R)
        pass

    def interpolate(self, h):
        return self.eim_interpolate(h, self.indices, self.B)

    def trim(self, num):
        """Trim zeros from remaining entries"""
        if num > 0:
            self.indices = self.indices[:num]
            self.R = self.R[:num]
            self.invV = self.invV[:num, :num]
            self.B = self.B[:num, :]
        pass


class EmpiricalInterpolant(object):
    """Class for managing the creation of empirical interpolants (EIs)"""

    def __init__(
        self, B_matrix: np.ndarray = np.ndarray([]), nodes: np.ndarray = np.ndarray([])
    ):
        """Initialise empirical interpolant"""

        self.B_matrix = B_matrix
        """Basis matrix"""

        self.nodes = nodes
        """Interpolant nodes"""

    def __str__(self):
        result = f"\nContents of instance of {type(self)}:\n"
        result = result + "B_matrix:\n"
        result = result + str(self.B_matrix) + "\n"
        result = result + "nodes: " + str(self.nodes)

        return result

    @log.callable("Computing empirical interpolant")
    def compute(self, reduced_basis: ReducedBasis) -> Self:
        """Compute empirical interpolant for a given reduced basis

        Parameters
        ----------
        reduced_basis : ReducedBasis
            Reduced basis used to compute the empirical interpolant

        Returns
        -------
        Self
            A reference to self, to allow for chaining of method calls
        """

        eim = _StandardEIM(
            reduced_basis.matrix_shape[0],
            reduced_basis.matrix_shape[1],
            dtype=reduced_basis.model.ordinate.dtype,
        )
        eim.make(reduced_basis.matrix)

        domain = reduced_basis.model.domain
        self.nodes = domain[eim.indices]

        nodes_sorted, B_matrix_sorted = zip(
            *sorted(zip(self.nodes, eim.B), key=lambda x: x[0])
        )
        self.nodes = np.asarray(nodes_sorted)
        self.B_matrix = np.asarray(B_matrix_sorted)   

        if reduced_basis.model.coordinate.resample:
            f = reduced_basis.model.coordinate
            Mc_min = reduced_basis.model.model_params['Mc'].min
            self.resample_B(f.min, Mc_min)
        
        return self

    @log.callable("Resampling the reduced basis")
    def resample_B(self, minimum_frequency, chirp_mass):
        """
        Interpolates the B matrix onto a uniformly sampled domain.
        """      
        T = time_to_freq(chirp_mass, minimum_frequency)
        
        interp_B = sp.interpolate.CubicSpline(self.nodes, self.B_matrix)
        
        self.nodes = np.linspace(f.min, f.max, (f.max-f.min)*T)
        
        self.B_matrix = interp_B(self.nodes)

        nodes_sorted, B_matrix_sorted = zip(
            *sorted(zip(self.nodes, self.B_matrix), key=lambda x: x[0])
        )
        self.nodes = np.asarray(nodes_sorted)
        self.B_matrix = np.asarray(B_matrix_sorted) 
    
    @log.callable("Writing empirical interpolant")
    def write(self, h5file: hdf5.File) -> None:
        """Write empirical interpolant to an open HDF5 file

        Parameters
        ----------
        h5file : hdf5.File
            Open HDF5 file
        """

        try:
            h5_group = h5file.create_group("empirical_interpolant")
            h5_group.create_dataset("B_matrix", data=self.B_matrix)
            h5_group.create_dataset("nodes", data=self.nodes)
        except IOError as e:
            log.handle_exception(e)

    @classmethod
    @log.callable("Instantiating empirical interpolant from file")
    def from_file(cls, file_in: hdf5.FileOrFilename) -> Self:
        """Instantiate an EmpiricalInterpolant from a Rombus file.

        Parameters
        ----------
        file_in : hdf5.FileOrFilename
            Rombus file (filename or open file) to read from

        Returns
        -------
        Self
            Returns a reference to self, so that method calls can be chained
        """

        try:
            h5file, close_file = hdf5.ensure_open(file_in)
            B_matrix = np.array(h5file["empirical_interpolant/B_matrix"])
            nodes = np.array(h5file["empirical_interpolant/nodes"])
            if close_file:
                h5file.close()
        except IOError as e:
            log.handle_exception(e)
        return cls(B_matrix=B_matrix, nodes=nodes)
