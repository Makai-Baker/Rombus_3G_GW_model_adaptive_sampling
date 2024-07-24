from typing import Any, Dict, List, Optional, Self, TypeAlias

import numpy as np

import rombus._core.mpi as mpi

from rombus.model import RombusModel
from rombus._core import hdf5

DEFAULT_TOLERANCE = 1e-14
DEFAULT_REFINE_N_RANDOM = 100

Sample: TypeAlias = np.ndarray

def _merge_dicts(x,y):
    return {k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y)}

class Samples(object):
    """Class for managing sets of parameter samples for Rombus."""

    def __init__(
        self, model: RombusModel, filename: Optional[str] = None, n_random: int = 0
    ):

        self.model: RombusModel = model
        """Rombus model for which the samples are computed"""

        self.n_random: int = n_random
        """Number of random points generated for this set"""

        # RNG current and initial state
        self._random: Optional[np.random._generator.Generator] = None
        """The numpy random number generator instance used for generating any random numbers in this Sample set"""
        self._random_starting_state: Optional[Dict[str, Any]] = None
        """The starting state of the generator"""

        # Initialise samples
        self.n_samples: np.int32 = np.int32(0)
        """Number of samples in this set."""
        self.samples: List[Sample] = []
        """List of samples in this set."""
        self.counts: Dict[int] = {}
        """Dict of num samples split onto mpi ranks"""

        if filename:
            self._add_from_file(filename)
        if self.n_random > 0:
            self._add_random_samples(self.n_random)

    def __str__(self):
        result = f"\nContents of instance of {type(self)}:\n"
        result = result + "   model:    " + str(self.model) + "\n"
        result = result + "   n_random: " + str(self.n_random) + "\n"
        result = result + "   n_samples:" + str(self.n_samples) + "\n"
        result = result + "   samples:  " + str(self.samples)

        return result

    @classmethod
    def from_file(cls, file_in: hdf5.FileOrFilename) -> Self:
        """Create an instance of a Sample set from a Rombus file on disk.

        Parameters
        ----------
        file_in : hdf5.FileOrFilename
            Rombus file (filename or opened file) to read from

        Returns
        -------
        Self
            Returns a reference to self so that method calls can be chained
        """

        h5file, close_file = hdf5.ensure_open(file_in)
        model_str = h5file["samples/model/model_str"].asstr()[()]
        model = RombusModel.load(model_str)
        samples = cls(model)
        samples.samples = [np.array(x) for x in h5file["samples/samples"]]
        samples.n_samples = np.int32(h5file["samples/n_samples"])
        if close_file:
            h5file.close()
        return samples

    def extend(self, new_samples: List[Sample]) -> None:
        """Add additional samples to the set.

        Parameters
        ----------
        new_samples : List[Sample]
            A list of new samples
        """
        new_samples, counts = self._decompose_samples(new_samples)
        n_new_samples = len(new_samples)

        self.samples.extend(new_samples)
        self.n_samples = self.n_samples + n_new_samples

        self.counts = _merge_dicts(self.counts,counts)
        #self.samples.extend(new_samples)
        #self.n_samples = self.n_samples + len(new_samples)

    def write(self, h5file: hdf5.File):
        """Save samples to an open HDF5 file.

        Parameters
        ----------
        h5file : hdf5.File
            An open HDF5 file
        """

        h5_group = h5file.create_group("samples")
        self.model.write(h5_group)
        h5_group.create_dataset("samples", data=self.samples)
        h5_group.create_dataset("n_samples", data=self.n_samples)

    def _add_from_file(self, filename_in: str) -> None:
        """Add samples from file to this set.  Accepts Numpy or CSV files.

        Parameters
        ----------
        filename_in : str
            Filename of a Numpy or CSV file.
        """

        # dividing greedypoints into chunks
        if mpi.RANK_IS_MAIN:
            if filename_in.endswith(".npy"):
                samples = np.load(filename_in, allow_pickle=True)
            elif filename_in.endswith(".csv"):
                samples = [
                    np.atleast_1d(x)
                    for x in np.genfromtxt(filename_in, delimiter=",", comments="#")
                ]
            else:
                raise Exception
        else:
            samples = None

        new_samples, counts = self._decompose_samples(samples)
        n_new_samples = len(new_samples)

        self.samples.extend(new_samples)
        self.n_samples = self.n_samples + n_new_samples

        self.counts = _merge_dicts(self.counts,counts)

    def _add_random_samples(self, n_samples: int) -> None:
        """Add randomly generated samples to the set.

        Parameters
        ----------
        n_samples : int
            Number of random samples to add
        """

        self._random = np.random.default_rng()
        self._random_starting_state = np.random.get_state()
        samples = []
        for _ in range(n_samples):
            new_sample = self.model.params.generate_random_sample(self._random)
            samples.append(new_sample)

        new_samples, counts = self._decompose_samples(samples)
        n_new_samples = len(new_samples)

        self.samples.extend(new_samples)
        self.n_samples = self.n_samples + n_new_samples

        self.counts = _merge_dicts(self.counts,counts)

    def _decompose_samples(
        self,
        samples: List[Sample],
    ) -> List[Sample]:
        """Split a list of samples accross MPI ranks.

        Parameters
        ----------
        samples : List[Sample]
            Set of samples to split

        Returns
        -------
        List[Sample]
            Set of samples selected for the local rank
        """

        chunks: List[List[Sample]] = [[]]
        chunk_counts = None                                                 #Added
        if mpi.RANK_IS_MAIN:
            chunks = [[] for _ in range(mpi.SIZE)]
            for i, chunk in enumerate(samples):
                chunks[i % mpi.SIZE].append(chunk)
            chunk_counts = {i: len(chunks[i]) for i in range(len(chunks))}  #Added

            chunks, chunk_counts = self._balance_chunks(chunks, chunk_counts)
            
        training_samples = mpi.COMM.scatter(chunks, root=mpi.MAIN_RANK)         #Added
        chunk_counts = mpi.COMM.bcast(chunk_counts, root=mpi.MAIN_RANK)         #Added

        return training_samples, chunk_counts
        #return mpi.COMM.scatter(chunks, root=mpi.MAIN_RANK)

    def _balance_chunks(
        self, 
        chunks, 
        counts: Dict
    ) -> List[Sample]:
        """Balance the number of samples on the cores

        Parameters
        ----------
        samples : Array
            Set of samples to split
        
        counts : Dict
            Counts of the number of samples on each core

        Returns
        -------
        List[Sample]
            Samples
        """
        max_samples = max(counts.values())
        for i, chunk in enumerate(chunks):
            while counts[i] < max_samples:
                new_sample = self.model.params.generate_random_sample(self._random)
                chunks[i].append(new_sample)
                counts[i] += 1 

        return chunks, counts