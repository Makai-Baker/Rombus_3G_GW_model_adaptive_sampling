import h5py  # type: ignore

import rombus.exceptions as exceptions
from typing import TypeVar, TypeAlias
from rombus._core.log import log

# The following are used for type hints
File: TypeAlias = h5py.File
Filename: TypeAlias = str
FileOrFilename = TypeVar("FileOrFilename", Filename, File)


def ensure_open(file_in: FileOrFilename) -> File:
    try:
        if type(file_in) == Filename:
            return h5py.File(file_in, "r"), True
        elif type(file_in) == File:
            return file_in, False
        else:
            raise exceptions.RombusHDF5Error(
                f"An attempt to open an invalid type ({type(file_in)}) as an HDF5 file was encountered."
            )
    except (IOError, exceptions.RombusException) as e:
        log.handle_exception(e)
