from mpi4py import MPI

MAIN_RANK = 0

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

RANK_IS_MAIN = RANK == MAIN_RANK
