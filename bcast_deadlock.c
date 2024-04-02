#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]) {
    int rank, size;

    int buf1[1];
    int buf2[1];

    MPI_Init(&argc, &argv); // Initializations
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        buf1[0] = 1;
        MPI_Bcast(buf1, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(buf2, 1, MPI_INT, 1, MPI_COMM_WORLD);
        std::cout << "rank0 : " << buf2[0];
    }
    else if (rank == 1) {
        buf2[0] = 2;
        MPI_Bcast(buf2, 1, MPI_INT, 1, MPI_COMM_WORLD);
        MPI_Bcast(buf1, 1, MPI_INT, 0, MPI_COMM_WORLD);
        std::cout << "rank1 : " << buf1[0];
    } //Ga deadlock kalo ada system buffer.

    MPI_Finalize();
    return 0;
}