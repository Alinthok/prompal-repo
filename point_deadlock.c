#include <mpi.h>

int main(int argc, char* argv[]) {
    int rank, size;

    int sendbuf[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    int recvbuf[10];

    MPI_Init(&argc, &argv); // Initializations
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Status status;

    if (rank == 0) {
        MPI_Recv(recvbuf, 10, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
        MPI_Send(sendbuf, 10, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }
    else if (rank == 1) {
        MPI_Recv(recvbuf, 10, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Send(sendbuf, 10, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
