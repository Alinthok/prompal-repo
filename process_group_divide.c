#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]) {
    // Get the rank and size in the original communicator
    int world_rank, world_size;
    MPI_Init(&argc, &argv); // Initializations
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the group of processes in MPI_COMM_WORLD
    MPI_Group world_group, subset_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    int n = 8;
    const int ranks1[4] = { 0, 1, 2, 3 };
    const int ranks2[4] = { 4, 5, 6, 7 };

    if (world_rank < 4) {
        MPI_Group_incl(world_group, 4, ranks1, &subset_group);
    }
    else {
        MPI_Group_incl(world_group, 4, ranks2, &subset_group);
    }

    int sendbuf, recvbuf;
    MPI_Comm new_comm;

    sendbuf = world_rank;

    // create new new communicator and then perform collective communications
    MPI_Comm_create(MPI_COMM_WORLD, subset_group, &new_comm);
    MPI_Allreduce(&sendbuf, &recvbuf, 1, MPI_INT, MPI_SUM, new_comm);


    int new_rank;

    // get rank in new group
    MPI_Group_rank(subset_group, &new_rank);
    printf("rank= %d newrank= %d recvbuf= %d\n", world_rank, new_rank, recvbuf);

    MPI_Finalize();
    return 0;
}