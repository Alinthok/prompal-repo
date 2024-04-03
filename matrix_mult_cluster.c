#include <iostream>
#include <ctime>
#include <cstdlib>
#include <mpi.h>

using namespace std;

void RowMatrixVectorMultiply(int dim, double* matrix_data, double* vector_data, double* result);

int main(int argc, char* argv[]) {
    int rank, size; // MPI RANK, SIZE

    MPI_Init(&argc, &argv); // Initializations
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int N_DIM = 128; // confirmed
    if (N_DIM % size) { // Valid Communicator Size?
        MPI_Finalize();
        return 0;
    }

    // Allocate memory for matrix, vector, and result
    double* matrix_data = new double[N_DIM * N_DIM];
    double* vector_data = new double[N_DIM];
    double* result = new double[N_DIM];

    // Initialize matrix and vector with random numbers on rank 0
    if (rank == 0) {
        srand(time(NULL)); // Seed for random numbers
        for (int i = 0; i < N_DIM; i++) {
            for (int j = 0; j < N_DIM; j++) {
                matrix_data[i * N_DIM + j] = rand() % 10; // Fill matrix with random numbers
            }
            vector_data[i] = rand() % 10; // Fill vector with random numbers
        }
    }

    RowMatrixVectorMultiply(N_DIM, matrix_data, vector_data, result);

    // Free dynamically allocated memory
    delete[] matrix_data;
    delete[] vector_data;
    delete[] result;

    MPI_Finalize();
    return 0;
}

void RowMatrixVectorMultiply(int dim, double* matrix_data, double* vector_data, double* result) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double* localresult = new double[dim / size]{};
    double* local_matrix = new double[(dim * dim) / size]; // local matrix

    double start = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatter(matrix_data, (dim * dim) / size, MPI_DOUBLE, local_matrix, (dim * dim) / size, MPI_DOUBLE, 0,
        MPI_COMM_WORLD); // Scatter the Matrix
    MPI_Bcast(vector_data, dim, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Broadcast the Vector

    double calc_timer = MPI_Wtime();

    // Calculate the results
    for (int i = 0; i < (dim / size); i++)
        for (int j = 0; j < dim; j++)
            localresult[i] += vector_data[j] * local_matrix[i * dim + j];

    calc_timer = MPI_Wtime() - calc_timer;
    cout << "Rank " << rank << " - CALC TIME = " << calc_timer << " seconds" << endl;

    MPI_Gather(localresult, (dim) / size, MPI_DOUBLE, result, (dim) / size, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Gather the results

    MPI_Barrier(MPI_COMM_WORLD);

    double end = MPI_Wtime();

    double max_time, min_start_time, max_calc_time;

    MPI_Reduce(&start, &min_start_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&end, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&calc_timer, &max_calc_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    // Print communication time on rank 0
    if (rank == 0) {
        double total_runtime = max_time - min_start_time;
        cout << "Time Needed for communication = " << total_runtime - max_calc_time << " seconds" << endl;
        cout << "Time Needed for computation = " << max_calc_time << " seconds" << endl;
        cout << "x/y = " << max_calc_time/(total_runtime - max_calc_time) << endl;
    }
    delete[] localresult;
    delete[] local_matrix;
}