#include <mpi.h>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    std::vector<int> sendbuf;  // Only allocated on rank 0
    int recvbuf;

    if (rank == 0) {
        sendbuf.resize(nproc);
        for (int i = 0; i < nproc; ++i) {
            sendbuf[i] = i + 1;
        }
    }

    // Use sendbuf.data() on rank 0, NULL on others
    MPI_Scatter(rank == 0 ? sendbuf.data() : nullptr, 1, MPI_INT,
                &recvbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::cout << "Process " << rank << " received value " << recvbuf << std::endl;

    MPI_Finalize();
    return 0;
}

