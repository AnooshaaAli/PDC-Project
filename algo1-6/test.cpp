#include <metis.h>
#include <iostream>

int main() {
    // idx_t nvtxs = 4;                         // Number of vertices
    // idx_t ncon = 1;                          // Number of balancing constraints
    // idx_t nparts = 4;                        // Number of partitions (you can change this)

    // idx_t xadj[5]    = {0, 2, 4, 6, 8};      // CSR row pointer
    // idx_t adjncy[8]  = {1, 2, 0, 3, 0, 3, 1, 2}; // CSR column indices

    // idx_t part[4];                           // Output array: which partition each vertex belongs to
    // idx_t objval;                            // Output: edge cut value

    idx_t nvtxs = 9;                         // Number of vertices
    idx_t ncon = 1;                          // Number of balancing constraints
    idx_t nparts = 9;                        // Change as needed (e.g. 2, 3...)

    // CSR representation of the undirected graph
    idx_t xadj[10] = {0, 2, 5, 7, 9, 10, 11, 14, 15, 16};
    idx_t adjncy[16] = {1, 2, 0, 2, 3, 1, 0, 1, 4, 3, 6, 5, 7, 8, 6, 6};

    idx_t part[9];       // Output partition array
    idx_t objval;        // Edge cut value

    // Call METIS
    // int result = METIS_PartGraphKway(
    //     &nvtxs, &ncon, xadj, adjncy, 
    //     NULL, NULL, NULL, 
    //     &nparts, NULL, NULL, NULL, 
    //     &objval, part
    // );

    int result = METIS_PartGraphRecursive(
        &nvtxs,             // Number of vertices
        &ncon,              // Balancing constraints
        xadj,               // CSR xadj
        adjncy,             // CSR adjncy
        NULL,               // Vertex weights (optional)
        NULL,               // Vertex sizes (optional)
        NULL,               // Edge weights (optional)
        &nparts,            // Number of partitions
        NULL,               // Desired weight for each partition (optional)
        NULL,               // Imbalance tolerance (optional)
        NULL,               // Options (optional)
        &objval,            // Edge cut value
        part                // Partition array (output)
    );

    if (result == METIS_OK) {
        std::cout << "Partitioning successful! Edge cut = " << objval << std::endl;
        for (int i = 0; i < nvtxs; ++i) {
            std::cout << "Vertex " << i << " in partition " << part[i] << std::endl;
        }
    } else {
        std::cout << "Partitioning failed!" << std::endl;
    }

    return 0;
}
