#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

#define MAX_NODES 1000
#define MAX_EDGES 5000
#define DAMPING 0.85
#define EPSILON 0.0001
#define MAX_ITER 100

int num_nodes, num_edges;
double IP[MAX_NODES], newIP[MAX_NODES];
int out_degree[MAX_NODES], in_degree[MAX_NODES];
double psi[MAX_NODES][MAX_NODES]; 
int followers[MAX_NODES][MAX_NODES];
int followers_count[MAX_NODES];

void initialize_graph() {
    FILE *f = fopen("edges.txt", "r");
    int u, v;
    while (fscanf(f, "%d %d", &u, &v) == 2) {
        psi[u][v] = 1.0; 
        out_degree[u]++;
        in_degree[v]++;
        followers[v][followers_count[v]++] = u;
    }
    fclose(f);

    for (int i = 0; i < num_nodes; i++) {
        IP[i] = 1.0 / num_nodes;
    }
}

void compute_influence_power() {
    int iter = 0;
    double delta;

    do {
        delta = 0.0;
        #pragma omp parallel for reduction(+:delta)
        for (int u = 0; u < num_nodes; u++) {
            double sum = 0.0;
            for (int i = 0; i < followers_count[u]; i++) {
                int f = followers[u][i];
                if (out_degree[f] > 0) {
                    sum += psi[f][u] * IP[f] / out_degree[f];
                }
            }
            newIP[u] = DAMPING * sum + (1.0 - DAMPING) / (in_degree[u] > 0 ? in_degree[u] : 1);
            delta += fabs(newIP[u] - IP[u]);
        }

        for (int i = 0; i < num_nodes; i++) {
            IP[i] = newIP[i];
        }

        iter++;
    } while (delta > EPSILON && iter < MAX_ITER);

    printf("Converged in %d iterations.\n", iter);
}

int main() {
    printf("Enter number of nodes and edges: ");
    scanf("%d %d", &num_nodes, &num_edges);

    initialize_graph();
    compute_influence_power();

    for (int i = 0; i < num_nodes; i++) {
        printf("Node %d Influence Power: %.6f\n", i, IP[i]);
    }

    return 0;
}
