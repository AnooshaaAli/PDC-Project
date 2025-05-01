#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <math.h>

#define MAX_NODES 1000
#define MAX_TOPICS 10
#define MAX_FOLLOWERS 100
#define DAMPING 0.85
#define EPSILON 0.0001
#define MAX_ITER 100

int num_nodes = 0;
double IP[MAX_NODES], newIP[MAX_NODES];
int out_degree[MAX_NODES], in_degree[MAX_NODES];
double psi[MAX_NODES][MAX_NODES];
int followers[MAX_NODES][MAX_FOLLOWERS];
int followers_count[MAX_NODES];
int interests[MAX_NODES][MAX_TOPICS];
int interest_len = 0;
int total_posts[MAX_NODES]; // N_{p_y}: posts by each user

// Action counters
typedef struct {
    int RT, RE, MT;
} ActionCount;
ActionCount action_matrix[MAX_NODES][MAX_NODES];

// Get α_i: action weights
double get_action_weight(const char *action) {
    if (strcmp(action, "RT") == 0) return 0.5;
    if (strcmp(action, "RE") == 0) return 0.35;
    if (strcmp(action, "MT") == 0) return 0.15;
    return 0.0;
}

// Jaccard similarity between two interest vectors
double jaccard(int *a, int *b, int len) {
    int intersect = 0, uni = 0;
    for (int i = 0; i < len; i++) {
        if (a[i] || b[i]) uni++;
        if (a[i] && b[i]) intersect++;
    }
    return (uni == 0) ? 0.0 : (double)intersect / uni;
}

void read_interest_vectors(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        perror("interest.txt");
        exit(1);
    }
    int user, val;
    while (!feof(f)) {
        int count = 0;
        if (fscanf(f, "%d", &user) != 1) break;
        while (count < MAX_TOPICS && fscanf(f, "%d", &val) == 1) {
            interests[user][count++] = val;
        }
        if (count > interest_len) interest_len = count;
        if (user >= num_nodes) num_nodes = user + 1;
    }
    fclose(f);
}

void read_edges(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        perror("edges.txt");
        exit(1);
    }
    int u, v;
    while (fscanf(f, "%d %d", &u, &v) == 2) {
        followers[v][followers_count[v]++] = u;
        out_degree[u]++;
        in_degree[v]++;
        if (u >= num_nodes) num_nodes = u + 1;
        if (v >= num_nodes) num_nodes = v + 1;
    }
    fclose(f);
}

void read_activities(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        perror("activity.txt");
        exit(1);
    }

    int u, v;
    long ts;
    char action[4];

    while (fscanf(f, "%d %d %ld %s", &u, &v, &ts, action) == 4) {
        if (strcmp(action, "RT") == 0) action_matrix[u][v].RT++;
        else if (strcmp(action, "RE") == 0) action_matrix[u][v].RE++;
        else if (strcmp(action, "MT") == 0) action_matrix[u][v].MT++;

        total_posts[v]++;  // Approximates number of posts by v

        if (u >= num_nodes) num_nodes = u + 1;
        if (v >= num_nodes) num_nodes = v + 1;
    }

    fclose(f);

    // Now compute ψ(u, v)
    for (int u = 0; u < num_nodes; u++) {
        for (int v = 0; v < num_nodes; v++) {
            int Npy = total_posts[v];
            if (Npy == 0) continue;

            double sim = jaccard(interests[u], interests[v], interest_len);
            double weighted_sum = 0.0;

            weighted_sum += 0.5 * sim * action_matrix[u][v].RT;
            weighted_sum += 0.35 * sim * action_matrix[u][v].RE;
            weighted_sum += 0.15 * sim * action_matrix[u][v].MT;

            psi[u][v] = weighted_sum / Npy;
        }
    }
}

void initialize_ip() {
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
                if (out_degree[f] > 0 && psi[f][u] > 0) {
                    sum += psi[f][u] * IP[f] / out_degree[f];
                }
            }
            int f_count = (in_degree[u] > 0) ? in_degree[u] : 1;
            newIP[u] = DAMPING * sum + (1.0 - DAMPING) / f_count;
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
    read_interest_vectors("interest.txt");
    read_edges("edges.txt");
    read_activities("activity.txt");
    initialize_ip();
    compute_influence_power();

    printf("\nFinal Influence Power Scores:\n");
    for (int i = 0; i < num_nodes; i++) {
        printf("Node %d: %.6f\n", i, IP[i]);
    }

    return 0;
}
