#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <omp.h>
#include <unordered_map>

#define MAX_NODES 1000
#define MAX_TOPICS 10
#define DAMPING 0.85
#define EPSILON 0.0001
#define MAX_ITER 100

using namespace std;

int num_nodes = 0;
double IP[MAX_NODES], newIP[MAX_NODES];
int out_degree[MAX_NODES], in_degree[MAX_NODES];
double psi[MAX_NODES][MAX_NODES];
int total_posts[MAX_NODES];
int interests[MAX_NODES][MAX_TOPICS];
int interest_len = 0;
vector<int> followers[MAX_NODES];

struct ActionCount {
    int RT = 0, RE = 0, MT = 0;
};

ActionCount action_matrix[MAX_NODES][MAX_NODES];

double get_action_weight(const string &action) {
    if (action == "RT") return 0.5;
    if (action == "RE") return 0.35;
    if (action == "MT") return 0.15;
    return 0.0;
}

double jaccard(int* a, int* b, int len) {
    int intersect = 0, uni = 0;
    for (int i = 0; i < len; i++) {
        if (a[i] || b[i]) uni++;
        if (a[i] && b[i]) intersect++;
    }
    return (uni == 0) ? 0.0 : (double)intersect / uni;
}

void read_interest_vectors(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Could not open interest.txt\n";
        exit(1);
    }

    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        int user, val;
        ss >> user;
        int count = 0;
        while (ss >> val && count < MAX_TOPICS) {
            interests[user][count++] = val;
        }
        if (count > interest_len) interest_len = count;
        if (user >= num_nodes) num_nodes = user + 1;
    }
}

void read_edges(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Could not open edges.txt\n";
        exit(1);
    }

    int u, v;
    while (file >> u >> v) {
        followers[v].push_back(u);
        out_degree[u]++;
        in_degree[v]++;
        if (u >= num_nodes) num_nodes = u + 1;
        if (v >= num_nodes) num_nodes = v + 1;
    }
}

void read_activities(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Could not open activity.txt\n";
        exit(1);
    }

    int u, v;
    long ts;
    string action;

    while (file >> u >> v >> ts >> action) {
        if (action == "RT") action_matrix[u][v].RT++;
        else if (action == "RE") action_matrix[u][v].RE++;
        else if (action == "MT") action_matrix[u][v].MT++;

        total_posts[v]++;
        if (u >= num_nodes) num_nodes = u + 1;
        if (v >= num_nodes) num_nodes = v + 1;
    }

    // Compute ψ(u, v)
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
            for (int f : followers[u]) {
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

    cout << "Converged in " << iter << " iterations." << endl;
}

int main() {
    read_interest_vectors("interest.txt");
    read_edges("edges.txt");
    read_activities("activity.txt");
    initialize_ip();
    compute_influence_power();

    cout << "\nFinal Influence Power Scores:\n";
    for (int i = 0; i < num_nodes; i++) {
        cout << "Node " << i << ": " << IP[i] << endl;
    }

    return 0;
}
