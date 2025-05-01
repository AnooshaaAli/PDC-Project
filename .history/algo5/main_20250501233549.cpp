#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <cmath>
#include <omp.h>
#include <set>

#define MAX_TOPICS 10
#define DAMPING 0.85
#define EPSILON 0.0001
#define MAX_ITER 100

using namespace std;

struct ActionCount {
    int RT = 0, RE = 0, MT = 0;
};

unordered_map<int, vector<int>> graph;
unordered_map<int, vector<int>> followers;
unordered_map<int, int> out_degree, in_degree;
unordered_map<int, unordered_map<int, ActionCount>> action_map;
unordered_map<int, vector<int>> interest_map;
unordered_map<int, int> total_posts;
unordered_map<int, double> IP, newIP;
unordered_map<int, double> psi;
unordered_map<int, int> node_to_component;
unordered_map<int, int> component_level;
unordered_map<int, vector<int>> component_nodes;
unordered_map<int, vector<int>> level_components;

int interest_len = 0;

double get_action_weight(const string &act) {
    if (act == "RT") return 0.5;
    if (act == "RE") return 0.35;
    if (act == "MT") return 0.15;
    return 0.0;
}

double jaccard(const vector<int> &a, const vector<int> &b) {
    int inter = 0, uni = 0;
    for (int i = 0; i < interest_len; i++) {
        if (a[i] || b[i]) uni++;
        if (a[i] && b[i]) inter++;
    }
    return (uni == 0) ? 0.0 : (double)inter / uni;
}

// Read interest vectors
void load_interest(const string &file) {
    ifstream f(file);
    string line;
    while (getline(f, line)) {
        stringstream ss(line);
        int user, bit;
        ss >> user;
        vector<int> vec;
        while (ss >> bit) vec.push_back(bit);
        interest_map[user] = vec;
        interest_len = vec.size();
    }
}

// Read edges
void load_graph(const string &file) {
    ifstream f(file);
    int u, v;
    while (f >> u >> v) {
        graph[u].push_back(v);
        followers[v].push_back(u);
        out_degree[u]++;
        in_degree[v]++;
    }
}

// Read activity
void load_activity(const string &file) {
    ifstream f(file);
    int u, v;
    long ts;
    string act;
    while (f >> u >> v >> ts >> act) {
        if (act == "RT") action_map[u][v].RT++;
        if (act == "RE") action_map[u][v].RE++;
        if (act == "MT") action_map[u][v].MT++;
        total_posts[v]++;
    }
}

void load_component_map(const string &file) {
    ifstream f(file);
    int node, comp;
    while (f >> node >> comp) {
        node_to_component[node] = comp;
        component_nodes[comp].push_back(node);
    }
}

void load_component_levels(const string &file) {
    ifstream f(file);
    int comp, level;
    while (f >> comp >> level) {
        component_level[comp] = level;
        level_components[level].push_back(comp);
    }
}

double compute_psi(int u, int v) {
    if (total_posts[v] == 0) return 0.0;
    if (!interest_map.count(u) || !interest_map.count(v)) return 0.0;

    double sim = jaccard(interest_map[u], interest_map[v]);
    cout << "U: " << u << "V: " << v << sim << endl;

    ActionCount &a = action_map[u][v];
    double weighted = 0.5 * sim * a.RT + 0.35 * sim * a.RE + 0.15 * sim * a.MT;
    return weighted / total_posts[v];
}

// Initialize IP
void initialize_ip() {
    for (const auto &[u, _] : graph) {
        IP[u] = 1.0 / graph.size();
    }
}

// PageRank within one component
void compute_influence_power(const vector<int> &nodes) {
    int iter = 0;
    double delta;
    do {
        delta = 0.0;
        #pragma omp parallel for reduction(+:delta)
        for (int i = 0; i < nodes.size(); ++i) {
            int u = nodes[i];
            double sum = 0.0;
            for (int f : followers[u]) {
                if (out_degree[f] > 0) {
                    double psi_val = compute_psi(f, u);
                    sum += psi_val * IP[f] / out_degree[f];
                }
            }
            double base = (in_degree[u] > 0) ? in_degree[u] : 1;
            newIP[u] = DAMPING * sum + (1.0 - DAMPING) / base;
            delta += fabs(newIP[u] - IP[u]);
        }
        for (int u : nodes) IP[u] = newIP[u];
        iter++;
    } while (delta > EPSILON && iter < MAX_ITER);
}

void run_algorithm_5() {
    for (const auto &[level, comps] : level_components) {
        #pragma omp parallel for
        for (int i = 0; i < comps.size(); ++i) {
            int comp = comps[i];
            compute_influence_power(component_nodes[comp]);
        }
    }
}

void output_IP() {
    for (const auto &[u, ip] : IP) {
        cout << "Node " << u << ": IP = " << ip << endl;
    }
}

int main() {
    load_interest("interest.txt");
    load_graph("edges.txt");
    load_activity("activity.txt");
    load_component_map("component.txt");
    load_component_levels("component_levels.txt");
    initialize_ip();
    run_algorithm_5();
    output_IP();
    return 0;
}
