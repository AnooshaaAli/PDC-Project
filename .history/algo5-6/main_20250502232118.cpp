#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <cmath>
#include <omp.h>
#include <set>
#include <algorithm>
#include <queue>

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
int num_nodes = 0;

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
    // cout << inter << endl;
    // cout << uni << endl;
    return (uni == 0) ? 0.0 : (double)inter / uni;
}

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
    // cout << "U: " << u << " V: " << v << " C: " << sim << endl;

    ActionCount &a = action_map[u][v];
    // cout << "Num Retweets: " << a.RT << ", Mentions: " << a.MT << ", Num Replies: " << a.RE << endl;
    double weighted = get_action_weight("RT") * sim * a.RT + get_action_weight("RE") * sim * a.RE + get_action_weight("MT") * sim * a.MT;
    return weighted / total_posts[v];
}

void initialize_ip() {
    set<int> all_nodes;

    for (const auto &[u, nbrs] : graph) {
        all_nodes.insert(u);
        for (int v : nbrs)
            all_nodes.insert(v);
    }

    for (const auto &[v, f_list] : followers) {
        all_nodes.insert(v);
        for (int f : f_list)
            all_nodes.insert(f);
    }

    num_nodes = all_nodes.size();

    for (int u : all_nodes) {
        IP[u] = 1.0 / num_nodes;
    }
}


void compute_influence_power(const vector<int> &nodes) {
    int iter = 0;
    double delta;

    do {
        delta = 0.0;
        for (int u : nodes) {
            double sum = 0.0;
            for (int f : followers[u]) {
                if (out_degree[f] > 0) {
                    double psi_val = compute_psi(f, u);
                    // cout << psi_val << endl;
                    // cout << f << out_degree[f] << endl;
                    sum += psi_val * IP[f] / out_degree[f];
                }
            }
            double base = (in_degree[u] > 0) ? in_degree[u] : 0.0;
            newIP[u] = DAMPING * sum + ((1.0 - DAMPING) * base) / num_nodes;
            // cout << newIP[u] << endl;
            delta += fabs(newIP[u] - IP[u]);
        }

        for (int u : nodes)
            IP[u] = newIP[u];

        iter++;

    } while (delta > EPSILON && iter < MAX_ITER);
}

void calculate_influence_power() {
    vector<int> levels;
    for (const auto &[level, _] : level_components)
        levels.push_back(level);
    sort(levels.begin(), levels.end());  

    for (int level : levels) {
        cout << "Processing Level " << level << "...\n";
        const auto &comps = level_components[level];

        #pragma omp parallel for
        for (int i = 0; i < comps.size(); ++i) {
            int thread_id = omp_get_thread_num();
            cout << "[OMP Thread " << thread_id << "] Processing Component " << comps[i] << endl;
            compute_influence_power(component_nodes[comps[i]]);
        }
    }
}

void output_IP() {
    for (const auto &[u, ip] : IP) {
        cout << "Node " << u << ": IP = " << ip << endl;
    }
}

unordered_map<int, double> compute_IL(int v) {
    unordered_map<int, double> IL;
    unordered_map<int, int> level_count;
    unordered_map<int, int> distance;
    queue<int> q;

    q.push(v);
    distance[v] = 0;

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        int lvl = distance[u];

        for (int nbr : followers[u]) { 
            if (!distance.count(nbr)) {
                distance[nbr] = lvl + 1;
                q.push(nbr);
            }
        }
    }

    for (const auto& [node, lvl] : distance) {
        if (node != v) { /
            IL[lvl] += IP[node];
            level_count[lvl]++;
        }
    }

    // Final average
    for (auto& [lvl, sum] : IL) {
        if (level_count[lvl] > 0)
            sum /= level_count[lvl];
    }

    return IL;
}

vector<int> find_seed_candidates() {
    vector<int> seeds;

    for (const auto& [v, _] : IP) {
        unordered_map<int, double> IL = compute_IL(v);

        vector<int> gamma;
        int max_lvl = 0;
        for (const auto& [lvl, _] : IL) {
            if (lvl > max_lvl) max_lvl = lvl;
        }

        for (int L = 1; L < max_lvl; ++L) {
            if (IL.count(L) && IL.count(L + 1)) {
                if (IL[L] > IL[L + 1]) {
                    gamma.push_back(L);
                }
            }
        }

        int L0 = gamma.empty() ? max_lvl : *min_element(gamma.begin(), gamma.end());

        double IL0 = IL.count(L0) ? IL[L0] : 0.0;

        cout << "\nNode " << v << ": IP = " << IP[v] << ", IL0 = " << IL0 << ", L0 = " << L0 << endl;
        if (IP[v] > IL0) {
            cout << "  ✔️ Node " << v << " selected as seed.\n";
            seeds.push_back(v);
        } else {
            cout << "  ❌ Node " << v << " NOT selected.\n";
        }
    }

    return seeds;
}

int main() {
    load_interest("interest.txt");
    load_graph("edges.txt");
    load_activity("activity.txt");
    load_component_map("components.txt");
    load_component_levels("component_levels.txt");
    initialize_ip();
    calculate_influence_power();
    output_IP();
    auto seeds = find_seed_candidates();
    cout << "\nSeed Candidates:\n";
    for (int s : seeds) {
        cout << "Node " << s << " (IP = " << IP[s] << ")\n";
    }
    return 0;
}
