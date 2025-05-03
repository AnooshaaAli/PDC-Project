#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <metis.h>
#include <stack>
#include <string>
#include <algorithm>

using namespace std;

struct DirectedEdge {
    int from, to;
    int retweet, reply, mention;
    int partition;
};

struct Vertex {
    int index = -1, lowlink = -1, level = -1, depth = -1;
    std::string type; // "scc" or "cac"
    int comp = -1;    // Component ID
};

struct Edge {
    int from, to;
};

void dfs(int v, const std::vector<DirectedEdge>& edges);
void discover(int v);
void explore(int v, const std::vector<DirectedEdge>& edges);
void finish(int v);

int index_counter = 0, comp_counter = 0;
std::stack<int> s;
std::vector<Vertex> vertices;
std::vector<std::vector<int>> adj;

void discover(int v) {
    vertices[v].index = index_counter;
    vertices[v].lowlink = index_counter;
    vertices[v].level = 1;
    vertices[v].depth = 1;
    index_counter++;
    s.push(v);
}

void explore(int v, const std::vector<DirectedEdge>& edges) {
    for (const auto& edge : edges) {
        if (edge.from != v) continue;

        int w = edge.to;
        if (vertices[w].index == -1) {
            dfs(w, edges);  // DFS wrapper that calls discover, explore, finish
        }

        if (!vertices[w].type.empty()) {
            vertices[v].level = std::max(vertices[w].level, vertices[v].level + 1);
        } else {
            vertices[v].level = std::max(vertices[w].level, vertices[v].level);
            vertices[v].lowlink = std::min(vertices[v].lowlink, vertices[w].lowlink);
        }
    }
}

void merge(int u, int v) {
    vertices[v].comp = comp_counter;
}

void finish(int v, const std::vector<DirectedEdge>& edges) {
    if (vertices[v].lowlink != vertices[v].index) return;

    int size = 0;
    int w;
    do {
        w = s.top(); s.pop();
        vertices[w].lowlink = vertices[v].lowlink;
        vertices[w].level = vertices[v].level;
        vertices[w].type = "scc";
        merge(v, w);
        size++;
    } while (w != v);

    if (size == 1) {
        vertices[v].type = "cac";
        bool m = false;
        std::vector<int> list;

        // ✅ Use the edge list instead of adj[v]
        for (const DirectedEdge& edge : edges) {
            if (edge.from != v) continue;
            int w = edge.to;

            if (vertices[w].type == "scc" && vertices[w].level == vertices[v].level - 1) {
                m = false;
                break;
            } else if (vertices[w].type == "cac" && vertices[w].level == vertices[v].level - 1) {
                list.push_back(w);
                m = true;
            }
        }

        if (m) {
            vertices[v].level -= 1;
            for (int w : list) {
                merge(w, v);
            }
        }
    }

    comp_counter++;
}


void dfs(int v, const std::vector<DirectedEdge>& edges) {
    discover(v);
    explore(v, edges);
    finish(v, edges);
}

void SCC_CAC_partition(const std::vector<std::vector<DirectedEdge>>& partitioned_edges) {
    for (const auto& edges : partitioned_edges) {
        for (const auto& edge : edges) {
            if (vertices[edge.from].comp == -1) {
                dfs(edge.from, edges);
            }
        }
    }
}


int main() {
    string filename = "dataset.txt";
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open file.\n";
        return 1;
    }
    cout << "Reading edges from file: " << filename << endl;
    cout << "----------------------------------------\n";
    set<int> nodes;
    vector<tuple<int, int, int, int, int>> raw_edges;

    // Step 1: Read edges and collect nodes
    int u, v, retweet, reply, mention;
    while (file >> u >> v >> retweet >> reply >> mention) {
        raw_edges.emplace_back(u, v, retweet, reply, mention);
        nodes.insert(u);
        nodes.insert(v);
    }

    cout << "Total edges: " << raw_edges.size() << endl;
    cout << "Total nodes: " << nodes.size() << endl;
    cout << "----------------------------------------\n";

    int nvtxs = *nodes.rbegin() + 1;

    // Step 2: Build undirected adjacency list
    vector<vector<int>> adjList(nvtxs);
    for (auto& [u, v, retweet, reply, mention] : raw_edges) {
        adjList[u].push_back(v);
        adjList[v].push_back(u); // undirected for METIS
    }

    cout << "Adjacency list built.\n";
    cout << "----------------------------------------\n";

    // Step 3: Convert to CSR format
    vector<idx_t> xadj(nvtxs + 1, 0);
    vector<idx_t> adjncy;
    for (int i = 0; i < nvtxs; ++i) {
        xadj[i + 1] = xadj[i] + adjList[i].size();
        for (auto& neighbor : adjList[i]) {
            adjncy.push_back(neighbor);
        }
    }

    cout << "CSR format built.\n";
    cout << "----------------------------------------\n";

    // Step 4: Set uniform node weights (optional: METIS can take nullptr here)
    vector<idx_t> vwgt(nvtxs, 1);

    cout << "Node weights set.\n";
    cout << "----------------------------------------\n";

    // Step 5: Partition using METIS
    idx_t nparts = 2;
    idx_t ncon = 1;
    vector<idx_t> part(nvtxs);
    idx_t objval;
    int result = METIS_PartGraphKway(
        &nvtxs, &ncon,
        xadj.data(), adjncy.data(), nullptr,
        nullptr, nullptr,
        &nparts, nullptr, nullptr, nullptr,
        &objval, part.data()
    );

    cout << "Partitioning...\n";
    cout << "----------------------------------------\n";

    if (result != METIS_OK) {
        cerr << "METIS partitioning failed!\n";
        return 1;
    }

    cout << "\nPartitioning successful!\nEdge cut: " << objval << "\n";
    for (int i = 0; i < nvtxs; ++i)
        cout << "Node " << i << " → Partition " << part[i] << endl;

    // Step 6: Group directed edges by source partition
    vector<vector<DirectedEdge>> partitioned_edges(nparts);
    for (auto& [u, v, retweet, reply, mention] : raw_edges) {
        if (part[u] == part[v]) { // Keep only intra-partition edges no  cross-partition edges
            DirectedEdge edge{u, v, retweet, reply, mention, part[u]};
            partitioned_edges[part[u]].push_back(edge);
        }
    }

    // Step 7: Output directed edges by partition
    for (int p = 0; p < nparts; ++p) {
        cout << "\nDirected Edges in Partition " << p << ":\n";
        for (const auto& edge : partitioned_edges[p]) {
            cout << "From " << edge.from << " → To " << edge.to
                 << " | Retweet = " << edge.retweet
                 << ", Reply = " << edge.reply
                 << ", Mention = " << edge.mention << endl;
        }
    }
    cout << "\n\n";
    
    vertices.resize(nvtxs); // assuming nvtxs is the number of nodes

    SCC_CAC_partition(partitioned_edges);

    for (size_t i = 0; i < vertices.size(); ++i) {
        cout << "Vertex " << i
             << " | Comp: " << vertices[i].comp
             << " | Type: " << vertices[i].type
             << " | Level: " << vertices[i].level
             << " | Lowlink: " << vertices[i].lowlink << endl;
    }
    

    return 0;
}
