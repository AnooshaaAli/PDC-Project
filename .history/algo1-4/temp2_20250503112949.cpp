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
#include <queue>
#include <unordered_set>


using namespace std;

// struct DirectedEdge {
//     int from, to;
//     int retweet, reply, mention;
//     int partition;
// };

// struct Vertex {
//     int index = -1, lowlink = -1, level = -1, depth = -1;
//     std::string type; // "scc" or "cac"
//     int comp = -1;    // Component ID
// };

// struct Edge {
//     int from, to;
// };

// void dfs(int v, const std::vector<DirectedEdge>& edges);
// void discover(int v);
// void explore(int v, const std::vector<DirectedEdge>& edges);
// void finish(int v);

// int index_counter = 0, comp_counter = 0;
// std::stack<int> s;
// std::vector<Vertex> vertices;
// std::vector<std::vector<int>> adj;

// void discover(int v) {
//     vertices[v].index = index_counter;
//     vertices[v].lowlink = index_counter;
//     vertices[v].level = 1;
//     vertices[v].depth = 1;
//     index_counter++;
//     s.push(v);
// }

// void explore(int v, const std::vector<DirectedEdge>& edges) {
//     for (const auto& edge : edges) {
//         if (edge.from != v) continue;

//         int w = edge.to;
//         if (vertices[w].index == -1) {
//             dfs(w, edges);  // DFS wrapper that calls discover, explore, finish
//         }

//         if (!vertices[w].type.empty()) {
//             vertices[v].level = std::max(vertices[w].level, vertices[v].level + 1);
//         } else {
//             vertices[v].level = std::max(vertices[w].level, vertices[v].level);
//             vertices[v].lowlink = std::min(vertices[v].lowlink, vertices[w].lowlink);
//         }
//     }
// }

// void merge(int u, int v) {
//     vertices[v].comp = comp_counter;
// }

// void finish(int v, const std::vector<DirectedEdge>& edges) {
//     if (vertices[v].lowlink != vertices[v].index) return;

//     int size = 0;
//     int w;
//     do {
//         w = s.top(); s.pop();
//         vertices[w].lowlink = vertices[v].lowlink;
//         vertices[w].level = vertices[v].level;
//         vertices[w].type = "scc";
//         merge(v, w);
//         size++;
//     } while (w != v);

//     if (size == 1) {
//         vertices[v].type = "cac";
//         bool m = false;
//         std::vector<int> list;

//         // Use the edge list instead of adj[v]
//         for (const DirectedEdge& edge : edges) {
//             if (edge.from != v) continue;
//             int w = edge.to;

//             if (vertices[w].type == "scc" && vertices[w].level == vertices[v].level - 1) {
//                 m = false;
//                 break;
//             } else if (vertices[w].type == "cac" && vertices[w].level == vertices[v].level - 1) {
//                 list.push_back(w);
//                 m = true;
//             }
//         }

//         if (m) {
//             vertices[v].level -= 1;
//             for (int w : list) {
//                 merge(w, v);
//             }
//         }
//     }

//     comp_counter++;
// }


// void dfs(int v, const std::vector<DirectedEdge>& edges) {
//     discover(v);
//     explore(v, edges);
//     finish(v, edges);
// }

// void SCC_CAC_partition(const std::vector<std::vector<DirectedEdge>>& partitioned_edges) {
//     for (const auto& edges : partitioned_edges) {
//         for (const auto& edge : edges) {
//             if (vertices[edge.from].comp == -1) {
//                 dfs(edge.from, edges);
//             }
//         }
//     }
// }

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

int index_counter = 0, comp_counter = 0;
std::stack<int> s;
std::vector<Vertex> vertices;
std::vector<std::vector<int>> adj;

void dfs(int v);
void discover(int v);
void explore(int v);
void finish(int v);

void discover(int v) {
    vertices[v].index = index_counter;
    vertices[v].lowlink = index_counter;
    index_counter++;
    s.push(v);
}

void explore(int v) {
    for (int w : adj[v]) {
        if (vertices[w].index == -1) {
            dfs(w);
            vertices[v].lowlink = std::min(vertices[v].lowlink, vertices[w].lowlink);
        } else if (vertices[w].comp == -1) {
            vertices[v].lowlink = std::min(vertices[v].lowlink, vertices[w].lowlink);
        }
    }
}

void merge(int u, int v) {
    vertices[v].comp = comp_counter;
}

void finish(int v) {
    if (vertices[v].lowlink != vertices[v].index) return;

    int size = 0;
    int w;
    do {
        w = s.top(); s.pop();
        vertices[w].lowlink = vertices[v].lowlink;
        vertices[w].type = "scc";
        merge(v, w);
        size++;
    } while (w != v);

    if (size == 1) {
        vertices[v].type = "cac";
    }

    comp_counter++;
}

void dfs(int v) {
    discover(v);
    explore(v);
    finish(v);
}

void SCC_CAC_partition(const std::vector<DirectedEdge>& edges, int n) {
    adj.assign(n, {});
    vertices.assign(n, Vertex());
    index_counter = 0;
    comp_counter = 0;
    while (!s.empty()) s.pop();

    for (const auto& edge : edges) {
        adj[edge.from].push_back(edge.to);
    }

    for (int v = 0; v < n; ++v) {
        if (vertices[v].comp == -1) {
            dfs(v);
        }
    }

    vector<std::unordered_set<int>> comp_adj(comp_counter);
    std::vector<int> in_degree(comp_counter, 0);
    std::vector<int> comp_level(comp_counter, 0);

    for (int u = 0; u < n; ++u) {
        for (int v : adj[u]) {
            int cu = vertices[u].comp;
            int cv = vertices[v].comp;
            if (cu != cv && !comp_adj[cu].count(cv)) {
                comp_adj[cu].insert(cv);
                in_degree[cv]++;
            }
        }
    }

    queue<int> q;
    for (int i = 0; i < comp_counter; ++i) {
        if (in_degree[i] == 0) {
            q.push(i);
            comp_level[i] = 0;
        }
    }

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : comp_adj[u]) {
            comp_level[v] = max(comp_level[v], comp_level[u] + 1);
            if (--in_degree[v] == 0) {
                q.push(v);
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        vertices[i].level = comp_level[vertices[i].comp];
    }
}



int main() {
    string filename = "data3.txt";
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open file.\n";
        return 1;
    }

    set<int> nodes;
    vector<tuple<int, int, int, int, int>> raw_edges;

    int u, v, retweet, reply, mention;
    while (file >> u >> v >> retweet >> reply >> mention) {
        raw_edges.emplace_back(u, v, retweet, reply, mention);
        nodes.insert(u);
        nodes.insert(v);
    }

    int nvtxs = *nodes.rbegin() + 1;

    vector<vector<int>> adjList(nvtxs);
    for (auto& [u, v, retweet, reply, mention] : raw_edges) {
        adjList[u].push_back(v);
        adjList[v].push_back(u); // for METIS only
    }

    vector<idx_t> xadj(nvtxs + 1, 0);
    vector<idx_t> adjncy;
    for (int i = 0; i < nvtxs; ++i) {
        xadj[i + 1] = xadj[i] + adjList[i].size();
        for (auto& neighbor : adjList[i]) {
            adjncy.push_back(neighbor);
        }
    }

    vector<idx_t> vwgt(nvtxs, 1);
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

    if (result != METIS_OK) {
        cerr << "METIS partitioning failed!\n";
        return 1;
    }

    vector<vector<DirectedEdge>> partitioned_edges(nparts);
    set<int> partition_nodes[nparts];

    for (auto& [u, v, retweet, reply, mention] : raw_edges) {
        if (part[u] == part[v]) {
            DirectedEdge edge{u, v, retweet, reply, mention, part[u]};
            partitioned_edges[part[u]].push_back(edge);
            partition_nodes[part[u]].insert(u);
            partition_nodes[part[u]].insert(v);
        }
    }

    for (int p = 0; p < nparts; ++p) {
        cout << "\n--- SCC/CAC for Partition " << p << " ---\n";

        // Map real node IDs to 0...n-1 for this partition
        map<int, int> node_to_index;
        map<int, int> index_to_node;
        int index = 0;
        for (int node : partition_nodes[p]) {
            node_to_index[node] = index;
            index_to_node[index] = node;
            index++;
        }

        vector<DirectedEdge> remapped_edges;
        for (auto& edge : partitioned_edges[p]) {
            remapped_edges.push_back({
                node_to_index[edge.from],
                node_to_index[edge.to],
                edge.retweet,
                edge.reply,
                edge.mention,
                edge.partition
            });
        }

        // Run SCC + CAC
        SCC_CAC_partition(remapped_edges, node_to_index.size());

        // Display results using real node IDs
    for (int i = 0; i < (int)node_to_index.size(); ++i) {
        int real_node = index_to_node[i];
        if (vertices[i].comp != -1) {
            cout << "Vertex " << real_node
                << " | Comp: " << vertices[i].comp
                << " | Type: " << vertices[i].type
                << " | Level: " << vertices[i].level
                << " | Lowlink: " << vertices[i].lowlink << endl;

            // Write components and component_levels to files
            string base = "partition_" + to_string(p) + "_";

            // Opening files to write the components and levels
            ofstream comp_out(base + "components.txt");
            ofstream level_out(base + "component_levels.txt");

            // Avoiding duplicates by keeping track of components already written
            map<int, int> component_written; // comp_id → level
            for (int i = 0; i < (int)node_to_index.size(); ++i) {
                int real_node = index_to_node[i];
                int comp_id = vertices[i].comp;
                comp_out << real_node << " " << comp_id << "\n";  // Writing component ID
                if (component_written.find(comp_id) == component_written.end()) {
                    component_written[comp_id] = vertices[i].level;
                }
            }

            // Writing component levels
            for (auto& [comp_id, level] : component_written) {
                level_out << comp_id << " " << level << "\n";  // Writing component level
            }

            // Closing the files
            comp_out.close();
            level_out.close();
        }
    }
    }


    return 0;
}

