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
vector<vector<DirectedEdge>> run_metis_partitioning(
    const string& filename,
    int nparts,
    set<int> partition_nodes[]
);
vector<DirectedEdge> remap_partition_edges(
    const vector<DirectedEdge>& edges,
    const vector<int>& nodes,
    map<int, int>& node_to_index,
    map<int, int>& index_to_node
);
void display_partition_results(
    const map<int, int>& index_to_node,
    const vector<Vertex>& vertices,
    const vector<DirectedEdge>& remapped_edges,
    const map<int, int>& node_to_index,
    int partition_id
);

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

void SCC_CAC_partition(const std::vector<DirectedEdge>& edges, int n){
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

vector<vector<DirectedEdge>> run_metis_partitioning(
    const string& filename,
    int nparts,
    set<int> partition_nodes[]
) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open file.\n";
        exit(1);
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
        adjList[v].push_back(u); // for METIS
    }

    vector<idx_t> xadj(nvtxs + 1, 0);
    vector<idx_t> adjncy;
    for (int i = 0; i < nvtxs; ++i) {
        xadj[i + 1] = xadj[i] + adjList[i].size();
        for (int neighbor : adjList[i]) {
            adjncy.push_back(neighbor);
        }
    }

    vector<idx_t> vwgt(nvtxs, 1);
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
        exit(1);
    }

    vector<vector<DirectedEdge>> partitioned_edges(nparts);
    for (auto& [u, v, retweet, reply, mention] : raw_edges) {
        if (part[u] == part[v]) {
            DirectedEdge edge{u, v, retweet, reply, mention, part[u]};
            partitioned_edges[part[u]].push_back(edge);
            partition_nodes[part[u]].insert(u);
            partition_nodes[part[u]].insert(v);
        }
    }

    return partitioned_edges;
}

vector<DirectedEdge> remap_partition_edges(
    const vector<DirectedEdge>& edges,
    const set<int>& node_ids,
    map<int, int>& node_to_index_out,
    map<int, int>& index_to_node_out
) {
    int index = 0;

    for (int node : node_ids) {
        node_to_index_out[node] = index;
        index_to_node_out[index] = node;
        ++index;
    }

    vector<DirectedEdge> remapped_edges;
    for (const auto& edge : edges) {
        remapped_edges.push_back({
            node_to_index_out[edge.from],
            node_to_index_out[edge.to],
            edge.retweet,
            edge.reply,
            edge.mention,
            edge.partition
        });
    }

    return remapped_edges;
}


void display_partition_results(
    const map<int, int>& index_to_node,
    const vector<Vertex>& vertices,
    const vector<DirectedEdge>& remapped_edges,
    const map<int, int>& node_to_index,
    int partition_id
) {
    cout << "\n--- SCC/CAC for Partition " << partition_id << " ---\n";

    string base = "partition_" + to_string(partition_id) + "_";
    ofstream comp_out(base + "components.txt");
    ofstream level_out(base + "component_levels.txt");
    ofstream edge_out(base + "edges.txt");

    map<int, int> component_written; // comp_id → level

    for (int i = 0; i < (int)index_to_node.size(); ++i) {
        int real_node = index_to_node.at(i);

        if (vertices[i].comp != -1) {
            cout << "Vertex " << real_node
                 << " | Comp: " << vertices[i].comp
                 << " | Type: " << vertices[i].type
                 << " | Level: " << vertices[i].level
                 << " | Lowlink: " << vertices[i].lowlink << endl;

            comp_out << real_node << " " << vertices[i].comp << "\n";

            if (component_written.find(vertices[i].comp) == component_written.end()) {
                component_written[vertices[i].comp] = vertices[i].level;
            }
        }
    }

    for (const auto& [comp_id, level] : component_written) {
        level_out << comp_id << " " << level << "\n";
    }

    // Write edges (in original node ID form) to edges.txt
    for (const auto& edge : remapped_edges) {
        int from_real = index_to_node.at(edge.from);
        int to_real = index_to_node.at(edge.to);

        edge_out << from_real << " " << to_real << "\n";
    }

    comp_out.close();
    level_out.close();
    edge_out.close();
}


int main() {
    const string filename = "data3.txt";
    const int nparts = 2;

    set<int> partition_nodes[nparts];
    vector<vector<DirectedEdge>> partitioned_edges =
        run_metis_partitioning(filename, nparts, partition_nodes);

        for (int p = 0; p < nparts; ++p)
        {
            cout << "\n--- SCC/CAC for Partition " << p << " ---\n";
        
            // Maps between real node IDs and 0...n-1 indices
            map<int, int> node_to_index;
            map<int, int> index_to_node;
        
            // Remap edges and populate maps
            vector<DirectedEdge> remapped_edges =
                remap_partition_edges(partitioned_edges[p], partition_nodes[p], node_to_index, index_to_node);
        
            // Run SCC + CAC
            SCC_CAC_partition(remapped_edges, node_to_index.size());
        
            // Display results using real node IDs
            display_partition_results(index_to_node, vertices, remapped_edges, node_to_index, p);
        }
        


    return 0;
}

