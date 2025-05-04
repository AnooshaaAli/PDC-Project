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
#include <unordered_map>
#include <cmath>
#include <omp.h>
#include <climits>
#include <mpi.h>
#include <chrono>

#define MAX_TOPICS 10
#define DAMPING 0.85
#define EPSILON 0.0001
#define MAX_ITER 100

using namespace std;

struct DirectedEdge
{
    int from, to;
    int retweet, reply, mention;
    int partition;
};

struct Vertex
{
    int index = -1, lowlink = -1, level = -1, depth = -1;
    std::string type; // "scc" or "cac"
    int comp = -1;    // Component ID
};

struct ActionCount
{
    int RT = 0, RE = 0, MT = 0;
};

struct BFS_Tree
{
    unordered_map<int, int> parent;   // Key: Node, Value: Parent Node
    unordered_map<int, int> distance; // Key: Node, Value: Distance from root
};

int index_counter = 0, comp_counter = 0;
std::stack<int> s;
std::vector<Vertex> vertices;
std::vector<std::vector<int>> adj;

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

void dfs(int v);
void discover(int v);
void explore(int v);
void finish(int v);

vector<vector<DirectedEdge>> run_metis_partitioning(const string &filename, int nparts, set<int> partition_nodes[]);
vector<DirectedEdge> remap_partition_edges(const vector<DirectedEdge> &edges, const vector<int> &nodes, map<int, int> &node_to_index, map<int, int> &index_to_node);
void display_partition_results(const map<int, int> &index_to_node, const vector<Vertex> &vertices, const vector<DirectedEdge> &remapped_edges, const map<int, int> &node_to_index, int partition_id);

// Algorithm 1-4 and METIS functions

void discover(int v)
{
    vertices[v].index = index_counter;
    vertices[v].lowlink = index_counter;
    index_counter++;
    s.push(v);
}

void explore(int v)
{
    for (int w : adj[v])
    {
        if (vertices[w].index == -1)
        {
            dfs(w);
            vertices[v].lowlink = std::min(vertices[v].lowlink, vertices[w].lowlink);
        }
        else if (vertices[w].comp == -1)
        {
            vertices[v].lowlink = std::min(vertices[v].lowlink, vertices[w].lowlink);
        }
    }
}

void merge(int u, int v)
{
    vertices[v].comp = comp_counter;
}

void finish(int v)
{
    if (vertices[v].lowlink != vertices[v].index)
        return;

    int size = 0;
    int w;
    do
    {
        w = s.top();
        s.pop();
        vertices[w].lowlink = vertices[v].lowlink;
        vertices[w].type = "scc";
        merge(v, w);
        size++;
    } while (w != v);

    if (size == 1)
    {
        vertices[v].type = "cac";
    }

    comp_counter++;
}

void dfs(int v)
{
    discover(v);
    explore(v);
    finish(v);
}

void SCC_CAC_partition(const std::vector<DirectedEdge> &edges, int n)
{
    adj.assign(n, {});
    vertices.assign(n, Vertex());
    index_counter = 0;
    comp_counter = 0;
    while (!s.empty())
        s.pop();

    for (const auto &edge : edges)
    {
        adj[edge.from].push_back(edge.to);
    }

    for (int v = 0; v < n; ++v)
    {
        if (vertices[v].comp == -1)
        {
            dfs(v);
        }
    }

    vector<std::unordered_set<int>> comp_adj(comp_counter);
    std::vector<int> in_degree(comp_counter, 0);
    std::vector<int> comp_level(comp_counter, 0);

    for (int u = 0; u < n; ++u)
    {
        for (int v : adj[u])
        {
            int cu = vertices[u].comp;
            int cv = vertices[v].comp;
            if (cu != cv && !comp_adj[cu].count(cv))
            {
                comp_adj[cu].insert(cv);
                in_degree[cv]++;
            }
        }
    }

    queue<int> q;
    for (int i = 0; i < comp_counter; ++i)
    {
        if (in_degree[i] == 0)
        {
            q.push(i);
            comp_level[i] = 0;
        }
    }

    while (!q.empty())
    {
        int u = q.front();
        q.pop();
        for (int v : comp_adj[u])
        {
            comp_level[v] = max(comp_level[v], comp_level[u] + 1);
            if (--in_degree[v] == 0)
            {
                q.push(v);
            }
        }
    }

    for (int i = 0; i < n; ++i)
    {
        vertices[i].level = comp_level[vertices[i].comp];
    }
}

vector<vector<DirectedEdge>> run_metis_partitioning(const string &filename, int nparts, set<int> partition_nodes[])
{
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Failed to open file.\n";
        exit(1);
    }

    set<int> nodes;
    vector<tuple<int, int, int, int, int>> raw_edges;
    int u, v, retweet, reply, mention;

    while (file >> u >> v >> retweet >> reply >> mention)
    {
        raw_edges.emplace_back(u, v, retweet, reply, mention);
        nodes.insert(u);
        nodes.insert(v);
    }

    int nvtxs = *nodes.rbegin() + 1;

    vector<vector<int>> adjList(nvtxs);
    for (auto &[u, v, retweet, reply, mention] : raw_edges)
    {
        adjList[u].push_back(v);
        adjList[v].push_back(u); // for METIS
    }

    vector<idx_t> xadj(nvtxs + 1, 0);
    vector<idx_t> adjncy;
    for (int i = 0; i < nvtxs; ++i)
    {
        xadj[i + 1] = xadj[i] + adjList[i].size();
        for (int neighbor : adjList[i])
        {
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
        &objval, part.data());

    if (result != METIS_OK)
    {
        cerr << "METIS partitioning failed!\n";
        exit(1);
    }

    vector<vector<DirectedEdge>> partitioned_edges(nparts);
    for (auto &[u, v, retweet, reply, mention] : raw_edges)
    {
        if (part[u] == part[v])
        {
            DirectedEdge edge{u, v, retweet, reply, mention, part[u]};
            partitioned_edges[part[u]].push_back(edge);
            partition_nodes[part[u]].insert(u);
            partition_nodes[part[u]].insert(v);
        }
    }

    return partitioned_edges;
}

vector<DirectedEdge> remap_partition_edges(const vector<DirectedEdge> &edges, const set<int> &node_ids, map<int, int> &node_to_index_out, map<int, int> &index_to_node_out)
{
    int index = 0;

    for (int node : node_ids)
    {
        node_to_index_out[node] = index;
        index_to_node_out[index] = node;
        ++index;
    }

    vector<DirectedEdge> remapped_edges;
    for (const auto &edge : edges)
    {
        remapped_edges.push_back({node_to_index_out[edge.from],
                                  node_to_index_out[edge.to],
                                  edge.retweet,
                                  edge.reply,
                                  edge.mention,
                                  edge.partition});
    }

    return remapped_edges;
}

void display_partition_results(const map<int, int> &index_to_node, const vector<Vertex> &vertices, const vector<DirectedEdge> &remapped_edges, const map<int, int> &node_to_index, int partition_id)
{
    cout << "                   SCC/CAC for Partition " << partition_id << " \n";

    string base = "partition_" + to_string(partition_id) + "_";
    ofstream comp_out(base + "components.txt");
    ofstream level_out(base + "component_levels.txt");
    ofstream edge_out(base + "edges.txt");

    map<int, int> component_written; // comp_id → level

    for (int i = 0; i < (int)index_to_node.size(); ++i)
    {
        int real_node = index_to_node.at(i);

        if (vertices[i].comp != -1)
        {
            cout << "Vertex " << real_node
                 << " | Comp: " << vertices[i].comp
                 << " | Type: " << vertices[i].type
                 << " | Level: " << vertices[i].level
                 << " | Lowlink: " << vertices[i].lowlink << endl;

            comp_out << real_node << " " << vertices[i].comp << "\n";

            if (component_written.find(vertices[i].comp) == component_written.end())
            {
                component_written[vertices[i].comp] = vertices[i].level;
            }
        }
    }
    cout << "---------------------------------------------------------------------" << endl;

    for (const auto &[comp_id, level] : component_written)
    {
        level_out << comp_id << " " << level << "\n";
    }

    // Write edges (in original node ID form) to edges.txt
    for (const auto &edge : remapped_edges)
    {
        int from_real = index_to_node.at(edge.from);
        int to_real = index_to_node.at(edge.to);

        edge_out << from_real << " " << to_real << "\n";
    }

    comp_out.close();
    level_out.close();
    edge_out.close();
}

void write_activity_file(const vector<DirectedEdge>& remapped_edges, const map<int, int>& index_to_node, const string& filename) {
    ofstream f(filename);
    for (const auto& edge : remapped_edges) {
        int from = index_to_node.at(edge.from);
        int to = index_to_node.at(edge.to);
        for (int i = 0; i < edge.retweet; ++i)
            f << from << " " << to << " RT" << endl;
        for (int i = 0; i < edge.reply; ++i)
            f << from << " " << to << " RE" << endl;
        for (int i = 0; i < edge.mention; ++i)
            f << from << " " << to << " MT" << endl;
    }
    f.close();
}

// Algorithm 5-6 functions

double get_action_weight(const string &act)
{
    if (act == "RT")
        return 0.5;
    if (act == "RE")
        return 0.35;
    if (act == "MT")
        return 0.15;
    return 0.0;
}

double jaccard(const vector<int> &a, const vector<int> &b)
{
    int inter = 0, uni = 0;
    for (int i = 0; i < interest_len; i++)
    {
        if (a[i] || b[i])
            uni++;
        if (a[i] && b[i])
            inter++;
    }
    return (uni == 0) ? 0.0 : (double)inter / uni;
}

void load_interest(const string &file)
{
    ifstream f(file);
    string line;
    while (getline(f, line))
    {
        stringstream ss(line);
        int user, bit;
        ss >> user;
        vector<int> vec;
        while (ss >> bit)
            vec.push_back(bit);
        interest_map[user] = vec;
        interest_len = vec.size();
    }
}

void load_graph(const string &file)
{
    ifstream f(file);
    int u, v;
    while (f >> u >> v)
    {
        graph[u].push_back(v);
        followers[v].push_back(u);
        out_degree[u]++;
        in_degree[v]++;
    }
}

void load_activity(const string &file)
{
    ifstream f(file);
    int u, v;
    long ts;
    string act;
    while (f >> u >> v >> act)
    {
        if (act == "RT")
            action_map[u][v].RT++;
        if (act == "RE")
            action_map[u][v].RE++;
        if (act == "MT")
            action_map[u][v].MT++;
        total_posts[v]++;
    }
}

void load_component_map(const string &file)
{
    ifstream f(file);
    int node, comp;
    while (f >> node >> comp)
    {
        node_to_component[node] = comp;
        component_nodes[comp].push_back(node);
    }
}

void load_component_levels(const string &file)
{
    ifstream f(file);
    int comp, level;
    while (f >> comp >> level)
    {
        component_level[comp] = level;
        level_components[level].push_back(comp);
    }
}

double compute_psi(int u, int v)
{
    if (total_posts[v] == 0)
        return 0.0;
    if (!interest_map.count(u) || !interest_map.count(v))
        return 0.0;

    double sim = jaccard(interest_map[u], interest_map[v]);

    ActionCount &a = action_map[u][v];
    double weighted = get_action_weight("RT") * sim * a.RT + get_action_weight("RE") * sim * a.RE + get_action_weight("MT") * sim * a.MT;
    return weighted / total_posts[v];
}

void initialize_ip()
{
    set<int> all_nodes;

    for (const auto &[u, nbrs] : graph)
    {
        all_nodes.insert(u);
        for (int v : nbrs)
            all_nodes.insert(v);
    }

    for (const auto &[v, f_list] : followers)
    {
        all_nodes.insert(v);
        for (int f : f_list)
            all_nodes.insert(f);
    }

    num_nodes = all_nodes.size();

    for (int u : all_nodes)
    {
        IP[u] = 1.0 / num_nodes;
    }
}

void compute_influence_power(const vector<int> &nodes)
{
    int iter = 0;
    double delta;

    do
    {
        delta = 0.0;
        for (int u : nodes)
        {
            double sum = 0.0;
            for (int f : followers[u])
            {
                if (out_degree[f] > 0)
                {
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

void calculate_influence_power()
{
    vector<int> levels;
    for (const auto &[level, _] : level_components)
        levels.push_back(level);
    sort(levels.begin(), levels.end());

    for (int level : levels)
    {
        cout << "Processing Level " << level << "...\n";
        const auto &comps = level_components[level];

        #pragma omp parallel for
        for (int i = 0; i < comps.size(); ++i)
        {
            int thread_id = omp_get_thread_num();
            cout << "[OMP Thread " << thread_id << "] Processing Component " << comps[i] << endl;
            compute_influence_power(component_nodes[comps[i]]);
        }
    }
    cout << "---------------------------------------------------------------------" << endl;
}

void output_IP()
{
    for (const auto &[u, ip] : IP)
    {
        cout << "Node " << u << ": IP = " << ip << endl;
    }
    cout << "---------------------------------------------------------------------" << endl;
}

unordered_map<int, double> compute_IL(int v)
{
    unordered_map<int, double> IL;
    unordered_map<int, int> level_count;
    unordered_map<int, int> distance;
    queue<int> q;

    q.push(v);
    distance[v] = 0;

    while (!q.empty())
    {
        int u = q.front();
        q.pop();
        int lvl = distance[u];

        for (int nbr : followers[u])
        {
            if (!distance.count(nbr))
            {
                distance[nbr] = lvl + 1;
                q.push(nbr);
            }
        }
    }

    for (const auto &[node, lvl] : distance)
    {
        if (node != v)
        {
            IL[lvl] += IP[node];
            level_count[lvl]++;
        }
    }

    for (auto &[lvl, sum] : IL)
    {
        if (level_count[lvl] > 0)
            sum /= level_count[lvl];
    }

    return IL;
}

vector<pair <int, double>> find_seed_candidates()
{
    vector<pair <int, double>> seeds_scores;

    for (const auto &[v, _] : IP)
    {
        unordered_map<int, double> IL = compute_IL(v);

        vector<int> gamma;
        int max_lvl = 0;
        for (const auto &[lvl, _] : IL)
        {
            if (lvl > max_lvl)
                max_lvl = lvl;
        }

        for (int L = 1; L < max_lvl; ++L)
        {
            if (IL.count(L) && IL.count(L + 1))
            {
                if (IL[L] > IL[L + 1])
                {
                    gamma.push_back(L);
                }
            }
        }

        int L0 = gamma.empty() ? max_lvl : *min_element(gamma.begin(), gamma.end());

        double IL0 = IL.count(L0) ? IL[L0] : 0.0;

        if (IP[v] > IL0)
        {
            seeds_scores.emplace_back(v, IP[v]);
        }
    }

    return seeds_scores;
}

void clear_maps()
{
    level_components.clear();
    component_nodes.clear();
    component_level.clear();
    IP.clear();
    newIP.clear();
    total_posts.clear();
    psi.clear();
    graph.clear();
    followers.clear();
}

// Algorithm 7 functions
BFS_Tree influence_bfs_tree(int seed, const unordered_set<int> &seed_set)
{
    BFS_Tree tree;
    queue<int> q;

    // Start with the seed node, mark it as its own parent, distance is 0
    tree.parent[seed] = -1; // root has no parent
    tree.distance[seed] = 0;
    q.push(seed);

    while (!q.empty())
    {
        int current = q.front();
        q.pop();
        if (graph.find(current) == graph.end())
            continue;
        // Only explore neighbors that are seed nodes (black nodes)
        for (int neighbor : graph.at(current))
        {
            // If the neighbor is not yet visited
            if (tree.parent.find(neighbor) == tree.parent.end())
            {
                // Add the neighbor to the tree
                tree.parent[neighbor] = current;
                tree.distance[neighbor] = tree.distance[current] + 1;
                if (seed_set.count(neighbor))
                {
                    q.push(neighbor);
                }
            }
        }
    }

    return tree;
}

double compute_rank(const BFS_Tree &tree)
{
    if (tree.distance.empty())
        return 0.0;

    double total_level = 0.0;
    for (const auto &[node, level] : tree.distance)
    {
        total_level += level;
    }

    return total_level / tree.distance.size(); // average level
}

unordered_set<int> extract_blackpath(const BFS_Tree &tree, const unordered_set<int> &seed_set)
{
    unordered_set<int> blackpath;

    for (const auto &[node, _] : tree.distance)
    {
        if (seed_set.count(node))
        {
            blackpath.insert(node);
        }
    }

    return blackpath;
}

vector<pair<int, double>> seed_selection_algorithm(const vector<pair<int, double>> &seeds)
{
    unordered_set<int> seed_set;
    vector<pair<int, double>> final_seeds_scores;

    // Step 1: build all initial influence-BFS trees
    unordered_map<int, BFS_Tree> bfs_trees;
    for (const auto &[seed, _] : seeds)
    {
        seed_set.insert(seed);
        bfs_trees[seed] = influence_bfs_tree(seed, seed_set);
    }

    // Step 2: repeat until no seeds remain
    while (!seed_set.empty())
    {
        // Step 3: find seed that has largest influence tree
        int max_seed = -1;
        size_t max_size = 0;

        for (const auto &[seed, tree] : bfs_trees)
        {
            size_t size = tree.parent.size(); // total nodes in the tree
            if (size > max_size && seed_set.count(seed))
            {
                max_size = size;
                max_seed = seed;
            }
        }

        if (max_seed == -1)
            break; // no valid seeds left

        const BFS_Tree &largest_tree = bfs_trees[max_seed];

        // Step 4: get blackpath of this tree
        unordered_set<int> blackpath = extract_blackpath(largest_tree, seed_set);

        // Step 5: find tree rooted at each blackpath node and compute rank
        int min_rank_seed = -1;
        float min_rank = std::numeric_limits<float>::max();
        BFS_Tree best_tree;

        for (int black_node : blackpath)
        {
            BFS_Tree tree = bfs_trees[black_node];
            float rank = compute_rank(tree);

            cout << "\nFor black_node: " << black_node << ", found rank: " << rank << endl;

            if (rank < min_rank)
            {
                min_rank = rank;
                min_rank_seed = black_node;
                best_tree = std::move(tree);
            }
        }

        if (min_rank_seed != -1)
        {
            final_seeds_scores.emplace_back(min_rank_seed, IP[min_rank_seed]);

            // Remove all blackpath nodes of this best_tree from seed_set
            for (int node : blackpath)
            {
                seed_set.erase(node);
                bfs_trees.erase(node); // also remove its BFS tree if exists
            }
        }
    }

    return final_seeds_scores;
}

vector<pair <int, double>> select_best_k_seeds(vector<pair <int, double>> &final_seeds, int k)
{
    vector<pair<int, double>> seed_scores;
    for (const auto &[seed, score] : final_seeds)
    {
        seed_scores.emplace_back(seed, score);
    }

    // Sort seeds by their scores in descending order
    sort(seed_scores.begin(), seed_scores.end(),
         [](const pair<int, double> &a, const pair<int, double> &b) {
             return a.second > b.second;
         });

    // Select top k seeds
    vector<pair <int, double>> best_k_seeds;
    for (int i = 0; i < k && i < seed_scores.size(); ++i)
    {
        best_k_seeds.emplace_back(seed_scores[i].first, seed_scores[i].second);
    }

    return best_k_seeds;
}

vector<pair<int, double>> process_partition(int p, const vector<DirectedEdge>& remapped_edges, map<int, int>& node_to_index, map<int, int>& index_to_node)
{
    cout << "-------------------------- Partition " << p << " ------------------------------" << endl;
    SCC_CAC_partition(remapped_edges, node_to_index.size());
    display_partition_results(index_to_node, vertices, remapped_edges, node_to_index, p);

    string base = "partition_" + to_string(p) + "_";
    load_interest("interest.txt");
    load_graph(base + "edges.txt");
    load_activity("activity.txt");
    load_component_map(base + "components.txt");
    load_component_levels(base + "component_levels.txt");
    initialize_ip();
    calculate_influence_power();
    output_IP();

    vector<pair<int, double>> seeds = find_seed_candidates();
    cout << "\nSeed Candidates:\n";
    for (const auto &[s, score] : seeds) {
    cout << "Node " << s << " (IP = " << score << ")\n";
    }

    vector<pair<int, double>> final_seeds = seed_selection_algorithm(seeds);
    cout << "Final seeds:\n";
    for (const auto &[seed, _] : final_seeds) {
    cout << seed << " ";
    }
    cout << "\n---------------------------------------------------------------------\n\n";

    clear_maps();
    return final_seeds;
}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    cout << "Rank: " << rank << ", Total Processes: " << size << endl;

    const int nparts = size;
    int k = 3;

    if (size < 2) {
        cerr << "This program requires at least 2 processes." << endl;
        MPI_Finalize();
        return 1;
    }

    const string filename = "initial_dataset.txt";
    vector<pair<int, double>> all_final_seeds;
    auto start_time;

    if (rank == 0) {
        auto start_time = chrono::high_resolution_clock::now();
        // Master node: partition and distribute
        set<int> partition_nodes[nparts];
        vector<vector<DirectedEdge>> partitioned_edges =
            run_metis_partitioning(filename, nparts, partition_nodes);

        // Send data to slaves
        for (int p = 1; p < nparts; ++p) {
            int edge_count = partitioned_edges[p].size();
            MPI_Send(&edge_count, 1, MPI_INT, p, 0, MPI_COMM_WORLD);
            MPI_Send(partitioned_edges[p].data(), edge_count * sizeof(DirectedEdge), MPI_BYTE, p, 1, MPI_COMM_WORLD);

            int node_count = partition_nodes[p].size();
            vector<int> nodes(partition_nodes[p].begin(), partition_nodes[p].end());
            MPI_Send(&node_count, 1, MPI_INT, p, 2, MPI_COMM_WORLD);
            MPI_Send(nodes.data(), node_count, MPI_INT, p, 3, MPI_COMM_WORLD);
        }

        // Master processes its own partition
        set<int> node_set = partition_nodes[0];
        map<int, int> node_to_index, index_to_node;
        vector<DirectedEdge> remapped_edges =
            remap_partition_edges(partitioned_edges[0], node_set, node_to_index, index_to_node);

        if (!remapped_edges.empty()) {
            vector<pair<int, double>> final_seeds =
                process_partition(0, remapped_edges, node_to_index, index_to_node);
            all_final_seeds.insert(all_final_seeds.end(), final_seeds.begin(), final_seeds.end());
        }
        else {
            cout << "Partition 0 has no edges.\n";
        }

        // Receive final seeds from each slave
        for (int p = 1; p < nparts; ++p) {
            int seed_count;
            MPI_Recv(&seed_count, 1, MPI_INT, p, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            vector<pair<int, double>> seeds(seed_count);
            MPI_Recv(seeds.data(), seed_count * sizeof(pair<int, double>), MPI_BYTE, p, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            all_final_seeds.insert(all_final_seeds.end(), seeds.begin(), seeds.end());
        }

        // Final output
        cout << "---------------- Final Influencers User Selected --------------------" << endl;
        if (all_final_seeds.empty()) {
            cout << "No influencers found." << endl;
        } else {
            vector<pair<int, double>> best_k_seeds = select_best_k_seeds(all_final_seeds, k);
            for (const auto &[seed, score] : best_k_seeds) {
                cout << "User " << seed << " with IP " << score << endl;
            }
        }
        cout << "---------------------------------------------------------------------" << endl;
    } else {
        // SLAVE node
        int edge_count;
        MPI_Recv(&edge_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        vector<DirectedEdge> edges(edge_count);
        MPI_Recv(edges.data(), edge_count * sizeof(DirectedEdge), MPI_BYTE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int node_count;
        MPI_Recv(&node_count, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        vector<int> nodes(node_count);
        MPI_Recv(nodes.data(), node_count, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        set<int> node_set(nodes.begin(), nodes.end());

        map<int, int> node_to_index, index_to_node;
        vector<DirectedEdge> remapped_edges =
            remap_partition_edges(edges, node_set, node_to_index, index_to_node);

        vector<pair<int, double>> final_seeds;
        if (!remapped_edges.empty()) {
            final_seeds = process_partition(rank, remapped_edges, node_to_index, index_to_node);
        }

        int seed_count = final_seeds.size();
        MPI_Send(&seed_count, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);
        MPI_Send(final_seeds.data(), seed_count * sizeof(pair<int, double>), MPI_BYTE, 0, 5, MPI_COMM_WORLD);
    }
    
    if (rank == 0)
    {
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end_time - start_time;
        cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;
        cout << "---------------------------------------------------------------------" << endl;
        cout << "---------------------------- End of Program -------------------------" << endl;
    }

    MPI_Finalize();
    return 0;
}
