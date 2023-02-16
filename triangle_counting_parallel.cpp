#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <vector>
#include <atomic>
#include <functional>
#include <utility>

using std::vector;
using std::thread;
using std::atomic;
using std::ref;
using std::pair;

typedef struct {
    long non_unique_triangles;
    long num_vertices;
    long num_edges;
    double time_taken;
} Result;

struct ThreadDistribution {
    long num_edges = 0;
    vector<uintV> vertices;
};

long countTriangles(uintV *array1, uintE len1, uintV *array2, uintE len2, intV u, uintV v)
{
    uintE i = 0, j = 0; // indexes for array1 and array2
    long count = 0;

    if (u == v)
        return count;

    while ((i < len1) && (j < len2)) {
        if (array1[i] == array2[j]) {
            if ((array1[i] != u) && (array1[i] != v)) {
                count++;
            } else {
                // triangle with self-referential edge -> ignore
            }
            i++;
            j++;
        } else if (array1[i] < array2[j]) {
            i++;
        } else {
            j++;
        }
    }
    return count;
}

void triangleCountSerial(Graph &g)
{
    uintV n = g.n_;
    long triangle_count = 0;
    double time_taken = 0.0;
    timer t1;
    t1.start();
    for (uintV u = 0; u < n; u++) {
        uintE out_degree = g.vertices_[u].getOutDegree();
        for (uintE i = 0; i < out_degree; i++) {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                             g.vertices_[u].getInDegree(),
                                             g.vertices_[v].getOutNeighbors(),
                                             g.vertices_[v].getOutDegree(),
                                             u,
                                             v);
        }
    }
    time_taken = t1.stop();
    std::cout << "Number of triangles : " << triangle_count << "\n";
    std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
}

void threadStrat1(Graph& g,
                 uintV first_vertex,
                 uintV final_vertex,
                 uint tid,
                 vector<Result>& thread_results,
                 atomic<long>& global_triangle_count)
{
    // Process each edge <u,v>
    timer t;
    t.start();
    long local_triangle_count = 0;
    long num_edges = 0;
    for (uintV u = first_vertex; u <= final_vertex; u++) {
        // For each outNeighbor v, find the intersection of inNeighbor(u) and
        // outNeighbor(v)
        uintE out_degree = g.vertices_[u].getOutDegree();
        num_edges += out_degree;
        for (uintE i = 0; i < out_degree; i++) {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            local_triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                                   g.vertices_[u].getInDegree(),
                                                   g.vertices_[v].getOutNeighbors(),
                                                   g.vertices_[v].getOutDegree(),
                                                   u,
                                                   v);
        }
    }

    // Atomically add to the global triangle count.
    long expected = global_triangle_count;
    long desired;
    do {
        desired = expected + local_triangle_count;
    } while(!global_triangle_count.compare_exchange_weak(expected, desired));

    Result res = {local_triangle_count, final_vertex - first_vertex + 1, num_edges, t.stop()};
    thread_results[tid] = res;
}

void triangleCountParallelStrat1(Graph &g, uint num_threads)
{
    uintV n = g.n_;
    atomic<long> triangle_count(0);
    double time_taken = 0.0;
    timer t1;
    vector<Result> thread_results;
    vector<thread> threads;

    double partitioning_time = 0.0;
    timer t2;
    t2.start();
    uintV quotient = n / num_threads;
    uintV remainder = n % num_threads;
    uintV step = quotient;
    uintV final_step = (remainder) ? step + remainder : step;
    uintV first_vertex = 0;
    uintV final_vertex = step - 1;
    partitioning_time += t2.stop();


    // The outNghs and inNghs for a given vertex are already sorted

    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    t1.start();
    for (uint tid = 0; tid < num_threads; tid++) {
        Result init_result = {0};
        thread_results.push_back(init_result);
        threads.push_back(
            thread(
                threadStrat1,
                ref(g),
                first_vertex,
                final_vertex,
                tid,
                ref(thread_results),
                ref(triangle_count)
            )
        );

        t2.start();
        first_vertex = final_vertex + 1;
        if (tid == num_threads - 2) {
            final_vertex += final_step;
        }
        else {
            final_vertex += step;
        }
        partitioning_time += t2.stop();
    }
    for (auto& lwp : threads) {
        lwp.join();
    }
    time_taken = t1.stop();
    // -------------------------------------------------------------------
    std::cout << "thread_id, num_vertices, num_edges, triangle_count, time_taken\n";
    for (uint tid = 0; tid < num_threads; tid++) {
        Result res = thread_results[tid];
        std::cout << tid << ", "
                  << res.num_vertices << ", "
                  << res.num_edges << ", "
                  << res.non_unique_triangles << ", "
                  << res.time_taken << "\n";
    }

    // Print the overall statistics
    std::cout << "Number of triangles : " << triangle_count << "\n";
    std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
    std::cout << "Partitioning time (in seconds) : " << std::setprecision(TIME_PRECISION) << partitioning_time << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
}

// TODO: Write this function to use strategy 2.
void threadStrat2(Graph& g,
                 const vector<long>& degrees,
                 const ThreadDistribution& thread_distro,
                 uint tid,
                 vector<Result>& thread_results,
                 atomic<long>& global_triangle_count)
{
    // Process each edge <u,v>
    timer t;
    t.start();
    long local_triangle_count = 0;
    for (uintV u : thread_distro.vertices) {
        // For each outNeighbor v, find the intersection of inNeighbor(u) and
        // outNeighbor(v)
        uintE out_degree = degrees[u];
        for (uintE i = 0; i < out_degree; i++) {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            local_triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                                   g.vertices_[u].getInDegree(),
                                                   g.vertices_[v].getOutNeighbors(),
                                                   g.vertices_[v].getOutDegree(),
                                                   u,
                                                   v);
        }
    }

    // Atomically add to the global triangle count.
    long expected = global_triangle_count;
    long desired;
    do {
        desired = expected + local_triangle_count;
    } while(!global_triangle_count.compare_exchange_weak(expected, desired));

    Result res = {
        local_triangle_count,
        static_cast<long>(thread_distro.vertices.size()),
        thread_distro.num_edges,
        t.stop()
    };
    thread_results[tid] = res;
}

int min(const vector<ThreadDistribution>& distros)
{
    // first = thread number who has min number of edges
    // sedond = min number of edges.
    pair<int, long> min_edges = {0, LONG_MAX};
    for (unsigned long t = 0; t < distros.size(); t++) {
        if (distros[t].num_edges < min_edges.second) {
            min_edges = {t, distros[t].num_edges};
        }
    }
    return min_edges.first;
}

// TODO: Write this function to use strategy 2.
void triangleCountParallelStrat2(Graph &g, uint num_threads)
{
    uintV n = g.n_;
    atomic<long> triangle_count(0);
    double time_taken = 0.0;
    timer t1;
    t1.start();
    vector<Result> thread_results;
    vector<thread> threads;

    // The outNghs and inNghs for a given vertex are already sorted

    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------

    // Distribute vertices to threads such that each thread processes roughly the same number of edges.
    // Step 1. Construct an associative array which holds the degree of each vertex.
    // We can use this array later on to save on processing, at the cost of some extra memory.
    timer t2;
    t2.start();
    vector<long> degrees(n);
    for (uintV v = 0; v < n; v++) {
        degrees[v] = g.vertices_[v].getOutDegree();
    }
    // Step 2. Hand out vertices to each thread, each time giving a vertex to whichever thread has the least
    // amount of edges currently allocated for it. Kind of a greedy strategy -- not sure if it's truly optimal, but
    // hopefully it's good enough for our purposes.

    vector<ThreadDistribution> thread_distros(num_threads, ThreadDistribution());
    // Assume cache block size is 64 bytes.
    // const int num_elems_per_block = 4 * 64 / sizeof(uintV);
    for (uintV v = 0; v < n; /* v += num_elems_per_block*/ v++) {
        int min_thread = min(thread_distros);
        // Assign the next two cache block's worth of vertices to min_thread to attempt to maintain a stride-1 access pattern
        // later in the program.
        // for (uintV u = v; (u < n) && (u < v + num_elems_per_block); u++) {
            // thread_distros[min_thread].vertices.push_back(u);
            // thread_distros[min_thread].num_edges += degrees[u];
        // }
            thread_distros[min_thread].vertices.push_back(v);
            thread_distros[min_thread].num_edges += degrees[v];
    }
    double partitioning_time = t2.stop();

    for (uint tid = 0; tid < num_threads; tid++) {
        Result init_result = {0};
        thread_results.push_back(init_result);
        threads.push_back(
            thread(
                threadStrat2,
                ref(g),
                ref(degrees),
                ref(thread_distros[tid]),
                tid,
                ref(thread_results),
                ref(triangle_count)
            )
        );
    }
    for (auto& lwp : threads) {
        lwp.join();
    }
    time_taken = t1.stop();
    // -------------------------------------------------------------------
    std::cout << "thread_id, num_vertices, num_edges, triangle_count, time_taken\n";
    for (uint tid = 0; tid < num_threads; tid++) {
        Result res = thread_results[tid];
        std::cout << tid << ", "
                  << 0 << ", "
                  << res.num_edges << ", "
                  << res.non_unique_triangles << ", "
                  << res.time_taken << "\n";
    }

    // Print the overall statistics
    std::cout << "Number of triangles : " << triangle_count << "\n";
    std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
    std::cout << "Partitioning time (in seconds) : " << std::setprecision(TIME_PRECISION) << partitioning_time << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
}

uintV get_next_vertex(uintV n)
{
    static atomic<uintV> next_vertex(0);

    if (next_vertex > n || next_vertex < 0) return -1;

    uintV expected = next_vertex;
    uintV desired;
    do {
        desired = expected + 1;
    } while(!next_vertex.compare_exchange_weak(expected, desired));

    return next_vertex;
}

void threadStrat3(Graph& g,
                 uint tid,
                 vector<Result>& thread_results,
                 atomic<long>& global_triangle_count,
                 double& global_partitioning_time)
{
    // Process each edge <u,v>
    timer t;
    t.start();
    long local_triangle_count = 0;
    long num_edges = 0;
    long num_vertices = 0;
    double partitioning_time = 0;
    timer t2;
    while (true) {
        // For each outNeighbor v, find the intersection of inNeighbor(u) and
        // outNeighbor(v)
        // NOTE: I consider this as part of the "partitioning" time. Wasn't too sure if that applies in dynamic task scenario,
        // so doing this just in case it does.
        if (tid == 0) t2.start();
        uintV u = get_next_vertex(g.n_);
        if (tid == 0) partitioning_time += t2.stop();
        if (u == -1) break;
        uintE out_degree = g.vertices_[u].getOutDegree();
        num_edges += out_degree;
        for (uintE i = 0; i < out_degree; i++) {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            local_triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                                   g.vertices_[u].getInDegree(),
                                                   g.vertices_[v].getOutNeighbors(),
                                                   g.vertices_[v].getOutDegree(),
                                                   u,
                                                   v);
        }
        num_vertices++;
    }

    // Atomically add to the global triangle count.
    long expected = global_triangle_count;
    long desired;
    do {
        desired = expected + local_triangle_count;
    } while(!global_triangle_count.compare_exchange_weak(expected, desired));

    if (tid == 0) global_partitioning_time = partitioning_time;

    Result res = {local_triangle_count, num_vertices, num_edges, t.stop()};
    thread_results[tid] = res;
}

void triangleCountParallelStrat3(Graph &g, uint num_threads)
{
    // uintV n = g.n_;
    atomic<long> triangle_count(0);
    double time_taken = 0.0;
    timer t1;
    vector<Result> thread_results;
    vector<thread> threads;

    double partitioning_time = 0.0;


    // The outNghs and inNghs for a given vertex are already sorted

    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    t1.start();
    for (uint tid = 0; tid < num_threads; tid++) {
        Result init_result = {0};
        thread_results.push_back(init_result);
        threads.push_back(
            thread(
                threadStrat3,
                ref(g),
                tid,
                ref(thread_results),
                ref(triangle_count),
                ref(partitioning_time)
            )
        );
    }
    for (auto& lwp : threads) {
        lwp.join();
    }
    time_taken = t1.stop();
    // -------------------------------------------------------------------
    std::cout << "thread_id, num_vertices, num_edges, triangle_count, time_taken\n";
    for (uint tid = 0; tid < num_threads; tid++) {
        Result res = thread_results[tid];
        std::cout << tid << ", "
                  << res.num_vertices << ", "
                  << res.num_edges << ", "
                  << res.non_unique_triangles << ", "
                  << res.time_taken << "\n";
    }

    // Print the overall statistics
    std::cout << "Number of triangles : " << triangle_count << "\n";
    std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
    std::cout << "Partitioning time (in seconds) : " << std::setprecision(TIME_PRECISION) << partitioning_time << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
}

int main(int argc, char *argv[])
{
    cxxopts::Options options(
            "triangle_counting_serial",
            "Count the number of triangles using serial and parallel execution");
    options.add_options(
            "custom",
            {
                    {"nWorkers", "Number of workers",
                     cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
                    {"strategy", "Strategy to be used",
                     cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                    {"inputFile", "Input graph file path",
                     cxxopts::value<std::string>()->default_value(
                    // TODO: Add back in         "/scratch/input_graphs/roadNet-CA")},
                    "/home/ben/Documents/SFU/cmpt-431/assignments/1/default_page_rank_graph/roadNet-CA")},
            });

    auto cl_options = options.parse(argc, argv);
    uint n_workers = cl_options["nWorkers"].as<uint>();
    uint strategy = cl_options["strategy"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();
    std::cout << std::fixed;
    std::cout << "Number of workers : " << n_workers << "\n";
    std::cout << "Task decomposition strategy : " << strategy << "\n";

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";

    switch (strategy) {
        case 0:
            std::cout << "\nSerial\n";
            triangleCountSerial(g);
            break;
        case 1:
            std::cout << "\nVertex-based work partitioning\n";
            triangleCountParallelStrat1(g, n_workers);
            break;
        case 2:
            std::cout << "\nEdge-based work partitioning\n";
            triangleCountParallelStrat2(g, n_workers);
            break;
        case 3:
            std::cout << "\nDynamic task mapping\n";
            triangleCountParallelStrat3(g, n_workers);
            break;
        default:
            break;
    }

    return 0;
}
