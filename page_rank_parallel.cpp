#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <vector>
#include <functional>
#include <atomic>
#include <utility>

using std::thread;
using std::vector;
using std::ref;
using std::atomic;
using std::pair;

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef float PageRankType;
#endif

typedef struct {
    long num_vertices;
    long num_edges;
    double barrier1_time;
    double barrier2_time;
    double getNextVertex_time;
    double time_taken;
} Result;

struct ThreadDistribution {
    long num_edges = 0;
    vector<uintV> vertices;
};

void pageRankSerial(Graph &g, int max_iters) {
    uintV n = g.n_;

    PageRankType *pr_curr = new PageRankType[n];
    PageRankType *pr_next = new PageRankType[n];

    for (uintV i = 0; i < n; i++) {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    // Push based pagerank
    timer t1;
    double time_taken = 0.0;
    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    t1.start();
    for (int iter = 0; iter < max_iters; iter++) {
        // for each vertex 'u', process all its outNeighbors 'v'
        for (uintV u = 0; u < n; u++) {
            uintE out_degree = g.vertices_[u].getOutDegree();
            for (uintE i = 0; i < out_degree; i++) {
                uintV v = g.vertices_[u].getOutNeighbor(i);
                pr_next[v] += (pr_curr[u] / out_degree);
            }
        }
        for (uintV v = 0; v < n; v++) {
            pr_next[v] = PAGE_RANK(pr_next[v]);

            // reset pr_curr for the next iteration
            pr_curr[v] = pr_next[v];
            pr_next[v] = 0.0;
        }
    }
    time_taken = t1.stop();
    // -------------------------------------------------------------------

    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++) {
        sum_of_page_ranks += pr_curr[u];
    }
    std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
    std::cout << "Time taken (in seconds) : " << time_taken << "\n";
    delete[] pr_curr;
    delete[] pr_next;
}

void threadStrat1(int tid,
                  vector<Result>& thread_results,
                  Graph& g,
                  const uintV first_vertex,
                  const uintV final_vertex,
                  int max_iters,
                  PageRankType* pr_curr,
                  atomic<PageRankType>* pr_next,
                  CustomBarrier& barrier)
{
    timer barrier1_timer, barrier2_timer, total_timer;
    double barrier1_time = 0;
    double barrier2_time =0;
    long num_edges = 0;
    total_timer.start();
    for (int iter = 0; iter < max_iters; iter++) {
        for (uintV u = first_vertex; u <= final_vertex; u++) {
            uintE out_degree = g.vertices_[u].getOutDegree();
            num_edges += out_degree;
            for (uintE i = 0; i < out_degree; i++) {
                uintV v = g.vertices_[u].getOutNeighbor(i);
                PageRankType expected = pr_next[v];
                PageRankType desired;
                do {
                    desired = expected + pr_curr[u] / out_degree;
                } while(!pr_next[v].compare_exchange_weak(expected, desired));
            }
        }
        barrier1_timer.start();
        barrier.wait();
        barrier1_time += barrier1_timer.stop();
        for (uintV v = first_vertex; v <= final_vertex; v++) {
            // No lock needed here, since v is only from this thread's subset of vertices.
            pr_next[v] = PAGE_RANK(pr_next[v]);
            pr_curr[v] = pr_next[v];
            pr_next[v] = 0.0;
        }
        barrier2_timer.start();
        barrier.wait();
        barrier2_time += barrier2_timer.stop();
    }
    long num_vertices = (final_vertex - first_vertex + 1) * max_iters;
    thread_results[tid] = {num_vertices, num_edges, barrier1_time, barrier2_time, 0, total_timer.stop()};
}

void pageRankParallelStrat1(Graph &g, int max_iters, int num_threads) {
    uintV n = g.n_;

    PageRankType* pr_curr = new PageRankType[n];
    atomic<PageRankType>* pr_next = new atomic<PageRankType>[n];

    for (uintV i = 0; i < n; i++) {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    // Push based pagerank
    timer t1, partition_timer;
    double time_taken = 0.0;
    double partitioning_time = 0;
    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    vector<thread> threads;
    vector<Result> thread_results(num_threads);
    CustomBarrier barrier(num_threads);

    t1.start();
    partition_timer.start();
    uintV quotient = n / num_threads;
    uintV remainder = n % num_threads;
    uintV step = quotient;
    uintV final_step = (remainder) ? step + remainder : step;
    uintV first_vertex = 0;
    uintV final_vertex = step - 1;
    partitioning_time += partition_timer.stop();
    for (int tid = 0; tid < num_threads; tid++) {
        threads.push_back(thread(threadStrat1,
                                 tid,
                                 ref(thread_results),
                                 ref(g),
                                 first_vertex,
                                 final_vertex,
                                 max_iters,
                                 pr_curr,
                                 pr_next,
                                 ref(barrier)));

        partition_timer.start();
        first_vertex = final_vertex + 1;
        if (tid == num_threads - 2) {
            final_vertex += final_step;
        }
        else {
            final_vertex += step;
        }
        partitioning_time += partition_timer.stop();
    }
    for (auto& thread : threads) {
        thread.join();
    }
    time_taken = t1.stop();

    std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, getNextVertex_time, total_time\n";
    for (int tid = 0; tid < num_threads; tid++) {
        std::cout << tid << ", "
                  << thread_results[tid].num_vertices << ", "
                  << thread_results[tid].num_edges << ", "
                  << thread_results[tid].barrier1_time << ", "
                  << thread_results[tid].barrier2_time << ", "
                  << thread_results[tid].getNextVertex_time << ", "
                  << thread_results[tid].time_taken << "\n";
    }

    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++) {
        sum_of_page_ranks += pr_curr[u];
    }
    std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
    std::cout << "Partitioning time (in seconds) : " << std::setprecision(TIME_PRECISION) << partitioning_time << "\n";
    std::cout << "Time taken (in seconds) : " << time_taken << "\n";
    delete[] pr_curr;
    delete[] pr_next;
}

int min(const vector<ThreadDistribution>& distros)
{
    // first = thread number who has min number of edges
    // second = min number of edges.
    pair<int, long> min_edges = {0, LONG_MAX};
    for (unsigned long t = 0; t < distros.size(); t++) {
        if (distros[t].num_edges < min_edges.second) {
            min_edges = {t, distros[t].num_edges};
        }
    }
    return min_edges.first;
}
void threadStrat2(int tid,
                  vector<Result>& thread_results,
                  Graph& g,
                  const vector<long>& degrees,
                  const ThreadDistribution& thread_distro,
                  int max_iters,
                  PageRankType* pr_curr,
                  atomic<PageRankType>* pr_next,
                  CustomBarrier& barrier)
{
    timer barrier1_timer, barrier2_timer, total_timer;
    double barrier1_time = 0;
    double barrier2_time = 0;
    total_timer.start();
    for (int iter = 0; iter < max_iters; iter++) {
        for (uintV u : thread_distro.vertices) {
            uintE out_degree = degrees[u];
            for (uintE i = 0; i < out_degree; i++) {
                uintV v = g.vertices_[u].getOutNeighbor(i);
                PageRankType expected = pr_next[v];
                PageRankType desired;
                do {
                    desired = expected + pr_curr[u] / out_degree;
                } while(!pr_next[v].compare_exchange_weak(expected, desired));
            }
        }
        barrier1_timer.start();
        barrier.wait();
        barrier1_time += barrier1_timer.stop();
        for (uintV v : thread_distro.vertices) {
            // No lock needed here, since v is only from this thread's subset of vertices.
            pr_next[v] = PAGE_RANK(pr_next[v]);
            pr_curr[v] = pr_next[v];
            pr_next[v] = 0.0;
        }
        barrier2_timer.start();
        barrier.wait();
        barrier2_time += barrier2_timer.stop();
    }
    long num_vertices = thread_distro.vertices.size() * max_iters;
    long num_edges = thread_distro.num_edges * max_iters;
    thread_results[tid] = {num_vertices, num_edges, barrier1_time, barrier2_time, 0, total_timer.stop()};
}


void pageRankParallelStrat2(Graph &g, int max_iters, int num_threads) {
    uintV n = g.n_;

    PageRankType* pr_curr = new PageRankType[n];
    atomic<PageRankType>* pr_next = new atomic<PageRankType>[n];

    for (uintV i = 0; i < n; i++) {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    // Push based pagerank
    timer t1, partition_timer;
    double time_taken = 0.0;
    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    vector<thread> threads;
    vector<Result> thread_results(num_threads);
    CustomBarrier barrier(num_threads);

    t1.start();

    partition_timer.start();
    vector<long> degrees(n);
    for (uintV v = 0; v < n; v++) {
        degrees[v] = g.vertices_[v].getOutDegree();
    }
    vector<ThreadDistribution> thread_distros(num_threads, ThreadDistribution());
    for (uintV v = 0; v < n; v++) {
        int min_thread = min(thread_distros);
        thread_distros[min_thread].vertices.push_back(v);
        thread_distros[min_thread].num_edges += degrees[v];
    }
    double partitioning_time = partition_timer.stop();

    for (int tid = 0; tid < num_threads; tid++) {
        threads.push_back(thread(threadStrat2,
                                 tid,
                                 ref(thread_results),
                                 ref(g),
                                 ref(degrees),
                                 ref(thread_distros[tid]),
                                 max_iters,
                                 pr_curr,
                                 pr_next,
                                 ref(barrier)));
    }
    for (auto& thread : threads) {
        thread.join();
    }

    time_taken = t1.stop();

    std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, getNextVertex_time, total_time\n";
    for (int tid = 0; tid < num_threads; tid++) {
        std::cout << tid << ", "
                  << thread_results[tid].num_vertices << ", "
                  << thread_results[tid].num_edges << ", "
                  << thread_results[tid].barrier1_time << ", "
                  << thread_results[tid].barrier2_time << ", "
                  << thread_results[tid].getNextVertex_time << ", "
                  << thread_results[tid].time_taken << "\n";
    }

    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++) {
        sum_of_page_ranks += pr_curr[u];
    }
    std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
    std::cout << "Partitioning time (in seconds) : " << std::setprecision(TIME_PRECISION) << partitioning_time << "\n";
    std::cout << "Time taken (in seconds) : " << time_taken << "\n";
    delete[] pr_curr;
    delete[] pr_next;
}

uintV get_next_vertex(atomic<uintV>& next_vertex, uintV n, uintV k)
{
    if (next_vertex >= n || next_vertex < -k) return -1;

    uintV expected = next_vertex;
    uintV desired;
    do {
        desired = expected + k;
    } while(!next_vertex.compare_exchange_weak(expected, desired));

    if (desired >= n) return -1;

    return desired;
}

void reset_next_vertex(atomic<uintV>& next_vertex, const uintV k)
{
    next_vertex = -k;
}

void threadStrat3(int tid,
                  vector<Result>& thread_results,
                  Graph& g,
                  atomic<uintV>& next_vertex_first,
                  atomic<uintV>& next_vertex_second,
                  const int max_iters,
                  const uintV k,
                  PageRankType* pr_curr,
                  atomic<PageRankType>* pr_next,
                  CustomBarrier& barrier)
{
    timer barrier1_timer, barrier2_timer, total_timer, get_next_vertex_timer;
    double barrier1_time = 0;
    double barrier2_time = 0;
    double get_next_vertex_time = 0;
    long num_vertices = 0;
    long num_edges = 0;
    total_timer.start();
    for (int iter = 0; iter < max_iters; iter++) {
        while (true) {
            get_next_vertex_timer.start();
            uintV u = get_next_vertex(next_vertex_first, g.n_, k);
            get_next_vertex_time += get_next_vertex_timer.stop();
            if (u == -1) {
                break;
            }
            for (uintV i = 0; i < k && u < g.n_; i++, u++) {
                uintE out_degree = g.vertices_[u].getOutDegree();
                for (uintE i = 0; i < out_degree; i++) {
                    uintV v = g.vertices_[u].getOutNeighbor(i);
                    PageRankType expected = pr_next[v];
                    PageRankType desired;
                    do {
                        desired = expected + pr_curr[u] / out_degree;
                    } while(!pr_next[v].compare_exchange_weak(expected, desired));
                }
                num_edges += out_degree;
            }
        }
        barrier1_timer.start();
        barrier.wait();
        barrier1_time += barrier1_timer.stop();
        while (true) {
            get_next_vertex_timer.start();
            uintV v = get_next_vertex(next_vertex_second, g.n_, k);
            get_next_vertex_time += get_next_vertex_timer.stop();
            if (v == -1) {
                break;
            }
            for (uintV i = 0; i < k && v < g.n_; i++, v++) {
                // No lock needed here, since v is only from this thread's subset of vertices.
                pr_next[v] = PAGE_RANK(pr_next[v]);
                pr_curr[v] = pr_next[v];
                pr_next[v] = 0.0;
                num_vertices++;
            }
        }
        barrier2_timer.start();
        barrier.wait();
        get_next_vertex_timer.start();
        reset_next_vertex(next_vertex_first, k);
        reset_next_vertex(next_vertex_second, k);
        get_next_vertex_time += get_next_vertex_timer.stop();
        barrier.wait();
        barrier2_time += barrier2_timer.stop();
    }

    thread_results[tid] = {num_vertices, num_edges, barrier1_time, barrier2_time, get_next_vertex_time, total_timer.stop()};
}

void pageRankParallelStrat3(Graph &g, int max_iters, int num_threads, uintV granularity) {
    uintV n = g.n_;

    PageRankType* pr_curr = new PageRankType[n];
    atomic<PageRankType>* pr_next = new atomic<PageRankType>[n];

    for (uintV i = 0; i < n; i++) {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    // Push based pagerank
    timer t1, partition_timer;
    double time_taken = 0.0;
    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    vector<thread> threads;
    vector<Result> thread_results(num_threads);
    CustomBarrier barrier(num_threads);
    atomic<uintV> next_vertex_first(-granularity);
    atomic<uintV> next_vertex_second(-granularity);

    t1.start();

    for (int tid = 0; tid < num_threads; tid++) {
        threads.push_back(thread(threadStrat3,
                                 tid,
                                 ref(thread_results),
                                 ref(g),
                                 ref(next_vertex_first),
                                 ref(next_vertex_second),
                                 max_iters,
                                 granularity,
                                 pr_curr,
                                 pr_next,
                                 ref(barrier)));
    }
    for (auto& thread : threads) {
        thread.join();
    }

    time_taken = t1.stop();

    std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, getNextVertex_time, total_time\n";
    for (int tid = 0; tid < num_threads; tid++) {
        std::cout << tid << ", "
                  << thread_results[tid].num_vertices << ", "
                  << thread_results[tid].num_edges << ", "
                  << thread_results[tid].barrier1_time << ", "
                  << thread_results[tid].barrier2_time << ", "
                  << thread_results[tid].getNextVertex_time << ", "
                  << thread_results[tid].time_taken << "\n";
    }

    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++) {
        sum_of_page_ranks += pr_curr[u];
    }
    std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
    std::cout << "Partitioning time (in seconds) : " << 0.0 << "\n";
    std::cout << "Time taken (in seconds) : " << time_taken << "\n";
    delete[] pr_curr;
    delete[] pr_next;
    pr_curr = nullptr;
    pr_next = nullptr;
}

int main(int argc, char *argv[]) {
    cxxopts::Options options(
            "page_rank_push",
            "Calculate page_rank using serial and parallel execution");
    options.add_options(
            "",
            {
                    {"nWorkers", "Number of workers",
                     cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
                    {"nIterations", "Maximum number of iterations",
                     cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
                    {"strategy", "Strategy to be used",
                     cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                    {"granularity", "Granularity to be used",
                     cxxopts::value<uint>()->default_value(DEFAULT_GRANULARITY)},
                    {"inputFile", "Input graph file path",
                     cxxopts::value<std::string>()->default_value(
                    "/scratch/input_graphs/roadNet-CA")},
            });

    auto cl_options = options.parse(argc, argv);
    uint n_workers = cl_options["nWorkers"].as<uint>();
    uint strategy = cl_options["strategy"].as<uint>();
    uint max_iterations = cl_options["nIterations"].as<uint>();
    uint granularity = cl_options["granularity"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
    std::cout << "Using INT\n";
#else
    std::cout << "Using FLOAT\n";
#endif
    std::cout << std::fixed;
    std::cout << "Number of workers : " << n_workers << "\n";
    std::cout << "Task decomposition strategy : " << strategy << "\n";
    std::cout << "Iterations : " << max_iterations << "\n";
    std::cout << "Granularity : " << granularity << "\n";

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";
    switch (strategy) {
        case 0:
            std::cout << "\nSerial\n";
            pageRankSerial(g, max_iterations);
            break;
        case 1:
            std::cout << "\nVertex-based work partitioning\n";
            pageRankParallelStrat1(g, max_iterations, n_workers);
            break;
        case 2:
            std::cout << "\nEdge-based work partitioning\n";
            pageRankParallelStrat2(g, max_iterations, n_workers);
            break;
        case 3:
            std::cout << "\nDynamic task mapping\n";
            pageRankParallelStrat3(g, max_iterations, n_workers, granularity);
            break;
        default:
            break;
    }

    return 0;
}
