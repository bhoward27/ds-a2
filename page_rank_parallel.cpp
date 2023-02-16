#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <vector>
#include <functional>
#include <atomic>

using std::thread;
using std::vector;
using std::ref;
using std::atomic;

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
                    // TODO: Add back in         "/scratch/input_graphs/roadNet-CA")},
                    "/home/ben/Documents/SFU/cmpt-431/assignments/1/default_page_rank_graph/roadNet-CA")},
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
        break;
    case 3:
        std::cout << "\nDynamic task mapping\n";
        break;
    default:
        break;
    }

    return 0;
}
