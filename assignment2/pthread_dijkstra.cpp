/* Name: Junxue ZHANG
 * ID: 20371613
 * Email: jzhangcs@ust.hk
 */

/*
 * This is code skeleton for COMP5112-17Spring assignment2
 * Compile: g++ -std=c++11 -lpthread -o pthread_dijkstra pthread_dijkstra.cpp
 * Run: ./pthread_dijkstra -n <number of threads> -i <input file>,
 * you will find the output in 'output.txt' file
 */

#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <climits>
#include <cstring>
#include <algorithm>
#include <sys/time.h>
#include <time.h>
#include <getopt.h>

#include <pthread.h>

using std::string;
using std::cout;
using std::endl;
using std::vector;

// #define MY_DEBUG

/*
 * utils is a namespace for utility functions
 * including I/O (read input file and print results) and one matrix dimension convert(2D->1D) function
 */
namespace utils {
    int num_threads; //number of thread
    int N; //number of vertices
    int *mat; // the adjacency matrix

    string filename; // input file name
    string outputfile; //output file name, default: 'output.txt'

    void print_usage() {
        cout << "Usage:\n" << "\tpthread_dijkstra -n <number of threads> -i <input file>" << endl;
        exit(0);
    }

    int parse_args(int argc, char **argv) {
        filename = "";
        outputfile = "output.txt";
        num_threads = 0;

        int opt;
        if (argc < 2) {
            print_usage();
        }
        while ((opt = getopt(argc, argv, "n:i:o:h")) != EOF) {
            switch (opt) {
                case 'n':
                    num_threads = atoi(optarg);
                    break;
                case 'i':
                    filename = optarg;
                    break;
                case 'o':
                    outputfile = optarg;
                    break;
                case 'h':
                case '?':
                default:
                    print_usage();
            }
        }
        if (filename.length() == 0 || num_threads == 0)
            print_usage();
        return 0;
    }

    /*
     * convert 2-dimension coordinate to 1-dimension
     */
    int convert_dimension_2D_1D(int x, int y) {
        return x * N + y;
    }

    int read_file(string filename) {
        std::ifstream inputf(filename, std::ifstream::in);
        inputf >> N;
        assert(N < (1024 * 1024 *
                    20)); // input matrix should be smaller than 20MB * 20MB (400MB, we don't have two much memory for multi-processors)
        mat = (int *) malloc(N * N * sizeof(int));
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                inputf >> mat[convert_dimension_2D_1D(i, j)];
            }

        return 0;
    }

    string format_path(int i, int *pred) {
        string out("");
        int current_vertex = i;
        while (current_vertex != 0) {
            string s = std::to_string(current_vertex);
            std::reverse(s.begin(), s.end());
            out = out + s + ">-";
            current_vertex = pred[current_vertex];
        }
        out = out + std::to_string(0);
        std::reverse(out.begin(), out.end());
        return out;
    }

    int print_result(int *dist, int *pred) {
        std::ofstream outputf(outputfile, std::ofstream::out);
        outputf << dist[0];
        for (int i = 1; i < N; i++) {
            outputf << " " << dist[i];
        }
        for (int i = 0; i < N; i++) {
            outputf << "\n";
            if (dist[i] >= 1000000) {
                outputf << "NO PATH";
            } else {
                outputf << format_path(i, pred);
            }
        }
        outputf << endl;
        return 0;
    }

    void find_local_minimal (int vertex_start, int vertex_end, int *dist, bool *visit, int &minValue, int &minIndex) {
        int min = INT_MAX;
        int u = -1;
        for (int i = vertex_start; i < vertex_end; ++i) {
            if (!visit[i]) {
                if (dist[i] < min) {
                    min = dist[i];
                    u = i;
                }
            }
        }
        minValue = min;
        minIndex = u;
    }


}//namespace utils

//Hint: use pthread condition variable or pthread barrier to do the synchronization
//------You may add helper functions and global variables here------

struct pthread_param {
    long thread_id;
    int num_threads;
    int *mat;
    int *global_minimal;
    int *global_minimal_index;
    long *global_minimal_thread;
    int *all_dist;
    bool *all_visit;
    int *all_pred;
    pthread_mutex_t *mutex;
    pthread_barrier_t *barrier;
};

void *thread_worker (void *worker_param);
void dijkstra(int N, int p, int *mat, int *all_dist, int *all_pred) {
    //------your code starts from here------
    long thread;
    pthread_t *thread_handles = (pthread_t*) calloc (p, sizeof (pthread_t));
    struct pthread_param **params = (struct pthread_param**) calloc (p, sizeof (struct pthread_param*));

    int *global_minimal = (int*) calloc (utils::N, sizeof (int));
    int *global_minimal_index = (int*) calloc (utils::N, sizeof (int));
    long *global_minimal_thread = (long *) calloc (utils::N, sizeof (long));

    bool *all_visit = (bool*) calloc (utils::N, sizeof (bool));

    pthread_mutex_t mutex;
    pthread_mutex_init (&mutex, nullptr);

    pthread_barrier_t barrier;
    pthread_barrier_init (&barrier, NULL, p);

    for (int i = 0; i < utils::N; ++i) {
        all_dist[i] = mat[utils::convert_dimension_2D_1D (0, i)];
        all_pred[i] = 0;
        all_visit[i] = false;
        global_minimal[i] = INT_MAX;
        global_minimal_index[i] = -1;
        global_minimal_thread[i] = 0l;
    }
    all_visit[0] = true;

    for (thread = 0; thread < p; ++ thread) {
        params[thread] = (struct pthread_param*) malloc (sizeof (struct pthread_param));
        params[thread]->thread_id = thread;
        params[thread]->global_minimal = global_minimal;
        params[thread]->global_minimal_index = global_minimal_index;
        params[thread]->global_minimal_thread = global_minimal_thread;
        params[thread]->num_threads = p;
        params[thread]->mat = mat;
        params[thread]->all_dist = all_dist;
        params[thread]->all_pred = all_pred;
        params[thread]->all_visit = all_visit;
        params[thread]->mutex = &mutex;
        params[thread]->barrier = &barrier;
        pthread_create (&thread_handles[thread], nullptr, thread_worker, (void*) params[thread]);
    }

    for (thread = 0; thread < p; ++ thread) {
        pthread_join (thread_handles[thread], nullptr);
    }

    pthread_barrier_destroy (&barrier);
    pthread_mutex_destroy (&mutex);

    for (thread = 0; thread < p; ++ thread) {
        free (params[thread]);
    }
    free (global_minimal);
    free (global_minimal_index);
    free (global_minimal_thread);
    free (all_visit);
    free (thread_handles);
    free (params);
    //------end of your code------
}

void *thread_worker (void *worker_param) {
    struct pthread_param *param = (struct pthread_param *) worker_param;
    long thread_id = param->thread_id;
    int num_threads = param->num_threads;

    int *mat = param->mat;
    int *all_dist = param->all_dist;
    bool *all_visit = param->all_visit;
    int *all_pred = param->all_pred;

    int *global_min = param->global_minimal;
    int *global_min_index = param->global_minimal_index;
    long *global_min_thread = param->global_minimal_thread;

    pthread_mutex_t *mutex = param->mutex;
    pthread_barrier_t *barrier = param->barrier;

    int local_n = utils::N / num_threads;
    if (thread_id == num_threads - 1) {
        local_n = utils::N - (num_threads - 1) * local_n;
    }
    int vertex_start = thread_id * (utils::N / num_threads);
    int vertex_end = vertex_start + local_n;

#ifdef MY_DEBUG
    std::cout << "Thread ID: " << thread_id << " handles vertex from: " << vertex_start
        << " to vertex: " << vertex_end << std::endl;
#endif

    for (int i = 1; i < utils::N; ++i) {

        int local_min = INT_MAX;
        int local_min_index = -1;
        utils::find_local_minimal (vertex_start, vertex_end, all_dist, all_visit, local_min, local_min_index);

        // Calculate the global minimal
        pthread_mutex_lock (mutex);
        if (local_min < global_min[i]) {
            global_min[i] = local_min;
            global_min_index[i] = local_min_index;
            global_min_thread[i] = thread_id;
        } else if (local_min == global_min[i] && thread_id < global_min_thread[i]) {
            global_min_index[i] = local_min_index;
            global_min_thread[i] = thread_id;
        }
        pthread_mutex_unlock (mutex);

        // Sync among threads to agree to the global minimal
        pthread_barrier_wait (barrier);

        // Copy the global variable to local
        int local_global_min = global_min[i];
        int u = global_min_index[i];

#ifdef MY_DEBUG
        std::cout << "Thread ID: " << thread_id << " copies global minimal: " << local_global_min << " at: " << u << std::endl;
#endif
        // Sync among threads after copying global variable
        pthread_barrier_wait (barrier);

        // Update the visited array
        if (u >= vertex_start && u < vertex_end) {
            all_visit[u] = true;
        }

        for (int v = vertex_start; v < vertex_end; ++v) {
            if (!all_visit[v]) {
                int new_dist = local_global_min + mat[utils::convert_dimension_2D_1D (u, v)];
                if (new_dist < all_dist[v]) {
                    all_dist[v] = new_dist;
                    all_pred[v] = u;
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    assert(utils::parse_args(argc, argv) == 0);
    assert(utils::read_file(utils::filename) == 0);

    assert(utils::num_threads <= utils::N);
    //`all_dist` stores the distances and `all_pred` stores the predecessors
    int *all_dist;
    int *all_pred;
    all_dist = (int *) calloc(utils::N, sizeof(int));
    all_pred = (int *) calloc(utils::N, sizeof(int));

    //time counter
    timeval start_wall_time_t, end_wall_time_t;
    float ms_wall;

    //start timer
    gettimeofday(&start_wall_time_t, nullptr);

    dijkstra(utils::N, utils::num_threads, utils::mat, all_dist, all_pred);

    //end timer
    gettimeofday(&end_wall_time_t, nullptr);
    ms_wall = ((end_wall_time_t.tv_sec - start_wall_time_t.tv_sec) * 1000 * 1000
               + end_wall_time_t.tv_usec - start_wall_time_t.tv_usec) / 1000.0;

    std::cerr << "Time(ms): " << ms_wall << endl;

    utils::print_result(all_dist, all_pred);

    free(utils::mat);
    free(all_dist);
    free(all_pred);

    return 0;
}
