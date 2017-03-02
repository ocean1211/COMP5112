/*
 * This is code skeleton for COMP5112 assignment1
 * Compile: mpic++ -o mpi_dijkstra mpi_dijkstra.cpp
 * Run: mpiexec -n <number of processes> mpi_dijkstra <input file>, you will find the output in 'output.txt' file
 */


#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <climits>
#include <cstring>
#include <algorithm>
#include "mpi.h"

using std::string;
using std::cout;
using std::endl;
using std::vector;

//#define MY_DEBUG

/*
 * utils is a namespace for utility functions
 * including I/O (read input file and print results) and one matrix dimension convert(2D->1D) function
 */
namespace utils {
    int N; //number of vertices
    int *mat; // the adjacency matrix

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
        std::ofstream outputf("output.txt", std::ofstream::out);
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
}//namespace utils
// you may add some helper functions here.

void find_local_minimum(int vertex_start, int vertex_end, int *dist, bool *visit, int *out) {
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
    out[0] = min;
    out[1] = u;
}

void dijkstra(int my_rank, int N, int p, MPI_Comm comm, int *mat, int *all_dist, int *all_pred) {

    //------your code starts from here------
    int loc_N; // I need a local copy for N
    int loc_n; //how many vertices I need to process.
    int *loc_mat; //local matrix
    int *loc_dist; //local distance
    int *loc_pred; //local predecessor
    bool *loc_visit; //local visit record array

    int local_min[2];
    int global_min[2];

    int *displs;
    int *rcounts;

    //step 1: broadcast N
    if (my_rank == 0) {
        loc_N = N;
    }
    MPI_Bcast (&loc_N, 1, MPI_INT, 0, comm);
#ifdef MY_DEBUG
    std::cout << "MPI process: " << my_rank << " received N: " << loc_N << std::endl;
#endif

    //step 2: find loc_n
    loc_n = loc_N / p;
    if (my_rank == p - 1) {
        loc_n = loc_N - (p - 1) * loc_n;
    }

#ifdef MY_DEBUG
    std::cout << "MPI process: " << my_rank << " find local n: " << loc_n << std::endl;
#endif

    //step 3: allocate local memory
    loc_mat = (int *) calloc(loc_N * loc_N, sizeof(int));
    loc_dist = (int *) calloc(loc_N, sizeof(int));
    loc_pred = (int *) calloc(loc_N, sizeof(int));
    loc_visit = (bool *) calloc(loc_N, sizeof(int));

    displs = (int *) calloc (p, sizeof(int));
    rcounts = (int *) calloc (p, sizeof (int));

    int displ = 0;
    for (int i = 0; i < p; ++i) {
        if (i != p -1) {
            rcounts[i] = loc_N / p;
            displs[i] = displ;
        } else {
            rcounts[i] = loc_N - (p - 1) * loc_n;
            displs[i] = displ;
        }
        displ += rcounts[i];
    }

    //step 4: broadcast matrix mat
    if (my_rank == 0) {
        memcpy (loc_mat, mat, loc_N * loc_N * sizeof (int));
    }
    MPI_Bcast (loc_mat, loc_N * loc_N, MPI_INT, 0, comm);

#ifdef MY_DEBUG
    std::cout << "MPI process: " << my_rank << " receive vertex matrix: [";
    for (int i = 0; i < loc_N * loc_N; ++i) {
        std::cout << loc_mat[i] << ",";
    }
    std::cout << "]" << std::endl;
#endif

    //step 4: dijkstra algorithm
    int vertex_start = my_rank * (loc_N / p);
    int vertex_end = vertex_start + loc_n;

#ifdef MY_DEBUG
    std::cout << "MPI process: " << my_rank << " will process vertex: [" << vertex_start
        << "," << vertex_end << ")" << std::endl;
#endif

    for (int i = 0; i < loc_N; ++i) {
        loc_dist[i] = loc_mat[utils::convert_dimension_2D_1D (0, i)];
        loc_pred[i] = 0;
        loc_visit[i] = false;
    }

    loc_visit[0] = true;

    for (int i = 1; i < loc_N; ++i) {
        find_local_minimum (vertex_start, vertex_end, loc_dist, loc_visit, local_min);
#ifdef MY_DEBUG
        std::cout << "MPI process: " << my_rank << " find local min: " << local_min[0] << " with vertex: "
            << local_min[1] << std::endl;
#endif

        MPI_Allreduce (local_min, global_min, 1, MPI_2INT, MPI_MINLOC, comm);
        int u = global_min[1];
        int globalMin = global_min[0];

#ifdef MY_DEBUG
        std::cout << "MPI process: " << my_rank << " receive global min: " << globalMin <<
            " with vertex: " << u << std::endl;
#endif
        loc_visit[u] = true;
        for (int v = vertex_start; v < vertex_end; ++v) {
            if (!loc_visit[v]) {
                int new_dist = globalMin + loc_mat[utils::convert_dimension_2D_1D(u, v)];
                if (new_dist < loc_dist[v]) {
                    loc_dist[v] = new_dist;
                    loc_pred[v] = u;
                }
            }
        }
    }

    //step 5: retrieve results back
    //Hint: use MPI_Gather(or MPI_Gatherv) function

    MPI_Gatherv (loc_dist + vertex_start, loc_n, MPI_INT, all_dist, rcounts, displs , MPI_INT, 0, comm);
    MPI_Gatherv (loc_pred + vertex_start, loc_n, MPI_INT, all_pred, rcounts, displs , MPI_INT, 0, comm);

    //step 6: remember to free memory
    free(loc_mat);
    free(loc_dist);
    free(loc_pred);
    free(loc_visit);

    //------end of your code------
}

int main(int argc, char **argv) {
    assert(argc > 1 && "Input file was not found!");
    string filename = argv[1];
    assert(utils::read_file(filename) == 0);

    //`all_dist` stores the distances and `all_pred` stores the predecessors
    int *all_dist;
    int *all_pred;
    all_dist = (int *) calloc(utils::N, sizeof(int));
    all_pred = (int *) calloc(utils::N, sizeof(int));

    //MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm comm;

    int p;//number of processors
    int my_rank;//my global rank

    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &my_rank);

    dijkstra(my_rank, utils::N, p, comm, utils::mat, all_dist, all_pred);

    if (my_rank == 0)
        utils::print_result(all_dist, all_pred);
    MPI_Finalize();

    free(utils::mat);
    free(all_dist);
    free(all_pred);

    return 0;
}
