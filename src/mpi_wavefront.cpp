#include <filesystem>
#include <cmath>
#include <iostream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "mpi.h"
#include <omp.h>

// MACROS
#define DEFAULT_INITIALIZATION 0.0 // default matrix value inzialization (double)s
#define DEFAULT_M_SIZE 10 // default size of the square matrix (NxN)
#define DEFAULT_GRANULARITY 3 // default size of the task (NxN)
#define DEFAULT_W 1 // default number of workers
#define DEFAULT_LOG_NAME "./log/mpi_wavefront.csv" // default log file name
#define DEFAULT_MACHINE 0 // default machine identifier

/*
Macros used to implement checks or debug printouts. 
Each debug value implies an increasing level of control, inspection and inefficiency of the code:
    1 - no additional printing, only consistency checks
    2 - very useful prints
    3 - Minor stdout printouts
    5 - Major stdout printouts
*/
#ifndef DEBUG
    #define DEBUG 0
#endif
// Execute code
#define DEBUG_CHECK(d,x) if (DEBUG >= d) {x}
// Print information
#define DEBUG_STDOUT(d,x) if (DEBUG >= d) { \
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()); \
    std::cout << "[" << std::put_time(std::localtime(&now), "%H:%M:%S") << "] " << x << std::endl; \
}

// To compute elapsed time
#define TIMERSTART(label)\
	std::chrono::time_point<std::chrono::system_clock> a##label, b##label;\
	a##label = std::chrono::system_clock::now();
#define TIMERSTOP(label, time_elapsed)\
	b##label = std::chrono::system_clock::now();\
	std::chrono::duration<double> delta##label = b##label-a##label;\
	time_elapsed = delta##label.count();



// CLASS AND STRUCT

// Define a template struct for efficiently operating on a triangular matrix
template <typename matrix_type>
struct TriangularMatrix {
    // Macro to compute the index in the internal array for the element (x, y)
    #define UM_INDEXING(x,y) y + x * matrix_dimension - (x * x + x)/2

    // Pointer to the internal array storing matrix elements
    matrix_type* matrix_array;

    // Boolean flag to determine if the matrix is column-major (true) or row-major (false)
    const bool column_matrix;

    // Size of the internal array, which is the number of elements in the upper triangular matrix
    const int array_size;

    // Dimension of the matrix (number of rows or columns)
    const int matrix_dimension;

    // Constructor for the TriangularMatrix class
    // Initializes the matrix with a given dimension, default value, and optional column-major flag
    TriangularMatrix(const int matrix_dimension, matrix_type default_value, const int column_matrix = 0):
        column_matrix(column_matrix),
        // Calculate the number of elements in the upper triangular matrix
        array_size((matrix_dimension*matrix_dimension + matrix_dimension )/2),
        // Initialize matrix dimension
        matrix_dimension(matrix_dimension) {
        
        // Allocate memory for the array of matrix elements
        matrix_array = new matrix_type[array_size];

        // Initialize all elements in the matrix to the default value
        for (int i = 0; i < array_size; i++) {
            matrix_array[i] = default_value;
        }
    }

    // Destructor to release the allocated memory
    ~TriangularMatrix() {
        delete[] matrix_array;
    }

    // Method to access an element at position (i, j)
    // Handles column-major order if specified
    matrix_type& element(int i, int j) {
        // Adjust indices if the matrix is in column-major order
        if (column_matrix) {
            int dummy_val = i;
            i = matrix_dimension - 1 - j;
            j = matrix_dimension - 1 - dummy_val;
        }

        // Debug checks to ensure indices are within bounds and element is in the upper triangular part
        DEBUG_CHECK(1,
            if (i < 0 || i >= matrix_dimension || j < 0 || j >= matrix_dimension ) {
                throw std::out_of_range("Index out of bounds. i = " + std::to_string(i) + ", j = " + std::to_string(j));
            }
            else if (j < i) {
                throw std::out_of_range("Element not in the triangular matrix. i = " + std::to_string(i) + ", j = " + std::to_string(j));
            }
        )
        // Return the matrix element at the computed index
        return matrix_array[UM_INDEXING(i, j)];
    }

    // Overload the stream insertion operator to print the matrix
    // Prints the matrix in a human-readable format
    friend std::ostream & operator<<(std::ostream& os, TriangularMatrix<matrix_type> & tc) {

        std::stringstream ret;
        std::string token(8, 'X'); // A placeholder string for elements below the diagonal

        // Iterate over matrix dimensions to print the matrix
        for (int i = 0; i < tc.matrix_dimension; i++) {
            for (int j = 0; j < tc.matrix_dimension; j++) {
                // Print 'X' for elements below the diagonal, otherwise print the matrix element
                if (i > j) 
                    ret << token << " ";
                else
                    ret << tc.element(i, j) << " ";
            }
            ret << std::endl; // New line for each row
        }
        return os << ret.str() << std::endl; // Output the matrix as a string
    }
};



// FUNCTIONS

// Function to print the usage instructions for the program
void print_usage() {
    std::cout << "usage: fastflow_wf [options]" << std::endl
        << "\t-h\t\tprints this message" << std::endl
        << "\t-N matrix_size\t\tsize of the square matrix [default=" << DEFAULT_M_SIZE << "]" << std::endl
        << "\t-G granularity\t\tsize of the subtask [default=" << DEFAULT_GRANULARITY << "]" << std::endl
        << "\t-W num_workers\t\tnumber of workers [default=" << DEFAULT_W << "]" << std::endl
        << "\t-l log_name\t\tlog file name [default=" << DEFAULT_LOG_NAME << "]" << std::endl
        << "\t-m machine_id\t\tinteger to diversify the machine used for testing [default=" << DEFAULT_MACHINE << "]" << std::endl;
}

// Function to initialize a file
// Creates necessary directories and files if they do not exist
bool init_file(const std::string& file_name) {

    // Get the parent directory of the file
    std::filesystem::path directory = std::filesystem::path(file_name).parent_path();
    
    // Check if the directory exists
    if (!std::filesystem::exists(directory)) {
        // Create the directory
        std::filesystem::create_directories(directory);
    }

    // Open the file in read mode
    std::ifstream file(file_name);
    
    // If the file exists and can be opened successfully, return true
    if (file.good()) {
        file.close(); // Close the file after checking
        return true;
    }
    file.close(); // Close the file if it could not be opened

    // Open the file in write mode to create a new file
    std::ofstream new_file(file_name);
    
    // If the new file could not be opened, return false
    if (!new_file.is_open()) {
        return false;
    }
    
    // Write a header line to the new file
    new_file << "num_process,matrix_size,num_workers,machine_identifier,granularity,debug_value,init_time,running_time" <<std::endl;
    new_file.close(); // Close the new file after writing the header
    return true;

}

// Compute one element of the wavefront computation that comprise the vectorial sum 
// and the cubic root
double compute(TriangularMatrix<double>* row_tm, TriangularMatrix<double>* col_tm, int i, int j){
    double partial_sum = 0;
    for (int k=0; k < (j-i); k++){
        // Compute partial sum for the row and column adjacent vector
        // Here the memory ordering significantly influence computation time 
        partial_sum += row_tm->element(i, i + k) * col_tm->element(j - k, j);
    }
    // Return the cubic root of the partial sum
    return std::cbrt(partial_sum);
}

// MAIN
int main (int argc, char *argv[]){
	// Initialize MPI
    // int prov;
	// MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &prov);
    MPI_Init(&argc, &argv);

    // store number of process and id of "this" process
	int numP;
	int myId;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);
    MPI_Comm_size(MPI_COMM_WORLD, &numP);

    // Set default values for matrix size, granularity, number of workers, machine identifier, and log file name
	int matrix_size = DEFAULT_M_SIZE;
    int granularity = DEFAULT_GRANULARITY;
	int num_workers = DEFAULT_W;
    int machine_identifier = DEFAULT_MACHINE;
	std::string log_file_name = DEFAULT_LOG_NAME;

    // Parse command-line arguments for configuration
	int opt;
	while ((opt = getopt(argc, argv, "hN:G:W:m:l:")) != -1) {
		switch (opt) {
            case 'm':
                // Set machine identifier
				machine_identifier = std::stoi(optarg);
				break;
            case 'G':
                // Set granularity of tasks
				granularity = std::stoi(optarg);
				break;
			case 'N':
                // Set matrix size
				matrix_size = std::stoi(optarg);
				break;
			case 'W':
                // Set number of workers
				num_workers = std::stoi(optarg);
				break;
			case 'l':
                // Set log file name
				log_file_name = optarg;
				break;
            case 'h':
                // Print usage instructions
                print_usage();
				MPI_Abort(MPI_COMM_WORLD, 0);
            case '?':       
            default:
                // Handle invalid options
				print_usage();
				MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

    // call to OpenMP to set a fixed number of threads per process
    omp_set_num_threads(num_workers);
    // Only if "this" is the main process we perform these checks and debug prints
    if(!myId){
        // Print initial debug message indicating the debug mode is active with the current DEBUG level
        DEBUG_STDOUT(2, "Debug mode is active - DEBUG=" << DEBUG);

        // Validate matrix size and granularity inputs
        if(matrix_size <= 1){
            std::cout << "matrix_size needs to be > 1" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        // Validate granularity
        if(granularity < 1 || granularity > matrix_size){
            std::cout << "granularity needs to be > 0 and <= matrix_size" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        // Output debug information with the final configuration
        DEBUG_STDOUT(3,
            "Starting with parameter: N=" << matrix_size << " G=" << granularity << " W=" << num_workers <<
            " m=" << machine_identifier << " l=" <<log_file_name
        )

        // Initialize the log file
        if(!init_file(log_file_name)){
            std::cout << "log file creation error" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }  
        // Print debug information about log file initialization
        DEBUG_STDOUT(3,
                    "Check or creation of the log file complete"
                )
    }

    // Record initialization time
    double init_time;
    /*we calculate the size of the task matrix, 
    in order to handle the problem at task level and not 
    at the level of individual elements of the original matrix, 
    but there is no need to create the actual matrix as the code 
    is  inherently thread and process safe*/
    int task_matrix_size = std::ceil(static_cast<double>(matrix_size - 1)/granularity);
    int task_work_load;
    TIMERSTART(mat_init)
    // Initialization of the row and column element matrices, 
    // one for each process because there is no shared memory
    TriangularMatrix<double> row_tm(matrix_size, DEFAULT_INITIALIZATION);
    TriangularMatrix<double> col_tm(matrix_size, DEFAULT_INITIALIZATION, 1);
    for (int k = 0; k < matrix_size; k++){
        row_tm.element(k,k) = (k+1.0)/matrix_size;
        col_tm.element(k,k) = (k+1.0)/matrix_size;
    }
    // Stop the initialization timer and print elapsed time
    TIMERSTOP(mat_init, init_time)
    if(!myId){
        DEBUG_STDOUT(3,
            "Matrix initialization complete - elapsed time: " << init_time
        )
    }
    
    // Record running time
    double running_time;

    TIMERSTART(run)
    // Core section of the code implementing the multi-process wavefront calculation.
    // Each process handles a portion of the diagonal of the task matrix, which represents 
    // a high-level subdivision of the original problem into smaller subtasks.
    // Within each process, threads are utilized exploiting data-driven parallelism to
    // efficiently compute these subtasks, ensuring that each diagonal segment of the matrix
    // is processed concurrently by multiple threads.

    // Iterate over the diagonals of the task matrix (virtual matrix)
    for (int diag = 0; diag < task_matrix_size; diag++){
        
        // the size of the subtask is computed
        if (diag == 0) task_work_load = (static_cast<double>(granularity) * (granularity + 1))/2;
        else  task_work_load = granularity * granularity;

        // the lenght of the current task matrix diagonal 
        int total_work_load = task_matrix_size - diag;
        // the number of task per process
        int proc_work_load = std::ceil(static_cast<double>(total_work_load) / numP);
        // the starting point of computation in the task matrix diagonal of each process
        int proc_start = proc_work_load*myId;
        // the ending point of computation in the task matrix diagonal of each process
        int proc_end = std::min(proc_work_load*(myId+1), total_work_load);
        // the real process work load, eliminating any overflow
        int real_proc_work_load = std::max(0, proc_end - proc_start);

        
        

        // All the values and arrays required for exchanging data between processes, some comment are from the official docmentation
        // iterator and variable for the computation
        int next_index = 0, size;
        // integer array (of length group size). 
        // Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from process i (integer)*/
        int* displs = new int[numP];
        // non-negative integer array (of length group size) 
        // containing the number of elements that are received from each process (non-negative integer)
        int* recvcounts = new int[numP];
        for (int k = 0;  k < numP; k++){
            size = std::max(0, (std::min(proc_work_load*(k+1), total_work_load) - proc_work_load*k)) * task_work_load;
            recvcounts[k] = size;
            displs[k] = next_index;
            next_index += size;
        }
        
        // print debug information
        DEBUG_STDOUT(5, "[" << myId << "] " << total_work_load << " - "  << proc_work_load << " - " << proc_start << " - " << proc_end << " (" << task_work_load << ")")
        // starting address of send buffer (choice)
        double* sendbuf = new double[real_proc_work_load * task_work_load];
        // Iterator used to correctly populate the buffer with computed values.
        int buff_it = 0;
        // Each process compute the respective section of the task matrix diagonal
        for (int k = proc_start; k < proc_end; k++){
            
            /* Each task is a similar subproblem to the original one and is solved in a similar way, 
               in which case the lower triangular matrix must also be taken into account 
               for all tasks that are not in the first daigonal*/
            int i = k, j = diag+k;
            int stop;

            // skipped in the first iteration, lower triangular part of the subtask
            if (diag != 0) {
                // iterate over the element of the diagonal
                for (int internal_diag = 0; internal_diag < granularity - 1; internal_diag++){      
                    
                    stop = std::max(-1,
                                std::max(i*granularity + granularity - 1 - matrix_size,
                                j*granularity + 1 + internal_diag - matrix_size));

                    // OpenMP is used to easily parallelise this cycle, 
                    // in which each iteration is independent of the others.
                    #pragma omp parallel for
                    for (int internal_k=internal_diag; internal_k>stop; internal_k--){

                        // Here, the coordinates relative to the element matrix are calculated and used
                        int abs_i = i*granularity + (granularity - internal_k - 1);
                        int abs_j = j*granularity + 1 + (internal_diag - internal_k);

                        // The send buffer is popolated
                        int atual_buff_it = buff_it + internal_diag - internal_k;
                        sendbuf[atual_buff_it] = compute(&row_tm, &col_tm, abs_i, abs_j);
                        row_tm.element(abs_i, abs_j) = sendbuf[atual_buff_it];
                        col_tm.element(abs_i, abs_j) = sendbuf[atual_buff_it];
                    }
                    // The send buffer iterator is increased
                    buff_it += internal_diag - stop;
                }
            }
            // upper triangular part of the subtask
            for (int internal_diag = 0; internal_diag < granularity; internal_diag++){
                
                stop = std::min(granularity - internal_diag,
                                std::min(matrix_size - 1 - internal_diag - j*granularity,
                                matrix_size - i*granularity));

                // OpenMP is used to easily parallelise this cycle, 
                // in which each iteration is independent of the others.
                #pragma omp parallel for
                for (int internal_k=0; internal_k < stop; internal_k++){

                    // Here, the coordinates relative to the element matrix are calculated and used
                    int abs_i = i*granularity + (internal_k);
                    int abs_j = j*granularity + 1 + (internal_diag + internal_k);

                    int atual_buff_it = buff_it + internal_k;

                    // The send buffer is popolated
                    sendbuf[atual_buff_it] = compute(&row_tm, &col_tm, abs_i, abs_j);
                    row_tm.element(abs_i, abs_j) = sendbuf[atual_buff_it];
                    col_tm.element(abs_i, abs_j) = sendbuf[atual_buff_it];
                }
                // The send buffer iterator is increased
                buff_it += stop;
            }
        }

        // print debug information if DEBUG > 5
        DEBUG_CHECK(5,
            std::cout << "SEND:";
            for (int p = 0; p < real_proc_work_load * task_work_load; p++){
                std::cout << " " << sendbuf[p];
                if (p%task_work_load == task_work_load - 1) std::cout << " |";
            }
            std::cout << std::endl;
        )

        /*Here, the function MPI_Allgatherv is used to enable the exchange of computed elements between processes 
        and ensure that each process has all the elements calculated up to this point. 
        Each subsequent iteration may require the presence of these elements.*/
        // Address of receive buffer (choice)
        double* recvbuf = new double[total_work_load*task_work_load];
        MPI_Allgatherv(sendbuf, real_proc_work_load * task_work_load , MPI_DOUBLE, 
                        recvbuf, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

        // print debug information if DEBUG > 5
        DEBUG_CHECK(5,
            std::cout << "REC:";
            for (int p = 0; p < total_work_load*task_work_load; p++){
                std::cout << " " << recvbuf[p];
                if (p%task_work_load == task_work_load - 1) std::cout << " |";
            }
            std::cout << std::endl;
        )

        // Now each process must popolate in its own element matrices in order to be able to 
        // use all elements received and calculated by the other processes
        buff_it = 0;
        for (int k = 0; k < total_work_load; k++){
            // The elements calculated by "this" process are already present
            if (k >= proc_start && k < proc_end){
                buff_it += task_work_load;
                continue;
            }
            int i = k, j = diag+k;
            int stop;

            /* In this section, the iterations follow the original order. 
               Note that parallelization is not used here because, 
               after testing, it was found that parallelizing this memory-intensive section actually 
               decreased performance. Thus, for efficiency, this section is executed in a single-threaded manner. */
            // for comments please refer to the previous similar section.
            if (diag != 0) {
                for (int internal_diag = 0; internal_diag < granularity - 1; internal_diag++){      
                    
                    stop = std::max(-1,
                                std::max(i*granularity + granularity - 1 - matrix_size,
                                j*granularity + 1 + internal_diag - matrix_size));

                    for (int internal_k=internal_diag; internal_k > stop; internal_k--){
                        int abs_i = i*granularity + (granularity - internal_k - 1);
                        int abs_j = j*granularity + 1 + (internal_diag - internal_k);
                        
                        int atual_buff_it = buff_it + internal_k - internal_diag;
                        row_tm.element(abs_i, abs_j) = recvbuf[atual_buff_it];
                        col_tm.element(abs_i, abs_j) = recvbuf[atual_buff_it];
                    }
                    buff_it += internal_diag - stop;
                }
            }
            for (int internal_diag = 0; internal_diag < granularity; internal_diag++){

                stop = std::min(granularity - internal_diag,
                                std::min(matrix_size - 1 - internal_diag - j*granularity,
                                matrix_size - i*granularity));

                for (int internal_k=0; internal_k < stop; internal_k++){
                    int abs_i = i*granularity + (internal_k);
                    int abs_j = j*granularity + 1 + (internal_diag + internal_k);
                    
                    int atual_buff_it = buff_it + internal_k;
                    row_tm.element(abs_i, abs_j) = recvbuf[atual_buff_it];
                    col_tm.element(abs_i, abs_j) = recvbuf[atual_buff_it];
                }
                buff_it += stop;
            }
        }

        // freeing allocated memory
        delete[] sendbuf;
        delete[] recvbuf;
        delete[] displs;
        delete[] recvcounts;
    }
    TIMERSTOP(run, running_time)

    // Only if "this" is the main process we perform these checks and debug prints
    if(!myId){
        // debug printing
        DEBUG_STDOUT(5,
            "Final Matrix row:" << std::endl << row_tm <<  "Final Matrix col:" << std::endl << col_tm
        )
        
        DEBUG_STDOUT(2,"Final value (matrix dim " << matrix_size << "): " << row_tm.element(0, row_tm.matrix_dimension - 1))
        DEBUG_STDOUT(2, "Computation complete - elapsed time: " << running_time*1000)
        

        // write the execution times to a file
        std::ofstream file;
        file.open(log_file_name, std::ios_base::app);
        file << numP << "," << matrix_size << "," << num_workers << "," << machine_identifier << "," << granularity << "," << DEBUG << ","
            << init_time*1000  << "," << running_time*1000 << std::endl;
        file.close();
    
    }
    MPI_Finalize(); 
}
