#include <vector>
#include <random>
#include <ff/ff.hpp>

#include <cmath>
#include <iostream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <filesystem>

using namespace ff;

// mutex only used if debug mode is used
#include <mutex>
std::mutex mtx;



// MACROS
#define DEFAULT_INITIALIZATION 0.0 // default matrix value inzialization (double)s
#define DEFAULT_M_SIZE 10 // default size of the square matrix (NxN)
#define DEFAULT_GRANULARITY 3 // default size of the task (NxN)
#define DEFAULT_W 2 // default number of workers
#define DEFAULT_LOG_NAME "./log/fastflow_stream_dm.csv" // default log file name
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
// execute code
#define DEBUG_CHECK(d,x) if (DEBUG >= d) {if (DEBUG >= 5) std::lock_guard<std::mutex> lock(mtx);x}
// print information
#define DEBUG_STDOUT(d,x) if (DEBUG >= d) { \
    if (DEBUG >= 5) std::lock_guard<std::mutex> lock(mtx); \
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()); \
    std::cout << "[" << std::put_time(std::localtime(&now), "%H:%M:%S") << "] " << x << std::endl; \
}

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



// Define types for task source and task sink
// task_source contains information to process a part of the matrix
using task_source = std::tuple<TriangularMatrix<double>*, TriangularMatrix<double>*, int, int, int, int>;
// task_sink is used to indicate completion of a processing task
using task_sink = std::tuple<int, int>;

struct SourceSinkDoubleMatrix: ff_monode_t<task_sink, task_source> {

    const double init_value = DEFAULT_INITIALIZATION;

    // 2 matrices where to store the elements row-major and column-major order
    TriangularMatrix<double> row_tm;
    TriangularMatrix<double> col_tm;
    // Matrix where to keep track of processed elements
    TriangularMatrix<bool> presence_um;

    // Dimensions and granularity of the matrix
    const int true_matrix_dimension;
	const int pres_matrix_dimension;
    // Used to determine the task dimension and to divide the global task
    const int granularity;

    // Constructor to initialize the SourceSinkDoubleMatrix
    SourceSinkDoubleMatrix(const int matrix_dimension, const int granularity, double* diagonal):
        // Initialisation real matrices, element matrices
        row_tm(matrix_dimension, init_value),
        col_tm(matrix_dimension, init_value, 1),
        // Initialisation virtual matrix, task matrix
        presence_um(std::ceil(static_cast<double>(matrix_dimension - 1)/granularity), 0),
        // Dimension (number of row or column) of the real element matrix
        true_matrix_dimension(matrix_dimension),
        // Dimension (number of row or column) of the presence matrix, virtual matrix, task matrix
        pres_matrix_dimension(std::ceil(static_cast<double>(matrix_dimension - 1)/granularity)),
        granularity(granularity)
        {
            // Output debug information abount sink initialization, level 3
            DEBUG_STDOUT( 3, "Sink node initialized with dim: " << matrix_dimension )
            
            // Initialize the diagonal elements of the matrices
            for (int k = 0; k < true_matrix_dimension; k++){
                row_tm.element(k,k) = diagonal[k];
                col_tm.element(k,k) = diagonal[k];
            }
            // Free the memory allocated for the diagonal array
            delete[] diagonal;
        }

    // Called when the node is initializated, send initial tasks to generate the "sub-task" diagonal
    int svc_init(){
        DEBUG_STDOUT( 3,
            "sink svc_init()"
        )

        // Send tasks for each row of the matrix
        for (int k = 0; k < pres_matrix_dimension; k++){
            ff_send_out(new task_source(&row_tm, &col_tm, k, k, granularity, 1));
        }
        return 0; 
    }

    // Process received (completed) tasks
    task_source* svc(task_sink* task) {

        // Handle null tasks sended by fastflow
        if (task == nullptr) {
            DEBUG_STDOUT(5, "SOURCE: null "<< task)

            return this->GO_ON;
        }

        // tuple unpacking
        int i = std::get<0>(*task);
        int j = std::get<1>(*task);
        DEBUG_STDOUT(5, "SOURCE: " << i << " - " << j)

        // Mark the received task as processed
        presence_um.element(i,j) = 1;

        // tasks are sent if they can be calculated, depending on the presence of the two adjacent tasks.
        if (i-1 >= 0  && j-1 >= 0 && presence_um.element(i-1, j-1))
            ff_send_out(new task_source(&row_tm, &col_tm, i-1 , j, granularity, 0));
        if (i+1 < pres_matrix_dimension  && j+1 < pres_matrix_dimension && presence_um.element(i+1, j+1))
            ff_send_out(new task_source(&row_tm, &col_tm, i , j+1, granularity, 0));

        // Broadcast end-of-stream signal if is reached the end of the matrix 
        if (i == 0 && j == pres_matrix_dimension - 1){
            broadcast_task(EOS);
            return EOS;
        }
        else return this->GO_ON;
    }

    // Printing debug information at the termination
    void svc_end(){
        DEBUG_STDOUT(3,
            "sink svc_end()"
        )
        DEBUG_STDOUT(5,
            "Final Matrix row:" << std::endl << row_tm <<  "Final Matrix col:" << std::endl << col_tm
        )
        DEBUG_STDOUT(2,
            "Final value (matrix dim " << true_matrix_dimension << "): " << row_tm.element(0, true_matrix_dimension-1)
        )
    } 
    
    // Printing debug information
    void eosnotify(ssize_t ch){
        DEBUG_STDOUT(3,
            "sink eosnotify()" 
        )
    } 
};

// Class that implements the worker nodes
struct Worker: ff_node_t<task_source, task_sink> {

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
    
    // Process incoming tasks sent by the source-sink
    task_sink* svc(task_source* task) {

        // tuple unpacking
        TriangularMatrix<double>* row_tm = std::get<0>(*task);
        TriangularMatrix<double>* col_tm = std::get<1>(*task);
        int i = std::get<2>(*task);
        int j = std::get<3>(*task);
        int granularity = std::get<4>(*task);
        int first_it = std::get<5>(*task);

        // Print debug information if DEBUG > 5
        DEBUG_STDOUT(5, "WORKER: " << i << " - " << j << " - " << first_it)

        // Variables for the coordinates relative to the element in the matrix, not the task matrix
        int abs_i;
        int abs_j;
        // Used to avoid repeating the computation
        double dummy_val;

        // If this is NOT the first iteration, calculate the lower triangular part of the task,
        // if this is the firts iteration, the lower triangular matrix is not needed and not computable
        if (!first_it) {
            // Iterate over diagonal
            for (int diag = 0; diag < granularity - 1; diag++){      
                // Iterate over diagonal element
                for (int k=diag; k>=0; k--){
                    // Compute the absolute (element matrix) coordinates
                    abs_i = i*granularity + (granularity - k - 1);
                    abs_j = j*granularity + 1 + (diag - k);
                    
                    // Break if out of bounds
                    if (abs_i >= row_tm->matrix_dimension || abs_j >= row_tm->matrix_dimension) break;

                    // Update element matrices with the computed value
                    dummy_val = compute(row_tm, col_tm, abs_i, abs_j);
                    row_tm->element(abs_i, abs_j) = dummy_val;
                    col_tm->element(abs_i, abs_j) = dummy_val;
                }
            }
        }

        // Compute upper triangular section of the task for all iteration 
        // (this use the lower part if is not the first iteration)
        // Iterate over diagonal
        for (int diag = 0; diag < granularity; diag++){
            // Iterate over diagonal element
            for (int k=0; k< (granularity - diag); k++){
                // Compute the absolute (element matrix) coordinates
                abs_i = i*granularity + (k);
                abs_j = j*granularity + 1 + (diag + k);
                
                // Break if out of bounds
                if (abs_i >= row_tm->matrix_dimension || abs_j >= row_tm->matrix_dimension) 
                    break;
                
                // Update element matrices with the computed value
                dummy_val = compute(row_tm, col_tm, abs_i, abs_j);
                row_tm->element(abs_i, abs_j) = dummy_val;
                col_tm->element(abs_i, abs_j) = dummy_val;
            }
        }

        // Send completed task to source-sink
		ff_send_out(new task_sink(i, j));
		delete task; // Clean up task memory
		return GO_ON; // Continue processing; do not terminate
    }

    // Printing debug information
    int svc_init(){
        DEBUG_STDOUT(3,
            "worker svc_init()"
        )
        return 0; 
    }

    // Printing debug information
    void svc_end(){
        DEBUG_STDOUT(3,
            "worker svc_end()"
        )
    } 

    // Printing debug information
    void eosnotify(ssize_t ch){
        DEBUG_STDOUT(3,
            "worker eosnotify()"
        )
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
    new_file << "matrix_size,num_workers,machine_identifier,granularity,debug_value,init_time,running_time" <<std::endl;
    new_file.close(); // Close the new file after writing the header
    return true;

}

// MAIN

int main(int argc, char *argv[]) {   
    // Print initial debug message indicating the debug mode is active with the current DEBUG level
    DEBUG_STDOUT(2, "Debug mode is active - DEBUG=" << DEBUG);

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
                return 0;
            case '?':       
            default:
                // Handle invalid options
                print_usage();
                return 1;
        }
    }

    // Validate matrix size and granularity inputs
    if(matrix_size <= 1){
        std::cout << "matrix_size needs to be > 1" << std::endl;
        return 1;
    }
    // Validate granularity
    if(granularity < 1 || granularity > matrix_size){
        std::cout << "granularity needs to be > 0 and <= matrix_size" << std::endl;
        return 1;
    }

    // Output debug information with the final configuration
    DEBUG_STDOUT(3,
        "Starting with parameter: N=" << matrix_size << " G=" << granularity << " W=" << num_workers <<
        " m=" << machine_identifier << " l=" <<log_file_name
    )

    // Initialize the log file
    if(!init_file(log_file_name)){
        std::cout << "log file creation error" << std::endl;
        return 1;
    }
    // Print debug information about log file initialization
    DEBUG_STDOUT(3,
                "Check or creation of the log file complete"
            )

    // Record initialization time
    double init_time;
    TIMERSTART(init)
    // Initialization of matrix values
    double* init = new double[matrix_size];
    for (int i = 0; i < matrix_size; i++){
        init[i] = (i + 1.0) / matrix_size;
    }

    // Create the Source-Sink object with initial values
    SourceSinkDoubleMatrix producer(matrix_size, granularity, init);

    // Configure the FastFlow farm with custon worker nodes
    ff_Farm<task_source, task_sink> farm(
                            [&]() {
                                std::vector<std::unique_ptr<ff_node> > W;
                                // Initialize worker nodes
                                for(auto i = 0; i < num_workers; ++i)
                                    W.push_back(make_unique<Worker>());
                                return W;
                            } (),
                            producer);
    // Remove the default collector from the farm as it's not needed
    farm.remove_collector();
    // Wrap the farm to enable feedback from Workers to the Emitter
    farm.wrap_around();

    // Stop the initialization timer and print elapsed time
    TIMERSTOP(init, init_time)
    DEBUG_STDOUT(3,
                "Initialization complete - elapsed time: " << init_time
            )


    // Record running time
    double running_time;
    TIMERSTART(run)
    // Run the FastFlow farm and wait for completion
    if (farm.run_and_wait_end() < 0) {
        error("running farm");
        return -1;
    }
    // Stop the computation timer and print elapsed time
    TIMERSTOP(run, running_time)
    DEBUG_STDOUT(2,
                "Computation complete - elapsed time: " << running_time * 1000
            )

    // Log the execution times and configuration parameters to the log file
    std::ofstream file;
    file.open(log_file_name, std::ios_base::app);
    file << matrix_size << "," << num_workers << "," << machine_identifier << "," << granularity << "," << DEBUG << ","
        << init_time * 1000  << "," << running_time * 1000 << std::endl;
    file.close();

    return 0;
}
