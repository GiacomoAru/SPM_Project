#include <vector>
#include <random>
#include <ff/ff.hpp>
#include <ff/parallel_for.hpp>

#include <iostream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <filesystem>

#include <thread>

using namespace ff;

// mutex only used if debug mode is used
#include <mutex>
std::mutex mtx;



// MACROS

#define DEFAULT_INITIALISATION 0.0 // default matrix value inzialization (double)s
#define DEFAULT_M_SIZE 10 // default size of the square matrix (NxN)
#define DEFAULT_W 2 // default number of workers
#define DEFAULT_LOG_NAME "./log/fastflow_data_sm.csv" // default log file name
#define DEFAULT_MACHINE 0 // default machine identifier


/*
Macros used to implement checks or debug printouts. 
Each debug value implies an increasing level of control, inspection and inefficiency of the code:
    1 - no additional printing, only consistency checks
    2 - very useful prints
    3 - Minor stdout printouts
    5 - Major stdout printouts + mutex
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

// to easily compute elapsed time
#define TIMERSTART(label)\
	std::chrono::time_point<std::chrono::system_clock> a##label, b##label;\
	a##label = std::chrono::system_clock::now();
#define TIMERSTOP(label, time_elapsed)\
	b##label = std::chrono::system_clock::now();\
	std::chrono::duration<double> delta##label = b##label-a##label;\
	time_elapsed = delta##label.count();



// CLASS AND STRUCT

// Define a template struct for efficiently operating on an upper triangular matrix
template <typename matrix_type>
struct UpperMatrix {
    
    // This macro calculates the index in the storage array for the element at position (x, y)
    #define UM_INDEXING(x,y) y + x * matrix_dimension - (x * x + x)/2

    // Pointer to an array that holds the matrix elements
    matrix_type* matrix_array;

    // Size of the data array, which is the number of elements in the upper triangular matrix
    const int array_size;

    // Dimension of the matrix (NxN)
    const int matrix_dimension;

    // Constructor for the UpperMatrix class
    // Initializes the matrix with a given dimension and default value
    UpperMatrix(const int matrix_dimension, matrix_type default_value):
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

    // Destructor for the UpperMatrix class
    // Frees the allocated memory for the matrix array
    ~UpperMatrix() {
        delete[] matrix_array;
    }

    // Accesses an element of the matrix at position (i, j)
    // Throws an exception if the indices are out of bounds or if (i, j) is below the diagonal
    matrix_type& element(int i, int j) {
        DEBUG_CHECK( 1,
            if (i < 0 || i >= matrix_dimension || j < 0 || j >= matrix_dimension ) {
                throw std::out_of_range("Index out of bounds");
            }
            else if (j < i) {
                throw std::out_of_range("Element not in the upper triangular");
            }
        )
        // Return the matrix element at the computed index
        return matrix_array[UM_INDEXING(i,j)];
    }

    // Overload the stream insertion operator to print the matrix
    // Prints the matrix in a human-readable format
    friend std::ostream & operator<<(std::ostream& os, UpperMatrix<matrix_type> & tc) {
        std::stringstream ret;
        std::string token(8, 'X');  // A placeholder string for elements below the diagonal, can be edited freely

        // Iterate over each element in the matrix
        for (int i = 0; i < tc.matrix_dimension; i++){
            for (int j = 0; j < tc.matrix_dimension; j++){
                if (i > j) 
                    ret << token << " ";  // Print placeholder for elements below the diagonal
                else
                    ret << tc.element(i,j) << " ";  // Print actual matrix element
                
            }
            ret << std::endl;  // Move to the next line after printing a row
        }
        return os << ret.str() << std::endl;  // Output the matrix as a string
    }
};



// FUNCTIONS

// Function to print the usage instructions for the program
void print_usage() {
    std::cout << "usage: fastflow_wf [options]" << std::endl
        << "\t-h\t\tprints this message" << std::endl
        << "\t-N matrix_size\t\tsize of the square matrix [default=" << DEFAULT_M_SIZE << "]" << std::endl
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
    new_file << "matrix_size,num_workers,machine_identifier,debug_value,init_time,running_time" <<std::endl;
    new_file.close(); // Close the new file after writing the header
    return true;
}



// MAIN

int main(int argc, char *argv[]) {   
    // Output debug information if debugging is enabled, level 2
    DEBUG_STDOUT(2, "Debug mode is active - DEBUG=" << DEBUG);

    // Default values for command-line arguments
    int matrix_size = DEFAULT_M_SIZE;   // Default matrix size
    int num_workers = DEFAULT_W;        // Default number of workers
    int machine_identifier = DEFAULT_MACHINE;  // Default machine identifier
    std::string log_file_name = DEFAULT_LOG_NAME; // Default log file name

    // Variable to store the option retrieved by getopt
    int opt;
    
    // Parse command-line options using getopt
    while ((opt = getopt(argc, argv, "hN:W:m:l:")) != -1) {
        switch (opt) {
            case 'm':
                // Set machine identifier
                machine_identifier = std::atoi(optarg);
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
                // Print usage instructions and exit
                print_usage();
                return 0;
            case '?':       
            default:
                // Handle unknown options or errors
                print_usage();
                return 1;
        }
    }

    // Validate matrix size
    if (matrix_size <= 1) {
        std::cout << "matrix_size needs to be > 1" << std::endl;
        return 1;
    }

    // Output debug information with level 3
    DEBUG_STDOUT(3,
                "Starting with parameter: N=" << matrix_size << " W=" << num_workers <<
                " m=" << machine_identifier << " l=" << log_file_name
            )

    // Initialize the log file
    if (!init_file(log_file_name)) {
        std::cout << "log file creation error" << std::endl;
        return 1;
    }
    // Output debug information with level 3 after file initialization
    DEBUG_STDOUT(3,
                "Check or creation of the log file complete"
            )

    // Variables to measure matrix initialization time
    double init_time;
    TIMERSTART(mat_init)
    
    // Initialize a single upper triangular matrix
    // um - represents an upper triangular matrix with column-major order
    UpperMatrix um(matrix_size, DEFAULT_INITIALISATION);
    for (int i = 0; i < matrix_size; i++) {
        um.element(i, i) = (i + 1.0) / matrix_size;  // Initialize diagonal elements
    }

    // Create the parallel computation object with the specified number of workers with fastflow
    ParallelForReduce<double> pf(num_workers);

    TIMERSTOP(mat_init, init_time)
    // Output debug information with level 3 after matrix initialization
    DEBUG_STDOUT(3,
                "Initialization complete - elapsed time: " << init_time
            )

    // Variables to measure computation time
    double running_time;
    TIMERSTART(run)
    
    // Perform the wavefront computation for each diagonal in parallel
    // The work is divided equally for each thread dividing the diagonal at each external iteration
    for (int diag = 1; diag < um.matrix_dimension; diag++) {
        long v = (um.matrix_dimension - diag);
        // There's no need to use a rediuction operator as we're going to write the results in memory. 
        // Also, the results on each diagonal are independent of each other.
        pf.parallel_for(0L, v, 1L,
                        0L,
                        [&](const long k) {
                            double partial_sum = 0.0;
                            // Output debug information with level 5 during summation
                            DEBUG_CHECK(5, std::cout << "[" << k << "," << diag + k << "]";)
                            for (int l = 0; l < diag; l++) {
                                DEBUG_CHECK(5, std::cout << " - (" << k << "," << l << ") + (" << diag + k - l << "," << diag + k << ")";)
                                // Compute partial sums using elements from the matrix
                                partial_sum += um.element(k, k + l) * um.element(diag + k - l, diag + k);
                            }
                            // Output debug information with level 5 after summation
                            DEBUG_CHECK(5, std::cout << " |" << std::endl;)
                            // Update matrix elements with the cube root of the computed partial sum
                            partial_sum = std::cbrt(partial_sum);
                            um.element(k, diag + k) = partial_sum;
                        });
    }
    TIMERSTOP(run, running_time)

    // Output debug information with level 5 showing the final matrix
    DEBUG_STDOUT(5,
        "Final Matrix row:" << std::endl << um
    )
    
    // Output debug information with level 2 showing the final value and computation time
    DEBUG_STDOUT(2, "Final value (matrix dim " << matrix_size << "): " << um.element(0, um.matrix_dimension - 1))
    DEBUG_STDOUT(2, "Computation complete - elapsed time: " << running_time * 1000)
    
    // Write execution times and other details to the specified log file
    std::ofstream file;
    file.open(log_file_name, std::ios_base::app);
    file << matrix_size << "," << num_workers << "," << machine_identifier << "," << DEBUG << ","
        << init_time * 1000 << "," << running_time * 1000 << std::endl;
    file.close();

    return 0;
}
