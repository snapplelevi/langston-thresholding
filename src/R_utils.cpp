#include "utils.h"

// incrementally write to an already stream
int output_results(std::ofstream outfile, std::string& message){
    outfile << message;
    return 0;
}

// https://stackoverflow.com/a/53230284/4996681
std::string get_str_pid(){
    int pid = getpid();
    char str_pid[20];   // ex. 34567
    snprintf(str_pid, sizeof(str_pif), "%d", pid);
    return str_pid;
}

void signal_handler( int signal_num ) {
	Rcpp::Rcout << "The interrupt signal is (" << signal_num << "). \n";
	// terminate program
    Rcpp::stop("exiting the process with signal: " + std::to_string(signal_num));
}

