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
    sprintf(str_pid, "%d", pid);
    return str_pid;
}

void signal_handler( int signal_num ) {
	Rcpp::Rcout << "The interrupt signal is (" << signal_num << "). \n";
	// terminate program
	exit(signal_num);
}

