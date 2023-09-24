#include <Rcpp.h>
#include <igraph.h>

#include "math_ext.h"
#include "spectral_methods.h"
#include "igraph_ext.h"

#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <getopt.h>

#include <cmath>
#include <algorithm>
#include <vector>
#include <set>




/////////////// update later after basic functionality of analysis program works //////////////
// Helper code to print out help information if inputs are incorrect

//' Display argument information on terminal for thresholding::thresholdAnalysis()
//' 
// [[Rcpp::export]]
void help(){
    Rcpp::Rcerr <<  "\n";
    Rcpp::Rcerr <<  "    Usage: \n";
    Rcpp::Rcerr <<  "    " << "thresholdAnaylsis"  << " [-OPTIONS]... <GRAPH FILE PATH> <OUTPUT FILE PATH> \n\n";
    Rcpp::Rcerr <<  "    Graph has to be in .ncol format. \n";
    Rcpp::Rcerr <<  "    Output file path is the prefix to the results files, which will be of the form: \n";
    Rcpp::Rcerr <<  "        <OUTPUT FILE PATH>.pid.<method_name>.txt\n\n";
    Rcpp::Rcerr <<  "    Options: \n";
    Rcpp::Rcerr <<  "      -l  --lower                  <value>     lower bound on thresholds to test (default 0.5)\n";
    Rcpp::Rcerr <<  "      -u  --upper                  <value>     upper bound on thresholds to test (default 0.99)\n";
    Rcpp::Rcerr <<  "      -i  --increment              <value>     threshold increment (default 0.01)\n";
    Rcpp::Rcerr <<  "      -w  --windowsize             <value>     sliding window size for spectral method (default 5)\n";
    Rcpp::Rcerr <<  "      -p  --minimumpartitionsize   <value>     minimum size of graph or subgraph after threshold (default 10)\n";
    Rcpp::Rcerr <<  "      -n  --num_samples            <value>     number of samples for significance and power calculations (default NULL)\n";
    Rcpp::Rcerr <<  "      -b  --bonferroni_correction              switch to perform bonferroni corrections in significance and power calculations (default FALSE)\n";
    Rcpp::Rcerr <<  "      -c  --minimum_cliquesize     <value>     minimum size of maximal cliques in maximal clique count (default 5)\n";
    Rcpp::Rcerr <<  "      -m  --methods                <value>     comma separated list of methods (defaults to none)\n";
    Rcpp::Rcerr <<  "                                                   0 - all\n";
    Rcpp::Rcerr <<  "                                                   1 - significance and power calculations (only valid for Pearson CC)\n";
    Rcpp::Rcerr <<  "                                                   2 - local-global\n";
    Rcpp::Rcerr <<  "                                                   3 - scale free\n";
    Rcpp::Rcerr <<  "                                                   4 - maximal cliques\n";
    Rcpp::Rcerr <<  "                                                   5 - spectral methods\n";
    Rcpp::Rcerr <<  "                                                   6 - random matrix theory\n";
    Rcpp::Rcerr <<  "                                                   7 - clustering coefficient\n";
    Rcpp::Rcerr <<  "                                                   8 - percolation\n";
    Rcpp::Rcerr <<  "      -h  --help                               print this help and exit\n";
    Rcpp::Rcerr <<  "\n";
    exit(0);
}

/*
int argument_parser(std::string &infile, 
                    std::string &outfile_prefix,
                    std::string methods, 
                    double &lower,
                    double &upper,
                    double &increment,
                    int &window_size,
                    int &min_partition_size,
                    int &min_clique_size,
                    double &min_alpha,
                    double &max_alpha,
                    double &alpha_increment,
                    int &num_samples,
                    double &significance_alpha,
                    bool &bonferroni_corrected)
{
    int next_option;

    const char* const short_options = "hl:u:i:w:p:n:bc:m:" ;
    const struct option long_options[] =
        {    //name,                    has_arg,    flag,        val
            { "help",                   0,          NULL,        'h'},
            { "lower",                  1,          NULL,        'l'},
            { "upper",                  1,          NULL,        'u'},
            { "increment",              1,          NULL,        'i'},
            { "windowsize",             1,          NULL,        'w'},
            { "minimumpartitionsize",   1,          NULL,        'p'},
            { "num_samples",            1,          NULL,        'n'},
            { "bonferroni_correction",  0,          NULL,        'b'},
            { "minimum_cliquesize",     1,          NULL,        'c'},
            { "methods",                1,          NULL,        'm'},
            { NULL, 0, NULL, 0 }
        };

    // Parse options
    while (1) {
        // Obtain an option
        next_option = getopt_long(argc, argv,
                                  short_options, long_options, NULL);

        if (next_option == -1)
            break; // No more options. Break loop.

        switch (next_option){

            case 'h' : // -h or --help
                help(argv[0]);
                break;

            case 'l' : // -l or --lower
                lower=atof(optarg);
                break;

             case 'u' : // -u or --upper
                upper=atof(optarg);
                break;

             case 'i' : // -i or --increment
                increment=atof(optarg);
                break;

            case 'w' : // -w or --windowsize
                windowsize=atoi(optarg);
                break;

            case 'p' : // -p or --minimumpartitionsize
                minimumpartitionsize=atoi(optarg);
                break;

            case 'n' : // -p or --minimumpartitionsize
                num_samples=atoi(optarg);
                break;

            case 'b' : // -p or --minimumpartitionsize
                bonferroni_corrected=1;
                break;

            case 'c' : // -p or --minimumpartitionsize
                minimum_cliquesize=atoi(optarg);
                break;

            case 'm' : // -m or --methods
                methods=optarg;
                break;

            case '?' : // Invalid option
                help(argv[0]); // Return help

            case -1 : // No more options
                break;

            default : // Something unexpected? Aborting
                return(1);
        }
    }

    // Mandatory arguments
    // Current index (optind) < than the total number of arguments
    //std::cout<<argc<<" number of arguments\n";
    //std::cout<<"optind: "<<optind<<"\n";
    if(optind == argc){
        Rcpp::Rcerr << "\n Mandatory argument(s) missing\n";
        help();
    }
    // Iterate over rest of the arguments (i.e. in argv[optind])
    if (argc - optind != 2){
         Rcpp::Rcerr << "\n Mandatory argument(s) missing\n";
        help();
    }
    infile = argv[optind];
    optind++;
    outfile = argv[optind];

    return 0;
} */

// Internal function -- add character or expanded keywords from a string (methods)
// and insert the method into a set (in case the user passes more than one of the 
// same argument)
//
// Invalid arguments will cause the program to exit with an error and alert the user
// of the incorrect argument passed and output information about the types of methods
// that are allowed to be passed to the function
int parse_string_methods(std::set<int> &analysis_methods, const std::string &str_methods){
    if(str_methods == "" ){
        analysis_methods.insert({-1});
    }
    else if (str_methods == "0" ){
        analysis_methods.insert({1, 2, 3, 4, 5, 6, 7});
    }
    else{
        std::istringstream sin(str_methods);
        std::vector<int> int_methods;
        std::string cur_method;

        while(std::getline(sin, cur_method, ',')){
            int_methods.push_back(std::stoi(cur_method));
        }
        analysis_methods.insert(int_methods.begin(), int_methods.end());
    }
    
}

// Manually exported in NAMESPACE
//' Main graph thresholding analysis function
//' @param infile: Name of .ncol graph file to read in for analysis
//' @param outfile_prefix: Prefix of output file in which analysis will be redirected to (Ex: <PREFIX>.iterative.txt )
//' @param methods: Comma separated list of analysis methods, listed if thresholding::help() is called (defaults to none)
// [[Rcpp::export]]
int thresholdAnalysis(std::string infile, 
                      std::string outfile_prefix,
                      std::string methods="", 
                      double lower=0.5,
                      double upper=0.99,
                      double increment=0.01,
                      int window_size=5,
                      int min_partition_size=10,
                      int min_clique_size=5,
                      double min_alpha=0,
                      double max_alpha=4,
                      double alpha_increment=0.1,
                      int num_samples=0,
                      double significance_alpha=0.01,
                      bool bonferroni_corrected=0)
{
    Rcpp::Rcout << "TESTING 1,2,3 THIS IS THE OPTIONAL PARAMETER lower: " << lower << '\n';

    // NOTE TO SELF: figure out a way to parse the string methods into something useful
    // similar to how it did inside of argument_parser, but without 
    /* check_args(infile, outfile_prefix, methods,
                lower, upper, increment, window_size,
                min_partition_size, min_clique_size, 
                min_alpha, max_alpha, alpha_increment, 
                num_sampes, significance_alpha, bonferroni_corrected); */

    // Ensure output file prefix exists
    if(outfile_prefix.empty()){
        Rcpp::Rcerr << "Error - No output prefix specified." << '\n';
        Rcpp::Rcerr << "Use thresholding::help() to get more information on thresholdAnalysis() inputs." << '\n';
        return 1;
    }

    // Ensure the window size is less than the minimum partition size
    if(min_partition_size <= window_size){
        Rcpp::Rcerr << "Warning: cannot have ";
        Rcpp::Rcerr << "min_partition_size <= windowsize. \n";
        Rcpp::Rcerr << "Proceeding using default values of ";
        Rcpp::Rcerr << "window_size = 5 and min_partition_size = 10";
        Rcpp::Rcerr << '\n';
        window_size = 5;
        min_partition_size = 10;
    }

    /*  Note to self:
            Only include this depending on the type of weights provided by the .wel
            0 < t < 1 is what is used for the correlation values. Look at the output of the 
            histogram program to see what the range of upper and lower can be. Delete these
            checks if the recommended upper and lower exceed 1.0.
    // Ensure thresholding bounds of upper and lower are within acceptable bounds
    if(lower < 0){
        Rcpp::Rcerr << "Error in thresholding limits:";
        Rcpp::Rcerr << "cannot have a lower bound below 0 (lower < 0).\n";
        Rcpp::Rcerr << "Please restart with a lower bound greater than or equal to 0.\n";
        Rcpp::Rcerr << "If a lower bound isn't specified to the function, ";
        Rcpp::Rcerr << "the default lower bound will be 0.5." << '\n';
        return 1;
    }
    else if(upper > 1.0){
        Rcpp::Rcerr << "Error in thresholding limits:";
        Rcpp::Rcerr << "cannot have an upper bound greater than 0 (upper > 1.0).\n";
        Rcpp::Rcerr << "Please restart with an upper bound less than or equal to 1.\n";
        Rcpp::Rcerr << "If an upper bound isn't specified to the function, ";
        Rcpp::Rcerr << "the default upper bound will be 0.99." << '\n';
        return 1;
    }

    */

    // Ensure lower does not exceed the value of upper
    if(upper <= lower){
        Rcpp::Rcerr << "Error in threshold limits: ";
        Rcpp::Rcerr << "cannot have lower >= upper.\n";
        Rcpp::Rcerr << "Please restart with corrected lower and upper bounds." << '\n';
        Rcpp::Rcerr << "Use thresholding::help() to get more information on thresholdAnalysis() inputs." << '\n';
        return 1;
    }


    // check that output file path exists
    std::ofstream out;
    out.open(outfile_prefix.c_str(), std::ofstream::out);
    // is it open
    if (out.fail()) {
        Rcpp::Rcerr << "Error opening file for writing: " << outfile_prefix << "\n";
        return 0;
    }

    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "------------------------------------------------\n";
    Rcpp::Rcout << "input graph file:      "  << infile << "\n";
    Rcpp::Rcout << "output file prefix:    "  << outfile_prefix << "\n";
    Rcpp::Rcout << "lower threshold:       "  << lower << "\n";
    Rcpp::Rcout << "upper threshold:       "  << upper << "\n";
    Rcpp::Rcout << "threshold increment:   "  << increment << "\n";
    Rcpp::Rcout << "------------------------------------------------\n";

    // Hold the result of the thresholding analysis process
    int status = 0;

    
    igraph_t G;

    // Ensure path to infile containing graph info is valid
    // read_graph checks to make sure that file opened is an existing file and passes
    // the file to the igraph_read_graph_ncol function

    // NOTE FOR FUTURE 9-20-23 work:
        // R session crashes likely due to the std::cout in igraph_ext.cpp.
        // resolve this so function doesn't crash when this is run
    Rcpp::Rcout << "Loading graph..." << '\n';
    read_graph(infile, G, IGRAPH_ADD_WEIGHTS_YES);
    Rcpp::Rcout << "Done! Graph has been loaded." << '\n';
    
    // Perform graph thresholding analysis
    igraph_integer_t E = igraph_ecount(&G); // number edges
    igraph_integer_t V = igraph_vcount(&G); // number vertices

    // Max number edges based on number of vertices
    // i.e. the number of edges in a complete graph with V vertices
    double orig_max_E = 0.5 * V * (V - 1.0);
    
    Rcpp::Rcout << "Number vertices:  " << V << "\n";
    Rcpp::Rcout << "Number edges:     " << E;
    Rcpp::Rcout << "  (maximum possible number edges " << int(orig_max_E) << ")";
    Rcpp::Rcout << '\n';
    Rcpp::Rcout << "------------------------------------------------\n\n";


    //std::set<int> &analysis_methods;
    /*
    parse_string_methods(analysis_methods, methods);

    if(methods.find(2)!=methods.end()){
        outfile_name = outfile_prefix + "local_global.txt";
        local_global_method(G,
                     min_alpha,
                     max_alpha,
                     alpha_increment,
                     windowsize,
                     minimumpartitionsize,
                     outfile_name);
        methods.erase(2);
    }

    if (methods.size() == 0){
        return 0;
    } */

    igraph_destroy(&G);
    return status;
}

