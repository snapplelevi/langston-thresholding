// RCpp wrapper of Carissa Bleker's main_threhsold C++ code
// made on dev
// Carissa BLeker
// cbleker@vols.utk.edu

#include <Rcpp.h>
#include <igraph.h>

#include <vector>     // std::vector
#include <iostream>   // std::cout, std::cerr, std::endl
#include <fstream>    // fopen, fclose (to read igraph), ofstream
#include <getopt.h>   // commandline argument parsing
#include <stdlib.h>   // atoi, atof
#include <sstream>    // stringstream
#include <string>     // getline

#include "utils.h"
#include "math_ext.h"
#include "igraph_ext.h"

#include "local_global.h"
#include "local_rank.h"

///////////////////////////////////////////////////////////////////////////////


// Manually exported in NAMESPACE
//' Strict thresholding for weighted graphs in .ncol format
//'
//' Description of threshold goes here
//'
//' @param infile DOCUMENT THIS
//' @param outfile DOCUMENT THIS
//' @param method DOCUMENT / CHANGE THIS
//' @param absolute DOCUMENT / CHANGE THIS
//' @param local_global_alpha DOCUMENT / CHANGE THIS
//' @param rank DOCUMENT / CHANGE THIS
//' @examples
//' print('NEED TO WRITE EXAMPLES FOR THIS LATER')
//' print('MULTIPLE IF POSSIBLE')
//' @returns Nothing. The thresholded graph is written to the file specified by outfile.
// [[Rcpp::export]]
int threshold(std::string& infile,
              std::string& outfile,
              std::string method,
              double absolute,
              double local_global_alpha,
              int rank){

    std::ifstream fin;
    fin.open(infile);
    if(fin.good() != 1){
        Rcpp::Rcerr << "The input file " << infile << " had a problem during opening.\n";
        Rcpp::Rcerr << "Please make sure the input file is spelled correctly and\n";
        Rcpp::Rcerr << "you have the proper permissions to read it.\n";
        Rcpp::stop(infile + " was not able to be opened");
    }

    // turn on attribute handling
    // for igraph to handle edge weights
    igraph_i_set_attribute_table(&igraph_cattribute_table);

    // Load graph (names = true)
    igraph_t G;
    Rcpp::Rcout << "Loading graph ... " << std::flush;
    read_graph(infile, G, IGRAPH_ADD_WEIGHTS_YES, true);
    Rcpp::Rcout << "Graph loaded." << "\n";
    Rcpp::Rcout << "--------------------------------------------------------------------\n";


    Rcpp::Rcout << "Original number vertices: " << igraph_vcount(&G) << "\n";
    Rcpp::Rcout << "Original number edges:    " << igraph_ecount(&G) << "\n";
    Rcpp::Rcout << "--------------------------------------------------------------------\n";


    if (method == "absolute"){
        threshold_graph(absolute, G);
        write_graph(outfile, G);
        Rcpp::Rcout << "Resulting number vertices: " << igraph_vcount(&G) << "\n";
        Rcpp::Rcout << "Resulting number edges:    " << igraph_ecount(&G) << "\n";
    }
    else if (method == "local-global"){
        igraph_t new_G;
        double mean_k;
        local_global_pruning(G, local_global_alpha, new_G, mean_k);
        write_graph(outfile, new_G);
        Rcpp::Rcout << "Resulting number vertices: " << igraph_vcount(&new_G) << "\n";
        Rcpp::Rcout << "Resulting number edges:    " << igraph_ecount(&new_G) << "\n";
    }
    else if (method=="rank"){
        igraph_t new_G;
        local_rank(G, rank, new_G);
        write_graph(outfile, new_G);
        Rcpp::Rcout << "Resulting number vertices: " << igraph_vcount(&new_G) << "\n";
        Rcpp::Rcout << "Resulting number edges:    " << igraph_ecount(&new_G) << "\n";
   }
   // TODO: refactor the the code to support a strict threshold instead of requiring that absolute
   // thresholding by default? Carissa's code requires that one of the above options are entered
   // to work properly. Is it meaningful to only threshold one side of a distribution?
    else {
        Rcpp::Rcerr << "The option: \"" << method << "\" is not supported.\n";
        Rcpp::Rcerr << "Exiting threshold() now...\n"; 
    }

    Rcpp::Rcout << "--------------------------------------------------------------------\n";
    Rcpp::Rcout << "Done! Thresholded graph saved to: " << outfile << "\n";
    Rcpp::Rcout << "--------------------------------------------------------------------\n";
    return 0;
}


///////////////////////////////////////////////////////////////////////////////
//     Command-line arguments                                                //
///////////////////////////////////////////////////////////////////////////////


// NOT EXPORTED
/*
void help(std::string prog_name){
    Rcpp::Rcerr <<  "\n";
    Rcpp::Rcerr <<  "    Usage: \n";
    Rcpp::Rcerr <<  "    " << prog_name     << " [-OPTIONS]... <GRAPH FILE PATH> <OUTPUT FILE PATH> \n\n";
    Rcpp::Rcerr <<  "    Graph has to be in .ncol format. \n";
    Rcpp::Rcerr <<  "    One of the following options have to be given: ";
    Rcpp::Rcerr <<  "   \n\n";
    Rcpp::Rcerr <<  "    Options: \n";
    Rcpp::Rcerr <<  "      -a  --absolute              <value>     Threshold graph at absolute of <value>\n";
    Rcpp::Rcerr <<  "      -l  --local-global          <value>     Use local-global method to threshold with alpha = <value>\n";
    Rcpp::Rcerr <<  "      -r  --rank                  <value>     Use top <value> ranked edges per vertex to threshold graph\n";
    Rcpp::Rcerr <<  "      -h  --help                              print this help and exit\n";
    Rcpp::Rcerr <<  "\n";
    exit(0);
}
*/

// NOT EXPORTED EITHER
/*
int argument_parser(int argc, char **argv,
    // Mandatory argument definitions
    std::string &infile,
    std::string &outfile,
    // Here flags (options without arguments) and arguments with defined type
    std::string &method,
    double &absolute,
    double &local_global_alpha,
    int &rank){


    const char* const short_options = "ha:l:u:r:" ;
    const struct option long_options[] =
        {    //name,                    has_arg,    flag,        val
            { "help",                   0,          NULL,        'h'},
            { "absolute",               1,          NULL,        'a'},
            { "local-global",           1,          NULL,        'l'},
            { "rank",                   1,          NULL,        'r'},
            { NULL, 0, NULL, 0 }
        };

    // Parse option
    int option;
    option = getopt_long(argc, argv,
                              short_options, long_options, NULL);

    switch (option){

        case 'h' : // -h or --help
            help(argv[0]);
            break;

        case 'a' : // -a or --absolute
            absolute=atof(optarg);
            method="absolute";
            break;

         case 'l' : // -l or --local-global
            local_global_alpha=atof(optarg);
            method="local-global";
            break;

         case 'r' : // -r or --rank
            rank=atoi(optarg);
            method="rank";
            break;

        case '?' : // Invalid option
            help(argv[0]); // Return help
            break;

        case -1 : // No more options
            break;

        default : // Something unexpected? Aborting
            return(1);
    }

    // Mandatory arguments
    // Current index (optind) < than the total number of arguments
    if(optind == argc){
        Rcpp::Rcerr << "\n Mandatory argument(s) missing\n";
        help(argv[0]);
    }
    // Iterate over rest of the arguments (i.e. in argv[optind])
    infile = argv[optind];
    optind++;
    outfile = argv[optind ];

    return 0;
}
*/
///////////////////////////////////////////////////////////////////////////////
//     Main                                                                  //
///////////////////////////////////////////////////////////////////////////////
/*
int main(int argc, char **argv){

    // Parse arguments
    // Mandatory argument definitions
    std::string infile;  //input file name
    std::string outfile;

    // Flags (options without arguments) and arguments with defined type
    std::string method;
    double absolute;
    double local_global_alpha;
    int    rank;

    argument_parser(argc, argv, infile, outfile,
                    method, absolute, local_global_alpha, rank);

    // check arguments
    if(outfile.empty()) {
        Rcpp::Rcerr << "No output file name specified. " << "\n";
        return 0;
    }

    ///////////////////////////////////////////////////////////////////////

    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "------------------------------------------------\n";
    Rcpp::Rcout << "input graph file:      "  << infile << "\n";
    Rcpp::Rcout << "output graph file      "  << outfile << "\n";
    Rcpp::Rcout << "threshold method       "  << method << "\n";
    Rcpp::Rcout << "------------------------------------------------\n";

    // turn on attribute handling
    // for igraph to handle edge weights
    igraph_i_set_attribute_table(&igraph_cattribute_table);

    // Load graph (names = true)
    igraph_t G;
    Rcpp::Rcout << "Loading graph ... " << std::flush;
    read_graph(infile, G, IGRAPH_ADD_WEIGHTS_YES, true);
    Rcpp::Rcout << "done." << "\n";

    int status;
    status = threshold(outfile, G, method, absolute, local_global_alpha, rank);
    return status;
} */

///////////////////////////////////////////////////////////////////////////////

