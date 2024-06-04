 #include <Rcpp.h>
 #include "igraph_ext.h"

///////////////////////////////////////////////////////////////////////////////
//     IO functions                                                          //
///////////////////////////////////////////////////////////////////////////////

// Read in graph
int read_graph(std::string& graph_file_path,
               igraph_t& G,
               igraph_add_weights_t is_weighted,
               igraph_bool_t names /*default false*/){

    FILE *graph_file;
    graph_file = fopen(graph_file_path.c_str(), "r");

    // Test if file exists
    if (graph_file == NULL){
        Rcpp::Rcerr << "Error - Unable to open file: " << graph_file_path << '\n';
        Rcpp::stop("Input file does not exist. Ending analysis early.");
        // exit(-1); --> R will crash when this is called. Rcpp::stop will handle program execution
    }

    // Read in file as graph
    igraph_read_graph_ncol(&G,
                           graph_file,
                           NULL,
                           names,
                           is_weighted,
                           IGRAPH_UNDIRECTED);

    fclose(graph_file);

    return 0;
}


// Write graph
int write_graph(std::string& graph_file_path,
                igraph_t& G){

    FILE *graph_file;
    graph_file = fopen(graph_file_path.c_str(), "w");

    if(graph_file == NULL){
        Rcpp::Rcerr << "Error - Unable to open graph output file: " << graph_file_path << '\n';
        Rcpp::stop("Output file does not exist. Ending operaion early.");
        //exit(-1); --> R will crash when this is called. Rcpp::stop will handle program execution
    }
    // write in file as weighted edge list
    igraph_write_graph_ncol(&G, graph_file, "name", "weight");

    fclose(graph_file);

    return 0;
}

