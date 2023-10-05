#include <Rcpp.h>
#include <igraph.h>

#include "math_ext.h"
#include "igraph_ext.h"

#include "random_matrix_theory.h"
#include "spectral_methods.h"
#include "scale_free.h"
#include "maximal_cliques.h"
#include "local_global.h"
#include "significance.h"
#include "local_rank.h"

#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <getopt.h>

#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <iomanip>


/////////////// update later after basic functionality of analysis program works //////////////
// Helper code to print out help information if inputs are incorrect

//' Display argument information on terminal for thresholding::thresholdAnalysis()
//' 
// [[Rcpp::export]]
void help(){
    /*
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
    */
    Rcpp::Rcout <<  "\n";
    Rcpp::Rcout <<  "--------------Langston Lab Thresholding Analysis Techniques (2023)---------------\n";
    Rcpp::Rcout <<  "                               thresholdAnalysis()                               \n";
    Rcpp::Rcout <<  "Synopsis:\n";
    Rcpp::Rcout <<  "    thresholdAnalysis(infile, outfile_prefix,\n";
    Rcpp::Rcout <<  "                      methods="",  lower=0.5,  upper=0.99, \n";
    Rcpp::Rcout <<  "                      increment=0.01,  window_size=5,  min_partition_size=10, \n";
    Rcpp::Rcout <<  "                      min_clique_size=5,  min_alpha=0,  max_alpha=4,\n";
    Rcpp::Rcout <<  "                      alpha_increment=0.1,  num_samples=0,\n";
    Rcpp::Rcout <<  "                      significance_alpha=0.01,  bonferroni_corrected=0)\n";
    Rcpp::Rcout <<  "\n";
    Rcpp::Rcout <<  "Arguments to thresholdAnalysis():\n";
    Rcpp::Rcout <<  "Required:\n";
    Rcpp::Rcout <<  "\t1.) infile: string input\n";
    Rcpp::Rcout <<  "\t        The weighted edge list (.wel) file input. This file is in the .ncol format as specified by\n";
    Rcpp::Rcout <<  "\t        the Large Graph Layout group: https://lgl.sourceforge.net/#FileFormat.\n\n";
    
    Rcpp::Rcout <<  "\t        In this application, the graph input file is simple, weighted, undirected. The vertices in the .wel\n";
    Rcpp::Rcout <<  "\t        file follow the following format where the arrow represents whitespace: \n";
    Rcpp::Rcout <<  "\t          vertex1⇥vertex2⇥weight1,2\n";
    Rcpp::Rcout <<  "\t          vertex1⇥vertex3⇥weight1,3\n";
    Rcpp::Rcout <<  "\t          ...\n\n";

    Rcpp::Rcout <<  "\t         NOTE: vertex names cannot contain whitespace.\n\n";

    Rcpp::Rcout <<  "\t2.) outfile_prefix: string input\n";
    Rcpp::Rcout <<  "\t        Prefix to the output files to the analysis file(s).\n";
    Rcpp::Rcout <<  "\t        Example: If the prefix is \"graph-output\", the output is graph-output.<method_name>.txt";
    Rcpp::Rcout <<  "\t                 where method name is the type of analysis performed, such as iterative or statistical_errors\n";
    Rcpp::Rcout <<  "\t                 The method(s) are controlled by the optional \"methods\" argument.\n\n";

    Rcpp::Rcout <<  "Optional:\n";
    Rcpp::Rcout <<  "\t3.) methods: string input  (defaults to empty string)\n";
    Rcpp::Rcout <<  "\t        Comma separated list of analysis operations to complete. These methods are represented by an integer";
    Rcpp::Rcout <<  "\t        which is mapped to its corresponding method internally. The following methods are currently available:\n\n";
    Rcpp::Rcout <<  "\t             0 - all (methods 1-7 will be performed)\n";
    Rcpp::Rcout <<  "\t             1 - significance and power calculations (only valid for Pearson CC)\n";
    Rcpp::Rcout <<  "\t             2 - local-global\n";
    Rcpp::Rcout <<  "\t             3 - scale free\n";
    Rcpp::Rcout <<  "\t             4 - maximal cliques\n";
    Rcpp::Rcout <<  "\t             5 - spectral methods\n";
    Rcpp::Rcout <<  "\t             6 - random matrix theory\n";
    Rcpp::Rcout <<  "\t             7 - clustering coefficient\n";
    Rcpp::Rcout <<  "\t             8 - percolation\n\n";
    Rcpp::Rcout <<  "\t        Example: methods=\"7, 5, 2\" (Note: methods will be performed in numerical order internally, but the order";
    Rcpp::Rcout <<  "\t                 which they are passed to the function doesn't matter.)\n\n";
    Rcpp::Rcout <<  "\t4.) lower: floating point input  (defaults to 0.5)\n";
    Rcpp::Rcout <<  "\t        Lower initial starting bound for the thresholding loop. The loop ends when the upper bound limit.\n";
    Rcpp::Rcout <<  "\t        Note: lower must be less than or equal to upper (lower <= upper) for function to continue.\n\n";
    Rcpp::Rcout <<  "\t        Example: lower=0.6";
}


// Internal function -- add character or expanded keywords from a string (methods)
// and insert the method into a set (in case the user passes more than one of the 
// same argument)
//
// Invalid arguments will cause the program to exit with an error and alert the user
// of the incorrect argument passed and output information about the types of methods
// that are allowed to be passed to the function

// TODO: strip white space from beginning and end of str_methods
void parse_string_methods(std::set<int> &analysis_methods, const std::string &str_methods){
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
void thresholdAnalysis(std::string infile, 
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
    // Stores the outfile name passed to analysis functions at multiple points throughout 
    // the analysis exeuction
    std::string outfile_name;

    // Ensure output file prefix exists
    // Return one for 
    if(outfile_prefix.empty()){
        Rcpp::Rcerr << "Error - No output file prefix specified." << '\n';
        Rcpp::Rcerr << "Use thresholding::help() to get more information on thresholdAnalysis() inputs." << '\n';
        Rcpp::stop("empty output file prefix. Ending analysis early.");
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

    // Ensure lower does not exceed the value of upper
    if(upper <= lower){
        Rcpp::Rcerr << "Error in threshold limits: ";
        Rcpp::Rcerr << "cannot have lower >= upper.\n";
        Rcpp::Rcerr << "Please restart with corrected lower and upper bounds." << '\n';
        Rcpp::Rcerr << "Use thresholding::help() to get more information on thresholdAnalysis() inputs." << '\n';
        Rcpp::stop("invalid upper and lower limits. Ending analysis early.");
    }

    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "Initial analysis parameters" << '\n';
    Rcpp::Rcout << "------------------------------------------------\n";
    Rcpp::Rcout << "input graph file:      "  << infile << "\n";
    Rcpp::Rcout << "output file prefix:    "  << outfile_prefix << "\n";
    Rcpp::Rcout << "lower threshold:       "  << lower << "\n";
    Rcpp::Rcout << "upper threshold:       "  << upper << "\n";
    Rcpp::Rcout << "threshold increment:   "  << increment << "\n";
    Rcpp::Rcout << "------------------------------------------------\n";


    // Ensure that the analysis_methods are valid before continuing the analysis process
    // and put the methods into the set for later analysis operations
    std::set<int> analysis_methods;
    parse_string_methods(analysis_methods, methods);

    /////////////////////////////////////////////////////////////////////////////////////////
    // 1 = Method for finding significance and power calculations (only valid for Pearson CC)
    // Type I error (false positive rate) and
    // Type II error (false negative rate) control
    // Have to have n - number of samples (not number of variables)
    if(analysis_methods.find(1) != analysis_methods.end()){
        outfile_name = outfile_prefix + ".statistical_errors.txt";
        control_statistical_errors(significance_alpha,
                                  num_samples,
                                  0, //E
                                  bonferroni_corrected,
                                  outfile_name);
        analysis_methods.erase(1);
    }

    // End program if the only method desired was the significance and power calculations
    if (analysis_methods.size() == 0){
        return;
    }

    // Turn on attribute handling
    // For igraph to handle edge weights
    igraph_set_attribute_table(&igraph_cattribute_table);

    // Hold the result of the thresholding analysis process
    int status = 0;
    
    // Load graph   
    // Ensure path to infile containing graph info is valid
    // read_graph checks to make sure that file opened is an existing file and passes
    // the file to the igraph_read_graph_ncol function
    igraph_t G;
    Rcpp::Rcout << "Loading graph..." << '\n';
    read_graph(infile, G, IGRAPH_ADD_WEIGHTS_YES);
    Rcpp::Rcout << "Done! Graph has been loaded." << '\n';
    
    igraph_integer_t E = igraph_ecount(&G); // number edges
    igraph_integer_t V = igraph_vcount(&G); // number vertices

    // Max number edges based on number of vertices
    // i.e. the number of edges in a complete graph with V vertices
    double orig_max_E = 0.5 * V * (V - 1.0);
    
    Rcpp::Rcout << "Number vertices:  " << V << "\n";
    Rcpp::Rcout << "Number edges:     " << E;
    Rcpp::Rcout << "  (maximum possible number edges " << int(orig_max_E) << ")\n";
    Rcpp::Rcout << "------------------------------------------------\n\n";

    ///////////////////////////////////////////////////////////////////////
    // Non-loop methods
    ///////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////
    // local-global (guzzi2014, rank)
    if(analysis_methods.find(2) != analysis_methods.end()){
        outfile_name = outfile_prefix + ".local_global.txt";
        local_global_method(G,
                     min_alpha,
                     max_alpha,
                     alpha_increment,
                     window_size,
                     min_partition_size,
                     outfile_name);
        analysis_methods.erase(2);
    }

    // Exit program if no only non-loop methods were requested
    if (analysis_methods.size() == 0){
        return;
    }

    ///////////////////////////////////////////////////////////////////////
    // Threshold loop
    // loop destroys the graph
    ///////////////////////////////////////////////////////////////////////

    // Ready the output file
    outfile_name = outfile_prefix + ".iterative.txt";
    std::ofstream out;
    out.open(outfile_name.c_str(), std::ofstream::out);
    // Is the iterative file open for writing?
    if (out.fail()) {
        Rcpp::Rcerr << "Error opening file for writing: " << outfile_name << "\n";
        Rcpp::stop("output file was unable to be opened for writing. Ending analysis early.");
    }

    // Output header for iterative file contents
    std::stringstream header;
    header << "threshold";
    header << "\tvertex-count\tedge-count";
    header << "\tconnected-component-count";
    header << "\tdensity\tdensity-orig-V";
    header << "\tlargest-cc-size\t2nd-largest-cc-size";
    header << "\tclustering-coefficient\trandom-clustering-coefficient";
    header << "\t2nd-eigenvalue\talmost-disconnected-component-count";
    header << "\tmaximal-clique-count\tclique-number";
    header << "\tpoisson-chi2\tpoisson-pvalue";
    header << "\tgoe-chi2\tgoe-pvalue";
    header << "\tscale-free-KS\tscale-free-KS-p-value\tscale-free-alpha";
    out << header.str();
    out << '\n';

    // Get the threshold increments - range() function from math_ext
    double t;
    static const std::vector<double> t_vector = range(lower, upper, increment);
    int num_increments = t_vector.size();

    Rcpp::Rcout << "Iterative thresholding\n";
    Rcpp::Rcout << "Number steps: " << num_increments << '\n';

    // Initialise necessary stuff
    int nearly_disconnected_components      = -1;
    igraph_real_t  second_eigenvalue        = std::nan("");

    igraph_integer_t    clique_count        = -1;   // number maximal cliques
    igraph_integer_t    clique_number       = -1;   // maximum clique size

    double  density                         = std::nan("");
    double  density_orig_V                  = std::nan("");

    double  poi_chi_sq_stat                 = std::nan("");
    double  goe_chi_sq_stat                 = std::nan("");

    double  poi_chi_sq_pvalue               = std::nan("");
    double  goe_chi_sq_pvalue               = std::nan("");

    igraph_integer_t     cc_count           = -1;
    igraph_integer_t     largest_cc_size    = -1;
    igraph_integer_t     largest2_cc_size   = -1;

    double  scale_free_pvalue               = std::nan("");
    double  scale_free_KS                   = std::nan("");
    double  scale_free_xmin                 = std::nan("");
    double  scale_free_alpha                = std::nan("");

    igraph_real_t clustering_coefficient    = std::nan("");
    igraph_real_t clustering_coefficient_r  = std::nan("");

    ///////////////////////////////////////////////////////////////////////
    // If methods includes clustering coefficient,
    // then need a copy of G to threshold
    // Rewiring at each threshold takes longer (?)
    igraph_t G_random;
    if (analysis_methods.find(7) != analysis_methods.end()){
        igraph_copy(&G_random, &G);

        // igraph rewire loses edge weights,
        // need to save them and reassign back ()
        igraph_rewire(&G_random, E, IGRAPH_REWIRING_SIMPLE);

        igraph_vector_t edge_weights;
        igraph_vector_init(&edge_weights, E);

        igraph_cattribute_EANV(&G, "weight", igraph_ess_all(IGRAPH_EDGEORDER_ID), &edge_weights);
        igraph_cattribute_EAN_setv(&G_random, "weight", &edge_weights);

        igraph_vector_destroy(&edge_weights);
    }
  
    // Main threshold loop
    for(int i_t = 0; i_t < num_increments; i_t++){
        t = t_vector[i_t];

        Rcpp::Rcout << "\nStep: " << i_t +1  << ", Threshold: " << t << '\n';

        // Threshold step
        int threshold_status = threshold_graph(t, G);
        V = igraph_vcount(&G);
        E = igraph_ecount(&G);


        if(threshold_status == 1){
            // 1 = all edges removed, stop
            Rcpp::Rcout <<" Graph is empty, finished. " << '\n';
            break;
        }
        else if( threshold_status == 3){
            // 3 = no edges removed, only skip if not first iteration
            if(i_t > 0){
                Rcpp::Rcout << " No edges removed, skipping. " << '\n';
                continue;
            }
            else{
                Rcpp::Rcout << "\t\tVertices: " << V << "\tEdges: " << E;
                Rcpp::Rcout << '\n';
            }
        }
        else if(threshold_status == 2){
            // 2 = some edges are removed, keep going
            // make sure graph is large enough to continue
            if ( (V < min_partition_size) || (E < min_partition_size) ){
                //not large enough
                Rcpp::Rcout <<" Graph too small, finished. " << '\n';
                break;
            }
            else{
                Rcpp::Rcout << "\t\tVertices: " << V << "\tEdges: " << E;
                Rcpp::Rcout << '\n';
            }
        }
        else{
            Rcpp::Rcerr << " Something went wrong during thresholding loop" << '\n';
            Rcpp::stop("weirdness occured during thresholding loop...");
        }

        ///////////////////////////////////////////////////////////////////////
        // Metrics to do by default
        ///////////////////////////////////////////////////////////////////////

        // Density
        density        = (double) E / ( 0.5 * V * (V -1) );
        density_orig_V = (double) E / orig_max_E;


        ///////////////////////////////////////////////////////////////////////
        // Metrics to only do if requested
        ///////////////////////////////////////////////////////////////////////

        for (auto& m : analysis_methods) {
                ///////////////////////////////////////////////////////////////
                // None
                if(m==-1){
                    break;
                }

                else if(m==5 || m==8){
                    ///////////////////////////////////////////////////////////////////////
                    //  Largest connected component sizes
                    // (basically percolation)
                    igraph_t G_cc;
                    largest_connected_component(G, G_cc, cc_count,
                        largest_cc_size, largest2_cc_size);

                    // Spectral Methods
                    if(m == 5){
                        if(largest_cc_size >= min_partition_size){
                            spectral_methods(G_cc,
                               window_size,
                               min_partition_size,
                               second_eigenvalue,
                               nearly_disconnected_components);
                        }
                    }
                    igraph_destroy(&G_cc);
                }

                ///////////////////////////////////////////////////////////////
                // Maximal Clique Number
                else if(m==4){
                    maximal_cliques(G, min_clique_size,
                                    clique_count, clique_number);
                }

                ///////////////////////////////////////////////////////////////
                // Scale free
                else if(m==3){
                    igraph_plfit_result_t scale_free_result;
                    scale_free_test(G, V, scale_free_result);
                    scale_free_pvalue = scale_free_result.p;
                    scale_free_KS = scale_free_result.D;
                    scale_free_xmin = scale_free_result.xmin;
                    scale_free_alpha = scale_free_result.alpha;
                }

                ///////////////////////////////////////////////////////////////
                // Random Matrix Theory
                else if(m==6){
                    random_matrix_theory(G,
                                         V,
                                         poi_chi_sq_stat,
                                         goe_chi_sq_stat,
                                         poi_chi_sq_pvalue,
                                         goe_chi_sq_pvalue);
                }

                ///////////////////////////////////////////////////////////////
                // Clustering coefficient
                else if(m==7){
                    // plain clustering coefficient
                    igraph_transitivity_undirected(&G, &clustering_coefficient, IGRAPH_TRANSITIVITY_NAN);

                    // threshold G_random an d get random clustering coefficient
                    int threshold_status2 = threshold_graph(t, G_random);
                    //std::cout << "\n" <<  threshold_status2 << "   " << igraph_vcount(&G_random) << "  " << igraph_ecount(&G_random) << std::endl;
                    igraph_transitivity_undirected(&G_random, &clustering_coefficient_r, IGRAPH_TRANSITIVITY_NAN);
                }

                ///////////////////////////////////////////////////////////////
                else{
                    //std::cerr << "Unknown method " << m << std::endl;
                    Rcpp::Rcerr << "Unknown method: " << m << ". Continuing analysis...\n";
                    continue;
                }
        }
        ///////////////////////////////////////////////////////////////////////
        // Make results into a string
        std::stringstream message;
        message << t;
        message << "\t" << V                 << "\t" << E;
        message << "\t" << cc_count;
        message << "\t" << density          << "\t" << density_orig_V;
        message << "\t" << largest_cc_size   << "\t" << largest2_cc_size;
        message << "\t" << clustering_coefficient;
        message << "\t" << clustering_coefficient_r;
        message << "\t" << second_eigenvalue;
        message << "\t" << nearly_disconnected_components;
        message << "\t" << clique_count      << "\t" << clique_number;
        message << "\t" << poi_chi_sq_stat   << "\t" << poi_chi_sq_pvalue;
        message << "\t" << goe_chi_sq_stat   << "\t" << goe_chi_sq_pvalue;
        message << "\t" << scale_free_KS     << "\t" << scale_free_pvalue;
        message << "\t" << scale_free_alpha;
        out << message.str();
        out << '\n';
    }

    if (analysis_methods.find(7) != analysis_methods.end()){
        igraph_destroy(&G_random);
    }
  
    out.close();
    Rcpp::Rcout << "Thresholding completed - all analysis methods completed.\n";

    igraph_destroy(&G);
}

