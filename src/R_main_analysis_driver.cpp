// R_main_analysis_driver.cpp
// Langston Lab Thresholding Methods
// September 2024

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
#include "utils.h"

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

// Internal function to convert list of methods in Rcpp::NumericVector into a
// std::set for use in the original analysis function
std::set<int> parse_methods_list(Rcpp::NumericVector methods){
    // set to return to the analysis() function
    std::set<int> retmethods;

    // Methods list was empty, so add -1 to the methods set to let
    // analysis know that no methods were requested (this is default behavior)
    // if the user doesn't specify any methods to analysis()
    if(methods.size() == 0){
        retmethods.insert(-1);
        return retmethods;
    }

    // Want to make sure method integers are in sorted order so if they enter
    methods = methods.sort();
    for(int i = 0; i < methods.size(); i++){

        // Currently only integers 2-8 are accepted as input for additional analysis methods
        // File output naming scheme attaches a concatenated string of all analysis methods
        // to the prefix (either default or specified by user).
        // Previously, adding a method int that wasn't valid would just add to the string
        // and nothing would happen in the main thresholding loop. 
        //
        // This means adding two invalid methods 23 and 45 would lead to the same file name
        // as specifying methods 2,3,4 and 5. However, it would be misleading as the analysis
        // call with 23 and 45 wouldn't have the desired outputs for methods 1-4. 
        // Therefore, this validates input before entering the main analysis loop.
        // 
        // If the user doesn't want any additinal analysis methods ran, then the methods list remains
        // empty and normal iterative thresholding is performed.
        int lower_method_bound = 2;
        int upper_method_bound = 8;

        // User wants all methods to be used (IF ONE IS SEEN: THEN ALL METHODS ARE USED)
        if(methods[i] == 1){
            retmethods.insert({2, 3, 4, 5, 6, 7, 8});
            break;
        }

        else if(methods[i] < lower_method_bound || methods[i] > upper_method_bound){
            Rcpp::Rcerr << "\nMethod #" << methods[i] << " is not a valid analysis method.\n";
            Rcpp::stop("invalid analysis method integer. Stopping analysis()");
        }
        
        // Insert method integer otherwise
        retmethods.insert(methods[i]);
    }

    return retmethods;
}

// Manually exported in NAMESPACE
//' Analysis methods for absolute thresholding of weighted graphs
//' 
//' Weighted graphs thresholding analysis. This function iteratively performs absolute thresholding the input graph at specified threshold values [\code{lower}, \code{upper}].
//' The function's methods parameter controls which graph analysis methods are performed at each threshold value.
//' Execution of \code{analysis()} ends after the '\code{upper}' threshold is reached, or if the graph becomes
//' too small to threshold further.
//'  
//' The results at each step are written output files. There is always at least one output file, but there can be 
//' up to two depending on the methods passed:
//' \itemize{
//'   \item \code{<outfile_prefix>.iterative.txt}:          (\strong{guaranteed})
//'   \item \code{<outfile_prefix>.statistical_errors.txt}: (method \strong{2})
//' }
//' 
//' Refer to Dr. Carissa Bleker's dissertation for more information about these analysis methods: \url{https://trace.tennessee.edu/utk_graddiss/5894/}
//' 
//' @param infile string. File path for .ncol graph file (\url{https://lgl.sourceforge.net/}) to read in for analysis. 
//' This file must be space delimited for this function to properly read in the graph's information.
//' @param outfile_prefix string. Prefix of output file in which analysis will be redirected to. If this is not specified,
//'        \code{thresholding::analysis()} will auto generate the output file prefix to include the input file's 
//'        prefix and the method numbers (in ascending order)
//'          For example if the user requests methods \code{4} and \code{7} with an input file named \code{"myfile.tsv"} would be:
//'             \code{myfile-47.<method_name>.txt}. 
//'        The input file prefix will be determined by the characters preceding the first period ('.') character.
//'         Otherwise, the user-provided value for \code{outfile_prefix} and the method numbers (if any) will be used in the output file's name.
//' @param methods Numeric vector of analysis method integers. Defaults to an empty list (no analysis methods). These will be performed on each thresholding step and recorded
//' into the output file(s). The number to method translation is given below:  
//' \itemize{
//'   \item 1 = all
//'   \item 2 = significance and power calculations (only valid for Pearson CC)
//'   \item 3 = scale free
//'   \item 4 = maximal cliques
//'   \item 5 = spectral methods
//'   \item 6 = random matrix theory
//'   \item 7 = clustering coefficient
//'   \item 8 = percolation
//' }
//'  Refer to the following papers for more detailed description of these methods: 
//' \itemize{
//' \item Dr. Carissa Bleker's thresholding dissertation: \url{https://trace.tennessee.edu/utk_graddiss/5894/} 
//' \item Dr. Langston, Grady, and Bleker's thresholding paper: \url{https://web.eecs.utk.edu/~mlangsto/JCB-Thresholding-Paper.pdf}
//' }
//'         The method name in the outputfile name may vary based on the methods used. The name will either be \code{iterative} or \code{statistical_errors}.
//' @param lower numeric.Lower bound to begin thresholding loop at (default = 0.5 ; lower >= 0)
//' @param upper numeric. Hard upper bound that ends thresholding  loop when \code{lower} value is greater than \code{upper} value (Default = 0.99)
//' @param increment numeric. Size of increment step in the thresholding loop
//' @param window_size numeric. Sliding window size for spectral method (Default = 5)
//' @param min_partition_size integer. minimum size of graph or subgraph after threshold (Default = 10)
//' @param min_clique_size integer. Minimum size of maximal cliques in maximal clique count (Default = 5)
//' @param num_samples integer. Number of samples in Pearson Correlation Coefficient data (only used for \strong{analysis method 2 - Power and Significance calculations}).  
//'        \emph{\strong{\code{num_samples} must be positive, non-zero, and match the number of samples from the original dataset for method 2 to work.}} \emph{Not used in any other methods.}
//' @param significance_alpha numeric. Probability of rejecting the null hypothesis when the null hypothesis is true.  (Default = 0.01)
//' @param overwrite boolean. Determines whether output file with given or generated prefix will be overwritten. 
//'        Set this to \code{TRUE} to force overwrite the output file. The default, \code{FALSE}, will display a menu asking
//'        whether or not you wish to overwrite the output file. The examples are defaulted to \code{TRUE} to avoid prompting this menu if the 
//'        example output files already exist.
//' @examples
//' #######    Variable Set-Up     #######
//' \dontrun{
//' library(thresholding)
//' data_file <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "thresholding")   
//' data_prefix <- './HumanCellCycleSubset-thresh-ex'  # prefix used for output file(s)
//' lower <- 0.6   
//'
//' #######    Example 1 - No methods #######
//' analysis(data_file, 
//'          outfile_prefix = data_prefix,
//'          lower = lower,
//'          overwrite = TRUE
//'          )
//' } 
//'
//' \dontrun{ 
//' #######    Example 2 - iterative and power/significance methods #######
//' data_file <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "thresholding") 
//' methods <- c(8, 2, 3)    # select the three desired analysis methods
//' lower <- 0.6             # choose lower bound thresholding value the thresholding loop begins at
//' num_samples <- 13        # ONLY FOR METHOD 2 - number of samples in data set
//' 
//' # Note: analysis() will autogenerate an output file name based 
//' #       on the input file path if a prefix is not passed.
//' analysis(data_file, 
//'          methods = methods,
//'          lower = lower,
//'          num_samples = num_samples,
//'          overwrite = TRUE
//'          )
//' }
//'
//' \dontrun{
//' #######    Example 3 - More Iterative methods #######
//' # WARNING FOR EXAMPLES: this code may take a few minutes to run.
//' data_file <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "thresholding") 
//' methods <- c(3,4,5,6)
//' lower <- 0.7
//' 
//' # Note: By supplying an output file prefix, the only modification to the output
//' #       file will be "-###" where ### represents the analysis methods passed
//' #       to analysis(). This is done to distinguish between output files generated by analysis()
//' #       that have the same prefix, but potentially different methods.
//' analysis(data_file,
//'          outfile_prefix = "./example_2_it",
//'          methods = methods,
//'          lower = lower,
//'          overwrite = TRUE
//'         )
//'
//' } 
//' # End "dontrun" block
//' @returns Nothing. \code{analysis()} writes all output to a file. The file path and file prefixes are printed on standard output when `analysis()` terminates.
// [[Rcpp::export]]
std::string analysis(std::string infile, 
              std::string outfile_prefix="",
              Rcpp::NumericVector methods=Rcpp::NumericVector::create(), 
              double lower=0.5,
              double upper=0.99,
              double increment=0.01,
              int window_size=5,
              int min_partition_size=10,
              int min_clique_size=5,
              int num_samples=0,
              double significance_alpha=0.01,
              bool overwrite = false
              )
{
    // Save a string to return in case of an early return
    const std::string fail_string = "analysis_failed";

    // Make sure the input file is able to be accessed
    std::ifstream exists;
    exists.open(infile);
    if(exists.good() != 1){
        Rcpp::Rcerr << "The input file " << infile << " had a problem during opening.\n";
        Rcpp::Rcerr << "Please make sure the input file is spelled correctly and\n";
        Rcpp::Rcerr << "you have the proper permissions to read it.\n";
        Rcpp::stop(infile + " was not able to be opened");
    }
    exists.close();
    
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
        Rcpp::Rcerr << "Use thresholding::help() to get more information on analysis() inputs." << '\n';
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
    analysis_methods = parse_methods_list(methods);


    // CARISSA'S METHOD:Get the pid to create unique file names if the same input file is run multiple times
    //std::string str_pid = get_str_pid();

    // STANDARD WAY: Create a string from the user's methods list (if non-empty)
    //               This way, if a user requests the same methods, then the file is just updated.
    //               Runs with different methods will still have different output file names.
    //               If no methods are requested, then the output prefix will not have the additional 
    //               hyphonated methods numbers in it.
    std::string str_methods;
    if(*analysis_methods.begin() == -1){
        str_methods = "";
    }
    else{
        std::set<int>::iterator sit;
        for(sit=analysis_methods.begin(); sit!=analysis_methods.end(); sit++){
            int current_method = *sit;
            str_methods += std::to_string(current_method);
        }
    }
    
    std::string stripped_prefix;
    
    if(outfile_prefix == ""){
        int findex;
        findex = infile.find_last_of(".");
        stripped_prefix = (findex != -1) ? infile.substr(0, findex) : infile;
        outfile_prefix = stripped_prefix + ((str_methods!="") ? "-" + str_methods : "");
    } 
    else{
      stripped_prefix = outfile_prefix;
      if(str_methods != "") outfile_prefix += "-" + str_methods;
    }
    
      
    // Present user with the option to overwrite the output based on the output file
    // and its prefix if the overwrite parameter is left as FALSE
    // Set the parameter to TRUE for unconditional overwriting
    Rcpp::Function r_glob("Sys.glob");   
    Rcpp::Function r_menu("menu");
    Rcpp::StringVector fileNames = r_glob(outfile_prefix + "." + "*");

    if(fileNames.length() > 0){
        for(auto file : fileNames){
            if(!overwrite && file != infile){
                Rcpp::Rcout << "\nYou are about to overwrite the output file:\n";
                Rcpp::Rcout << "    " << file << "\n";
                Rcpp::Rcout << "Continue with graph threshold analysis?\n";
                Rcpp::CharacterVector options = Rcpp::CharacterVector::create("Yes", "No");
                int response = Rcpp::as<int>(Rcpp::wrap(r_menu(options)));
                if(response == 2){
                    Rcpp::Rcout << "--analysis() will not overwrite " << file << ".\n";
                    Rcpp::stop("\rLeaving analysis() early.");
                }
            }        
            
        }
        Rcpp::Rcout << "\n";
        Rcpp::Rcout << "----Continuing with analysis().\n";
        Rcpp::Rcout << "----Overwriting output files with prefix of " << outfile_prefix << "\n\n";
    }

    // Stores the outfile name passed to analysis functions at multiple points throughout 
    // the analysis execution
    std::string outfile_name;

    /////////////////////////////////////////////////////////////////////////////////////////
    // 1 = Method for finding significance and power calculations (only valid for Pearson CC)
    // Type I error (false positive rate) and
    // Type II error (false negative rate) control
    // Have to have n - number of samples (not number of variables)
    if(analysis_methods.find(2) != analysis_methods.end()){
        outfile_name = outfile_prefix + ".statistical_errors.txt";
        // TODO: Fix setjmp() failure in dependency to restore old formal argument: 
        // bool bonferroni_corrected=0,  
        // This was was in OG packge, but  setjmp() to fail in lglib::invstudenttdistribution 
        // in src/specialfuntions.cpp
        control_statistical_errors(significance_alpha,
                                  num_samples,
                                  0, //E
                                  false,   // was bonferroni_corrected previously. Turned off so doesn't break when true
                                  outfile_name);
        analysis_methods.erase(2);
    }

    // End program if the only method desired was the significance and power calculations
    if (analysis_methods.size() == 0){
        return fail_string;
    }

    // Turn on attribute handling
    // For igraph to handle edge weights
    igraph_set_attribute_table(&igraph_cattribute_table);

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
    Rcpp::Rcout << "  (maximum possible number of edges " << int(orig_max_E) << ")\n";
    Rcpp::Rcout << "-----------------------------------------------------\n\n";

    ///////////////////////////////////////////////////////////////////////
    // Non-loop methods
    ///////////////////////////////////////////////////////////////////////

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
    header << "\tvertex.count\tedge.count";
    header << "\tconnected.component.count";
    header << "\tdensity\tdensity.orig.V";
    header << "\tlargest.cc.size\tsecond.largest.cc.size";
    header << "\tclustering.coefficient\trandom.clustering.coefficient";
    header << "\tsecond.eigenvalue\talmost.disconnected.component.count";
    header << "\tmaximal.clique.count\tclique.number";
    header << "\tpoisson.chi2\tpoisson.pvalue";
    header << "\tgoe.chi2\tgoe.pvalue";
    header << "\tscale.free.KS\tscale.free.KS.p.value\tscale.free.alpha";
    out << header.str();
    out << '\n';

    // Get the threshold increments - range() function from math_ext
    double t;
    const std::vector<double> t_vector = range(lower, upper, increment);
    int num_increments = t_vector.size();

    Rcpp::Rcout << "Iterative thresholding\n";
    Rcpp::Rcout << "Number steps: " << num_increments << '\n';
    Rcpp::Rcout << "Iterative analysis methods requested:" << '\n';

    // Remind user of iterative methods for analysis they have selected
    // and are performed on the graph at each thresholding step
    std::string mname;
    for(const int &m : analysis_methods){

        if(m == -1){
            Rcpp::Rcout << "\tNo additional iterative analysis methods requested\n";
            Rcpp::Rcout << "\tContinuing with default iterative analysis...\n";
            break;
        }
        
        switch(m){
            case 3:
                mname = "Scale Free";
                break;
            case 4:
                mname = "Maximal Clique";
                break;
            case 5:
                mname = "Spectral Thresholding";
                break;
            case 6:
                mname = "Random Matrix Theory";
                break;
            case 7:
                mname = "Clustering Coefficient";
                break;
            case 8:
                mname = "Percolation";
                break;
            default:
                Rcpp::Rcout << "\t" << m << ": Method not recognized. Continuing with analysis...\n";
                continue;
        }

        Rcpp::Rcout << "\t" << m << ": " << mname << " Methods\n";
    }

    Rcpp::Rcout << "-----------------------------------------------------\n\n"; // formatting

    // Initialise necessary stuff for iterative analysis methods to use
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
    // Thresholds graph and runs requested additional analysis methods
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
                // (No additional analysis methods requested)
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
                    // int threshold_status2 = threshold_graph(t, G_random);
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
    Rcpp::Rcout << "Thresholding and analysis completed!!!\n";
    Rcpp::Rcout << "\nUse the thresholding::get_results() function along with the\n";
    Rcpp::Rcout << "output file PREFIX (which excludes the analysis method integers\n";
    Rcpp::Rcout << "and extension from the file path) to retrieve the recommended\n";
    Rcpp::Rcout << "threshold for each requested method.\n";
    Rcpp::Rcout << "\nExample: \n";
    Rcpp::Rcout << "\tFile name/path: " << outfile_name << '\n';
    Rcpp::Rcout << "\tPrefix:         " << stripped_prefix << '\n';

    igraph_destroy(&G);

    return stripped_prefix;
}

 