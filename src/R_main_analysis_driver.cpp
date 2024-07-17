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

/////////////// update later after basic functionality of analysis program works //////////////
// Helper code to print out help information if inputs are incorrect

// Display argument information on terminal for thresholding::analysis()
// 
// NOT EXPORTED YET AS OF VERSION 1.0.0
// void help(){
//     /*
//     Rcpp::Rcerr <<  "\n";
//     Rcpp::Rcerr <<  "    Usage: \n";
//     Rcpp::Rcerr <<  "    " << "thresholdAnaylsis"  << " [-OPTIONS]... <GRAPH FILE PATH> <OUTPUT FILE PATH> \n\n";
//     Rcpp::Rcerr <<  "    Graph has to be in .ncol format. \n";
//     Rcpp::Rcerr <<  "    Output file path is the prefix to the results files, which will be of the form: \n";
//     Rcpp::Rcerr <<  "        <OUTPUT FILE PATH>.pid.<method_name>.txt\n\n";
//     Rcpp::Rcerr <<  "    Options: \n";
//     Rcpp::Rcerr <<  "      -l  --lower                  <value>     lower bound on thresholds to test (default 0.5)\n";
//     Rcpp::Rcerr <<  "      -u  --upper                  <value>     upper bound on thresholds to test (default 0.99)\n";
//     Rcpp::Rcerr <<  "      -i  --increment              <value>     threshold increment (default 0.01)\n";
//     Rcpp::Rcerr <<  "      -w  --windowsize             <value>     sliding window size for spectral method (default 5)\n";
//     Rcpp::Rcerr <<  "      -p  --minimumpartitionsize   <value>     minimum size of graph or subgraph after threshold (default 10)\n";
//     Rcpp::Rcerr <<  "      -n  --num_samples            <value>     number of samples for significance and power calculations (default NULL)\n";
//     Rcpp::Rcerr <<  "      -b  --bonferroni_correction              switch to perform bonferroni corrections in significance and power calculations (default FALSE)\n";
//     Rcpp::Rcerr <<  "      -c  --minimum_cliquesize     <value>     minimum size of maximal cliques in maximal clique count (default 5)\n";
//     Rcpp::Rcerr <<  "      -m  --methods                <value>     comma separated list of methods (defaults to none)\n";
//     Rcpp::Rcerr <<  "                                                   0 - all\n";
//     Rcpp::Rcerr <<  "                                                   1 - significance and power calculations (only valid for Pearson CC)\n";
//     Rcpp::Rcerr <<  "                                                   2 - local-global\n";
//     Rcpp::Rcerr <<  "                                                   3 - scale free\n";
//     Rcpp::Rcerr <<  "                                                   4 - maximal cliques\n";
//     Rcpp::Rcerr <<  "                                                   5 - spectral methods\n";
//     Rcpp::Rcerr <<  "                                                   6 - random matrix theory\n";
//     Rcpp::Rcerr <<  "                                                   7 - clustering coefficient\n";
//     Rcpp::Rcerr <<  "                                                   8 - percolation\n";
//     Rcpp::Rcerr <<  "      -h  --help                               print this help and exit\n";
//     Rcpp::Rcerr <<  "\n";
//     */
//     Rcpp::Rcout <<  "\n";
//     Rcpp::Rcout <<  "--------------Langston Lab Thresholding Analysis Techniques (2023)---------------\n";
//     Rcpp::Rcout <<  "                               analysis()                               \n";
//     Rcpp::Rcout <<  "Synopsis:\n";
//     Rcpp::Rcout <<  "    analysis(infile, outfile_prefix,\n";
//     Rcpp::Rcout <<  "                      methods="",  lower=0.5,  upper=0.99, \n";
//     Rcpp::Rcout <<  "                      increment=0.01,  window_size=5,  min_partition_size=10, \n";
//     Rcpp::Rcout <<  "                      min_clique_size=5,  min_alpha=0,  max_alpha=4,\n";
//     Rcpp::Rcout <<  "                      alpha_increment=0.1,  num_samples=0,\n";
//     Rcpp::Rcout <<  "                      significance_alpha=0.01,  bonferroni_corrected=0)\n";
//     Rcpp::Rcout <<  "\n";
//     Rcpp::Rcout <<  "Arguments to analysis():\n";
//     Rcpp::Rcout <<  "Required:\n";
//     Rcpp::Rcout <<  "\t1.) infile: string input\n";
//     Rcpp::Rcout <<  "\t        The weighted edge list (.wel) file input. This file is in the .ncol format as specified by\n";
//     Rcpp::Rcout <<  "\t        the Large Graph Layout group: https://lgl.sourceforge.net/#FileFormat.\n\n";
    
//     Rcpp::Rcout <<  "\t        In this application, the graph input file is simple, weighted, undirected. The vertices in the .wel\n";
//     Rcpp::Rcout <<  "\t        file follow the following format where the arrow represents whitespace: \n";
//     Rcpp::Rcout <<  "\t          vertex1⇥vertex2⇥weight1,2\n";
//     Rcpp::Rcout <<  "\t          vertex1⇥vertex3⇥weight1,3\n";
//     Rcpp::Rcout <<  "\t          ...\n\n";

//     Rcpp::Rcout <<  "\t         NOTE: vertex names cannot contain whitespace.\n\n";

//     Rcpp::Rcout <<  "\t2.) outfile_prefix: string input\n";
//     Rcpp::Rcout <<  "\t        Prefix to the output files to the analysis file(s).\n";
//     Rcpp::Rcout <<  "\t        Example: If the prefix is \"graph-output\", the output is graph-output.<method_name>.txt";
//     Rcpp::Rcout <<  "\t                 where method name is the type of analysis performed, such as iterative or statistical_errors\n";
//     Rcpp::Rcout <<  "\t                 The method(s) are controlled by the optional \"methods\" argument.\n\n";

//     Rcpp::Rcout <<  "Optional:\n";
//     Rcpp::Rcout <<  "\t3.) methods: string input  (defaults to empty string)\n";
//     Rcpp::Rcout <<  "\t        Comma separated list of analysis operations to complete. These methods are represented by an integer";
//     Rcpp::Rcout <<  "\t        which is mapped to its corresponding method internally. The following methods are currently available:\n\n";
//     Rcpp::Rcout <<  "\t             0 - all (methods 1-7 will be performed)\n";
//     Rcpp::Rcout <<  "\t             1 - significance and power calculations (only valid for Pearson CC)\n";
//     Rcpp::Rcout <<  "\t             2 - local-global\n";
//     Rcpp::Rcout <<  "\t             3 - scale free\n";
//     Rcpp::Rcout <<  "\t             4 - maximal cliques\n";
//     Rcpp::Rcout <<  "\t             5 - spectral methods\n";
//     Rcpp::Rcout <<  "\t             6 - random matrix theory\n";
//     Rcpp::Rcout <<  "\t             7 - clustering coefficient\n";
//     Rcpp::Rcout <<  "\t             8 - percolation\n\n";
//     Rcpp::Rcout <<  "\t        Example: methods=\"2, 5\""; 
//     Rcpp::Rcout <<  "\t                 methods=\"6, 1\"  (Note: methods will be performed in numerical order internally, but the order";
//     Rcpp::Rcout <<  "\t                 which they are passed to the function doesn't matter.)\n\n";

//     Rcpp::Rcout <<  "\t4.) lower: floating point input  (defaults to 0.5)\n";
//     Rcpp::Rcout <<  "\t        Initial lower bound for  thresholding loop. The loop ends when the current threshold value surpasses the upper bound limit.\n";
//     Rcpp::Rcout <<  "\t        NOTe: lower must be less than or equal to upper (lower <= upper) for function to continue.\n\n";
//     Rcpp::Rcout <<  "\t        Example: lower=0.6\n\n";

//     Rcpp::Rcout <<  "\t5.) upper: floating point input (defaults to 0.99)\n";
//     Rcpp::Rcout <<  "\t        Upper bound for the thresholding loop; Thresholding ends when the current threshold value is greater than";
//     Rcpp::Rcout <<  "\t        this parameter.\n";
//     Rcpp::Rcout <<  "\t        NOTE: the upper must be greater than or equal to the lower input (lower <= upper) for the function to continue.\n\n";
//     Rcpp::Rcout <<  "\t        Example: upper=0.93\n\n";

//     Rcpp::Rcout <<  "\t6.) increment: floating point input (defaults to 0.01)\n";
//     Rcpp::Rcout <<  "\t        This value controls the step of the thresholding loop. On each pass, the graph is thresholded at the current";
//     Rcpp::Rcout <<  "\t        thresholding value.\n";
//     Rcpp::Rcout <<  "\t        After the thresholding step, the value is incremented by the increment parameter (which is 0.01 by default).\n";
//     Rcpp::Rcout <<  "\t        The increment parameter gives finer control to which thresholding values are used in the analysis. This parameter can also\n";
//     Rcpp::Rcout <<  "\t        work alongside the lower and upper parameters to limit the scope and depth of thresholding userd.\n\n";
//     Rcpp::Rcout <<  "\t        Example: increment=0.0001  - finer grain thresholding\n";
//     Rcpp::Rcout <<  "\t                increment=0.05    - coarser thresholding\n\n";

//     Rcpp::Rcout <<  "\t7.) window_size: integer input (defaults to 5)\n";
//     Rcpp::Rcout <<  "\t         NOTE: this parameter is only used for spectral graph methods, which are used in the local-global (#2) and spectral (#5) analysis methods\n";
//     Rcpp::Rcout <<  "\t         Used in spectral methods to create a differences vector with the specified sliding window width.\n";
//     Rcpp::Rcout <<  "\t         window_size controls the size of the difference vector and the distance of the window between each difference pair.\n";
//     Rcpp::Rcout <<  "\t         For example, a vector of size 10 with elements [0,1,2,3,4,5,6,7,8,9] exists.\n";
//     Rcpp::Rcout <<  "\t         IF window_size=7, the output vector will have a size of 3. The first elements compared are 0 [ind = 0] and 7 [ind = 7].\n";
//     Rcpp::Rcout <<  "\t         This difference is then stored. Next, the window is shifted by one. The next elements compared are\n";
//     Rcpp::Rcout <<  "\t         1 [ind = 1] and 8 [ind = 8]. This process repeats until the window extends past the end of the vector.\n";
//     Rcpp::Rcout <<  "\t         The vector that this internal difference function makes is then [7, 7, 7].\n\n";
//     Rcpp::Rcout <<  "\t         Example: window_size=10   (NOTE: window_size should be greater than the min_partition_size. If this is not true, then\n";
//     Rcpp::Rcout <<  "\t         the default values for each parameter will be used.)\n\n";

//     Rcpp::Rcout <<  "\t8.) min_partition_size: integer input (defaults to 10)\n";
//     Rcpp::Rcout <<  "\t";        
// }

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
        // User wants 
        if(methods[i] == 0){
            retmethods.insert({1, 2, 3, 4, 5, 6, 7, 8});
            break;
        }
        // Currently only integers 0-8 are accepted as input for additional analysis methods
        // File output naming scheme attaches a concatenated string of all analysis methods
        // to the prefix (either default or specified by user).
        // Previously, adding a method int that wasn't valid would just add to the string
        // and nothing would happen in the main thresholding loop. 
        //
        // This means adding two invalid methods 12 and 34 would lead to the same file name
        // as specifying methods 1,2,3, and 4. However, it would be misleading as the analysis
        // call with 12 and 34 wouldn't have the desired outputs for methods 1-4. 
        // Therefore, this validates input before entering the main analysis loop.
        // 
        // If the user doesn't want any additinal analysis methods ran, then the methods list remains
        // empty and normal iterative thresholding is performed.
        else if(methods[i] < 0 || methods[i] > 8){
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
//' up to three depending on the methods passed:
//' \itemize{
//'   \item \code{<outfile_prefix>.iterative.txt}:          (\strong{guaranteed})
//'   \item \code{<outfile_prefix>.statistical_errors.txt}: (method \strong{1})
//'   \item \code{<outfile_prefix>.local_global.txt}:        (method \strong{2})
//' }
//' 
//' Refer to Dr. Carissa Bleker's dissertation for more information about these analysis methods: \link{https://trace.tennessee.edu/utk_graddiss/5894/}
//' 
//' @param infile File path for .ncol graph file (\link{https://lgl.sourceforge.net/}) to read in for analysis. 
//' This file must be space delimited for this function to properly read in the graph's information.
//' @param methods Numeric vector of method integers. Defaults to an empty list. The number to method translation is given below:
//' \itemize{
//'   \item 0 = all
//'   \item 1 = significance and power calculations (only valid for Pearson CC)
//'   \item 2 = local-global
//'   \item 3 = scale free
//'   \item 4 = maximal cliques
//'   \item 5 = spectral methods
//'   \item 6 = random matrix theory
//'   \item 7 = clustering coefficient
//'   \item 8 = percolation
//' }
//'  Refer to the following papers for more detailed description of these methods: 
//' \itemize{
//' \item Carissa Bleker's thresholding dissertation: \link{https://trace.tennessee.edu/utk_graddiss/5894/} 
//' \item Dr. Langston, Grady, and Bleker's thresholding paper: \link{https://web.eecs.utk.edu/~mlangsto/JCB-Thresholding-Paper.pdf}
//' }
//' @param outfile_prefix Prefix of output file in which analysis will be redirected to. If this is not specified,
//'        \code{thresholding::analysis()} will auto generate the output file prefix to include the input file's 
//'        prefix and the ascending method numbers. 
//'        The input file prefix will be determined by the characters preceding the first period ('.') character.
//'          An example if the user requests methods \code{4} and \code{7} with an input file named \code{"myfile.tsv"} would be:
//'             \code{myfile-47.<method_name>.txt}
//'         Method name can vary based on the methods used. This will either be \code{iterative}, \code{local_global}, or \code{statistical_errors}.
//' @param lower Lower bound to begin thresholding loop at (default = 0.5 ; lower >= 0)
//' @param upper Hard upper bound that ends thresholding  loop when \code{lower} value is greater than \code{upper} value (Default = 0.99)
//' @param increment Size of increment step in the thresholding loop
//' @param window_size Sliding window size for spectral method (Default = 5)
//' @param min_partition_size minimum size of graph or subgraph after threshold (Default = 10)
//' @param min_clique_size Minimum size of maximal cliques in maximal clique count (Default = 5)
//' @param min_alpha Starting alpha value used in \strong{method 2 - local global pruning}.  \emph{Not used in any other methods.}  (Default = 0.0) 
//' @param max_alpha Ending alpha value used in \strong{method 2 - local global pruning}.  \emph{Not used in any other methods.}   (Default = 4.0)
//' @param alpha_increment Size of increment of alpha value in \strong{method 2 - local global pruning}'s main loop. \emph{Not used in any other methods.}  (Default = 0.1)
//' @param num_samples Number of samples in Pearson Correlation Coefficient data (only used for \strong{analysis method 1 - Power and Significance calculations}).  
//'        \emph{\strong{\code{num_samples} must be positive, non-zero, and match the number of samples from the original dataset for method 1 to work.}} \emph{Not used in any other methods.}
//' @param significance_alpha Probability of rejecting the null hypothesis when the null hypothesis is true.  (Default = 0.01)
//' @param bonferroni_corrected Option to perform Bonferroni correction in \strong{method 1 - significance and power calculations}. 
//'        Applies Bonferroni correction to the value of significance_alpha if set to \code{TRUE}. \emph{Not used in any other methods.} (default \code{FALSE})
//' @param overwrite Determines whether output file with given or generated prefix will be overwritten. 
//'        Set this to \code{TRUE} to force overwrite the output file. The default, \code{FALSE}, will display a menu asking
//'        whether or not you wish to overwrite the output file. 
//' @examples
//' #######    Variable Set-Up     #######
//' data_file <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "thresholding")   # .ncol weighted edge list
//' data_prefix <- './example/HumanCellCycleSubset-thresh'  # prefix used for output file(s)
//' lower <- 0.6   
//'
//' #######    Example 1 - No methods #######
//' analysis(data_file, 
//'          data_prefix,
//'          lower = lower,
//'          )
//'
//' #######    Example 2 - Iterative methods #######
//' # WARNING FOR EXAMPLES: this code may take a few minutes to run.
//' data_file <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "thresholding") 
//' methods <- c(3,4,5,6)
//' lower <- 0.7
//' 
//' # Note: By supplying an output file prefix, the only modification to the output
//'         file will be "-###" where ### represents the analysis methods passed
//'         to analysis(). This is done to distinguish between output files generated by analysis()
//'         that have the same prefix, but potentially different methods.
//' analysis(data_file,
//'          outfile_prefix = "./example_2_it",
//'          methods = methods,
//'          lower = lower,
//'         )
//'
//'
//' #######    Example 3 - iterative and power/significance methods #######
//' data_file <- system.file('extdata', 'HumanCellCycleSubset.ncol', package = "thresholding") 
//' methods <- c(8, 1, 3)    # select the three desired analysis methods
//' lower <- 0.6             # choose lower bound thresholding value the thresholding loop begins at
//' num_samples <- 13        # ONLY FOR METHOD 1 - number of samples in data set
//' 
//' # Note: analysis() will autogenerate an output file name based 
//' #       on the input file path if a prefix is not passed.
//' analysis(data_file, 
//'          methods = methods,
//'          lower = lower,
//'          num_samples = num_samples,
//'          )
//' @returns Nothing. \code{analysis()} writes all output to a file. The file path and file prefixes are printed on standard output when `analysis()` terminates.
// [[Rcpp::export]]
void analysis(std::string infile, 
              Rcpp::NumericVector methods=Rcpp::NumericVector::create(), 
              std::string outfile_prefix="",
              double lower=0.5,
              double upper=0.99,
              double increment=0.01,
              int window_size=5,
              int min_partition_size=10,
              int min_clique_size=5,
              double min_alpha=0.0,
              double max_alpha=4,
              double alpha_increment=0.1,
              int num_samples=0,
              double significance_alpha=0.01,
              bool bonferroni_corrected=0,
              bool overwrite = false
              )
{

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
      outfile_prefix += "-" + str_methods;
    }
    
      
    // Present user with the option to overwrite the output based on the output file
    // and its prefix
    // if the overwrite parameter is left as FALSE
    // Set the parameter to TRUE for unconditional overwriting
    Rcpp::Function r_glob("Sys.glob");   
    Rcpp::Function r_menu("menu");
    Rcpp::StringVector fileNames = r_glob(outfile_prefix + "." + "*");

    if(fileNames.length() > 0){
        for(auto file : fileNames){
            if(!overwrite){
                Rcpp::Rcout << "\nYou are about to overwrite the output file:\n";
                Rcpp::Rcout << "    " << file << "\n";
                Rcpp::Rcout << "Continue with graph threshold analysis?\n";
                Rcpp::CharacterVector options = Rcpp::CharacterVector::create("Yes", "No");
                int response = Rcpp::as<int>(Rcpp::wrap(r_menu(options)));
                if(response == 2){
                    Rcpp::Rcout << "--analysis() will not overwrite " << file << ".\n";
                    Rcpp::stop("\rLeaving analysis()f early.");
                }
            }        
            
        }
        Rcpp::Rcout << "\n----Continuing with analysis().\n";
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
    Rcpp::Rcout << "  (maximum possible number of edges " << int(orig_max_E) << ")\n";
    Rcpp::Rcout << "-----------------------------------------------------\n\n";

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
                Rcpp::Rcout << "\t" << m << ": Method not recognized.\n";
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
    Rcpp::Rcout << "Thresholding and analysis completed!!!\n";
    Rcpp::Rcout << "\nUse the thresholding::get_results() function along with the\n";
    Rcpp::Rcout << "output file PREFIX (which excludes the analysis method integers\n";
    Rcpp::Rcout << "and extension from the file path) to retrieve the recommended\n";
    Rcpp::Rcout << "threshold for each requested method.\n";
    Rcpp::Rcout << "\nExample: \n";
    Rcpp::Rcout << "\tFile name/path:" << outfile_name << '\n';
    Rcpp::Rcout << "\tPrefix:        " << stripped_prefix << '\n';

    igraph_destroy(&G);
}


/*

Manual usage statement for analysis().
//' @usage analysis(infile, outfile_prefix, [methods], [lower], [upper], \
//'                 [increment], [window_size], [min_partition_size], \
//'                 [min_clique_size], [min_alpha], [max_alpha], \
//'                 [alpha_increment], [num_samples], [significance_alpha], \
//'                 [bonferroni_corrected])

*/