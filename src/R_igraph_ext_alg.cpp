#include <Rcpp.h>
#include "igraph_ext.h"
// Threshold graph
// by removing edges with abs weight less than "t"
// and subsequently vertices with no neighbours
int threshold_graph(double t, igraph_t &G){
    // 1 = all edges removed
    // 2 = some edges are removed
    // 3 = no edges removed

    // identify edges to remove
    igraph_vector_t edge_indices;
    igraph_vector_init(&edge_indices, 0);

    igraph_real_t w;
    igraph_integer_t E = igraph_ecount(&G);
  
    for (int i=0; i<E; i++){
        w = igraph_cattribute_EAN(&G, "weight", i);
        if(w < t && w > -t){
            // add this edge index to list of edges to be deleted
            igraph_vector_push_back(&edge_indices, i);
        }
    }

    if(igraph_ecount(&G) <= igraph_vector_size(&edge_indices)){
        // removing all edges, so could just delete all vertices
        //igraph_delete_vertices(&G, igraph_vss_all());

        // clean up
        igraph_vector_destroy(&edge_indices);
        return 1;
    }
    else if(igraph_vector_size(&edge_indices) > 0){

        // remove edges
        igraph_delete_edges(&G, igraph_ess_vector(&edge_indices));
        //std::cout << " Removed " << igraph_vector_size(&edge_indices) << " edges and ";

        // clean up
        igraph_vector_destroy(&edge_indices);

        // identify and remove degree 0 vertices
        igraph_vector_t vertex_degrees;
        igraph_vector_init(&vertex_degrees, 0);

        igraph_vector_t vertex_indices;
        igraph_vector_init(&vertex_indices, 0);

        igraph_degree(&G, &vertex_degrees, igraph_vss_all(), IGRAPH_ALL, false);

        for(long int i=0; i<igraph_vcount(&G); i++){
            if(VECTOR_IG(vertex_degrees)[i] == 0){
                igraph_vector_push_back(&vertex_indices, i);
            }
        }

        if(igraph_vector_size(&vertex_indices) > 0){
            // remove them
            igraph_delete_vertices(&G, igraph_vss_vector(&vertex_indices));
        }
        //std::cout << igraph_vector_size(&vertex_indices) << " vertices. " << std::flush;

        // clean up
        igraph_vector_destroy(&vertex_degrees);
        igraph_vector_destroy(&vertex_indices);


        return 2;
    }
    // last option is no edges to remove, don't change G, so do nothing
    else if(igraph_vector_size(&edge_indices) == 0){
        //std::cout << " Removed 0 edges and 0 vertices. ";

        // clean up
        igraph_vector_destroy(&edge_indices);

        return 3;
    }

    else{
        return -1;
    }
}

// Identify largest connected component of the graph and induce
int largest_connected_component(igraph_t &G, igraph_t &G_cc,
                                igraph_integer_t &cc_count,
                                igraph_integer_t &V_cc,
                                igraph_integer_t &V2_cc){
    // See also igraph_decompose, but since we only need
    // the largest CC, there is no point inducing all CCs.
    igraph_vector_t membership;
    igraph_vector_init(&membership, 0);
    igraph_vector_t csize;
    igraph_vector_init(&csize, 0);

    igraph_clusters(&G, &membership, &csize, &cc_count, IGRAPH_STRONG);

    // iterate over csize to find largest CC
    V_cc = 0;
    int max_cc_index;
    int this_cc_size;
    for(int i =0; i<cc_count; i++){
        this_cc_size = VECTOR_IG(csize)[i];
        if(this_cc_size > V_cc){
            max_cc_index = i;
            V_cc = this_cc_size;
        }
     }

    if(V_cc < 2){
        return 0;
    }

    // check if there is more than on CC with max_cc_size
    // TODO what to do if there is more that 1?
    //int num_max_cc = 0;
    //for(int i=0; i<cc_count; i++){
    //    this_cc_size = VECTOR_IG(csize)[i];
    //    if(this_cc_size == max_cc_size){
    //        num_max_cc++;
    //    }
    //}

    // get the size of the second largest CC
    V2_cc = 0;
    for(int i =0; i<cc_count; i++){
        this_cc_size = VECTOR_IG(csize)[i];
        if( (this_cc_size > V2_cc) && (this_cc_size < V_cc) ){
            V2_cc = this_cc_size;
        }
    }

    // identified largest CC, now collect its vertices
    igraph_vector_t vertices_in_cc;
    igraph_vector_init(&vertices_in_cc, V_cc);

    int j = 0; // index in vertices_in_cc
    for(int i=0; j<V_cc; i++){
        if(VECTOR_IG(membership)[i] == max_cc_index){
            VECTOR_IG(vertices_in_cc)[j] = i;
            j += 1;
        }
    }

    // induce the subgraph of the largest CC
    igraph_induced_subgraph(&G, &G_cc, igraph_vss_vector(&vertices_in_cc), IGRAPH_SUBGRAPH_AUTO);

    igraph_vector_destroy(&membership);
    igraph_vector_destroy(&csize);
    igraph_vector_destroy(&vertices_in_cc);

    return 0;
}

// Fiedler vector: eigen-vector corresponding to first non-zero eigen-value
// Assume connected graph -> 2nd eigenvector
int Fiedler_vector(igraph_t &G,
                   igraph_vector_t &eigenvector,
                   igraph_real_t &eigenvalue){

    //std::cout << " Attempting to get the Fiedler vector and value. " << std::flush;
    // make sure G has edges and that G is connected
    if(igraph_ecount(&G) < 1){
        Rcpp::stop("Fielder failed - no edges in graph.");
    }

    igraph_bool_t is_connected;
    igraph_is_connected(&G, &is_connected, IGRAPH_STRONG);
    if(!is_connected){
        Rcpp::stop(" Fielder failed - graph not connected.");
    }

    // dimension of laplacian = num vertices
    igraph_integer_t V = igraph_vcount(&G);

    igraph_matrix_t laplacian;
    igraph_matrix_init(&laplacian, V, V);

    // init eigen values and vectors
    igraph_vector_t values;
    igraph_matrix_t vectors;

    igraph_vector_init(&values, V); // first two eigenvalues will go in here.
    igraph_matrix_init(&vectors, V, 2); // number vertices by 2 eigenvectors

    // Weighted Laplacian of G
    //igraph_vector_t weights;
    //igraph_vector_init(&weights, 0);
    //
    //igraph_cattribute_EANV(&G, "weight", igraph_ess_all(IGRAPH_EDGEORDER_ID), &weights);
    //igraph_laplacian(&G, &laplacian, NULL, false, &weights);

    // Unweighted Laplacian of G
    igraph_laplacian(&G, &laplacian, NULL, false, NULL);


    // Eigen decomposition for symmetric matrices using LAPACK
    igraph_lapack_dsyevr(&laplacian, IGRAPH_LAPACK_DSYEV_SELECT, 0, 0, 0, 1, 2, 1e-8, &values, &vectors, 0);

    // should be 0.0 (or there abouts)
    //std::cout << " 1st eigenvalue: " << VECTOR_IG(values)[0] << std::endl;
    //std::cout << " 2nd eigenvalue: " << VECTOR_IG(values)[1] << std::endl;

    // set eigenvalue and eigenvector of interest
    eigenvalue = VECTOR_IG(values)[1];
    igraph_matrix_get_col(&vectors, &eigenvector, 1);

    // remove laplacian, vectors and values
    igraph_matrix_destroy(&laplacian);
    igraph_vector_destroy(&values);
    //igraph_vector_destroy(&weights);
    igraph_matrix_destroy(&vectors);

    return 0;
}

// See igraph_get_adjacency for unweighted graph
int get_weighted_adjacency(igraph_t &G, igraph_matrix_t &Adj){
    igraph_eit_t edgeit;                             // interate over edges
    long int V = igraph_vcount(&G);               // number vertices / matrix size

    long int from, to;
    igraph_integer_t ffrom, fto;

    igraph_matrix_resize(&Adj, V, V);
    igraph_matrix_null(&Adj); // set all entries to zero

    // create the edge iterator
    igraph_eit_create(&G, igraph_ess_all(IGRAPH_EDGEORDER_ID), &edgeit);

    // edge weights
    igraph_vector_t weights;
    igraph_vector_init(&weights, 0);
    igraph_cattribute_EANV(&G, "weight", igraph_ess_all(IGRAPH_EDGEORDER_ID), &weights);
    igraph_real_t w;


    while (!IGRAPH_EIT_END(edgeit)){
        long int edge=IGRAPH_EIT_GET(edgeit);
        igraph_edge(&G, (igraph_integer_t) edge, &ffrom, &fto);
        from=ffrom;
        to=fto;
        w = VECTOR_IG(weights)[edge];

        MATRIX_IG(Adj, from, to) = w;

        if (from != to){
            MATRIX_IG(Adj, to, from) = w;
        }

        IGRAPH_EIT_NEXT(edgeit);
     }

    // clean
    igraph_eit_destroy(&edgeit);


    return 0;
}
