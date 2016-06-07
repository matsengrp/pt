#include <cstdarg>
#include <search.h>
#include "partition.hpp"

/// @file partition.cpp
/// @brief Implementation of the Partition class.
/// Much of this was directly copied from libpll examples, so (c) Thomas Flouri.

namespace pt {

#define STATES    4
#define RATE_CATS 4


// if you remove this, you can remove the include for cstdarg above.
static void fatal(const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}


unsigned int ParseFasta(std::string path, unsigned int seq_count, char *** headers_out, char *** seqdata_out) {
  /* open FASTA file */
  pll_fasta_t * fp = pll_fasta_open(&path[0], pll_map_fasta);
  if (!fp) {
    std::cout << "Error opening file " << path << std::endl;
    exit(EXIT_FAILURE);
  }

  char * seq = NULL;
  char * hdr = NULL;
  long seqlen;
  long hdrlen;
  long seqno;

  /* allocate arrays to store FASTA headers and sequences */
  char ** headers = (char **)calloc(seq_count, sizeof(char *));
  char ** seqdata = (char **)calloc(seq_count, sizeof(char *));

  /* read FASTA sequences and make sure they are all of the same length */
  unsigned int i;
  int sites = -1;
  for (i = 0; pll_fasta_getnext(fp,&hdr,&hdrlen,&seq,&seqlen,&seqno); ++i)
  {
    if (i >= seq_count){
      fatal("FASTA file contains more sequences than expected.\n");
      exit(EXIT_FAILURE);
    }

    if (sites != -1 && sites != seqlen)
      fatal("FASTA file does not contain equal size sequences\n");

    if (sites == -1) sites = seqlen;

    headers[i] = hdr;
    seqdata[i] = seq;
  }

  /* did we stop reading the file because we reached EOF? */
  if (pll_errno != PLL_ERROR_FILE_EOF)
    fatal("Error while reading file %s", &path[0]);

  /* close FASTA file */
  pll_fasta_close(fp);

  if (sites < 0)
    fatal("Unable to read alignment");

  if (i != seq_count)
    fatal("Some taxa are missing from FASTA file");

  *headers_out = headers;
  *seqdata_out = seqdata;
  return (unsigned int) sites;
}


void EquipPartitionWithData(
  pll_partition_t * partition,
  pll_utree_t * tree,
  unsigned int tip_nodes_count,
  char ** headers,
  char ** seqdata) {

/*  obtain an array of pointers to tip nodes */
  pll_utree_t ** tipnodes = (pll_utree_t  **)calloc(tip_nodes_count,
                                                    sizeof(pll_utree_t *));
  pll_utree_query_tipnodes(tree, tipnodes);

  /* create a libc hash table of size tip_nodes_count */
  hcreate(tip_nodes_count);

  /* populate a libc hash table with tree tip labels */

  unsigned int * data = (unsigned int *)malloc(tip_nodes_count *
                                               sizeof(unsigned int));
  unsigned int i;
  for (i = 0; i < tip_nodes_count; ++i)
  {
    data[i] = i;
    ENTRY entry;
    entry.key = tipnodes[i]->label;
    entry.data = (void *)(data+i);
    hsearch(entry, ENTER);
  }

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < tip_nodes_count; ++i)
  {
    ENTRY query;
    query.key = headers[i];
    ENTRY * found = NULL;

    found = hsearch(query,FIND);

    if (!found)
      fatal("Sequence with header %s does not appear in the tree", headers[i]);

    unsigned int tip_clv_index = *((unsigned int *)(found->data));

    pll_set_tip_states(partition, tip_clv_index, pll_map_nt, seqdata[i]);
  }

  /* destroy hash table */
  hdestroy();

  /* we no longer need these two arrays (keys and values of hash table... */
  free(data);
  free(tipnodes);
}


void SetModelParameters(pll_partition_t * partition) {
  /* initialize the array of base frequencies */
  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };

  /* substitution rates for the 4x4 GTR model. This means we need exactly
     (4*4-4)/2 = 6 values, i.e. the number of elements above the diagonal */
  double subst_params[6] = {1,1,1,1,1,1};

  /* we'll use 4 rate categories, and currently initialize them to 0 */
  double rate_cats[4] = {0};

  /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
  pll_compute_gamma_cats(1, 4, rate_cats);

  /* set frequencies at model with index 0 (we currently have only one model) */
  pll_set_frequencies(partition, 0, frequencies);

  /* set 6 substitution parameters at model with index 0 */
  pll_set_subst_params(partition, 0, subst_params);

  /* set rate categories */
  pll_set_category_rates(partition, rate_cats);
}


/* a callback function for performing a full traversal */
static int cb_full_traversal(pll_utree_t * node)
{
  return 1;
}




/// @brief Constructor for Partition from a tree and an alignment.
/// @param[in] newick_path
/// Path to a file with a Newick-format tree string.
/// @param[in] fasta_path
/// Path to a FASTA-formatted alignment.
Partition::Partition(
  std::string newick_path,
  std::string fasta_path) {

  // Parse the unrooted binary tree in newick format, and store the number of
  // tip nodes in tip_nodes_count.
  tree_ = pll_utree_parse_newick(&newick_path[0], &tip_nodes_count_);
  if(!tree_) {
    std::cout << "Parsing failure: tree must be an unrooted binary tree.\n";
    exit(EXIT_FAILURE);
  };

  /*
  if(!TreeHealthy(tree_)) {
    std::cout << "Missing branch lengths in tree.\n";
    exit(EXIT_FAILURE);
  }
  */

  char ** headers = NULL;
  char ** seqdata = NULL;
  unsigned int sites_count =
    ParseFasta(fasta_path, tip_nodes_count(), &headers, &seqdata);

  partition_ = pll_partition_create(
    tip_nodes_count(),      // Number of tip sequences we want to have.
    inner_nodes_count(),    // Number of CLV buffers to be allocated for inner nodes.
    STATES,                 // Number of states that our data have.
    sites_count,            // Number of sites in the alignment.
    1,                      // Number of different substitution models (or eigen decomposition)
                            //     to use concurrently (i.e. 4 for LG4).
    branch_count(),         // Number of probability matrices to be allocated.
    RATE_CATS,              // Number of rate categories we will use.
    inner_nodes_count(),    // How many scale buffers to use.
    PLL_ATTRIB_ARCH_AVX);   // List of flags for hardware acceleration.

  SetModelParameters(partition_);

  EquipPartitionWithData(
    partition_, tree_, tip_nodes_count(), headers, seqdata);

  for(unsigned int i = 0; i < tip_nodes_count(); ++i) {
    free(seqdata[i]);
    free(headers[i]);
  }
  free(seqdata);
  free(headers);


  params_indices_ = (unsigned int *)malloc(RATE_CATS * sizeof(unsigned int));
  for(unsigned int i = 0; i < RATE_CATS; ++i) {
    params_indices_[i] = 0;
  }
  travbuffer_ = (pll_utree_t **)malloc(nodes_count() * sizeof(pll_utree_t *));
  branch_lengths_ = (double *)malloc(branch_count() * sizeof(double));
  matrix_indices_ = (unsigned int *)malloc(branch_count() * sizeof(unsigned int));
  operations_ = (pll_operation_t *)malloc(inner_nodes_count() *
                                                sizeof(pll_operation_t));

}


/// @brief Destructor for Partition.
Partition::~Partition() {
  pll_partition_destroy(partition_);
  free(params_indices_);
  free(travbuffer_);
  free(branch_lengths_);
  free(matrix_indices_);
  free(operations_);
  pll_utree_destroy(tree_);
}


/// @brief Returns a the tree as a Newick string.
/// @return Newick string.
std::string Partition::ToNewick() {
  return std::string(pll_utree_export_newick(tree_));
}


/// @brief Perform a full tree traversal and return the likelihood.
/// Eventually, we will want to break this up for efficiency gains
/// so that we aren't always doing a full traversal every time we
/// want to calculate the likelihood, but we'll do that carefully!
/// @return Log likelihood.
double Partition::FullTraversalLogLikelihood() {

  /* Perform a full postorder traversal of the unrooted tree. */
  unsigned int traversal_size;
  unsigned int matrix_count, ops_count;
  if(!pll_utree_traverse(
      tree_,
      cb_full_traversal,
      travbuffer_,
      &traversal_size)) {
    fatal("Function pll_utree_traverse() requires inner nodes as parameters");
  }

  /* Given the computed traversal descriptor, generate the operations
     structure, and the corresponding probability matrix indices that
     may need recomputing. */
  pll_utree_create_operations(
    travbuffer_,
    traversal_size,
    branch_lengths_,
    matrix_indices_,
    operations_,
    &matrix_count,
    &ops_count);

  /* Update matrix_count probability matrices for model with index 0. The i-th
     matrix (i ranges from 0 to matrix_count - 1) is generated using branch
     length branch_lengths[i] and can be referred to with index
     matrix_indices[i]. */
  pll_update_prob_matrices(
    partition_,
    params_indices_,
    matrix_indices_,
    branch_lengths_,
    matrix_count);

  /* Use the operations array to compute all ops_count inner CLVs. Operations
     will be carried out sequentially starting from operation 0 towrds ops_count-1. */
  pll_update_partials(partition_, operations_, ops_count);

  double logl = pll_compute_edge_loglikelihood(
    partition_,
    tree_->clv_index,
    tree_->scaler_index,
    tree_->back->clv_index,
    tree_->back->scaler_index,
    tree_->pmatrix_index,
    params_indices_,
    NULL); // Can supply a persite_lnl parameter.

  return logl;
}


}
