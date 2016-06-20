#include "pll-utils.hpp"
#include <cstdarg>
#include <search.h>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

/// @file pll-utils.hpp
/// @brief Utilities for interfacing with libpll.
/// Much of this was directly copied from libpll examples, so parts (c) Thomas
/// Flouri.
/// Note that here we have a confluence of the pt C++ and the libpll C.
/// Thus this is the only place where we have mixed the two style conventions.


namespace pt {

// if you remove this, you can remove the include for cstdarg above.
void fatal(const char *format, ...) {
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}


/* a callback function for performing a full traversal */
int cb_full_traversal(pll_utree_t *node) { return 1; }

static int utree_traverse_check(pll_utree_t * node,
                                 int (*cbtrav)(pll_utree_t *))
{
  if (!node->next)
  {
    return cbtrav(node);
  }
  if (!cbtrav(node))
    return 0;

  return utree_traverse_check(node->next->back, cbtrav) &&
         utree_traverse_check(node->next->next->back, cbtrav);
}


/// @brief Determine if a tree is healthy, i.e. has branch lengths.
/// @param[in] tree
/// A pll_utree_t to check.
int cb_branch_healthy(pll_utree_t *tree)
{
if (!tree->length) return 0;

  if (tree->next) {
    if (!tree->next->length || !tree->next->next->length) return 0;

    return (cb_branch_healthy(tree->next->back) &&
            cb_branch_healthy(tree->next->next->back));
  }

  return true;

}
bool TreeHealthy(pll_utree_t *tree) {
  return (utree_traverse_check(tree,cb_branch_healthy) && utree_traverse_check(tree,cb_branch_healthy));
}


/// @brief Parse a Fasta file.
/// @param[in] path
/// A Fasta path.
/// @param[in] seq_count
/// How many sequences are expected.
/// @param[out] headers_out
/// A pointer to a char**, which will be filled with an array of header strings
/// that will need to be freed.
/// @param[out] seqdata_out
/// A pointer to a char**, which will be filled with an array of sequence
/// strings that will need to be freed.
unsigned int ParseFasta(std::string path, unsigned int seq_count,
                        char ***headers_out, char ***seqdata_out) {
  pll_fasta_t *fp = pll_fasta_open(&path[0], pll_map_fasta);
  if (!fp) fatal("Error opening file %s", &path[0]);

  char *seq = NULL;
  char *hdr = NULL;
  long seqlen;
  long hdrlen;
  long seqno;

  /* allocate arrays to store FASTA headers and sequences */
  char **headers = (char **)calloc(seq_count, sizeof(char *));
  char **seqdata = (char **)calloc(seq_count, sizeof(char *));

  /* read FASTA sequences and make sure they are all of the same length */
  unsigned int i;
  int sites = -1;
  for (i = 0; pll_fasta_getnext(fp, &hdr, &hdrlen, &seq, &seqlen, &seqno);
       ++i) {
    if (i >= seq_count)
      fatal("FASTA file contains more sequences than expected.\n");

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

  if (sites < 0) fatal("Unable to read alignment");

  if (i != seq_count) fatal("Some taxa are missing from FASTA file");

  *headers_out = headers;
  *seqdata_out = seqdata;
  return (unsigned int)sites;
}


/// @brief Add tip data to a partition.
/// Look up every sequence by its name and set the corresponding CLV in the
/// partition.
/// @param[in] partition
/// A partition to equip.
/// @param[in] tree
/// The tree to use for equipping the partition.
/// @param[in] tip_nodes_count
/// The number of tip nodes.
/// @param[out] headers
/// An array of header strings.
/// @param[out] seqdata
/// An array of sequence strings.
void EquipPartitionWithData(pll_partition_t *partition, pll_utree_t *tree,
                            unsigned int tip_nodes_count, char **headers,
                            char **seqdata) {
  /*  obtain an array of pointers to tip nodes */
  pll_utree_t **tipnodes =
      (pll_utree_t **)calloc(tip_nodes_count, sizeof(pll_utree_t *));
  pll_utree_query_tipnodes(tree, tipnodes);

  /* create a libc hash table of size tip_nodes_count */
  hcreate(tip_nodes_count);

  /* populate a libc hash table with tree tip labels */

  unsigned int *data =
      (unsigned int *)malloc(tip_nodes_count * sizeof(unsigned int));
  unsigned int i;
  for (i = 0; i < tip_nodes_count; ++i) {
    data[i] = i;
    ENTRY entry;
    entry.key = tipnodes[i]->label;
    entry.data = (void *)(data + i);
    hsearch(entry, ENTER);
  }

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < tip_nodes_count; ++i) {
    ENTRY query;
    query.key = headers[i];
    ENTRY *found = NULL;

    found = hsearch(query, FIND);

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

std::vector<std::string> ssplit(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> tokens;
    while (getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}


/// @brief Set up the model parameters of the given partition.
/// @todo This should take a RAxML_info file, parse it, and insert those values
/// here.
void SetModelParameters(pll_partition_t *partition, std::string path) {
  /* initialize the array of base frequencies */
  /// @todo `Base frequencies` should go here
  std::ifstream file(path);
  std::string read;
  std::string contents;

  while(std::getline(file,read))
  {
    contents+=read;
    contents.push_back('\n');
   }

  std::size_t pos1= contents.find("frequencies: ");
  std::size_t pos2= contents.find('\n',pos1);

  std::string sstr=contents.substr(pos1+13,pos2-pos1-13);

  std::vector<std::string> freqvector=ssplit(sstr,' ');
  double frequencies[freqvector.size()];

  for(int i=0; i<freqvector.size(); i++)
    frequencies[i]=std::stod(freqvector.at(i));

  ///double frequencies[4] = {0.17, 0.19, 0.25, 0.39}; OLD

  /* substitution rates for the 4x4 GTR model. This means we need exactly
     (4*4-4)/2 = 6 values, i.e. the number of elements above the diagonal */
  /// @todo `rates` should go here
  pos1= contents.find("ac ag at cg ct gt:");
  pos2= contents.find('\n',pos1);
  sstr=contents.substr(pos1+19,pos2-pos1-19);
  std::vector<std::string> ratevector=ssplit(sstr,' ');
  double subst_params[ratevector.size()];

  for(int i=0; i<ratevector.size(); i++)
    subst_params[i]=std::stod(ratevector.at(i));
  if(ratevector.size()!=(((freqvector.size())*freqvector.size() -4)/2))
    fatal("Wrong number of rate values.");


/*for(int i=0;i<6;i++)
  subst_params[i] = 1;*/

  /* we'll use RATE_CATS rate categories, and currently initialize them to 0 */
  double rate_cats[RATE_CATS] = {0};

  /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
  /// @todo `alpha` should go here
  pos1= contents.find("alpha[0]: ");
  pos2= contents.find(' ', pos1+10);
  sstr=contents.substr(pos1+10,pos2-pos1-10);
  double alpha= stod(sstr);
  pll_compute_gamma_cats(alpha, RATE_CATS, rate_cats);


  file.close();

  /* set frequencies at model with index 0 (we currently have only one model) */
  pll_set_frequencies(partition, 0, frequencies);

  /* set 6 substitution parameters at model with index 0 */
  pll_set_subst_params(partition, 0, subst_params);

  /* set rate categories */
  pll_set_category_rates(partition, rate_cats);
}


}
