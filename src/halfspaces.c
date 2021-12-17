#define qh_QHimport
#include "qhull_ra.h"

double** intersections(
  double*   halfspaces,
  double*   interiorpoint,
  unsigned  dim,
  unsigned  n,
  unsigned* nintersections,
  unsigned* exitcode,
  unsigned  print
)
{
  char opts[250];
  sprintf(opts, "qhull s Fp FF H H%f", interiorpoint[0]); //, interiorpoint[1] , interiorpoint[2]);
  for(unsigned i=1; i < dim; i++){
    char x[20];
    sprintf(x, ",%f", interiorpoint[i]);
    strcat(opts, x);
  }
  printf(opts); printf("\n");

  qhT qh_qh; /* Qhull's data structure */
  qhT* qh= &qh_qh;
  QHULL_LIB_CHECK
  qh_meminit(qh, stderr);
  boolT ismalloc  = False; /* True if qhull should free points in qh_freeqhull() or reallocation */
  FILE *errfile   = NULL;
  FILE* outfile = print ? stdout : NULL;
  qh_zero(qh, errfile);
  *exitcode = qh_new_qhull(qh, dim+1, n, halfspaces, ismalloc, opts, outfile,
                           errfile);
  printf("exitcode: %u\n", *exitcode);

  double** out;
  if(!(*exitcode)){
    *nintersections = qh->num_facets;
    out = malloc(*nintersections * sizeof(double*));
    facetT *facet;
    unsigned i_facet = 0;
    FORALLfacets{
      if(facet->offset != 0){
        out[i_facet] = malloc(dim * sizeof(double));
        for(unsigned i=0; i < dim; i++){
          out[i_facet][i] = - facet->normal[i] / facet->offset +
                            qh->feasible_point[i]; // = interiorpoint ? yes
        }
        i_facet++;
      }else{
        (*nintersections)--;
      }
    }

  }

  /* Do cleanup regardless of whether there is an error */
  int curlong, totlong;
  qh_freeqhull(qh, !qh_ALL);                /* free long memory */
  qh_memfreeshort(qh, &curlong, &totlong);  /* free short memory and memory allocator */

  if(*exitcode){
    free(out);
    return 0;
  }else{
    return out;
  }

}
