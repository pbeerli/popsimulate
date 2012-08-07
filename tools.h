#ifndef TOOLS_INCLUDE
#define TOOLS_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
H E L P E R     R O U T I N E S 

 

 Peter Beerli 1996, Seattle
 beerli@genetics.washington.edu
 $Id: tools.h,v 1.1.1.1 1998/06/06 05:57:37 beerli Exp $
-------------------------------------------------------*/

#include "migration.h"

extern double lengthof (node * p);
extern node *crawlback (const node * theNode);
extern node *crawl (node * theNode);
extern node *showtop (node * theNode);
extern void adjust_time (node * theNode, double tyme);
extern void insert_migr_node (world_fmt * world, node * up, node * down,
		     migr_table_fmt * migr_table, long *migr_table_counter);
extern void children (node * mother, node ** brother, node ** sister);
/* math tools */
extern double incompletegamma (double x, double alpha);
extern double polygamma (long n, double z);
extern void invert_matrix (double **a, long nsize);
extern boolean nrcheck (double **m, double **tm, double *v, long nrows, double *r1, double *r2, boolean do_newton);
extern double rannor (double mean, double sd);
extern double find_chi (long df, double prob);
extern double probchi (long df, double chi);
#ifndef HAVE_LGAMMA
extern double lgamma (double z);
#endif
extern double sum (double *vector, long n);
extern char lowercase (char c);
extern char uppercase (char c);

/*filemanipulation */
extern void init_files (world_fmt * world, data_fmt * data, option_fmt * options);
extern void exit_files (world_fmt * world, data_fmt * data, option_fmt * options);
extern void openfile (FILE ** fp, char *filename, char *mode, const char *application, char *perm);
extern void read_savesum(world_fmt *world, option_fmt *options, data_fmt *data);
extern void write_savesum(world_fmt *world);

/* string manipulation */
extern void translate (char *text, char from, char to);
/* time reporting */
extern void get_time (char *nowstr, char ts[]);
/* printing aid */
extern void print_llike (double llike, char *strllike);


#endif /*TOOLS_INCLUDE */
