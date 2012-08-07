



/* --------------------------------------------*
 * sort.h             
 * comparison routines for the different sorts
 * part of the lamarc package   
 * --------------------------------------------*
 * Peter Beerli 1994                            
 */
#ifndef MIGRATION_SORT
#define MIGRATION_SORT
/*-------------------------------------------
 * comparison between characters
 * used in qsort()
 */
int charcmp (const void *v1, const void *v2);

/*-------------------------------------------
 * comparison between strings
 * used in qsort()
 */
int stringcmp (const void *v1, const void *v2);


/*-------------------------------------------
 * comparison between numbers (doubles,longs,
 *  and ints) used in qsort()
 */
int numcmp (const void *v1, const void *v2);
int longcmp (const void *v1, const void *v2);
int intcmp (const void *v1, const void *v2);

/*-------------------------------------------
 * comparison between the times in the
 * the timevector (vtlist) and a given time
 * used in qsort()
 */
int agecmp (const void *x, const void *y);

/*-------------------------------------------
 * comparison in nodep->v 
 * used only for the tipnodes, sorts the
 * tips=2=having '?' alleles
 * to the end of nodelist.
 * used in qsort()
 */
int delcmp (const void *x, const void *y);


/*------------------------------------------
 * compares migr_table time entries
 * used in qsort of the migr_table
 */
int migr_time_cmp (const void *x, const void *y);


/*-------------------------------------------
 * comparison between the times in the
 * the timevector (vtlist) and a given time
 * used in bsearch()
 */
int searchagecmp (const void *x, const void *y);
#endif
