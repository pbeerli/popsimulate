#ifndef MIGRATION_RANDOM
#define MIGRATION_RANDOM
/* -------------------------------------------------------                        
   R A N D O M   G E N E R A T O R   R O U T I N E S 

   creates options structures,
   reads options from parmfile if present

   prints options,
   and finally helps to destroy itself.

   Peter Beerli 1996, Seattle
   beerli@genetics.washington.edu
   $Id: random.h,v 1.2 2002/05/16 05:03:03 beerli Exp $
   ------------------------------------------------------- */

/*-----------------------------------------------------*
 * calculates the start seed for randum
 * picks up the global variable iseed
 * PB 94
 *-----------------------------------------------------*/
extern void getseed (option_fmt * options);

/*-----------------------------------------------------*
 * RANDUM 
 * generates an uniform random number between o..1
 * using seed (generated by getseed)
 * JF <93 (6bit) and Mary K. Kuhner 1997 (32bit)
 * RANDINT
 * generates an uniform random integer number
 * between a..b
 * using seed (generated by getseed)
 *-----------------------------------------------------*/
extern double mk_randum (void);
#define RANDUM() mk_randum()
#define RANDINT(a,b) (long)((a)+(mk_randum()* (((b)-(a))+1.)))

#endif
