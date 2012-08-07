/// generates a coalescent genealogy and additional
// parameter settings
//
// needs parmfile piped into
// syntax parmfile is
//
// datatype mutation_rate_type
// numpop numloci simseed1 <simseed2 simseed3>
// sample1 <sample2 ....>
// popsize1 <popsize2 ....>
// migration rate matrix, diagonals are -
// for sequences:
// mutation_rate <locus_alpha> <site_alpha>
// <numsites ttratio> <type of tree output (with-mig | no-mig (default))>
// for msat:
// flanking_regio1 flanking_regio2 ....
// replength1 replength2 ....
//
// datatype = S_equence R_ratevariation_among_sites
//      and many others like M and N etc
#undef HIGHBITS
#include "definitions.h"
#include "migration.h"
#include "random.h"
#include <string.h>
#include <stdio.h>
#include <sys/time.h>
//#ifdef SGI
#include "ranlib.h"
//#else
//#include <ranlib.h>
//#endif
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

typedef struct param_fmt
  {
    char type;
    char dummy;
    char gamma[2];
    double ttratio;
    long *sites;
    double *param;
    double mu;
    double *allparam;
    long numpop;
    long numpop2;
    long sumtips;
    long *savesample;
    long *sample;
    long loci;
    long unique;
    long simseed;
    long gsimseed1;
    long gsimseed2;
    node **tips;
    double tyme;
    double oldtyme;
    long mcount1;
    long mcount2;
    long ccount;
    boolean withmig;
    long locus;
  }
param_fmt;
long r1, r2, r3, r4;
long shape_mut;
long bug;
boolean het;

long **migmat;
long simseed[6];

double *inheritance;
long *flanking;
long *replength;

boolean has_msatlength;
char input[LINESIZE];

void jumble (long *s, long n);
node *showtop (node * theNode);
void dtree (node * theNode);
double lengthof (node * p);
node *crawlback (const node * theNode);
void initparam (param_fmt * param);
void resetparam (param_fmt * param);
void freetree (node * p);
void createtips (param_fmt * param);
void maketree (param_fmt * param);
double event (param_fmt * param);
double migrations (long pop, long offset, param_fmt * param);
double coalescences (long pop, param_fmt * param);
double probmigration (param_fmt * param);
long  migrate (long pop, param_fmt * param);
long choosepopmig (param_fmt * param);
long chooseindm (long pop, param_fmt * param);
long lookup (long ind, long pop, param_fmt * param);
void create_migrnode (long ind, long pop, long from, param_fmt * param);
node *initmigrnode (long pop, long from,param_fmt * param);
void coalesce (long pop, param_fmt * param);
long choosepopc (param_fmt * param);
void chooseindc (long *i1, long *i2, long pop, param_fmt * param);
void create_coalnode (long ind1, long ind2, long pop, param_fmt * param);
node *initinternalnode (long pop, param_fmt * param);
double mk_randum (void);
void getsimseed2 (long simseed);
void decideandchange (param_fmt * param);
void treeout (node * joint, node * p, long s, boolean withmig);
double denom (param_fmt * param);
void translate (char *text, char from, char to);
void mutation (param_fmt * param);
long fromwhichpop(long pop, param_fmt *param);
void read_item(double *ptr);

int
main (int argc, char *argv[])
{
  long i, j,nr, pop;
  long locus;
  float alpha=1.;

  param_fmt *param;
  //double cv;
  boolean fullout=FALSE;
  het=FALSE; //globally defined
  param = (param_fmt *) calloc (1, sizeof (param_fmt));

  initparam (param);
  migmat = (long **) calloc (param->numpop, sizeof (long *));
  for(i=0;i<param->numpop;i++)
    migmat[i] = (long *) calloc (param->numpop, sizeof (long));
  if(toupper(param->type)=='J')
    {
      param->type='I';
      fullout=TRUE;
    }
  if(toupper(param->type)=='P')
    {
      param->type = 'L';
      fullout=TRUE;
    }
  if(toupper(param->type)=='R')
    {
      param->type = 'S';
      het = TRUE;
      fullout=TRUE;
    }
  getsimseed2 (param->simseed);

  if(fullout && param->type=='L')
    fprintf (stdout, "#%c%c\n", 'P', toupper (param->gamma[0]));
  else
    fprintf (stdout, "#%c%c\n", fullout ? 'J' : toupper (param->type),
             toupper (param->gamma[0]));

  if (toupper (param->gamma[0]) == 'G' || het)
    {
      setall (param->gsimseed1, param->gsimseed2);// from ranlib sets rand-generator
      for(i=0;i<100;i++)
        gengam(1.,1.); //clean generator

      fprintf (stdout, "#%li %li %li\n", /*RANDINT (0, MAXLONG)*/ param->simseed,
               param->gsimseed1, param->gsimseed2);
    }
  else
    {
      fprintf (stdout, "#%li\n", param->simseed);
    }
  fprintf (stdout, "#%li\n#", param->numpop);
  for (pop = 0; pop < param->numpop; pop++)
    {
      if (toupper (param->type) != 'S' && toupper (param->type) != 'I' && toupper(param->type) != 'N' &&  toupper(param->type) != 'L' && toupper(param->type) != 'P')
        fprintf (stdout, "%li ", param->savesample[pop]/2);
      else
        fprintf (stdout, "%li ", param->savesample[pop]);
    }
  fprintf (stdout, "\n");
  if (toupper (param->type) == 'S' ||toupper (param->type) == 'I'  ||  toupper(param->type) == 'N' || toupper(param->type) == 'L' || toupper(param->type) == 'P')
    {
      fprintf (stdout, "#%li ", param->loci);

      /*  for(locus=0;locus<param->loci;locus++) */
      fprintf (stdout, "%li ", param->sites[0]);
      fprintf (stdout, "%f\n", param->ttratio);
      if (het && toupper (param->gamma[0]) == 'G')
        {
          alpha = (float) param->allparam[param->numpop2 + 2];
        }
      else
        {
          if(het)
            alpha = (float) param->allparam[param->numpop2 + 1];
        }

      for(locus=0;locus<param->loci;locus++)
        {
          fprintf (stdout, "# rate among sites for locus %li (%f)\n# ",locus,alpha);
          for(i=0;i<param->sites[0];i++)
            {
              if(het)
                {
                  fprintf(stdout,"%f ",gengam(alpha,alpha));
                }
              else
                {
                  fprintf(stdout,"1. ");
                }
            }
          fprintf(stdout,"\n");
        }
    }
  else
    {
      fprintf (stdout, "#%li %i\n", param->loci, has_msatlength);
      if(has_msatlength)
	{
	  fprintf (stdout, "# ");
	  for(locus=0;locus<param->loci;locus++)
	    {
	      fprintf(stdout,"%li %li ",flanking[locus], replength[locus]);
	    }
	  fprintf(stdout,"\n");
	}
    }
  // coefficient of variation
  //cv = param->type=='G' ? param->allparam[param->numpop2 + 1] : 1.;
  for (locus = 0; locus < param->loci; locus++)
    {
      param->locus = locus;
      mutation (param);
      createtips (param);
      maketree (param);
      fprintf (stdout, "#$ locus %li\n",locus + 1);
      fprintf (stdout, "#$ %20.20lf\n",param->mu);
      treeout (param->tips[0], param->tips[0], 1, param->withmig);
      nr = 1;
      //freetree (param->tips[0]);
      fprintf(stderr,"Genealogies with migration for %li POPULATION\n\nNumber of migration events [columns=from, rows=to]:\n",param->numpop);
      for(i=0;i<param->numpop;i++)
        {
          for(j=0;j<param->numpop;j++)
            {
              fprintf(stderr,"%li ",migmat[i][j]);
            }
          fprintf(stderr,"\n");
          memset(migmat[i],0,sizeof(long)*param->numpop);
        }

      resetparam (param);
    }
  return 0;
}

void
mutation (param_fmt * param)
{
  /*gammadeviate */
  /*gengam(mean,shape) */
  static long count=0;
  long i,pop;
  float beta, alpha, cv, mean, mu;
  mean = (float) param->allparam[param->numpop2];
  fprintf(stderr,"Mean mutation rate is %f\n",mean);
  if (toupper (param->gamma[0]) == 'G')
    {
      alpha = (float) param->allparam[param->numpop2 + 1];
      //      alpha = 1. / (cv * cv);
      beta = (float) alpha/mean;// CAVEAT: beta for ranlib is 1/beta in
      // mathematica!!!!!!! There are two different interpretation
      // of beta depending on what formula one uses.
      mu = (float) gengam (beta, alpha);
    }
  else
    {
      //cv = 99999999;
      alpha = 1e37;
      beta = 1e37;
      mu = mean;
    }
  //  fprintf(stdout,"# mutation rate = %f (alpha=%f)\n", mu, alpha);
  /* create parameter */
  for (pop = 0; pop < param->numpop; pop++)
    {
      param->param[pop] = inheritance[param->locus] * param->allparam[pop] * mu;
    }
  for(i=param->numpop;i<param->numpop2;i++)
    {
      param->param[i] = param->allparam[i] / mu;
    }
  param->mu = mu;
}

void
freetree (node * p)
{
  if (p->type == 'm')
    {
      freetree (p->next->back);
      free (p->next);
      free (p);
      return;
    }
  if (p->type == 'i')
    {
      freetree (p->next->back);
      freetree (p->next->next->back);
      free (p->next->next);
      free (p->next);
      free (p);
      return;
    }
  if (p->type == 't')
    {
      free (p->nayme);
      free (p);
      return;
    }
}


///
/// better random_seed function taken from a discussion on 
/// http://sourceware.org/ml/gsl-discuss/2004-q1/msg00071.html
/// March 30, 2007, Robert G. Brown (http://www.phy.duke.edu/~rgb/)
/// at Duke University Dept. of Physics, Box 90305  Durham, N.C. 27708-030
/// suggested the code 
unsigned long int random_seed()
{
  
  unsigned int seed;
  struct timeval tv;
  FILE *devrandom;
  
  if ((devrandom = fopen("/dev/urandom","r")) == NULL) 
    {
      gettimeofday(&tv,0);
      seed = tv.tv_sec + tv.tv_usec;
    }
  else 
    {
      fread(&seed,sizeof(seed),1,devrandom);
      fclose(devrandom);
    }
  return(seed);
}


void
initparam (param_fmt * param)
{
  long pop;
  char tmp;
  param->sumtips = 0;
  param->unique = 1;
  param->type = toupper(getchar ());
  if(param->type=='R')
    {
      het=TRUE;
      param->type='S';
    }
  scanf ("%1s%*[^\n]", param->gamma);
  if (het || toupper (param->gamma[0]) == 'G')
    {
      scanf ("%li %li %li %li %li%*[^\n]", &param->numpop, &param->loci, &param->simseed, &param->gsimseed1, &param->gsimseed2);
      shape_mut = 1+(het ? 1 : 0)+ (toupper (param->gamma[0]) == 'G' ? 1 : 0);
      if(param->simseed <=0)
	param->simseed = random_seed();
      if(param->gsimseed1 <=0)
	param->gsimseed1 = random_seed();
      if(param->gsimseed2 <=0)
	param->gsimseed2 = random_seed();
    }
  else
    {
      scanf ("%li %li %li%*[^\n]", &param->numpop, &param->loci,
             &param->simseed);
      shape_mut = 1;
      if(param->simseed <=0)
	param->simseed = random_seed();
    }
  param->numpop2 = param->numpop * param->numpop;
  param->sites = (long *) calloc (1, sizeof (long) * param->loci);
  param->sample = (long *) calloc (1, sizeof (long) * param->numpop);
  param->savesample = (long *) calloc (1, sizeof (long) * param->numpop);
  param->param = (double *) calloc (1, sizeof (double) * param->numpop2);
  param->allparam = (double *) calloc (1, sizeof (double) * (param->numpop2
                                       + shape_mut));
  for (pop = 0; pop < param->numpop; pop++)
    {
      scanf ("%li", &param->savesample[pop]);
      if(param->type!='S' &&param->type!='I' && param->type!='J' && toupper(param->type)!='N' &&  toupper(param->type) != 'L' && toupper(param->type) != 'P')
        param->savesample[pop] *= 2;
      param->sumtips += param->savesample[pop];
    }
  // this code assumes that the migration diagonal is - or another character that
  // is not a number or . so there are numpop slots for sizes and numpop(numpop-1) for migrations
  for (pop = 0; pop < param->numpop2 + shape_mut; pop++)
    {
      read_item (&param->allparam[pop]);
    }
  long l;
  //inheritance is default 4, mtdna=1 etc
  inheritance = (double*) calloc(param->loci,sizeof(double));
  for (l=0;l<param->loci;l++)
    {
      inheritance[l]=4.0;
    }

  if (param->type == 'S' ||param->type == 'J' ||param->type == 'I' || param->type == 'N' || toupper(param->type) == 'L' || toupper(param->type) == 'P')
    {
      /*       for (pop=0;pop<param->loci;pop++){ */
      scanf ("%li", &param->sites[0]);
      for (pop = 1; pop < param->loci; pop++)
        {
          param->sites[pop] = param->sites[0];
        }
      scanf ("%lf", &param->ttratio);
      /*       param->ttratio=2.0; */
    }
  else
    {
      char opt='_';
      scanf("%c",&opt);
      if(opt=='f')
	{
	  tmp ='\0';
	  int z=0;
	  has_msatlength=TRUE;
	  flanking = (long*) calloc(param->loci,sizeof(long));
	  for(z=0; z < param->loci; z++)
	    {
	      int ret = scanf("%li", &flanking[z]); // read text this should be the flanking region
	      fprintf(stderr,"%li\n",flanking[z]);
	      if (ret<1)
		{
		  fprintf(stderr,"Missing flanking region!\n");
		  exit(-1);
		}
	    }
	  replength= (long*) calloc(param->loci,sizeof(long));
	  for(z=0; z < param->loci; z++)
	    {
	      int ret = scanf("%li", &replength[z]); // read text this should be the flanking re
	      if (ret<1)
		{
		  fprintf(stderr,"Missing repeatlength!\n");
		  exit(-1);
		}
	    }
	}
      else
	{
	  ungetc(opt,stdin);
	}
    }
  tmp =' ';
  while(tmp==' ')
    tmp = fgetc(stdin); // read blank
  //tmp is now something: this is not necessary tmp = fgetc(stdin); // read ??
  param->withmig=FALSE;
  int locus;
  //  while(tmp!= '\0')
  while(tmp > 0)
    {
      fprintf(stderr,"%c\n",tmp);
      switch(tmp){
      case 'w':
	param->withmig=TRUE;
	fscanf(stdin,"%*[^\n]");
	tmp = fgetc(stdin); //reads end of line character
	tmp = fgetc(stdin);
	break;
      case 'i':
	for (locus=0;locus<param->loci;locus++)
	  {
	    fscanf(stdin,"%lf",&inheritance[locus]);
	  }
	fgets(input,LINESIZE,stdin);
	tmp='\0';
	break;
      case 'n': /*not with migrate nodes*/
	param->withmig=FALSE;
	fscanf(stdin,"%*[^\n]");
	tmp = fgetc(stdin); //reads end of line character
	tmp = fgetc(stdin);
	break;
      case 'e':
      default:
	tmp = '\0';
      }
    }
  
  memcpy (param->sample, param->savesample, sizeof (long) * param->numpop);
  param->tips = (node **) calloc (1, sizeof (node *) * (param->sumtips * 2 + 1));
}

void read_item(double *ptr)
{
  long z=0;
  char ch;
  char val[1024];
  boolean reading=TRUE;
  ch=getchar();
  while(reading)
    {
      if(strchr("0123456789.",ch))
        {
          while(strchr("0123456789.",ch))
            {
              val[z++] = ch;
              if(z>1024)
                exit(-1);
              ch=getchar();
            }
          reading =FALSE;
        }
      else
        {
          ch=getchar();
        }
    }
  val[z] = '\0';
  fprintf(stderr,"%s=%f\n",val,atof(val));
  *ptr = atof(val);
}

void
resetparam (param_fmt * param)
{
  long pop;
  param->sumtips = 0;
  param->unique = 1;
  param->tyme = 0;
  memcpy (param->sample, param->savesample, sizeof (long) * param->numpop);
  for (pop = 0; pop < param->numpop; pop++)
    {
      param->sumtips += param->sample[pop];
    }
  param->tips = (node **) realloc (param->tips, sizeof (node *) * (param->sumtips * 2 + 1));
  memset (param->tips, 0, sizeof (node *) * (param->sumtips * 2 + 1));
}


void
createtips (param_fmt * param)
{
  char achar = 'A';
  char bchar = 'A';
  char cchar = 'A';
  char sumchar[4];
  long asum = 1;
  long bsum = 1;
  node *atip;
  long i, pop, ind, sumind = 0;
  long *look;
  look = (long *) calloc (1, sizeof (long) * param->numpop);
  for (pop = 0; pop < param->numpop; pop++)
    look[pop] = pop;
  jumble (look, param->numpop);
  param->ccount = 0;
  param->mcount1 = 0;
  param->mcount2 = 0;
  for (i = 0; i < param->numpop; i++)
    {
      achar = bchar = cchar = 'A';
      asum=bsum=0;

      pop = look[i];
      for (ind = 0; ind < param->sample[pop]; ind++)
        {
          atip = (node *) calloc (1, sizeof (node));
          atip->actualpop = pop;
          atip->pop = pop;
          atip->type = 't';
	  atip->top = 1;
          atip->id = param->unique++;
          atip->nayme = (char *) calloc (1, sizeof (char) * 5);
          sprintf (atip->nayme, "%li", pop);
          if (asum++ % 576 == 0)
            {
              achar += 1;
              asum = 1;
              bchar = 'A';
              bsum = 1;
              cchar = 'A';
            }
          else
            {
              if (bsum++ % 24 == 0)
                {
                  bchar += 1;
                  bsum = 1;
                  cchar = 'A';
                }
            }

          sumchar[0] = achar;
          sumchar[1] = bchar;
          sumchar[2] = cchar++;
          sumchar[3] = '\0';
          strcat (atip->nayme, sumchar);
          param->tips[sumind++] = atip;
        }
    }
  free (look);
}

void
maketree (param_fmt * param)
{

  while (param->sumtips > 1)
    {
      param->oldtyme = param->tyme;
      param->tyme += event (param);
      //fprintf(stderr,"%f %li ",param->tyme, param->sumtips);
      decideandchange (param);
    }
  fprintf(stderr,"Last time in tree=%f, Total (?) N=%f\n", param->tyme, param->tyme/(param->allparam[param->numpop2]*4.));
}

double
event (param_fmt * param)
{
  return (-log (RANDUM ()) / denom (param));
  /*  return (1./(denom(param))); */
}


double
denom (param_fmt * param)
{
  long i,pop,offset;
  double sm, sum = 0;
  //fprintf(stderr,"INHERITENCE: %f %li\n", inheritance[param->locus],param->locus);
  for (pop = 0; pop < param->numpop; pop++)
    {
      sum += ((double) param->sample[pop] * (param->sample[pop] - 1)) / (param->param[pop]);
      offset = param->numpop + pop *(param->numpop-1);
      sm=0.0;
      for(i=offset;i<offset+param->numpop-1;i++)
        sm +=  param->param[i];
      sum += param->sample[pop] * sm;
    }
  return sum;
}


void
decideandchange (param_fmt * param)
{
  long from= -1;
  double *r;
  long *s;
  long pop,j, offset;
  double rsum;
  double rr = RANDUM ();
  double dd = denom (param);
  r = (double *) calloc (1, sizeof (double) * (1+param->numpop2));
  for (pop = 0; pop < param->numpop; pop++)
    {
      r[pop] = coalescences (pop, param) / dd;
      offset = param->numpop + (param->numpop-1)*pop;
      r[param->numpop+pop] = migrations(pop,offset,param) /dd;
    }

  for (pop = 1; pop < (2*param->numpop); pop++)
    {
      r[pop] += r[pop - 1];
    }
  rsum = r[pop-1];
  pop = 0;
  while (rr >= r[pop]/rsum)
    {
      //    fprintf(stderr,"pop=%li %f %f\n",pop,r[pop]/rsum,rr);
      pop++;
    }
  if (pop >= param->numpop)
    {
      from = migrate (pop-param->numpop, param);
      fprintf(stderr, "m %20.20f %li %li %li %010li\n", param->tyme, from, pop-param->numpop, param->sumtips, param->locus);
    }
  else
    {
      fprintf(stderr,"c %20.20f %li %li %li %010li (%f)\n", param->tyme, pop, pop, param->sumtips, param->locus,
 param->sumtips * (param->sumtips-1) * (param->tyme - param->oldtyme));
      coalesce (pop, param);
    }
  free(r);

}


double
migrations (long pop, long offset, param_fmt * param)
{
  long i;
  double sm=0.0;
  for(i=offset;i<offset+param->numpop-1;i++)
    sm += param->param[i];
  return ((double) param->sample[pop]) * sm;
}

double
coalescences (long pop, param_fmt * param)
{
  return ((double) param->sample[pop] * (param->sample[pop] - 1)) / (param->param[pop]);
}

long
migrate (long pop, param_fmt * param)
{
  long ind, from;
  ind = chooseindm (pop, param);
  from = fromwhichpop(pop,param);
  migmat[pop][from] +=1;
  create_migrnode (ind, pop, from, param);
  return from;
}

// known
// M_ij
//
//
//
//
//
long fromwhichpop(long pop, param_fmt *param)
{
  long i,z=1;
  double *r, rr, rsum=0.0;
  long offset = param->numpop + pop* (param->numpop-1);
  r = (double*) calloc(1,sizeof(double)*(param->numpop));
  r[0] = param->param[offset];
  for(i=offset +1;i<offset+param->numpop-1;i++)
    {
      r[z] = r[z-1] + param->param[i];
      z++;
    }
  rr = RANDUM();
  for(i=0;i<z;i++)
    {
      r[i] /= r[z-1];
      if(r[i]>rr)
        {
          free(r);
          return ((i<pop) ? i : 1 + i);
        }
    }
  free(r);

  return-1; // never do this ((z<=pop) ? z-1 :  z);
}

long
chooseindm (long pop, param_fmt * param)
{
  long ind = RANDINT (0, param->sample[pop] - 1);
  return lookup (ind, pop, param);
}


long
lookup (long ind, long pop, param_fmt * param)
{
  long i;

  if (ind == 0)
    {
      for (i = 0; i < param->sumtips; i++)
        {
          if (param->tips[i]->pop == pop)
            return i;
        }
    }
  else
    {
      for (i = 0; i < param->sumtips; i++)
        {
          if (param->tips[i]->pop == pop)
            {
              if (ind == 0)
                return i;
              else
                ind--;
            }
        }
    }
  fprintf (stderr, "An error occured in lookup, %li in pop %li not found\n", ind, pop);
  return (-1);
}

void
create_migrnode (long ind, long pop, long from, param_fmt * param)
{
  node *mnode = initmigrnode (pop, from, param);
  mnode->next->back = param->tips[ind];
  param->tips[ind]->back = mnode->next;
  param->tips[ind] = mnode;
  param->sample[pop] -= 1;
  param->sample[from] += 1;
}

node *
initmigrnode (long pop, long from, param_fmt * param)
{
  node *n1, *n2;
  n1 = (node *) calloc (1, sizeof (node));
  n2 = (node *) calloc (1, sizeof (node));
  n1->id = param->unique;
  n2->id = param->unique++;
  n1->next = n2;
  n2->next = n1;
  n1->type = 'm';
  n2->type = 'm';
  n1->tyme = param->tyme;
  n2->tyme = param->tyme;
  n1->top = 1;
  n2->top = 0;
  n1->actualpop = pop;
  n1->pop = from;
  n2->actualpop = pop;
  n2->pop = n1->pop;
  return n1;
}

void
coalesce (long pop, param_fmt * param)
{
  long ind, ind1, ind2;
  chooseindc (&ind1, &ind2, pop, param);
  param->ccount++;
  if (ind1 > ind2)
    {
      ind = ind2;
      ind2 = ind1;
      ind1 = ind;
    }
  create_coalnode (ind1, ind2, pop, param);
}

void
chooseindc (long *i1, long *i2, long pop, param_fmt * param)
{
  long ind2, ind1 = RANDINT (0, param->sample[pop] - 1);
  ind2 = RANDINT (0, param->sample[pop] - 1);
  while (ind2 == ind1)
    ind2 = RANDINT (0, param->sample[pop] - 1);
  *i1 = lookup (ind1, pop, param);
  *i2 = lookup (ind2, pop, param);
}

void
create_coalnode (long ind1, long ind2, long pop, param_fmt * param)
{
  node *n1, *n2;
  node *cnode = initinternalnode (pop, param);
  cnode->next->back = param->tips[ind1];
  param->tips[ind1]->back = cnode->next;
  cnode->next->next->back = param->tips[ind2];
  param->tips[ind2]->back = cnode->next->next;
  param->tips[ind1] = cnode;
  param->tips[ind2] = param->tips[param->sumtips - 1];
  n1 = crawlback (param->tips[ind1]->next);
  n2 = crawlback (param->tips[ind1]->next->next);
  n1->length = lengthof (n1);
  n2->length = lengthof (n2);
  param->sumtips -= 1;
  param->sample[pop] -= 1;
}

node *
initinternalnode (long pop, param_fmt * param)
{
  node *n1, *n2, *n3;
  n1 = (node *) calloc (1, sizeof (node));
  n2 = (node *) calloc (1, sizeof (node));
  n3 = (node *) calloc (1, sizeof (node));
  n1->id = param->unique;
  n2->id = param->unique;
  n3->id = param->unique++;
  n1->next = n2;
  n2->next = n3;
  n3->next = n1;
  n1->type = 'i';
  n2->type = 'i';
  n3->type = 'i';
  n1->tyme = param->tyme;
  n2->tyme = param->tyme;
  n3->tyme = param->tyme;
  n1->top = 1;
  n2->top = 0;
  n3->top = 0;
  n1->actualpop = pop;
  n2->actualpop = pop;
  n3->actualpop = pop;
  n1->pop = pop;
  n2->pop = pop;
  n3->pop = pop;
  return n1;
}


void
getsimseed2 (long insimseed)
{
  long i;
  double clearsimseed;
  for (i = 1; i <= 1000; i++)	/* clear the random numbergenerator */
    clearsimseed = RANDUM ();
  for (i = 0; i <= 5; i++)
    simseed[i] = 0;
  i = 0;
  do
    {
      simseed[i] = insimseed & 63;
      insimseed /= 64;
      i++;
    }
  while (insimseed != 0);
}


double
mk_randum (void)
{				/* randum -- slow but machine independent */
  /* random number generator -- slow but machine independent */
  long i, j, k, sum;
  long mult[6], newsimseed[6];
  double x;

  mult[0] = 13;
  mult[1] = 24;
  mult[2] = 22;
  mult[3] = 6;
  for (i = 0; i <= 5; i++)
    newsimseed[i] = 0;
  for (i = 0; i <= 5; i++)
    {
      sum = newsimseed[i];
      k = i;
      if (i > 3)
        k = 3;
      for (j = 0; j <= k; j++)
        sum += mult[j] * simseed[i - j];
      newsimseed[i] = sum;
      for (j = i; j <= 4; j++)
        {
          newsimseed[j + 1] += newsimseed[j] / 64;
          newsimseed[j] &= 63;
        }
    }
  memcpy (simseed, newsimseed, sizeof (long)*6);
  simseed[5] &= 3;
  x = 0.0;
  for (i = 0; i <= 5; i++)
    x = x / 64.0 + simseed[i];
  x /= 4.0;
  return x;
}				/* randum */

void
translate (char *text, char from, char to)
{
  int i, j, gap = 0;
  while (text[gap] == from)
    gap++;
  for (i = gap, j = 0; text[i] != '\0'; i++)
    {
      if (text[i] != from)
        {
          text[j++] = text[i];
        }
      else
        {
          if (text[i - 1] != from)
            {
              text[j++] = to;
            }
        }
    }
  text[j] = '\0';
}




void
treeout (node * joint, node * p, long s, boolean withmig)
{
  /* write out file with representation of final tree */
  long w;
  double x;
  char migstring[30];
  if (p->type == 't')
    {
      translate (p->nayme, ' ', '_');
      fprintf (stdout, "%s", p->nayme);
    }
  else
    {
      putc ('(', stdout);
      treeout (joint, crawlback (p->next), s, withmig);
      putc (',', stdout);
      treeout (joint, crawlback (p->next->next), s, withmig);
      putc (')', stdout);
    }
  if (p == joint)
    {
      x = 0.0;
    }
  else
    {
      x = crawlback(p)->tyme - p->tyme;
    }
  if (x > 0.0)
    {
      w = (long) (0.4343 * log (x));
    }
  else
    {
      if (x == 0.0)
        w = 0;
      else
        w = (long) (0.4343 * log (-x)) + 1;
    }
  if (w < 0)
    w = 0;
  fprintf (stdout, ":%*.10f", (int) (w + 7), x /*< 10000. ? x : 10000.*/);
  if(p!=joint)
    {
      p = showtop(p->back);
      while(p->type == 'm')
        {
	  if(withmig)
	    {
	      sprintf (migstring, " [&M %li %li:%g]",
		       p->pop, p->actualpop, p->tyme - showtop(p->next->back)->tyme);
	      fprintf(stdout,"%s",migstring);
	    }
          p = showtop(p->back);
        }
    }
  else
    {
      fprintf (stdout,";\n");
    }
}				/* treeout */

double
lengthof (node * p)
{
  if (p->type == 'm')
    fprintf (stderr, "a migration node was feed into lengthof");
  return fabs (p->tyme - crawlback (p)->tyme);
}				/* length */

node *
crawlback (const node * theNode)
{
  node *tmp = theNode->back;

  while (tmp->type == 'm')
    {
      tmp = tmp->next->back;
    }
  return tmp;
}



void
dtree (node * theNode)
{
  /* pass atr (*tree).root->back for normal results */
  node *tmp;
  if (bug)
    {
      if ((theNode->type == 'i' || theNode->type == 'm') && theNode->top)
        {
          if (theNode->type == 'm')
            fprintf (stdout, "-M%li-", theNode->id);
          else
            {
              tmp = showtop (theNode);
              while (tmp->type == 'm' && tmp->next->back != NULL)
                {
                  tmp = tmp->next->back;
                }
              if (tmp != NULL && tmp->type == 'i')
                fprintf (stdout, "-%li-", showtop (theNode)->id);
              else
                fprintf (stdout, "-%li*-", showtop (theNode)->id);
            }
        }
      else
        {
          if (theNode->type == 't')
            fprintf (stdout, "-{%li:%s}\n", theNode->id, theNode->nayme);
          else
            fprintf (stdout, "ROOT=%li-", theNode->id);
        }
      if (theNode->tyme > 0)
        {
          if (theNode->next->back != NULL)
            dtree (theNode->next->back);
          if (theNode->type != 'm' && theNode->next->next->back != NULL)
            dtree (theNode->next->next->back);
        }
    }
}

node *
showtop (node * theNode)
{
  if (theNode == NULL)
    return NULL;
  else
    {
      if (theNode->top)
        {
          return theNode;
        }
      else
        {
          if (theNode->next->top)
            {
              return theNode->next;
            }
          else
            {
              return theNode->next->next;
            }
        }
    }

}

void
jumble (long *s, long n)
{
  long *temp, i, rr, tn = n;

  temp = (long *) calloc (1, sizeof (long) * n);
  memcpy (temp, s, sizeof (long) * n);
  for (i = 0; i < n && tn > 0; i++)
    {
      s[i] = temp[rr = RANDINT (0, tn - 1)];
      temp[rr] = temp[tn - 1];
      tn--;
    }
  free (temp);
}

