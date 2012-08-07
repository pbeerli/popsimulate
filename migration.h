
#ifndef MIGRATION_HEADER	/*migration.h */
#define MIGRATION_HEADER

#include "definitions.h"
/*-----------------------------------------------------------------
  Maximum likelihood estimation of migration rates 
  using coalescent trees

  Peter Beerli
  Genetics 357360
  University of Washington
  Seattle, WA 98195-7360, USA
  beerli@genetics.washington.edu
 
  With help of Joe Felsenstein (joe@agenetics.washington.edu),
  Mary Kuhner and Jon Yamato (mkkuhner@genetics.washington.edu)

  *----------------------------------------------------------------
 
  Time, 'tyme', in this tree is measured from the tips to the root.
  I.e. the tips are at tyme '0', and the root node has the largest
  value for 'tyme'. 
  

  */

typedef char allele_type[DEFAULT_ALLELENMLENGTH];

typedef double sitelike[5];
typedef sitelike *ratelike;
typedef ratelike *phenotype;

typedef double contribarr[MAXCATEGS];
typedef short val[MAXCATEGS];

typedef char allele_fmt[DEFAULT_ALLELENMLENGTH];

typedef struct lr_data_fmt
  {
    short type;
    long elem;
    char *value;
  }
lr_data_fmt;

typedef struct lratio_fmt
  {
    long alloccounter;
    long counter;
    lr_data_fmt *data;
  }
lratio_fmt;

typedef struct valrec
  {
    double rat, ratxi, ratxv, zz, z1, y1, ww1, zz1, ww2, zz2, z1zz, z1yy,
      xiz1, xiy1xv, ww1zz1, vv1zz1, ww2zz2, vv2zz2;
  }
valrec;

typedef valrec ***tbl_fmt;

typedef union xarray_fmt
  {
    double *a;
    phenotype s;
  }
xarray_fmt;

typedef struct seqmodel_fmt
  {
    double freqa, freqt, freqg, freqc;
    double freqr, freqy;
    double freqar, freqcy, freqgr, freqty;
    double aa, bb;
    long endsite;
    double xi, xv, ttratio, fracchange;
    long *sites;
    tbl_fmt tbl;
    long *alias;
    long *ally;
    long *category;
    short *weight;
    long weightsum;
    long *aliasweight;
    long *location;
    double **term;
    double **slopeterm;
    double **curveterm;
    val *mp;
    contribarr *contribution;
    contribarr like;
    contribarr nulike;
    contribarr clai;


  }
seqmodel_fmt;

/* defines the data structure read from infile */
typedef struct _data
  {
    FILE *infile;
    FILE *utreefile;
    FILE *weightfile;
    FILE *catfile;
    FILE *sumfile;
    char *****yy;
    allele_fmt **allele;
    long *maxalleles;
    boolean *skiploci;
    char **popnames;
    long **numind;
    char ***indnames;
    long numpop;
    long loci;
    seqmodel_fmt *seq;
    double freq;
    double freqlast;
    char dlm;
  }
data_fmt;


typedef struct _option
  {
    FILE *parmfile;
    FILE *seedfile;
    /*general options */
    long nmlength;		/* length of individual names */
    long popnmlength;		/* length of population names */
    long allelenmlength;	/* length of allele names */
    /*input/output options */
    int menu;
    boolean progress;
    boolean verbose;
    boolean movingsteps;
    double acceptfreq;
    boolean printdata;
    boolean usertree;
    short treeprint;
    boolean plot;
    boolean plotnow;
    short plotmethod;
    boolean simulation;
    char infilename[256];
    char outfilename[256];
    char mathfilename[256];
    char treefilename[256];
    char utreefilename[256];
    char catfilename[256];
    char weightfilename[256];
    char sumfilename[256];
    char title[81];
    lratio_fmt *lratio;
    long fsttype;
    boolean printfst;
    short profile;
    char profilemethod;
    boolean qdprofile;
    boolean printprofsummary;
    boolean printprofile;
    /* data options */
    char datatype;
    short migration_model;
    char dlm;
    long micro_stepnum;
    long micro_threshold;
    double **steps;
    boolean interleaved;
    double *ttratio;
    boolean freqsfrom;
    long categs;
    double *rate;
    long rcategs;
    double *rrate;
    double *probcat;
    double probsum;
    boolean autocorr;
    double lambda;
    boolean weights;
    double freqa;
    double freqc;
    double freqg;
    double freqt;
    /* random number options */
    short autoseed;
    long inseed;
    long saveseed;
    /* mcmc options */
    short thetaguess;
    short migrguess;
    boolean gamma;
    short burn_in;
    double *thetag;
    long numthetag;
    double *mg;
    long nummg;
    long schains;
    long sincrement;
    long ssteps;
    long lchains;
    long lincrement;
    long lsteps;
    double lcepsilon;
    long numpop;
    /* save genealogy summary options*/
    boolean readsum;
    boolean writesum;
  }
option_fmt;


/* used in the tree structure */
typedef struct _node
  {
    struct _node *next, *back;
    boolean tip;
    char type;
    long number;
    long pop;
    long actualpop;
    long id;
    xarray_fmt x;
    double *s;
    double lxmax;
    char *nayme;
    boolean top;
    boolean dirty;
    double v, tyme, length, xcoord;
    short ycoord, ymin, ymax;
	long species;
	long aspecies;
  }
node;

typedef struct vtlist
  {
    node *eventnode;		/* node with age=tyme */
    double age;			/* tyme from top nodelet */
    double interval;		/* interval t[i+1].age - t[i].age */
    long *lineages;
    long from;
    long to;
    /*  long pop; */
    long slice;
  }
vtlist;

typedef struct tree_fmt
  {
    node **nodep;
    node *root;
    long pop;
    long tips;
  }
tree_fmt;

typedef struct timelist_fmt
  {
    long numpop;
    long copies;
    long allocT;
    long T;
    long oldT;
    vtlist *tl;
  }
timelist_fmt;

typedef struct tarchive_fmt
  {
    long copies;
    double *km;
    double *kt;
    long *p;
    long *l;
  }
tarchive_fmt;

typedef struct timearchive_fmt
  {
    long allocT;
    long T;
    long numpop;
    long sumtips;
    double param_like;
    double thb;
    double alpha;
    tarchive_fmt *tl;
    double *param;
    double *param0;
    double *likelihood;
    long trials;
    double normd;
  }
timearchive_fmt;

typedef struct _plotmax
  {
    double x1;
    double y1;
    double l1;
    double x2;
    double y2;
    double l2;
  }
plotmax_fmt;

typedef struct _quantile_fmt
  {
    char *name;
    double *param;
  }
quantile_fmt;

typedef struct _world
  {
    /* generalities */
    option_fmt *options;
    data_fmt *data;
    long loci;
    long skipped;		/*loci with no data */
    long locus;			/* the current locus */
    long numpop;
    long numpop2;
    long sumtips;
    /* time archives, contains the data/results for summarizing */
    timearchive_fmt *atl;
    /*tree material */
    node **nodep;
    node *root;
    long unique_id;

    /* parameter */
    double *param0;
    double *param00;
    double **fstparam;

    /* mcmc related */
    long *lineages;
    timelist_fmt *treetimes;

    /* reporting */
    double *likelihood;
    plotmax_fmt **plotmax;
    quantile_fmt *quantiles;
    double maxdatallike;	/* the maximum log likleihood of a chain */
    double allikemax;		/* the maximum log likelihood  the best tree */
    boolean in_last_chain;	/*for last CHAIN option in print-tree */
    double param_like;
    double ***cov;
    long migration_counts;
    char ****plane;
    FILE *outfile;
    FILE *treefile;
    FILE *mathfile;
  }
world_fmt;




typedef struct _nr_fmt
  {
    long partsize;		/*number of part-variables, fixed per model */
    double *parts;		/* parts of the first and second derivatives */
    double *d;			/* first derivates */

    long numg;
    double *dv;
    double *gdelta;
    double *od;
    double *delta;
    tarchive_fmt *tl;
    world_fmt *world;
    double **dd;		/* second derivatives */
    double *param;		/* changed values of param */
    double *oparam;		/* saved old parameters */
    double llike;		/* parameter LOGlikelihood */
    double lastllike;		/* parameter LOGlikelihood */
    double ollike;		/* old "   " */
    double *datalike;		/*P(D|G) */
    double *locilikes;
    double **apg0;		/* part-loglikelihoods of param0 */
    double *apg;		/* part-loglikelihoods of param */
    double apg_max;		/* maxvalue of apg */
    double *gamma;
    long gammaI;
    double alpha;
    long numpop;		/* number of populations */
    long numpop2;		/* 2*numpop */
    double PGC;			/* uncorrected llike */
    double oPGC;		/* uncorrected ollike */
    long copy_nr;		/* number of genealogies */
    boolean *skiploci;
    long profilenum;            /* number of profile parameters*/
    long *profiles;             /* which profile parameters*/
    long *indeks;             /* which noprofile parameters*/
  }
nr_fmt;


typedef struct helper_fmt
  {
    long locus;
    nr_fmt *nr;
    timearchive_fmt *atl;
    double *param;
    long which;
    double weight;
    double ll;
  }
helper_fmt;



typedef struct _migr_table_fmt
  {
    long from;
    long to;
    double time;
  }
migr_table_fmt;

typedef struct proposal_fmt
  {
    world_fmt *world;
    boolean mig_removed;
    char datatype;
    short migration_model;
    long sumtips;
    long numpop;
    long endsite;
    double fracchange;
    double *param0;
    node *root;
    node *origin;
    node *target;
    node *realtarget;
    node *tsister;
    node *realtsister;
    node *osister;
    node *realosister;
    node *ocousin;
    node *realocousin;
    node *oback;
    node *realoback;
    node **line_f;
    node **line_t;
    node *connect;
    double likelihood;
    double time;
    double v;
    double vs;
    xarray_fmt xt;
    xarray_fmt xf;
    double mf;
    double mt;
    node **aboveorigin;
    node **bordernodes;
    long listsize;
    migr_table_fmt *migr_table;
    migr_table_fmt *migr_table2;
    long migr_table_counter;
    long migr_table_counter2;
    long old_migr_table_counter;
    long old_migr_table_counter2;
    long timeslice;
  }
proposal_fmt;

/* here a some definitions which should be remove in a future version
   (a) FLOC : place holder for first locus strcuture should finally deal
   with different number of individuals per locus */
#define FLOC 0


/* global variables (this should be not used at all) */
extern char appl[8];


#include "sighandler.h"
#include "tools.h"
#include "sort.h"

#endif
