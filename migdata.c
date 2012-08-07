/* reads a migration-tree and creates data acoording to the type
   (A, S, M).

   Peter Beerli
   beerli@fsu.edu
   migdata version 2.0
   September 2010

the program can take two arguments 
migdata first_pop max_snips

*/
#include <stdio.h>
#include <strings.h>
#include "migdata.h"

#undef MAXALLELES
#define MAXALLELES 500
#define STARTALLELE 100
#define RANDINT(a,b) (long)((a)+(randum(seed)* (((b)-(a))+1.)))
#define MAXLOCI 20*10000
FILE *plotfile;
char type, is_gamma;
long ntips, nextnode;
boolean haslengths, uselengths;
int has_msatlength;
char ch;
double mu;
longer seed;
long inseed, maxnch;
long loci;
long sites[10000];
long bug;
node *root;
FILE *outfile;
FILE *ancestor;
FILE *datefile;
short ****data;
long popnum;
long *maxind;
long *indnum;
long *indnumc;
char ***names;
double **ages;
double maxtyme;
long mc, cc;
double aa, bb, ttratio;
long mutated;
long *oldsite, *newsite, **polymorph;
node **whereptr;
boolean fullout;
double **siterates;
long *popseq;
long maxsnips_per_locus;
long *flanking;
long *replength;

short newstate (short oldstate, double len, short start);
void evolve (node * p, short oldstate, long locus, long *ind, long site);
void printdata (void);
double randum (long *seed);
void getseed (long seed);
void inputdata (void);
double lengthof (node * p);
node *crawlback (node * theNode);
void dtree (node * theNode);
node *showtop (node * theNode);
void pbprocesslength (node * p);
char processlength (FILE *file, node **p);
int eoln (FILE * f);
void probtrav (node * r);
void treeread (FILE *file, node **p, node *q);
char processbracket(FILE *file, node **p);
node * add_migration(node *p,long to,long from,double utime);
node *allocate_nodelet(long num, char type);
node * create_interior_node(node **q);
node * create_root_node(node **q);
node * create_tip_node(FILE *file, node **q, char *ch);
void treeout (node * joint, node * p, long s);
void length_to_times(node **p);
void findnames(node *r, long **indcount, double mu);
void evolve_allele(node *p, short oldstate, long locus, long *ind);
void evolve_seq (node * p, short oldstate, long locus, long *ind, long site);
void printdata_seq (void);
void printdata_allele(void);
void printdata_micro(void);
long ran_poisson(float mean);
long ran_normal(float mean);

void printdata_snp (void);

void translate (char *text, char from, char to);
	      
void
evolve_infsites (node * p, short oldstate, long locus, 
long *ind, long site);

void print_pairwise_fst(void);	      
void print_pairwise_ia_fst(void);	      
void
fill_tree (node * p, node *q, short state, short newstate, long locus, long *ind, long site);

char swiss[][40]={"Romanshorn","Arbon","Kreuzlingen","Frauenfeld","Guendelhart","Homburg","Aarau","L'Abbaye","Aigle","Alpnach","Appenzell","Apples","Arni","Aubonne","Bad","Baden","Ballens","Basel","Chessel","Chexbres","Chur","Corbeyrier","Crissier","Cunter","Dallenwil","Davos","Ddingen","Egg","Emmetten","Engelberg","Ennetmoos","Epalinges","Erlenbach","Esslingen","Etoy","Fchy","Feldbach","Fribourg","Gais","Geneva","Gimel","Giswil","Gonten","Grub","Gryon","Habsburg","Heiden","Hergiswil","Herisau","Hundwil","Interlaken","Kerns","Lausanne","Lauterbrunnen","Lavey-Morcles","Leysin","Locarno","Longirod","Lucerne","Lugano","Lungern","Lutzenberg","Meilen","Montreux","Noville","Oberdorf","Oberegg","Ollon","Ormont-Dessous","Ormont-Dessus","Paudex","Prilly","Pully","Rapperswil","Rehetobel","Reichenau","Renens","Rennaz","Reute","Roche","Rolle","Sachseln","Sarnen","Schaffhausen","Schlatt-Haslen","Schleinikon","Schnengrund","Schwellbrunn","Schwyz","Schwende","Sion","Solothurn","Speicher","St.Gallen","Stans","Stansstad","Stein","Teufen","Thal","Thalwil","Thun","Trogen","Urnaesch","Uzwil","Vevey","Villeneuve","Wald","Waldstatt","Walzenhausen","Weiach","Wildhaus","Winterthur","Wolfenschiessen","Wolfhalden","Yvorne","Zug","Zurich"};



double watterson_a(long n)
{
  long k;
  double sum=0.;
  for(k=1; k<n; k++)
    {
      sum += 1. / k;
    }
  return sum;
}
double watterson_b(long n)
{
  long k;
  double sum=0.;
  for(k=1; k<n; k++)
    {
      sum += 1. / (k*k);
    }
  return sum;
}

double watterson(long segreg, long n)
{
  return segreg / watterson_a(n);
}

//wattersonvar[s_, n_] :=
//((an = (Sum[1./k, {k, 1,
// n}]))(thetaw = watterson[s, n]) +  (bn = (Sum[
// 1./(k^2), {k, 1, n}]))thetaw^2 )/(an^2)
double wattersonvar(double  thetaw, long n)
{
  double an = watterson_a(n);
  double bn = watterson_b(n);
  return (an * thetaw + bn * thetaw * thetaw) / (an*an);
}

double wattersonstd(double  thetaw, long n)
{
  return sqrt(wattersonvar(thetaw,n));
}


double  factorial(double n)
{
  long nn=(long) n;
  switch(nn){
  case  0: return 1.;
  case  1: return 1.;
  case  2: return 2.;
  case  3: return 6.;
  case  4: return 24.;
  case  5: return 120.;
  case  6: return 720.;
  case  7: return 5040.;
  case  8: return 40320.;
  case  9: return 362880.;
  case 10: return 3628800.;
  case 11: return 39916800.;
  case 12: return 479001600.;
  default:
    return exp(lgamma(n+1.));
  }
}


void
openfile (FILE ** fp, char *filename, char *mode, char *application, char *perm)
{
  FILE *of;
  char file[100];
  strcpy (file, filename);
  while (1) {
    of = fopen (file, mode);
    if (of)
      break;
    else {
      switch (*mode) {
      case 'r':
	printf ("%s:  can't read %s\n", application, file);
	file[0] = '\0';
	while (file[0] == '\0') {
	  printf ("Please enter a new filename>");
	  fgets (file,99, stdin);
	}
	break;
      case 'w':
	printf ("%s: can't write %s\n", application, file);
	file[0] = '\0';
	while (file[0] == '\0') {
	  printf ("Please enter a new filename>");
          fgets (file,99, stdin);
	}
	break;
      }
    }
  }
  *fp = of;
  if (perm != NULL)
    strcpy (perm, file);
}

void
uppercase (ch)
     char *ch;
{				/* make ch upper-case */
  *ch = (islower (*ch) ? toupper (*ch) : (*ch));
}				/* uppercase */

void
getch (c)
     char *c;
{				/* get next nonblank character */
  do {
    *c = getc (stdin);
  }
  while ((*c == ' ') || (*c == '\n') || (*c == '\t'));
}				/* getch */

char
processlength (FILE *file, node **p)
{
  char ch;
  long digit, ordzero;
  double valyew, divisor;
  boolean pointread, minusread;

  ordzero = '0';
  pointread = FALSE;
  minusread = FALSE;
  valyew = 0.0;
  divisor = 1.0;
  ch = getc(file);
  digit = ch - ordzero;
  /*
  while (((unsigned long) digit <= 9) | (ch == '.') || (ch == '-')) {
    if (ch == '.')
      pointread = TRUE;
    else if (ch == '-')
      minusread = TRUE;
    else {
      valyew = valyew * 10.0 + digit;
      if (pointread)
	divisor *= 10.0;
    }
    ch = getc (file);
    digit = ch - ordzero;
  }
  if (!minusread)
    (*p)->length = valyew / divisor;
  else
    (*p)->length = 0.0;
  printf("%f\n",(*p)->length);
  return ch;
  */
  fscanf(file, "%lf%c",&(*p)->length,&ch);
  printf("%f\n",(*p)->length);
  return ch; 
}

void
treeread(FILE *file, node **pp, node *q)
{
  node *p;
  char ch=getc(file);
  while(ch!=';'){
    switch(ch){
    case '(': p = create_interior_node(&q);
      q = p->next;
      ch = getc(file);
      break;
    case ',': q = q->next;
      if (q->top){
	fprintf(stderr,"Multifurcation handling not yet installed");
	exit(-1);
      }
      ch = getc(file);
      break;
    case ')': p = showtop(q);
      q = p->back;
      ch = getc(file);
      break;
    case ' ':
    case '\n':
    case '\t': ch= getc(file);
      break;
    case ':': ch = processlength(file,&p);
      break;
    case '[': ch = processbracket(file, &p);
      q->back = p;
      p->back = q;
      break;
    default: p = create_tip_node(file, &q, &ch);
      break;
    }
  }
  (*pp) = p;
//  fscanf(file, "%*[^\n]");
  getc(file);
}

void time_from_root(node **p);
void time_from_tip(node **p);
int
main (int argc, char *argv[])
{
  long pop, oldstate=0,oldmutated=0;
  char line[20000];
  int locus, site;
  long polycount;
  long startpop=0;
  fullout=FALSE;
  loci = 1;
  maxnch = 10;
  openfile (&outfile, "infile", "w", argv[0], NULL);
  openfile (&datefile, "datefile", "w", argv[0], NULL);
  inputdata ();
  maxsnips_per_locus=LONG_MAX;
  if(argc>1)
  {
      startpop = atoi(argv[1]);
      if(argc>2)
	maxsnips_per_locus=atol(argv[2]);
  }
  popseq=(long *) calloc(1,sizeof(long) * popnum);
  for(pop=0;pop<popnum;pop++)
  {
      popseq[pop] = startpop++;
      if(startpop==popnum)
          startpop=0;
      fprintf(stderr,"%li ", popseq[pop]);
  }
  getseed (inseed);
  polymorph = (long **) calloc(1,sizeof(long*) * loci);
  if(type=='I')
    {
      oldsite = (long *) calloc(1,sizeof(long) * 2);
      newsite = (long *) calloc(1,sizeof(long) * 2);
      whereptr = (node **) calloc(1,sizeof(node*) * 2);
    }	  

  if(type=='I')
    {
      openfile (&ancestor, "ancestor", "w", argv[0], NULL);
//      polymorph = (long **) realloc(polymorph,sizeof(long *) * loci);
    }
  for (locus = 0; locus < loci; locus++) {
    maxtyme=0;
    memset(indnumc,0,sizeof(long)*popnum);
    fgets (line,MAXLOCI,stdin);// read locus tree header line
    if(strchr(line,'$')!=NULL)
      {
	fgets (line,MAXLOCI,stdin);// read locus tree mutation rate
	sscanf(line,"#$%lf",&mu);
      }
      else
	mu=1.0;
    treeread (stdin,&root,NULL);
    root = showtop(root);
    time_from_root(&root);
    time_from_tip(&root);
    //length_to_times(&root);
    // I infinite sites?
    // S sequence
    // N linked snps
    // L panel known (is first population?
    // P panel unknown 
    if(type=='I' || type=='S' || type=='N' || type=='L'){
            //  probtrav (root);
      findnames(root,&indnumc, mu);
      root->length= 10000.;
      root->back->next->next->tyme = root->tyme + 10000.;
      if(type=='I')
	{
	  mutated=0; oldmutated=0;
	  oldsite = (long *) realloc(oldsite,sizeof(long) * 2*sites[locus]);
	  polymorph[locus] = (long *) calloc(1,
				       sizeof(long) * 2*sites[locus]);
	  memset(oldsite,0,sizeof(long)*2*sites[locus]);
	  memset(polymorph[locus],0,sizeof(long)*2*sites[locus]);
	  newsite = (long *) realloc(newsite,sizeof(long) * 2*sites[locus]);
	  memset(oldsite,0,sizeof(long)*2*sites[locus]);
	  whereptr = (node **) realloc(whereptr,sizeof(node*) * 2*sites[locus]);
	  memset(whereptr,0,sizeof(node*)*2*sites[locus]);
	}	  
      else
      {
	mutated=sites[locus];
        polymorph[locus] = (long *) calloc(1,
                                           sizeof(long) * 2*sites[locus]);
        memset(polymorph[locus],0,sizeof(long)*2*sites[locus]);
      }
      
      for (site = 0; site < sites[locus] || mutated < sites[locus]; site++) {
	for (pop = 0; pop < popnum; pop++) {
	  maxind[pop] = 0;
	}
	oldstate = RANDINT(0,3);
	oldmutated=mutated;
	evolve (root, oldstate, locus, maxind, site);
	if(type=='I' && mutated==oldmutated)
	  {
	    whereptr[mutated] = root;
	    oldsite[mutated] = oldstate;
	    newsite[mutated] = oldstate;
	    mutated++;
	  }
      }
      if(type=='I')
	{
	  polycount=0;
	  for(site=0;site < sites[locus] || site < mutated; site++) 
	    {
	      for (pop = 0; pop < popnum; pop++) {
		maxind[pop] = 0;
	      }
	      fill_tree(root, whereptr[site],oldsite[site],newsite[site],locus, maxind,site);
	      if(newsite[site]!=oldsite[site]){
		polymorph[locus][site] = 1;	    
		polycount++;
		fprintf(ancestor,"%c ",oldsite[site]==0 ? 'A' : \
			(oldsite[site]==1 ? 'G' : (oldsite[site]==2 ? 'C' : 'T')));
		fprintf(stderr,"%i ",site);
	      } 
	    }
		fprintf(stderr,"\nTotal polymorphic sites: %li\n",polycount);
	  fprintf(ancestor,"\n");
	  fprintf(stdout,"\n");
	}
    }
    else {
      findnames(root,&indnumc, mu);

      root->length=10000.;
      root->back->next->next->tyme = root->tyme + 10000.;
      if(type=='A')
	evolve (root, newstate('A',0,0), locus, maxind,1);
      else
	evolve (root, STARTALLELE, locus, maxind,1);
    }
  }
  printdata ();
  fclose (outfile);
  fclose (datefile);

  if(type=='I')
    {
      fclose(ancestor);
    }
  return 0;
}

void
inputdata (void)
{
  long locus, ii, i;
  long pop, ind;
  char input[1024], *temp;
  fgets (input,1024,stdin);
  sscanf (input, "#%c%c", &type, &is_gamma);
  if(type=='J')
    {
      fullout=TRUE;
      type='I';
    }
  if(type=='P') //panel 
    {
      fullout=TRUE;
      type='L';
    }
  fgets (input,1023,stdin);
  sscanf (input, "#%li", &inseed);
    fgets (input,1023,stdin);
  sscanf (input, "#%li", &popnum);
    fgets (input,1023,stdin);
  indnum = calloc (1, sizeof (long) * popnum);
  indnumc = calloc (1, sizeof (long) * popnum);
  maxind = calloc (1, sizeof (long) * popnum);
  if (popnum == 2)
    sscanf (input, "#%li%li", &indnum[0], &indnum[1]);
  else {
    temp = calloc (1, sizeof (char) * 100);
    temp = strtok (input, " ");
    indnum[0] = atoi (temp + 1);
    for (pop = 1; pop < popnum; pop++) {
      temp = strtok (NULL, " ");
      if (temp == NULL)
	break;
      indnum[pop] = atoi (temp);
    }
  }
  fgets (input,1023,stdin);
  if(type=='I' || type=='S' || type == 'N' || type=='L' || type == 'P'){
    sscanf (input, "#%li%li%lf", &loci, &sites[0], &ttratio);
    siterates = (double**) calloc(1,sizeof(double*)*loci);
    for(ii=0;ii<loci;ii++)
    {
        siterates[ii] = (double*) calloc(1,sizeof(double)*sites[0]);
        fgets(input,1024,stdin); // reads # rate ... for locus x
	while ( fgetc(stdin) != '#' );
        for(i=0;i<sites[0];i++)
        {
            fscanf(stdin,"%lf ",&siterates[ii][i]);
        }
    }
    
    for (locus = 1; locus < loci; locus++) {
      sites[locus] = sites[0];
    }
    aa = 2.0 / (1.0 + ttratio);
    bb = (2.0 * ttratio - 1.0) / (ttratio + 1.0);
  }
  else {
    sscanf(input,"#%li %i",&loci,&has_msatlength);
    flanking = (long *) calloc(loci,sizeof(long));
    replength = (long *) calloc(loci,sizeof(long));
    for (pop = 0; pop < popnum; pop++) 
      indnum[pop] *= 2;
    if(has_msatlength==1)
      {
	char dummy;
	fscanf(stdin,"%c",&dummy);
	for (locus = 0; locus < loci; locus++) {
	  sites[locus] = 1;
	  fscanf(stdin,"%li%li", &flanking[locus], &replength[locus]);
	}
	fgets(input,LINESIZE,stdin);
      }
    else 
      {
	for (locus = 0; locus < loci; locus++) {
	  sites[locus] = 1;
	  flanking[locus]=0;
	  replength[locus]=1.0;
	}
      }
  }

  data = (short****) calloc (1, sizeof (short ***) * popnum);
  names= (char***)calloc(1,sizeof(char **)* popnum);
  ages= (double**)calloc(1,sizeof(double *)* popnum);
  for (pop = 0; pop < popnum; pop++) {
    names[pop] = (char **) calloc(1,sizeof(char**) *  indnum[pop]);
    ages[pop] = (double *) calloc(1,sizeof(double) *  indnum[pop]);
    data[pop] = calloc (1, sizeof (short **) * indnum[pop]);

    for (ind = 0; ind < indnum[pop]; ind++) {
      names[pop][ind] = (char *) calloc(1,sizeof(char)*11);
      data[pop][ind] = calloc (1, sizeof (short *) * loci);
      for (locus = 0; locus < loci; locus++) {
	data[pop][ind][locus] = calloc (1, sizeof (short) * 2*sites[locus]);
      }
    }
  }
}


long  newstate_brownian(long oldstate, double len, long  start)
{
  long temp;
  if(start==0)
    return oldstate;
  else {
    temp = ran_normal(len);
    //printf("%li ", temp);
    return oldstate + temp;
  }
}

long ran_normal(float mean)
{
  const double tpi = 6.2831853071795864769;
  long value=0;
  double u = -2.0 * log(randum(seed));
  double v = tpi * randum(seed);
  value = lround(mean * sqrt(u) * cos(v));
  return value;
}


long  newstate_micro(long oldstate, double len, long  start)
{
  long temp;
  if(start==0)
    return oldstate;
  else {
    temp = ran_poisson(len);
    //printf("%li ", temp);
    return oldstate + temp;
  }
}


long ran_poisson(float mean)
{
   long value = 0;
   double prob, test = randum(seed);
   prob = exp(-mean);
   while (test > prob) {
	  test *= randum(seed);
	  value += ((randum(seed) < 0.5) ? -1. : 1.);
   }
   return value;
}



long  newstate_micro_old(long oldstate, double len, long  start )
{
  static double *pis, **steps;
  long i;
  long k;
  long state;
  double temp,kk, pisum=0,sum;
  if(start==0){
    pis = (double*) realloc(pis,sizeof(double)*MAXALLELES);
    steps = (double**) calloc(1,sizeof(double*)*MAXALLELES);
    steps[0] = (double*) calloc(1,sizeof(double)*MAXALLELES*MAXALLELES);
    for(i=1;i<MAXALLELES;i++) {
      steps[i] = steps[0] + MAXALLELES*i;
    }
    for(i=1;i<MAXALLELES;i++) {
      for (k=i;k<MAXALLELES && k<i+20;k += 2){
	steps[i-1][k-i] += pow(0.5,k) / (factorial((k-i+1)/2)* factorial((k+i-1)/2));
      }
    }
  }
  if (start==0) {
    return oldstate;
  }
  memset(pis,0,sizeof(double)*MAXALLELES);
  kk = randum(seed);
  for(i=1;i<MAXALLELES;i++) {
    temp=0.0;
    for (k=0;k<20 && i+k<MAXALLELES;k += 2){
      temp += pow(len,k) * steps[i-1][k];
    }
    pis[i-1] = temp * exp(-len);
    pisum += pis[i-1];
  }
  sum=0;
  if(pisum!=0.0){
    for(i=1;i<MAXALLELES;i++) {
      sum+=pis[i-1]/pisum;
      if(kk<sum)
	break;
    }
  }
  else
    i=1;
  state = oldstate + (randum(seed) > 0.5 ? -1 : 1) * (i-1);
  if (state<1) {
    fprintf(stderr,"state was smaller than 1\n");
    state=1;
  }
  if (state>200) {
    fprintf(stderr,"state was bigger than 200\n");
    state=200;
  }
  return state;
}



short newstate(short oldstate, double len, short  start )
{
    static short alleles[93];/* for chars from ! to ~*/
    static short num;
    static boolean wrapped=FALSE;
    long i;
    if (type=='M') {
      num=newstate_micro(oldstate,len,start);
      return num;
    }
    if (type=='B') {
      num=newstate_brownian(oldstate,len,start);
      return num;
    }
    if (start==0) {
	for(i=0;i<93;i++)
	    alleles[i]='\0';
	alleles[0]=oldstate;
	num=1;
	return oldstate;
    }
    else {
	alleles[num]=alleles[num-1]+1;
	if(alleles[num]>176){
	  alleles[num]=41;
	  wrapped=TRUE;
	}
	if(wrapped){
	  if(alleles[num]>='A')
	    fprintf(stderr,"Wrapped around all possible allele states!\n");
	}  
	if (alleles[num]=='?') {
	    num++;
	    alleles[num]=alleles[num-1]+1;
	    num++;
	    return alleles[num-1];
	}
	else {
	    num++;
	    return alleles[num-1];
	}
    }
}

void
evolve(node * p, short oldstate, long locus, long *ind, long site)
{
    memset(ind,0,sizeof(long)*popnum);
    switch(type){
        case 'L':
        case 'N':
        case 'P':
        case 'S': evolve_seq (p,oldstate,locus,ind,site);
            break;
        case 'A': evolve_allele(p,oldstate,locus,ind);
            break;
        case 'B':
        case 'M': evolve_allele(p,oldstate,locus,ind);
            break;
        case 'I':
            evolve_infsites(p,oldstate,locus,ind,site);
            break;
    }
}

void evolve_allele(node *p, short oldstate, long locus, long *ind)
{
    char state = oldstate;
    double rr;
    double rrr;
    if(type=='A'){
      rr= randum(seed);
      rrr=1-exp(-(showtop(crawlback(p))->tyme-p->tyme));
      if (rr<rrr)
	state = (char) newstate(oldstate,0,1);
    }
    else {
      if(crawlback(p)->type!='r')
	state = newstate(oldstate,showtop(crawlback(p))->tyme - p->tyme, 1);
      else
	state = newstate(oldstate,showtop(crawlback(p))->tyme - p->tyme, 0);
    }
    if (p->tip) {
	data[atoi(p->nayme)][ind[atoi(p->nayme)]][locus][0]=state;
	ind[atoi(p->nayme)] += 1;
    }
    else {
	evolve_allele(crawlback(p->next),state, locus, ind);
	evolve_allele(crawlback(p->next->next),state, locus, ind);
    }
}

	      
void
evolve_seq (node * p, short oldstate, long locus, long *ind, long site)
{
  long state = oldstate;
  double rr = randum (seed);
  double x;
  
  if (p != root)
    x = (showtop(crawlback(p))->tyme - p->tyme)*siterates[locus][site];
  else
    x = 0.0;

  p->probxi = 0.25 + 0.25 * exp (-aa * x) - 0.5 * exp (-(aa + bb) * x);
  p->probxv = 0.5 * (1.0 - exp (-aa * x));
  if (rr < p->probxi) {
    if (oldstate < 2)
      state = 1 - oldstate;
    else
      state = 5 - oldstate;
  }
  else if (rr < p->probxi + p->probxv) {
    if (oldstate < 2)
      state = (long) (2 * randum (seed)) + 2;
    else
      state = (long) (2 * randum (seed));
  }
  else
    state = oldstate;
  if (p->tip) {
    data[atoi (p->nayme)][ind[atoi (p->nayme)]][locus][site] = (short) state;
    ind[atoi (p->nayme)] += 1;
  }
  else {
    evolve_seq (crawlback(p->next), state, locus, ind, site);
    evolve_seq (crawlback(p->next->next), state, locus, ind, site);
  }
}

	      
void
fill_tree (node * p, node *q, short state, short newstate, long locus, long *ind, long site)
{
  if(q==p)
    {
      state=newstate;
    }
  if (p->tip) {
    data[atoi (p->nayme)][ind[atoi (p->nayme)]][locus][site] = (short) state;
    ind[atoi (p->nayme)] += 1;
    //    fprintf(stdout,"%10.10s %li %li\n",p->nayme,ind[atoi(p->nayme)],site);   
  }
  else {
    fill_tree (crawlback(p->next), q, state, newstate,locus, ind, site);
    fill_tree (crawlback(p->next->next), q, state, newstate, locus, ind, site);
  }
}

	      
void
evolve_infsites (node * p, short oldstate, long locus, long *ind, long site)
{
  long state = oldstate;
  double rr = randum (seed);
  if (rr < p->probxi) {
    if (oldstate < 2)
      state = 1 - oldstate;
    else
      state = 5 - oldstate;
    whereptr[mutated] = p;
    oldsite[mutated] = oldstate;
    newsite[mutated] = state;
    mutated++;
  }
  else 
    {
      if (rr < p->probxi + p->probxv) {
	if (oldstate < 2)
	  state = (long) (2 * randum (seed)) + 2;
	else
	  state = (long) (2 * randum (seed));
	whereptr[mutated] = p;
	oldsite[mutated] = oldstate;
	newsite[mutated] = state;
	mutated++;
      }
      else
	state = oldstate;
    }
  if (p->tip) {
    //    data[atoi (p->nayme)][ind[atoi (p->nayme)]][locus][site] = (short) state;
    // ind[atoi (p->nayme)] += 1;
  }
  else {
    evolve_infsites (crawlback(p->next), state, locus, ind, site);
    evolve_infsites (crawlback(p->next->next), state, locus, ind, site);
  }
}


void
printdata(void)
{
    switch(type){
        case 'I':
        case 'S': printdata_seq ();
            break;
        case 'L':
        case 'P':
        case 'N': printdata_snp ();
            break;
        case 'A': printdata_allele();
            break;
        case 'B':
        case 'M': printdata_micro();
            break;
    }
}


void
printdata_seq (void)
{
  double thetaw;
  double hw, hb;
  long countw, countb;
  long ind1,ind2, pop2;
  long i, pop, locus, ind, site, sum, indsum=0;
  char ch=' ';
  fprintf (outfile, "   %li %li\n", popnum, loci);
  fprintf (datefile, "   %li %li\n", popnum, loci);
  for (locus = 0; locus < loci; locus++) {
    fprintf (outfile, "%li ", sites[locus]);
  }
  fprintf (outfile, "\n");
//  for (pop = 0; pop < popnum; pop++) {
  for (i = 0; i < popnum; i++)
  {
      pop = popseq[i];
      indsum += maxind[pop];
      fprintf (outfile, "%-3li  %s%3li\n", maxind[pop], swiss[pop], pop);
      fprintf (datefile, "%-3li  %s%3li\n", maxind[pop], swiss[pop], pop);
      for (locus = 0; locus < loci; locus++) {
          for (ind = 0; ind < maxind[pop]; ind++) {
              fprintf (outfile, "%-10.10s", names[pop][ind]);
              fprintf (datefile, "%-10.10s", names[pop][ind]);
              for (site = 0; site < sites[locus]; site++) {
                  if(type=='I' && polymorph[locus][site]==0 && !fullout)
                  {
                      continue;
                  }
                  if(type=='S')
                  {
                      polymorph[locus][site] += data[pop][ind][locus][site]!=
                          data[0][0][locus][site];
                  }
                  switch (data[pop][ind][locus][site]) {
                      case 0:
                          ch = 'A';
                          break;
                      case 1:
                          ch = 'G';
                          break;
                      case 2:
                          ch = 'C';
                          break;
                      case 3:
                          ch = 'T';
                          break;
                      case 4:
                          ch = 'T';
                          break;
		  }
                  fprintf (outfile, "%c", ch);
              }
#ifdef NEWFORMAT
              fprintf (outfile, "@");
#else
              fprintf (outfile, "\n");
#endif
              fprintf (datefile, "%f\n", ages[pop][ind]);
          }
	  fprintf(outfile,"#\n");
      }
  }
  if(type=='S')
  {    
      for (locus = 0; locus < loci; locus++) {
          sum =0;
          for (site = 0; site < sites[locus]; site++) {
              sum += polymorph[locus][site]>0;
          }
          fprintf(stderr,"\nLocus %li> %li variable sites out of %li (Theta_Watterson = %f, std=%f)\n",
                  locus+1,sum,sites[locus], (thetaw=watterson(sum, indsum))/sites[locus], wattersonstd(thetaw,indsum)/sites[locus]);
      }
      /*
      // print FST as 1 - Hw/Hb
      hw = hb = 0.0;
      countw = countb = 0;
      for (locus = 0; locus < loci; locus++) 
	{
          for (site = 0; site < sites[locus]; site++) 
	    {
	      if(polymorph[locus][site]>0)
		for(pop=0;pop<popnum;pop++)
		  {
		    //fb calculation
		    for (ind1 = 0; ind1 < maxind[pop]; ind1++) 
		      {
			for(pop2=0;pop2<pop;pop2++)
			  {
			    for (ind2 = 0; ind2 < maxind[pop2]; ind2++) 
			      {
				if(data[pop][ind1][locus][site]!=data[pop2][ind2][locus][site])
				  {
				    hb += 2.0;
				    countb += 2;
				  }
				else
				  {
				    countb +=2;
				  }
			      }
			  }
			
			//fw calculation
			for (ind2 = 0; ind2 < ind1; ind2++) 
			  {
			    if(data[pop][ind1][locus][site]!=data[pop][ind2][locus][site])
			      {
				hw += 2.0;
				countw += 2;
			      }
			    else
			      {
				countw +=2;
			      }
			  }
		      }
		  }
	    }
	}
      hw /= countw;
      hb /= countb;
      fprintf(stderr,"FST = 1 - Hw/Hb = 1 - %f/%f = %f\n",hw,hb,1.-hw/hb);
      */
      print_pairwise_fst();
  }
  //  else
  //  {
  //    print_pairwise_ia_fst();
  //  }
}

#define MYREAL double
void print_pairwise_fst(void)
{
  // print FST as 1 - Hw/Hb
  double hw, hb, mm;
  long site;
  long countw, countb;
  long locus, pop, pop2, ind1, ind2;
  MYREAL *hwp;
  MYREAL **hbp;
  long **countbp;
  long *countwp;
  hwp = (MYREAL *) calloc(popnum,sizeof(MYREAL));
  countwp = (long *) calloc(popnum,sizeof(long));
  hbp = (MYREAL**) calloc(popnum, sizeof(MYREAL *));
  countbp = (long **) calloc(popnum,sizeof(long *));
  for(pop=0; pop < popnum; pop++)
    {
      hbp[pop] = (MYREAL*) calloc(popnum, sizeof(MYREAL));
      countbp[pop] = (long*) calloc(popnum, sizeof(long));
    }
  
  hw = hb = 0.0;
  countw = countb = 0;
  for (locus = 0; locus < loci; locus++) 
    {
      for (site = 0; site < sites[locus]; site++) 
	{
	  if(polymorph[locus][site]>0)
	    for(pop=0;pop<popnum;pop++)
	      {
		//fb calculation
		for (ind1 = 0; ind1 < maxind[pop]; ind1++) 
		  {
		    for(pop2=0;pop2<pop;pop2++)
		      {
			for (ind2 = 0; ind2 < maxind[pop2]; ind2++) 
			  {
			    if(data[pop][ind1][locus][site]!=data[pop2][ind2][locus][site])
			      {
				hbp[pop][pop2] += 1.0;
				hbp[pop2][pop] += 1.0;
				countbp[pop][pop2] += 1;
				countbp[pop2][pop] += 1;
				hb += 2.0;
				countb += 2;
			      }
			    else
			      {
				countb +=2;
				countbp[pop][pop2] += 1;
				countbp[pop2][pop] += 1;
			      }
			  }
		      }
		    
		    //fw calculation
		    for (ind2 = 0; ind2 < ind1; ind2++) 
		      {
			if(data[pop][ind1][locus][site]!=data[pop][ind2][locus][site])
			  {
			    hwp[pop] += 1.0;
			    countwp[pop] += 1;
			    hw += 2.0;
			    countw += 2;
			  }
			else
			  {
			    countwp[pop] += 1;
			    countw +=2;
			  }
		      }
		  }
	      }
	}
    }
  hw /= countw;
  hb /= countb;
  fprintf(stderr,"FST = 1 - Hw/Hb = 1 - %f/%f = %f\n",hw,hb,1.-hw/hb);
  for(pop=0;pop < popnum; pop++)
    {
      fprintf(stderr,"# ");
      for(pop2=0; pop2 < popnum; pop2++)
	{
	  if(pop==pop2)
	    {
	      fprintf(stderr,"*%g ",hwp[pop]/countwp[pop]);
	    }
	  else
	    {
	      mm = hwp[pop]/countwp[pop] + hwp[pop2]/countwp[pop2];
	      mm /= 2.;
	      fprintf(stderr,"%g ",1.- mm /(hbp[pop][pop2] / countbp[pop][pop2]));
	    }
	}
      fprintf(stderr,"\n");
    }
}


void print_pairwise_ia_fst(void)
{
  // print FST as Ht - Hw / Ht
  double hw, hb, mm;
  long site;
  long countw, countb;
  long locus, pop, pop2, ind1, ind2;
  MYREAL *hwp;
  MYREAL **hbp;
  long **countbp;
  long *countwp;
  hwp = (MYREAL *) calloc(popnum,sizeof(MYREAL));
  countwp = (long *) calloc(popnum,sizeof(long));
  hbp = (MYREAL**) calloc(popnum, sizeof(MYREAL *));
  countbp = (long **) calloc(popnum,sizeof(long *));
  for(pop=0; pop < popnum; pop++)
    {
      hbp[pop] = (MYREAL*) calloc(popnum, sizeof(MYREAL));
      countbp[pop] = (long*) calloc(popnum, sizeof(long));
    }
  
  hw = hb = 0.0;
  countw = countb = 0;
  for (locus = 0; locus < loci; locus++) 
    {
      for (site = 0; site < sites[locus]; site++) 
	{
	  if(polymorph[locus][site]>0)
	    for(pop=0;pop<popnum;pop++)
	      {
		//fb calculation
		for (ind1 = 0; ind1 < maxind[pop]; ind1++) 
		  {
		    for(pop2=0;pop2<pop;pop2++)
		      {
			for (ind2 = 0; ind2 < maxind[pop2]; ind2++) 
			  {
			    if(data[pop][ind1][locus][site]!=data[pop2][ind2][locus][site])
			      {
				hbp[pop][pop2] += 1.0;
				hbp[pop2][pop] += 1.0;
				countbp[pop][pop2] += 1;
				countbp[pop2][pop] += 1;
				hb += 2.0;
				countb += 2;
			      }
			    else
			      {
				countb +=2;
				countbp[pop][pop2] += 1;
				countbp[pop2][pop] += 1;
			      }
			  }
		      }
		    
		    //fw calculation
		    for (ind2 = 0; ind2 < ind1; ind2++) 
		      {
			if(data[pop][ind1][locus][site]!=data[pop][ind2][locus][site])
			  {
			    hwp[pop] += 1.0;
			    countwp[pop] += 1;
			    hw += 2.0;
			    countw += 2;
			  }
			else
			  {
			    countwp[pop] += 1;
			    countw +=2;
			  }
		      }
		  }
	      }
	}
    }
  hw /= countw;
  hb /= countb;
  fprintf(stderr,"FST = 1 - Hw/Hb = 1 - %f/%f = %f\n",hw,hb,1.-hw/hb);
  for(pop=0;pop < popnum; pop++)
    {
      fprintf(stderr,"# ");
      for(pop2=0; pop2 < popnum; pop2++)
	{
	  if(pop==pop2)
	    {
	      fprintf(stderr,"*%g ",hwp[pop]/countwp[pop]);
	    }
	  else
	    {
	      mm = hwp[pop]/countwp[pop] + hwp[pop2]/countwp[pop2];
	      mm /= 2.;
	      fprintf(stderr,"%g ",1.- mm /(hbp[pop][pop2] / countbp[pop][pop2]));
	    }
	}
      fprintf(stderr,"\n");
    }
}

void
printdata_snp (void)
{
  long snips = 0;
  long maxsnips = maxsnips_per_locus;
  long pop, test, locus, ind, site;
  char ch=' ';
  long *numb;
  short **weight;
  weight = (short **) calloc(1,sizeof(short*)*loci);
  numb = (long *) calloc(1,sizeof(long)*loci);
  fprintf (outfile, "   %li %li\n", popnum, loci);
  for (locus = 0; locus < loci; locus++) {
    weight[locus] = calloc(1,sizeof(short)*sites[locus]);
    for (site = 0; site < sites[locus]; site++) {
      for (pop = 0; pop < (type!='L' ? popnum : 1); pop++) {
	for (ind = 0; ind < maxind[pop]; ind++) {
	  if(ind==0 && pop==0)
	    test=data[pop][ind][locus][site];
	  else {
	    {
                if(type!='L')
                {
                    if(data[pop][ind][locus][site] != test)
                    {
                        weight[locus][site]=1.;
                        numb[locus] += 1;
                        goto sitedone;
                    }
                }
                else
                {
                    if(data[0][ind][locus][site] != test)
                    {
                        weight[locus][site]=1.;
                        numb[locus] += 1;
                        goto sitedone;
                    }                    
                }
                
                
	    }
	  }
	}
      }
    sitedone:
      /* DO NOTHING HERE*/
      ;
    }
  }

  for (locus = 0; locus < loci; locus++) 
    {  
      snips=0;
      for (site = 0; site < sites[locus]; site++) 
	{
	  if(snips >= maxsnips)
	    weight[locus][site]=0;
	  else
	    {
	      if(weight[locus][site]>0)
		snips++;
	    }
	}
    }

  for (locus = 0; locus < loci; locus++) {
      fprintf(stderr,"Locus %li> Number of variable sites in %s = %li\n",
              locus+1,(type=='L' ? "panel population" : "whole data set"),
              numb[locus]);
      if(maxsnips_per_locus<numb[locus])
	fprintf (outfile, "%li ", maxsnips_per_locus);
      else
	fprintf (outfile, "%li ", numb[locus]);
      //      if(type!='L')
      //fprintf (outfile, "%li ", numb[locus]);
      //else
      //{
      //        //numb[locus]=1;
      //	fprintf (outfile, "%li ", numb[locus]);
      // }
  }
  fprintf (outfile, "\n");
  for (pop = 0; pop < popnum; pop++) {
    fprintf (outfile, "%-3li  %s%3li\n", maxind[pop], swiss[pop], pop);
    for (locus = 0; locus < loci; locus++) 
      {	
	for (ind = 0; ind < maxind[pop]; ind++) {
	  fprintf (outfile, "%-10.10s", names[pop][ind]);
	  for (site = 0; site < sites[locus]; site++) {
	    if(weight[locus][site]==0)
	      continue;

	    if(fullout && type=='L' && pop==0)
              ch = '?';
	    else
	      {
              switch (data[pop][ind][locus][site]) {
                  case 0:
                      ch = 'A';
                      break;
                  case 1:
                      ch = 'G';
                      break;
                  case 2:
                      ch = 'C';
                      break;
                  case 3:
                      ch = 'T';
                      break;
                  case 4:
                      ch = 'T';
                      break;
              }
          }
	  fprintf (outfile, "%c", ch);
              //if(type=='L')
              //break;
	}
	fprintf (outfile, "\n");
      }
    }
  }
}



void printdata_allele(void)
{
    /*for twopop*/
    long pop, locus, ind;
    char *tmp;
    char delimiter='/';
    tmp = (char *) calloc(1,sizeof(char)*100);
    fprintf(outfile,"   %li %li %c \n",popnum, loci, delimiter);
    for(pop=0;pop<popnum;pop++) {
		fprintf(outfile,"%-3li  %s%3li\n",maxind[pop]/2,swiss[pop], pop);			
	    for (ind=0;ind<maxind[pop]/2;ind++){
	      sprintf(tmp,"%s%s", 
		      names[pop][ind],names[pop][ind+(maxind[pop]/2)]);
	      fprintf (outfile, "%-10.10s", tmp);
		for(locus=0;locus<loci;locus++){
		    fprintf(outfile,"%c%c%c ",
			    data[pop][ind][locus][0],delimiter,
			    data[pop][ind+(maxind[pop]/2)][locus][0]);
		}
		fprintf(outfile,"\n");
	    }
    }
}

void printdata_micro(void)
{
    long pop, locus, ind;
    long mini;
    char * indname = (char *) calloc(1024,sizeof(char));
    fprintf(outfile,"   %li %li . \n",popnum,loci);
    if(has_msatlength==1)
      {
	fprintf(outfile,"#@Microsat-repeatlength ");
	for(locus=0;locus<loci;locus++)
	  {
	    fprintf(outfile,"%li ",replength[locus]);
	  }
	fprintf(outfile,"\n");
      }
  /* check  if no allele is below 0 and if so add minimum to all in this
     locus*/
    for(locus=0;locus<loci;locus++){
      mini = 1; 
      for(pop=0;pop<popnum;pop++) {
	for (ind=0;ind<maxind[pop];ind++){
	  if(data[pop][ind][locus][0]<mini)
	    mini = data[pop][ind][locus][0];
	  
	}
      }
      if(mini<=0)
	{
	  mini--;
	  for(pop=0;pop<popnum;pop++) {
	    for (ind=0;ind<maxind[pop];ind++){
	      data[pop][ind][locus][0] -= mini;
	    }
	  }
	}
    }
    /* printout corrected values*/
    for(pop=0;pop<popnum;pop++) {
		fprintf(outfile,"%-3li  %s%3li\n",maxind[pop]/2,swiss[pop],pop);			
	    for (ind=0;ind<maxind[pop]/2;ind++){
	      //	      fprintf (outfile, "%-5.5s%-5.5s", names[pop][ind],names[pop][ind+(maxind[pop]/2)]);
	      sprintf(indname, "%-5.5s%-5.5s", names[pop][ind],names[pop][ind+(maxind[pop]/2)]);
	      translate(indname,' ','_');
	      fprintf(outfile,"%-10.10s",indname);
	      for(locus=0;locus<loci;locus++){
		    fprintf(outfile,"%04li.%04li ",
			    flanking[locus]+replength[locus]*data[pop][ind][locus][0],
			    flanking[locus]+replength[locus]*data[pop][ind+(maxind[pop]/2)][locus][0]);
		}
		fprintf(outfile,"\n");
	    }
    }
    //    print_pairwise_ia_fst();
    free(indname);
}


int
eof (f)
     FILE *f;
{
  register long ch;

  if (feof (f))
    return 1;
  if (f == stdin)
    return 0;
  ch = getc (f);
  if (ch == EOF)
    return 1;
  ungetc (ch, f);
  return 0;
}


int
eoln (FILE * f)
{
  register long ch;

  ch = getc (f);
  if (ch == EOF)
    return 1;
  ungetc (ch, f);
  return (ch == '\n');
}

void
memerror ()
{
  printf ("Error allocating memory\n");
  exit (-1);
}



void
getseed (long inseed)
{
  long i;
  double clearseed;
  for (i = 1; i <= 1000; i++)	/* clear the random numbergenerator */
    clearseed = randum (seed);
  for (i = 0; i <= 5; i++)
    seed[i] = 0;
  i = 0;
  do {
    seed[i] = inseed & 63;
    inseed /= 64;
    i++;
  } while (inseed != 0);
}


double
randum (long *seed)
{				/* randum -- slow but machine independent */
  /* random number generator -- slow but machine independent */
  long i, j, k, sum;
  longer mult, newseed;
  double x;

  mult[0] = 13;
  mult[1] = 24;
  mult[2] = 22;
  mult[3] = 6;
  for (i = 0; i <= 5; i++)
    newseed[i] = 0;
  for (i = 0; i <= 5; i++) {
    sum = newseed[i];
    k = i;
    if (i > 3)
      k = 3;
    for (j = 0; j <= k; j++)
      sum += mult[j] * seed[i - j];
    newseed[i] = sum;
    for (j = i; j <= 4; j++) {
      newseed[j + 1] += newseed[j] / 64;
      newseed[j] &= 63;
    }
  }
  memcpy (seed, newseed, sizeof (longer));
  seed[5] &= 3;
  x = 0.0;
  for (i = 0; i <= 5; i++)
    x = x / 64.0 + seed[i];
  x /= 4.0;
  return x;
}				/* randum */

double
lengthof (node * p)
{
  if (p->type == 'm')
    fprintf (stderr, "a migration node was feed into lengthof");
  return fabs (p->tyme - showtop(crawlback (p))->tyme);
}				/* length */

node *
crawlback (node * theNode)
{
  node *tmp = theNode->back;

  while (tmp->type == 'm') {
    tmp = tmp->next->back;
  }
  return tmp;
}

void
dtree (node * theNode)
{
/* pass atr (*tree).root->back for normal results */
  node *tmp;
  if (bug) {
    if ((theNode->type == 'i' || theNode->type == 'm') && theNode->top) {
      if (theNode->type == 'm')
	fprintf (stdout, "-M%li-", theNode->id);
      else {
	tmp = showtop (theNode);
	while (tmp->type == 'm' && tmp->next->back != NULL) {
	  tmp = tmp->next->back;
	}
	if (tmp != NULL && tmp->type == 'i')
	  fprintf (stdout, "-%li-", showtop (theNode)->id);
	else
	  fprintf (stdout, "-%li*-", showtop (theNode)->id);
      }
    }
    else {
      if (theNode->type == 't')
	fprintf (stdout, "-{%li:%s}\n", theNode->id, theNode->nayme);
      else
	fprintf (stdout, "ROOT=%li-", theNode->id);
    }
    if (theNode->type != 't') {
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
  else {
    if (theNode->top) {
      return theNode;
    }
    else {
      if (theNode->next->top) {
	return theNode->next;
      }
      else {
	return theNode->next->next;
      }
    }
  }

}

void
probtrav (node * r)
{
  /* traverses the tree, filling in rates of change of a branch
     in the node at the top of the branch */
  double x;

  if (r != root)
    x = (showtop(crawlback(r))->tyme - r->tyme);
  else
    x = 0.0;
//  r->probxi = 0.25 + 0.25 * exp (-aa * x) - 0.5 * exp (-(aa + bb) * x);
//  r->probxv = 0.5 * (1.0 - exp (-aa * x));
  if (!r->tip) {
    probtrav (crawlback(r->next));
    probtrav (crawlback(r->next->next));
  }
}				/* probtrav */


char processbracket(FILE *file, node **p)
{
    long pop1, pop2;
    double utime;
    char c;
    c = getc(file);
    if(c=='&'){
	c = getc(file);   	
	switch(c){
	case 'T':
	  fscanf(file,"%li %li:%lf",&pop1,&pop2,&utime);
	  c = getc(file);
	  break;
	case 'M':
	  fscanf(file,"%li %li:%lf",&pop1,&pop2,&utime);
	  c = getc(file);
	  (*p) = add_migration(*p,pop1,pop2,utime);
	  break;
	default:
	  while(c!=']')
	    c=getc(file);
	  break;
	}
    }
    else {
	while(c!=']')
	    c=getc(file);
    }
    c=getc(file);
    return c;
}


node * add_migration(node *p,long to,long from,double utime)
{
    node *tmp;
    tmp = allocate_nodelet(2,'m');
    tmp->top = TRUE;
    tmp->next->back = p;
    p->back=tmp->next;        
    tmp->length = p->length - utime;
    p->length = utime;    
    tmp->pop = tmp->next->pop = from;
    tmp->actualpop = tmp->next->actualpop = to;
    return tmp;
}	



node *allocate_nodelet(long num, char type)
{
  static long unique_id=0;
  boolean isfirst=TRUE;
  long j;
  node *p, *q = NULL, *pfirst=NULL;
  long temp;
  temp = unique_id;
  for (j = 0; j < num; j++) {
    p = (node *) malloc(sizeof(node));
    p->tip = FALSE;
    p->pop = 0;
    p->actualpop = 0;
    p->type = type;
    p->top = FALSE;
    p->dirty = TRUE;
    p->next = q;
    p->x.s = NULL;	
    p->x.a = NULL;
    p->lxmax=0.0;
    p->back = NULL;
    p->nayme = NULL;
    p->v = 0.0;
    p->tyme = 0.0;
    p->length = 0.0;
    p->number = temp;
    p->id = unique_id++;
    if(isfirst){
      isfirst=FALSE;
      pfirst = p;
    }
    q = p;
  }
  pfirst->next = q;
  return pfirst;
}

node * create_interior_node(node **q)
{
  node *p;
  p = allocate_nodelet(3,'i');
  p->top = TRUE;
  p->back = *q;
  if((*q)==NULL)
    create_root_node(&p);
  else
    (*q)->back = p;
  return p;
}

node * create_root_node(node **q)
{
  node *p;
  p = allocate_nodelet(3,'r');
  p->top = TRUE;
  p->next->back = *q;
  (*q)->back = p->next;
  return p;
}


node * create_tip_node(FILE *file, node **q, char *ch)
{
  node *p;
  char c;
  char *nayme;
  long nl;
  long i=1;
  nayme = (char *) calloc(1,sizeof(char) * 200);
  nayme[0] = (*ch);
  while(strchr("[):;,\t\n",c=getc(file))==NULL)
    nayme[i++]=c;
  nayme[i]='\0';
  p = allocate_nodelet(1,'t');
  nl = strlen(nayme);
  p->nayme=(char *) calloc(1,sizeof(char) *(nl+1));
  p->top = TRUE;
  p->tip = TRUE;
  strcpy(p->nayme,nayme);
  p->back = *q;
  (*q)->back = p;
  free(nayme);
  (*ch) = c;
  return p;
}

void
translate (char *text, char from, char to)
{
  int i, j, gap = 0;
  while (text[gap] == from)
    gap++;
  for (i = gap, j = 0; text[i] != '\0'; i++) {
    if (text[i] != from) {
      text[j++] = text[i];
    }
    else {
      if (text[i - 1] != from) {
	text[j++] = to;
      }
    }
  }
  text[j] = '\0';
}




void
treeout (node * joint, node * p, long s)
{
  /* write out file with representation of final tree */
  static long col=0;
  long w;
  double x;
  char migstring[30];
  if (p->type == 't') {
    translate (p->nayme, ' ', '_');
    fprintf (stdout, "%s", p->nayme);
    col += strlen (p->nayme);
  }
  else {
      putc ('(', stdout);
      col++;
      treeout (joint, crawlback (p->next), s);
      putc (',', stdout);
      col++;
      if (col > 80) {
	putc ('\n', stdout);
	col = 0;
      }
      treeout (joint, crawlback (p->next->next), s);
      putc (')', stdout);
      col++;
  }
  if (p == joint) {
    x = 0.0;
  }
  else {
    x = showtop(crawlback(p))->tyme - p->tyme;
  }
  if (x > 0.0) {
    w = (long) (0.4343 * log (x));
  }
  else {
	if (x == 0.0)
	    w = 0;
	else
	w = (long) (0.4343 * log (-x)) + 1;
  }
  if (w < 0)
    w = 0;
  fprintf (stdout, ":%*.10f", (int) (w + 7), x < 10000. ? x : 10000.);
  col += w + 8;
  if (col > 80) {
	    putc ('\n', stdout);
	    col = 0;
	    }
  if(p!=joint){
	p = showtop(p->back);
	while(p->type == 'm') {
	    sprintf (migstring, " [&M %li %li:%g]",
		    p->pop, p->actualpop,p->tyme - p->next->back->tyme );
	    fprintf(stdout,"%s",migstring);
	    col += strlen(migstring)+1;
	    if (col > 80) {
	    putc ('\n', stdout);
	    col = 0;
	    }
	    p = showtop(p->back);      
	}
    }
    else{
	fprintf (stdout,";\n");
	col=0;
    }
}				/* treeout */

void length_to_times(node **p)
{
  node *q;
  if(*p==NULL)
    return;
  if(!(*p)->tip){
    length_to_times(&(*p)->next->back);
    if((*p)->type=='i')
      length_to_times(&(*p)->next->next->back);
  }
  q=showtop((*p)->back);
  if(q!=NULL)
     q->tyme = q->next->tyme = q->next->next->tyme = (*p)->tyme + (*p)->length;
}

void time_from_root(node **p)
{
  node *q;
  if(*p==NULL)
    return;
  q = showtop(*p);
  if(q->dirty)
    {
      if(q->type!='r')
	q->tyme = q->length + showtop(q->back)->tyme;
      else
	q->tyme = 10000;
      if(maxtyme < q->tyme)
	maxtyme = q->tyme;
      q->dirty = FALSE;
      if(q->tip)
	printf("@@%s %lf\n",q->nayme,q->tyme);
    }
  if(!(*p)->tip){
    time_from_root(&(*p)->next->back);
    if((*p)->type=='i')
      time_from_root(&(*p)->next->next->back);
  }
}

void time_from_tip(node **p)
{
  node *q;
  if(*p==NULL)
    return;
  (*p)->tyme  =fabs((*p)->tyme - maxtyme);
  if(!(*p)->tip){
    time_from_tip(&(*p)->next->back);
    if((*p)->type=='i')
      time_from_tip(&(*p)->next->next->back);
  }
}

void findnames(node *r, long **indcount, double mu)
{
  long ii,jj;
  if(r==NULL)
    return;
  if(!r->tip){

    findnames(crawlback(r->next),indcount, mu);
    findnames(crawlback(r->next->next),indcount, mu);
  }
  else {
    ii = atoi(r->nayme);
    jj = (*indcount)[ii];
    strncpy(names[ii][jj],r->nayme,10);
    ages[ii][jj] = r->tyme / mu;
    (*indcount)[ii] += 1;
  }
}









