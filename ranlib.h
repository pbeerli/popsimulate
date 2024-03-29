/* Prototypes for all user accessible RANLIB routines */
#ifndef _MYRANLIB_H__
extern void advnst(long k);
extern float genbet(float aa,float bb);
extern float genchi(float df);
extern float genexp(float av);
extern float genf(float dfn, float dfd);
extern float gengam(float a,float r);
extern void genmn(float *parm,float *x,float *work);
extern void genmul(long n,float *p,long ncat,long *ix);
extern float gennch(float df,float xnonc);
extern float gennf(float dfn, float dfd, float xnonc);
extern float gennor(float av,float sd);
extern void genprm(long *iarray,int larray);
extern float genunf(float low,float high);
extern void getsd(long *iseed1,long *iseed2);
extern void gscgn(long getset,long *g);
extern long ignbin(long n,float pp);
extern long ignnbn(long n,float p);
extern long ignlgi(void);
extern long ignpoi(float mu);
extern long ignuin(long low,long high);
extern void initgn(long isdtyp);
extern long mltmod(long a,long s,long m);
extern void phrtsd(char* phrase,long* seed1,long* seed2);
extern float ranf(void);
extern void setall(long iseed1,long iseed2);
extern void setant(long qvalue);
extern void setgmn(float *meanv,float *covm,long p,float *parm);
extern void setsd(long iseed1,long iseed2);
extern float sexpo(void);
extern float sgamma(float a);
extern float snorm(void);
#endif /*MYRANLIB*/

