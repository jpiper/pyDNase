// Some of the below are courtesy of the Cephes library
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

double gammln(double xx){
    //Returns the value ln[Γ(xx)] for xx > 0.
    //Internal arithmetic will be done in double precision, a nicety that you can omit if ﬁve-ﬁgure accuracy is good enough.
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677, 24.01409824083091,-1.231739572450155, 0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

double betacf(int a, int b, double x)
{
    //Used by betai: Evaluates continued fraction for incomplete beta function by modiﬁed Lentz’smethod (§5.2).
    int m,m2;
    double aa,c,d,del,h,qab,qam,qap;
    qab=a+b; //These q’s will be used in factors that occur
    qap=a+1.0; //in the coeﬃcients (6.4.6).
    qam=a-1.0;
    c=1.0; //First step of Lentz’s method.
    d=1.0-qab*x/qap;
    if (fabs(d) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;
    for (m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d; //One step (the even one) of the recurrence.
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d;// Next step of the recurrence (the odd one).
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;// Are we done?
    }
    return h;
}


double betai(int a, int b, double x)
{

    //Returns the incomplete beta function I
    //x(a,b).
    double bt;
    if (x < 0.0 || x > 1.0)
        {
            printf("Bad x in routine betai");
        }
    if (x == 0.0 || x == 1.0)
    {
        bt=0.0;
    }

    else //Factors in front of the continued fraction.

    bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));

    if (x < (a+1.0)/(a+b+2.0)) //Use continued fraction directly.
    {
    return bt*betacf(a,b,x)/a;
    }
    else //Use continued fraction after making the symmetry transformation.
        return 1.0-bt*betacf(b,a,1.0-x)/(double)b;

}

double bdtrc( int k, int n, double p )
{
double dk, dn;

if( (p < 0.0) || (p > 1.0) )
	goto domerr;
if( k < 0 )
	return( 1.0 );

if( n < k )
	{
domerr:
	return( -0);
	}

if( k == n )
	return( 0.0 );
dn = n - k;
if( k == 0 )
	{
	if( p < .01 )
		dk = -expm1( dn * log1p(-p) );
	else
		dk = 1.0 - pow( 1.0-p, dn );
	}
else
	{
	dk = k + 1;
	dk = betai( dk, dn, p );
	}
return( log(dk) );
}

double bdtr( int k, int n, double p )
{
double dk;

if( (p < 0.0) || (p > 1.0) )
	goto domerr;
if( k < 0 )
	return( 1.0 );

if( n < k )
	{
domerr:
	return( -0);
	}

if( k == n )
	return( 0.0 );
	dk = k + 1;
	dk = betai( n-k, k+1, 1-p );
return( log(dk) );
}

unsigned int * rolling_window(unsigned int * ist, unsigned int w_size, unsigned int size)
{
	//unsigned int * totals  = malloc(sizeof(unsigned int) * size + 1 - w_size);;
	unsigned int * totals = malloc(sizeof(unsigned int) * (size + 1 - w_size));
	unsigned int rollingadd = 0;
	for (unsigned int i = 0; i <w_size; i++)
	{
		rollingadd += ist[i];
	}
	totals[0] = rollingadd;
	for (unsigned int i = w_size; i < size; i++) {
		rollingadd -= ist[i-w_size];
		rollingadd += ist[i];
		totals[i-w_size+1] = rollingadd;
	}
	return totals;
}

struct tuple2{
	float * fpscores;
	unsigned int * mles;
} ;

struct tuple2 * wellington(unsigned int * f,  unsigned int * r, unsigned int length, unsigned int * shoulder_sizes, unsigned int shoulders, unsigned int * fp_sizes, unsigned int fps)
{
	float * scores = calloc(length, sizeof(float));
	unsigned int * mle = calloc(length, sizeof(unsigned int));
	for (unsigned int i = 0; i < shoulders; i++)
	{
		unsigned const int shoulder = shoulder_sizes[i];
		unsigned int * const f_bindingArray = rolling_window(f,shoulder, length);
		unsigned int * const b_bindingArray = rolling_window(r,shoulder, length);

		for(unsigned int j = 0 ;j < fps ;j++)
		{

			unsigned const int fp_size = fp_sizes[j];
			unsigned int * const fw_fpscores = rolling_window(f,fp_size, length);
			unsigned int * const rv_fpscores = rolling_window(r,fp_size, length);
			unsigned const int halffpround = floor((fp_size-1)/2);

			for (unsigned int i = halffpround + shoulder ; i < length-shoulder-halffpround; i++)
			{

				unsigned const int xForward  = f_bindingArray[i-halffpround-1-shoulder+1];
				unsigned const int nForward  = xForward + fw_fpscores[i-halffpround];
				unsigned const int xBackward = b_bindingArray[i+halffpround+1];
				unsigned const int nBackward = xBackward + rv_fpscores[i-halffpround];

				if (xForward > 0 & xBackward > 0)
				{
					float const p = (float)shoulder / (shoulder + fp_size);
					float const score = bdtrc(xForward - 1, nForward, p) + bdtrc(xBackward - 1, nBackward, p);

					if(score < scores[i])
					{
						scores[i] = score;
						mle[i] = fp_size;
					}
				}
			}
			free(fw_fpscores);
			free(rv_fpscores);
		}
		free(f_bindingArray);
		free(b_bindingArray);
	}
	struct tuple2 * const retarr = malloc(sizeof(struct tuple2));
	retarr->fpscores = scores;
	retarr->mles = mle;
	return retarr;
}
