 #include <stdio.h>
 #include <math.h>
 #include <gsl/gsl_multifit.h>
 #include "fit.h"
 //quadratically extrapolates a window of width 'width' and gives
 //a)an array 'weights' of width 'width', which is the log of the next extrapolated window
 //b)the coefficiants of the quadratic that is the next extrapolated window, 'coeffs'
double fit(int width, int l, double ** m_bin, double * weights, FILE *fw, double * coeffs)
 {
   int i, count = 0;
   double chisq, min_weight = 0, max_weight = 0, med_weight;
   gsl_matrix *cov, *X;
   X = gsl_matrix_alloc (width+1, 3);
   cov = gsl_matrix_alloc (3, 3);
   
   gsl_vector *c, *y;
   c = gsl_vector_alloc (3);
   y = gsl_vector_alloc(width+1);

   //start fitting; first, weights[] is the log of the first window, which we then put into
   //gsl_multifit to get the coefficients for a quadratic fit:
		for(i = 0;i < width+1;i ++){
			weights[i] = (m_bin[l][i]>=1)?log(m_bin[l][i]):(-1);
			}
		if(weights[0] == -1)
			weights[0] = weights[1];
		gsl_vector_set (y, 0, weights[0]);
		for(i = 1;i < width+1;i ++){
			if(weights[i] == -1){
				if(i != width){
					count ++;
					weights[i] = (weights[i+1]+weights[i-1])/2;		//we want a smooth curve
					}
				else
					weights[i] = weights[i-1];
				}
			gsl_vector_set (y, i, weights[i]);
			}
		
	//printf("\n\ncount of weights = -1 is %d\n\n", count);
   
   for(i = 0; i < width+1; i ++){
   gsl_matrix_set (X, i, 0, 1.0);
   gsl_matrix_set (X, i, 1, i);
   gsl_matrix_set (X, i, 2, i*i);
   }
 
   {
     gsl_multifit_linear_workspace * work 
       = gsl_multifit_linear_alloc (width+1, 3);
     gsl_multifit_linear (X, y, c, cov,
                           &chisq, work);
     gsl_multifit_linear_free (work);
   }
 
 #define C(i) (gsl_vector_get(c,(i)))
 #define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
 /*
   {
     printf ("# best fit: Y = %g + %g X + %g X^2\n", 
             C(0), C(1), C(2));
 
     printf ("# covariance matrix:\n");
     printf ("[ %+.5e, %+.5e, %+.5e  \n",
                COV(0,0), COV(0,1), COV(0,2));
     printf ("  %+.5e, %+.5e, %+.5e  \n", 
                COV(1,0), COV(1,1), COV(1,2));
     printf ("  %+.5e, %+.5e, %+.5e ]\n", 
                COV(2,0), COV(2,1), COV(2,2));
     printf ("# chisq = %g\n", chisq);
   }
 */
for(i=0;i<3;i++)
coeffs[i] = C(i);

for(i = 0;i < 2*width+1;i ++){
	//printf("%lf\t", (i<width+1)?weights[i]:(0));
	fprintf(fw,"%d\t%lf\t%lf\n", l*width + i, (i<width+1)?weights[i]:(0), coeffs[0]+i*coeffs[1]+i*i*coeffs[2]);
	}
	fprintf(fw,"\n\n");
min_weight = max_weight = coeffs[0]+(width+1)*coeffs[1]+(width+1)*(width+1)*coeffs[2];
//extrapolate from the fit: second time around, weights[] are the extrapolated ln(histogram[m]) for the next bin
for(i = width+1;i < 2*width+2;i ++){
	weights[i-width-1] = coeffs[0]+i*coeffs[1]+i*i*coeffs[2];
	if(min_weight > weights[i-width-1])
		min_weight = weights[i-width-1];
	if(max_weight < weights[i-width-1])
		max_weight = weights[i-width-1];
	}
//to keep the range of weights[] reasonable:
med_weight = (min_weight+max_weight)/2.0;
for(i = 0;i < width+1;i ++)
	weights[i] -= med_weight;

//end fitting
   gsl_matrix_free (cov);
   gsl_matrix_free (X);
   gsl_vector_free (y);
   gsl_vector_free (c);
   return med_weight;
 }

