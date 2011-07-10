#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_multifit.h>
#include "func-umbrella-nematic.h"
#include "mt64.h"
#include "fit.h"
/* generates a random number on [0,1)-real-interval */
//double genrand64_real2(void);
#define r genrand64_real2()
#define twopi 6.2831853
/*this macro checks whether s is within the correct window, 
  then discretizes the variables s and e so as to bin them
  so that the respective histograms can be built.*/
#define increment(s, e, l)	\
	{	\
	int_s = (int)(s);	\
	int_e = (int)(e);	\
	if((s_l <= int_s)&&(s_r >= int_s)){	\
		s_bin[l][int_s-s_l] += 1;	\
			e_bin[l][int_e - min_e] += 1;\
			if(s_r > int_s){	\
				esbin[int_s][int_e - min_e] += 1;	\
				count ++;	\
				}			\
		}	\
	}	\

//gcc -g -lm mt19937-64.c fit.c func-umbrella-nematic.c -I. -lgsl -lgslcblas; gdb ./a.out
//r 8 1000 P_E.dat P_S.dat P_ES.dat weights_etc.dat 2000
//gcc -lm mt19937-64.c fit.c func-umbrella-nematic.c -I. -lgsl -lgslcblas; ./a.out 8 1000 P_E.dat P_S.dat P_ES.dat weights_etc.dat 2000 20 10; echo 'p "P_S.dat" u 1:3' | gnuplot -persist

double ener4(double H0, double spin, int i, double* a, int**c, double p);
double total_ener(double H0, int N2, double* a, int**c, double p);
int mstate(double* S00, double*S01, double* a, double* old_a, int N, int s_l, int s_r);
int main(int argc, char *argv[])
{

time_t t1,t2; (void) time(&t1); int N_BIN; double p; int n = 0;
p = atoi(argv[9]);
//int STEPS = 30000; int ISEQ = 3000;	//~1 min for N=8, n_bin = 8 => 256000 MCs/s
//int STEPS = 720000; int ISEQ = 72000;
FILE * fp4, * fp3, * fp5, * fp6, * fw;
fp3 = fopen(argv[3],"w+");			//P(E)
fp4 = fopen(argv[4],"w+"); 			//P(S)
//fp5 = fopen(argv[5],"w+");			//esbin
fp6 = fopen(argv[6],"w+");			//i, raw data - s_bin[][], weights[]
fw = fopen("weights.dat","w+");
int i, j, k, k_max, l,
s_l, s_r, int_s,
int_e, 
N, N2, s_size, 				//N^2; we have a square lattice of N^2 spins
width, 				//width of window
accept, noaccept, 	//acceptance ratios
steps, iseq, 		//total number of spin flips
min_e, e_size;
double epsilon, H0,
e, delta_e, 	
old_e, new_e, 
s, new_s, factor_s = 1, 			//the order parameter
min_expweights, max, 
factor_for_energy, 	//the factor we multiply the energy in window n+1 by
med_weight = 1,		//not used
theta, old_theta, new_theta, delta_theta, var_theta = 6.,
fit_s_bin_l_i, area = 0,
denom, acceptance_ratio,
time_in_seconds;
unsigned long int count;
N = atoi(argv[1]); epsilon = ((double)(atoi(argv[2]))/1000.);
N_BIN = N;
N2 = N*N;
time_in_seconds = atoi(argv[8]);
steps = (253186*time_in_seconds)/N_BIN;
iseq = steps/10;
s_size = (1*N2);
factor_s = s_size/(double)N2;	//need not be an integer
width = s_size/N_BIN;			//assuming that (s_size)%N_BIN = 0
H0 = ((double)(atoi(argv[7])))/1000.;
printf("H0 = %lf\n",H0);
e_size = (int)(2*N2+ceil((H0/2)*N2));
min_e = -e_size;
e = s = 0; 

init_genrand64((unsigned)(time(0)));

double S[2][2]; //traceless symmetric tensor, the largest eigenvalue of which is the nematic order parameter s;
				//won't be used, just for show
				/*
				it looks like
								
				              |2*(cos_i(\theta))^2 - 1			2*sin_i(\theta)*cos_i(\theta)|
				sum over i of |																 |
				              |2*cos_i(\theta)*sin_i(\theta)	2*(sin_i(\theta))^2 - 1		 |
				which is of the form
				              |a	 b|
				              |		  |
				              |b	-a|
				we define our order parameter to be the larger eigenvalue
				s = sqrt(a^2 + b^2)
				*/

double s00, s01, S00, S01, S00_temp, S01_temp;
				/*
				with reference to the martrix above, S00 is a and S01 is b
				instead of summing over the entire lattice at each spin flip to calculate s, we do something clever: 
				s00 and s01 are the change in a and b respectively:
				s00 = 2*pow(cos(new_theta),2) - 2*pow(cos(old_theta),2);
				s01 = 2*cos(new_theta)*sin(new_theta) - 2*cos(old_theta)*sin(old_theta);
				S00_temp = S00 + s00;
				S01_temp = S01 + s01;
				new_s = sqrt(S00_temp*S00_temp + S01_temp*S01_temp);	//largest eigenvalue
				*/

double* weights = (double*)calloc((width+1),sizeof(double));		//weight for each s - see Virnau&Muller
double* expweights = (double*)calloc((width+1),sizeof(double));		//exp of weights
double* coeffs = (double*)calloc((3),sizeof(double));				//coeffs to calculate weights; not used 

/*the weight factor for each bin; initially set to 1 for the first window l = 0, 
  from which it can be extrapolated through to the second and so on*/
for(i = 0;i < width+1;i ++)
	expweights[i] = 1;

/*the histogram for each window is built up separately.
  first the raw data is collected in s_bin[]:*/
double** s_bin;	
s_bin = (double**)calloc(N_BIN,sizeof(double*));
for(i=0;i<N_BIN;i++) {
	s_bin[i]=(double*)calloc((width+1),sizeof(double));
}

/*then it is unweighted using expweights[] and its logarithm is stored in:*/
double** fit_s_bin;	
fit_s_bin = (double**)calloc(N_BIN,sizeof(double*));
for(i=0;i<N_BIN;i++) {
	fit_s_bin[i]=(double*)calloc((width+1),sizeof(double));	//width+1 because we need an overlap for fitting purposes
}

/*after which the different windows are fit together using ratio_of_fit_s_bin[]
  to give the ln of P(s) plus some constant:*/
double* ln_s_bin;
ln_s_bin = (double*)calloc((s_size),sizeof(double));

//the scale factor for each window
double* ratio_of_fit_s_bin, max_ratio = 0;
ratio_of_fit_s_bin = (double*)calloc((N_BIN),sizeof(double));

/*2-D bin for the energy and order parameter together; initially the raw histogram is built up
  which is then unweighted using the ratio_of_fit_s_bin[] and expweights[] obtianed from the previous quantities.
  this is actually all we need; the next few arrays are actually superfluous*/
double** esbin;	
esbin = (double**)calloc(s_size,sizeof(double*));			//row indices for s
for(i=0;i<s_size;i++) {
	esbin[i]=(double*)calloc((e_size),sizeof(double));	//column indices for e
}

/*energy histogram; built up and unweighted similar to esbin[], 
  for comparison and as a test of correctness*/
double** e_bin;
e_bin = (double**)calloc(N_BIN,sizeof(double*));
for(i=0;i<N_BIN;i++) {
	e_bin[i]=(double*)calloc((e_size),sizeof(double));
}

//built up by summing over s for esbin[]
double* final_e_bin;
final_e_bin = (double*)calloc((e_size),sizeof(double));

//simply ln(final_e_bin[])
double* ln_e_bin, max_ln_s = -10000;
ln_e_bin = (double*)calloc((e_size),sizeof(double));

//simply exp(ln_s_bin[])
double* final_s_bin;
final_s_bin = (double*)calloc((s_size),sizeof(double));

if(!(ln_e_bin&&ln_s_bin&&final_e_bin&&final_s_bin)){
	fputs("Error allocating memory\n",stderr);
	exit(1);
 }

double* a;	//lattice
double* old_a;	//old_a of lattice
a=(double*)calloc(N2,sizeof(double));
old_a=(double*)calloc(N2,sizeof(double));
if(!(a&&old_a)) {
	fputs("Error allocating memory\n",stderr);
	exit(1);
 }

//neighbors[i][j] gives the jth neighbour of the spin at the ith position in the lattice
int** neighbors;
neighbors = (int**)calloc(N2,sizeof(int*));
for(i=0;i<N2;i++) {
	neighbors[i]=(int*)calloc((4),sizeof(int));
}

for(i = 0;i < N2;i ++){
	neighbors[i][0] = 	   ((i < N)?(i + N2 - N):(i - N));	//above
	neighbors[i][1] =  ((i >= N2-N)?(i - N2 + N):(i + N));	//below
	neighbors[i][2] =   ((i%N == 0)?(i -  1 + N):(i - 1));	//left
	neighbors[i][3] = ((i%N == N-1)?(i +  1 - N):(i + 1));	//right
	}

//initialized for calculating acceptance ratios
for(i=0;i<N2;i++){
	a[i] = r;
	}

printf("temperature is %lf\np is %lf\nN is %d\n", epsilon, p, N);

//START check for acceptance ratio
double min_accept, max_accept;
min_accept = 0.45;	//min accepted acceptance ratio
max_accept = 0.55;	//max accepted acceptance ratio
count = 0;
do{
	accept = noaccept = 0;
	for(k = 0;k < 1000;k ++){
		i = (int)(N2*r); 	//we use N because integer casting is in effect flooring the float
		//assert(i<N2);
		new_theta = fmod(a[i] + var_theta*(r - 0.5),twopi);
		old_e = ener4(H0,a[i],i,a,neighbors,p);		//old energy
		new_e = ener4(H0,new_theta,i,a,neighbors,p);	//new energy
		delta_e = new_e - old_e;
		if((delta_e < 0.0)||(exp(-delta_e*epsilon) > r)){	//Metropolis algorithm
			a[i] = new_theta;
			e += delta_e;			//energy
			if(k%N2 == 0)
				e = total_ener(H0,N2,a,neighbors,p);
			accept ++;
			continue;			//new config works
			}
		//new config doesn't work
		noaccept ++;
		count ++;
		}//end k loop
acceptance_ratio = ((double)accept)/((double)(accept + noaccept));
//printf("var_theta = %lf\tacceptance_ratio = %lf\n", var_theta, acceptance_ratio);
if(min_accept > acceptance_ratio)
	var_theta *= 0.99;
else if(acceptance_ratio > max_accept)
	if(var_theta < 6.)
		var_theta *= 1.01;
	else
		break;
} while ( (min_accept > acceptance_ratio) || (acceptance_ratio > max_accept) );
printf("var_theta = %lf\tacceptance_ratio = %lf\n", var_theta, acceptance_ratio);
accept = noaccept = 0;

//END check for acceptance ratio

//START loop for window l
for(l = 0;l < N_BIN;l ++){
	second_iteration:;
	count = S00 = S01 = 0;
	s_l = l*s_size/N_BIN;		//left  bound of s
	s_r = (l+1)*s_size/N_BIN;	//right bound of s
	s = factor_s*mstate(&S00, &S01, a, old_a, N, l*N2/N_BIN, (l+1)*N2/N_BIN);//choose a state with s_l < s < s_r
	e = total_ener(H0,N2,a,neighbors,p);
	
	//START loop to explore phase space in window l
	for(k = 0;k < steps;k ++){
		i = (int)((N2)*r);
		//assert(i<N2);
		old_theta = a[i];	//old spin
		delta_theta = var_theta*(r - 0.5);
		new_theta = fmod(old_theta + delta_theta,twopi);
		/*first we check if s is within the window,
		  then we put e into the metropolis algorithm.*/
		//update S[i][j]: 
		//s00 = 2*pow(cos(new_theta),2) - 2*pow(cos(old_theta),2);
		s00 = cos(2*new_theta) - cos(2*old_theta);
		//s01 = 2*cos(new_theta)*sin(new_theta) - 2*cos(old_theta)*sin(old_theta);
		s01 = sin(2*new_theta) - sin(2*old_theta);
		S00_temp = S00 + s00;
		S01_temp = S01 + s01;
		new_s = factor_s*sqrt(S00_temp*S00_temp + S01_temp*S01_temp);	//largest eigenvalue
		if((s_l > new_s)||(s_r < new_s-1)){//rejection step
			if(k > iseq){
				increment(s,e,l);
				}
			continue;
			}
		old_e = ener4(H0,a[i],i,a,neighbors,p);			//old energy
		new_e = ener4(H0,new_theta,i,a,neighbors,p);	//new energy
		delta_e = new_e - old_e;
		//Metropolis algorithm
		if((delta_e < 0.0)||(exp(-delta_e*epsilon + weights[(int)s-s_l] - weights[(int)new_s-s_l]) > r)){	
			a[i] = new_theta;
			S00 = S00_temp;
			S01 = S01_temp;
			e += delta_e;			//energy
			if(k%N2 == 0)
				e = total_ener(H0,N2,a,neighbors,p);	//because errors in e can add up
			s = new_s;
			if(k > iseq){
				accept ++;
				increment(s,e,l);
				}
				continue;			//new config works		
			}
		//new config doesn't work
		if(k > iseq){
			noaccept ++;
			increment(s,e,l);
			}
		}//end k loop
	
	//this is just used to normalise the graph stored for (visual) checking to see if all's well:
	min_expweights = 0;
	for(i = 0;i < width+1;i ++)
		min_expweights += expweights[i];
	
	for(i = 0;i < width+1;i ++){
		fprintf(fp6,"%d\t%lf\t%lf\t%lf\n",l*width+i,s_bin[l][i]/(double)k_max,expweights[i]/min_expweights,s_bin[l][i]/(double)k_max*expweights[i]/min_expweights);
		}fprintf(fp6,"\n\n");
	
	//for the zeroth window:
	if(l == 0){
		for(i = 0;i < width+1;i ++){
			//first we unweight the raw data
			s_bin[l][i] *= expweights[i];
			//if the value is non-zero, we can take the log, otherwise use a special value,
			//for example -1, to denote that there is no data there.
			fit_s_bin[l][i] = (s_bin[l][i]>0)?log(s_bin[l][i]):(-1);
			ln_s_bin[i] = fit_s_bin[l][i];
			}
		}
	else{
	//for every window after that, we fit the histogram to the one to its left
		for(i = 0;i < width+1;i ++)
			{
			//first we unweight the raw data
			s_bin[l][i] *= expweights[i];
			//if the value is non-zero, we can take the log, fit it, and go on
			if(s_bin[l][i] > 0)
				fit_s_bin_l_i = log(s_bin[l][i]);
			else
				{
				//if we're at the first histogram element of the window,
				//take the value to be equal to the one on it's left for smoothness
				//which happens to be the last element of the previous window
				if(i == 0)
					fit_s_bin_l_i  = fit_s_bin[l-1][width];
				else
				//if we're not at the first histogram element of the window,
				//take the value to be equal to the one on it's left for smoothness
				//which happens to be the last element of the same window
					fit_s_bin_l_i = fit_s_bin[l][i-1];
				}
			fit_s_bin[l][i] = fit_s_bin_l_i;
			}
		//we look at the difference between the first element of this window and the last element 
		//of the previous window and scale every element of the current window accordingly;
		//since we take the log first, instead of dividing the two quantities, we subtract.
		ratio_of_fit_s_bin[l] = fit_s_bin[l][0] - fit_s_bin[l-1][width];
		if(isnan(ratio_of_fit_s_bin[l]))
			printf("\nERROR: left = %lf\tright = %lf\tratio_of_fit_s_bin[%d] = left - right = %lf\n", fit_s_bin[l][0], fit_s_bin[l-1][width], 						l,ratio_of_fit_s_bin[l]);
		for(i = 0;i < width;i ++){
			factor_for_energy = exp(-ratio_of_fit_s_bin[l] + weights[i]);
			for(j = 0;j < e_size;j ++)
				esbin[l*width + i][j] *= factor_for_energy;
			}
		/*we now extrapolate the current window to get a rough estimate of the next window
		  which we can use to weight the transition probabilities and 'flatten' the
		  resulting histogram obtained during the simulation, so as to get better accuracy
		  on the side of the histogram that would otherwise have had less entries;
		  we input the width of your window 'width', the window index l,
		  the pointer to arrays s_bin and weights, the file pointer fw, and 
		  the pointer to the array of cofficients 'coeffs'. The only output of 'fit' 
		  that we're actually using in the code is the array 'weights' --
		  see 'fit.c' for further documentation. 
		*/
		med_weight = fit(width,l,s_bin,weights,fw,coeffs);
		//the windows have to overlap to be fit together, therefore width+1:
		for(i = 0;i < width+1;i ++){
			if(fit_s_bin[l][i] != -1)
				fit_s_bin[l][i] -= ratio_of_fit_s_bin[l];
			//but ln_s_bin[] is only of size s_size:
			if(l*width+i < s_size){
				ln_s_bin[l*width+i] = fit_s_bin[l][i];
				if(max_ln_s < ln_s_bin[l*width+i])
					max_ln_s = ln_s_bin[l*width+i];
				}
			expweights[i] = exp(weights[i]);
			}
		}
	}//end l loop

for(j = 0;j < e_size;j ++){
	for(i = 0;i < s_size;i ++){
		final_e_bin[j] += esbin[i][j];
		//fprintf(fp5,"%d\t%d\t%lf\n",i,j,esbin[i][j]);
		}
	ln_e_bin[j] = (final_e_bin[j]>0)?log(final_e_bin[j]):(-1);
	fprintf(fp3,"%d\t%lf\t%lf\n", j, ln_e_bin[j], final_e_bin[j]);
	}

area = 0;
for(i = 0;i < s_size;i ++){
	final_s_bin[i] = exp(ln_s_bin[i]);
	area += final_s_bin[i]/s_size;
	}

for(i = 0;i < s_size;i ++){
	fprintf(fp4,"%d\t%lf\t%lf\n",i,ln_s_bin[i],final_s_bin[i]/area);
	}
/*
fclose(fp4);
fclose(fw);
free(ratio);
free(s_bin);
free(weights);
free(expweights);
free(fit_s_bin);
*/
printf("acceptance ratio with var_theta = %lf is acceptance_ratio = %f, noaccept_ratio = %f\n", var_theta, ((float)accept)/((float)(accept + noaccept)), ((float)noaccept)/((float)(accept + noaccept)));
(void) time(&t2); printf("\nSuccess! Time to taken to run = %d seconds\nMC steps per second = %d\n", (int)(t2-t1), (int)(steps/(double)(t2-t1))); 
return 0;

}

//*****************************************************************************************************************

int mstate(double* S00, double* S01, double* a, double* old_a, int N, int s_l, int s_r){
	//choose a state with s within the correct window
	int i, j, N2, count = 0, count_n_bin, greater_or_less = 1;
	N2 = N*N;
	double theta, s, delta_theta, denom, s00, s01, tS00, tS01, diff = 1;
	count_n_bin = theta = 0; 
	if(s_r == N2)
		delta_theta = 0;
	else
		delta_theta = twopi/(N2 - (float)(s_r));
	for(i = 0;i < N2;i ++)
		{
		count_n_bin ++;
		if(count_n_bin < ((float)s_r))	//starts off with s inside the specified range
			old_a[i] = a[i] = 0;		
		else{
			theta += delta_theta;
			old_a[i] = a[i] = theta;
			}
		}
	for(i = 0;i < N2;i ++)
		a[i] = old_a[i] = 0;
	denom = 5.;
	TRYAGAIN:;
	for(i = 0;i < N2;i ++){
		if(greater_or_less==1)
			old_a[i] = a[i];
		a[i] = fmod(old_a[i] + (r-0.5)/denom,6.28318531);
			}
	tS00 = tS01 = 0;
	for(i = 0;i < N2;i ++)
		{
		tS00 += 2*pow(cos(a[i]),2) - 1;
		tS01 += 2*cos(a[i])*sin(a[i]);
		}
	s = sqrt(tS00*tS00 + tS01*tS01);
	if(fabs(s-(s_r+s_l)/2.)>diff)
	{
	if((s-(s_r+s_l)/2.)>diff)
		{
		greater_or_less = 1;
		}
	else if((s-(s_r+s_l)/2.)<-diff)
		{
		denom *= 1.5;
		greater_or_less = -1;
		}
	count ++;
	if(count%100==0)
		printf("count = %d\twrong: s = %lf\ts_l = %d\ts_r = %d\n", count, s, s_l, s_r);
	if(count > 1000){
		printf("there's a problem in finding the state with the reqired s-value\n");
		return 0;
		}
	goto TRYAGAIN;
	}
	*S00 = tS00; *S01 = tS01;
	printf("s = %lf\ts_l = %d\ts_r = %d\n", s, s_l, s_r);
	return s;
}

double ener4(double H0, double new_theta, int i, double* a, int**neighbors, double p) //"""Returns the contribution to the energy from a spin at position i in the lattice"""
{
	return-(pow(cos(new_theta - a[neighbors[i][0]]),p) + 
			pow(cos(new_theta - a[neighbors[i][1]]),p) + 
			pow(cos(new_theta - a[neighbors[i][2]]),p) + 
			pow(cos(new_theta - a[neighbors[i][3]]),p) +
			H0*(cos(4*new_theta)+ 1)/4.);
		/*H0*(pow(cos(new_theta),4)
		   +pow(sin(new_theta),4)));*/
}

double total_ener(double H0, int N2, double* a, int**c, double p){	//"""returns the total energy of the lattice configuration"""	
	double summ = 0.0;
	int i;
	for(i = 0;i < N2;i ++)
		summ += 0.5*ener4(H0,a[i],i,a,c,p);
	return summ;
	}
