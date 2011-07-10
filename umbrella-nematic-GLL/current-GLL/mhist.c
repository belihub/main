//gcc -lm mhist.c; ./a.out 8 1600 3 30 allebin_8_.dat allembin_8_.dat E1_8_.dat
//gcc -lm mhist.c; ./a.out 20 1600 30 30 all_P_E_20_.dat all_P_ES_20_.dat
//echo 'p "E.dat" u 1:2 index 0, "E.dat" u 1:2 index 1, "E.dat" u 1:2 index 2' | gnuplot -persist
//echo 'p "E1.dat" u 1:2, "E2.dat" u 1:2' | gnuplot -persist
//echo 'p "M1.dat" u 1:2 tit "data", "M2.dat" u 1:2 tit "interpolated"' | gnuplot -persist 
//printf("errorcount = %d\n", errorcount); errorcount ++;
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
int main(int argc, char *argv[])
{
double maxZ, minZ, A, epsilon, T0, dT, tmp1, tmp2 = 0, ener, E, M, M2, M4, M2_vink, M4_vink, norm_p_m;
int e, m, i, j, k, l, s, sum0, N, N2, d1 = 0, d2 = 0, multihist = 1, notdone = 1, RANGE, RANGE1, errorcount = 0, dummy1;
N = atoi(argv[1]);
N2 = N*N;
epsilon = ((double)(atoi(argv[2]))/1000.);
dT = 	  ((double)(atoi(argv[3]))/1000.);
RANGE = (atoi(argv[4]));
RANGE1 = RANGE*2;
FILE *fte1, *ftm1, *fte2, *ftm2, *fe, *fem, *prob, *prob2;
fe =  fopen(argv[5],"r");	//allebin
fem = fopen(argv[6],"r");	//allembin
	fte1 = fopen("E1.dat","w+");	//directly from simulation
	ftm1 = fopen("M1.dat","w+");	//directly from simulation
	fte2 = fopen("E2.dat","w+");	//interpolated
	ftm2 = fopen("M2.dat","w+");	//interpolated
if(notdone){
	prob = fopen("prob.dat","w+");
	}
else{
	prob = fopen("prob.dat","r");
	prob2= fopen("prob2.dat","w+");
	}

double *** embin, *** embin_trans, ** ebin, ** mbin, * numer, * numerator, * n, * n_old, * rho, * sum, * beta, * sum1, * beta1, * Z, * Z1, * Zold, * U, * S, * area, * area1, ** p_e, ** p_m;
embin = (double ***)calloc(RANGE,sizeof(double **));
for(i=0;i<RANGE;i++){
	embin[i] = (double **)calloc(N2,sizeof(double *));		//row indices for M
	for(j=0;j<N2;j++) {
		embin[i][j]=(double *)calloc((2*N2),sizeof(double));//column indices for E
		}
	}
embin_trans = (double ***)calloc(RANGE,sizeof(double **));
for(i=0;i<RANGE;i++){
	embin_trans[i] = (double **)calloc(2*N2,sizeof(double *));		//row indices for M
	for(j=0;j<2*N2;j++) {
		embin_trans[i][j]=(double *)calloc((N2),sizeof(double));//column indices for E
		}
	}
ebin = (double **)calloc(RANGE,sizeof(double *));
for(i = 0;i < RANGE;i ++){
	ebin[i] = (double *)calloc(2*N2,sizeof(double));
	}
mbin = (double **)calloc(RANGE,sizeof(double *));
for(i = 0;i < RANGE;i ++){
	mbin[i] = (double *)calloc(2*N2,sizeof(double));
	}
numer = (double *)calloc(2*N2,sizeof(double));
numerator = (double *)calloc(2*N2,sizeof(double));
n = (double *)calloc(RANGE,sizeof(double));
n_old = (double *)calloc(RANGE,sizeof(double));
rho = (double *)calloc(RANGE,sizeof(double));
sum = (double *)calloc(RANGE,sizeof(double));
beta = (double *)calloc(RANGE,sizeof(double));
sum1 = (double *)calloc(RANGE1,sizeof(double));
beta1 = (double *)calloc(RANGE1,sizeof(double));
Z = (double *)calloc(RANGE,sizeof(double));
Zold = (double *)calloc(RANGE,sizeof(double));
Z1 = (double *)calloc(RANGE1,sizeof(double));
U = (double *)calloc(RANGE1,sizeof(double));
S = (double *)calloc(RANGE1,sizeof(double));
area = (double *)calloc(RANGE,sizeof(double));
area1 = (double *)calloc(RANGE1,sizeof(double));
p_e = (double **)calloc(RANGE1,sizeof(double *));
for(i = 0;i < RANGE1;i ++)
	p_e[i] = (double *)calloc(2*N2,sizeof(double));
p_m = (double **)calloc(RANGE1,sizeof(double *));
for(i = 0;i < RANGE1;i ++)
	p_m[i] = (double *)calloc(  N2,sizeof(double));

for(k = 0;k < RANGE;k ++)
	beta[k] = epsilon + k*dT;

//beta1[k] falls between beta[k] and beta[k+1]
for(k = 0;k < RANGE1;k ++)
	beta1[k] = epsilon+k*dT*RANGE/((float)RANGE1);	

//begin input of data to the respective arrays
//111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
double num, denom, diff;
for(k = 0;k < RANGE;k ++){
	ener = 0;
	for(i = 0;i < 2*N2;i ++){
		fscanf(fe,"%d\t%lf\t%lf\n",&d1,&tmp2,&ebin[k][i]);
		area[k] += ebin[k][i];
		ener += i*ebin[k][i];
		}
	//printf("area[%d] = %lf\n", k, area[k]);
	if(area[k] != 0)
		ener /= area[k];
	else 
		ener = 0; 
	fprintf(fte1,"%lf\t%lf\n",beta[k],ener);
	}fprintf(fte1,"\n\n");
for(k = 0;k < RANGE;k ++){
	for(j = 0;j < 2*N2;j ++){
		ebin[k][i] = 0;
		for(i = 0;i < N2;i ++){
			fscanf(fem,"%d\t%d\t%lf\n",&d1,&d2,&embin[k][i][j]);
			ebin[k][j] += embin[k][i][j];	//for checking
			}
		}
	}
/*
//this is where we check:
for(k = 0;k < RANGE;k ++){
	ener = 0;
	area[k] = 0;
	for(i = 0;i < 2*N2;i ++){
		area[k] += ebin[k][i];
		ener += i*ebin[k][i];
		}if(area[k] != 0) ener /= area[k]; else ener = 0; fprintf(fte1,"%lf\t%lf\n",beta[k],ener); 
	}fprintf(fte1,"\n\n");
//end check
*/
for(k = 0;k < RANGE;k ++){
	M = M2 = M4 = M2_vink = M4_vink = area[k] = 0;
	for(i = 0;i < N2;i ++){
		for(j = 0;j < 2*N2;j ++){
			mbin[k][i] += embin[k][i][j];
			}
		area[k] += mbin[k][i];
		M += i*mbin[k][i];
		M2 += i*i*mbin[k][i];
		M4 += pow(i,4)*mbin[k][i];
		}
	if(area[k] > 0){ M /= area[k]; M2 /= area[k]; M4 /= area[k];}
	else M = M2 = M4 = 0; 
	for(i = 0;i < N2;i ++){
		M2_vink += (i-M)*(i-M)*mbin[k][i];
		M4_vink += pow((i-M),4)*mbin[k][i];
		}
	if(area[k] > 0){M2_vink /= area[k]; M4_vink /= area[k];}
	fprintf(ftm1,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",beta[k],M,M2,M4,M4/(M2*M2),M4_vink/(M2_vink*M2_vink)); 
	}
	
	printf("\ndone!\n");

//22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
if(notdone){ 	//do the multiple histogram analysis, calculate stuff and print p_m to file
for(k = 0;k < RANGE;k ++){
Zold[k] = 1000;
Z[k] = 0;
for(e = 0;e < 2*N2;e ++)
	n[k] += ebin[k][e];	//ebin[k][e] is N_k(E)
}

for(e = 0;e < 2*N2;e ++)
	for(i = 0;i < RANGE;i ++)
		numerator[e] += ebin[i][e];

do{
maxZ = 0;
for(k = 0;k < RANGE;k ++){
	Z[k] = 0;
	for(e = 0;e < 2*N2;e ++){
		num = denom = 0;
		//for(k = 0;k < RANGE;k ++)
		//	num += ebin[k][e];
		for(j = 0;j < RANGE;j ++)
			denom += n[j]*exp((beta[k]-beta[j])*(e))/Zold[j];
		Z[k] += numerator[e]/denom;
		}
	if(maxZ < Z[k])
		maxZ = Z[k];
	if(Z[k] == 0){
		printf("stopped at k=%d\n", k);
		goto HERE;
		}
	}
minZ = -Z[0];
for(k = 0;k < RANGE;k ++)
	if(minZ < -Z[k])
		minZ = -Z[k];
//"normalization:"
A = 1/(sqrt(-minZ*maxZ));
diff = 0;
for(k = 0;k < RANGE;k ++){	
	diff += pow(((Z[k]-Zold[k])/Zold[k]),2);
	}
printf("diff = %lf\n",diff);
for(k = 0;k < RANGE;k ++){
	Zold[k] = A*Z[k];
	}
} while (diff > 0.0000000001);	
for(k = 0;k < RANGE;k ++)
printf("Z[%d]=%e\n",k,Z[k]);

//interpolated values:
printf("\n<U>:\t");
for(k = 0;k < RANGE1;k ++){
	Z1[k] = 0;
	for(e = 0;e < 2*N2;e ++){
		num = denom = 0;
		for(j = 0;j < RANGE;j ++)
			num += ebin[j][e];
		for(j = 0;j < RANGE;j ++)
			denom += n[j]*exp((beta1[k]-beta[j])*(e))/Z[j];
		Z1[k] += num/denom;
		p_e[k][e] = num/denom;
		}
	for(e = 0;e < 2*N2;e ++){
		p_e[k][e] /= Z1[k];
		fprintf(prob,"%d\t%lf\n",e,p_e[k][e]);
		U[k] += e*p_e[k][e];
		}
	printf("k=%d, Z=%e, U=%lf\n",k,Z1[k],U[k]);
	fprintf(fte2,"%lf\t%lf\n",beta1[k],U[k]);
	}fprintf(fte2,"\n\n");fprintf(prob,"\n\n");
for(k=0;k<RANGE;k++)
	printf("%lf\n",n[k]);
printf("\ndone\n");

for(k=0;k<RANGE;k++){
	sum1[k] = 0;
	//Z1[j] = 0;
	}
//3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333

printf("Begin multihist!\n");
//first, just to check, the energy:
/*
for(k = 0;k < RANGE1;k ++){
	Z1[k] = 0;
	for(e = 0;e < 2*N2;e ++){
		num = denom = 0;
		for(j = 0;j < RANGE;j ++)
			for(m = 0;m < N2;m ++)
				num += embin[j][m][e];
		for(j = 0;j < RANGE;j ++)
			denom += n[j]*exp((beta1[k]-beta[j])*(e))/Z[j];
		Z1[k] += num/denom;		//we don't get a match using the old Z1[k], but
		p_e[k][e] = num/denom;	//we get a perfect match using a re-calculated Z1[k] and dividing by 2
		//U[k] += e*num/denom/Z1[k];	//...due to normalisation?
		}
	for(e = 0;e < 2*N2;e ++){
		p_e[k][e] /= Z1[k];
		U[k] += e*p_e[k][e];
		}
	printf("j=%d\tU=%lf\n",k,U[k]);
	fprintf(fte2,"%lf\t%lf\n",beta1[k],U[k]/2.);
	}*/
//then the magnetization
for(k = 0;k < RANGE1;k ++){
	area1[k] = 0;
	for(m = 0;m < N2;m ++){
		for(e = 0;e < 2*N2;e ++){
			num = 0;
			for(i = 0;i < RANGE;i ++)
				num += embin[i][m][e];
			denom = 0;
			for(j = 0;j < RANGE;j ++)
				denom += n[j]*exp((beta1[k]-beta[j])*(e))/Z[j];
			p_m[k][m] += num/denom;
			}
		p_m[k][m] /= Z1[k];
		area1[k] += p_m[k][m];
		S[k] += m*p_m[k][m];
		}//printf("area1[%d] = %lf\n", area1[k]);
	for(m = 0;m < N2;m ++){
		p_m[k][m] /= area1[k];
		}
	printf("k=%d, Z=%e, S=%lf\n",k,Z1[k],S[k]);
	}
for(k = 0;k < RANGE1;k ++){
	for(m = 0;m < N2;m ++){
		fprintf(prob,"%d\t%lf\n",m,p_m[k][m]);
		}
	}
//end multihist 

for(k = 0;k < RANGE1;k ++){
	M = M2 = M4 = M2_vink = M4_vink = norm_p_m = 0;
	for(m = 0;m < N2;m ++){
		norm_p_m += p_m[k][m];
		M += m*p_m[k][m];
		M2 += pow(m,2)*p_m[k][m];
		M4 += pow(m,4)*p_m[k][m];
		}
//printf("norm_p_m at %d = %lf\n", k, norm_p_m);
	for(m = 0;m < N2;m ++){
		M2_vink += (m-M)*(m-M)*p_m[k][m];
		M4_vink += pow((m-M),4)*p_m[k][m];
		}
	fprintf(ftm2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", beta1[k], M, M2, M4, M4/(M2*M2), M4_vink/(M2_vink*M2_vink));
	}

}//end notdone

else{	//read p_m off file and calculate the quantities we're interested in
for(k = 0;k < RANGE1;k ++)
	for(e = 0;e < 2*N2;e ++){
		fscanf(prob,"%d\t%lf\n",&dummy1,&p_e[k][e]);
		fprintf(prob2,"%d\t%lf\n",dummy1,p_e[k][e]);
		}
fscanf(prob,"\n\n");
fprintf(prob2,"\n\n");
for(k = 0;k < RANGE1;k ++)
	for(m = 0;m < N2;m ++){
		fscanf(prob,"%d\t%lf\n",&dummy1,&p_m[k][m]);
		fprintf(prob2,"%d\t%lf\n",dummy1,p_m[k][m]);
		}
for(k = 0;k < RANGE1;k ++){
	M = M2 = M4 = M2_vink = M4_vink = norm_p_m = 0;
	for(m = 0;m < N2;m ++){
		norm_p_m += p_m[k][m];
		M += m*p_m[k][m];
		M2 += pow(m,2)*p_m[k][m];
		M4 += pow(m,4)*p_m[k][m];
		}
//printf("norm_p_m at %d = %lf\n", k, norm_p_m);
	for(m = 0;m < N2;m ++){
		M2_vink += (m-M)*(m-M)*p_m[k][m];
		M4_vink += pow((m-M),4)*p_m[k][m];
		}
	fprintf(ftm2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", beta1[k], M, M2, M4, (M2*M2)/M4, (M2_vink*M2_vink)/M4_vink);
	}

for(k = 0;k < RANGE1;k ++){
	E = 0;
	for(e = 0;e < 2*N2;e ++){
		E += e*p_e[k][e];
		}
	fprintf(fte2,"%lf\t%lf\n", beta1[k], E);
}

}

HERE: ;
return 0;
//printf("\nError: divide by zero - aborted\n");
}
