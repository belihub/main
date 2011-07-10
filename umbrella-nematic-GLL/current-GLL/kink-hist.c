/* this transforms the graph into one with two peaks of equal area */
//gcc -lm kink-hist.c; ./a.out 10 140 P_M_10_140.dat kink_10_140.dat binder_M_10_140.dat; echo 'p "kink_10_140..dat" u 1:2' | gnuplot -persist
//printf("steps = %d\n", steps); steps ++;
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
int main(int argc, char *argv[])
{
	FILE *ftest = fopen("test.dat","w+"),
	*fp3 = fopen(argv[3],"r"),		//>P_M_N_epsilon.dat
	*fp4 = fopen(argv[4],"w+"),		//>kink_N_epsilon.dat
	*fp5 = fopen(argv[5],"w+");		//>binder_M_N_epsilon.dat
	int i, i1 = 2, i2 = 91, i_min, j, k, N, N2, steps = 0, counter;
	double epsilon, S, S0, M, M2, M4, E, binder_M, area, slope, delta, min = 1000000, area_l, area_r;
	N = atoi(argv[1]);
	N2 = N*N;
	epsilon = ((double)(atoi(argv[2]))/1000.);
	double * ln_s_bin;
	ln_s_bin = (double *)calloc((N2),sizeof(double));
	double * final_s_bin;
	final_s_bin = (double *)calloc((N2),sizeof(double));
	double * test_s_bin;
	test_s_bin = (double *)calloc((N2),sizeof(double));

	for(i = 0;i < N2;i ++){
		fscanf(fp3,"%d\t%lf\t%lf\n",&j,&ln_s_bin[i],&final_s_bin[i]);
		if(j != i)
			{
			printf("j != i");
			goto HERE;
			}
		}

/* create two peaks of equal areas once you have the averaged distribution*/
counter = 0;
slope = 0;
delta = 0.01;
int sign = 1, old_sign = 1;	//>shows which lobe is/was greater; if the signs are different then we've overshot

do{
	min = 100000;
	area_l = area_r = 0;
	for(i = 0;i < N2;i ++){
		test_s_bin[i] = ln_s_bin[i] - slope*i;
		if((i > N2/2-N2/4)&&(i < N2/2+N2/4)&&(min > test_s_bin[i])){
			min = test_s_bin[i];	//>the value of the minimum
			i_min = i;				//>the position of the minimum
			}
		}
	for(i = 0;i < i_min;i ++)
		area_l += test_s_bin[i] - min;	//>area of the left lobe
	for(i = i_min;i < N2;i ++)
		area_r += test_s_bin[i] - min;	//>area of the right lobe
	sign = ((area_r > area_l)?1:(-1));
	if((counter > 0)&&((sign*old_sign) == -1)){
		delta /= 2.;					//>go the other way with smaller delta
		printf("slope = %lf\n",slope);	
		}
	slope += delta*sign;				//>if the right peak is larger, increase the slope, else vice versa
	counter ++;
	if(counter > 100000){
		printf("counter limit reached\n");
		break;
		}
	old_sign = sign;
	} while (fabs(area_l - area_r) > 0.01);
	
printf("area_r - area_l = %lf\n",area_r-area_l);
for(i = 0;i < N2;i ++)
	fprintf(fp4,"%d\t%lf\n",i,test_s_bin[i] - min);
HERE:;
return 0;
}
