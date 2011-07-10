//./a.out 32 29 all_Q_M_32_.dat 200 100 all_P_E_32_.dat all_P_E_32_plot.dat all_Q_M_32_plot.dat H0
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#define twopi 6.2831853
int main(int argc, char *argv[])
{
	FILE 
	//INPUT:
	*fm = fopen(argv[3],"r"), 		//i,ln_s_bin[i],final_s_bin[i]/area
	*fe = fopen(argv[6],"r"), 		//j, ln_e_bin[j], final_e_bin[j]
	//OUTPUT:
	*fm_plot = fopen(argv[8],"w+"), //(double)i/(double)N2, final_s_bin[i]/area (separated by \n\n for plotting)
	*fe_plot = fopen(argv[7],"w+"), //(double)i/(double)N2, final_e_bin[i]/area (separated by \n\n for plotting)
	*fsusc = fopen(argv[10],"w+"),	//start_epsilon+j*delta_epsilon,M,susc,BCM
	*fspec = fopen(argv[11],"w+");//start_epsilon+j*delta_epsilon,E,spec,BCE
	int i, j, N, N2, number_of_inverse_temperatures, e_size, dummy_int;
	double dummy, dummy2, area, M, M2,M4,E,E2,E4,susc,spec,BCM,BCE,ln_mbin_i,start_epsilon,delta_epsilon, H0;
	N = atoi(argv[1]);
	number_of_inverse_temperatures = atoi(argv[2]);
	start_epsilon = (double)(atoi(argv[4]))/1000.;
	delta_epsilon = (double)(atoi(argv[5]))/1000.;
	N2 = N*N;
	H0 = ((double)(atoi(argv[9])))/1000.;
	e_size = (int)(2*N2+ceil((H0/2.)*N2));	//how large the energy range is depends upon the hamiltonian
	double * mbin, *ebin;
	mbin = (double *)calloc((N2),sizeof(double));
	ebin  = (double *)calloc((e_size),sizeof(double));
	for(j = 0;j < number_of_inverse_temperatures;j ++){
		area = 0;
		for(i = 0;i < N2;i ++){
			fscanf(fm,"%d\t%lf\t%lf\n",&dummy_int,&ln_mbin_i,&mbin[i]);
			area += mbin[i];
			}
		for(i = 0;i < N2;i ++)
			fprintf(fm_plot,"%lf\t%lf\n",(double)i/(double)N2,mbin[i]/area*N2);
		fprintf(fm_plot,"\n\n");
		M = M2 = M4 = 0;	//M2 is M^2 etc.
		for(i = 0;i < N2;i ++){
			M += (i+0.5)*(mbin[i]/area);
			M2 += pow((i+0.5),2)*(mbin[i]/area);
			M4 += pow((i+0.5),4)*(mbin[i]/area);
			}
		susc = (M2 - M*M)*(start_epsilon+j*delta_epsilon)/N2;	//note that \epsilon is the inverse temperature
		BCM = M4/(M2*M2);
		printf("M = %lf\tM2 = %lf\tsusc = %lf\t",M,M2,susc);
		fprintf(fsusc,"%lf\t%lf\t%lf\t%lf\n",start_epsilon+j*delta_epsilon,M,susc,BCM);
		
		area = 0;
		for(i = 0;i < e_size;i ++){
			fscanf(fe,"%d\t%lf\t%lf\n",&dummy_int,&dummy,&ebin[i]);
			area += ebin[i];
			}
		for(i = 0;i < e_size;i ++)
			fprintf(fe_plot,"%lf\t%lf\n",(double)i/(double)N2,ebin[i]/area);
		fprintf(fe_plot,"\n\n");
		
		E = E2 = E4 = 0;
		for(i = 0;i < e_size;i ++){
			E += (i+0.5)*(ebin[i]/area);
			E2 += pow((i+0.5),2)*(ebin[i]/area);
			E4 += pow((i+0.5),4)*(ebin[i]/area);
			}
		spec = (E2 - E*E)*pow((start_epsilon+j*delta_epsilon),2)/(double)N2;
		BCE = E4/(E2*E2);
		printf("E = %lf\tE2 = %lf\tspec = %lf\n",E,E2,spec);
		fprintf(fspec,"%lf\t%lf\t%lf\t%lf\n",start_epsilon+j*delta_epsilon,E,spec,BCE);
		
		}
	return 0;
}
