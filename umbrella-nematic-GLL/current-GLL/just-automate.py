#N0 = 64		#linear system dimension L = N0
T0 = 900	#starting \epsilon times 1000
dT = 20		#the increment in inverse temperature \epsilon
NUM_TEMP=30	#the number of different inverse temperatures to be simulated
P = 10
#H0 = 3000	#the factor in the hamiltonian H times 1000
dH = 3000	#the increment in H times 1000
NUM_H = 1	#the number of different H's to bu simulated
hours = 0
minutes = 2
seconds = 0
time_in_seconds = hours*3600 + minutes*60 + seconds

import sys
import os

N0 = int(sys.argv[2])
H0 = int(sys.argv[3])
#time_in_seconds = int(sys.argv[3])

if sys.argv[1] == '1':
	os.system('rm -rf *.out')

	#gcc -lm mt19937-64.c fit.c umbrella-magnetic.c -I. -lgsl -lgslcblas; ./a.out 8 1000 P_S.dat P_E.dat P_ES.dat weights_etc.dat 2000 ; echo 'p "P_S.dat" u 1:3' | gnuplot -persist
	#key: ./a.out N0 T0+j*dT P_S.dat P_E.dat P_ES.dat weights_etc.dat H0+k*dH

	os.system('gcc -lm mt19937-64.c fit.c func-umbrella-nematic.c -I. -lgsl -lgslcblas')

	for k in range(0,NUM_H):
		for j in range(0,NUM_TEMP):
			label_P_H_N_T = str(P) + '_' + str(H0+k*dH) + '_' + str(N0) + "_" + str(T0+j*dT) + '.dat'
			cmd = './a.out ' + str(N0) + " " + str(T0+j*dT) + ' P_E_' + label_P_H_N_T + ' P_S_' + label_P_H_N_T + ' P_ES_' + label_P_H_N_T + ' s_bin_and_weights_' + label_P_H_N_T + " " + str(H0+k*dH) + ' ' + str(time_in_seconds) + ' ' + str(P) + ' &'
			print cmd
			os.system(cmd)

elif sys.argv[1] == '2':
	os.system('rm -rf *.out')
	os.system('rm -rf *.ps')
	os.system('rm -rf susc.dat')
	os.system('rm -rf spec.dat')
	j = 0
	k = 0
	notdone = 1
	#notdone is 1 if we're running second-automate.py for the first time after a first-automate.py run,
	#the individual files generated for each combination of parameters T and H by func-umbrella-nematic.c 
	#are concatenated so that they can be analysed together.
	#otherwise, notdone=0 so that we won't waste time re-concatenating the files.

	if notdone == 1:
		os.system('rm -rf all_P_S_*')
		os.system('rm -rf all_P_M_*')
		os.system('rm -rf all_P_E_*')
		for k in range(0,NUM_H):
			for j in range(0,NUM_TEMP):
				label_P_H_N = str(P) + '_' + str(H0+k*dH) + '_' + str(N0) + '.dat'
				label_P_H_N_T = str(P) + '_' + str(H0+k*dH) + '_' + str(N0) + "_" + str(T0+j*dT) + '.dat'
				cmd = 'cat' + ' P_S_' + label_P_H_N_T + ' >> all_P_S_' + label_P_H_N
				print cmd
				os.system(cmd)
				cmd = 'cat' + ' P_E_' + label_P_H_N_T + ' >> all_P_E_' + label_P_H_N
				print cmd
				os.system(cmd)
				#cmd = 'cat' + ' P_ES_' + label_P_H_N_T + ' >> all_P_ES_' + label_P_H_N
				#print cmd
				#os.system(cmd)

	#if we're simulating more than one value of H, the following code should be put inside an appropriate for loop.
	#anal-hist.c analyses the data to give susceptibilities, specific heats, binder cumulants etc.
	label_P_H_N = str(P) + '_' + str(H0+k*dH) + '_' + str(N0) + '.dat'
	cmd  = 'gcc -lm anal-hist.c'
	print cmd
	os.system(cmd)
	cmd  = "./a.out" + ' ' + str(N0) + ' ' + str(NUM_TEMP) + ' all_P_S_' + label_P_H_N + ' ' + str(T0) + ' ' + str(dT) + ' ' + ' all_P_E_' + label_P_H_N + ' all_P_E_plot' + label_P_H_N + ' all_P_S_plot' + label_P_H_N + ' ' + str(H0) + ' susc' + label_P_H_N + ' spec' + label_P_H_N
	print cmd
	os.system(cmd)

	#the files are plotted for immediate checking
	cmd = 'echo \'p "susc' + label_P_H_N + '" u 1:3\' | gnuplot -persist'
	print cmd
	os.system(cmd)

	cmd = 'echo \'p "spec' + label_P_H_N + '" u 1:3\' | gnuplot -persist'
	print cmd
	os.system(cmd)

	cmd = 'echo \'p "susc' + label_P_H_N + '" u 1:4 tit "BCM"\' | gnuplot -persist'
	print cmd
	os.system(cmd)

	cmd = 'echo \'p "spec' + label_P_H_N + '" u 1:4 tit "BCE"\' | gnuplot -persist'
	print cmd
	os.system(cmd)
	
	#the titles are self-explanatory
	plotfile = '"all_P_S_plot' + label_P_H_N + '"'
	pcmd = 'echo \'set term postscript enhanced color; set xlabel "S"; set ylabel "P(S)"; set title "P(S) for different epsilons, for the 2-D GLL model with L = ' + str(N0) + ' and H = ' + str(H0/1000+k*dH) + '"; p '
	for j in range(0,NUM_TEMP):
		pcmd = pcmd + plotfile + ' u 1:2 index ' + str(j) + ' tit "' + str((T0+j*dT)/1000.) + '"'
		if j < NUM_TEMP-1:
			pcmd = pcmd + ', '
	pcmd = pcmd + '\' | gnuplot -persist >> P_M.ps'
	print pcmd
	os.system(pcmd)

	plotfile = '"all_P_E_plot' + label_P_H_N + '"'
	pcmd = 'echo \'set term postscript enhanced color; set xlabel "E"; set ylabel "P(E)"; set title "P(E) for different epsilons, for the 2-D GLL model with L = ' + str(N0) + ' and H = ' + str(H0/1000+k*dH) + '"; p '
	for j in range(0,NUM_TEMP):
		pcmd = pcmd + plotfile + ' u 1:2 index ' + str(j) + ' tit "' + str((T0+j*dT)/1000.) + '"'
		#pcmd = pcmd + '"P_E_' + label_P_H_N_T + '" u 1:3 tit "' + str((T0+j*dT)/1000.) + '"'
		if j < NUM_TEMP-1:
			pcmd = pcmd + ', '
	pcmd = pcmd + '\' | gnuplot -persist >> P_E_.ps'
	print pcmd
	os.system(pcmd)

	plotfile = '"all_P_S_plot' + label_P_H_N + '"'
	pcmd = 'echo \'set xlabel "S"; set ylabel "P(S)"; set title "P(S) for different epsilons, for the 2-D GLL model with L = ' + str(N0) + ' and H = ' + str(H0/1000+k*dH) + '"; p '
	for j in range(0,NUM_TEMP):
		pcmd = pcmd + plotfile + ' u 1:2 index ' + str(j) + ' tit "' + str((T0+j*dT)/1000.) + '"'
		if j < NUM_TEMP-1:
			pcmd = pcmd + ', '
	pcmd = pcmd + '\' | gnuplot -persist'
	print pcmd
	os.system(pcmd)

	plotfile = '"all_P_E_plot' + label_P_H_N + '"'
	pcmd = 'echo \'set xlabel "E"; set ylabel "P(E)"; set title "P(E) for different epsilons, for the 2-D GLL model with L = ' + str(N0) + ' and H = ' + str(H0/1000+k*dH) + '"; p '
	for j in range(0,NUM_TEMP):
		pcmd = pcmd + plotfile + ' u 1:2 index ' + str(j) + ' tit "' + str((T0+j*dT)/1000.) + '"'
		#pcmd = pcmd + '"P_E_' + label_P_H_N_T + '" u 1:3 tit "' + str((T0+j*dT)/1000.) + '"'
		if j < NUM_TEMP-1:
			pcmd = pcmd + ', '
	pcmd = pcmd + '\' | gnuplot -persist'
	print pcmd
	os.system(pcmd)

	cmd = 'echo \'set term postscript enhanced color; set xlabel "S"; set ylabel "P(S)"; set title "GLL model with L = ' + str(N0) + ' and H = ' + str(H0) + '"; p "all_P_S_' + label_P_H_N + '" u 1:3 notitle\' | gnuplot -persist >> P_M.ps'
	print cmd
	os.system(cmd)

	cmd = 'echo \'set term postscript enhanced color; set xlabel "T"; set ylabel "(<S^2> - <S><S>)/N2/T"; set title "susceptibility per spin for the GLL model with L = ' + str(N0) + 'and H = ' + str(H0) + '"; p "susc' + label_P_H_N + '" u 1:3 notitle\' | gnuplot -persist >> susc.ps'
	print cmd
	os.system(cmd)

	cmd = 'echo \'set term postscript enhanced color; set xlabel "T"; set ylabel "(<E^2> - <E><E>)/N2/(T^2)"; set title "specific heat per spin for the GLL model with L = ' + str(N0) + 'and H = ' + str(H0) + '"; p "spec' + label_P_H_N + '" u 1:3 notitle\' | gnuplot -persist >> spec.ps'
	print cmd
	os.system(cmd)

	cmd = 'echo \'set term postscript enhanced color; set xlabel "T"; set ylabel "<S^4> /<S^2><S^2>"; set title "BC for M for the GLL model with L = ' + str(N0) + 'and H = ' + str(H0) + '"; p "susc' + label_P_H_N + '" u 1:4 notitle\' | gnuplot -persist >> BCM.ps'
	print cmd
	os.system(cmd)

	cmd = 'echo \'set term postscript enhanced color; set xlabel "T"; set ylabel "<E^4> /<E^2><E^2>"; set title "BC for E for the GLL model with L = ' + str(N0) + 'and H = ' + str(H0) + '"; p "spec' + label_P_H_N + '" u 1:4 notitle\' | gnuplot -persist >> BCE.ps'
	print cmd
	os.system(cmd)
	
	#the data for different temperatures can be used to calculate the density of states in between those temperatures, if they're close enough; this is done by mhist.c.
	'''
	#gcc -lm mhist.c; ./a.out 20 1600 30 30 all_P_E_20_.dat all_P_ES_20_.dat
	cmd  = 'gcc -lm mhist.c'
	print cmd
	os.system(cmd)
	cmd  = "./a.out" + ' ' + str(N0) + ' ' + str(T0) + ' ' + str(dT) + ' ' + str(NUM_TEMP) + ' ' + ' all_P_E_' + label_P_H_N + ' all_P_ES_' + label_P_H_N
	print cmd
	os.system(cmd)

	cmd = 'echo \'p "M1.dat" u 1:2 tit "data", "M2.dat" u 1:2 tit "interpolated"\' | gnuplot -persist'
	print cmd
	os.system(cmd)
	'''
