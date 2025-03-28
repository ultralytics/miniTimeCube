% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

%Function taken from the NIST website:  
%http://physics.nist.gov/cgi-bin/Xcom/xcom2
%Determining the throughput of gammas in lead.  
%ONLY valid for gamma radiation in Lead.

%Input:
%MeV_in:[MeV]   energy of the incident gamma photon.  Linear interpolation
%               will be used to get the value of mu for energies not in the 
%               given array.  valid ranges are 0.001 MeV to 1000 MeV
%x:[cm]         thickness of the lead shield in cm.

%Output:
%I_fract: [-]   the FRACTIONAL transmission of gamma flux at given energy
%               and shield thickness.


function I_fract = getGammaTransInLead(MeV_in, x)

MeV = [0.001
0.0015
0.002
0.002484
0.002534
0.002586
0.003
0.003066
0.003301
0.003554
0.003699
0.003851
0.004
0.005
0.006
0.008
0.01
0.01304
0.015
0.0152
0.01553
0.01586
0.02
0.03
0.04
0.05
0.06
0.08
0.088
0.1
0.15
0.2
0.3
0.4
0.5
0.6
0.8
1
1.022
1.25
1.5
2
2.044
3
4
5
6
7
8
9
10
11
12
13
14
15
16
18
20
22
24
26
28
30
40
50
60
80
100
150
200
300
400
500
600
800
1000
];

mu = [58933.98
26580.96
14447.16
12332.25
18552.24
24789.24
22169.7
22583.61
20207.88
17355.87
16238.88
15082.2
14084.28
8190.882
5214.132
2525.418
1425.438
1255.6215
1228.122
1416.366
1566.054
1604.61
953.0136
328.1796
152.4096
83.76858
51.3702
23.95008
51.41556
60.51024
21.6594
10.61424
4.233222
2.435832
1.699866
1.323378
0.9534672
0.7714602
0.7570584
0.6443388
0.5769792
0.513702
0.5107536
0.47628
0.4737852
0.483084
0.4969188
0.5127948
0.529578
0.5464746
0.5634846
0.5799276
0.5961438
0.6119064
0.627102
0.6413904
0.6549984
0.6804
0.703647
0.7249662
0.7449246
0.763182
0.7803054
0.7962948
0.862974
0.9135504
0.9534672
1.0131156
1.055754
1.1236806
1.164618
1.212246
1.239462
1.257606
1.271214
1.288224
1.299564
];

this_mu = interp1(MeV, mu, MeV_in, 'linear');

I_fract = exp(-this_mu*x);


return, I_fract
