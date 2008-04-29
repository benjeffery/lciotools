pdgnames = {1  :  "d"  ,
2  :  "u"  ,
3  :  "s"  ,
4  :  "c"  ,
5  :  "b"  ,
6  :  "t"  ,
7  :  "b'"  ,
8  :  "t'"  ,
21  :  "g"  ,
22  :  "gamma"  ,
23  :  "Z0"  ,
24  :  "W+"  ,
25  :  "h0"  ,
2101  :  "ud_0"  ,
3101  :  "sd_0"  ,
3201  :  "su_0"  ,
211  :  "pi+"  ,
311  :  "K0"  ,
321  :  "K+"  ,
411  :  "D+"  ,
421  :  "D0"  ,
431  :  "D_s+"  ,
511  :  "B0"  ,
521  :  "B+"  ,
531  :  "B_s0"  ,
541  :  "B_c+"  ,
111  :  "pi0"  ,
221  :  "eta"  ,
331  :  "eta'"  ,
441  :  "eta_c"  ,
551  :  "eta_b"  ,
130  :  "K_L0"  ,
310  :  "K_S0"  ,
10213  :  "b_1+"  ,
10313  :  "K_10"  ,
10323  :  "K_1+"  ,
10413  :  "D_1+"  ,
10423  :  "D_10"  ,
10433  :  "D_1s+"  ,
10113  :  "b_10"  ,
10223  :  "h_10"  ,
10333  :  "h'_10"  ,
10443  :  "h_1c0"  ,
20213  :  "a_1+"  ,
20313  :  "K*_10"  ,
20323  :  "K*_1+"  ,
20413  :  "D*_1+"  ,
20423  :  "D*_10"  ,
20433  :  "D*_1s+"  ,
20113  :  "a_10"  ,
20223  :  "f_10"  ,
20333  :  "f'_10"  ,
20443  :  "chi_1c0"  ,
100443  :  "psi'"  ,
100553  :  "Upsilon'"  ,
2112  :  "n0"  ,
2212  :  "p+"  ,
3112  :  "Sigma-"  ,
3122  :  "Lambda0"  ,
3212  :  "Sigma0"  ,
3222  :  "Sigma+"  ,
3312  :  "Xi-"  ,
3322  :  "Xi0"  ,
4112  :  "Sigma_c0"  ,
4122  :  "Lambda_c+"  ,
4212  :  "Sigma_c+"  ,
4222  :  "Sigma_c++"  ,
4132  :  "Xi_c0"  ,
4312  :  "Xi'_c0"  ,
4232  :  "Xi_c+"  ,
4322  :  "Xi'_c+"  ,
4332  :  "Omega_c0"  ,
5112  :  "Sigma_b-"  ,
5122  :  "Lambda_b0"  ,
5212  :  "Sigma_b0"  ,
5222  :  "Sigma_b+"  ,
11  :  "e-"  ,
12  :  "nu_e"  ,
13  :  "mu-"  ,
14  :  "nu_mu"  ,
15  :  "tau-"  ,
16  :  "nu_tau"  ,
17  :  "tau'"  ,
18  :  "nu'_tau"  ,
32  :  "Z'0"  ,
33  :  'Z"0'  ,
34  :  "W'+"  ,
35  :  "H0"  ,
36  :  "A0"  ,
37  :  "H+"  ,
39  :  "Graviton"  ,
41  :  "R0"  ,
42  :  "LQ"  ,
1103  :  "dd_1"  ,
2103  :  "ud_1"  ,
2203  :  "uu_1"  ,
3103  :  "sd_1"  ,
3203  :  "su_1"  ,
3303  :  "ss_1"  ,
213  :  "rho+"  ,
313  :  "K*0"  ,
323  :  "K*+"  ,
413  :  "D*+"  ,
423  :  "D*0"  ,
433  :  "D*_s+"  ,
513  :  "B*0"  ,
523  :  "B*+"  ,
533  :  "B*_s0"  ,
543  :  "B*_c+"  ,
113  :  "rho0"  ,
223  :  "omega"  ,
333  :  "phi"  ,
443  :  "J/psi"  ,
553  :  "Upsilon"  ,
10211  :  "a_0+"  ,
10311  :  "K*_00"  ,
10321  :  "K*_0+"  ,
10411  :  "D*_0+"  ,
10421  :  "D*_00"  ,
10431  :  "D*_0s+"  ,
10111  :  "a_00"  ,
10221  :  "f_00"  ,
10331  :  "f'_00"  ,
10441  :  "chi_0c0"  ,
215  :  "a_2+"  ,
315  :  "K*_20"  ,
325  :  "K*_2+"  ,
415  :  "D*_2+"  ,
425  :  "D*_20"  ,
435  :  "D*_2s+"  ,
115  :  "a_20"  ,
225  :  "f_20"  ,
335  :  "f'_20"  ,
445  :  "chi_2c0"  ,
1114  :  "Delta-"  ,
2114  :  "Delta0"  ,
2214  :  "Delta+"  ,
2224  :  "Delta++"  ,
3114  :  "Sigma*-"  ,
3214  :  "Sigma*0"  ,
3224  :  "Sigma*+"  ,
3314  :  "Xi*-"  ,
3324  :  "Xi*0"  ,
3334  :  "Omega-"  ,
4114  :  "Sigma*_c0"  ,
4214  :  "Sigma*_c+"  ,
4224  :  "Sigma*_c++"  ,
4314  :  "Xi*_c0"  ,
4324  :  "Xi*_c+"  ,
4334  :  "Omega*_c0"  ,
5114  :  "Sigma*_b-"  ,
5214  :  "Sigma*_b0"  ,
5224  :  "Sigma*_b+"  ,
}

def pdgToName(pdg):
    try:
        n = name[abs(pdg)]
        if pdg < 0: n += "bar"
        return n
    except KeyError:
        return str(pdg)
