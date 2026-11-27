\r kplus_setup.gp
{ my(t0=getabstime()); for (i=1,3, my(ell=[11,29,37][i]); my(B = setupBoth(ell)); printf("ell=%d ok deg K=%d  %d ms\n", ell, poldegree(B[5]), getabstime()-t0)); }
