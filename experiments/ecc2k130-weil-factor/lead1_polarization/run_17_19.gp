\\ usage: gp -q < run_17_19.gp     (ell = 17 is the small analogue of ell = 131: (-7/ell) = -1 and U_ell = 271 is of inert type)
\r decide_ppav.gp
L = [17, 19];
for (i = 1, #L, t = getabstime(); run(L[i]); printf("  [time %d ms]\n\n", getabstime() - t));
