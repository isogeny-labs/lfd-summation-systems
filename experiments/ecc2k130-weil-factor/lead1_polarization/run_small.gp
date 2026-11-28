\\ usage: gp -q < run_small.gp
\r decide_ppav.gp
L = [3, 5, 11, 13];
for (i = 1, #L, t = getabstime(); run(L[i]); printf("  [time %d ms]\n\n", getabstime() - t));
