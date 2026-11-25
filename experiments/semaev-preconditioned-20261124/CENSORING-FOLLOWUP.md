# Secondary interval analysis of censored baselines

Added after the first held-out split timeouts were observed, while the
prespecified run was still in progress. The frozen protocol and candidate
are unchanged. The strict complete-pair gate remains separately reported;
this is an explicitly secondary analysis, not a relabeling of that gate.

The split and original baselines can time out while the primary candidate
finishes. Discarding those cases would hide potentially important gains,
but replacing a timeout by its cap would invent an exact runtime. Use
intervals instead:

- Successful repetition with charged time t: [t,t].
- GB timeout: [10,infinity], because the charged GB phase alone reaches
  its ten-second real-time alarm. Do not use the partial worker wall time,
  which also contains imports not charged to the primary metric.
- Error, process timeout or unrun row: [0,infinity]. The forty-second
  process cap is not a charged-solver-time lower bound, since imports and
  startup are outside that metric.

For three repetitions, monotonicity of the median gives the valid interval
[median(lower endpoints), median(upper endpoints)]. For the faster of two
baseline medians, take the minimum of their lower and upper endpoints.
For a candidate interval [cL,cU] and best-baseline interval [bL,bU], its
speedup lies in [bL/cU,bU/cL], using zero/infinity for uninformative bounds.
Geometric means are monotone too, so the product of the lower bounds gives
a conservative geometric-mean improvement bound over every case, including
censored cases. No case is removed from that calculation.

Example: original median 6.9 seconds, split median at least 10 seconds,
candidate median 1.1 seconds. The faster baseline is known to be original,
even though the exact split time is unknown. Its gain is determined without
inventing the split runtime. If both baselines are censored at 10 seconds,
the same candidate has only a lower-bound gain of 10/1.1, not an exact gain.

The practical target is supported by this secondary analysis only if the
primary candidate has complete, exact-output-verified repetitions for every
case and the conservative geometric-mean lower bound is at least 2x in each
cohort. Keep the original complete-pair results and target status alongside
the bound analysis. No statistical or asymptotic conclusion follows.
