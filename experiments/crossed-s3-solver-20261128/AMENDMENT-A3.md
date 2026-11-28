# A3: separate Gröbner basis stage replay

28 November 2026, after the stopped initial attempt and before any GB-only
replay measurements. The preserved stderr shows repeated `Changing variable
order for possibly more generic staircase` messages. The original msolve
command tried to produce a rational parametrization and cycled through input
and timing headers thousands of times. A diagnostic on the same n=3 tree
input with `-c 0` exited after the GB step with a request for a random linear
form. A second diagnostic with `-g 2 -c 0` completed and printed a reduced
grevlex basis. These diagnostics ran while another host workload was active
and are not timing evidence.

The prospective `measure_gb.py` replay uses the identical 27 frozen `.ms`
inputs, target/base selection, 20-second cap, seed, thread count and three
repetitions, but adds `-g 2 -c 0`. It writes to `raw-gb/` and `runs-gb.json`;
it never appends to the stopped `runs.json`. This isolates F4/GB computation
and reports observed matrix work and degree. It does **not** measure rational
parametrization, solution extraction, or full relation-query cost. Timing
starts only after the unrelated msolve process on this host ends. Any
completed GB-stage advantage requires a separate extraction and total-cost
comparison before promotion.
