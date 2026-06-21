# psychnet validation results

Generated 2026-06-21. Reference packages: qgraph 1.9.8, IsingFit 0.4, mgm 1.2.15.

## Part A: EBICglasso vs qgraph (real questionnaire data)

19 datasets; all structure agreement = 1.000, max edge delta <= 0.00801.

## Part B: Ising vs IsingFit (real binary data)

2 datasets; min structure agreement 0.992.

## Part C: synthetic ground truth

psychnet and qgraph agree, and recover the true graph at the same F-measure.

See results_realdata.csv, results_ising.csv, results_synthetic.csv.
