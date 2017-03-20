set table "3.pgf-plot.table"; set format "%.5f"
set format "%.7e";; plot '../data/exp3b2.csv' using 2:($3+$2) index 1
