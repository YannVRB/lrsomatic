#!/bin/bash -ue
# Stop early, with a readable message, when this bedMethyl would make `modkit dmr pair` fail.
#
# modkit dmr pair rejects rows with N_other > 0, and SUP basecalls call both 5mC and 5hmC, so a pileup that has
# other modifications (bedMethyl column 14, Nother_mod) but no 5hmC rows (code `h` in column 4) was made without
# 5hmC. HAC data with no other modifications is fine. Only the first 200000 lines are read; no `pipefail` on
# purpose, because `head` closing the pipe would make zcat fail with SIGPIPE.
file="${bedmethyl}"
verdict=\$(zcat "\$file" | head -n 200000 | awk '
    \$4 == "h" { h = 1 }
    \$14 > 0   { other = 1 }
    END { if (NR == 0) print "empty"; else if (other && !h) print "no5hmc"; else print "ok" }')

case "\$verdict" in
    ok) ;;
    empty)
        echo "\$file: no rows in the bedMethyl file" >&2
        exit 1 ;;
    *)
        echo "\$file: rows with other modifications (N_other > 0) and no 5hmC rows; re-run lrsomatic with 'Phased methylation pileup' ticked" >&2
        exit 1 ;;
esac
