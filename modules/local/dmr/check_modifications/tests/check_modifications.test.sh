#!/bin/bash
# Tests for templates/check_modifications.sh. Renders the Nextflow template the way Nextflow does
# (${bedmethyl} -> the fixture path, \$ -> $) and runs it on small bedMethyl fixtures.
set -u
here=$(cd "$(dirname "$0")" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
failed=0

# bedMethyl row with 18 columns: code in column 4, Nother_mod in column 14
row() { printf 'chr1\t100\t101\t%s\t14\t.\t100\t101\t255,0,0\t14\t50.00\t2\t2\t%s\t0\t0\t0\t0\n' "$1" "$2"; }

run_case() { # name, fixture path, expected exit code, expected stderr fragment ("" = none)
    local name=$1 file=$2 want=$3 fragment=$4
    if ! sed -e "s|\${bedmethyl}|$file|g" -e 's|\\\$|$|g' "$here/../templates/check_modifications.sh" > "$tmp/run.sh"; then
        echo "FAIL $name: template could not be rendered"
        failed=1
        return
    fi
    bash "$tmp/run.sh" 2> "$tmp/err"
    local got=$?
    if [ "$got" -eq "$want" ] && { [ -z "$fragment" ] || grep -q "$fragment" "$tmp/err"; }; then
        echo "PASS $name"
    else
        echo "FAIL $name: exit $got (wanted $want), stderr: $(cat "$tmp/err")"
        failed=1
    fi
}

{ row m 10; row h 10; } | gzip > "$tmp/with_5hmc.bed.gz"
{ row m 3; row m 4; } | gzip > "$tmp/no_5hmc_other_mods.bed.gz"
{ row m 0; row m 0; } | gzip > "$tmp/hac_5mc_only.bed.gz"
: | gzip > "$tmp/empty.bed.gz"

run_case "5hmC rows present, N_other > 0: usable" "$tmp/with_5hmc.bed.gz" 0 ""
run_case "no 5hmC rows and N_other > 0: stop" "$tmp/no_5hmc_other_mods.bed.gz" 1 "DMR not skipped"
run_case "HAC 5mC only and N_other = 0: usable" "$tmp/hac_5mc_only.bed.gz" 0 ""
run_case "no rows at all: stop" "$tmp/empty.bed.gz" 1 "no rows"
exit $failed
