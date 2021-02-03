# convert LEPHARE .spec files into pdz and specz folders

mkdir -p SPEC_OUT_BROADENED/
for f in SPEC_OUT/*; do
	grep -v '^#' "$f"|head -n1|while read ID zspec zphot
	do
		grep '^  0.00000  ' -B 10000000 -A 0 "$f" -m1|grep -v '^  0.00000  '
		cat smoothened/$ID
		grep '^  0.00000  ' -A 10000000 "$f" |grep -v '^  '
	done > ${f/SPEC_OUT/SPEC_OUT_BROADENED}
done
