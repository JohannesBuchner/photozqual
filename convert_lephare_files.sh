# convert LEPHARE .spec files into pdz and specz folders

for f in SPEC_OUT/*; do
	grep -v '^#' "$f"|head -n1|while read ID zspec zphot
	do
		grep '^  0.00000  ' -A 10000000 "$f" | grep '^  '|head -n -2 > pdz/$ID
		if [[ $zspec != '-99.00000' ]]; then
			echo $zspec > specz/$ID
		fi
	done
done
