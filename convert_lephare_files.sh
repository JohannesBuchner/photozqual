# convert LEPHARE .spec files into pdz and specz folders

mkdir -p specz pdz

for f in SPEC_OUT/*.spec; do
	grep -v '^#' "$f"|head -n1|while read ID zspec zphot
	do
		sed s,0\.212200-313,0.212200E-313,g "$f" | grep '^  0.00000  ' -A 10000000 | grep '^  '|head -n -2 > pdz/$ID
		if [[ $zspec != '-99.00000' ]]; then
			echo $zspec > specz/$ID
		fi
	done
done
