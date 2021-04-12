
if [[ "$1" == "run" ]]
then

if [ ! -e SPEC_OUT ]; then
	echo "ERROR: lephare output directory SPEC_OUT does not exist"
	exit 1
fi

bash convert_lephare_files.sh &&
python3 photozquality.py $(ls pdz) &&
bash convert_back_lephare_files.sh

fi

if [[ "$1" == "clean" ]]
then

rm -rf *.png *.pdf pdz/  prev/  smoothened/ specz/ SPEC_OUT_BROADENED/

fi

if [[ "$1" == "pack" ]]
then

tar -czvf SPEC_OUT_corrected.tar.gz *.pdf SPEC_OUT_BROADENED/* smoothened/ specz/

fi
