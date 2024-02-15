#! /bin/bash

INPUT_METABARCODING=""
INPUT_QPCR=""
SHOULD_RENAME="false"
METADATA_SCRIPT="$(pwd)/metadata_extractor.py"
METABARCODING_INITIALIZER="$(pwd)/eDNAExplorer_Metabarcoding_Metadata_Initializer.R"
SERVICE_CREDENTIALS="$(pwd)/google_earth_engine_service.json"
CONDASH=f"/usr/local/miniconda/etc/profile.d/conda.sh"
while getopts "1:2:p:c:" opt; do
    case $opt in
        1) INPUT_METABARCODING="$OPTARG"
        ;;
        2) INPUT_QPCR="$OPTARG"
        ;;
        p) PROJECTID="$OPTARG"
        ;;
    esac
done

JOBFOLDER="${PROJECTID}_$(openssl rand -hex 5)"

# Cleanup steps
cleanup() {
    rm -r "$JOBFOLDER"
}

# Set a trap to execute the cleanup function on EXIT
trap cleanup EXIT

# download input metadata 1
if [ -n "$INPUT_METABARCODING" ]; then
    mkdir -p $JOBFOLDER/tmp1
    aws s3 cp s3://$AWS_BUCKET/projects/$PROJECTID/$INPUT_METABARCODING $JOBFOLDER/ --no-progress

    #de-duplicate rows
    awk -F, 'BEGIN {OFS=","} NR == 1 {for (i=1; i<=NF; i++) col[$i] = i} !seen[$col["Sample ID"], $col["Latitude"], $col["Longitude"], $col["Sample Date"], $col["Spatial Uncertainty"]]++' "$JOBFOLDER/$INPUT_METABARCODING" > "$FOLDER/temp_file.csv"

    mv "$JOBFOLDER/temp_file.csv" "$JOBFOLDER/$INPUT_METABARCODING"

    header=$(head -n 1 $JOBFOLDER/$INPUT_METABARCODING)
    echo $header

    # Input 1 is metabarcoding file
    echo "$INPUT_METABARCODING is metabarcoding."
    OUTPUT_METADATA="MetadataOutput_Metabarcoding.csv"
    # split input file
    tail -n +2 $JOBFOLDER/$INPUT_METABARCODING | split -l 1 - $JOBFOLDER/tmp1/temp_
    for file in $JOBFOLDER/tmp1/temp_*; do
        if [[ ! "$file" =~ \.csv$ ]]; then
            (echo "$header"; cat "$file") > "${file}.csv"
        fi
    done

    # run metadata_extractor.py in parallel
    parallel -j 50 "output_file=\$(mktemp); echo \$output_file; cat {}; retries=3; while ((retries > 0)); do python3 $METADATA_SCRIPT --input {} --output \$output_file --account $SERVICE_ACCOUNT --credentials $SERVICE_CREDENTIALS && if [[ -s \$output_file ]]; then tail -n +2 \$output_file | cat >> $JOBFOLDER/output.txt; head -n 1 \"\$output_file\" | cat > $JOBFOLDER/headers.txt; break; else ((retries--)); fi; done; rm \$output_file" ::: $JOBFOLDER/tmp1/temp_*.csv


    # combine headers and outfile
    cat "$JOBFOLDER/headers.txt" "$JOBFOLDER/output.txt" > $JOBFOLDER/$OUTPUT_METADATA

    # upload output metadata
    aws s3 cp $JOBFOLDER/$OUTPUT_METADATA s3://$AWS_BUCKET/projects/$PROJECTID/$OUTPUT_METADATA --no-progress

    #cleanup
    cleanup

    # trigger initializer script
    sudo -u ubuntu bash -c "source $CONDASH && conda activate reports && Rscript --vanilla $METABARCODING_INITIALIZER $PROJECTID"
fi

# download input qPCR file
if [ -n "$INPUT_QPCR" ]; then
    mkdir -p $JOBFOLDER/tmp2
    aws s3 cp s3://$AWS_BUCKET/projects/$PROJECTID/$INPUT_QPCR $JOBFOLDER/ --no-progress

    #de-duplicate rows
    awk -F, 'BEGIN {OFS=","} NR == 1 {for (i=1; i<=NF; i++) col[$i] = i} !seen[$col["Sample ID"], $col["Latitude"], $col["Longitude"], $col["Sample Date"], $col["Spatial Uncertainty"]]++' "$JOBFOLDER/$INPUT_QPCR" > "$JOBFOLDER/temp_file.csv"

    mv "$JOBFOLDER/temp_file.csv" "$JOBFOLDER/$INPUT_QPCR"

    header=$(head -n 1 $JOBFOLDER/$INPUT_QPCR)
    echo $header

    # Input 2 is qPCR file
    echo "$INPUT_QPCR is qPCR."
    OUTPUT_METADATA="MetadataOutput_qPCR.csv"
    # split input file
    tail -n +2 $JOBFOLDER/$INPUT_QPCR | split -l 1 - $JOBFOLDER/tmp2/temp_
    for file in $JOBFOLDER/tmp2/temp_*; do
        if [[ ! "$file" =~ \.csv$ ]]; then
            (echo "$header"; cat "$file") > "${file}.csv"
        fi
    done

    # run metadata_extractor.py in parallel
    parallel -j 50 "output_file=\$(mktemp); echo \$output_file; cat {}; retries=3; while ((retries > 0)); do python3 $METADATA_SCRIPT --input {} --output \$output_file --account $SERVICE_ACCOUNT --credentials $SERVICE_CREDENTIALS && if [[ -s \$output_file ]]; then tail -n +2 \$output_file | cat >> $JOBFOLDER/output.txt; head -n 1 \"\$output_file\" | cat > $JOBFOLDER/headers.txt; break; else ((retries--)); fi; done; rm \$output_file" ::: $JOBFOLDER/tmp2/temp_*.csv

    # combine headers and outfile
    cat "$JOBFOLDER/headers.txt" "$JOBFOLDER/output.txt" > $JOBFOLDER/$OUTPUT_METADATA

    # upload output metadata
    aws s3 cp $JOBFOLDER/$OUTPUT_METADATA s3://$AWS_BUCKET/projects/$PROJECTID/$OUTPUT_METADATA --no-progress

    #cleanup
    cleanup
fi
