#!/bin/bash

get_closest_numbers() {
    number=$1
    lower=$(( (number / 100) * 100 ))
    upper=$(( lower + 100 ))
    echo "$lower $upper"
}

TEMPLATE_FILE="condor-data-seed-template.job"

#RUN_NUMBERS=(53018 53046 53079 53080 53081 53194 53195 53196 53197 53494 53513 53517 53530 53531 53532 53571 53578 53579 53580 53581 53586 53587 53590 53686)
RUN_NUMBERS=(53877)

if [ ! -f "$TEMPLATE_FILE" ]; then
    echo "errors: template file $TEMPLATE_FILE not exist!"
    exit 1
fi

for RUN in "${RUN_NUMBERS[@]}"; do

    OUTPUT_FILE="condor-data-seed-${RUN}.job"

    cp "$TEMPLATE_FILE" "$OUTPUT_FILE"

    sed -i "s/RUNNUMBER/$RUN/g" "$OUTPUT_FILE"

    result=$(get_closest_numbers ${RUN})
    result_arr=($result)
    numlower=${result_arr[0]}
    numupper=${result_arr[1]}
    runrange=run_000${numlower}_000${numupper}

    sed -i "s/RUNRANGE/${runrange}/g" "$OUTPUT_FILE"

    echo "Producing: $OUTPUT_FILE"

    condor_submit ${OUTPUT_FILE}
done

echo "All done!"
