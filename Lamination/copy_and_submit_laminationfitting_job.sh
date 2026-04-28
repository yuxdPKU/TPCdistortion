#!/bin/bash

TEMPLATE_FILE="condor-laminationfitting-template.job"

RUN_NUMBERS=(79516)

if [ ! -f "$TEMPLATE_FILE" ]; then
    echo "errors: template file $TEMPLATE_FILE not exist!"
    exit 1
fi

for RUN in "${RUN_NUMBERS[@]}"; do
    OUTPUT_FILE="condor-laminationfitting-${RUN}.job"

    cp "$TEMPLATE_FILE" "$OUTPUT_FILE"

    sed -i "s/RUNNUMBER/$RUN/g" "$OUTPUT_FILE"

    echo "Producing: $OUTPUT_FILE"

    condor_submit ${OUTPUT_FILE}
done

echo "All done!"
