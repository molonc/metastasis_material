#!/usr/bin/env bash
# Find all files ending with "metrics.csv" under the current directory
# and compress each into a corresponding "<name>metrics.csv.gz".
set -euo pipefail

find . -type f -name '*metrics.csv' -print0 |
while IFS= read -r -d '' file; do
    echo "Compressing: $file"
    gzip -f "$file"
done

echo "Done."
