#!/usr/bin/env bash
set -euo pipefail

study_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
raw_dir="${study_dir}/data/raw"
mkdir -p "${raw_dir}"

fetch() {
    local url="$1"
    local output="$2"
    local destination="${raw_dir}/${output}"
    local partial="${destination}.partial"

    if [[ -s "${destination}" ]]; then
        echo "already present: ${destination}"
        return
    fi

    curl --fail --location --retry 3 --retry-delay 2 \
        --output "${partial}" "${url}"
    mv "${partial}" "${destination}"
    echo "downloaded: ${destination}"
}

fetch "https://eprints.soton.ac.uk/462776/1/457354.pdf" \
    "insel-1990-thesis.pdf"
fetch "https://eprints.soton.ac.uk/46442/1/071.pdf" \
    "ship-science-71.pdf"
fetch "https://eprints.soton.ac.uk/46441/1/072.pdf" \
    "ship-science-72.pdf"
