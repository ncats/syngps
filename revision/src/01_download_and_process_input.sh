#!/bin/bash


# Downloads buyables data used in ASKCOS
# File modified from https://gitlab.com/mlpds_mit/askcosv2/askcos2_core/-/blob/main/scripts/download_data.sh?ref_type=heads#L29

# Set up destination directory relative to this script's location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGET_DIR="$SCRIPT_DIR/../data/input"
BUYABLES_DIR="$TARGET_DIR/asckos_buyables_raw"
mkdir -p "$BUYABLES_DIR"

if [ ! -f "$BUYABLES_DIR/buyables.json.gz" ]; then
    echo "$BUYABLES_DIR/buyables.json.gz not found. Downloading.."
    curl -L --progress-bar -o "$BUYABLES_DIR/buyables.json.gz" \
      "https://www.dropbox.com/scl/fi/jaqo5r8jji7n19ijrohj1/buyables_ori_prop_new2.json.gz?rlkey=lnbcj7t1ygrjgqi3g7a7vshz2&st=krp61fnf&dl=1"
    echo "buyables.json.gz Downloaded."
fi

if [ ! -f "$BUYABLES_DIR/chembridge_buyables.json.gz" ]; then
    echo "$BUYABLES_DIR/chembridge_buyables.json.gz not found. Downloading.."
    curl -L --progress-bar -o "$BUYABLES_DIR/chembridge_buyables.json.gz" \
      "https://www.dropbox.com/scl/fi/v5jayotulkdueiebck39g/chembridge_buyables_prop_new.json.gz?rlkey=zj3ncebc29ohzmhenj8fmeeo2&dl=1"
    echo "chembridge_buyables.json.gz Downloaded."
fi

if [ ! -f "$BUYABLES_DIR/mcule_buyables_fd2.json.gz" ]; then
    echo "$BUYABLES_DIR/mcule_buyables_fd2.json.gz not found. Downloading.."
    curl -L --progress-bar -o "$BUYABLES_DIR/mcule_buyables_fd2.json.gz" \
      "https://www.dropbox.com/scl/fi/1ikitf6n1lylyq3eroqxp/mcule_buyables_fd2_prop_new.json.gz?rlkey=2rjslc4n7gvvplsznx05vom4i&dl=1"
    echo "mcule_buyables_fd2.json.gz Downloaded."
fi

if [ ! -f "$BUYABLES_DIR/chemspace_buyables_2026Apr.json.gz" ]; then
    echo "$BUYABLES_DIR/chemspace_buyables_2026Apr.json.gz not found. Downloading.."
    curl -L --progress-bar -o "$BUYABLES_DIR/chemspace_buyables_2026Apr.json.gz" \
      "https://www.dropbox.com/scl/fi/iau2q04rt7woiof70up5e/chemspace_buyables_2026Apr.json.gz?rlkey=jcrzv1s3o2dpfrtdzljlm0ngz&st=8rhguy74&dl=1"
    echo "chemspace_buyables_2026Apr.json.gz Downloaded."
fi

# Downloads SimpRetro synthesis targets

if [ ! -f "$TARGET_DIR/retrostar_raw_test.csv" ]; then
    echo "$TARGET_DIR/retrostar_raw_test.csv not found. Downloading.."
    curl -L --progress-bar -o "$TARGET_DIR/retrostar_raw_test.csv" \
      "https://www.dropbox.com/scl/fo/ef2688z253mc5hq9oki4s/AHHjWUKhnCk1_Y_Oc1VDC6I/raw_test.csv?rlkey=y06bcu97k2x0zo31s7pd2u2bc&dl=0"
    echo "retrostar_raw_test.csv Downloaded."
fi


# Downloads RetroStar synthesis targets

if [ ! -f "$TARGET_DIR/SMILES.txt" ]; then
    echo "${TARGET_DIR}/SMILES.txt not found. Downloading.."
    curl -L --progress-bar -o "$TARGET_DIR/SMILES.txt" \
      "https://raw.githubusercontent.com/catalystforyou/SimpRetro/master/SMILES.txt"
    echo "SMILES.txt Downloaded."
fi

# Process ASKCOS buyables into a single inventory TSV
echo "Processing ASKCOS buyables..."
conda run -n syngps_rev python "$SCRIPT_DIR/process_askcos_inv.py"
echo "ASKCOS inventory processing complete."

