#!/bin/bash

# Download emolecules inventory file from Google Drive to revision/data/emolecules_inv.txt
# Requires: curl (pre-installed on macOS/Linux)

FILE_ID="1x33LmAizIdA5Dgw7IJp7_k7gRdNVv5cT"
DEST_DIR="$(dirname "$0")/data"
DEST_FILE="$DEST_DIR/emolecules_inv.txt"
COOKIE_FILE="$(mktemp)"

# Create destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

BASE_URL="https://drive.usercontent.google.com/download"

# Step 1: Fetch the page to get Google's virus-scan confirmation cookie
curl -c "$COOKIE_FILE" -s "${BASE_URL}?id=${FILE_ID}&export=download" > /dev/null

# Step 2: Download using the confirmation cookie (bypasses large-file warning)
curl -L -b "$COOKIE_FILE" \
  "${BASE_URL}?id=${FILE_ID}&export=download&confirm=t" \
  -o "$DEST_FILE"

rm -f "$COOKIE_FILE"
echo "Downloaded to $DEST_FILE"