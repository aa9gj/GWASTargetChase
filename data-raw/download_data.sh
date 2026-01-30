#!/usr/bin/env bash
#
# download_data.sh
# ================
# Downloads real data from IMPC, OpenTargets Platform, OpenTargets Genetics,
# and GENCODE for use with GWASTargetChase.
#
# Usage:
#   chmod +x data-raw/download_data.sh
#   cd data-raw
#   bash download_data.sh
#
# After downloading, use prepare_real_data.R to process and save as .rda files.
#
# NOTE: These downloads are large (several GB total). Ensure sufficient disk space.

set -euo pipefail

DOWNLOAD_DIR="downloaded_data"
mkdir -p "$DOWNLOAD_DIR"

echo "=============================================="
echo "  GWASTargetChase - Data Download Script"
echo "=============================================="
echo ""
echo "Download directory: $DOWNLOAD_DIR"
echo ""

# ------------------------------------------------------------------------------
# 1. IMPC Data (International Mouse Phenotyping Consortium)
# ------------------------------------------------------------------------------
echo "----------------------------------------------"
echo "1. Downloading IMPC phenotype hits per gene..."
echo "----------------------------------------------"
# Check https://www.mousephenotype.org/data/release for the latest release
IMPC_URL="https://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/phenotypeHitsPerGene.csv.gz"
if [ ! -f "$DOWNLOAD_DIR/phenotypeHitsPerGene.csv" ]; then
    wget -P "$DOWNLOAD_DIR" "$IMPC_URL" || \
        curl -L -o "$DOWNLOAD_DIR/phenotypeHitsPerGene.csv.gz" "$IMPC_URL"
    gunzip "$DOWNLOAD_DIR/phenotypeHitsPerGene.csv.gz"
    echo "IMPC data downloaded and extracted."
else
    echo "IMPC data already exists, skipping."
fi
echo ""

# ------------------------------------------------------------------------------
# 2. OpenTargets Platform - Genetic Association Evidence + Disease Annotations
# ------------------------------------------------------------------------------
# OpenTargets Platform has moved to Google Cloud Storage for recent releases.
# Check https://platform.opentargets.org/downloads for the latest release version.
# Update the version number below as needed.
OT_PLATFORM_VERSION="24.06"

echo "----------------------------------------------"
echo "2a. Downloading OpenTargets Platform evidence"
echo "    (associationByDatasourceDirect)..."
echo "    Version: $OT_PLATFORM_VERSION"
echo "----------------------------------------------"
OT_ASSOC_URL="http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/${OT_PLATFORM_VERSION}/output/etl/json/associationByDatasourceDirect"
ASSOC_DIR="$DOWNLOAD_DIR/associationByDatasourceDirect"
if [ ! -d "$ASSOC_DIR" ]; then
    mkdir -p "$ASSOC_DIR"
    wget -r -np -nH --cut-dirs=8 -P "$ASSOC_DIR" "$OT_ASSOC_URL/" 2>&1 || {
        echo "wget recursive download failed. Trying alternative approach..."
        echo "You may need to download manually from:"
        echo "  $OT_ASSOC_URL"
        echo "Or use gsutil:"
        echo "  gsutil -m cp -r gs://open-targets-data-releases/${OT_PLATFORM_VERSION}/output/etl/json/associationByDatasourceDirect $ASSOC_DIR"
    }
    # Concatenate all JSON parts into a single file
    if ls "$ASSOC_DIR"/*.json 1>/dev/null 2>&1; then
        cat "$ASSOC_DIR"/*.json > "$DOWNLOAD_DIR/associationByDatasourceDirect.json"
        echo "Association data downloaded and concatenated."
    else
        echo "WARNING: No JSON files found. Check download."
    fi
else
    echo "Association data directory already exists, skipping."
fi
echo ""

echo "----------------------------------------------"
echo "2b. Downloading OpenTargets Platform diseases..."
echo "----------------------------------------------"
OT_DISEASE_URL="http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/${OT_PLATFORM_VERSION}/output/etl/json/diseases"
DISEASE_DIR="$DOWNLOAD_DIR/diseases"
if [ ! -d "$DISEASE_DIR" ]; then
    mkdir -p "$DISEASE_DIR"
    wget -r -np -nH --cut-dirs=8 -P "$DISEASE_DIR" "$OT_DISEASE_URL/" 2>&1 || {
        echo "wget recursive download failed. Trying alternative approach..."
        echo "You may need to download manually from:"
        echo "  $OT_DISEASE_URL"
        echo "Or use gsutil:"
        echo "  gsutil -m cp -r gs://open-targets-data-releases/${OT_PLATFORM_VERSION}/output/etl/json/diseases $DISEASE_DIR"
    }
    # Concatenate all JSON parts
    if ls "$DISEASE_DIR"/*.json 1>/dev/null 2>&1; then
        cat "$DISEASE_DIR"/*.json > "$DOWNLOAD_DIR/diseases.json"
        echo "Disease data downloaded and concatenated."
    else
        echo "WARNING: No JSON files found. Check download."
    fi
else
    echo "Disease data directory already exists, skipping."
fi
echo ""

# ------------------------------------------------------------------------------
# 3. OpenTargets Genetics - Locus-to-Gene (L2G) + Study Index
# ------------------------------------------------------------------------------
# OpenTargets Genetics data is available via FTP/Google Cloud Storage.
# The L2G data is in Parquet format (read by the arrow R package).
# Check https://genetics-docs.opentargets.org/data-access for latest info.
OT_GENETICS_VERSION="latest"

echo "----------------------------------------------"
echo "3a. Downloading OpenTargets Genetics L2G data..."
echo "----------------------------------------------"
L2G_URL="http://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/${OT_GENETICS_VERSION}/l2g"
L2G_DIR="$DOWNLOAD_DIR/l2g"
if [ ! -d "$L2G_DIR" ]; then
    mkdir -p "$L2G_DIR"
    wget -r -np -nH --cut-dirs=6 -P "$L2G_DIR" "$L2G_URL/" 2>&1 || {
        echo "wget recursive download failed. Try:"
        echo "  gsutil -m cp -r gs://genetics-portal-data/${OT_GENETICS_VERSION}/l2g/ $L2G_DIR"
    }
    echo "L2G data downloaded."
else
    echo "L2G data directory already exists, skipping."
fi
echo ""

echo "----------------------------------------------"
echo "3b. Downloading OpenTargets Genetics study index..."
echo "----------------------------------------------"
STUDY_URL="http://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/${OT_GENETICS_VERSION}/lut/study-index"
STUDY_DIR="$DOWNLOAD_DIR/study-index"
if [ ! -d "$STUDY_DIR" ]; then
    mkdir -p "$STUDY_DIR"
    wget -r -np -nH --cut-dirs=7 -P "$STUDY_DIR" "$STUDY_URL/" 2>&1 || {
        echo "wget recursive download failed. Try:"
        echo "  gsutil -m cp -r gs://genetics-portal-data/${OT_GENETICS_VERSION}/lut/study-index/ $STUDY_DIR"
    }
    # Concatenate study index JSON files
    if ls "$STUDY_DIR"/*.json 1>/dev/null 2>&1; then
        cat "$STUDY_DIR"/*.json > "$DOWNLOAD_DIR/study_index.json"
        echo "Study index downloaded and concatenated."
    else
        echo "WARNING: No JSON files found. Check download."
    fi
else
    echo "Study index directory already exists, skipping."
fi
echo ""

# ------------------------------------------------------------------------------
# 4. GENCODE Human GTF (needed for geneticAssocPrep and l2gPrep)
# ------------------------------------------------------------------------------
echo "----------------------------------------------"
echo "4. Downloading GENCODE human GTF file..."
echo "----------------------------------------------"
GENCODE_VERSION="44"
GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.annotation.gtf.gz"
if [ ! -f "$DOWNLOAD_DIR/gencode.v${GENCODE_VERSION}.annotation.gtf" ]; then
    wget -P "$DOWNLOAD_DIR" "$GENCODE_URL" || \
        curl -L -o "$DOWNLOAD_DIR/gencode.v${GENCODE_VERSION}.annotation.gtf.gz" "$GENCODE_URL"
    gunzip "$DOWNLOAD_DIR/gencode.v${GENCODE_VERSION}.annotation.gtf.gz"
    echo "GENCODE GTF downloaded and extracted."
else
    echo "GENCODE GTF already exists, skipping."
fi
echo ""

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
echo "=============================================="
echo "  Download Summary"
echo "=============================================="
echo ""
echo "Files in $DOWNLOAD_DIR:"
ls -lh "$DOWNLOAD_DIR/" 2>/dev/null || echo "(directory listing failed)"
echo ""
echo "Next steps:"
echo "  1. Open R and run: source('data-raw/prepare_real_data.R')"
echo "  2. This will process the downloads and save .rda files to data/"
echo "  3. Rebuild the package to include the new data"
echo ""
echo "Done!"
