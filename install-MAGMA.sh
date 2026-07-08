# Create directories
USER_ID="<USER>"
MAGMA_DIR="/myriadfs/home/${USER_ID}/SOFTWARE/MAGMA"

mkdir -p \
  "${MAGMA_DIR}/v1.10" \
  "${MAGMA_DIR}/reference/NCBI37.3" \
  "${MAGMA_DIR}/reference/g1000_eur" \
  "${MAGMA_DIR}/reference/synonyms" \
  "${MAGMA_DIR}/docs"

# Download MAGMA v1.10
cd "${MAGMA_DIR}/v1.10"
wget -O magma_v1.10.zip "https://vu.data.surf.nl/index.php/s/zkKbNeNOZAhFXZB/download"
unzip -o magma_v1.10.zip
chmod +x magma

# Download manual
cd "${MAGMA_DIR}/docs"
wget -O manual_v1.10.pdf "https://vu.data.surf.nl/index.php/s/MUiv3y1SFRePnyG/download"

# Download GRCh37 / Build 37 gene locations
cd "${MAGMA_DIR}/reference/NCBI37.3"
wget -O NCBI37.3.zip "https://vu.data.surf.nl/index.php/s/Pj2orwuF2JYyKxq/download"
unzip -o NCBI37.3.zip

# Download 1000 Genomes European reference data
cd "${MAGMA_DIR}/reference/g1000_eur"
wget -O g1000_eur.zip "https://vu.data.surf.nl/index.php/s/VZNByNwpD8qqINe/download"
unzip -o g1000_eur.zip

# Optional: download dbSNP synonym file
cd "${MAGMA_DIR}/reference/synonyms"
wget -O dbsnp151.synonyms.zip "https://vu.data.surfsara.nl/index.php/s/MSeFJuAVKJ4HLHv/download"
unzip -o dbsnp151.synonyms.zip

# Add MAGMA to PATH (note: allows you to execute easily - no need to specify location of executable)
echo "export PATH=\"${MAGMA_DIR}/v1.10:\$PATH\"" >> ~/.bashrc
source ~/.bashrc

# Check that runs
which magma
magma --version