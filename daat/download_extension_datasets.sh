# Dataset 2 - Human Glioblastoma Visium
cd "/projectnb/ds596/projects/Team_7/data"

# Create a folder for it
mkdir humanGlioblastoma
cd humanGlioblastoma

# Download the filtered count matrix
wget https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5

# Download the spatial folder
wget https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma_spatial.tar.gz

# Extract spatial folder
tar -xzf Parent_Visium_Human_Glioblastoma_spatial.tar.gz

# Confirm what you have
ls -lh


# Dataset 3 - Xenium Mouse Brain
cd "/projectnb/ds596/projects/Team_7/data"

mkdir xeniumMouseBrain
cd xeniumMouseBrain

# Start with the SMALL subset for parameter testing — much faster to download and run
wget https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip

# Unzip
unzip Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip

# Confirm contents
ls -lh outs/
