import cell2location
import scanpy as sc

# These are the exact files used in the paper tutorials
# They download automatically to a local cache

# snRNA-seq reference
adata_snrna = sc.read(
    filename="snRNA_mousebrainSection1.h5ad",
    backup_url="https://cell2location.cog.sanger.ac.uk/paper/mouse_brain_snrna_seq/snRNA_mousebrainSection1.h5ad"
)

# Visium spatial
adata_vis = sc.read(
    filename="visium_mousebrain.h5ad",  
    backup_url="https://cell2location.cog.sanger.ac.uk/paper/mouse_brain_visium/visium_mousebrain.h5ad"
)

# Check what you have
print(adata_snrna)
print(adata_snrna.obs["annotation_1_print"].value_counts())  # cell type labels

