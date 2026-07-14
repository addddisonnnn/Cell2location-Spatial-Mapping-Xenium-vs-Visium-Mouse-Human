"""
00_gpu_setup.py
---------------
GPU configuration for cell2location v0.3 (PyMC3/Theano backend).

MUST be imported before any other cell2location, theano, or pymc3 imports.
Import this at the very top of any script that uses cell2location.

Usage:
    import scripts.00_gpu_setup  # at the top of any script
    # or run standalone to verify GPU is detected:
    python scripts/00_gpu_setup.py
"""

import os

# Force Theano to use CUDA GPU.
# This must be set before theano is imported anywhere in the Python process.
os.environ["THEANO_FLAGS"] = "device=cuda,floatX=float32"

# ── Verify GPU is accessible ─────────────────────────────────────────────────
if __name__ == "__main__":
    import subprocess

    print("=" * 60)
    print("GPU CONFIGURATION CHECK")
    print("=" * 60)

    # Check nvidia-smi
    result = subprocess.run(["nvidia-smi"], capture_output=True, text=True)
    if result.returncode == 0:
        print("nvidia-smi output:")
        print(result.stdout)
    else:
        print("WARNING: nvidia-smi not found or no GPU available.")
        print("Cell2location will run on CPU — expect very slow training.")

    # Check theano device
    try:
        import theano
        print(f"Theano device: {theano.config.device}")
    except ImportError:
        print("Theano not importable — check Singularity container is active.")

    # Check torch
    try:
        import torch
        print(f"PyTorch CUDA available: {torch.cuda.is_available()}")
        if torch.cuda.is_available():
            print(f"GPU name: {torch.cuda.get_device_name(0)}")
    except ImportError:
        print("PyTorch not importable.")

    print("=" * 60)
    print("If GPU is shown above, you are ready to run the pipeline.")
    print("If no GPU is shown, request a V100 desktop session before proceeding.")
    print("=" * 60)