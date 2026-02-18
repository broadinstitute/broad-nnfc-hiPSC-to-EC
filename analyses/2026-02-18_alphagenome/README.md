## Install alphagnenome on Sherlock from source

The major problem is that the default GCC version on Sherlock is too old to compile alphagenome. To get around this, we can install a newer version of GCC from conda-forge and make sure that it is the compiler used when installing alphagenome from source.

```bash
# Create a new conda environment
conda create -n alphagenome_conda -c conda-forge -c bioconda --solver=libmamba python=3.12 gcc=9 libmamba -y
# Activate the newly create environment
conda activate alphagenome_conda
# Install the latest GCC
conda install -c conda-forge gcc_linux-64 gxx_linux-64 cmake make
# Create new variables pointing to the compilers
export CC="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gcc"
export CXX="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++"
# Verify that $HOME/miniforge3/envs/alphagenome_conda/bin
# is the first entry in the PATH
echo $PATH | tr ':' '\n' | head
# If it is and you still get the error, install compilers from conda-forge and try again
conda install -c conda-forge compilers
hash -r
# Verify that points to the miniconda compilers
which g++
g++ --version


# At this point, you should be able to install alphagenome from source without
pip install -e ./alphagenome_research
``````