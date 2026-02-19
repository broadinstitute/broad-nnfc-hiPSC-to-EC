## Install alphagnenome from source on Sherlock

The easiest way to install alphagenome is via pip or conda.

```bash
# Install via pip
pip install -U alphagenome
# Or install via conda
conda install bioconda::alphagenome
```

```bash
If for some reason you need to install alphagenome from source on Sherlock, you can follow the instructions below.

$ pip install -U alphagenome


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

# Install other dependencies from binaries
conda install pandas numpy cython meson-python ninja

# At this point, you should be able to install alphagenome from source without
pip install --no-build-isolation -e ./alphagenome_research
```




## Extra notes

Depending on your setup, you may need to set additional environment variables to ensure that the correct compiler and libraries are used when installing alphagenome from source. For example, you may need to set the following variables:

```bash
export CC="$CONDA_PREFIX/bin/gcc"
export CXX="$CONDA_PREFIX/bin/g++"
export CFLAGS="-O2"
export CXXFLAGS="-O2 -std=c++17"
export LDFLAGS="-L$CONDA_PREFIX/lib"
export CPPFLAGS="-I$CONDA_PREFIX/include"
hash -r
```