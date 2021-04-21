# Install megaman package

- Software requirements: 
    
    Anaconda/Miniconda
    python versions 2.7, 3.5 or 3.6
    a C++ compiler such as gcc/g++
    
Following the instructions from the [source](https://github.com/mmp2/megaman), the megaman package can be installed on MacOS or Linux system by running the following conda commands in the terminal. 

```
$ conda create -n manifold_env python=3.5 -y
# can also use python=2.7 or python=3.6

$ conda activate manifold_env
$ conda install -channel=conda-forge -y pip nose coverage cython numpy scipy scikit-learn pyflann pyamg h5py plotly

$ cd /tmp/
$ git clone https://github.com/mmp2/megaman.git
$ cd megaman
$ python setup.py install
$ make test
```

# Use megaman package

```
$ conda activate manifold_env
$ pip install jupyter # install for the first time
$ jupyter notebook
```
