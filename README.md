# create a causal environment
```bash
conda create -n causal python=3.8

conda activate causal
conda install ipykernel
python -m ipykernel install --user --name  causal
```
jupyterlab has the causal kernel, create a new jupyter file

# load zhangtao's package
```bash
git clone https://github.com/eexcalibur/my_utils.git
git clone https://github.com/eexcalibur/climate_causal.git
```
The example case is at climate_causal/ENSO_example/causal_ENSO_index_v1.ipynb, open it and select the causal kernel 



run the import part of the causal_ENSO_index_v1.ipynb, if there are loading errors, install the following packages
# install packages
```bash
conda install numpy
conda install pandas
conda install netcdf4
conda install scipy
conda install matplotlib
conda install cartopy
conda install scikit-learn
conda install statsmodels
conda install -c conda-forge cmocean
conda install -c conda-forge metpy
conda install -c conda-forge rpy2
conda install  -c conda-forge -c cdat/label/v8.2.1 cdat
conda install -c conda-forge python-igraph
conda install -c pytorch pytorch torchvision 
conda install -c conda-forge loguru
```


## install R and  rpy2

R lib is at /home/tzhang/soft/miniconda3/envs/causal3/lib/R

```r
R
install.packages("pcalg",repos='http://cran.us.r-project.org')
ERROR: dependencies ‘graph’, ‘RBGL’, ‘ggm’ are not available for package ‘pcalg’

install the error package using BiocManager
install.packages("BiocManager", repos='http://cran.us.r-project.org')
BiocManager::install('graph')
BiocManager::install('RBGL')
BiocManager::install('ggm')
install.packages("pcalg",repos='http://cran.us.r-project.org')
```

## test  causal package
```python
import os
os.environ['R_HOME'] = '/home/tzhang/soft/miniconda3/envs/causal3/lib/R'
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
pcalg = importr('pcalg')
sys.path.append("/home/tzhang/climate_causal/notears/notears/")
import linear
import utils
import nonlinear
```


