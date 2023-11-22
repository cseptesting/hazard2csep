# oq2pycsep


## Installation

```shell
conda env create -f environment.yml
pip install -e .
```

or

```shell
conda create -n hazard2csep
conda install -c conda-forge pycsep
conda install -c conda-forge openquake.engine
pip install -e .
```

## Download data

Download ESHM20

```shell
cd eshm_test
git clone https://gitlab.seismo.ethz.ch/efehr/eshm20 --depth=1
```

Download ESHM13
```shell
cd eshm_test
wget http://hazard.efehr.org/export/sites/efehr/.galleries/dwl_europe2013/SHARE_OQ_input_20140807.zip_2063069299.zip -O temp.zip && unzip temp.zip -d eshm13 && rm temp.zip 
unzip eshm13/SHARE_OQ_input_20140807/source_models.zip -d eshm13/SHARE_OQ_input_20140807/

```

## Get Regions

```shell
cd eshm_test
python get_regions.py
```


Region files will be output to `eshm_test/regions/`


