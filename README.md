# hazard2csep

Interface tool between [pyCSEP](https://github.com/SCECcode/pycsep) and [OpenQuake](https://github.com/gem/oq-engine) hazard models.


## Overview

This application tool reads the files of Source Model used in Seismic Hazard Analysis, in the OpenQuake format, and creates (pseudo)equivalent CSEP forecasts for testing purposes. The tool is based on the [pyCSEP](https://github.com/SCECcode/pycsep) and [OpenQuake](https://github.com/gem/oq-engine) libraries. 

The tool has the following main functionalities:

* Read, create and merge the spatial-magnitude region of the source model components (multiple).
* Intersect regions of different source models, so testing can be performed on a common region.
* Project the source model rates to the CSEP-style spatial-magnitude grid.

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

## Get Regions

Gets the closed (convex hull) region of the source model, intersected to a CSEP-style grid

```shell
cd examples/a_eshm13_faultareasources
hazard2csep region --plot fs_bg_source_model.xml
```

## Get Forecast

Projects the source model rates to a CSEP-style grid

```shell
cd examples/a_eshm13_faultareasources
hazard2csep project --plot fs_bg_source_model.xml
```

or alternatively, project the model onto a predetermined region:

```shell
cd examples/a_eshm13_faultareasources
hazard2csep project --plot -r region.txt fs_bg_source_model.xml
```



## Download related data

Download ESHM20 (Danciu et al., 2021)

```shell
git clone https://gitlab.seismo.ethz.ch/efehr/eshm20 --depth=1
```

Download ESHM13 (Woessner et al., 2015)
```shell
wget http://hazard.efehr.org/export/sites/efehr/.galleries/dwl_europe2013/SHARE_OQ_input_20140807.zip_2063069299.zip -O temp.zip && unzip temp.zip -d eshm13 && rm temp.zip 
unzip eshm13/SHARE_OQ_input_20140807/source_models.zip -d eshm13/SHARE_OQ_input_20140807/

```


## References and Links

* pyCSEP: https://github.com/SCECcode/pycsep


Savran, W. H., Bayona, J. A., Iturrieta, P., Asim, K. M., Bao, H., Bayliss, K., Herrmann, M., Schorlemmer, D., Maechling, P. J. & Werner, M. J. (2022). pyCSEP: A python toolkit for earthquake forecast developers. Seismological Research Letters, 93(5), 2858-2870. https://doi.org/10.1785/0220220033

* OpenQuake: https://github.com/gem/oq-engine

Pagani, M., Monelli, D., Weatherill, G., Danciu, L., Crowley, H., Silva, V., ... & Vigano, D. (2014). OpenQuake engine: An open hazard (and risk) software for the global earthquake model. Seismological Research Letters, 85(3), 692-702.


* The 2013 European Seismic Hazard Model

Woessner, J., Danciu L., D. Giardini and the SHARE consortium (2015), The 2013 European Seismic Hazard Model: key components and results, Bull. Earthq. Eng., doi:10.1007/s10518-015-9795-1.

* The 2020 European Seismic Hazard Model:

Danciu L., Nandan S., Reyes C., Basili R., Weatherill G., Beauval C., Rovida A., Vilanova S., Sesetyan K., Bard P-Y., Cotton F., Wiemer S., Giardini D. (2021) - The 2020 update of the European Seismic Hazard Model: Model Overview. EFEHR Technical Report 001, v1.0.0, https://doi.org/10.12686/a15

