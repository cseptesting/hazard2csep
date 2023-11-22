# hazard2csep

## Example A

### Description

This example shows how to use the `hazard2csep` command line tool to convert a hazard file to the CSEP2 format. The example uses the Faults and Background Seismicity source model of the European Seismic Hazard Model 2013 [[EFEHR13 hazard repository](http://hazard.efehr.org/en/Documentation/specific-hazard-models/europe/overview/)]. 

### Usage

The following command will convert the OpenqQuake hazard file to the CSEP2 format. The output will be written to the current working directory. 

```shell
hazard2csep project --plot fs_bg_source_model.xml
```

or can be alternatively be run through the python script

```shell
python run.py
```
