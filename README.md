# pymage
Python module to Measure AGE bias in SN Cosmology

# Installation
```
cd WHERE YOU STORE CODES (like ~/Libraries)
git clone https://github.com/MickaelRigault/pymage.git 
```
then 
```
cd pymage
python setup.py install
```

# Usage

## Data query
First define the global variable "DATAPATH" (in your .bash_profile)
Let's image you want to download galex data for "SNF20070712-003" which is at ra=250.667625, dec=50.5975166667:

```python
from pymage import query
gquery = query.GALEXQuery()
# Grab the metadata from MAST archive:
gquery.download_target_metadata("SNF20070712-003", 250.667625, 50.5975166667, dl=False)
```
this updates the metadata attribute:
```python
print(gquery.metadata)
```
```
	  name	        ra	      dec	  filters	project	dataURL	basename	baseurl
0	SNF20070712-003	250.667625	50.597517	FUV	AIS	NaN	AIS_6_1_15	http://galex.stsci.edu/data/GR6/pipe/02-vsn/50...
1	SNF20070712-003	250.667625	50.597517	NUV	AIS	NaN	AIS_6_1_15	http://galex.stsci.edu/data/GR6/pipe/02-vsn/50...
```
Now to download the actual galex images on your computer, simply do:
```python
gquery.download_target_data("SNF20070712-003")
```

Once this is done, to retreive the fullpath of the local existing data:
```python 
data_fullpath = gquery.get_target_data("SNF20070712-003")
```

### Quick versions

If you have downloaded the data on your computer. You can directly retreive the data:
```python
from pymage import query
gquery = query.GALEXQuery()
data_fullpath = gquery.get_target_data("SNF20070712-003")
```

You can also download images right after downloading the metadata by setting the `dl` option to true:
```python
from pymage import query
gquery = query.GALEXQuery()
gquery.download_target_metadata("SNF20070712-003", 250.667625, 50.5975166667, dl=False)
data_fullpath = gquery.get_target_data("SNF20070712-003")
```

