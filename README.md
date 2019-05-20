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


## Target Photometry

You can get the photometry assocated to a given target using the `PhotoTarget` class inside `pymage.phototarget.py`.
This module relies on `pymage.query.py` for accessing  the data associated to the target of interest. 

Let's consider SN2005ir ra=19.18233, dec=0.79456, zcmb=0.07517:
```python
from pymage import phototarget
pt = phototarget.PhotoTarget(ra=19.18233, dec=0.79456, zcmb=0.07517, name="SN2005ir")
# Let's load sdss images and attach them to the current object
pt.load_instrument("sdss")
# Let's do the same for galex (which is slightly slower)
pt.load_instrument("galex")
# Instruments (astrobject's Instrument objects) are stored as pt.instruments
pt.instruments
"""
{'galex': {'fuv': [<astrobject.instruments.galex.GALEX at 0x11314a358>,
   <astrobject.instruments.galex.GALEX at 0x11324b0f0>,
   <astrobject.instruments.galex.GALEX at 0x116cc27f0>],
  'nuv': [<astrobject.instruments.galex.GALEX at 0x113148ef0>,
   <astrobject.instruments.galex.GALEX at 0x11315e7b8>,
   <astrobject.instruments.galex.GALEX at 0x116cecfd0>]},
 'sdss': {'sdssg': [<astrobject.instruments.sdss.SDSS at 0x117408390>],
  'sdssi': [<astrobject.instruments.sdss.SDSS at 0x2356cdf60>],
  'sdssr': [<astrobject.instruments.sdss.SDSS at 0x117b04978>],
  'sdssu': [<astrobject.instruments.sdss.SDSS at 0x1173a7fd0>],
  'sdssz': [<astrobject.instruments.sdss.SDSS at 0x1173b2198>]}}
"""
```
If you now want the local target photometry (within say 3 'kpc' or 2 'arcsec', you can get a photopoints associated to the loaded instruments (astrobject's object containing the magnitude, flux, error etc).

For instance, let's get the photopoint for the sdss 'g' band:
```python
sdssg_ppoint = pt.get_photopoints("sdss.sdssg", 3, 'arcsec')[0]
print(sdssg_ppoint.flux, sdssg_ppoint.mag, sdssg_ppoint.mag_err, sdssg_ppoint.magabs)
"""
1.617e-16 18.7031 (0.0096, 0.00956) -19.025
"""
```
Remark that you need `[0]` because the code could have returned several photopoints if your target had several images at the given band. Also, because the object has a redshift, the code can do the mag<->mag_abs conversion. 

You could also ask the local photometry in "kpc" (in that case galex nuv band. there are several images for SN2005ir):

```python
nuv_ppoint = pt.get_photopoints("galex.nuv", 4, 'kpc')
print([[p_.magabs, p_.mag_err] for p_ in nuv_ppoint])
"""
[[-16.831, (0.032, 0.032)],
 [-16.835, (0.125, 0.109)],
 [-16.913, (0.151, 0.136)]]
"""
```
