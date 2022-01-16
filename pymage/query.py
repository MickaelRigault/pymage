#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pandas
import numpy as np

from astroquery.mast import Observations

from . import io

KNOWN_INSTRUMENTS  = ["galex", "sdss", "panstarrs"]

def query_metadata(instrument, ra, dec):
    """ query metadata for the given intrument and coordinates
    the metadata are information about filename and url
    """
    return eval("query_%s_metadata(ra, dec)"%_test_instrument_(instrument))

def load_metadata(instrument, load_empty=False):
    """ Load existing metadata for the given instrument 
    Returns
    -------
    Pandas.DataFrame 
    //(columns are (at least) "name","ra","dec","filters", "project", "basename", "baseurl")
    """
    empty = pandas.DataFrame(columns=["name","ra","dec","filters", "project",
                                      "basename", "baseurl"])
    if load_empty:
        return empty
    
    metadata_file = _get_metadata_file_(_test_instrument_(instrument))
    return pandas.read_csv(metadata_file) if os.path.exists(metadata_file) else empty    

def get_directory(instrument):
    """ The default instrument directory on your computer """
    return  io.DATAPATH+"%s/"%instrument.upper()

def metadata_to_url(instrument, baseurl, basename, bands, **kwargs):
    """ build url (or local fullpath) from metadata information """
    return eval("_%s_info_to_urlpath_(baseurl=baseurl, basename=basename, bands=bands,**kwargs)"%_test_instrument_(instrument) )

def _get_metadata_file_(instrument):
    """ returns the metadata associated to the given instruement """
    return eval("_%s_METADATA_FILE"%instrument.upper())

def _test_instrument_(instrument):
    """ """
    if instrument.lower() not in KNOWN_INSTRUMENTS:
        raise NotImplementedError("unknown instrument: %s"%instrument +"\n"+"known instruments: "+",".join(KNOWN_INSTRUMENTS))
    return instrument.lower()

##########################
#                        #
#    GENERAL TOOLS       #
#                        #
##########################
class _Query_( object ):
    """ Simply Class to manage data IO """
    INSTRUMENT = "to_be_defined"
    
    def __init__(self, empty=False):
        """ initialize query. It loads the known `instrument` meta data. """
        self.metadata = load_metadata(self.INSTRUMENT, load_empty=empty)

    # --------- #
    #  GETTER   #
    # --------- #
    def get_target_data(self, targetname, must_exists=True, fromdir="default",
                            filters="*", **kwargs):
        """ returns the full path of data on your computer. 
        
        """
        if fromdir is None or fromdir in ["default"]:
            fromdir = self._default_dldir

        urls, localpaths = self._build_target_data_url_and_path_(targetname, fromdir,
                                                                     filters=filters, **kwargs)
        return [l_ for l_ in localpaths if os.path.exists(l_) or not must_exists]

    def get_target_coords(self, targetname):
        """ """
        if not self.is_target_known(targetname):
            raise AttributeError("unknown target. Please run download_target_metadata()")
        # This assumes all entry at the same name have the same coordinate as it should.
        return np.asarray(self.metadata[self.metadata["name"]==targetname
                                        ].iloc[0].get(["ra","dec"]).values, dtype="float")

    def get_target(self, targetname):
        """ Loads an astrobject target (name, ra, dec) and returns it """
        from astrobject import get_target
        ra,dec = self.get_target_coords(targetname)
        return get_target(name=targetname, ra=ra, dec=dec)

    
    def get_target_instruments(self, targetname, cachedl=False, filters="*"):
        """ Return a list of Astrobject's Instrument for each entry coresponding to the given target """
        if not self.is_target_known(targetname):
            raise AttributeError("unknown target. Please run download_target_metadata(), and download the associated files")

        from astrobject import instruments
        target_data = self.get_target(targetname)
        
        # Cache Download
        load_prop = dict(target=target_data, instrument=self.INSTRUMENT.lower())
        if cachedl:
            sourcefile = self.download_target_data(targetname, "default", filters=filters, dl=False)[0]
            load_prop["cache"]=False
        else:
            sourcefile = self.get_target_data(targetname,filters=filters)
            
        # load files
        return [instruments.get_instrument(f_, **load_prop)
                for f_ in sourcefile]
    
    
    # ------------- #
    #  Downloader   #
    # ------------- #
    def query_metadata(self, ra, dec):
        """ """
        df_ = query_metadata(self.INSTRUMENT, ra, dec)
        df_["ra"]   = ra
        df_["dec"]  = dec
        return df_
    
    def download_target_metadata(self, targetname, ra, dec, 
                                 update=True, store=True, 
                                 overwrite=False, dl=True):
        """ Look for metadata online archive, save them and download the corresponding  files. 
        (options enables to store or not and download or not)

        Parameters:
        -----------
        targetname: [str]
            instrument metadata will be associated to this name.

        ra,dec: [float, float] // in degree
            Coordinates, in degree, of the target.

        // options //
        update: [bool] -optional-
            Shall the downloaded metadata be inserted in the object's self.metadata ?
            => If False, the downloaded metadata is returned

        store: [bool] -optional-
            If the downloaded metadata has been inserted to the object's self.metadata
            shall this file stored in your computer be updated too ? [you should !]

        overwrite: [bool] -optional-
            If the object's metadata already contains the target, shall we overwrite it ?
            If not, this entire function is skipped 
            (see self.download_target_data(targetname) do download data of already known targets)

        dl: [bool] -optional-
            Shall the data associated to the target's metadata be downloaded ?

        Returns
        -------
        None (or DataFrame if update=False)
        """
        if self.is_target_known(targetname) and not overwrite:
            print("no need")
            return
        
        df_ = self.query_metadata(ra,dec)
        df_["name"] = targetname
        
        if not update:
            return df_
        
        # merge with current
        if self.is_target_known(targetname):
            self.metadata.drop(index=self.metadata.index[self.metadata["name"]==targetname], inplace=True)
            
        self.metadata = pandas.concat([self.metadata, df_], sort=False)
        if store:
            self.store()
            
        # Downloading
        if dl:
            self.download_target_data(targetname)
            
        return 0 # 0 means no problem

    def store(self):
        """ """
        fileout = _get_metadata_file_(self.INSTRUMENT)
        if not os.path.isdir(os.path.dirname(fileout)):
            if io.DATAPATH == "_notdefined_":
                raise AttributeError("You must define the global variable DATAPATH to bve able to download/store data")
                
            os.mkdir(os.path.dirname(fileout))
                
        self.metadata.to_csv(fileout, index=False)
        
    def download_target_data(self, targetname, dirout="default",
                                 overwrite=False,  dl=True, verbose=True, **kwargs):
        """ Download the target photometric data. 
        
        Parameters
        ----------
        targetname: 
            Name of a target known by the class

        dirout: [string] -optional-
            Where shall the data be downloaded
            - "default": the default local structure | use this if unsure
            - "PATH": provide any path, the data will be downloaded here
            - "StringIO": download the data inside StringIO files [they will be returned]

        overwrite: [bool] -optional-
            If the file already exists where you want to download them, should this overwrite them?
            
        dl: [bool] -optional-
            Should the download be actually launched ?
            If False: the returns the urls to be downloaded and where they will be.

        **kwargs goes to _build_target_data_url_and_path_

        Returns
        -------
        None (or list of StringIO if dirout='StringIO')
    
        """
        if dirout is None or dirout in ["default"]:
            dirout = self._default_dldir
            
        if io.DATAPATH == "_notdefined_":
            raise AttributeError("You must define the global variable DATAPATH to bve able to download/store data")

        if dirout in ["StringIO", "stringio", "iostring", "io", "BytesIO","BytesIO","bytes"]:
            urls = self._build_target_data_url_and_path_(targetname, "default", **kwargs)[0]
            # Bytes IO are more suitable for internet requests
            localpaths = ["BytesIO" for i in range(len(urls))]
            is_stringio=True
            overwrite=True
        else:
            urls, localpaths = self._build_target_data_url_and_path_(targetname, dirout, **kwargs)
            is_stringio=False

        if not dl:
            return urls, localpaths
        
        return [io.download_single_url(url_, localpath_,
                    overwrite=overwrite, verbose=verbose)
                    for url_, localpath_ in zip(urls, localpaths)]
        
    def _build_target_data_url_and_path_(self, targetname, dirout, filters=None, **kwargs):
        """ Returns the URL and download location of data """
        if not self.is_target_known(targetname):
            raise AttributeError("unknown target. Please run download_target_metadata()")

        url_ = np.asarray([metadata_to_url(self.INSTRUMENT, row["baseurl"], row["basename"], bands=row["filters"], **kwargs)
               for index_, row in self.metadata[self.metadata["name"]==targetname].iterrows()
                               if (filters is None or filters in ["all","*"]) or row["filters"] in filters]).flatten()
        localpath_ = np.asarray([metadata_to_url(self.INSTRUMENT, dirout, row["basename"], bands=row["filters"], **kwargs)
                for index_, row in self.metadata[self.metadata["name"]==targetname].iterrows()
                               if (filters is None or filters in ["all","*"]) or row["filters"] in filters]).flatten()
        
        return url_, localpath_
        
    def is_target_known(self, targetname):
        """ Test if the given target has known metadata. """
        return targetname in self.known_targets

    ##################
    #   Properties   #
    ##################
    @property
    def known_targets(self):
        """ list of targets inside metadata """
        return self.metadata["name"].values

    @property
    def _default_dldir(self):
        """ """
        return get_directory(self.INSTRUMENT)

# ====================== #
#                        #
#    GALEX               #
#                        #
# ====================== #
GALEX_DIR            = os.path.join(io.DATAPATH, "GALEX")
_GALEX_METADATA_FILE = os.path.join(GALEX_DIR,"target_source.csv")

def query_mast(ra, dec, instrument=None, radius="10 arcsec"):
    """ returns a table containing the existing data information """
    from astroquery.mast import Observations
    t = Observations.query_region("%.5f %+.5f"%(ra,dec), radius=radius)
    if instrument is not None:
        return t[t["obs_collection"]==instrument]
    
    return t

def _galex_info_to_urlpath_(baseurl, basename, which=["int", "skybg"], bands=["NUV","FUV"]):
    """ Build URL or fullpath for the given data """
    
    return [[baseurl+"/"+basename+"-%sd-%s.fits.gz"%(band[0].lower(),todl_) for band in np.atleast_1d(bands)]
           for todl_ in np.atleast_1d(which)]
    
def query_galex_metadata(ra, dec):
    """ look for GALEX meta data inside MAST archive """        
    t = query_mast(ra, dec, instrument="GALEX")
    df = pandas.DataFrame(dict(t[["filters", "project","target_name"]]))
    df["basename"] = [t_.replace("_1_","_sg") if t_.startswith("AIS") else t_ for t_ in df.pop("target_name")]
    df["basename"] = [b_.replace("_sg","_sg0") if "_sg" in b_ and len(b_.split("_sg")[-1])==1 else b_ for b_ in df["basename"]]
    df = df.assign(baseurl= ["/".join(url_.split("/")[:-1]) for url_ in t["dataURL"].data])
    return df

#
# CLASS 
#
class GALEXQuery( _Query_ ):
    """ Simply Class to manage the GALEX data IO """
    INSTRUMENT = "GALEX"
    def get_target_instruments(self, targetname, contains=None, buffer_safe_width=0.05):
        """ """
        if not self.is_target_known(targetname):
            raise AttributeError("unknown target. Please run download_target_metadata(), and download the associated files")

        from astrobject import instruments
        # Which data to use
        target_data = self.get_target_data(targetname)
        all_data_int    = [f for f in target_data if "int" in f and (contains is None or contains in f)]
        # Which data to use
        target = self.get_target(targetname)
        instru = []
        for fullpath in [f for f in all_data_int if f.replace("int","skybg") in target_data]:
            inst_ = instruments.get_instrument(fullpath)
            if not inst_.is_target_in(target,buffer_safe_width=buffer_safe_width):
                print("Given target not inside GALEX FoV for %s - skipped"%fullpath)
                continue
            
            inst_.set_sky(fullpath.replace("int","skybg"))
            inst_.set_target(target)
            instru.append(inst_)
            
        return instru     
    
# ====================== #
#                        #
#    SDSS                #
#                        #
# ====================== #
SDSS_DIR            = io.DATAPATH+"SDSS/"
_SDSS_METADATA_FILE = SDSS_DIR+"target_source.csv"

SDSS_BASEURL = "https://dr12.sdss.org"

def query_sdss_metadata(ra, dec):
    """ look for SDSS meta data inside MAST archive """
    import requests
    url_ = SDSS_BASEURL+'/fields/raDec?ra=%.5f&dec=%+.5f'%(ra,dec)
    r = requests.post(url_)
    r.raise_for_status() # raise a status if issue, like wrong auth

    fitsband = [l for l in r.text.splitlines() if "FITS" in l]
    if len(fitsband)==0:
        raise AttributeError("Was not able to find sdss data while searching for \n %s"%url_)

    filters  = []
    basename = []
    baseurl  = []
    for lband in fitsband:
        url = lband.split('"')[1]
        filters.append(lband.split("-band")[0][-1])
        baseurl.append("/".join([SDSS_BASEURL]+url.split("/")[1:-1]))
        basename.append(url.split("/")[-1].split(".")[0])
        
    df =  pandas.DataFrame(np.asarray([filters,baseurl,basename]).T, columns=["filters","baseurl","basename"])
    df["project"] = "dr12"
    return df
    

def _sdss_info_to_urlpath_(baseurl, basename, bands=None):
    """ band name inside the sdss url. *bands is not used.* 
    simply does: 
        "return baseurl+'/'+basename+'.fits.bz2'"
    """
    return baseurl+"/"+basename+".fits.bz2"

#
# CLASS 
#
class SDSSQuery( _Query_ ):
    """ """
    INSTRUMENT = "SDSS"



    

# ====================== #
#                        #
#  PanStarrs             #
#                        #
# ====================== #
PANSTARRS_DIR            = io.DATAPATH+"PanSTARRS/"
_PANSTARRS_METADATA_FILE = PANSTARRS_DIR+"target_source.csv"

def query_panstarrs_metadata(ra, dec, size=240, filters="grizy"):
    
    """ Query ps1filenames.py service to get a list of images
    
    Parameters
    ----------
    ra, dec: [floats] 
        position in degrees
    size: [float] 
        image size in pixels (0.25 arcsec/pixel)
    filters: [strings]
        string with filters to include
        
    Returns
    --------
    Table (a table with the results)
    """
    

    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    d_ = pandas.read_csv(url, sep=" ")
    d_["basename"] = d_.pop("shortname")
    d_["baseurl"]  = d_.pop("filename")
    d_["project"]  = "ps1"
    d_["filters"] = d_["filter"]
    return d_

def _ps_pstodata_(dataframe, directory):
    """ """
    return [directory+"%s_"%target_+"%s"%basename_
                for target_,basename_ in zip(dataframe["name"],dataframe["basename"])]
    
def _ps_pstourl_(dataframe, format="fits", size=240,  output_size=None, filters="grizy"):
    """ Get URL for images in the table
    
    Parameters
    ----------
    format: [string] 
        data format (options are "jpg", "png" or "fits")
                
    Returns
    -------
    String (a string with the URL)
    """
    if format not in ["jpg","png","fits"]:
        raise ValueError("format must be one of jpg, png, fits (%s given)"%format)
    
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}&filters={filters}").format(
               **{"ra":dataframe["ra"].values[0],"dec":dataframe["dec"].values[0],
                      "size":size,"format":format,"filters":filters})
    if output_size:
        url = url + "&output_size={}".format(output_size)
        
    # sort filters from red to blue
    urlbase = url + "&red="
    return [urlbase+filename for filename in dataframe['baseurl'].values]

def get_ps_url(ra, dec, size=240, output_size=None, filters="grizy", format="fits", color=False):
    
    """ Get URL for images in the table
    
    Parameters
    ----------
    ra, dec: [floats] 
        position in degrees

    size: [float] 
        image size in pixels (0.25 arcsec/pixel)

    filters: [strings]
        string with filters to include

    format: [string] 
        data format (options are "jpg", "png" or "fits")
        
    color: [bool]
        if True, creates a color image (only for jpg or png format).
        Default is return a list of URLs for single-filter grayscale images.
        
    Returns
    -------
    String (a string with the URL)
    """
        
    df = query_panstarrs_metadata(ra, dec, size=size, filters=filters)
    return _ps_pstourl_(df, output_size=output_size, size=size, format=format, color=color)
#
# Class
#
class PanSTARRSQuery( _Query_ ):
    """ """
    INSTRUMENT = "PanSTARRS"

    def _build_target_data_url_and_path_(self, targetname, dirout, filters=None, **kwargs):
        """ Returns the URL and download location of data """
        if not self.is_target_known(targetname):
            raise AttributeError("unknown target. Please run download_target_metadata()")

        dataframe = self.metadata[self.metadata["name"]==targetname]
        if filters is not None:
            dataframe = dataframe[dataframe["filters"].isin(filters)]
        url_ = _ps_pstourl_(dataframe)
        localpath_ = _ps_pstodata_(dataframe, dirout)
        
        return url_, localpath_
