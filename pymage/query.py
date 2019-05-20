#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pandas
import numpy as np

from astroquery.mast import Observations

from . import io



KNOWN_INSTRUMENTS  = ["galex", "sdss"]

def query_metadata(instrument, ra, dec):
    """ query metadata for the given intrument and coordinates
    the metadata are information about filename and url
    """
    return eval("query_%s_metadata(ra, dec)"%_test_instrument_(instrument))

def load_metadata(instrument, load_empty=False):
    """ Load existing metadata for the given instrument 
    Returns
    -------
    Pandas.DataFrame #  columns are (at least) "name","ra","dec","filters", "project", "basename", "baseurl"
    """
    empty = pandas.DataFrame(columns=["name","ra","dec","filters", "project", "basename", "baseurl"])
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

    
    def get_target_instruments(self, targetname, contains=None):
        """ Return a list of Astrobject's Instrument for each entry coresponding to the given target """
        if not self.is_target_known(targetname):
            raise AttributeError("unknown target. Please run download_target_metadata(), and download the associated files")

        from astrobject import instruments
        target_data = self.get_target(targetname)
        return [instruments.get_instrument(f_, astrotarget=target_data)
                for f_ in self.get_target_data(targetname) if contains is None or contains in f_]
    
    
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
            fileout = _get_metadata_file_(self.INSTRUMENT)
            if not os.path.dirname(fileout):
                if io.DATAPATH == "_notdefined_":
                    raise AttributeError("You must define the global variable DATAPATH to bve able to download/store data")
                
                os.mkdir(os.path.dirname(fileout))
                
            self.metadata.to_csv(fileout, index=False)
            
        # Downloading
        if dl:
            self.download_target_data(targetname)

    def download_target_data(self, targetname, dirout="default",
                                 overwrite=False, load_metadata=True, **kwargs):
        """ Download the target photometric data. """
        if dirout is None or dirout in ["default"]:
            dirout = self._default_dldir
        if io.DATAPATH == "_notdefined_":
            raise AttributeError("You must define the global variable DATAPATH to bve able to download/store data")
        
        urls, localpaths = self._build_target_data_url_and_path_(targetname, dirout, **kwargs)
        for url_, localpath_ in zip(urls, localpaths):
            io.download_single_url(url_, localpath_, overwrite=overwrite)
        
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
GALEX_DIR            = io.DATAPATH+"GALEX/"
_GALEX_METADATA_FILE = GALEX_DIR+"target_source.csv"

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
    df = df.assign(baseurl= ["/".join(url_.split("/")[:-1]) for url_ in t["dataURL"].data])
    return df

#
# CLASS 
#
class GALEXQuery( _Query_ ):
    """ Simply Class to manage the GALEX data IO """
    INSTRUMENT = "GALEX"
    def get_target_instruments(self, targetname, contains=None):
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
            if not inst_.is_target_in(target):
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
