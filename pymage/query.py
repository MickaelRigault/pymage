#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pandas
import numpy as np

from astroquery.mast import Observations

from . import io

GALEX_DIR            = io.DATAPATH+"GALEX/"
_GALEX_METADATA_FILE = GALEX_DIR+"target_source.csv"
GALEX_METADATA       = pandas.read_csv(_GALEX_METADATA_FILE) if os.path.exists(_GALEX_METADATA_FILE) else \
  pandas.DataFrame(columns=["name","ra","dec","filters", "project","dataURL", "basename", "baseurl"])

    
##########################
#                        #
#    GENERAL TOOLS       #
#                        #
##########################
def query_mast(ra, dec, instrument=None, radius="10 arcsec"):
    """ returns a table containing the existing data information """
    from astroquery.mast import Observations
    t = Observations.query_region("%.5f %+.5f"%(ra,dec), radius="10 arcsec")
    if instrument is not None:
        return t[t["obs_collection"]==instrument]
    
    return t


# ====================== #
#                        #
#    GENERAL TOOLS       #
#                        #
# ====================== #
def _galex_info_to_urlpath_(baseurl, basename, todl=["int", "skybg"], bands=["NUV","FUV"]):
    """ Build URL or fullpath for the given data """
    
    return [[baseurl+"/"+basename+"-%sd-%s.fits.gz"%(band[0].lower(),todl_) for band in np.atleast_1d(bands)]
           for todl_ in np.atleast_1d(todl)]
    
def query_galex_metadata(ra, dec, dirout="default"):
    """ look for GALEX meta data inside MAST archive """
    if dirout == "default":
        dirout = GALEX_DIR
        
    t = query_mast(ra, dec, instrument="GALEX")
    df = pandas.DataFrame(dict(t[["filters", "project","target_name"]]))
    df["basename"] = df.pop("target_name")
    df = df.assign(baseurl= ["/".join(url_.split("/")[:-1]) for url_ in t["dataURL"].data])
    return df

                        
class GALEXQuery( object ):
    """ Simply Class to manage the GALEX data IO """
    def __init__(self):
        """ initialize Galex query. It loads the known galex meta data. """
        self.metadata = GALEX_METADATA
    
    def download_target_metadata(self, targetname, ra, dec, 
                                 update=True, store=True, 
                                 overwrite=False, dl=True):
        """ Look for GALEX metadata inside MAST archive, save them and download the corresponding galex files. 
        (options enables to store or not and download or not)

        Parameters:
        -----------
        targetname: [str]
            galex metadata will be associated to this name.

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
        
        df_ = query_galex_metadata(ra,dec)
        df_["name"] = targetname
        df_["ra"] = ra
        df_["dec"] = dec
        if not update:
            return df_
        
        # Downloading
        
        # merge with current
        self.metadata.drop(index=self.metadata.index[self.metadata["name"]==targetname], inplace=True)
        self.metadata = pandas.concat([self.metadata, df_], sort=False)
        if store:
            self.metadata.to_csv(_GALEX_METADATA_FILE, index=False)
        
        if dl:
            self.download_target_data(targetname)
            
    def get_target_data(self, targetname, must_exists=True, which=["int","skygb"]):
        """ returns the full path of galex data on your computer. """
        urls, localpaths = self._build_target_data_url_and_path_(targetname, GALEX_DIR, which)
        return [l_ for l_ in localpaths if os.path.exists(l_) or not must_exists]
             
    def download_target_data(self, targetname, todl=["int","skygb"], dirout="default", overwrite=False, 
                            load_metadata=True):
        """ Download the target galex data. """            
        if dirout is None or dirout in ["default"]:
            dirout = GALEX_DIR
            
        urls, localpaths = self._build_target_data_url_and_path_(targetname, dirout, todl)
        for url_, localpath_ in zip(urls, localpaths):
            io.download_single_url(url_, localpath_, overwrite=overwrite)
        
    def _build_target_data_url_and_path_(self, targetname, dirout, todl=["int","skygb"]):
        """ Returns the URL and download location of GALEX data """
        if not self.is_target_known(targetname):
            raise AttributeError("unknown target. Please run download_target_metadata()")

        
        url_ = np.asarray([_galex_info_to_urlpath_(row["baseurl"], row["basename"], bands=row["filters"])
               for index_, row in self.metadata[self.metadata["name"]==targetname].iterrows()]).flatten()
        localpath_ = np.asarray([_galex_info_to_urlpath_(dirout, row["basename"], bands=row["filters"])
                for index_, row in self.metadata[self.metadata["name"]==targetname].iterrows()]).flatten()
        
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
