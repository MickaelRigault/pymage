#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from astrobject.baseobject import AstroTarget

from . import query
GALEXQuery = query.GALEXQuery()
SDSSQuery  = query.SDSSQuery()

class PhotoTarget( AstroTarget ):
    """ """
    PROPERTIES = ["instruments"]

    # -------- #
    #  LOADER  #
    # -------- #
    def load_instrument(self, which=["galex", "sdss"]):
        """ """
        which = np.atleast_1d(which)
        if "galex" in which:
            self._load_galex_()
        if "sdss" in which:
            self._load_sdss_()
        
    def _load_galex_(self):
        """ """
        if "galex" not in self.instruments:
            self.instruments["galex"] = {"nuv":[], "fuv":[]}
            
        for inst_ in GALEXQuery.get_target_instruments(self.name):
            inst_.set_target(self, test_inclusion=False)
            self.instruments["galex"][inst_.bandname].append(inst_)

    def _load_sdss_(self):
        """ """
        if "sdss" not in self.instruments:
            self.instruments["sdss"] = {"sdssu":[], "sdssg":[], "sdssr":[], "sdssi":[], "sdssz":[]}
            
        for inst_ in SDSSQuery.get_target_instruments(self.name):
            inst_.set_target(self, test_inclusion=False)
            self.instruments["sdss"][inst_.bandname].append(inst_)
        
    # -------- #
    #  GETTER  #
    # -------- #
    def get_photopoints(self, instruband, radius, runits="arcsec", **kwargs):
        """ instruband = instrument.band i.e. galex.nuv or sdss.sdssu """
        instru, band = instruband.split(".")
        return [inst_.get_target_photopoint(radius, runits=runits, **kwargs) 
                for inst_ in self.instruments[instru][band]]
    
    # ================ #
    #   Properties     #
    # ================ #
    @property
    def instruments(self):
        """ """
        if self._properties["instruments"] is None:
           self._properties["instruments"] = {}
        return self._properties["instruments"]
    
