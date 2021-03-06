#! /usr/bin/env python
# -*- coding: utf-8 -*-


import os
import io
import warnings
DATAPATH = os.getenv("DATAPATH","_notdefined_")
if DATAPATH == "_notdefined_":
    warnings.warn("You don't have a global variable named 'DATAPATH'. You need one to be able to download data. They will be stored in $DATAPATH/{INSTRUNAME}/bla")

def download_single_url(url, fileout=None, mkdir=True,
                        overwrite=False, verbose=True, chunk=1024, **kwargs):
    """ Download the url target using requests.get.
    the data is returned (if fileout is None) or stored in `fileout`
    Pa
    """
    import requests
    
    if fileout is not None and not overwrite and os.path.isfile( fileout ):
        if verbose:
            print("%s already exists: skipped"%fileout)
        return
    else:
        if verbose and fileout:
            print("downloading %s -> %s"%(url,fileout))
    
    request_fnc = "get" if not "data" in kwargs else "post"
    response = getattr(requests,request_fnc)(url, **kwargs)
    if response.status_code == 200:
        if fileout in ["BytesIO", "StringIO"]:
            return getattr(io, fileout)(response.content)
        
        with open(fileout, 'wb') as f:
            for data in response.iter_content(chunk):
                f.write(data)
    else:
        print("Issue downloading")
        print("response.status_code: ", response.status_code)
