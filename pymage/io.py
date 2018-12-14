#! /usr/bin/env python
# -*- coding: utf-8 -*-


import os
DATAPATH = os.getenv("DATAPATH")


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
            print("downloading %s"%fileout)
    
    request_fnc = "get" if not "data" in kwargs else "post"
    cout = getattr(requests,request_fnc)(url, **kwargs)
    with open(fileout, 'wb') as f:
        for data in cout.iter_content(chunk):
            f.write(data)
