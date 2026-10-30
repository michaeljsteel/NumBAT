# functions for python notebooks

import time
import IPython.display as ipydisp
from pathlib import Path


# Define a function to display a single image with a caption
def img_single(fn1, cap='', height=300):
    for fn in (fn1,):
        try:
            fp = Path(fn)
        except:
            print(f"File arguments must be strings or path objects.\n Received: {fn}.")
            return
        if not fp.exists():
            print(f"Can't find image file {fn}.")
            return

    cache_buster=int(time.time() * 1000)

    html_string = f"""
        <div style="margin-right: 10px;">
            <img src="{fn1}?{cache_buster}" height="{height}px">
            <p style="text-align: left;">{cap}</p>
        </div>"""

    ipydisp.display(ipydisp.HTML(html_string))

# Define a function to display a pair of images side by side with captions
def img_pair(fn1, fn2, cap1='', cap2='', height1=300, height2=300):

    for fn in (fn1, fn2):
        try:
            fp = Path(fn)
        except:
            print(f"File arguments must be strings or path objects.\n Received: {fn}.")
            return
        if not fp.exists():
            print(f"Can't find image file {fn}.")
            return

    cache_buster=int(time.time() * 1000)
    html_string = f"""
    <div style="display: flex; justify-content: center; align-items: flex-start;">
        <div style="margin-right: 10px;">
            <img src="{fn1}?{cache_buster}" height="{height1}px">
            <p style="text-align: center;">{cap1}</p>
        </div>
        <div>
            <img src="{fn2}?{cache_buster}" height="{height2}px">
            <p style="text-align: center;">{cap2}</p>
        </div>
    </div>"""

    ipydisp.display(ipydisp.HTML(html_string))


    # Define a function to display a pair of images side by side with captions
def img_triple(fn1, fn2, fn3, cap1='', cap2='', cap3='', height1=300, height2=300, height3=300):
    for fn in (fn1, fn2, fn3):
        try:
            fp = Path(fn)
        except:
            print(f"File arguments must be strings or path objects.\n Received: {fn}.")
            return
        if not fp.exists():
            print(f"Can't find image file {fn}.")
            return

    cache_buster=int(time.time() * 1000)
    html_string = f"""
    <div style="display: flex; justify-content: center; align-items: flex-start;">
        <div style="margin-right: 10px;">
            <img src="{fn1}?{cache_buster}" height="{height1}px">
            <p style="text-align: center;">{cap1}</p>
        </div>
        <div>
            <img src="{fn2}?{cache_buster}" height="{height2}px">
            <p style="text-align: center;">{cap2}</p>
        </div>
        <div>
            <img src="{fn3}?{cache_buster}" height="{height3}px">
            <p style="text-align: center;">{cap3}</p>
        </div>
    </div>"""

    ipydisp.display(ipydisp.HTML(html_string))
