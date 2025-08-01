# functions for python notebooks

import time
import IPython.display as ipydisp


# Define a function to display a single image with a caption
def img_single(fn, cap='', height=300):
    cache_buster=int(time.time() * 1000)

    html_string = f"""
        <div style="margin-right: 10px;">
            <img src="{fn}?{cache_buster}" height="{height}px">
            <p style="text-align: left;">{cap}</p>
        </div>"""

    ipydisp.display(ipydisp.HTML(html_string))

# Define a function to display a pair of images side by side with captions
def img_pair(fn1, fn2, cap1='', cap2='', height1=300, height2=300):
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