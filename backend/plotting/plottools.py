
"""
plottools.py
-------------
Utility functions for figure saving, image manipulation, and progress display in NumBAT plotting routines.

Includes:
- RGB color calculation for polarization vectors
- Figure saving with project-specific preferences
- Image cropping, joining, and border trimming
- Terminal progress bar for loops

Copyright (C) 2017-2025  Michael Steel.
"""

# NumBAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


from PIL import Image, ImageOps
from pathlib import Path

import math
import matplotlib.pyplot as plt
import numpy as np

def get_rgb_for_poln(px, py, pz):
    """
    Calculate normalized RGB values for a polarization vector.

    Parameters
    ----------
    px, py, pz : float or complex
        Components of the polarization vector.

    Returns
    -------
    rgb : ndarray
        Normalized RGB values (0-1) as a numpy array.
    """
    vp = np.array([px, py, pz])
    pmod = math.sqrt(px*px.conjugate() + py*py.conjugate() + pz.pz.conjugate())
    # RGB values are in range 0-1
    rgb = np.abs(vp) / pmod
    return rgb


def save_and_close_figure(fig, fig_fname):
    """
    Save a matplotlib figure to file with project-specific DPI and close it.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    fig_fname : str
        Output filename. PNG files are saved with default DPI, others with tight bounding box.
    """
    import numbat  # here to avoid circular import

    plot_prefs = numbat.NumBATPlotPrefs()

    if plot_prefs is None:  # If NumBATApp() has not been called this will be None
        dpi = 300
    else:
        dpi = plot_prefs.plot_output_resolution_dpi

    if fig_fname[-3:-1] == 'png':
        fig.savefig(fig_fname, dpi=dpi)
    else:
        fig.savefig(fig_fname, bbox_inches='tight', dpi=dpi)

    plt.close(fig)


def fig_trim_border(im, clip=None, trimwhite=False, border=20):
    """
    Crop and trim borders from a PIL image, optionally removing white space.

    Parameters
    ----------
    im : PIL.Image
        The image to process.
    clip : tuple or None
        (left, top, right, bottom) pixels to crop from each side.
    trimwhite : bool
        If True, trim white border from the image.
    border : int
        Extra border to add after trimming.

    Returns
    -------
    im : PIL.Image
        The processed image.
    """
    if clip:
        sz = im.size
        area = (clip[0], clip[1], sz[0]-clip[2], sz[1]-clip[3])
        im = im.crop(area)
    if trimwhite:
        if im.mode != 'RGB':
            im = im.convert('RGB')
        imrev = ImageOps.invert(im)
        bbox = np.array(imrev.getbbox())
        if border:
            bbox2 = bbox + np.array([-border, -border, border, border])
        im = im.crop(bbox2)
    return im


def join_figs(l_fns, fnout, clip=None, trimwhite=False, border=20, delete_inputs=False):
    """
    Join two images horizontally, matching heights, and save the result.

    Parameters
    ----------
    l_fns : list of str
        List of input image filenames (should be two).
    fnout : str
        Output filename.
    clip : tuple or None
        (left, top, right, bottom) pixels to crop from each side.
    trimwhite : bool
        If True, trim white border from the image.
    border : int
        Extra border to add after trimming.
    delete_inputs : bool
        If True, delete the input files after joining.
    """
    images = []
    for fn in l_fns:
        im = Image.open(fn)
        images.append(im)

    # Match heights of first two images
    nsz0 = [images[1].size[0], int(images[0].size[1] * images[1].size[0] / images[0].size[0])]
    images[0] = images[0].resize(nsz0)

    widths, heights = zip(*(i.size for i in images))
    total_width = sum(widths)
    max_height = max(heights)

    new_im = Image.new('RGB', (total_width, max_height))

    x_offset = 0
    yoffs = [0, 0]
    if images[0].size[1] > images[1].size[1]:
        yoffs = [0, int((images[0].size[1] - images[1].size[1]) / 2)]
    else:
        yoffs = [int((images[1].size[1] - images[0].size[1]) / 2), 0]

    for i, im in enumerate(images):
        new_im.paste(im, (x_offset, yoffs[i]))
        x_offset += im.size[0]

    new_im.save(fnout)

    if delete_inputs:
        for fn in l_fns:
            Path(fn).unlink()

#fill = '█',   # TODO: this messes with pdflatex in docs. Fix

def progressBar(sequence, prefix='', suffix='', decimals=1, length=100, fill='x', printEnd="\r"):
    """
    Display a terminal progress bar for an iterable sequence.

    Parameters
    ----------
    sequence : iterable
        The sequence to iterate over.
    prefix : str
        Prefix string for the progress bar.
    suffix : str
        Suffix string for the progress bar.
    decimals : int
        Number of decimals in percent complete.
    length : int
        Character length of the bar.
    fill : str
        Bar fill character.
    printEnd : str
        End character (e.g., "\r", "\r\n").

    Yields
    ------
    item : object
        Items from the input sequence, with progress bar display.
    """
    total = len(sequence)
    def printProgressBar(iteration):
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=printEnd)

    if len(sequence) > 1:
        printProgressBar(0)  # Initial Call
        for i, item in enumerate(sequence):   # Update Progress Bar
            yield item
            printProgressBar(i + 1)
    else:
        yield sequence[0]
    print()  # Print New Line on Complete

