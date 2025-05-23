# Copyright (C) 2017-2025  Michael Steel.

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


import matplotlib.pyplot as plt
from PIL import Image

def save_and_close_figure(fig, fig_fname):

    if fig_fname[-3:-1] == 'png':
        fig.savefig(fig_fname)
    else:
        fig.savefig(fig_fname, bbox_inches='tight')

    plt.close(fig)


def join_figs(l_fns, fnout, clip=None):

    images = []
    for fn in l_fns:
        im = Image.open(fn)

        if clip is not None:
            sz = im.size
            area = (clip[0], clip[1], sz[0]-clip[2], sz[1]-clip[3])
            im = im.crop(area)
        images.append(im)

    #matches heights of first two images
    nsz0 = [images[1].size[0],
            int(images[0].size[1]*
            images[1].size[0]/
            images[0].size[0])]
    images[0]=images[0].resize(nsz0)

    widths, heights = zip(*(i.size for i in images))

    total_width = sum(widths)
    max_height = max(heights)

    new_im = Image.new('RGB', (total_width, max_height))

    x_offset = 0
    yoffs  = [0,0]
    if images[0].size[1]>images[1].size[1]:
        yoffs=[0, int((images[0].size[1]-images[1].size[1])/2)]
    else:
        yoffs=[int((images[1].size[1]-images[0].size[1])/2), 0]

    for i, im in enumerate(images):
        new_im.paste(im, (x_offset,yoffs[i]))
        x_offset += im.size[0]

    new_im.save(fnout)

#fill = '█',   # TODO: this messes with pdflatex in docs. Fix
def progressBar(sequence, prefix = '', suffix = '',
                decimals = 1, length = 100,
                fill = 'x',
                printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
        @params:
            iterable    - Required  : iterable object (Iterable)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
            printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)

    """
    total = len(sequence)
    # Progress Bar Printing Function
    def printProgressBar (iteration):
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)

    if len(sequence) > 1:
        printProgressBar(0)  # Initial Call
        for i, item in enumerate(sequence):   # Update Progress Bar
            yield item
            printProgressBar(i + 1)
    else:
        yield sequence[0]

    print() # Print New Line on Complete