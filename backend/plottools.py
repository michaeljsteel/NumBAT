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
