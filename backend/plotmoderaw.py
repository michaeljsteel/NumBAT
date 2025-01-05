
# Copyright (C) 2017-2025  Michael Steel, Bjorn Sturmberg, Kokou Dossou.

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


import numpy as np
import numpy.random
import matplotlib.pyplot as plt

from nbtypes import SI_um
import plottools

def do_raw_fem_mode_plot(l_comps, mode_helper, fem_mesh, fem_evecs, mode_num):

    prefix = 'ttraw'

    fig, ax = plt.subplots(dpi=300)

    v_F = None
    v_x, v_y = fem_mesh.get_fullmesh_nodes_xy()
    v_x *= 1/SI_um
    v_y *= 1/SI_um

    sigx = (np.max(v_x)-np.min(v_x))/400
    sigy = (np.max(v_y)-np.min(v_y))/400

    rng = numpy.random.default_rng()
    v_x += rng.uniform(-0.5,0.5,len(v_x))*sigx
    v_y += rng.uniform(-0.5,0.5,len(v_y))*sigy

    for comp in l_comps:
        icomp = {'x':0, 'y':1, 'z': 2}[comp]

        # now an array [6,n_el], field at 6 nodes for each elemet
        v_F = fem_evecs[icomp, :6, mode_num, :]
        if comp == 'z':
            v_Fr = np.imag(v_F)
        else:
            v_Fr = np.real(v_F)


        # flatten to an array of 6 x n_msh_el
        v_Fr_flat= v_Fr.flatten('F')

        #sizes =np.abs(v_Fr_flat)*20
        sizes=.2
        colors =v_Fr_flat
        vm = np.max(np.abs(v_Fr_flat))

        ax.scatter(v_x, v_y, s=sizes, c=colors, marker='o', vmin=-vm/5, vmax=vm/5, cmap='seismic')
        ax.set_xlabel(r'$x$ [μm]')
        ax.set_ylabel(r'$y$ [μm]')


        plottools.save_and_close_figure(fig, prefix+f'_{comp}')
