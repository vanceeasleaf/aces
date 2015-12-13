#!/usr/bin/env python

# Copyright (C) 2011 Atsushi Togo
# All rights reserved.
#
# This file is part of phonopy.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# * Neither the name of the phonopy project nor the names of its
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import sys
import numpy as np

try:
    import yaml
except ImportError:
    print "You need to install python-yaml."
    exit(1)
    
try:
    from yaml import CLoader as Loader
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from phonopy.units import VaspToTHz



class object:
    def __init__(self, **entries):
        self.__dict__.update(entries)
 
def get_plot_data(data):
    segment_positions = []
    frequencies = []
    distances = []
    npoints = data['nqpoint'] / data['npath']

    for j, v in enumerate(data['phonon']):
        frequencies.append([f['frequency'] for f in v['band']])
        distances.append(v['distance'])
        
        if j % npoints == 0:
            segment_positions.append(v['distance'])

    return distances, frequencies, segment_positions


def plotband(factor=1.0,
                        f_max=None, 
                        f_min=None,
                        is_gnuplot=False,
                        is_points=False,
                        is_vertial_line=True,
                        output_filename='band.png',
                        labels=r'Gamma X Gamma',
                        show_legend=False,
                        title=None,filename0='band.yaml'
                        ):
    options=object(factor=factor,
                        f_max=f_max, 
                        f_min=f_min,
                        is_gnuplot=is_gnuplot,
                        is_points=is_points,
                        is_vertial_line=is_vertial_line,
                        output_filename=output_filename,
                        labels=labels,
                        show_legend=show_legend,
                        title=title
                        )
    args=[]
    if not options.is_gnuplot:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')
        if options.labels:
            
            from matplotlib import rc
            rc('text', usetex=False)

    colors = ['r-', 'g-', 'g-', 'c-', 'm-', 'y-', 'k-', 'b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--']
    if options.is_points:
        colors = [x + 'o' for x in colors]

    count = 0


    if len(args) == 0:
        filenames = [filename0]
    else:
        filenames = args

    if options.is_gnuplot:
        print "# distance  frequency (bands are separated by blank lines)"
    min1=100000
    max1=-100000
    for i, filename in enumerate(filenames):
        string = open(filename).read()
        data = yaml.load(string, Loader=Loader)
        distances, frequencies, segment_positions = get_plot_data(data)
        min1=min(min1,min(min(frequencies)))
        max1=max(max1,max(max(frequencies)))
        
        if options.is_gnuplot:
            print "# segments:",
            for v in segment_positions:
                print "%10.8f" % v,
            print "%10.8f" % distances[-1]
            
        elif options.is_vertial_line and len(filenames) == 1:
            for v in segment_positions[1:]:
                plt.axvline(x=v,linestyle='--',  linewidth=.5, color='black')

        for j, freqs in enumerate(np.array(frequencies).T):
            print freqs
            if (np.abs(freqs)<1e-0).all():continue
            if options.is_gnuplot:
                for d, f in zip(distances, freqs * options.factor):
                    print d,f
                print
            else:
                if j==0:
                    plt.plot(distances, freqs * options.factor, colors[i],
                             label=filename)
                else:
                    plt.plot(distances, freqs * options.factor, colors[i])

        if options.is_gnuplot:
            print

        
    if not options.is_gnuplot:
        plt.ylabel('Frequency(THz)')
        plt.xlabel('Wave vector')
        plt.xlim(distances[0], distances[-1])
        if not options.f_max == None:
            plt.ylim(ymax = options.f_max)
        #else: plt.ylim(ymax=max1* options.factor)
        if not options.f_min == None:
            plt.ylim(ymin = options.f_min)
        else: plt.ylim(ymin=min1* options.factor)
        plt.axhline(y=0, linestyle=':', linewidth=0.5, color='g')
        if len(filenames) == 1:
            xticks = segment_positions + [distances[-1]]
            if options.labels:
                labels = [x.replace('Gamma',r'$\Gamma$') for x in options.labels.split()]
                
                if len(labels)==len(xticks):
                    plt.xticks(xticks, labels)
                else:
                    print "Numbers of labels and band segments don't match."
                    sys.exit(1)
            else:
                plt.xticks(xticks, [''] * len(xticks))
        else:
            plt.xticks([])
        
        if not options.title == None:
            plt.title(options.title)

        if options.show_legend:
            plt.legend()
        
        if not options.output_filename == None:
            plt.rcParams['pdf.fonttype'] = 42
            plt.rcParams['font.family'] = 'serif'
            plt.savefig(options.output_filename,bbox_inches='tight',transparent=True)
        else:
            plt.show()
        plt.close()
def plotbanddos(freq,dos,factor=1.0,
                        f_max=None, 
                        f_min=None,
                        is_gnuplot=False,
                        is_points=False,
                        is_vertial_line=True,
                        output_filename='banddos.png',
                        labels=r'Gamma X Gamma',
                        show_legend=False,
                        title=None,filename0='band.yaml'):
    args=[]
    options=object(factor=factor,
                        f_max=f_max, 
                        f_min=f_min,
                        is_gnuplot=is_gnuplot,
                        is_points=is_points,
                        is_vertial_line=is_vertial_line,
                        output_filename=output_filename,
                        labels=labels,
                        show_legend=show_legend,
                        title=title
                        )
    if not options.is_gnuplot:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')
        if options.labels:
            
            from matplotlib import rc
            rc('text', usetex=False)

    colors = ['r-', 'g-', 'g-', 'c-', 'm-', 'y-', 'k-', 'b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--']
    if options.is_points:
        colors = [x + 'o' for x in colors]

    count = 0


    if len(args) == 0:
        filenames = [filename0]
    else:
        filenames = args

    if options.is_gnuplot:
        print "# distance  frequency (bands are separated by blank lines)"
    min1=100000
    max1=-100000
    plt.figure(1, (8, 6))
    plt.axes([.1, .07, .69, .85])
    for i, filename in enumerate(filenames):
        string = open(filename).read()
        data = yaml.load(string, Loader=Loader)
        distances, frequencies, segment_positions = get_plot_data(data)
        min1=min(min1,min(min(frequencies)))
        max1=max(max1,max(max(frequencies)))
        
        if options.is_gnuplot:
            print "# segments:",
            for v in segment_positions:
                print "%10.8f" % v,
            print "%10.8f" % distances[-1]
            
        elif options.is_vertial_line and len(filenames) == 1:
            for v in segment_positions[1:]:
                plt.axvline(x=v,linestyle='--',  linewidth=.5, color='black')

        for j, freqs in enumerate(np.array(frequencies).T):
            if (np.abs(freqs)<1).all():continue
            if options.is_gnuplot:
                for d, f in zip(distances, freqs * options.factor):
                    print d,f
                print
            else:
                if j==0:
                    plt.plot(distances, freqs * options.factor, colors[i],
                             label=filename)
                else:
                    plt.plot(distances, freqs * options.factor, colors[i])

        if options.is_gnuplot:
            print
    min1*=options.factor    
    if not options.is_gnuplot:
        plt.ylabel('Frequency(THz)')
        plt.xlabel('Wave vector')
        plt.xlim(distances[0], distances[-1])
        if not options.f_max == None:
            plt.ylim(ymax = options.f_max)
        #else: plt.ylim(ymax=max1* options.factor)
        if not options.f_min == None:
            plt.ylim(ymin = options.f_min)
        else: plt.ylim(ymin=min1)
        if abs(min1)>1:
            plt.axhline(y=0, linestyle=':', linewidth=0.5, color='k')
        if len(filenames) == 1:
            xticks = segment_positions + [distances[-1]]
            if options.labels:
                labels = [x.replace('Gamma',r'$\Gamma$') for x in options.labels.split()]
                
                if len(labels)==len(xticks):
                    plt.xticks(xticks, labels)
                else:
                    print "Numbers of labels and band segments don't match."
                    sys.exit(1)
            else:
                plt.xticks(xticks, [''] * len(xticks))
        else:
            plt.xticks([])
        
        if not options.title == None:
            plt.title(options.title)

        if options.show_legend:
            plt.legend()
        #plt.grid('on')
        
        aa=plt.ylim()
        # Plot the band structure and DOS
        
        plt.axes([.8, .07, .2, .85])
        #plt.fill_between(dos, freq,  color='lightgrey', edgecolor='k', lw=1)
        plt.plot(dos, freq,'k',linewidth=1)
        if abs(min1)>1:
            plt.axhline(y=0, linestyle=':', linewidth=0.5, color='k')
        plt.ylim(aa)
        #plt.ylim(ymin=min1* options.factor)
        plt.xticks([], [])
        plt.yticks([], [])
        plt.xlabel("DOS")
        if not options.output_filename == None:
            plt.rcParams['pdf.fonttype'] = 42
            plt.rcParams['font.family'] = 'serif'
            plt.savefig(options.output_filename,bbox_inches='tight',transparent=True)
        else:
            plt.show()
        plt.close()

def band3d(imgname='band3d.png',filename='life.txt'):
    from aces.graph import surf,scatter3d
    import pandas as pd
    import numpy as np
    data=pd.read_csv(filename,sep='\t')
    scatter3d(data['kx'],data['ky'],data['w0'],filename=imgname)
    """
    surf(np.array([data['kx']]).reshape(-1,8),
        np.array([data['ky']]).reshape(-1,8),
        np.array([data['freq']]).reshape(-1,8),filename=imgname)
"""

