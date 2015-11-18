__author__ = "Jason Piper"

"""
Dear god this is the ugliest piece of shit source I've ever written.
Here's a cat to cheer you up until I fix it.

                \`*-.
                 )  _`-.
                .  : `. .
                : _   '  \
                ; *` _.   `*-._
                `-.-'          `-.
                  ;       `       `.
                  :.       .        \
                  . \  .   :   .-'   .
                  '  `+.;  ;  '      :
                  :  '  |    ;       ;-.
                  ; '   : :`-:     _.`* ;
        [meow] .*' /  .*' ; .*`- +'  `*'
               `*-*   `*-*  `*-*'

"""

from matplotlib import rcParams
import pyDNase
import pylab as plt
import numpy as np
from clint.textui import progress
import argparse

parser = argparse.ArgumentParser(description='Plots a report of cuts, footprint scores, and conservation scores')
parser.add_argument("-w", "--window_size", help="Size of flanking area around centre of the regions to plot (default: 50)",default=50,type=int)
parser.add_argument("-y", help="ymax (default: auto)",default=0,type=int)
parser.add_argument("regions", help="BED file of the regions you want to generate the average profile for")
parser.add_argument("reads", help="The BAM file containing the DNase-seq data")
parser.add_argument("output", help="filename to write the output to (use .pdf or .png)")
args  = parser.parse_args()

xsize   = args.window_size
reads   = pyDNase.BAMHandler(args.reads)
regions = pyDNase.GenomicIntervalSet(args.regions)

#Set all strands to positive
#for each in regions:
#    each.strand = "+"

regions.resizeRegions(xsize)

fw = []
rv = []
#TODO: Make this memory efficient - we don't need to store all the fw and rvs
plt.figure(num=None,figsize=(4,12))
#plt.subplot(211)
plt.subplot2grid((4,1),(0, 0))
print("Plotting cut counts...")
for each in progress.bar(regions):
    if sum(reads[each]["+"]) and sum(reads[each]["-"]):
        fw.append(reads[each]["+"])
        rv.append(reads[each]["-"])


fw = [map(float,i) for i in fw]
rv = [map(float,i) for i in rv]
fw = [np.divide(np.subtract(i, min(i)), np.subtract(max(i) , min(i))) for i in fw]
rv = [np.divide(np.subtract(i, min(i)), np.subtract(max(i) , min(i))) for i in rv]


plt.plot(np.mean(fw,axis=0),c="red")
plt.plot(np.mean(rv,axis=0),c="blue")

del(fw)
del(rv)

#Pad the axis out reads bit
rcParams['xtick.major.pad'] = 20
rcParams['ytick.major.pad'] = 20

#Sort out the X axis
ticks = [0,xsize,xsize*2]
labels = [-xsize,0,xsize]
plt.xticks(ticks, labels)

#Makes ticks only appear on the left hand side
plt.gca().yaxis.set_ticks_position('left')

#Remove top, right, and bottom borders
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
#plt.gca().spines['bottom'].set_visible(False)

#plt.gca().tick_params(axis='both', which='major', labelsize=20)
if args.y:
    plt.gca().set_ylim(0,args.y)


plt.gca().set_ylabel('Average DNase\n Activity', multialignment='center')

plt.subplot2grid((4,1),(1, 0), rowspan=3)

matrix = []
for i in progress.bar(sorted(regions,key=lambda x: x.importorder)):
    cuts = reads[i]
    cuts["+"] = np.array(cuts["+"],dtype="float")
    cuts["-"] = np.array(cuts["-"],dtype="float")

    pcuts = np.divide(np.subtract(cuts["+"], min(cuts["+"])), np.subtract(max(cuts["+"]) , min(cuts["+"])))
    ncuts = np.divide(np.subtract(cuts["-"], min(cuts["-"])), np.subtract(max(cuts["-"]) , min(cuts["-"])))

   # pcuts = (cuts["+"] - cuts["+"].min()) / (cuts["+"].max() - cuts["+"].min())
    #ncuts = (cuts["-"] - cuts["-"].min()) / (cuts["-"].max() - cuts["-"].min())
    matrix.append((pcuts - ncuts).tolist())
    #matrix.append((cuts["+"] - cuts["-"]).tolist())


cdict1 = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.1),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 1.0),
                   (0.5, 0.1, 0.0),
                   (1.0, 0.0, 0.0))
        }
import matplotlib.colors
from matplotlib.colors import LinearSegmentedColormap
blue_red1 = LinearSegmentedColormap('BlueRed1', cdict1)
plt.pcolormesh(np.array(matrix),cmap = blue_red1,norm = matplotlib.colors.Normalize(-0.5,0.5))

ticks = [0,xsize,xsize*2]
labels = [-xsize,0,xsize]
plt.xticks(ticks, labels)

plt.tight_layout()
plt.savefig(args.output)
#plt.show()