# Copyright (C) 2013 Jason Piper - j.piper@warwick.ac.uk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import matplotlib.pyplot as plt
import pyDNase
from pyDNase.footprinting import wellington

#Load test data
reads = pyDNase.BAMHandler("example.bam")
regions = pyDNase.GenomicIntervalSet("example.bed")

#Plot cuts data
plt.plot(reads[regions[0]]["+"],c="red")
plt.plot(-reads[regions[0]]["-"],c="blue")

#Footprint and plot the results
footprinter = wellington(regions[0],reads)
plt.plot(footprinter.scores,c="black")

plt.show()