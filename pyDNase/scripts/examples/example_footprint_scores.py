import matplotlib.pyplot as plt
import pyDNase
from pyDNase.footprinting import wellington

#Load test data
reads = pyDNase.BAMHandler(pyDNase.example_reads())
regions = pyDNase.GenomicIntervalSet(pyDNase.example_regions())

#Plot cuts data
plt.plot(reads[regions[0]]["+"],c="red")
plt.plot([-i for i in reads[regions[0]]["-"]],c="blue")

#Footprint and plot the results
footprinter = wellington(regions[0],reads)
plt.plot(footprinter.scores,c="black")

plt.show()