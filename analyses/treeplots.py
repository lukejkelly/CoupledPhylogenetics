"""Plotting newick trees to see which clades they support."""
# Not currently working as I cannot get Bio/biopython to install
from os.path import expanduser
from Bio import Phylo
from io import StringIO
from matplotlib import pyplot as plt

b = Phylo.read(StringIO("(1:1000,(2:500,3:500):500)"), "newick")
Phylo.draw_ascii(b)

fig, ax = plt.subplots()
Phylo.draw(b, do_show=False, axes=ax3)
# fig3.savefig(expanduser("~/Workspace/TraitLabSDLT-coupled/data/SIM_B3.pdf"))

b_x = Phylo.read(StringIO("((2:15.419,((4:4.0413,5:4.0413):3.9196,3:7.9609):7.4581):6.4127,1:21.8317)"), "newick")
b_y = Phylo.read(StringIO("((2:30.1801,3:30.1801):25.7177,((4:27.3402,5:27.3402):19.0021,1:46.3423):9.5556)"), "newick")
[Phylo.draw_ascii(lol, column_width = 20) for lol in [b_x, b_y]]
