from os.path import expanduser
from Bio import Phylo
from io import StringIO
from matplotlib import pyplot as plt


b3 = Phylo.read(StringIO("(1:1000,(2:500,3:500):500)"), "newick")
b4 = Phylo.read(StringIO("((1:487.593,2:487.593):512.407,(3:393.0037,4:393.0037):606.9963)"), "newick")
b5 = Phylo.read(StringIO("((((4:88.3817,5:88.3817):34.4521,3:122.8338):503.1649,2:625.9987):374.0013,1:1000)"), "newick")
b10 = Phylo.read(StringIO("((3:422.9429,((6:215.9649,7:215.9649):126.9294,(8:66.4352,(9:55.0157,10:55.0157):111.4194):176.4591):80.0486):577.0571,(1:465.4652,((4:382.6481,5:382.6481):224.7189,2:607.367):358.0982):34.5348)"), "newick")

Phylo.draw_ascii(b3)
Phylo.draw_ascii(b4)
Phylo.draw_ascii(b5)
Phylo.draw_ascii(b10)

fig3, ax3 = plt.subplots()
Phylo.draw(b3, do_show=False, axes=ax3)
fig3.savefig(expanduser("~/Workspace/TraitLabSDLT-coupled/data/SIM_B3.pdf"))

fig4, ax4 = plt.subplots()
Phylo.draw(b4, do_show=False, axes=ax4)
fig4.savefig(expanduser("~/Workspace/TraitLabSDLT-coupled/data/SIM_B4.pdf"))

fig5, ax5 = plt.subplots()
Phylo.draw(b5, do_show=False, axes=ax5)
fig5.savefig(expanduser("~/Workspace/TraitLabSDLT-coupled/data/SIM_B5.pdf"))

fig10, ax10 = plt.subplots()
Phylo.draw(b10, do_show=False, axes=ax10)
fig10.savefig(expanduser("~/Workspace/TraitLabSDLT-coupled/data/SIM_B10.pdf"))


b5x = Phylo.read(StringIO("((2:15.419,((4:4.0413,5:4.0413):3.9196,3:7.9609):7.4581):6.4127,1:21.8317)"), "newick")
b5y = Phylo.read(StringIO("((2:30.1801,3:30.1801):25.7177,((4:27.3402,5:27.3402):19.0021,1:46.3423):9.5556)"), "newick")
[Phylo.draw_ascii(lol, column_width = 20) for lol in [b5x, b5y]]
