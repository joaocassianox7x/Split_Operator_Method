import numpy as np
import seaborn as sea
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import os
import glob


os.system('rm figs/*.png')
z0=np.loadtxt("out/data0.txt")
norm = np.linalg.norm(z0)
pot = np.loadtxt("potencial.dat")

txt_files = glob.glob("out/*.txt")

def plotar(i):
    z = np.loadtxt("out/data"+str(i)+".txt")
    sea.heatmap(np.amax(z)*(pot>0)+z,cbar=False)
    plt.savefig("figs/figura%05d.png"%(i))
    plt.clf()
    plt.close()
    

Parallel(n_jobs=2)(delayed(plotar)(x) for x in range(len(txt_files)))


os.system('ffmpeg -r:v 24 -i "figs/figura%05d.png" -codec:v libx264 -preset veryslow -codec:v libx264 -preset veryslow -pix_fmt yuv420p -crf 28 -an "teste.mp4"')