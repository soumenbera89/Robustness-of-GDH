# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 13:49:08 2016

@author: SOUMEN
"""

import numpy as np
import matplotlib.pyplot as plt
NH4, KG, Gibbs = np.loadtxt("gibbs.csv", skiprows=1, unpack = True, delimiter= ',')
NH41, KG1, Gibbs1 = np.loadtxt("gibbs_normal.csv", skiprows=1, unpack = True, delimiter= ',')
x = np.linspace(0,100,200)
y = -0.0000001*x**3+0.00005*x**2-0.0071*x+0.0995
fig = plt.figure(figsize=(12,6))
ax1 = fig.add_subplot(121)
#plt.figure(1)
ax1.plot(x, y, 'b-', label='theory')
ax1.plot(Gibbs,'bs')
ax1.axhline(linewidth=1, linestyle='--',color='k')
#plt.axvspan(0.0, 19.0,facecolor='0.6', alpha=0.1)
#plt.axvspan(19.0, 100.0,facecolor='0.6', alpha=0.2)
ax1.set_xlabel("Extent of Reaction")
ax1.set_ylabel("Gibbs function")
#ax.annotate('Backward', xy=(10, 0.15), xytext=(20, 0.05),
#           arrowprops=dict(arrowstyle="->"),
#            )
ax1.annotate('Amination Reaction', xy=(80, -0.35), xytext=(20, -0.35),
            arrowprops=dict(arrowstyle="->"),
            )
# Hide the right and top spines            
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
ax1.text(20, 0.11, r'$\mathit{y=-1E^{-07}x^{3}+5E^{-05}x^{2}-0.0071x+0.0995}$')
ax1.text(40, 0.075, r'$\mathit{R^{2}=0.7453}$')
ax1.grid(which='major', axis='both', linestyle='-', alpha=0.075)   
ax2 = fig.add_subplot(122)
ax2.plot(Gibbs1,'b-')
ax2.set_xlabel("Extent of Reaction")
ax2.set_ylabel("Gibbs function")
# Hide the right and top spines            
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')
ax2.grid(which='major', axis='both', linestyle='-', alpha=0.075) 
ax2.annotate('Amination Reaction', xy=(80, -0.35), xytext=(20, -0.35),
            arrowprops=dict(arrowstyle="->"),
            )
ax2.axhline(linewidth=1, linestyle='--',color='k')
plt.show()
fig.savefig('gibbs1.png', format='png', dpi=300)
plt.close(fig)
#ax.view_init(elev=elevation_angle, azim=azimuthal_angle)
