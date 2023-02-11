import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors
import matplotlib.cbook as cbook
from matplotlib.ticker import LogLocator, AutoLocator
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
fig=plt.figure(figsize=(7,8))
plt.subplots_adjust(wspace=0.01,hspace=0.01)#bottom=0.01,left=0.01,top=0.99,right=0.99,

ncolor = 9
cmap = plt.get_cmap("viridis", ncolor)
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
ax11=plt.subplot2grid((3,9), (0,0), colspan=8, rowspan=1)
ax11.tick_params(axis="x", direction="in", which='both', labelsize=13)
ax11.tick_params(axis="y", direction="in",               labelsize=13)

plt.xscale('log')

C2C04_rotacf_300 = np.loadtxt('300/rotacf_P2_ilht.xvg')
C2C04_rotacf_325 = np.loadtxt('325/rotacf_P2_ilht.xvg')
C2C04_rotacf_350 = np.loadtxt('350/rotacf_P2_ilht.xvg')
C2C04_rotacf_375 = np.loadtxt('375/rotacf_P2_ilht.xvg')
C2C04_rotacf_400 = np.loadtxt('400/rotacf_P2_ilht.xvg')
C2C04_rotacf_425 = np.loadtxt('425/rotacf_P2_ilht.xvg')
C2C04_rotacf_450 = np.loadtxt('450/rotacf_P2_ilht.xvg')
C2C04_rotacf_475 = np.loadtxt('475/rotacf_P2_ilht.xvg')
C2C04_rotacf_500 = np.loadtxt('500/rotacf_P2_ilht.xvg')

ax11.plot(C2C04_rotacf_300[:,0], C2C04_rotacf_300[:,1], '-', label="300", color=cmap(0))
ax11.plot(C2C04_rotacf_325[:,0], C2C04_rotacf_325[:,1], '-', label="325", color=cmap(1))
ax11.plot(C2C04_rotacf_350[:,0], C2C04_rotacf_350[:,1], '-', label="350", color=cmap(2))
ax11.plot(C2C04_rotacf_375[:,0], C2C04_rotacf_375[:,1], '-', label="375", color=cmap(3))
ax11.plot(C2C04_rotacf_400[:,0], C2C04_rotacf_400[:,1], '-', label="400", color=cmap(4))
ax11.plot(C2C04_rotacf_425[:,0], C2C04_rotacf_425[:,1], '-', label="425", color=cmap(5))
ax11.plot(C2C04_rotacf_450[:,0], C2C04_rotacf_450[:,1], '-', label="450", color=cmap(6))
ax11.plot(C2C04_rotacf_475[:,0], C2C04_rotacf_475[:,1], '-', label="475", color=cmap(7))
ax11.plot(C2C04_rotacf_500[:,0], C2C04_rotacf_500[:,1], '-', label="500", color=cmap(8))

ax11.set_xlim(0.1, 10000)
ax11.set_ylim(0, 1)
#ax11.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
#ax11.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax11.text(0.05, 0.11,'(a)',ha='center',va='center',size=16,color='black',transform=ax11.transAxes)
#ax11.leg=plt.legend(loc='lower right', bbox_to_anchor=[1.04, 0.18], fontsize=12, labelspacing=0.3, shadow=False, fancybox=False, frameon=False)
#ax11.leg.get_frame().set_alpha(0.2)
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
ax21 = plt.subplot2grid((3,9), (1,0), colspan=8, rowspan=1)
ax21.set_ylabel("Rotational correlation functions",      fontsize=13)
ax21.tick_params(axis="x", direction="in", which='both', labelsize=13)
ax21.tick_params(axis="y", direction="in",               labelsize=13)

plt.xscale('log')

C2C12_rotacf_300 = np.loadtxt('300/rotacf_P2_ilht.xvg')
C2C12_rotacf_325 = np.loadtxt('325/rotacf_P2_ilht.xvg')
C2C12_rotacf_350 = np.loadtxt('350/rotacf_P2_ilht.xvg')
C2C12_rotacf_375 = np.loadtxt('375/rotacf_P2_ilht.xvg')
C2C12_rotacf_400 = np.loadtxt('400/rotacf_P2_ilht.xvg')
C2C12_rotacf_425 = np.loadtxt('425/rotacf_P2_ilht.xvg')
C2C12_rotacf_450 = np.loadtxt('450/rotacf_P2_ilht.xvg')
C2C12_rotacf_475 = np.loadtxt('475/rotacf_P2_ilht.xvg')
C2C12_rotacf_500 = np.loadtxt('500/rotacf_P2_ilht.xvg')

ax21.plot(C2C12_rotacf_300[:,0], C2C12_rotacf_300[:,1], '-', label="300", color=cmap(0))
ax21.plot(C2C12_rotacf_325[:,0], C2C12_rotacf_325[:,1], '-', label="325", color=cmap(1))
ax21.plot(C2C12_rotacf_350[:,0], C2C12_rotacf_350[:,1], '-', label="350", color=cmap(2))
ax21.plot(C2C12_rotacf_375[:,0], C2C12_rotacf_375[:,1], '-', label="375", color=cmap(3))
ax21.plot(C2C12_rotacf_400[:,0], C2C12_rotacf_400[:,1], '-', label="400", color=cmap(4))
ax21.plot(C2C12_rotacf_425[:,0], C2C12_rotacf_425[:,1], '-', label="425", color=cmap(5))
ax21.plot(C2C12_rotacf_450[:,0], C2C12_rotacf_450[:,1], '-', label="450", color=cmap(6))
ax21.plot(C2C12_rotacf_475[:,0], C2C12_rotacf_475[:,1], '-', label="475", color=cmap(7))
ax21.plot(C2C12_rotacf_500[:,0], C2C12_rotacf_500[:,1], '-', label="500", color=cmap(8))

ax21.set_xlim(0.1, 10000)
ax21.set_ylim(0, 1)
#ax21.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
#ax21.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax21.text(0.05, 0.11,'(b)',ha='center',va='center',size=16,color='black',transform=ax21.transAxes)
#ax21.leg=plt.legend(loc='upper right', fontsize=15, labelspacing=0.2, frameon=False, shadow=False, fancybox=False) 
#ax21.leg.get_frame().set_alpha(0.2)
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
ax31=plt.subplot2grid((3,9), (2,0), colspan=8, rowspan=1)
ax31.set_xlabel("Time (ps)",                             fontsize=13)
ax31.tick_params(axis="x", direction="in", which='both', labelsize=13)
ax31.tick_params(axis="y", direction="in",               labelsize=13)

plt.xscale('log')

C2C20_rotacf_300 = np.loadtxt('300/rotacf_P2_ilht.xvg')
C2C20_rotacf_325 = np.loadtxt('325/rotacf_P2_ilht.xvg')
C2C20_rotacf_350 = np.loadtxt('350/rotacf_P2_ilht.xvg')
C2C20_rotacf_375 = np.loadtxt('375/rotacf_P2_ilht.xvg')
C2C20_rotacf_400 = np.loadtxt('400/rotacf_P2_ilht.xvg')
C2C20_rotacf_425 = np.loadtxt('425/rotacf_P2_ilht.xvg')
C2C20_rotacf_450 = np.loadtxt('450/rotacf_P2_ilht.xvg')
C2C20_rotacf_475 = np.loadtxt('475/rotacf_P2_ilht.xvg')
C2C20_rotacf_500 = np.loadtxt('500/rotacf_P2_ilht.xvg')

ax31.plot(C2C20_rotacf_300[:,0], C2C20_rotacf_300[:,1], '-', label="300", color=cmap(0))
ax31.plot(C2C20_rotacf_325[:,0], C2C20_rotacf_325[:,1], '-', label="325", color=cmap(1))
ax31.plot(C2C20_rotacf_350[:,0], C2C20_rotacf_350[:,1], '-', label="350", color=cmap(2))
ax31.plot(C2C20_rotacf_375[:,0], C2C20_rotacf_375[:,1], '-', label="375", color=cmap(3))
ax31.plot(C2C20_rotacf_400[:,0], C2C20_rotacf_400[:,1], '-', label="400", color=cmap(4))
ax31.plot(C2C20_rotacf_425[:,0], C2C20_rotacf_425[:,1], '-', label="425", color=cmap(5))
ax31.plot(C2C20_rotacf_450[:,0], C2C20_rotacf_450[:,1], '-', label="450", color=cmap(6))
ax31.plot(C2C20_rotacf_475[:,0], C2C20_rotacf_475[:,1], '-', label="475", color=cmap(7))
ax31.plot(C2C20_rotacf_500[:,0], C2C20_rotacf_500[:,1], '-', label="500", color=cmap(8))

ax31.set_xlim(0.1, 10000)
ax31.set_ylim(0, 1)
#ax31.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
#ax31.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax31.text(0.05, 0.11,'(c)',ha='center',va='center',size=16,color='black',transform=ax31.transAxes)
#ax31.leg=plt.legend(loc='upper right', fontsize=15, labelspacing=0.2, frameon=False, shadow=False, fancybox=False) 
#ax31.leg.get_frame().set_alpha(0.2)
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
ax4x=plt.subplot2grid((3,9), (0,8), rowspan=9)
axins=inset_axes(ax4x, width="25%", height="63%", loc='lower left',
                   bbox_to_anchor=(1.05, 0.225, 1, 1),
                   bbox_transform=ax4x.transAxes, borderpad=0)
im=ax4x.imshow([[300, 400], [400, 500]])
fig.colorbar(im, cax=axins, ticks=[300, 325, 350, 375, 400, 425, 450, 475, 500])
axins.tick_params(direction="in", labelsize=13)
ax4x.remove()
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
plt.tight_layout()
plt.savefig("rotacf_sample.pdf", dpi=300)
plt.savefig("rotacf_sample.png", dpi=300)
plt.show()
plt.clf()
plt.close('all')
#-----------------------------------------------------------------------
