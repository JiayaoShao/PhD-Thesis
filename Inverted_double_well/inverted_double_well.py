import numpy as np
import scipy
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D


#plt.rcParams['text.usetex'] = False
fig_width_pt = 500  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
#ratio = 1.1
ratio = 4./3                            # Sane ratio
#ratio = 2./(pylab.sqrt(5)-1.0)                # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width/ratio            # height in inches
fig_size =  [fig_width,fig_height]

sns.set_style("ticks")
params = {'backend': 'ps',
          'axes.labelsize': 10,
          'font.family': 'serif',
 #         'font.serif': 'Computer Modern Roman',
          'font.weight': 'normal',
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': False,
          'figure.figsize': fig_size}
plt.rcParams.update(params)

data = scipy.io.loadmat('inverted_doublewell.mat')
globals().update((k, v) for k, v in data.items())
T = T[0,0];xv = xv[0,0]; c = c[0,0];
def V(x,y):
    return 0.5*c*x**2-0.25*c*(x**4)/(xv**2)+0.5*y**2


cm = sns.color_palette('Paired',6)

saddlelevel = V(saddle[0,0],saddle[0,1])

fig, axes = plt.subplots(2,2); axes = axes.flatten()

#axes[0].contour(x_, u_, V(x_, u_), np.arange(0,3*saddlelevel,0.25*saddlelevel), colors='k', linewidths=0.3)
axes[0].contourf(x_,y_, V(x_,y_), np.arange(0,3*saddlelevel,0.25*saddlelevel), cmap='Blues_r')
#axes[0].plot(saddle[0,0], saddle[0,1],'o', color=cm[4])
#axes[0].plot(-saddle[0,0], saddle[0,1],'o', color=cm[4])
axes[0].plot(phix_sink[0,:], phiy_sink[0,:], color=cm[4], lw=1)
axes[0].plot(phix_saddle[0,:], phiy_saddle[0,:], color=cm[4], lw=1)
axes[0].plot(phinull[0,0], phinull[0,1],'o',markersize=3,color=cm[5])
axes[0].plot(pend[0,0], pend[0,1],'o', markersize=3,color=cm[5])
axes[0].plot(phix[0,:], phiy[0,:], color=cm[5])
axes[0].set_xlabel(r'$x$')
axes[0].set_ylabel(r'$y$')
axes[0].set_xlim((-1.2,1.2))
axes[0].set_ylim((-0.8,0.8))
axes[0].text(-0.2, 1.1, '(a)', transform=axes[0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right')


axes[1].text(-10,0.2,'linear', ha='center',
             bbox=dict(facecolor='white', edgecolor='none', pad=2.0, alpha=0.8))
axes[1].text(15,0.2,'nonlinear', ha='center',
             bbox=dict(facecolor='white', edgecolor='none', pad=2.0, alpha=0.8))
axes[1].text(40,0.2,'linear', ha='center',
             bbox=dict(facecolor='white', edgecolor='none', pad=2.0, alpha=0.8))
axes[1].plot([0,0],[-100,100], 'k--',lw=1.0)
axes[1].plot([T,T],[-100,100], 'k--',lw=1.0)
axes[1].plot(t[0,:], phiz[0,:], color=cm[3], label=r'$z$')
axes[1].plot(sinkt[0,:], phiz_sink[0,:], color=cm[2])
axes[1].plot(saddlet[:,0], phiz_saddle[0,:], color=cm[2])
axes[1].set_xlim((-20,50))
axes[1].set_ylim((-0.15,0.25))
axes[1].grid()
lgd=axes[1].legend(loc='lower left', bbox_to_anchor=(0.2, 0.00), ncol=2, frameon=True)
lgd.get_frame().set_linewidth(0.0)
axes[1].set_xlabel(r'$t$')
axes[1].text(-0.2, 1.1, '(b)', transform=axes[1].transAxes, fontsize=12, fontweight='bold', va='top', ha='right')


fig.delaxes(axes[2])
axes[2] = fig.add_subplot(2, 2, 3, projection='3d')  # reassign as 3D
surf = axes[2].plot_surface(Xrot, Yrot, Zrot, cmap='viridis', edgecolor='none')
axes[2].plot3D(phix[0,:],phiy[0,:],phiz[0,:],color=cm[5])
axes[2].plot3D(phix_sink[0,:],phiy_sink[0,:],phiz_sink[0,:],color=cm[4], lw=1)
axes[2].plot3D(phix_saddle[0,:],phiy_saddle[0,:],phiz_saddle[0,:],'--',color=cm[4], lw=1)
axes[2].view_init(elev=40, azim=30)
axes[2].set_xlabel(r'$x$')
axes[2].set_ylabel(r'$y$')
axes[2].set_zlabel(r'$z$')
axes[2].annotate('(c)', xy=(-0.25, 1.1), xycoords='axes fraction',
                 fontsize=12, fontweight='bold', va='top', ha='right')

axes[3].text(-10,0.25,'linear', ha='center',
             bbox=dict(facecolor='white', edgecolor='none', pad=2.0, alpha=0.8))
axes[3].text(15,0.25,'nonlinear', ha='center',
             bbox=dict(facecolor='white', edgecolor='none', pad=2.0, alpha=0.8))
axes[3].text(40,0.25,'linear', ha='center',
             bbox=dict(facecolor='white', edgecolor='none', pad=2.0, alpha=0.8))
axes[3].plot([0,0],[-100,100], 'k--',lw=1.0)
axes[3].plot([T,T],[-100,100], 'k--',lw=1.0)
axes[3].plot(t[0,:-1], thz[0,:], color=cm[3], label='white noise')
axes[3].plot(sinkt[0,:], etasink[0,:], color=cm[2])
axes[3].plot(saddlet[:,0], etasaddle[0,:], color=cm[2])
axes[3].set_xlim((-20,50))
axes[3].set_ylim((-0.2,0.3))
axes[3].grid()
lgd=axes[3].legend(loc='lower left', ncol=1, frameon=True)
lgd.get_frame().set_linewidth(0.0)
axes[3].set_xlabel(r'$t$')
axes[3].text(-0.2, 1.1, '(d)', transform=axes[3].transAxes, fontsize=12, fontweight='bold', va='top', ha='right')


plt.tight_layout(pad=0.1)
plt.savefig('inverted_double_well.pdf',bbox_inches="tight")

