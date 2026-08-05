import numpy as np
import scipy
import matplotlib.pyplot as plt
import seaborn as sns


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

data = scipy.io.loadmat('c2h4.mat')
globals().update((k, v) for k, v in data.items())


cm = sns.color_palette('Paired',6)


fig, axes = plt.subplots(1, 3, figsize=(12, 4))
plt.tight_layout()


axes[0].plot(x[0,:], u[0,:],'x-',color=cm[4], lw=1)
axes[0].set_xlabel(r'roll configuration $\theta$', labelpad=2)
axes[0].set_ylabel(r'heave configuration $z$', labelpad=-2)
axes[0].set_xlim((-0.4,0.85))
axes[0].set_ylim((-0.3,0.7))


axes[1].plot(u[0,:],v[0,:],'x-', color=cm[4], lw=1)
axes[1].set_xlim((-0.3,0.7))
axes[1].set_ylim((-0.4,0.6))
axes[1].set_xlabel(r'heave configuration',labelpad=2)
axes[1].set_ylabel(r'heave momentum',labelpad=-2)



axes[2].plot(x[0,:],y[0,:],'x-',color=cm[4], lw=1)
axes[2].set_xlim((-0.4,0.85))
axes[2].set_ylim((-0.8,1.0))
axes[2].set_xlabel(r'roll configuration',labelpad=2)
axes[2].set_ylabel(r'roll momentum',labelpad=-2)


plt.tight_layout(pad=0.1)
plt.savefig('c2h4.pdf')

