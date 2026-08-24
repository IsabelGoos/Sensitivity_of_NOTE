from obspy.clients.syngine import Client as SyngineClient
import obspy
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert
#import scipy.signal

c_s      = SyngineClient()
st_synth = obspy.Stream()

for k in range(0, 181):
   if (np.mod(k, 10)==0):
      print(k)
   st_synth += c_s.get_waveforms(model="ak135f_5s", receiverlatitude=0, receiverlongitude=k,sourcelatitude=0, sourcelongitude=0, sourcedepthinmeters=30000,sourcedoublecouple=[145,43,61,3.51e+24],dt="0.1", units="displacement", components="Z")

mat = np.zeros([len(st_synth)-1,len(st_synth[0].times())])
x = st_synth[0].times()/60
y = np.zeros(len(st_synth)-1)
for k in range(1-1,len(st_synth)-1):
   st = st_synth[k]-np.mean(st_synth[k])
   analytic_signal = hilbert(st.data)
   env = np.abs(analytic_signal)
   # plt.plot(st_synth[k].times(),k+st_synth[k].data/abs(np.max(st_synth[k].data)))
   # plt.plot(st_synth[k].times(),k+10*env/abs(np.max(env)),'k',lw=0.5)
   mat[k,:] = env#/abs(np.max(env))
   # plt.plot(st_synth[k].times(),np.cos(np.pi/2+k*np.pi/180/2)+env/abs(np.max(env)),'k',lw=0.5)
   y[k] = np.cos(np.pi/2+k*np.pi/180/2)

print(np.shape(st_synth))
k_plot1 = 60  # choose angle index
k_plot2 = 128
trace1 = st_synth[k_plot1]
trace2 = st_synth[k_plot2]
trace1 = trace1 - np.mean(trace1.data)
trace2 = trace2 - np.mean(trace2.data)
analytic_signal1 = hilbert(trace1.data)
analytic_signal2 = hilbert(trace2.data)
env1 = np.abs(analytic_signal1)
env2 = np.abs(analytic_signal2)
t = x  # minutes
plt.rcParams.update({'font.size': 18})
fig, ax1 = plt.subplots(figsize=(6.5, 3.), constrained_layout=True)
plt.plot(t, np.array(trace1.data), '-',  linewidth=1, label=r"$\Delta=60^\circ$")
plt.plot(t, np.array(trace2.data), '--', linewidth=1, label=r"$\Delta=128.3^\circ$")
plt.xlim(0, 60)
plt.xlabel("Time [min]")
plt.ylabel("Amplitude")
plt.legend(fontsize=16)
#fig.text(0.75, 0.35, r"$\Delta = 60^\circ$")
plt.savefig('filename3.png', dpi=300, bbox_inches='tight', transparent=True)
plt.show()

   
X,Y = np.meshgrid(x,y)

plt.rcParams.update({'font.size': 18})
fig, ax1 = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)

h = ax1.pcolor(X,Y,mat*5000/mat.max(),vmin=0,vmax=1,cmap="cividis")
ax1.set_xlim([0,60])
ax1.set_ylim([-1,0])
ax1.set_xlabel(r'Time [min]')
ax1.set_ylabel(r'$\cos\theta_z$')

ax2 = ax1.twinx()
h = ax2.pcolor(X,Y,mat*5000/mat.max(),vmin=0,vmax=1,cmap="cividis")
vec = np.flipud([0,np.cos(105*np.pi/180),np.cos(120*np.pi/180),np.cos(135*np.pi/180),np.cos(150*np.pi/180),-1])
ax2.set_yticks(vec)  # Aligner les ticks avec la pression
ax2.set_yticklabels([f"{int(d)}" for d in [180,120,90,60,30,0]])  # Affichage en km
ax2.set_ylabel(r'$\Delta$ [$^\circ$]')
#ax2.plot([0,100],[-0.84,-0.84],'--w',lw=2)

cbar = fig.colorbar(h, ax=ax1, label="Seismic energy [arb. units]", orientation='horizontal', location='top') 
cbar.set_label("Seismic energy [arb. units]", labelpad=10)

plt.savefig('filename2.png', dpi=300, bbox_inches='tight', transparent=True)

for i in range(len(y)):
    print(i, y[i], (np.arccos(y[i]) * 2.0 - np.pi) * (180.0 / np.pi))

fig, ax1 = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
plt.plot(x, mat[60, :]*5000/mat[60, :].max())
#plt.yscale('log')
plt.xlabel(r'Time [min]')
plt.ylabel("Seismic energy [arb. units]")
plt.savefig('filename1.png', dpi=300, bbox_inches='tight', transparent=True)
plt.show()

