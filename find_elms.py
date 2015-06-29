import numpy as np
import MDSplus as mds
import matplotlib.pyplot as plt

def find_spikes(shot,t1,t2,an=1,shift=0.2):

   tree = mds.Tree('spectroscopy', shot, 'ReadOnly')
   node = tree.getNode('\FS04F')
   tms = node.dim_of().getData().data()
   fs04f = node.getData().data()

   ind = np.where(((tms>t1)&(tms<t2)))
   u = fs04f[ind]/1e15
   t = tms[ind]
   #t is the time vector in ms
   #u is the D_alpha signal
 
   #implementing 2nd order first derivative on signal, discard 2 points at edge
   #works also by implementing 1st order derivative at the edge and keeping all the points
   du = 0.5*(u[2:]-u[0:-2])
   u = u[1:-1]
   t = t[1:-1]

   #Smoothing first derivative over N points to minimize spurious detections...not sure it should be implemented 
   #as sometimes it causes the algorithm to miss tiny spikes
   #Nc = 5
   #yc = np.ones([1,N],'int')
   #du = np.convolve(du,yc)/Nc;
   #Discarding boundaries
   #du = du[Nc-1:-(Nc-1)];

   #setting parameters
   vdu = du.std()
   M = len(du)

   #Setting automatic vs manual threshold
   # an = 1 #1 = automatic, 0 =  manual
   # manual value of threshold
   sfac = 0.5 

   #Initial step
   facdu = vdu*np.sqrt(2*np.log(M))*(an*0.92*(1-np.exp(-(M/0.05)**(0.107)))+(1-an)*sfac)
   facdu_old = 0
   ind1 = np.array([],'int')
   i = 0

   while (i<100) and (facdu!=facdu_old):
      i = i+1
      facdu_old = facdu
      dump2 = np.where(abs(du)>facdu)[0]
      dump1 = np.where(du>facdu)[0]#avoids spurious detections in negative time derivatives
      ind1=np.concatenate((ind1,dump1),axis=1)
      du[dump2] = np.nan
      dump = np.where(~np.isnan(du))[0]
      vdu = np.std(du[dump])
      M = len(dump)
      #The last term allows sfac to change as M decreases or to keep it fixed
      facdu = vdu*np.sqrt(2*np.log(M))*(an*0.92*(1-np.exp(-(M/0.05)**(0.107)))+(1-an)*sfac)
	  
   if len(ind1)==0:
      print("No ELM detected, try lowering s-fac")
   
   #sorting indices
   ind1 = np.sort(ind1)

   #Selecting only the first index, i.e. when the ELM is triggered
   #Additionally, selecting ELM at least shift ms apart...to avoid spurious detections
   shift = np.ceil(shift/np.mean(np.diff(t)))
   dump = np.where(np.diff(ind1)>shift)[0]
   dump = dump+np.ones([1,len(dump)],'int')
   ind1 = np.concatenate(([ind1[0,]],ind1[dump[0,]]),axis=1)
   plt.plot(t,u,t[ind1],u[ind1],'r*')
   plt.show()

t1 = input('Enter initial time [ms] ')
t2 = input('Enter final time [ms] ')
shot = input('Enter shot number ')
find_spikes(shot,t1,t2)
