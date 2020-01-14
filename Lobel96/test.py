import numpy as np
import copy as cp

# import fdtd2d as slv
# import domain as dm
# import boundary as bd
# import waveform as wv
# import source as sc
# import probe as pb

# dx, dy = 5e-3, 5e-3
# epsrb, sigb = 1., .0
# I, J = 50, 50
# epsr = np.ones((I,J))
# sig = np.zeros((I,J))
# sampled_frequencies = np.array([8003e6])
# wv_frequency = 800e6
# ls_x, ls_y = 4*dx, 4*dy
# magnitude = 1e1
# time_window = 5e-8
# Rs = 14e-2
# dtheta = 12
# Rx = pb.get_rx_array(Rs,dtheta,sampled_frequencies.size,save_signal=True)
# Tx = sc.Source(sc.get_tx_position(Rs,dtheta,0),wv.GaussianSignal(dx,dy,wv_frequency),magnitude)
# mymodel = slv.FDTD2D(dm.Domain(dx,dy,I,J),epsrb,sigb,sampled_frequencies,ls_x,ls_y,time_window,probes=Rx)
# epsr[np.ix_(range(round(I/4)-round(.1*I),round(I/4)+round(.1*I)),
#             range(round(J/4)-round(.1*J),round(J/4)+round(.1*J)))] = 5.0
# mymodel.run(Tx,epsr,sig)
# mymodel.plot_probes()

a = np.zeros((3,4,5))
a[:,:,0] = np.arange(12).reshape((3,4))
a[:,:,1] = np.arange(12,24).reshape((3,4))
a[:,:,2] = np.arange(24,36).reshape((3,4))
a[:,:,3] = np.arange(36,48).reshape((3,4))
a[:,:,4] = np.arange(48,60).reshape((3,4))
b = a.reshape((12,5)).copy()
print(a[:,:,0])
print(b)