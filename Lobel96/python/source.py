import copy as cp
import numpy as np
import waveform as wv

class Source:
    
    center_position = (float(),float())
    min_x, max_x = float(), float()
    min_y, max_y = float(), float()
    direction = ''
    magnitude = float()
    signal = wv.Waveform()
    i_s, j_s = int(), int()
    i_e, j_e = int(), int()
    
    def __init__(self, center_position, waveform, magnitude):
        self.center_position = center_position
        self.signal = cp.deepcopy(waveform)
        self.magnitude = magnitude
        
def get_tx_position(Rs,dtheta,s_indx):
    theta = np.arange(0, 360, dtheta)
    pos_x = Rs*np.cos(np.radians(theta[s_indx]))
    pos_y = Rs*np.sin(np.radians(theta[s_indx]))    
    return (pos_x,pos_y)