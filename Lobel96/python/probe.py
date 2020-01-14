import numpy as np
import matplotlib.pyplot as plt

class Probe:
    
    signal = np.array([])
    position = (float(),float())
    fourier = np.array([],dtype=complex)
    save_signal = bool()
    
    def __init__(self,position,number_frequencies,save_signal=False):
        self.position = position
        self.fourier = np.zeros(number_frequencies,dtype=complex)
        self.save_signal = save_signal
        
    def allocate_signal(self,time):
        self.signal = np.zeros(time.size)
        self.idx = 0
        
    def append(self,value):
        self.signal[self.idx] = value
        self.idx +=1
        
    def set_field_freq(self,value,frequency_indx=None):
        if frequency_indx is None:
            self.fourier = np.copy(value)
        else:
            self.fourier[frequency_indx] = value

    def get_field_freq(self,frequency_indx=None):
        if frequency_indx is None:
            return np.copy(self.fourier)
        else:
            return self.fourier[frequency_indx]
    
    def plot(self,time=None,filename=None,fileformat=None):
        if time is None:
            plt.plot(self.signal)
            plt.xlabel('Iterations')
        else:
            plt.plot(time,self.signal)
            plt.xlabel('Time [sec]')
        plt.ylabel('Electric Field Intensity [V/m]')
        plt.title('Probe at x = %.2e' %self.position[0] 
                  + ' [m], y = %.2e' % self.position[1] + ' [m]')
        plt.grid()
        if filename is None:
            plt.show()
        else:
            if fileformat is None:
                plt.savefig(filename,'eps')
            else:
                plt.savefig(filename,fileformat)
        plt.close()
        
def get_rx_array(Rs,dtheta,number_frequencies,save_signal=False):
    theta = np.arange(0, 360, dtheta)
    Rx = []
    for i in range(theta.size):
        pos_x = Rs*np.cos(np.radians(theta[i]))
        pos_y = Rs*np.sin(np.radians(theta[i]))
        Rx.append(Probe((pos_x,pos_y),number_frequencies,save_signal))
    return Rx