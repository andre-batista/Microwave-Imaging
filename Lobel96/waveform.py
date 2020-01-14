import numpy as np

c = 299792458

class Waveform:
    waveform = np.array([])
    def __init__(self):
        pass
    def compute_waveform(self,time):
        pass
        
class GaussianSignal(Waveform):
    
    nc = int()
    frequency_waveform = float()
    dx, dy = float(), float()
    maximum_frequency = float()
    
    def __init__(self,dx,dy,frequency_waveform,number_of_cells_per_wavelength=20):
        self.nc = number_of_cells_per_wavelength
        self.frequency_waveform = frequency_waveform
        self.dx, self.dy = dx, dy
    
    def compute_waveform(self,time):
        self.maximum_frequency = c/(self.nc*max([self.dx,self.dy]))
        tau = np.sqrt(2.3)/np.pi/self.frequency_waveform
        t_0 = 4.5*tau
        self.waveform = np.exp(-((time-t_0)/tau)**2)
        
class Sinusoidal(Waveform):
    
    frequency = float()
    
    def __init__(self,frequency_waveform):
        self.frequency = frequency_waveform
        
    def compute_waveform(self,time):
        self.waveform = np.sin(2*np.pi*self.frequency*time)

class UnitStep(Waveform):
    
    start_time_step = float()
    
    def __init__(self,start_time_step):
        self.start_time_step = start_time_step
    
    def compute_waveform(self,time):
        self.waveform = np.zeros(time.size)
        self.waveform[time>=self.start_time_step] = 1.

class DerivativeGaussian(Waveform):
    
    nc = int()
    dx, dy = float(), float()
    maximum_frequency = float()
    
    def __init__(self,dx,dy,number_of_cells_per_wavelength=20):
        self.nc = number_of_cells_per_wavelength
        self.dx, self.dy = dx, dy
    
    def compute_waveform(self,time):
        self.maximum_frequency = c/(self.nc*max([self.dx,self.dy]))
        tau = (self.nc*max([self.dx,self.dy]))/(2*c)
        t_0 = 4.5*tau
        self.waveform = -((np.sqrt(2*np.exp(1))/tau)*(time-t_0)
                          * np.exp(-((time-t_0)/tau)**2))

class CosineModulatedGaussian(Waveform):
    
    modulation_frequency = float()
    bandwith = float()
    
    def __init__(self,modulation_frequency,bandwith):
        self.modulation_frequency = modulation_frequency
        self.bandwith = bandwith
        
    def compute_waveform(self,time):
        frequency = self.modulation_frequency
        tau = .966/self.bandwith
        t_0 = 4.5*tau
        self.waveform = (np.cos(2*np.pi*frequency*(time-t_0))
                         *np.exp(-((time-t_0)/tau)**2))
