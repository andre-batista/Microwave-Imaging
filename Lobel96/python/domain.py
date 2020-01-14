import numpy as np

class Domain:
    
    dx, dy = float(), float()
    min_x, max_x = float(), float()
    min_y, max_y = float(), float()
    size_x, size_y = float(), float()
    cell_center_coordinates_x = np.array([])
    cell_center_coordinates_y = np.array([])
    I, J = int(), int()
    
    def __init__(self, dx, dy, I, J):
        self.dx, self.dy = dx, dy
        self.I, self.J = I, J