#!/usr/bin/env python
from scipy.optimize import minimize

class Optimizer:
    """
    
    """

    def __init__(self):
        pass


    def add_noise(self):
        pass

   
    @staticmethod
    def linear_combination(weights: List[float], X: List[np.ndarray]) -> int:
        w0, w1 = weights
        x0, x1 = X
        return w0*x0 + w1*x1


    def optimize_linear_combination(self, func) -> :
        minimize(func, np.array(1, 1), 

