import numpy as np


class HMFNum(object):

    def __init__(self, hmasses : np.ndarray, **props) -> None:
        self.mass = hmasses
        self.props = props
        return
    
    def dndlnM(self, mbins : np.ndarray):
        hist, edges = np.histogram(self.mass, bins = mbins)
