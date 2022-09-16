import csv
import numpy as np


class StreamingInstabilityData:

    def __init__(self, rho_ice, rho_sil, file_path="./data/si-data"):
        self.rho_ice = rho_ice
        self.rho_sil = rho_sil
        self.file_path = file_path
