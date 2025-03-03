
import warnings
from .eco import EcoFile

class HGSEcoFileParser(EcoFile):
    def __init__(self, *args, **kwargs):
        warnings.warn('This class is deprecated and will be removed soon. '\
            'Use EcoFile instead')
        super().__init__(*args, **kwargs)

