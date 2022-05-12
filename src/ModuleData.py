from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import pickle

from .ModelData import ModelData

@dataclass
class ModuleData:
    """
        Base class for prediction modules
    """

    models: dict
    modelsPath: dict
    modelsNames: list

    def loadModels( self, modelsPath : dict ) -> ModuleData :
        """" Loads each model binary file into the object """

        models = { }
        modelsNames = []
        for modelName in modelsPath :
            model = ModelData( name=modelName )
            models[ modelName ] = model.loadmodel( modelsPath[ modelName ] )
            modelsNames.apend( modelName )

        return ModuleData( modelsPath=modelsPath, models=models, modelsNames=modelsNames )



