from __future__ import annotations
from dataclasses import dataclass
from .ModelData import ModelData

@dataclass
class ModuleData:
    """
        Base class for prediction modules
    """

    predictors: dict
    modelsPath: dict
    predictorsNames: list

    def loadModels( modelsPath : dict ) -> ModuleData :
        """" Loads each model binary file into the object """

        predictors = { }
        predictorNames = []
        for predictorName in modelsPath :
            predictor = ModelData( name=predictorName, model='', modelPath=modelsPath[ predictorName ], params=[])
            predictors[ predictorName ] = predictor.loadModel( modelsPath[ predictorName ] )
            predictorNames.append( predictorName )


        return ModuleData( modelsPath=modelsPath, predictors=predictors, predictorsNames=predictorNames )



