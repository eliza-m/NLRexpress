from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from sklearn.neural_network import MLPClassifier
import pickle



@dataclass
class ModelData:
    """
        Base class for prediction models
    """

    model: MLPClassifier
    modelPath: Path
    params: dict
    name: str

    def loadModel(self, path: Path) -> ModelData :
        """" Loads model binary file into the object """

        model = pickle.load(open(path, 'rb'))
        return ModelData(name=self.name, modelPath=path, model=model, params=self.params)


    def train(self, X, Y, params, outFile:Path, printToFile=False) -> ModelData :
        """ trains model """
        model = MLPClassifier()
        for param in params:
            model.setattr(param, params[param])

        model.fit(X, Y)

        if printToFile:
            pickle.dump(model, open(outFile, 'wb'))

        return ModelData(name=self.name, modelPath=outFile, models=model, params=params)

