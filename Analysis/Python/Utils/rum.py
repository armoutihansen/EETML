from sklearn.base import BaseEstimator, ClassifierMixin
from functools import reduce
import pandas as pd
from pystata import config
import contextlib
with contextlib.redirect_stdout(None):
    config.init('se')
from pystata import stata

class RUM(BaseEstimator, ClassifierMixin):
    
    def __init__(self, utility='linear', verbose=False) -> None:
        self.verbose = verbose
        self.utility = utility

    def fit(self, X: pd.DataFrame, y: pd.Series) -> BaseEstimator:
        
        verbose = self.verbose
        _data = pd.concat([y, X], axis=1)
        _cmd = "clogit " + reduce(lambda x, y: x + ' ' + y if y not in ['ObsID','sid'] else x, _data.columns.to_list()) + ', group(ObsID) vce(cluster sid)'
        stata.pdataframe_to_data(_data, force=True)
        stata.run(_cmd, quietly=verbose)
        self.coef_ = self._get_params(X, y)
        
        return self

    def _get_params(self, X: pd.DataFrame, y: pd.Series) -> pd.DataFrame:
        
        _data = pd.concat([y, X], axis=1)
        stata.pdataframe_to_data(_data, force=True)
        cols = [col for col in X.columns if col not in ['ObsID','sid']]
        temp_cols = ["("+x+": _b["+x+"]/_b[sigma])" for x in X.columns if x not in ['ObsID', 'sid', 'sigma']]
        _cmd = "nlcom (sigma:_b[sigma])"+reduce(lambda x, y: x+y, temp_cols)
        stata.run(_cmd, quietly=True)
        table = stata.get_return()['r(table)']
        coefs = table[0]
        se = table[1]
        
        return pd.DataFrame([coefs, se], columns=cols, index=['coefs', 'se'])
    
    def predict_proba(self, X: pd.DataFrame) -> pd.Series:
        
        _data = X
        stata.pdataframe_to_data(_data, force=True)
        stata.run("predict p, pc1")
        proba = stata.pdataframe_from_data(selectvar=-1)["p"]
        
        return proba
    
    def predict(self, X: pd.DataFrame) -> pd.Series:
        
        proba = self.predict_proba(X)
        pred = proba.apply(lambda x: 1 if x > .5 else 0)
        
        return pred


