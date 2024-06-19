from sklearn.base import BaseEstimator, ClassifierMixin
from functools import reduce
import numpy as np
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
        
        quietly = 1 - self.verbose
      #  y['choice_b'] = 1 - y['choice_a']
        _data = pd.concat([y, X], axis=1)
        _data = pd.wide_to_long(_data,
                                stubnames=['choice'] + list(set([x[:-2] for x in X.columns if x not in ['ObsID','sid']])),
                                sep='_',
                                i=['sid','ObsID'],
                                j='alt',
                                suffix=r'\w+').reset_index()
        choice_col = _data.pop('choice')
        _data.insert(0, 'choice', choice_col)
        print(_data.columns)
        _data.fillna(0, inplace=True)
        # _data = pd.wide_to_long(_data,
        #                         stubnames=[x for x in X.columns if x not in ['ObsID','sid']],
        #                         sep='_',
        #                         i=['sid','ObsID'],
        #                         j='alt',
        #                         suffix=r'\w+').reset_index()
        _cmd = "clogit " + reduce(lambda x, y: x + ' ' + y if y not in ['ObsID','sid','alt'] else x, _data.columns.to_list()) + ', group(ObsID) vce(cluster sid)'
        print(_cmd)
        stata.pdataframe_to_data(_data, force=True)
        stata.run(_cmd, quietly=False)
        self.coef_ = self._get_params(X, y)
        
        return self

    def _get_params(self, X: pd.DataFrame, y: pd.Series) -> pd.DataFrame:
        
      #  y['choice_b'] = 1 - y['choice_a']
        _data = pd.concat([y, X], axis=1)
        quietly = 1 - self.verbose
        _data = pd.wide_to_long(_data,
                                stubnames=['choice'] + list(set([x[:-2] for x in X.columns if x not in ['ObsID','sid']])),
                                sep='_',
                                i=['sid','ObsID'],
                                j='alt',
                                suffix=r'\w+').reset_index()
        choice_col = _data.pop('choice')
        _data.insert(0, 'choice', choice_col)
        _data.fillna(0, inplace=True)
        # _data = pd.wide_to_long(_data,
        #                         stubnames=[x for x in X.columns if x not in ['ObsID','sid']],
        #                         sep='_',
        #                         i=['sid','ObsID'],
        #                         j='alt',
        #                         suffix=r'\w+').reset_index()
        stata.pdataframe_to_data(_data, force=True)
        cols = [col for col in X.columns if col not in ['ObsID','sid']]
        temp_cols = ["("+x+": _b["+x+"]/_b[sigma])" for x in X.columns if x not in ['ObsID', 'sid', 'sigma']]
        _cmd = "nlcom (sigma:_b[sigma])"+reduce(lambda x, y: x+y, temp_cols)
        stata.run(_cmd, quietly=quietly)
        # stat_ret = stata.get_return()
        # coefs = stat_ret['r(b)'][0]
        # se = np.sqrt(np.diag(stat_ret['r(V)']))
        
        # return pd.DataFrame([coefs, se], columns=cols, index=['coefs', 'se'])
        return stata.get_return()
    
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

