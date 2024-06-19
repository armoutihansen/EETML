import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split

def _rename(df: pd.DataFrame) -> pd.DataFrame:
    print("\tRenaming columns...")
    df = df.rename(columns={
         's_x': 'disadvantage_a',
         's_y': 'disadvantage_b',
         'r_x': 'advantage_a',
         'r_y': 'advantage_b',
         'q': 'kind',
         'v': 'unkind',
         'self_x': 'piA_a',
         'other_x': 'piB_a',
         'self_y': 'piA_b',
         'other_y': 'piB_b',
         'self_z': 'piA_c',
         'other_z': 'piB_c',
         'choice_x': 'choice_a'
    })
    df['sid'],_ = pd.factorize(df['sid'])
    df['ObsID'] = np.arange(len(df))
#    df['choice_a'] = df['choice_x']
    df['choice_b'] = 1-df['choice_a']
    df['lambdaS_a'] = df['piB_a']-df['piA_a']
    df['lambdaS_b'] = df['piB_b']-df['piA_b']
    df['lambdaD_a'] = df['disadvantage_a']*(df['piB_a']-df['piA_a'])
    df['lambdaD_b'] = df['disadvantage_b']*(df['piB_b']-df['piA_b'])
    df['lambdaA_a'] = df['advantage_a']*(df['piB_a']-df['piA_a'])
    df['lambdaA_b'] = df['advantage_b']*(df['piB_b']-df['piA_b'])
    df['lambdaK_a'] = df['kind']*(df['piB_a']-df['piA_a'])
    df['lambdaK_b'] = df['kind']*(df['piB_b']-df['piA_b'])
    df['lambdaU_a'] = df['unkind']*(df['piB_a']-df['piA_a'])
    df['lambdaU_b'] = df['unkind']*(df['piB_b']-df['piA_b'])
    df['sigma_a'] = df['piA_a']
    df['sigma_b'] = df['piA_b']
    # df = df.rename(columns={
    #     's_x': 'disadvantage_a',
    #     's_y': 'disadvantage_b',
    #     'r_x': 'advantage_a',
    #     'r_y': 'advantage_b',
    #     'q': 'kind',
    #     'v': 'unkind',
    #     'self_x': 'piA_a',
    #     'other_x': 'piB_a',
    #     'self_y': 'piA_b',
    #     'other_y': 'piB_b',
    #     'self_z': 'piA_c',
    #     'other_z': 'piB_c',
    #     'choice_x': 'choice_a'
    # })
    # df['choice_b'] = 1-df['choice_a']
    df = df.fillna(0)
    print(df.columns)
    return df

def _return_long(df: pd.DataFrame) -> pd.DataFrame:
    # df['ObsID'] = np.arange(len(df))
    # df['choice_b'] = 1-df['choice_a']
    # df['lambdaS_a'] = df['piB_a']-df['piA_a']
    # df['lambdaS_b'] = df['piB_b']-df['piA_b']
    # df['lambdaD_a'] = df['disadvantage_a']*(df['piB_a']-df['piA_a'])
    # df['lambdaD_b'] = df['disadvantage_b']*(df['piB_b']-df['piA_b'])
    # df['lambdaA_a'] = df['advantage_a']*(df['piB_a']-df['piA_a'])
    # df['lambdaA_b'] = df['advantage_b']*(df['piB_b']-df['piA_b'])
    # df['lambdaK_a'] = df['kind']*(df['piB_a']-df['piA_a'])
    # df['lambdaK_b'] = df['kind']*(df['piB_b']-df['piA_b'])
    # df['lambdaU_a'] = df['unkind']*(df['piB_a']-df['piA_a'])
    # df['lambdaU_b'] = df['unkind']*(df['piB_b']-df['piA_b'])
    # df['sigma_a'] = df['piA_a']
    # df['sigma_b'] = df['piA_b']
    df_long = pd.wide_to_long(df,
                            stubnames=[
                                'choice',
                                'lambdaS',
                                'lambdaD',
                                'lambdaA',
                                'lambdaK',
                                'lambdaU',
                                'sigma'
                                ],
                            sep="_",
                            i=['sid','gid'],
                            j='alt',
                            suffix='\w+').reset_index()
    return df_long

def _remove_cols(df: pd.DataFrame) -> pd.DataFrame:
    print("\tRemoving columns...")
    df = df[['sid', 'ObsID', 'choice_a', 'choice_b',
            'lambdaS_a', 'lambdaS_b', 'lambdaD_a', 'lambdaD_b',
            'lambdaA_a', 'lambdaA_b', 'lambdaK_a', 'lambdaK_b',
            'lambdaU_a', 'lambdaU_b', 'sigma_a', 'sigma_b']]
    return df
    
def preprocess(df: pd.DataFrame,
               test_size=0.2,
               shuffle=True,
               random_state=1) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    
    # print("Preprocessing dataframe...")
    
    df = _rename(df)
    df = _remove_cols(df)
    
    train, test = train_test_split(df,
                                test_size=test_size,
                                random_state=random_state,
                                stratify=df['sid'])
    # if long:
    #     print("\tCreating long dataframe")
    #     df = _return_long(df)
    #     train = _return_long(train)
    #     test = _return_long(test)
    # print("\tDone!")

    return df, train, test