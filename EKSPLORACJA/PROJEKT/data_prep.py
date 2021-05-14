import pandas as pd
import numpy as np

df = pd.read_csv('DATA.csv')

to_drop = ['customerID']
df.drop(to_drop, inplace=True, axis=1)

df = df.set_index('lp')
print(df.columns)