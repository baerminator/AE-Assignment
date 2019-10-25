import pandas as pd
import matplotlib.pyplot as plt
df  = pd.read_csv("Hull.csv")
df2 = pd.read_csv("points.csv")

df2.plot(kind='scatter',x='x',y='y') # scatter plot
df.plot(kind='scatter',x='x',y='y') # scatter plot
plt.show()