#print('|------------------------------------------------------|')
#print('| *          |                             |   *       |')
#print('|      *     |         3D PCA plot         |      *    |')
#print('|            |                             |  *        |')
#print('|  *         |   Siavash Salek Ardestani   |       *   |')
#print('|          * |            Contact:         |    *      |')
#print('|   *    *   |     siasia6650@gmail.com    |*          |')
#print('|      *     |                             |       *   |')
#print('|------------------------------------------------------|')

import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse

# Argument parser setup
P = argparse.ArgumentParser()
P.add_argument('--evec', help='eigenvec file', required=True)
P.add_argument('--eval', help='eigenval file', required=True)
P.add_argument('--s', help='Size of scatter points', type=int, required=False, default=10)
P.add_argument('--x', help='Width of figure', type=int, required=False, default=10)
P.add_argument('--y', help='Height of figure', type=int, required=False, default=10)
P.add_argument('--o', help='Output file prefix', type=str, required=True)
args = P.parse_args()

# Load data
data = pd.read_csv(args.evec, delimiter=r"\s+", header=None)
eigenval = pd.read_csv(args.eval, delimiter=r"\s+", header=None, nrows=3).round(2)
df = pd.DataFrame([('PC1'), ('PC2'), ('PC3')])
df2 = pd.DataFrame([['Breed'], ['ID']])
eigenval = (df + '(' + eigenval.astype(str) + '%)').astype(str)
keep_col1 = [0, 1, 2, 3, 4]
PCA = data[keep_col1]
header = df2.append(eigenval, ignore_index=True)
PCA = PCA.rename(columns=header[0])

# Define custom color and shape palettes
color_palette = {
    "PR_CR": "#01665E", "PR_Ph": "#80CDC1", "PR_PR325": "#35978F",
    "MK": "black", "GD": "#4D9221", "GD2": "#7FBC41"
}

shape_palette = {
    "MK": "+", "GD": "o", "GD2": "o",
    "PR_CR": "s", "PR_Ph": "s", "PR_PR325": "s"
}

# Create 3D PCA plot
fig = plt.figure(figsize=(args.x, args.y))
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)

# Ensure Breed is categorical
PCA['Breed'] = pd.Categorical(PCA['Breed'])
labels = np.unique(PCA['Breed'])

for label in labels:
    PCA1 = PCA[PCA['Breed'] == label]
    color = color_palette.get(label, "gray")  # Default color if not found
    marker = shape_palette.get(label, "o")   # Default shape if not found

    ax.scatter(
        PCA1.iloc[:, 2:3], PCA1.iloc[:, 4:5], PCA1.iloc[:, 3:4], 
        s=args.s, color=color, alpha=1, edgecolor='k', linewidth=0.5, label=label, marker=marker
    )

# Set axis labels
ax.set_xlabel(PCA.columns[2])
ax.set_ylabel(PCA.columns[4])
ax.set_zlabel(PCA.columns[3])

# Add legend and save the plot
#plt.legend()
ax.view_init(20, 65)
plt.savefig(args.o + '3dPCA.pdf', dpi=300, format="pdf")
plt.show()
