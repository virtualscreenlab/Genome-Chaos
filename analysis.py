import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read data from the ".dat" file
file_path = 'wt.dat'  # Replace with the actual path to your file
data = pd.read_csv(file_path, sep='\s+', header=None, names=['x', 'y'])

# Create a scatter plot with a density estimate
sns.scatterplot(x='x', y='y', data=data, color='white', alpha=0.1)
kde = sns.kdeplot(x='x', y='y', data=data, fill=True, cmap='Reds', levels=3, alpha=1)

# Find the vertices of the triangle
triangle_vertices = [(1.6, -220), (1.75, -220), (1.6, -230)]

# Set x and y axis limits
plt.xlim(-0.5, 2.0)
plt.ylim(-280, -160)

# Add x and y axis labels
plt.xlabel('RMSD (Å)', fontsize=12)
plt.ylabel('Energy Score (REU)', fontsize=12)

# Add a title
#plt.title('Triangle and KDE Plot', fontsize=14)

# Show the plot
plt.show()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read data from the ".dat" file
file_path = 'mut.dat'  # Replace with the actual path to your file
data = pd.read_csv(file_path, sep='\s+', header=None, names=['x', 'y'])

# Create a scatter plot with a density estimate
sns.scatterplot(x='x', y='y', data=data, color='white', alpha=0.1)
kde = sns.kdeplot(x='x', y='y', data=data, fill=True, cmap='Reds', levels=3, alpha=1)

# Find the vertices of the triangle
triangle_vertices = [(1.6, -220), (1.75, -220), (1.6, -230)]

# Set x and y axis limits
plt.xlim(-0.5, 2.0)
plt.ylim(-280, -160)

# Add x and y axis labels
plt.xlabel('RMSD (Å)', fontsize=12)
plt.ylabel('Energy Score (REU)', fontsize=12)

# Add a title
#plt.title('Triangle and KDE Plot', fontsize=14)

# Show the plot
plt.show()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read data from the ".dat" file
file_path = 'wt.dat'  # Replace with the actual path to your file
data = pd.read_csv(file_path, sep='\s+', header=None, names=['x', 'y'])

# Create a scatter plot with a density estimate
sns.scatterplot(x='x', y='y', data=data, color='black', alpha=1)
#kde = sns.kdeplot(x='x', y='y', data=data, fill=True, cmap='Reds', levels=3, alpha=1)
sns.scatterplot(x='x', y='y', data=data, color='black', alpha=1, s=20, edgecolor='none')


# Find the vertices of the triangle
triangle_vertices = [(1.6, -220), (1.75, -220), (1.6, -230)]

# Set x and y axis limits
plt.xlim(0, 2.0)
plt.ylim(-300, -160)

plt.yticks(fontsize=12)
plt.xticks([0, 0.5, 1, 1.5, 2], fontsize=12)

# Add a dashed line starting from x = 0.5
plt.axvline(x=0.5, color='red', linestyle='--', linewidth=2)
plt.axhline(y=-280, color='red', linestyle='--', linewidth=2)

# Add x and y axis labels
plt.xlabel('RMSD (Å)', fontsize=12)
plt.ylabel('Energy Score (REU)', fontsize=12)

# Add a title
#plt.title('Triangle and KDE Plot', fontsize=14)

# Show the plot
plt.show()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read data from the ".dat" file
file_path = 'mut.dat'  # Replace with the actual path to your file
data = pd.read_csv(file_path, sep='\s+', header=None, names=['x', 'y'])

# Create a scatter plot with a density estimate
sns.scatterplot(x='x', y='y', data=data, color='black', alpha=1)
#kde = sns.kdeplot(x='x', y='y', data=data, fill=True, cmap='Reds', levels=3, alpha=1)
sns.scatterplot(x='x', y='y', data=data, color='black', alpha=1, s=20, edgecolor='none')


# Find the vertices of the triangle
triangle_vertices = [(1.6, -220), (1.75, -220), (1.6, -230)]

# Set x and y axis limits
plt.xlim(0, 2.0)
plt.ylim(-300, -160)

plt.yticks(fontsize=12)
plt.xticks([0, 0.5, 1, 1.5, 2], fontsize=12)

# Add a dashed line starting from x = 0.5
plt.axvline(x=0.5, color='red', linestyle='--', linewidth=2)
plt.axhline(y=-280, color='red', linestyle='--', linewidth=2)

# Add x and y axis labels
plt.xlabel('RMSD (Å)', fontsize=12)
plt.ylabel('Energy Score (REU)', fontsize=12)

# Add a title
#plt.title('Triangle and KDE Plot', fontsize=14)

# Show the plot
plt.show()
