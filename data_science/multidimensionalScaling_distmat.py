"""
==================================
multidimensionalScaling_distmat.py
==================================

PT-BR:
Script para calcular escalonamento multidimensional a partir de uma matriz upper triangular obtida pelo algoritmo distmat do EMBOSS.
Antes de executar esse script, deve-se editar manualmente o arquivo distmat, removendo as linhas inciais (8 linhas) e coluna inicial e final que ficam em branco. Nomear a coluna com os IDs das sequncias (ultima coluna) com o header 'labels'
Esse script necesita do modulo lower_triangular_matrix.py (preenche a matriz lower triangular, criando uma matriz simétrica)

EN:
Script to calculate multidimensional scaling from a triangular upper matrix obtained using the EMBOSS distmat algorithm.
Before executing this script, you must manually edit the distmat file, removing the starting lines (8 lines) and the starting and ending columns that are blank. Name the column with the sequence IDs (last column) with the 'labels' header
This script needs the lower_triangular_matrix.py module (fills the lower triangular matrix, creating a symmetric matrix)

"""

print(__doc__)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from sklearn.preprocessing import MinMaxScaler

# Carregar os dados. Costumo executar o distmat (EMBOSS) usando diferentes modelos de correcao de distancias para sequencias mais divergentes 
#dados = pd.read_csv('matriz_teste.csv', sep=',')
dados = pd.read_csv('mafft.fasta_trim.distmat_JinNei-gamma', sep='\t')
#dados = pd.read_csv('mafft.fasta_trim.distmat_TajimaNei', sep=',')
#dados = pd.read_csv('mafft.fasta_trim.distmat_K2P', sep=',')
#dados = pd.read_csv('mafft.fasta_trim.distmat_JC', sep=',')

import lower_triangular_matrix
dados, labels = lower_triangular_matrix.lower_tringular(dados)

# transforma a lista de classes em vetor.
vetor_classes = np.array(labels)

#suppress exponential notation, define an appropriate float formatter. Tive que fazer isso porque o numpy estava convertendo os numeros float para notacao cientifica
np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.3f}'.format})
dados_log = np.log10(dados)

#Then, we apply the MDS procedure to get a 2-dimensional dataset. The random_state is set in order to make every plot reproducible.
mds = MDS(n_components=2, random_state=0)
dados_2d = mds.fit_transform(dados_log) # coordenadas escalonadas para 2 dimensoes

#Finally, we can plot the dataset.
plt.rcParams['figure.figsize'] = [7, 7]
plt.rc('font', size=14) # 'rc' stands for 'runtime configuration' and contains the default styles for every plot you create.
plt.scatter(dados_2d[:,0], dados_2d[:,1])
plt.show()

# Essa fucao gera um MDS classico, igual da funcao cmdscale() do R. O metodo anteior (MDS) faz um escalonamento ligeiramente diferente.
# carregar o modulo cmdScale.py
import cmdScale
dados_2d, eigenvalues = cmdScale.cmdscale(dados_log)

plt.rcParams['figure.figsize'] = [7, 7]
plt.rc('font', size=14) # 'rc' stands for 'runtime configuration' and contains the default styles for every plot you create.
plt.scatter(dados_2d[:,0], dados_2d[:,1])
plt.show()
