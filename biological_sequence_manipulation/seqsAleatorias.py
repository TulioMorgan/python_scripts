"""
=================
seqsAleatorias.py
=================

PT-BR:
Script para selecionar sequencias aleatórias de um arquivo multifasta.

EN:
Script to select random sequences from a multifasta file.

usage: python3 seqsAleatorias.py <multifasta file> <number of sequences to retrieve>

"""

print(__doc__)

import sys

file = open(sys.argv[1],'r')
proteinas = file.readlines()
file.close()

numero_proteinas = int(str(sys.argv[2]))

# formata o arquivo fasta. Cria um dicionario, ou seja, cada vez que o programa roda, as proteinas ficam aleatorizadas.
import formataFasta
fasta_formatado = formataFasta.SeqOneLine(proteinas)

file = open('proteinas_selecionadas.fasta','w')
for chave in fasta_formatado.keys():
	seqSelecionadas = chave+'\n'+fasta_formatado[chave]+'\n'
	file.writelines(seqSelecionadas)

	numero_proteinas -= 1
	if numero_proteinas > 0:
		continue
	else:
		break

file.close()


