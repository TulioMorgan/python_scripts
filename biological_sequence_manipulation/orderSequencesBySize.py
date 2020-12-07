"""
=========================
orderSequencesBySize.py
=========================

PT-BR:
Script para ordernar sequencias em um arquivo multifasta de acordo com seus tamanhos, da maior para a menor.

EN:
Script to order sequences in a multifasta file according to their sizes, from highest to lowest.

Type python3 orderSequencesByLength.py -help to see the instructions.
"""

print(__doc__)

import sys
from collections import Counter

usage = "\norderSequencesBySize.py\n\n\
		Script para ordernar sequencias em um arquivo multifasta de acordo com seus tamanhos.\n\n\
		*** Argumento obrigatório\n\
		--input <arquivo multifasta> (arquivo multifasta)\n\n\
		*** Argumento opcional\n\
		--reverse (indica para reverter a ordem de classificacao, da menor para a maior)\n\
		**** Um arquivo de saida eh gerado com o nome <input>_filtrado_<comprimento> contendo apenas as sequencas que atingiram o cutoff determinado no argumento --comprimento\n"

def orderBySize(arquivo_fasta, reverse):
	
	import formataFasta
	fasta_formatado = formataFasta.SeqOneLine(arquivo_fasta)

	dict_lengths = {} # chave = id das sequencias. Valor = comprimento das sequencias.
	for chave in fasta_formatado.keys():
		dict_lengths[chave] = len(fasta_formatado[chave])

	dict_lengths_ordenado = sorted(dict_lengths.items(), key=lambda kv: kv[1]) # ordena o dicionario de acordo com o valor de cada chave (menor para maior)
	if not reverse:
		dict_lengths_ordenado = dict_lengths_ordenado[::-1] # reverter a ordem do dicionario, pois os pares (chave,valor) estao em ordem crescente e quero decrescente.
	fasta_formatado_ordenado = []
	for ID in dict_lengths_ordenado:
		for chave in fasta_formatado.keys():
			if ID[0] == chave:
				fasta_formatado_ordenado.append(chave+'\n'+fasta_formatado[chave])
		
	return fasta_formatado_ordenado

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print(usage)
		exit()
	else:
		if '--input' in sys.argv:
			argINPUT = sys.argv.index('--input')+1
			file = open(sys.argv[argINPUT],'r')
			fileINPUT = file.readlines()
			file.close()
		else:
			print(usage)
			exit()

		if '--reverse' in sys.argv:
			reverse = True
		else:
			reverse = False

		fasta_formatado_ordenado = orderBySize(fileINPUT, reverse)

		file = open(str(sys.argv[argINPUT])+'_ordered','w')
		for sequencia in fasta_formatado_ordenado:
			file.writelines(sequencia+'\n')
		file.close()
