"""
=======================
removeSmallSequences.py
=======================

PT-BR:
Script para remover sequencias menores que um valor definido pelo usuario (--comprimento <int>).

EN:
Script to remove sequences with length less than a user-defined value (--comprimento <int>).

Type python3 removeSmallSequences.py -help to see the instructions.
"""

print(__doc__)

import sys, re

usage = "\nremoveSmallSequences.py\n\n\
		Script para remover sequencias menores que determinado valor (--comprimento <int>) um arquivo multifasta.\n\n\
		*** Argumento obrigatório\n\
		--input <arquivo multifasta> (arquivo de montagem do genoma)\n\n\
		--comprimento <int> (valor minimo de tamanho da sequencia fasta para ser mantida no data set)\n\
		**** Um arquivo de saida eh gerado com o nome <input>_filtrado_<comprimento> contendo apenas as sequencas que atingiram o cutoff determinado no argumento --comprimento\n"

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

	if '--comprimento' in sys.argv:
		argCOMPRIMENTO = (sys.argv.index('--comprimento')+1)
		comprimento = int(sys.argv[argCOMPRIMENTO])
	else:
		print(usage)
		exit()

multifastaFiltrado = []
header = 'None'
seq = ''
for i,linha in enumerate(fileINPUT):
	linha = linha.rstrip()
	if linha.startswith('>') and header != linha:
		#print(linha)
		if seq:
			if comprimento <= int(len(seq)): # manter a sequencia no data set, pois eh maior que cutoff.
				multifastaFiltrado.append(header+'\n')
				multifastaFiltrado.append(seq+'\n')
				
		header = linha
		seq = ''

	elif len(fileINPUT) == i+1:
		seq = seq + linha

		if seq:
			if comprimento <= int(len(seq)): # manter a sequencia no data set, pois eh maior que cutoff.
				multifastaFiltrado.append(header+'\n')
				multifastaFiltrado.append(seq+'\n')

	else:
		seq = seq + linha
file = open(sys.argv[argINPUT]+'_filtrado_'+str(comprimento), 'w')
for linha in multifastaFiltrado:
	file.writelines(linha)
file.close()

