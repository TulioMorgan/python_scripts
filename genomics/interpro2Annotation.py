"""
======================
interpro2Annotation.py
======================

PT-BR:
Esse script analisa os dominios estruturais e familias de cada proteina fornecidos de uma analise com InterProScan. Salva o identificador da proteina e a seguir todos as descricoes dos dominios (sem redundancias estritas) com as coordendas entre parenteses.
Nao salva informacoes a respeito de peptideo sinal e helice transmembrana!
Usar as opcoes -f tsv --iprlookup ao utilizar o interpro. Com isso, teremos a saida no formato tabular e os dominios IPR e descricao correspondentes para cada match (se houver a correspondencia).

EN:
This script analyzes and filter the output of InterProScan. Save the protein identifier followed by all domain descriptions (without strict redundancies) with the coordinates between parentheses.
It does not save information about signal peptide and transmembrane helix!
Use the -f tsv --iprlookup options when using the interpretation. With that, the output will be in tabular format and the corresponding IPR domains and descriptions for each match, if any.

Type python3 interpro2Annotation.py -help to see the instructions

"""

usage = "\n\tinterpro2Annotation.py\n\n\
	Script para recuperar dominios estruturais e familias de cada proteina fornecidos de uma analise com InterPro.\n\
	Salva o identificador da proteina e todos as descricoes dos dominios (sem redundancias estritas) com as coordendas entre parenteses\n\
	AVISO: A saida do InterProScan deve ser no formato tabular (-f tsv)\n\n\
	USAGE: python3 interpro2Annotation.py --input <arquivo de saida do InterProScan> --output <nome do arquivo de saida>\n"

print(__doc__)

import sys, re, operator
from collections import Counter

if len(sys.argv) < 5:
	print(usage)
	exit()
else:
	if '--input' in sys.argv:
		argInput = sys.argv.index('--input')+1
		file = open(sys.argv[argInput],'r')
		linhas = file.readlines()
		file.close()
	else:
		print(usage)
		exit()

	if '--output' in sys.argv:
		argOutput = sys.argv.index('--output')+1
		nome_output= str(sys.argv[argOutput])
	else:
		print(usage)
		exit()


####################################COMPARACOES ENTRE DOMINIOS DA MESMA PROTEINA####################################

file = open(nome_output, 'w')

secrecao = ''
transmembrana = ''
dominio = []
coordenada_refinar = []
dominios_final = []

padrao_dominios = re.compile(r'Pfam|PIRSF|PANTHER|Gene3D|PRINTS|SUPERFAMILY|ProSitePatterns|ProSiteProfiles|CDD|SMART|Hamap|ProDom|SFLD|TIGRFAM')

for i,elemento in enumerate(linhas):
	elemento = elemento.rstrip('\n')
	
	elemento_split = elemento.split('\t')

	if not i or protein_id == elemento_split[0] and not i+1 == len(linhas): # se estivermos na 1a linha do arquivo do interpro ou estivermos analisando a mesma proteina.
		protein_id = elemento_split[0] # id da proteina
	
		presente_dominios = padrao_dominios.search(elemento)

		if presente_dominios:
			if len(elemento_split) > 11:
				dominio.append(elemento_split[12]+' ('+elemento_split[11]+')\t'+'('+elemento_split[6]+'-'+elemento_split[7]+')')
			else:
				dominio.append(elemento_split[5]+' ('+elemento_split[4]+')\t'+'('+elemento_split[6]+'-'+elemento_split[7]+')')

			
	elif protein_id != elemento_split[0] or i+1 == len(linhas): # quando isso acontecer chegamos em uma nova proteina ou chegamos ao fim do arquivo. Teremos que analisar os dominios salvos da proteina anterior.
		
		if i+1 == len(linhas) and presente_dominios: # esse 'if' eh para imprimir a ultima linha do arquivo, quando esse for o caso.
			if len(elemento_split) > 11:
				dominio.append(elemento_split[12]+' ('+elemento_split[11]+')\t'+'('+elemento_split[6]+'-'+elemento_split[7]+')')
			else:
				dominio.append(elemento_split[5]+' ('+elemento_split[4]+')\t'+'('+elemento_split[6]+'-'+elemento_split[7]+')')
		
		if len(dominio) == 0:
			file.writelines(protein_id+'\n') # iniciar a escrita do arquivo com o identificador da proteina. Como existem dominios para esse ID, usar tabulacao para imprimir os dominios na frente.
		
		else:
			file.writelines(protein_id+'\t') # iniciar a escrita do arquivo com o identificador da proteina. Como existem dominios para esse ID, usar tabulacao para imprimir os dominios na frente.

		dominio.sort() # usei a funcao lista.sort(), que modifica a lista 'in-place'. É util quando nao temos interesse na lista original (desordenada). Se qusier manter a lista original, usar a função sorted(lista) e salvar o resultado em uma nova variavel (em uma nova lista).
		for k, each in enumerate(dominio):
			each_split = each.split('\t')

			if not k:
				dominio_refinar = each_split[0]
				coordenada_refinar.append(each_split[1])

			elif each_split[0] == dominio_refinar:
				#nao salva o nome do dominio, pois eh repetido. Salva apenas a coordenada.
				coordenada_refinar.append(each_split[1])

			elif each_split[0] != dominio_refinar:
				#chegou em um novo dominio dessa proteina. Devemos finalizar o dominio antigo, imprimindo seu nome e as coordenadas
				coordenadas_final = ''.join(coordenada_refinar)
				dominios_final.append(dominio_refinar+'\t'+coordenadas_final+'\t|\t')

				coordenada_refinar.clear()
				coordenada_refinar.append(each_split[1]) # salva a coordenada do proximo dominio.
				dominio_refinar = each_split[0] # alterar o valor da variavel dominio_refinar para o ID da proxima proteina.

			if k+1 == len(dominio):
				# esse elif eh para o caso de chegar no fim da lista de dominios de determinada proteina.
				coordenadas_final = ''.join(coordenada_refinar)
				dominios_final.append(dominio_refinar+'\t'+coordenadas_final+'\n')
				coordenada_refinar.clear()

		file.writelines(dominios_final)
		#limpar todas as listas para salvar os resultados da proxima proteina.
		dominio.clear()
		dominios_final.clear()

		protein_id = elemento_split[0]
		if presente_dominios:
			if len(elemento_split) > 11:
				dominio.append(elemento_split[12]+' ('+elemento_split[11]+')\t'+'('+elemento_split[6]+'-'+elemento_split[7]+')')
			else:
				dominio.append(elemento_split[5]+' ('+elemento_split[4]+')\t'+'('+elemento_split[6]+'-'+elemento_split[7]+')')
		
file.close()