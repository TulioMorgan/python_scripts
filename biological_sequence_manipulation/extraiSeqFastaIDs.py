"""
====================
extraiSeqFastaIDs.py
====================

PT-BR:
Script para recuperar sequencias fasta de um arquivo multifasta, de acordo com uma lista de IDs fornecidos

EN:
Script to retrieve fasta sequences from a multifasta file, according to a list of IDs.

usage: python3 extraiSeqFastaIDs.py <IDs list> <multifasta file>
"""

print(__doc__)

import re, sys, formataFasta

def ExtraiSeqFastaIDs(lista_ids, fasta_formatado):
	
	fasta_selected = {}
	final_list = []

	### --- corrigir lista de ids --- ###
	if isinstance(lista_ids,str): #testa se a variavel 'lista' eh uma string, ou seja, provavelmente existe apenas 1 ID. Nao vou quebrar a string (usando split como tinha feito antes). Se qusier submeter mais de 1 ID, colocar apenas 1 em cada linha.
		final_list = [lista_ids] # final_list deve ser uma lista

		if '' in final_list: #remove elemento vazio que pode ser criado pelo 're.split'
			final_list.remove('')

	else:
		#print('## Cada ID da lista foi considerado como estando em uma linha do arquivo:', str(sys.argv[1]))
		final_list = lista_ids
	#### ---- recueperar sequencias fasta de acordo com IDs ----####
	#print('## IDs extraidos')
	for ID in final_list:
		ID = ID.rstrip()
		#print(ID)
		if ID == '':
			continue
		found = False
		for chave in fasta_formatado.keys():
			if ('>'+ID) == chave or ID == chave:
				fasta_selected[chave] = fasta_formatado[chave]
				found = True
				break
		if found:
			continue
		else: # se nao encontrou o ID completo, vamos quebrar o ID do aquivo fasta para ver se encontramos. Talvez o IDs da lista foram obtidos de um arquivo de saida do BLAST que faz truncagem do ID no primeiro espaco em branco.
			for chave in fasta_formatado.keys():
				chaveSplit = re.split(' ',chave) # quabra no espaco em branco.
				for palavra in chaveSplit: # percorre cada palavra do identificador contido na lista chaveSplit.
					if ('>'+ID) == palavra or ID == palavra:
						fasta_selected[chave] = fasta_formatado[chave]
						found = True
						break
		if not found:
			print('Nao foi possivel detectar o ID', ID, 'no arquivo fasta:', str(sys.argv[2]))

	return fasta_selected #retorna a variavel final criada dentro da funcao 

if __name__ == '__main__': #esta utilizando o script de forma independente (sem ser como modulo de outro).
	import sys
	if len(sys.argv) != 3 or sys.argv[1] == '-help' or sys.argv[1] == '--help' or sys.argv[1] == '-h' or sys.argv[1] == '--help':
		print("usage: python3 extraiSeqFastaIDs.py <IDs list> <multifasta file>")
		sys.exit()

	else:	
		file = open(sys.argv[2],'r') #multifasta
		arquivo = file.readlines()
		file.close()

		file = open(sys.argv[1],'r') #ids
		lista_ids = file.readlines()
		file.close()

		fasta_formatado = formataFasta.SeqOneLine(arquivo) # a variavel arquivo pode ser inserida nesse programa

		sequencias_selecionadas = extraiSeqFastaIDs(lista_ids, fasta_formatado) #chama a funcao acima que recuepra sequencias fasta de acordo com os ids.
		
		# definir nome_output
		if '/' in str(sys.argv[1]):
			nome1 = str(sys.argv[1]).split('/')[-1]
		else:
			nome1 = str(sys.argv[1])
		if '/' in str(sys.argv[2]):
			nome2 = str(sys.argv[2]).split('/')[-1]
		else:
			nome2 = str(sys.argv[2])

		# salvar arquivo de saida.	
		file = open(nome1+'_vs_'+nome2,'w')
		for chave in sequencias_selecionadas.keys():
			file.writelines(chave+'\n'+sequencias_selecionadas[chave]+'\n')
		file.close()

