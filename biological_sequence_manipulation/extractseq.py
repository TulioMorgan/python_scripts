"""
=============
extractseq.py
=============

PT-BR:
Extrai sequencia/subsequencia de um arquivo (multi)fasta.

EN:
Extract sequence/subsequence from a (multi)fasta file.

USAGE: python3 extractseq.py <(multi)fasta file> <sequence ID>

"""

print(__doc__)

import sys,re,FormataFasta
 
if len(sys.argv) != 3:
	print("USAGE: python3 extractseq.py <arquivo multifasta formatado/desformatado> <nome do contig para extrair sequencia/subsequencia>")
	sys.exit()

else:

	file = open(sys.argv[1],'r')
	arquivo = file.readlines()
	file.close()

	print ('Corrigindo arquivo fasta...')
	fasta_formatado = FormataFasta.SeqOneLine(arquivo) #chama a funcao SeqOneLine do script 'FormataFasta.py'
	print ('Arquivo fasta corrigido.')

	counter = 0
	idsFound = []
	for chave in fasta_formatado:
		padrao = re.compile((str(sys.argv[2])).lower())
		match = padrao.search((chave).lower())

		if match:
			counter += 1
			idsFound.append(chave[1:])
	if counter > 1:
		print('\nERRO: Mais de um match foi encontrado no arquivo fasta para a busca com "'+str(sys.argv[2])+'"\nIDs encontrados: '+', '.join(idsFound)+'\nExecute novamente o script com o nome completo do ID de interesse\n')
		exit()

	found = False
	for chave in fasta_formatado:
		padrao = re.compile((str(sys.argv[2])).lower())
		match = padrao.search((chave).lower())

		if match:
			found = True
			intervalo = input("Região de "+str(sys.argv[2])+" para extrair (e.g. 31-450) [1-"+str(len(fasta_formatado[chave]))+"]: ")
			coordenadas = intervalo.split('-') 
			sequencia = fasta_formatado[chave][int(coordenadas[0])-1:int(coordenadas[1])]
				
			if sequencia and not int(coordenadas[1]) > len(fasta_formatado[chave]):
				nome = input("Nome do arquivo de saida: ")
				file = open(nome,'w')
				file.writelines(chave+'\t[Intervalo extraido: '+intervalo+']\n'+sequencia)
				file.close()
				
			else:
				print('Intervalo invalido.')
				exit()

	if not found:
		print('\nERRO: Não foi possível encontrar o ID "'+str(sys.argv[2])+'" no arquivo "'+str(sys.argv[1])+'"\n')


