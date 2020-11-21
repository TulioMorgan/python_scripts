"""
================
formataFasta.py
================

PT-BR:
Formata um arquivo multifasta, coloca cada sequencia em apenas uma linha (remove quebra de linhas).
Retorna um dicionario, onde a chave é o header da sequencia e o valor é a própria sequencia.

EN:
Formats a multifasta file, placing each sequence in only one line (removes line breaks).
Returns a dictionary, where the key is the sequence header and the value is the sequence.
"""

print(__doc__)

def SeqOneLine(arq_fasta):
	import re
	dicionario = {}
	fasta_formatado = {}
	salva_gene = False
	contig = []

	for seq in arq_fasta: # vamos buscar o gene dentro do genoma. Esse for eh para recuperar o contig (sem quebra de linha) onde esta o gene. Se ja temos esse contig salvo, nao preciso entrar nesse for.
		if seq == '\n':
			pass
		seq = seq.rstrip()
		notSeq = re.findall('[:_.-]|\(|\)',seq) # nao permite salvar linhas que tenham caracteres estranhos a sequencias de nucleotideos ou proteinas.
		if seq.startswith('>'): # para encontrar o identificador do contig no arquivo de genoma desformatado (com quebra de linhas dentro dos contigs)
	
			if salva_gene:
				contig = ''.join(contig)
				dicionario[nome_contig] = contig
				contig = []

			nome_contig = seq # salva o nome desse contig. Nao quero novo loop for no genoma desformatado para recuperar o mesmo contig.
			#print(nome_contig)
			#nome_contig = '>'+(''.join(re.findall(r'(?<=>)[a-zA-Z0-9._-|]+(?<!:| )', seq))) # adicionei essa expressao regular para recuperar apenas a parte incial dos nomes das sequencias. Fiz isso pois estava dando problema com o alinhador needle. Estava pegando parte do ID da sequencia diferente do que eu tinha configurado no arquivo de ids para filtaragem (com o script ExtraiSeqFastaIDs.py)
			#print(nome_contig)
			salva_gene = True

		elif salva_gene and not notSeq:
			contig.append(seq)
		
	else:
		if salva_gene and not notSeq:
			contig = ''.join(contig)
			dicionario[nome_contig] = contig

	return dicionario


if __name__ == '__main__':
	import sys
	if len(sys.argv) == 1 or sys.argv[1] == '-help' or sys.argv[1] == '--help' or sys.argv[1] == '-h' or sys.argv[1] == '--help':
		print("\nUSAGE: python3 FormataFasta.py <arquivo fasta desformatado>\n")
		sys.exit()

	elif sys.argv[1] != '-help' or sys.argv[1] != '--help' or sys.argv[1] != '-h' or sys.argv[1] != '--help':
		file = open(sys.argv[1],'r')
		arquivo = file.readlines()
		file.close()

		fasta_formatado = SeqOneLine(arquivo)
		#print (fasta_formatado)
		for chave in fasta_formatado.keys():
			print(chave+'@'+fasta_formatado[chave])
