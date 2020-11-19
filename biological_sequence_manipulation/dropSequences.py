"""
================
dropSequences.py
================

PT-BR:
Script para remover sequencias com nucleotideos ambiguos (R, S, W, Y, ...) ou indeterminados (N) do arquivo multifasta. Eh necessario indicar uma porcentagem (float) no argumento --limit para ser o limite maximo de nucleotideos indeterminados permitidos na sequencia.
Sao fornecidos dois arquivos fasta de saida: um com as sequencias de "boa qualidade" e outro com as sequencias descartadas (sequencias que possuem nucleotideos ambiguos ou indeterminados compondo porcentagem supeior aquela indicada por --limit).

EN:
Script to drop sequences with  indetermined or ambiguous nucleotides (N, R, S, W, Y, ...) from a multifasta file. It is necessary to indicate a percentage (float) in the --limit argument to be the maximum limit of indeterminate nucleotides allowed in the sequence.
Two output files are provided: one with the "good quality" sequences and the other with the discarded sequences (sequences that have ambiguous or indeterminate nucleotides making up a percentage higher than indicated by --limit).

Type dropSequences.py -help to see the usage information.

"""

print(__doc__)

usage = "\tdropSequences.py\n\n\
	Script para remover sequencias com nucleotideos ambiguos (R, S, W, Y, ...) ou indeterminados (N) do arquivo multifasta.\n\n\
	***Argumentos obrigatorios:\n\
	--sequences <arquivo> (arquivo fasta para ser analisado)\n\
	--limit <float> (valor indicando o limite máximo de nucleotideos ambiguos/indeterminados permitidos na sequencia. Acima desse valor, a sequencia sera removida. Coloque 0 para remover sequencias com 1 ou mais desses nucleotideos)\n\n"

import sys, formataFasta_17Mar2020

def analizeSequences(fasta, limite):

	fasta_formatado = formataFasta_17Mar2020.SeqOneLine(fasta)

	fasta_formatado_selecionadas = {} # contera as sequencias que possuem somente ACTG
	fasta_formatado_descartadas = {} # contera as sequencias que possuem nucleotideos ambiguos ou indeterminados, segundo o limite permitido em '--limit'

	for chave in fasta_formatado.keys():
		sequencia = list(fasta_formatado[chave].upper())
		nucl_ambiguos_indeterminados = 0
		for nucl in sequencia: # percorrer a sequenica em busca de nucleotideos ambiguos ou indeterminados.
			if nucl != 'A':
				if nucl != 'C':
					if nucl != 'T':
						if nucl != 'G':
							nucl_ambiguos_indeterminados += 1
		# testar se a sequencia atual teve mais nucleotideos ambiguos/ideterminados do que o limite permitido
		percentegem_ambiguos_indeterminados = float(nucl_ambiguos_indeterminados/len(sequencia))
		if percentegem_ambiguos_indeterminados <= limite:
			fasta_formatado_selecionadas[chave] = ''.join(sequencia)
		else:
			fasta_formatado_descartadas[chave] = ''.join(sequencia)


	return fasta_formatado_selecionadas, fasta_formatado_descartadas

if __name__ == '__main__': #esta utilizando o script de forma independente (sem ser como modulo de outro).
	import sys
	if len(sys.argv) < 5:
		print(usage)
		exit()

	else:	
		if '--sequences' in sys.argv:
			argSEQS = sys.argv.index('--sequences')+1
			file = open(sys.argv[argSEQS],'r')
			fasta = file.readlines()
			file.close()
		else:
			print(usage)
			exit()

		if '--limit' in sys.argv:
			argLIMIT = sys.argv.index('--limit')+1
			limite = float(sys.argv[argLIMIT])
		else:
			print(usage)
			exit()

		resultado_selecionadas, resultado_descartadas = analizeSequences(fasta, limite)

		# criar arquivos de saida para sequencias selecionadas
		file = open(str(sys.argv[argSEQS])+'.selected','w')
		for chave in resultado_selecionadas.keys():
			file.writelines(chave+'\n'+resultado_selecionadas[chave]+'\n')
		file.close()

		# criar arquivos de saida para sequencias descartadas.
		file = open(str(sys.argv[argSEQS])+'.discarted','w')
		for chave in resultado_descartadas.keys():
			file.writelines(chave+'\n'+resultado_descartadas[chave]+'\n')
		file.close()
