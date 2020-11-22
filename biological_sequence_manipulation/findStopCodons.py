"""
=================
findStopCodons.py
=================

PT-BR:
Retorna sequencias (aminoacidos ou nucleotideos) sem stop codons internos ou reporta todas as sequencias até atingir o primeiro stop codon inframe.

EN:
Returns sequences (amino acids or nucleotides) without internal stop codons or reports all sequences until reaching the first in-frame stop codon.

Type python3 findStopCodons.py -help to see the instructions

"""

print(__doc__)

usage = '\nfindStopCodons.py\n\n\
		Script para detectar stop codons in-frame em sequencias de aminoacidos ou nucleotideos\n\
		*** Argumentos obrigatorios:\n\
		--input <arquivo> (arquivo fasta de sequencias)\n\
		--seq_type <prot|nuc>\n\n\
		*** Argumentos de opcionais\n\
		--report_seq <True|False> (True: reporta todas as sequencias, removendo o stop codon e o restante da sequencia apos ele). Default: False. Retorna apenas sequencias sem stop codons inframe\n\n'

def findStopCodons(arquivoFasta):
	import re, formataFasta

	fastaFormatado = formataFasta.SeqOneLine(arquivoFasta)
	 
	if sequencia == 'prot': # vai buscar sinais de stop codons, i.e, '*' ou '.'
		file = open(str(sys.argv[argFasta])+'ids_proteinas_stop_codons_internos.txt', 'w')
		for chave in fastaFormatado.keys():

			stopCodons = re.findall(r'[\*\.]*',fastaFormatado[chave])

			if stopCodons.count('*') > 1:
				file.writelines(chave+'\n')

			elif stopCodons.count('.') > 1:
				file.writelines(chave+'\n')

		file.close()

	elif sequencia == 'nuc': # busca stop codons internos nos transcritos (TGA, TAG, TAA)
		output_sequencias = {} # dicionario de saida
		for chave in fastaFormatado.keys():
			stop_codon_interno = False
			seqUpper = fastaFormatado[chave].upper()
			codons = re.findall(r'[A-Z]{3}', seqUpper)
			save_codon = [] # lista de codons, sem stop codons internos 
			for i,codon in enumerate(codons):
				if codon == 'TGA' or codon == 'TAG' or codon == 'TAA':
					if i+1 != len(codons): # verifica se o stop codon NAO eh o ultimo na sequencia
						if report_seq: # vamos reportar a sequencia sem o stop codon. Usar isso somente na pipeline de dN/dS usando o pal2nal.pl
							print('\nArgumento --report_seq ativo. Reportar a sequencia sem o Stop codon.')
							print(chave)
							print('nao salvou o', codon, 'da posicao', codons.index(codon)+1, '\nO fim da sequencia eh o codon numero', len(codons), '\n')
							pass #nao salva o stop codon na lista de codons.
						else: # caso nao deseje reportar a sequencia com stop codon, vamos pular para a proxima sequencia.
							stop_codon_interno = True # -report_seq eh False, logo, nao reportar a sequencia que tenha stop codon interno
							break

					else: # o stop codon eh o ultimo na sequencia
						save_codon.append(codon)
				else: # nao eh um stop codon
					save_codon.append(codon)

			if not stop_codon_interno:
				output_sequencias[chave] = ''.join(save_codon)


	return output_sequencias

if __name__ == '__main__':
	import sys
	if len(sys.argv) < 3:
		print(usage)
		exit()

	else:
		if '--input' in sys.argv:
			argFasta = sys.argv.index('--input')+1
			file = open(sys.argv[argFasta],'r')
			arquivoFasta = file.readlines()
			file.close()
		else:
			print(usage)
			exit()

		if '--seq_type' in sys.argv:
			argSeq = sys.argv.index('--seq_type')+1
			sequencia = str(sys.argv[argSeq])
			if sequencia == 'prot':
				pass
			elif sequencia == 'nuc':
				pass
			else:
				print('ERRO: Argumento --seq_type inválido. Utilize --seq_type prot ou --seq_type nuc\n')
				exit()
		else:
			print(usage)
			exit()

		if '--report_seq' in sys.argv: # vai reportar a sequencia sem o stop codon
			argReport = sys.argv.index('--report_seq')+1 
			report_seq = str(sys.argv[argReport])
			if report_seq == 'True':
				report_seq = True
			elif report_seq == 'False':
				report_seq = False
			else:
				print('ERRO: Argumento --report_seq inválido. Utilize --report_seq True ou --report_seq False\n')
				exit()
		else:
			report_seq = False

		output = findStopCodons(arquivoFasta)

		file = open('output_findStopCodons.fasta', 'w')
		for chave in output.keys():
			file.writelines(chave+'\n'+output[chave]+'\n')
		file.close()