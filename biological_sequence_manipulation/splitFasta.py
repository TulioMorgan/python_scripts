"""
=============
splitFasta.py
=============

PT-BR: divide aleatoriamente um arquivo multifasta em arquivos fasta indepedentes. O numero de sequencias em cada arquivo fasta eh definido pelo usuario.

EN: randomly divides a multifasta file into independent fasta files. The number of sequences in each file is defined by the user.

Type python3 splitFasta.py -help to see the instructions.
"""

print(__doc__)


usage = "\n	splitFasta.py\n\n\
	Script para particionar arquivos multifasta de acordo com o numero de sequencias para cada arquivo.\n\n\
	*** Argumentos obrigatórios\n\
		--input <arquivo multifasta>\n\
		--num_seqs <int> (numero de sequencias em cada arquivo)\n\
		--output <nome arquivo de saida>\n\n"

def divideFasta(arquivoFasta, numSequencias):
	
	import FormataFasta
	fasta_formatado = FormataFasta.SeqOneLine(arquivoFasta)

	# numero de sequencias no arquivo multifasta
	numSeqTotal = len(fasta_formatado.keys())

	# converter dicionario para lista
	seqList = []
	for chave in fasta_formatado.keys():
		seqList.append(chave+'\n'+fasta_formatado[chave]+'\n')

	iteracao = 1
	chunck = 1
	while iteracao <= numSeqTotal:
		
		file = open(output+'_'+str(chunck),'a')

		file.writelines(seqList[iteracao-1])

		if iteracao != 0 and (iteracao % numSequencias) == 0:
			chunck += 1
			file.close()

		iteracao += 1


if __name__ == '__main__':
	import sys
	if len(sys.argv) < 7:
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

		if '--num_seqs' in sys.argv:
			argNum = sys.argv.index('--num_seqs')+1
			numSequencias = int(str(sys.argv[argNum]))
		else:
			print(usage)
			exit()

		if '--output' in sys.argv:
			argOutput = sys.argv.index('--output')+1
			output = str(sys.argv[argOutput])
		else:
			print(usage)
			exit()

		divideFasta(arquivoFasta, numSequencias)