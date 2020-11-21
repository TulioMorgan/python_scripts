"""
===================
concatAlignemnts.py
===================

PT-BR:
Script para concatenar alinhamentos de sequencias no formato fasta. Útil para reconstrucao filogenetica/filogenomica utilizando "super-genes".
AVISO: a versao atual deste script exige que pelo menos um dos arquivos de alinhamento contenha sequencias de todas as especies em analise.

EN:
Script to concatenate sequence alignments in fasta format. Useful for phylogenetic/phylogenomic reconstruction using "super-genes".
WARNING: the current version of this script requires that at least one of the alignment files contain sequences from all species analyzed.

Type python3 oncatAlignemnts.py -help to see the usage information.
"""

print(__doc__)

usage = '\nconcatAlignemnts.py\n\n\
Script para concatenar alinhamentos de sequencias para reconstrucao filogenetica/filogenomica.\n\
AVISO: a versao atual deste script exige que pelo menos um dos arquivos de alinhamento contenha sequencias de todas as especies em analise.\n\n\
	***Argumentos obrigatorios:\n\
		--lista_arquivos <arquivo contendo nomes dos arquivos de alinhamento> (Dica: utilize o comamndo `ls` e redirecione o output para um arquivo)\n\
		--campo <int> (campo do ID com nome da especie. Ex.: se o nome da especie estiver no primeiro campo, configure --campo 1)\n\n\
	***Argumentos opcionais:\n\
		--separador <"caracter separador de campo no ID"> (indique o separador de campo nao-default, caso necessario. Ex.: se o separador entre os campo for virgula, utilize --separador ",")\n\n\
Separadores de campo default:\npipe |\nespacos em branco (\s,\\t)\n'

def concatenaAlinhamento(lista_arquivos, campo):
	import subprocess, FormataFasta, re
	numSeqsArquivos = []
	alinhamentosConcatenados = []

	# contar quantas sequencias existem em cada arquivo. Isso eh importante, pois o arquivo guia deve conter sequencias para todas as especies em estudo.
	for arquivo in listaArquivos:
		arquivo = arquivo.rstrip()
		file = open(arquivo, 'r')
		linhas = file.readlines()
		file.close()
		contarSeq = 0
		for linha in linhas:
			contarSeq += linha.count('>')
		numSeqsArquivos.append(str(contarSeq)+':'+arquivo)
		#print(arquivo, 'possui',contarSeq, 'sequencias.')
	numSeqsArquivos.sort()
	#print(numSeqsArquivos)

	# o ultimo arquivo da lista de arquivos ordenados em ordem crescente de numero de sequencias sera o guia, pois eh provavel que possua todas as sequencias que constam nos demais (mas isso pode nao ser verdade... ver 'usage'). Vou iniciar com ele e buscar a mesma especie dentro dos outros arquivos.
	arquivoAbrir = (numSeqsArquivos[-1].split(':'))[1] # recuperar o nome do arquivo com maior numero de sequencias. 
	file = open(arquivoAbrir, 'r')
	#print(arquivoAbrir)
	arquivoGuia = file.readlines()
	file.close()

	#formatar a sequencia fasta de alinhamento
	arquivoGuiaFormatado = FormataFasta.SeqOneLine(arquivoGuia)
	
	for chave in arquivoGuiaFormatado.keys():
		
		if separador: # separador nao-default foi fornecido
			especie = chave.split(separador)
		else: # utilize separadores default
			especie = re.split(r'[|\t ]+',chave)
		especie = especie[campo-1] # queremos apenas o nome da especie.
		
		if especie.startswith('>'): # if else para nao repetir o sinal '>'
			alinhamentosConcatenados.append(especie+'\n'+arquivoGuiaFormatado[chave]) # adiciona a sequencia da especie atual na lista final de alinhamentos concatenados. Vamos procurar as sequencias dessa especie nos outros arquivos e adciona-las na lista...
		else:
			alinhamentosConcatenados.append('>'+especie+'\n'+arquivoGuiaFormatado[chave])

		for arquivo in numSeqsArquivos[:-1]: # desconsidera o ultimo arquivo da lista de arquivos ordenada, pois ele eh o arquivoGuia. Ja salvamos sua sequencia na lista 'alinhamentosConcatenados'
			arquivo = (arquivo.split(':'))[1]
			file = open(arquivo,'r') # abre o arquivo de alinhamento atual
			conteudoArquivo = file.readlines()
			file.close()
			conteudoArquivoFormatado = FormataFasta.SeqOneLine(conteudoArquivo) # formata o arquivo de alinhamento atual.

			# vamos buscar exatamente o nome da especie dentro do arquivo de sequencias atual. Usar match literal.
			find = False # controla se encontrou a sequencia da especie em determiando arquivo (conteudoArquivoFormatado). Importante para adicionar gaps no alinhamento concatenado, caso nao entre a sequencia no determinado arquivo.
			for chave in conteudoArquivoFormatado.keys():
				comprimentoGene = len(conteudoArquivoFormatado[chave]) # salva o comprimento dos genes do arquiov de alinhamento atual. Poderia ter salvo apenas 1 vez, pois sao todos iguais.
				if separador: # separador nao-default foi fornecido
					chaveSplit = chave.split(separador)
				else: # utilize separadores default
					chaveSplit = re.split(r'[|\t ]+',chave)
				#print(chaveSplit)
				for palavra in chaveSplit: # percorre cada palavra obtida com split. Devemos encontrar um match literal com a especie do arquivo guia que estamos trabalhando no momento. OBS: ja poderia ter seleciona o campo da especies e testar o if diretamente nele...
				 	if especie == palavra: # se verdadeiro encontramos a especie dentro do arquivo de sequencia atual. Vamos salvar sua sequencia na lista 'alinhamentosConcatenados'
				 		#print ('match de', especie,'do arquivo', arquivoAbrir,'com', palavra, 'do arquivo', arquivo)
				 		alinhamentosConcatenados.append(conteudoArquivoFormatado[chave])
				 		find = True
				 		break
			if not find: # se entrar nesse 'if' significa que a especie nao esta presente em 'conteudoArquivoFormatado'. Devemos adicionar gaps no lugar desse 'buraco' no alinhamento concatenado.
				gaps = '-'*comprimentoGene # cria uma string de 'gaps' do tamanho da sequencia alinhada. O comprimento dessa string de gaps deve ser o mesmo dos genes dentro do arquivo de alinhamento atual 'conteudoArquivoFormatado'.
				alinhamentosConcatenados.append(gaps)
		# apos salvar das seqeuncias da mesma especies (ou adiconar gaps, caso  ela nao exista em dado arquivo), adiconamos um '\n' para separar das proximas sequencias concatenadas.
		alinhamentosConcatenados.append('\n')
	
	alinhamentosConcatenados = ''.join(alinhamentosConcatenados)
	file = open('alinhamentos_concatenados.fasta','w')
	file.writelines(alinhamentosConcatenados)
	file.close()

if __name__ == '__main__':
	import sys
	if len(sys.argv) < 5:
		print(usage)
		exit()

	else:
		if '--lista_arquivos' in sys.argv:
			argInput= sys.argv.index('--lista_arquivos')+1
			file = open(sys.argv[argInput],'r')
			listaArquivos = file.readlines()
			file.close()
		else:
			print(usage)
			exit()

		if '--campo' in sys.argv:
			argCampo = sys.argv.index('--campo')+1
			campo = int(sys.argv[argCampo])
		else:
			print(usage)
			exit()

		if '--separador' in sys.argv:
			argSep = sys.argv.index('--separador')+1
			separador = str(sys.argv[argSep])
		else:
			separador = None

		concatenaAlinhamento(listaArquivos, campo)