"""
=============================================
pairwiseGlobalAlignment_sequenceIdentities.py
=============================================

PT-BR:
Script para realizar alinhamento global par-a-par. O input consiste de um arquivo multifasta. O argumento --identity configura o limite mínimo de identidade para que duas sequencias sejam consideradas redundantes. As saidas consistem de arquivo tabulares indicando a identidade de sequencias entre cada par e um arquivo fasta de sequencias não-redundantes. O programa utiliza EMBOSS (needle) para realizar os alinhamentos.

EN:
Script to perform global alignments. The input consists of a multifasta file. The --identity argument sets the minimum identity for two strings be considered redundant. The outputs consist of tabular files indicating the sequence identity between each pair and a fasta file with non-redundant sequences. The script uses EMBOSS (needle) to perform the alignments.

Type python3 pairwiseGlobalAlignment_sequenceIdentities.py -help to see the instructions.

"""

print(__doc__)

usage = "\n\tUSAGE: python3 pairwiseGlobalAlignment_sequenceIdentities.py [options]\n\n\
	*** Argumentos obrigatorios:\n\
		--sequences <arquivo fasta> (arquivo de sequencias para alinhamentos)\n\
		--cpus <integer> (numero de CPUs para executar alinhamentos)\n\n\
		--identity <float> (limite de identidade. Pares de sequencias com valor de identidade maior ou igual ao especificado serao consideradas redundantes.)\n\n\
	*** Argumento opcionais:\n\
		--ids <lista de IDs> (para realizar alinhamento somente das sequencias especificadas na lista de identificadores)\n\
		--dna2protein (quando possuir sequencia de nucleotideos e quiser alinha-las como proteína)\n\n\
	*** Modulos obrigatorios:\n\
		formataFasta_17Mar2020.py\n\
		extraiSeqFastaIDs_24Mar2020.py\n\n"


#realizar alinhamento global par a par de sequencias para obter identidade
#Vou usar o EMBOSS (ferramenta needle) para isso.

#export PATH=$PATH:/home/tuliomorgan/doutorado/programas/EMBOSS-6.6.0/emboss

def loopAlinhamento(cpus, fasta_selected, cutoff_identidade):
	import sys, re, subprocess

	## -- teste para saber o se needle esta configurado no sistema -- ##
	needle = subprocess.Popen(['which', 'needle'],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	case = needle.communicate() # communicate() returns a tuple (stdoutdata, stderrdata)
	if case[0]:
		print ('\nFerramenta Needle do pacote EMBOSS encontrada com sucesso.\n')
	else:
		print ("ERRO: Ferramenta Needle não encontrada. Obtenha o pacote EMBOSS em http://emboss.open-bio.org/html/adm/ch01s01.html ou configure a variavel de ambiente PATH para o diretorio que contem o executável needle.\n")
		exit()

	### -- iniciar subrotina de alinhamentos -- ##

	# lista com os valores de identidade de cada par, mas abaixo do vlaores de cutoff de identidade. Ja adicionei o cabecalho
	resultados_sequencias_nao_redundantes = ['query\tsubject\tidentidade(%)']
	# lista com os valores de identidade de cada par, mas maior ou igual do vlaores de cutoff de identidade. Ja adicionei o cabecalho
	resultados_sequencias_Redundantes = ['query\tsubject\tidentidade(%)']

	# criar um dicionario com todas as sequenicas do dataset. Vai ser util quando for recuperar as sequencias qeu devem ser mantidas no dataset, pois vou removendo-as da lista fasta_selected a medida que geramos os arquivos chunkID
	list_all_sequences = []
	for elemento in fasta_selected:
		list_all_sequences.append(elemento)

	# didionario de proteinas main. Como vou remove-las do dataset (para ela nao ser alinahda novamente com outras em outros loops), vamos criar um dicionario para salva-las e adiciona-las novamente na saida script.
	dict_proteinas_main = {}

	#dicionario com os tamanhos das sequencias. Vai ser util quando for calcular a %identidade
	dict_lengths = {}

	# enquanto houver 2 ou mais sequencias no arquivo fasta, alinhamentos devem ser conduzidos.
	while len(fasta_selected) >= 2:
		
		# definir proteina principal (vai ser alinhada contra todas as demais).
		proteinMain = fasta_selected[0] # primeiro elemento de 'fasta_selected'. Eh a proteina mais longa do dataset (uma vez que foram organizadas por tamanho com o 'orderSequencesBySize'.
		
		file = open('temp_proteinMain_alinhar.fasta', 'w')
		# salva o ID e a sequencia da proteina Main
		file.writelines(proteinMain)
		file.close()
		#print(fasta_selected)
		# remove a proteina principal do dicionario 'fasta_selected'. Nao devo alinhar ela com ela mesma.
		fasta_selected.remove(proteinMain)
		#print(fasta_selected)
		# adiciona-la no dicionario de proteinas Main
		proteinMain = proteinMain.split('\n')
		dict_proteinas_main[proteinMain[0]] = proteinMain[1]
		#print(dict_proteinas_main)

		# quebrar a lista 'fasta_selected' em numero de arquivos de acordo com o numero de processadores requisitados.
		# primeiro avalia se o numero de sequencias é maior ou igual ao numero de cpus requisitados. Caso contrario, adequa o numero de cpus ao numero de sequencias.
		if len(fasta_selected) <= int(cpus):
			cpus = len(fasta_selected)

		# elimina arquivos gerados do loop de alinhamento anterior.
		subprocess.run(['rm temp_*_chunkID*'], shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)

		# criacao dos arquivos de sequencias para alinhamento com a proteina main. Serao criados numero de aquivos igual ao numero de cpus requisitados, de modo que os alinhamentos sejam executados em paralelo.
		nSplit = 1
		while nSplit <= int(cpus):
			#print('iteracao:',str(nSplit))
			#print('cpus:',str(cpus))

			file = open('temp_sequences_chunkID'+str(nSplit), 'w')

			# se necessário, o primeiro arquivo (nSplit == 1) ira receber mais sequencias. Isso ocorre quando o numero de sequencias e/ou numero de cpus é impar.
			if nSplit == 1:
				numSeqEachFile = (len(fasta_selected)//int(cpus)) + (len(fasta_selected)%int(cpus))
			else:
				numSeqEachFile = len(fasta_selected)//int(cpus)

			# a variavel "numSeqEachFile" é o numero de sequencias em cada arquivo fasta
			while numSeqEachFile:
				file.writelines(fasta_selected[0]+'\n')
				seq_id = fasta_selected[0].split('\n')[0].lstrip('>')
				tamanho = len(fasta_selected[0].split('\n')[1])
				dict_lengths[seq_id] = tamanho
				del fasta_selected[0]
				numSeqEachFile -= 1
			file.close()
			nSplit += 1
			#print(dict_lengths)

		# chamada da funcao de alinhamento global par a par. A proteina main sera alinhada com cada proteina de cada arquivo fasta gerado no bloco anterior.
		chunkid = 1
		childProcess = []
		while chunkid <= cpus:

			nomeOutput = 'temp_align_sequences_chunkID'+str(chunkid)
			nomeSequencesAlinhar = 'temp_sequences_chunkID'+str(chunkid) # nomes dos aruqivos gerados no while anterior

			with open(nomeOutput, 'w') as exonerateOutput:

				needle_cline = subprocess.Popen(['needle', '-asequence', 'temp_proteinMain_alinhar.fasta', '-bsequence',nomeSequencesAlinhar, '-gapopen', str(10), '-gapextend', str(0.5), '-outfile',nomeOutput],stdout = subprocess.PIPE, stderr=subprocess.PIPE)
				childProcess.append(needle_cline)  # ... e adiciona o comando em uma lista de processos
				chunkid += 1

		for needleCall in childProcess: # para cada processo da lista 'childProcess' em execucao...
			needleCall.wait() # ... informa para aguardar sua finalizacao. Dessa forma, os processos serao executados paralelamente.

		subprocess.run(['cat temp_align_sequences_chunkID* > temp_all_alignments'], shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # concatena os resultados de alinhamentos.

		# traneferir os resultados de alinhamentos para uma variavel do programa.
		file = open('temp_all_alignments', 'r')
		alinhamentos = file.readlines()
		file.close()

		# analisar a identidade entre as sequencias.
		idsProteinasManter = []
		for linha in alinhamentos:
			# Nesse for, tambem recuperar o comprimento das sequencias alinhadas (tamanho_seq_1 & tamanho_seq_2). Calcular a %identidade considerando o tamanho da menor sequencia. Pode ser que ela seja completamente englobado pela outra (com 100% de identidade). Nesse caso, se pegar somente o tamanho da maior sequencia, a %identidade nao da 100% (mesmo a sequencia menor sendo totalmente englobada pela maior).
			if linha.startswith('# 1:'): # id da sequencia proteinMain.
				linhaSplit = linha.split()
				proteinaPrincipal = linhaSplit[2]
				tamanho_seq_1 = len(proteinMain[1]) # sequencia 'main'

			elif linha.startswith('# 2:'): # id da sequencia que alinhou com proteinMain.
				linhaSplit = linha.split()
				proteinaAlinhada = linhaSplit[2] # esse ID fica truncado pelo Needle (quebra no espaco em branco e pega o primeiro campo)
				
				# para recuperar o tamanho da sequencia 2, devemos busca-la do dicionario dict_lengths
				for chave in dict_lengths.keys():
					chaveSplit = chave.split(' ')
					if chave == proteinaAlinhada: # o ID proteinaAlinhada fica truncado pelo Needle (quebra no espaco em branco e pega o primeiro campo)
						tamanho_seq_2 = dict_lengths[chave]
						break
					elif chaveSplit[0] == proteinaAlinhada:
						tamanho_seq_2 = dict_lengths[chave]
						break

			# selecionar os IDs de proteinas nao redundantes. Essas irao para um novo loop de alinhamento.
			elif linha.startswith('# Identity:'):
				identidade = ''.join(re.findall(r'(?<=\()[ 0-9]+.[0-9]+(?<!%)', linha))
				numero_matches = int(''.join(re.findall(r'(?<=:)[ 0-9]+(?<!/)', linha)).lstrip(' ')) # recupera o numero de matches.

				# calcular a %identidade contra a menor sequencia. Pode ser que ela tenha sido englobada pela maior, com 100% de identidade.
				if tamanho_seq_1 > tamanho_seq_2:
					identidade_menor_seq = (numero_matches/tamanho_seq_2)*100
				elif tamanho_seq_1 < tamanho_seq_2:
					identidade_menor_seq = (numero_matches/tamanho_seq_1)*100
				#nao precisa testar tamanho_seq_1 == tamanho_seq_2 pois se isso ocorrer a identidade fornecida pela needle esta eh a certa.
			
				# se verdadeiro, significa que a sequencia menor foi englobado pela maior. Podemos remove-lo do dataset.
				if identidade_menor_seq == 100:
					#print('proteina descartar (englobada)',proteinaAlinhada)
					resultados_sequencias_Redundantes.append(proteinaPrincipal+'\t'+proteinaAlinhada+'\t'+str("{:.2f}".format(identidade_menor_seq)))

				elif float(identidade) < cutoff_identidade*100: # se verdadeiro, as sequencias NAO sao consideradas redundantes.
					#print('proteina manter',proteinaAlinhada)
					idsProteinasManter.append(proteinaAlinhada)
					# salvar o resultado de indentidade para as duas sequencias alinhadas
					resultados_sequencias_nao_redundantes.append(proteinaPrincipal+'\t'+proteinaAlinhada+'\t'+identidade)
				else:
					#print('proteina descartar',proteinaAlinhada)
					resultados_sequencias_Redundantes.append(proteinaPrincipal+'\t'+proteinaAlinhada+'\t'+identidade)
		
		# salva as sequencias nao redundantes em um novo dicionario
		new_fasta_selected = []
		for ids in idsProteinasManter:
			for elemento in list_all_sequences:
				ID = elemento.split('\n')[0].lstrip('>')
				if ids == ID:
					new_fasta_selected.append(elemento)
		#print('new_fasta_selected',new_fasta_selected)
		# substitui o novo fasta_selected (com as sequencias nao redundantes comparadas com a proteinaMain) com o nome fasta_select, e volta para o while do inicio do script para selecionar nova proteinMain.
		fasta_selected = new_fasta_selected

		if len(fasta_selected) == 1: #so sobrou uma sequencia no 'fasta_selected' e essa sequencia eh boa, pois seu ID foi mantido ate aqui (nao eh redundante)
			# adiciona-la no dicionario de proteinas Main
			dict_proteinas_main[fasta_selected.split('\n')[0]] = fasta_selected.split('\n')[1]

	return dict_proteinas_main, resultados_sequencias_nao_redundantes, resultados_sequencias_Redundantes
	
### -- uso do script de forma independente, nao como modulo de outro script -- ###

if __name__ == '__main__':
	import sys
	if len(sys.argv) < 4:
		print(usage)
		exit()

	else:
		traduz = False

		if 'dna2protein' in sys.argv:
			sys.argv.remove('dna2protein')
			traduz = True

		if '--sequences' in sys.argv:	
			argSequences = sys.argv.index('--sequences')+1
			file = open(sys.argv[argSequences])
			arq_fasta = file.readlines()
			file.close()
			
			import orderSequencesBySize #formata o arquivo multifasta e ordena por tamanho
			reverse = False # ordenar as sequencias da maior para menor.
			fasta_selected = orderSequencesBySize.orderBySize(arq_fasta, reverse)

		else:
			print(usage)
			exit()

		if '--cpus' in sys.argv: # numero de divisoes que serao aplicadas ao dataset. Isso pertimite executar alinhamentos em paralelo.
			argcpus = sys.argv.index('--cpus')+1
			cpus = int(sys.argv[argcpus])
		else:
			print(usage)
			exit()

		if '--identity' in sys.argv: # limite de identidade (inclusive) para considerar duas sequencias redundantes.
			argIdentity = sys.argv.index('--identity')+1
			cutoff_identidade = float(sys.argv[argIdentity])
		else:
			print(usage)
			exit()

		if '--ids' in sys.argv:
			argIDS = sys.argv.index('--ids')+1
			file = open(sys.argv[argIDS], 'r')
			lista_ids = file.readlines()
			file.close()
			
			import extraiSeqFastaIDs_24Mar2020
			fasta_selected = extraiSeqFastaIDs_24Mar2020.ExtraiSeqFastaIDs(lista_ids, fasta_selected)

		else:
			lista_ids = None

		if traduz:#se o argumento 'dna2protein' foi passado para o programa, fazemos a traducao da sequencia de DNA para depois realizar o alinhemnto global par a par.

			from Bio.Seq import Seq
			translated_fasta_selected = {}
			for chave in fasta_selected.keys():
				sequence = fasta_selected[chave]
				sequence = Seq(sequence) #transformar a sequencia j (uma string) em um objeto 'Seq' e salvar na variavel seq.
				protein_raw = str(sequence.translate()) #a variavel seq eh um objeto Seq e podemos aplicar a funcao translate. Apos, a traducao em aminoacidos, obtem-se resulado com a sequencia de aminoacidos e outras coisas (p.e. [Seq('TLYIKMVCLSLSEERGYIEALSVW*', HasStopCodon(ExtendedIUPACProtein(), '*'))]). Mas como quero imprimir somente a sequencia de aminoacidos, transformamos a variavel em string. Muito facil
				translated_fasta_selected[chave+' [translated]'] = protein_raw

			fasta_selected = translated_fasta_selected

		#chamada da funcao loopAlinhamento.
		dict_proteinas_main, resultados_sequencias_nao_redundantes, resultados_sequencias_Redundantes = loopAlinhamento(cpus, fasta_selected, cutoff_identidade)
		
		#### ----- para imprimir a chave e o valor (sequencia) devemos fazer um loop ----- #####
		file = open('fasta_selected_sem_redundanciaMP.fasta','w')
		for chave in dict_proteinas_main.keys():
			file.writelines(chave+'\n'+dict_proteinas_main[chave]+'\n')		
		file.close()

		### --- imprimir os valores de identidade para sequenicas consideradas NAO redundantes --- ###
		file = open('identidade_sequencias_nao_redundantes.txt','w')
		for elemento in resultados_sequencias_nao_redundantes:
			file.writelines(elemento+'\n')			
		file.close()

		### --- imprimir os valores de identidade para sequenicas consideradas Redundantes --- ###
		file = open('identidade_sequencias_redundantes.txt','w')
		for elemento in resultados_sequencias_Redundantes:
			file.writelines(elemento+'\n')			
		file.close()