"""
=========================
runExonerate.v2.py
=========================

PT-BR:
Script para executar para executar jobs do Exonerate em paralelo. Utilizo principalmente para predição de genes em novas montagens de genomas.

EN:
Script to run to run Exonerate in parallel. I use it mainly to predict genes in new genome assemblies.

Type python3 runExonerate.v2.py -help to see the instructions.

"""
usage = ('\n	runExonerate.v2.py\n\n\
	Script para executar jobs do Exonerate em paralelo.\n\
	A atualizacao desse script foi feita em 23 de Marco de 2020 para incluir a opcao `--refine`\n\n\
# -------------------------------------------------------------------------------- #\n\n\
 *** Argumentos obrigatorios:\n\
		--genome <arquivo multifasta de genoma>\n\
		--sequences <arquivo multifasta de sequencias de nucleotideos ou aminoacidos>\n\
		--chunks <numero de particoes do database (geralmente o arquivo de genoma)>\n\
			`quanto mais particoes, mais rapido tende a ser a execucao do programa`\n\
		--querytype <dna | protein>\n\
		--targettype <dna | protein>\n\
		--model <modelo para alinhamentos>\n\
			`para alinhamento de transcritos contra genoma utilizo est2genome.`\n\
			`para alinhamento de proteinas contra genoma utilizo protein2genome.`\n\
		--bestn <int>\n\
 *** Argumentos opcionais:\n\
		--fsmmemory <int>. Default: 10000\n\
		--showtargetgff <yes|no>. Default: no\n\
		--exhaustive <yes|no>. Default: no\n\
		--dpmemory <int>. Default: 15000\n\n')

import sys
if len(sys.argv) < 9:
	print (usage)
	exit()

else:
	import subprocess

	if '--genome' in sys.argv:
		argGenome = sys.argv.index('--genome')+1
		genome = str(sys.argv[argGenome]) # nome do arquivo de genoma para passar ao exonerate. Geralmente eh o target
	else:
		print(usage)
		exit()

	if '--sequences' in sys.argv:	
		argSequences = sys.argv.index('--sequences')+1
		sequences = str(sys.argv[argSequences]) # nome do arquivo de sequences de genes para passar ao exonerate. Geralmente eh a query
	else:
		print(usage)
		exit()

	if '--chunks' in sys.argv: # numero de divisoes que serao aplicadas ao TARGET. A query nao sera dividida.
		argChunks = sys.argv.index('--chunks')+1
		chunks = int(sys.argv[argChunks])
	else:
		print(usage)
		exit()

	if '--querytype' in sys.argv:
		argQueryType = sys.argv.index('--querytype')+1
		queryType = str(sys.argv[argQueryType])
		if queryType != 'protein':
			if queryType != 'dna':
				print('Argumento',queryType ,'nao eh valido para --querytype. Utilize protein ou dna.\n')
				print (usage)
				exit()
	else:
		print(usage)
		exit()

	if '--targettype' in sys.argv:
		argTargetType = sys.argv.index('--targettype')+1
		targetType = str(sys.argv[argTargetType])
		if targetType != 'protein':
			if targetType != 'dna':
				print('Argumento',targetType ,'nao eh valido para --targettype. Utilize protein ou dna.\n')
				print (usage)
				exit()
	else:
		print(usage)
		exit()
	
	if '--model' in sys.argv:
		argModel = sys.argv.index('--model')+1
		model = str(sys.argv[argModel])
	else:
		print(usage)
		exit()

	if '--bestn' in sys.argv:
		argBestn = sys.argv.index('--bestn')+1
		bestn = str(sys.argv[argBestn])
	else:
		print(usage)
		exit()

	if '--percent' in sys.argv:
		argPercent = sys.argv.index('--percent')+1
		percent = str(sys.argv[argPercent])
	else:
		print(usage)
		exit()

	if '--score' in sys.argv:
		argScore = sys.argv.index('--score')+1
		score = str(sys.argv[argScore])
	else:
		print(usage)
		exit()

	if '--output' in sys.argv:
		argOutput = sys.argv.index('--output')+1
		outputDir = str(sys.argv[argOutput])+'/'
	else:
		testaDir = subprocess.Popen(['cd output_exonerate/'], shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
		dirExistente = testaDir.communicate()[1]
		if dirExistente:
			print('Existe um diretorio chamado output_exonerate, provavelmente com resultados de uma execucao do Exonerate. Remova esse diretorio, altere seu nome ou utilize o argumento --output na chamada do script para indicar um nome diferente.\n')
			exit()
		else:
			outputDir = 'output_exonerate/'

	if '--fsmmemory' in sys.argv:
		argFsmmemory = sys.argv.index('--fsmmemory')+1
		fsmmemory = str(sys.argv[argFsmmemory])
	else:
		fsmmemory = '10000'

	if '--showtargetgff' in sys.argv:
		argShowtargetgff = sys.argv.index('--showtargetgff')+1
		showtargetgff = str(sys.argv[argShowtargetgff])
	else:
		showtargetgff = 'no'
	
	if '--showalignment' in sys.argv:
		argShowalignment = sys.argv.index('--showalignment')+1
		showalignment = str(sys.argv[argShowalignment])
	else:
		showalignment = 'no'

	if '--ryo' in sys.argv:
		argRyo = sys.argv.index('--ryo')+1
		ryo = str(sys.argv[argRyo]) # geralmente uso: --ryo ">%ti\t%td\t%tl\n%tas\n>%qi\t%qd\t%ql\n%qas"  Para imprimir no output a sequencia fasta do target que alinhou e tambem a sequencia fasta da query que alinhou.
	else:
		ryo = 'none'

	if '--exhaustive' in sys.argv:
		argExhaustive = sys.argv.index('--exhaustive')+1
		exhaustive = str(sys.argv[argExhaustive]) # usar o alinahamento exaustivo, porem demora muito mais que o heuristico (default)

	if '--dpmemory' in sys.argv: # dynamic programming memory usage
		argDPmemory = sys.argv.index('--dpmemory')+1
		dpmemory = str(sys.argv[argDPmemory]) # util quando usar o alinahamento exaustivo. Aumenta a velocidade de execucao.
	else:
		dpmemory = '15000'

	if '--refine' in sys.argv: # Force exonerate to refine alignments generated by heuristics using dynamic programming over larger regions. This takes more time, but improves the quality of the final alignments.
		argRefine = sys.argv.index('--refine')+1
		refine = str(sys.argv[argRefine]) # util quando usar o alinahamento exaustivo. Aumenta a velocidade de execucao.
	else:
		refine = 'none'	
	
	nomeGenoma = str((genome.split('/'))[-1]) # pegar somente o nome do arquivo de genoma, sem o caminho relativo, para adicionar no nome do arquivo de saida.
	nomeSequences = str((sequences.split('/'))[-1])

	subprocess.Popen(['mkdir', outputDir])
	
	chunkid = 1
	childProcess = []
	while chunkid <= chunks:

		nomeOutput = nomeGenoma+'_vs_'+nomeSequences+'_output_chunkID'+str(chunkid) # nome do arquivo de saida ira conter de qual particao do target foi oriundo.
		
		with open(nomeOutput, 'w') as exonerateOutput:

			if '--exhaustive' in sys.argv: # a chamada do Exonerate muda dependendo da configuracao de alinhamento. Alinhamento exsaustivo é mais acurado, porem mais demorado para executar.
				exonerate = subprocess.Popen(['exonerate', '--query', sequences, '--target', genome, '--querytype', queryType, '--targettype', targetType, '--targetchunkid', str(chunkid), '--targetchunktotal', str(chunks), '--model', model, '--bestn', bestn, '--percent', percent, '--score', score, '--fsmmemory', fsmmemory, '--showalignment', showalignment, '--showtargetgff', showtargetgff, '--ryo', ryo, '--exhaustive', exhaustive, '--dpmemory', dpmemory, '--refine', refine], stdout = exonerateOutput, stderr = subprocess.PIPE) # inicia um processo...
				childProcess.append(exonerate)  # ... e adiciona o comando em uma lista de processos


			else: # aqui vai usar alinhamento heuristico (mais rapido, provavelmente menos acurado)
				exonerate = subprocess.Popen(['exonerate', '--query', sequences, '--target', genome, '--querytype', queryType, '--targettype', targetType, '--targetchunkid', str(chunkid), '--targetchunktotal', str(chunks), '--model', model, '--bestn', bestn, '--percent', percent, '--score', score, '--fsmmemory', fsmmemory, '--showalignment', showalignment, '--showtargetgff', showtargetgff, '--ryo', ryo, '--dpmemory', dpmemory, '--refine', refine], stdout = exonerateOutput, stderr = subprocess.PIPE) # inicia um processo...
				childProcess.append(exonerate)  # ... e adiciona o comando em uma lista de processos

		chunkid += 1

	for exonerateCall in childProcess: # para cada processo da lista 'childProcess' em execucao...
		exonerateCall.wait() # ... informa para aguardar sua finalizacao. Dessa forma, os processos serao executados paralelamente.

	subprocess.run(['cat *_output_chunkID* > temp'], shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # concatena os resultados do exonerate.
	output = nomeGenoma+'_vs_'+nomeSequences+'_output_complete.txt' # nome do arquivo de resultados concatenados.
	subprocess.Popen(['mv', 'temp', output], stdout = subprocess.PIPE, stderr = subprocess.PIPE) # renomeia o arquivo de resultados concatenados.

	nMover = 1
	listaNomesArquivos = []
	while nMover <= chunks: # loop para mover arquivos de saida do exonerate para o diretorio 'outputDir' criado no inicio do script.
		nomeArquivoMover = nomeGenoma+'_vs_'+nomeSequences+'_output_chunkID'+str(nMover) # vai acessar cada nome de arquivo, como se fosse a variavel 'nomeOutput'
		subprocess.Popen(['mv', nomeArquivoMover, outputDir], stdout = subprocess.PIPE, stderr = subprocess.PIPE) # move cada arquivo para o diretorio de output
		nMover += 1

	subprocess.Popen(['mv', output, outputDir+output], stdout = subprocess.PIPE, stderr = subprocess.PIPE) # move o arquivo de resultados concatenados para o diretorio de output.

	import exonerate2gff # converter arquivo de saida do exonerate em arquivo GFF (estilo NCBI).
	output = nomeGenoma+'_vs_'+nomeSequences+'_output_complete.txt' # nome do arquivo de resultados concatenados.
	
	file = open(outputDir+output, 'r')
	outputExonerateConcantenado = file.readlines()
	file.close()

	exonerate2gff.exonerate2GFF(outputExonerateConcantenado, output+'.gff3')

	subprocess.Popen(['mv', output+'.gff3', outputDir+output+'.gff3'])