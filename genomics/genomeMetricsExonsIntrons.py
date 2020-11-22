"""
============================
genomeMetricsExonsIntrons.py
============================

PT-BR:
Analisa a distribuicao de tamanho dos introns/exons e o numero de introns/exons de genes em um arquivo GFF3.

EN:
Analyzes the distribution of size of introns/exons and the number of introns/exons in a GFF3 file.

Type python3 genomeMetricsExonsIntrons.py -help to see the instructions

"""

print(__doc__)

usage = "\n\tgenomeMetricsExonsIntrons.py\n\n\
		Script para analisar a distribuicao de tamanho dos introns/exons e o numero de introns/exons de genes em um arquivo GFF3.\n\n\
		*** Argumento obrigatório\n\
		--input <arquivo GFF3>\n\
		--output <nome arquivo de saida>\n\n"

import re, sys, statistics 

if len(sys.argv) < 5:
	print(usage)
	exit()
else:
	if '--input' in sys.argv:
		argINPUT = sys.argv.index('--input')+1
		file = open(sys.argv[argINPUT],'r')
		fileINPUT = file.readlines()
		file.close()
	else:
		print(usage)
		exit()

	if '--output' in sys.argv:
		argOutput = sys.argv.index('--output')+1
		nomeOutput= str(sys.argv[argOutput])
	else:
		print(usage)
		exit()

	allIntronLengths = []
	allExonLengths = []
	#singleExonLengths = []
	geneExonCoordinates = []
	exonCount = 0
	exonCountGenome = 0
	intronCountGenome = 0
	numberOfGenes = 0

	for i,linha in enumerate(fileINPUT):
		if linha == '\n': # existem quebras de linha entre os genes no arquivo GFF3 de saida do EVM.
			continue # jumps to next iteration. Ignore the code below.
		linha = linha.rstrip()

		if '\tgene\t' in linha:
			numberOfGenes += 1

		if i+1==len(fileINPUT) or '\tgene\t' in linha and i:
			
			geneExonCoordinates = '\t'.join(geneExonCoordinates)
			geneExonCoordinates = geneExonCoordinates.split('\t')
	
			toCountExons = geneExonCoordinates.copy() # criei essa variavel para obter o tamnho dos exons. A variavel original 'geneExonCoordinates' vai sendo deletada enquanto calculo o tamanho dos introns no loop while abaixo.
			#print(toCountExons)
			#print(geneExonCoordinates)
			# obter tamnaho de introns de cada gene.
			itera = 0
			if exonCount > 1: # somente genes com mais de 1 exon possuem introns...
				while itera < exonCount-1: # com 2 exons faco apenas 1 iteracao, com 3 exons sao 2 iteracoes e assim por diante.
					allIntronLengths.append(int(geneExonCoordinates[2])-int(geneExonCoordinates[1]))
					geneExonCoordinates.remove(geneExonCoordinates[0])
					geneExonCoordinates.remove(geneExonCoordinates[0])
					#print(geneExonCoordinates)
					itera += 1
			
			#print(toCountExons)
			# obter tamanho de exons de cada gene (mais simples!)
			itera = 0
			#print(exonCount)
			if exonCount > 1:
				while itera < exonCount-1: # com 2 exons faco apenas 1 iteracao, com 3 exons sao 2 iteracoes e assim por diante.
					
					allExonLengths.append(int(toCountExons[1])-int(toCountExons[0]))
					toCountExons.remove(toCountExons[0])
					toCountExons.remove(toCountExons[0])
					#print(toCountExons)
					itera += 1
			
			elif exonCount == 1:
				allExonLengths.append(int(toCountExons[1])-int(toCountExons[0]))
			
			intronCountGenome += (exonCount-1)
			geneExonCoordinates = []
			exonCount = 0

		elif '\texon\t' in linha:
			linhaSplit = linha.split('\t')
			fita = linhaSplit[6]
			
			if fita == '+': 
				geneExonCoordinates.append(linhaSplit[3]+'\t'+linhaSplit[4])
			else:
				geneExonCoordinates.insert(0,linhaSplit[3]+'\t'+linhaSplit[4]) # se o gene estiver an fita reversa, altero a posicao que os os exons aparecem (originalmente estao escritos na ordem inversa, isto eh, o exon na porcao 3' do gene esta escrito primeiro)

			exonCount += 1 # vai sendo restaurado a zero em cada gene.
			exonCountGenome += 1

	meanIntronSize = (sum(allIntronLengths))/int(intronCountGenome)
	meanIntronPerGene = int(intronCountGenome)/int(numberOfGenes)

	meanExonSize = (sum(allExonLengths))/int(exonCountGenome)
	meanExonPerGene = int(exonCountGenome)/int(numberOfGenes)
	
	"""
	print("Numero de genes: "+str(numberOfGenes))
	print("Numero total de exons: "+str(exonCountGenome))
	#print("Tamanhos de exons: "+str(allExonLengths))
	"""
	#print("Numero total de introns: "+str(intronCountGenome))
	#print("Media de tamanho de introns: "+str(meanIntronSize))
	#print("Media de introns por gene: "+str(meanIntronPerGene))
	#print("Tamanhos de introns: "+str(allIntronLengths))
	#print("Desvio padrao amostral tamanho de introns: ", statistics.stdev(allIntronLengths))

	file = open(nomeOutput,'w')
	file.writelines('# Genes (protein-coding)\tTotal exons\tMean exons/gene\tExon mean length\tTotal introns\tMean introns/gene\tIntron mean length\tStandard Deviation (Intron mean length)\n')
	file.writelines(str(numberOfGenes)+'\t'+str(exonCountGenome)+'\t'+str(meanExonPerGene)+'\t'+str(meanExonSize)+'\t'+str(intronCountGenome)+'\t'+str(meanIntronPerGene)+'\t'+str(meanIntronSize)+'\t'+str(statistics.stdev(allIntronLengths)))
	file.close()
