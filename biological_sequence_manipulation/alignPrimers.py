"""
===============
alignPrimers.py
===============

PT-BR:
Script para realizar alinhamentos entre primers (foward, reverse, probe) e sequencias alvo. Esse processo necessita do Exonerate v.2.2.0 instalado e com arquivo executavel acessivel pelo sistema. As sequencias de entrada devem estar no formato fasta.
Os arquivos de saida consistem de uma tabela, onde o header eh composto pela sequencia do primer e a primeira coluna sao os IDs das sequencias alvo. A tabela eh preenchida com zero ou um para cada posicao do alinhamento. Zero significa match e 1 significa mismatch. Em caso de mismatch, o nucleotideo alterado na sequencia fica em parenteses. Alem disso, fornece arquivos fasta das regioes amplificadas.

EN:
Script to perform alignments between primers (foward, reverse, probe) and target sequences. This process requires Exonerate v.2.2.0 installed and with the executable file accessible by the system. The input strings must be in the fasta format.
The output files consist of a table, where the first row is the primer sequence and the first column are the target sequence IDs. The table is filled with 0`s or 1`s for each position in the alignment. Zero means match and 1 means mismatch. In the case of mismatch, the altered nucleotide is depicted in parentheses. In addition, it provides fasta files from the amplified regions.

Type python3 alignPrimers.py -help to see the usage information.

"""

print(__doc__)

usage = "\n\tUSAGE: python3 alignPrimers.py [options]\n\n\
	Script para alinhar primers em sequenicas alvo. Se houver nucleotideos diferentes de ACTG ou gaps, nao ira salvar o alinhamento no arquivo de saida.\n\n\
	*** Argumentos obrigatorios:\n\
		--sequences <arquivo fasta> (arquivo de sequencias para alinhamentos)\n\
		--primers <arquivo fasta> (a primeira sequencia é considerada o primer foward, a segunda o primer reverse, e a terceira (se houver) é considerada a probe. As sequencias dos primers e da probe (essa é opcional) devem estar no memo arquivo fasta e na direção 5'3')\n"

import sys, re, subprocess

def loopAlinhamento(fasta_selected):
	
	### -- iniciar subrotina de alinhamentos -- ##
	
	# gerar variaveis de saida.
	alinhamentos_primer_foward = [] # reporta matches e mismatches das sequencias com o primer foward
	alinhamentos_primer_reverse = [] # reporta matches e mismatches das sequencias com o primer reverse	
	alinhamentos_probe = [] # reporta matches e mismatches das sequencias com o probe
	sequencias_amplificacao = [] # sequencias fasta da regiao amplificada pelos primers.

	# percorrer cada sequencia para alinhar com o primers/probes passados ao script.
	for chave in fasta_selected.keys():
		# remover e substituir caracteres que atraplham salvar o arquivo.
		seqID = chave.split('\t')[0] # tem uma tabulacao no ID das sequencias. VOu pegar so o primeiro campo.
		seqID = seqID.lstrip('>') # remove o sinal de maior '>' no inicio do ID.
		seqID = re.sub(r'/|\|','_',seqID)
		seqID = re.sub(r'\"','',seqID)
		
		# salvar a sequencia alvo em um arquivo fasta
		targetSeq = seqID+'.fasta'
		file = open(targetSeq, 'w')
		file.writelines(chave+'\n'+fasta_selected[chave])
		file.close()

		# tambem salvar a sequencia complementar reversa apra alinhar no primer reverse. Se nao fizer isso, o e exonerate faz a revcomp do primer, o que nao esta certo.
		revComp = ''
		for nuc in fasta_selected[chave][::-1]:
			if nuc.lower() == 'a':
				revComp += 't'
			elif nuc.lower() == 't':
				revComp += 'a'
			elif nuc.lower() == 'c':
				revComp += 'g'
			elif nuc.lower() == 'g':
				revComp += 'c'
			elif nuc.lower() == 'r':
				revComp += 'y'
			elif nuc.lower() == 'y':
				revComp += 'r'
			elif nuc.lower() == 'k':
				revComp += 'm'
			elif nuc.lower() == 'm':
				revComp += 'k'
			elif nuc.lower() == 'b':
				revComp += 'v'
			elif nuc.lower() == 'v':
				revComp += 'b'
			elif nuc.lower() == 'd':
				revComp += 'h'
			elif nuc.lower() == 'h':
				revComp += 'd'
			else:
				revComp += nuc
		revComp = revComp.upper()
		targetSeq_revComp = seqID+'_revComp.fasta'
		file = open(targetSeq_revComp, 'w')
		file.writelines(chave+'\n'+revComp)
		file.close()

		# --- alinhar PRIMER FOWARD --- #
		nomeOutput = seqID+'_foward_primer_aligned'
		with open(nomeOutput, 'w') as exonerateOutput:
			exonerate = subprocess.Popen(['exonerate', '--query', 'primer_foward.fasta', '--target', targetSeq, '--querytype', 'dna', '--targettype', 'dna', '--model', 'affine:bestfit', '--bestn', '1', '--percent', '50', '--score', '30','--dnawordlen', '3', '--showalignment', 'yes', '--exhaustive', 'yes', '--subopt', 'yes', '--ryo', '#Resultado\t%qi\t%ti\tidentidade:%ei\tmismatches:%em\n'], stdout = exonerateOutput, stderr = subprocess.PIPE)

		stderr = exonerate.communicate()
		#if stderr:
			#print('ERRO reportado pelo exonerate do alinhamento do primer_foward.fasta com '+chave+'\n'+str(stderr))

		# --- alinhar PRIMER REVERSE --- #
		nomeOutput = seqID+'_reverse_primer_aligned'
		with open(nomeOutput, 'w') as exonerateOutput:
			exonerate = subprocess.Popen(['exonerate', '--query', 'primer_reverse.fasta', '--target', targetSeq_revComp, '--querytype', 'dna', '--targettype', 'dna', '--model', 'affine:bestfit', '--bestn', '1', '--percent', '50', '--score', '30','--dnawordlen', '3', '--showalignment', 'yes', '--exhaustive', 'yes', '--subopt', 'yes', '--ryo', '#Resultado\t%qi\t%ti\tidentidade:%ei\tmismatches:%em\n'], stdout = exonerateOutput, stderr = subprocess.PIPE)

		stderr = exonerate.communicate()
		#if stderr:
			#print('ERRO reportado pelo exonerate do alinhamento do primer_reverse.fasta com '+chave+'\n'+str(stderr))

		# --- alinhar PRIMER PROBE, se houver --- #
		if probe:
			nomeOutput = seqID+'_probe_primer_aligned'
			with open(nomeOutput, 'w') as exonerateOutput:
				exonerate = subprocess.Popen(['exonerate', '--query', 'probe.fasta', '--target', targetSeq, '--querytype', 'dna', '--targettype', 'dna', '--model', 'affine:bestfit', '--bestn', '1', '--percent', '50', '--score', '30','--dnawordlen', '3', '--showalignment', 'yes', '--exhaustive', 'yes', '--subopt', 'yes', '--ryo', '#Resultado\t%qi\t%ti\tidentidade:%ei\tmismatches:%em\n'], stdout = exonerateOutput, stderr = subprocess.PIPE)

			stderr = exonerate.communicate()
			#if stderr:
				#print('ERRO reportado pelo exonerate do alinhamento do probe.fasta com '+chave+'\n'+str(stderr))

		# remover as sequencias fasta, pois ja efetuamos os alinhamentos.
		subprocess.run(['rm',targetSeq])
		subprocess.run(['rm', targetSeq_revComp])

		#### ---- recuperar a linha com o resultado no formato 'vulgar', que tem as coordenadas do alinhamento na target ---- ####
		
		############### ---------------- alinhamento com primer foward ------------------ ###################
		vulgar_raw = subprocess.Popen(['grep','vulgar', seqID+'_foward_primer_aligned'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = vulgar_raw.communicate()
		if stderr:
			print('ERRO ao executar o comando `grep` para recuperar a linha `vulgar:` de saida do exonerate.\n')
		elif not stdout: # nao gerou alinhamento.
			print('AVISO: a sequencia '+chave.lstrip('>')+' NAO gerou alinhamento com o primer foward '+primers[0])
			analisar_sequencia_foward = False
		else:
			coordenadas = (str(stdout).split(' '))[6:8] # recupera as coordendas inicial e final do bloco 'vulgar'.
			regiao_alinhada_foward = (fasta_selected[chave][int(coordenadas[0]):int(coordenadas[1])]).upper()
			regiao_alinhada_foward = list(regiao_alinhada_foward)

			# salvar a coordenada inicial de alinhamneto do primer foward para recuperar a regiao de amplificacao. Vai ser util para construir as figuras de alinhaentos multiplos.
			coordenada_inicial_regiao_amplificacao = int((str(stdout).split(' '))[6])

			# verificar se os tamanhos das regiao alinhada e do primer coincidem antes de chamar a funcao para identificar as posicoes de match e mismatch (count_matches_mismatches)
			if len(regiao_alinhada_foward) != len(seq_foward):
				print('AVISO: Os comprimentos da regiao alinhada no gene e o tamanho do primer_foward.fasta NAO coincidem.\nPossivelmente a região alinhada possui gaps na sequencia.\n Recomenda-se verificar a sequencia',chave,'e o primer foward.\n')
				# buscar onde esta o gap na sequencia.
				vulgar_block = str(stdout).split(' ')
				if 'G' in vulgar_block: # se verdadeiro, existiu um gap no alinhamento, por isso as sequencias estao com tamanho diferente.
					# posicao do gap eh dada pelo elemento anteior ao 'G' no vulgar block. Mas nao precisamos subtrair 1, pois o insert o elemento desejado na posicao anterior ao especificado.
					posicao_gap = vulgar_block.index('G')
					regiao_alinhada_foward.insert(posicao_gap, '-')

			# verificar existencia de nucleotideos diferentes de ACTG ou gaps. Se encontrar nao salva o alinhamento.
			analisar_sequencia_foward = True
			for nucl in regiao_alinhada_foward:
				if nucl == 'A' or nucl == 'C' or nucl == 'G' or nucl == 'T' or nucl == '-':
					continue
				else: # se entra nesse else significa que existem nucleotideos ambiguos ou nao determinados na sequencia.
					analisar_sequencia_foward = False
					print("AVISO: foi detectado nucleotideos ambiguos (Y,W,S,M,...) ou indeterminados (N) na sequencia "+chave+" ao alinhar com o primer foward.\nEla nao sera analisada pela funcao 'count_matches_mismatches', logo, nao sera incluida na tabela final de matches e mismatches.\n")
					break
		
		############### ---------------- alinhamento com primer reverse ------------------ ###################
		vulgar_raw = subprocess.Popen(['grep','vulgar', seqID+'_reverse_primer_aligned'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = vulgar_raw.communicate()
		if stderr:
			print('ERRO ao executar o comando `grep` para recuperar a linha `vulgar:` de saida do exonerate.\n')
		elif not stdout: # nao gerou alinhamento.
			print('AVISO: a sequencia '+chave.lstrip('>')+' NAO gerou alinhamento com o primer reverse '+primers[2])
			analisar_sequencia_reverse = False
		else:
			coordenadas = (str(stdout).split(' '))[6:8] # recupera as coordendas inicial e final do bloco 'vulgar'.
			# pegar a complementar reversa (variavel "revComp"). Foi criada antes dos alinhamentos para gerar o arquivo fasta
			regiao_alinhada_reverse = (revComp[int(coordenadas[0]):int(coordenadas[1])]).upper()
			regiao_alinhada_reverse = list(regiao_alinhada_reverse)

			# salvar a coordenada INICIAL de alinhamneto do primer reverse para recuperar a regiao de amplificacao. Nao eh a coordenada final pois a sequencia alinhada foi a complementar reversa. Basta subtrair essa coordenada inicial do comprimento da sequencia apra encontrar a coordenada final com relação a fita senso.
			coordenada_final_regiao_amplificacao = len(revComp) - int((str(stdout).split(' '))[6])

			# verificar se os tamanhos das regiao alinhada e do primer coincidem antes de chamar a funcao para identificar as posicoes de match e mismatch (count_matches_mismatches)
			if len(regiao_alinhada_reverse) != len(seq_reverse):
				print('AVISO: Os comprimentos da regiao alinhada no gene e o tamanho do primer_reverse.fasta NAO coincidem.\nRevisar a sequencia',chave,'e o primer reverse. A região alinhada possui gaps na sequencia.\n')
				# buscar onde esta o gap na sequencia.
				vulgar_block = str(stdout).split(' ')
				if 'G' in vulgar_block: # se verdadeiro, existiu um gap no alinhamento, por isso as sequencias estao com tamanho diferente.
					# posicao do gap eh dada pelo elemento anteior ao 'G' no vulgar block. Mas nao precisamos subtrair 1, pois o insert o elemento desejado na posicao anterior ao especificado.
					posicao_gap = vulgar_block.index('G')
					regiao_alinhada_reverse.insert(posicao_gap, '-')

			# verificar existencia de nucleotideos diferentes de ACTG ou gaps. Se encontrar nao salva o alinhamento.
			analisar_sequencia_reverse = True
			for nucl in regiao_alinhada_reverse:
				if nucl == 'A' or nucl == 'C' or nucl == 'G' or nucl == 'T' or nucl == '-':
					continue
				else: # se entra nesse else significa que existem nucleotideos ambiguos ou nao determinados na sequencia.
					analisar_sequencia_reverse = False
					print("AVISO: foi detectado nucleotideos ambiguos (Y,W,S,M,...) ou indeterminados (N) na sequencia "+chave+" ao alinhar com o primer reverse.\nEla nao sera analisada pela funcao 'count_matches_mismatches', logo, nao sera incluida na tabela final de matches e mismatches.\n")
					break


		############### ---------------- alinhamento com primer probe ------------------ ###################
		if probe: # primeiro avaliar se existe a probe para o set de primers.
			vulgar_raw = subprocess.Popen(['grep','vulgar', seqID+'_probe_primer_aligned'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			stdout, stderr = vulgar_raw.communicate()
			if stderr:
				print('ERRO ao executar o comando `grep` para recuperar a linha `vulgar:` de saida do exonerate.\n')
			elif not stdout: # nao gerou alinhamento.
				print('AVISO: a sequencia '+chave.lstrip('>')+' NAO gerou alinhamento com a probe '+primers[4])
				analisar_sequencia_probe = False
			else:
				coordenadas = (str(stdout).split(' '))[6:8] # recupera as coordendas inicial e final do bloco 'vulgar'.
				regiao_alinhada_probe = (fasta_selected[chave][int(coordenadas[0]):int(coordenadas[1])]).upper()
				regiao_alinhada_probe = list(regiao_alinhada_probe)

				# verificar se os tamanhos das regiao alinhada e do primer coincidem antes de chamar a funcao para identificar as posicoes de match e mismatch (count_matches_mismatches)
				if len(regiao_alinhada_probe) != len(seq_probe):
					print('AVISO: Os comprimentos da regiao alinhada no gene e o tamanho do primer_probe.fasta NAO coincidem.\nRevisar a sequencia',chave,'e o primer probe. A região alinhada possui gaps na sequencia.\n')
					# buscar onde esta o gap na sequencia.
					vulgar_block = str(stdout).split(' ')
					if 'G' in vulgar_block: # se verdadeiro, existiu um gap no alinhamento, por isso as sequencias estao com tamanho diferente.
						# posicao do gap eh dada pelo elemento anteior ao 'G' no vulgar block. Mas nao precisamos subtrair 1, pois o insert o elemento desejado na posicao anterior ao especificado.
						posicao_gap = vulgar_block.index('G')
						regiao_alinhada_probe.insert(posicao_gap, '-')
					
				# verificar existencia de nucleotideos diferentes de ACTG ou gaps. Se encontrar nao salva o alinhamento.
				analisar_sequencia_probe = True
				for nucl in regiao_alinhada_probe:
					if nucl == 'A' or nucl == 'C' or nucl == 'G' or nucl == 'T' or nucl == '-':
						continue
					else: # se entra nesse else significa que existem nucleotideos ambiguos ou nao determinados na sequencia.
						analisar_sequencia_probe = False
						print("AVISO: foi detectado nucleotideos ambiguos (Y,W,S,M,...) ou indeterminados (N) na sequencia "+chave+" ao alinhar com o primer probe.\nEla nao sera analisada pela funcao 'count_matches_mismatches', logo, nao sera incluida na tabela final de matches e mismatches.\n")
						break

		################-----------------------------------################################------------------############################
		# chama a funcao para contar e indicar as posicoes de match e mismatch entre seqeuncia e primer APENAS se nao houver nucleotideos ambiguos e indeterminados nas sequenicas alinhadas.
		if analisar_sequencia_foward and analisar_sequencia_reverse and analisar_sequencia_probe:
			alinhamentos_primer_foward.append(count_matches_mismatches(chave, regiao_alinhada_foward, seq_foward))
			alinhamentos_primer_reverse.append(count_matches_mismatches(chave, regiao_alinhada_reverse, seq_reverse))
			alinhamentos_probe.append(count_matches_mismatches(chave, regiao_alinhada_probe, seq_probe))

			# salvar a sequencia fasta da regiao amplificada.
			sequencias_amplificacao.append(chave+'\n'+fasta_selected[chave][coordenada_inicial_regiao_amplificacao:coordenada_final_regiao_amplificacao]) # adicionei -5 e +5 as coordenadas para recuperar com margem de sequencia onde os primers anelam.

	return alinhamentos_primer_foward, alinhamentos_primer_reverse, alinhamentos_probe, sequencias_amplificacao

def count_matches_mismatches(regiao_alinhada_ID, regiao_alinhada_sequencia, primer_sequencia):
	saida_temp = [] # gerar uma lista temporaria para salvar dados. depois tranferir para a lista de saida 'alinhamentos_primers'
	regiao_alinhada_ID = regiao_alinhada_ID.lstrip('>') # remove o sinal de maior para adicionar o ID da sequencia.
	regiao_alinhada_ID = re.sub(r'\t+','_',regiao_alinhada_ID) # substitui eventuais tabulacao no ID por underline, pois o separador de campo no arquivo final sera tabulação.
	saida_temp.append(regiao_alinhada_ID)
	for i in range(0,len(primer_sequencia),1):
		
		if regiao_alinhada_sequencia[i] == primer_sequencia[i]:
			saida_temp.append('0') # zero significa match

		elif primer_sequencia[i] == 'R':
			if regiao_alinhada_sequencia[i] == 'A' or regiao_alinhada_sequencia[i] == 'G':
				saida_temp.append('0') # zero significa match
			else:
				saida_temp.append('1('+regiao_alinhada_sequencia[i]+')') # um significa mismatch. ALem disso, reporta qual o nucleotideo na sequencia dessa posicao que nao deu match

		elif primer_sequencia[i] == 'Y':
			if regiao_alinhada_sequencia[i] == 'C' or regiao_alinhada_sequencia[i] == 'T':
				saida_temp.append('0') # zero significa match
			else:
				saida_temp.append('1('+regiao_alinhada_sequencia[i]+')') # um significa mismatch. ALem disso, reporta qual o nucleotideo na sequencia dessa posicao que nao deu match


		elif primer_sequencia[i] == 'S':
			if regiao_alinhada_sequencia[i] == 'G' or regiao_alinhada_sequencia[i] == 'C':
				saida_temp.append('0') # zero significa match
			else:
				saida_temp.append('1('+regiao_alinhada_sequencia[i]+')') # um significa mismatch. ALem disso, reporta qual o nucleotideo na sequencia dessa posicao que nao deu match


		elif primer_sequencia[i] == 'W':
			if 	regiao_alinhada_sequencia[i] == 'A' or regiao_alinhada_sequencia[i] == 'T':
				saida_temp.append('0') # zero significa match
			else:
				saida_temp.append('1('+regiao_alinhada_sequencia[i]+')') # um significa mismatch. ALem disso, reporta qual o nucleotideo na sequencia dessa posicao que nao deu match


		elif primer_sequencia[i] == 'K':
			if regiao_alinhada_sequencia[i] == 'G' or regiao_alinhada_sequencia[i] == 'T':
				saida_temp.append('0') # zero significa match
			else:
				saida_temp.append('1('+regiao_alinhada_sequencia[i]+')') # um significa mismatch. ALem disso, reporta qual o nucleotideo na sequencia dessa posicao que nao deu match


		elif primer_sequencia[i] == 'M':
			if regiao_alinhada_sequencia[i] == 'A' or regiao_alinhada_sequencia[i] == 'C':
				saida_temp.append('0') # zero significa match
			else:
				saida_temp.append('1('+regiao_alinhada_sequencia[i]+')') # um significa mismatch. ALem disso, reporta qual o nucleotideo na sequencia dessa posicao que nao deu match


		elif primer_sequencia[i] == 'B':
			if regiao_alinhada_sequencia[i] == 'G' or regiao_alinhada_sequencia[i] == 'C' or regiao_alinhada_sequencia[i] == 'T': 
				saida_temp.append('0') # zero significa match
			else:
				saida_temp.append('1('+regiao_alinhada_sequencia[i]+')') # um significa mismatch. ALem disso, reporta qual o nucleotideo na sequencia dessa posicao que nao deu match


		elif primer_sequencia[i] == 'D':
			if regiao_alinhada_sequencia[i] == 'A' or regiao_alinhada_sequencia[i] == 'G' or regiao_alinhada_sequencia[i] == 'T': 
				saida_temp.append('0') # zero significa match
			else:
				saida_temp.append('1('+regiao_alinhada_sequencia[i]+')') # um significa mismatch. ALem disso, reporta qual o nucleotideo na sequencia dessa posicao que nao deu match


		elif primer_sequencia[i] == 'H':
			if regiao_alinhada_sequencia[i] == 'A' or regiao_alinhada_sequencia[i] == 'C' or regiao_alinhada_sequencia[i] == 'T': 
				saida_temp.append('0') # zero significa match
			else:
				saida_temp.append('1('+regiao_alinhada_sequencia[i]+')') # um significa mismatch. ALem disso, reporta qual o nucleotideo na sequencia dessa posicao que nao deu match


		elif primer_sequencia[i] == 'V':
			if regiao_alinhada_sequencia[i] == 'A' or regiao_alinhada_sequencia[i] == 'C' or regiao_alinhada_sequencia[i] == 'G': 
				saida_temp.append('0') # zero significa match
			else:
				saida_temp.append('1('+regiao_alinhada_sequencia[i]+')') # um significa mismatch. ALem disso, reporta qual o nucleotideo na sequencia dessa posicao que nao deu match

		else:
			saida_temp.append('1('+regiao_alinhada_sequencia[i]+')') # um significa mismatch. ALem disso, reporta qual o nucleotideo na sequencia dessa posicao que nao deu match

	return saida_temp

### -- uso do script de forma independente, nao como modulo de outro script -- ###

if __name__ == '__main__':
	import sys
	if len(sys.argv) < 5:
		print(usage)
		exit()
	else:
		if '--sequences' in sys.argv:	
			argSequences = sys.argv.index('--sequences')+1
			file = open(sys.argv[argSequences],'r')
			arq_fasta = file.readlines()
			file.close()
			
			import formataFasta_17Mar2020 #formata o arquivo multifasta.
			fasta_selected = formataFasta_17Mar2020.SeqOneLine(arq_fasta)

		else:
			print(usage)
			exit()

		if '--primers' in sys.argv:
			argPrimers = sys.argv.index('--primers')+1
			file = open(sys.argv[argPrimers],'r')
			primers = file.readlines()
			file.close()
			
			primer_foward = primers[0]+primers[1] # o primer foward é considerado a sempre a primeira sequencia do arquivos fasta de primers. Já tem o \n apos o ID da sequencia, logo, nao preciso incluir outro.
			seq_foward = list((primers[1].rstrip()).upper()) # sequencia do primer foward. Vai ser usado para analisar quais as poscioes de match/mismatches/gaps entre o primer e a sequencia alvo
			file = open('primer_foward.fasta', 'w')
			file.writelines(primer_foward)
			file.close()

			if not primers[2].startswith('>'):
				print('\nERRO: a sequencia dos primers NAO podem ser separadas por quebras de linha!\n')
				print(usage)
				exit()
			else:
				primer_reverse = primers[2]+primers[3] # não precisa fazer a complementar reversa, pois o exonerate alinha a RevComp também.
				seq_reverse = list((primers[3].rstrip()).upper()) # sequencia do primer reverse eh o elemento primer[3]. Vai ser usado para analisar quais as poscioes de match/mismatches/gaps entre o primer e a sequencia alvo
				file = open('primer_reverse.fasta', 'w')
				file.writelines(primer_reverse)
				file.close()

			if len(primers) > 4: # significa que tambem temos a probe para alinhar.
				probe = True
				probe = primers[4]+primers[5]
				seq_probe = list((primers[5].rstrip()).upper()) # sequencia da probe eh o elemento primer[5]. Vai ser usado para analisar quais as poscioes de match/mismatches/gaps entre o primer e a sequencia alvo
				file = open('probe.fasta', 'w')
				file.writelines(probe)
				file.close()
			else:
				probe = False # contra se deve ou tentar alinhamento da probe com a sequencia.

			#print(primer_foward)
			#print(primer_reverse)
			#print(probe)

		else:
			print(usage)
			exit()

		### --- chamada da funcao loopAlinhamento --- ###

		alinhamentos_primers, alinhamentos_primer_reverse, alinhamentos_probe, sequencias_amplificacao = loopAlinhamento(fasta_selected)

		# remove arquivos fasta dos primers gerados anteriormente
		subprocess.run(['rm','primer_reverse.fasta'])
		subprocess.run(['rm','primer_foward.fasta'])
		if probe:
			subprocess.run(['rm','probe.fasta'])

		## -- imprimir resultados -- ##

		# pegar o nome do arquivo fasta de primers.
		if '/' in str(sys.argv[argPrimers]):
			nomeArqFastaPrimers = str(sys.argv[argPrimers]).split('/')[-1]
		else:
			nomeArqFastaPrimers = str(sys.argv[argPrimers])

		# pegar o nome do arquivo fasta das sequencias.
		if '/' in str(sys.argv[argSequences]):
			nomeArqFastaSequencias = str(sys.argv[argSequences]).split('/')[-1]
		else:
			nomeArqFastaSequencias = str(sys.argv[argSequences])

		file = open('alinhamentos_primer_foward_'+nomeArqFastaPrimers+'_vs_'+nomeArqFastaSequencias+'.aligned','w')
		# salvar a sequencia do primer foward
		seq_foward = '\t'.join(seq_foward) # converter a lista em string, com os elementos separados por tabulação.
		ID_primer_foward = primers[0].lstrip('>').rstrip() # ID do primer foward.
		file.writelines(ID_primer_foward+'\t'+seq_foward+'\n')
		for linha in alinhamentos_primers:
			# linha eh uma lista iniciada com o ID da sequencia que alinhou no primer, seguida de 0 e 1 para match e mismatches, repectivamente.
			linha = '\t'.join(linha)
			file.writelines(linha+'\n')
		file.close()

		file = open('alinhamentos_primer_reverse_'+nomeArqFastaPrimers+'_vs_'+nomeArqFastaSequencias+'.aligned','w')
		# salvar a sequencia do primer reverse
		seq_reverse = '\t'.join(seq_reverse) # converter a lista em string, com os elementos separados por tabulação.
		ID_primer_reverse = primers[2].lstrip('>').rstrip() # ID do primer foward.
		file.writelines(ID_primer_reverse+'\t'+seq_reverse+'\n')
		for linha in alinhamentos_primer_reverse:
			# linha eh uma lista iniciada com o ID da sequencia que alinhou no primer, seguida de 0 e 1 para match e mismatches, repectivamente.
			linha = '\t'.join(linha)
			file.writelines(linha+'\n')
		file.close()

		if probe:
			file = open('alinhamentos_primer_probe_'+nomeArqFastaPrimers+'_vs_'+nomeArqFastaSequencias+'.aligned','w')
			# salvar a sequencia do primer probe
			seq_probe = '\t'.join(seq_probe) # converter a lista em string, com os elementos separados por tabulação.
			ID_primer_probe = primers[4].lstrip('>').rstrip() # ID do primer foward.
			file.writelines(ID_primer_probe+'\t'+seq_probe+'\n')
			for linha in alinhamentos_probe:
				# linha eh uma lista iniciada com o ID da sequencia que alinhou no primer, seguida de 0 e 1 para match e mismatches, repectivamente.
				linha = '\t'.join(linha)
				file.writelines(linha+'\n')
			file.close()

		# salvar seqeuncias fasta da regiao amplificada.
		file = open('sequencias_regiao_amplificacao_'+nomeArqFastaPrimers+'_vs_'+nomeArqFastaSequencias+'.fasta', 'w')
		for linha in sequencias_amplificacao:
			file.writelines(linha+'\n')
		file.close()

		## zipar arquivos de saida do exonerate
		nome_resultados_foward = 'saida_exonerate_primer_foward_'+nomeArqFastaPrimers+'_vs_'+nomeArqFastaSequencias+'.zip'
		subprocess.run(['zip temp.zip *foward_primer_aligned'], shell=True)
		subprocess.run(['mv', 'temp.zip', nome_resultados_foward])
		subprocess.run(['rm *foward_primer_aligned'], shell=True)

		nome_resultados_reverse = 'saida_exonerate_primer_reverse_'+nomeArqFastaPrimers+'_vs_'+nomeArqFastaSequencias+'.zip'
		subprocess.run(['zip temp.zip *reverse_primer_aligned'], shell=True)
		subprocess.run(['mv', 'temp.zip', nome_resultados_reverse])
		subprocess.run(['rm *reverse_primer_aligned'], shell=True)

		if probe:
			nome_resultados_probe = 'saida_exonerate_probe_'+nomeArqFastaPrimers+'_vs_'+nomeArqFastaSequencias+'.zip'
			subprocess.run(['zip temp.zip *probe_primer_aligned'], shell=True)
			subprocess.run(['mv', 'temp.zip', nome_resultados_probe])
			subprocess.run(['rm *probe_primer_aligned'], shell=True)

