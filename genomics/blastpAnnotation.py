"""
===================
blastpAnnotation.py
===================

PT-BR: Script para anotacao funcional de proteinas preditas. Analisa todos os hits de uma analise com BLASTp ou DIAMOND e recupera a descricao dos subjects que esta disponivel ao utilizar parametro "stitle" no argumento -outfmt/--outfmt de execucao do BLAST/DIAMOND. Realiza concatenacao (merge) de todas as descricoes de subject para cada query de forma a definir a anotacao mais provavel (majority consensus).

EN: Script for functional annotation of predicted proteins. It analyzes all the hits of an analysis with BLASTp or DIAMOND and retrieves the description of the subjects that is available when using parameter "stitle" in the argument -outfmt /--outfmt  of BLAST/DIAMOND. Performs concatenation (merge) of all subject descriptions for each query in order to define the most likely annotation (majority consensus).

Type python3 blastpAnnotation.py -help to see the instructions.

"""

print(__doc__)

usage = '\n\tEste programa analisa todos os hits de uma analise com BLASTp ou DIAMOND. Recupera a descricao dos subjects que esta disponivel ao utilizar parametro "stitle" no argumento -outfmt/--outfmt de execucao do BLAST. Realiza concatenacao (merge) de todas as descricoes de subject para cada query de forma a definir a anotacao mais provavel (majority consensus)\n\n\
	Usar uma linha de comando do BLASTp semelhante a: blastp -query proteinas.faa -db database.faa -outfmt "6 qseqid sseqid qstart qend pident length slen evalue sstart send stitle qcovs" -evalue 1e-10 -num_threads 10 > output.txt\n\n\
	Usar uma linha de comando do DIAMOND semelhante a: diamond blastp --query proteinas.faa --db database.faa --outfmt 6 qseqid sseqid qstart qend pident length slen evalue sstart send stitle qlen --evalue 1e-10 --query-cover 60 --id 40 --max-target-seqs 200 --threads 10 > output.txt\n\n\
		AVISO: Usar preferencialmente o DIAMOND, pois eh mais rapido que o BLASTp e podemos reportar resultados filtrados com relacao a cobertura da query (--query-cover) e porcentagem de identidade (--id)\n\n\
USAGE: python3 blastpAnnotation.py --input <arquivo_saida_blast.txt> --max_descriptions <integer> --output <nome_arquivo_output>\n\n\
	*** Argumentos:\n\
		--input (arquivo de saida do BLASTp/DIAMOND no formato outfmt 6 com "stitle" na coluna 11)\n\
		--max_descriptions (numero maximo de anotacoes reportadas para cada proteina). Default: 5\n\
		--output (nome do arquivo para salvar os resultados de anotacao)\n\
	'

# definir palavaras e descrições não importantes para anotação. A funcao abaixo recupera termos muito abundantes e que por isso nao sao interessantes para anotacao (termos genericos, como "precusor" ou "component")
trash = ['hypothetical', 'protein', 'variant', 'precursor', 'component', 'partial', 'putative', 'uncharacterized', 'containing', 'predicted', 'like', 'conserved', 'altname:', 'related', 'to', 'type', 'unnamed', 'product', 'repeat', 'general', 'and', 'associated', 'probable', 'possible', 'contig', 'domain', 'superfamily', 'family', 'binding', 'module','factor'] # algumas descricoes/palavras genericas de subjects.

def genericWords(saida_blast):

	import re, operator
	from collections import Counter

	descriptions = []
	for i,linha in enumerate(saida_blast):

		if linha.startswith('sp|'): # descricao referente a proteina do swiss-prot
			match = re.search(r'(?<=Full=).+(?=;|\\n)',linha) # match com todos os caracteres entre "Full=" e ";" ou "\n" (quebra de linha, pois encontrei descricoes do swiss-prot nao tinham ';' depois da descricao)
			if match: # tive que colocar essa condicao, pois econtrei IDs sem descricoes! Dava erro pois tentanva a funcao group() onde nao existiu match com a funcao re.search()
				descriptions.append(match.group())

		elif linha.startswith('tr|'): # descricao referente ao TrEMBL
			match = re.search(r'(?<= ).+(?= OS=)',linha) # match com todos os caracteres entre " " e " OS="
			if match:
				descriptions.append(match.group())

		elif re.search(r'(?<= )[A-Z]{2}_',linha): # descricao referente ao RefSeq
			match = re.search(r'(?<= ).+(?= \[)',linha) # match com todos os caracteres entre " " e " ["
			if match:
				descriptions.append(match.group()) 
			
		elif linha.startswith('>pdb'): # descricao referente ao PDB
			print('A recuperacao de descricoes do PDB nao esta configurada. Ignorando a descricao: '+linha[1:])
	
	words = (' '.join(descriptions).lower()) # string de descricoes (eh enorme). Converti tudo em minusculo para ser case insensitive.
	words = re.sub(r'(?<= )[0-9]+[,\-]*[0-9]*(?= )','',str(words)) # substitui palavras que sao numeros por 'nada'.  Nao deve contabilizar numeros como palavaras genericas (remove por exemplo o numero 3 e 3,3 e 33 e 33,3 e 33-3 e assim por diante)
	words = re.split(r'\s|\-|/|=|,|;|\)|\(', words) # transforma a string em lista de cada palavra (nao sao mais descricoes, apenas palavras independentes)
	
	words = dict(Counter(words)) # transforma a lista em um dicionario, onde a chave sao as palavras e o valor eh o numero de vezes que apareceu na lista.

	for i,chave in enumerate(words.keys()): #com esse for, vou buscar chave 'vazia' e remove-la do dicionario.
		if chave == '':
			del words[chave]
			iteracoes = i
			break

	wordsOrdenado = sorted(words.items(), key=lambda kv: kv[1]) # ordena o dicionario de acordo com o valor de cada chave (menor para maior)
	wordsOrdenado = wordsOrdenado[::-1] # reverter a ordem do dicionario, pois os pares (chave,valor) estao em ordem crescente e quero decrescente.

	# inspecao manual de quais palavras sao ou podem ser consideradas genericas. Com isso, posso estabelecer um corte e usa-lo automaticamente em proximas analises.
	file = open('termosDescricoesProteinas.txt','w')
	for linha in wordsOrdenado:
		file.writelines(str(linha)+'\n')
	file.close()

	# remover termos muito especificos ou termos muito genericos.
	# os termos genericos estao listados na variavel 'trash' do inicio do script.
	# os termos muitos especificos aparecem 1 ou 2 vezes na variavel 'wordsOrdenado'
	file = open('termosDescricoesProteinasFiltrado.txt','w') # esse arquivo de saida contem as descricoes "boas" para anotacao funcional das proteinas preditas.
	wordsOrdenadoFiltrado = []
	for elemento in wordsOrdenado: # percorre cada tupla de 'wordsOrdenado'. O segundo elemento de cada tupla eh quantas vezes determinada palavra aparece.
		#print(elemento)
		if elemento[1] <= 2 or (elemento[0].lower()) in trash: # remover termos muito especificos. Esses termos aparecem apenas 1 ou 2 vezes no banco de dados inteiro. Nao sao interessantes par anotacao (geralmente nao sao informativos de funcionalidade). Alem disso, remove termos que aparecem na lista 'trash'
			continue
		elif len(elemento[0]) > 2 and (elemento[0].lower()) not in trash: # o comprimento da palavra deve ser maior que 2. Menos que isso eh pouco informativo (inclusive pode ser algum 'numero' (e.g., 3)). Alem disso, nao fazer uso de termos que estao na lista 'trash', pois sao muito genericos.
			wordsOrdenadoFiltrado.append(elemento)
			file.writelines(str(elemento)+'\n')
	file.close()

# funcao para recuperar descricoes mais abundantes das query (proteinas de interesse).
def blastpAnnotation(blastpResult, maxDespcritions, nomeOutput):

	import re, operator
	from collections import Counter

	stitles = [] # salva as descricoes de subjects.

	for i,elemento in enumerate(linhas):

		if elemento == '\n': # nao trabalhar com linhas vazias.
			continue

		elemento = elemento.rstrip('\n')
		elemento = re.split('\t',elemento)
		
		if not i: #comecar a salvar ID_query e stitles. So entra aqui na primeira linha.
			ID_query = elemento[0] # identificador do query
			stitles.append(elemento[10]) # descricao do subject do arquivo de saida do blast (contem tambem o ID do subject)
			
		elif elemento[0] == ID_query and not i+1 == len(linhas): #esse elif eh para continuar salvando 'stitles' que tem a mesma query. Caso esteja no fim do arquivo, nao entra nesse elif, deve-se entrar no proxima para imprimir os resultados dessa query.
			ID_query = elemento[0] # identificador do query
			stitles.append(elemento[10]) # descricao do subject do arquivo de saida do blast (contem tambem o ID do subject)

		elif elemento[0] != ID_query or i+1 == len(linhas): # se mudar a query ou chegar no fim do arquivo de linhas (BLASTp), vamos analisar as descricoes da query que estavamos trabalhando.

			if i+1 == len(linhas): #caso entre no elif acima pois chegou no fim do arquivo, devemos salvar o 'stitles' dessa ultima linha.
				stitles.append(elemento[10])
				
			file = open(nomeOutput, 'a')
			file.writelines(ID_query+'\t') # imprimir o identificador da proteina query no arquivo de saida. Depois colocar a descricao do BLAST na frente (vai ser separado por tabulacao)
			file.close()
			
			ID_query = elemento[0] # se entrou no 'if' pela condicao elemento[0] diferent de ID_query, devemos fazer a operacao descrita para salvar a nova query e pegar seus subjects.
			descriptions = []
			for k,subject in enumerate(stitles): # vamos percorrer a lista 'stitles' que tem as descricoes cruas dos subjects. Vamos refina-las.
				if subject.startswith('sp'): # descricao referente a proteina do swiss-prot
					match = re.search(r'(?<=Full=).+(?=;|\\n)', subject) # match com todos os caracteres entre "Full=" e ";" ou "\n" (quebra de linha, pois encontrei descricoes do swiss-prot nao tinham ';' depois da descricao)
					if match: # tive que colocar essa condicao, pois econtrei IDs sem descricoes! Dava erro pois tentanva a funcao group() onde nao existiu match com a funcao re.search()
						descriptions.append(match.group())

				elif subject.startswith('tr|'): # descricao referente ao TrEMBL
					match = re.search(r'(?<= ).+(?= OS=)', subject) # match com todos os caracteres entre " " e " OS="
					if match:
						descriptions.append(match.group())
			
				elif re.search(r'(?=)[A-Z]{2}_', subject): # descricao referente ao RefSeq
					match = re.search(r'(?<= ).+(?= \[)', subject) # match com todos os caracteres entre " " e " ["
					if match:
						descriptions.append(match.group())


			words = (' '.join(descriptions).lower()) # string de descricoes (eh enorme). Converti tudo em minusculo para ser case insensitive.
			words = re.sub(r'[\]\[\)\()\+]+','', words) # algumas descricoes de subject continham metacaracteres ][)(+ que estavam dando problema no modulo de regex do python3.6 (o erro ocorria na linha 168 do meu script)
			words = re.sub(r'(?<= )[0-9]+[,\-]*[0-9]*(?= )','', str(words)) # substitui palavras que sao numeros por 'nada'.  Nao deve contabilizar numeros como palavaras genericas (remove por exemplo o numero 3 e 3,3 e 33 e 33,3 e 33-3 e assim por diante)
			words = re.split(r'\s|\-|/|=|,|;|\)|\(', words) # transforma a string em lista de cada palavra (nao sao mais descricoes, apenas palavras independentes)
			words = dict(Counter(words)) # transofrma a lista em um dicionario, onde a chave sao as palavras e o valor eh o numero de vezes que apareceu na lista.
			#print(words)
			for i,chave in enumerate(words.keys()): #com esse for, vou buscar chave 'vazia' e remove-la do dicionario.
				if chave == '':
					del words[chave]
					break

			wordsOrdenado = sorted(words.items(), key=lambda kv: kv[1]) # ordena o dicionario de acordo com o valor de cada chave (menor para maior)
			wordsOrdenado = wordsOrdenado[::-1] # reverter a ordem do dicionario, pois os pares (chave,valor) estao em ordem crescente e quero decrescente.
			
			# remover termos muito especificos ou termos muito genericos.
			# os termos genericos estao listados na variavel 'trash' do inicio do script.
			# os termos muitos especificos aparecem 1 ou 2 vezes na variavel 'wordsOrdenado'
			wordsOrdenadoFiltrado = []
			for tupla in wordsOrdenado: # percorre cada tupla de 'wordsOrdenado'. O segundo elemento de cada tupla eh quantas vezes determinada palavra aparece.
				#print(tupla)
				if tupla[1] <= 2 or (tupla[0].lower()) in trash or 'hypothetical' in (tupla[0].lower()): # remover termos muito especificos. Esses termos aparecem apenas 1 ou 2 vezes no banco de dados inteiro. Nao sao interessantes par anotacao (geralmente nao sao informativos de funcionalidade). Alem disso, remove termos que aparecem na lista 'trash'
					#print(tupla)
					continue
				elif len(tupla[0]) > 2 and (tupla[0].lower()) not in trash: # o comprimento da palavra deve ser maior que 2. Menos que isso eh pouco informativo (inclusive pode ser algum 'numero' (e.g., 3)). Alem disso, nao fazer uso de termos que estao na lista 'trash', pois sao muito genericos.
					wordsOrdenadoFiltrado.append(tupla)
		
			stitles.clear() # ja usei o que precisava da lista 'stitles'

			descricao_final = [] # dicionario final de anotacoes
			
			for n,valor in enumerate(wordsOrdenadoFiltrado): # a tupla 'wordsOrdenadoFiltrado' contem cada palavra e o numero de vezes que ela aparece (no formato chave, valor do dicionario). Esta em ordem decrescente de numero de vezes que aparece

				if n < int(maxDespcritions) and not valor[0].isnumeric(): # quero analisar as 5 palavras mais abundantes das descricoes, logo pego os 5 primeiros elementos da tupla 'wordsOrdenadoFiltrado'. Nao sei se o valor 5 é bom.

					for k, subject in enumerate(descriptions): # pegar cada frase dentro de uma lista de frases de descricoes

						if valor[0] in subject.lower() and not re.search(valor[0], (str(descricao_final).lower())): #analisa se uma palavra esta dentro da descricao maior (subject) e se essa palavra (valor[0]) nao ira induzir redundancias para lista 'descricao_final'. Isso previne redundancias.				 
							descricao_final.append(subject) # adiciona a descricao (subject) como chave do dicionario e o numero de vezes que a palavra que deu origem (valor[0]) a essa acao como o valor.
							break #paralisa o 'for' para nao ficar repetindo o salvamento de descricoes a partir de uma palavra
			
			if not descricao_final: # caso nenhuma descricao seja encontrada...
				descricao_final.append('hypothetical protein')

			descricao_final_string = '; '.join(descricao_final) # cnoverte a lista em string.

			file = open(nomeOutput, 'a')
			file.writelines(descricao_final_string)
			file.writelines('\n')
			file.close()

			stitles.append(elemento[10]) # descricao do subject do ID atual

if __name__ == '__main__':
	import sys
	if len(sys.argv) < 3:
		print(usage)
		exit()

	else:
		if '--input' in sys.argv:
			argInput = sys.argv.index('--input')+1
			file = open(sys.argv[argInput],'r')
			linhas = file.readlines()
			file.close()
		else:
			print(usage)
			exit()

		if '--output' in sys.argv:
			nomeOutput = str(sys.argv[sys.argv.index('--output')+1])
		else:
			print(usage)
			exit()

		if '--max_descriptions' in sys.argv:
			maxDespcritions = int(sys.argv[sys.argv.index('--max_descriptions')+1])
		else:
			maxDespcritions = 5

		genericWords(linhas) # chamada da funcao para contar descricoes de subject mais abundantes no arquivo de saida. Como sao muito comuns ou genericas, nao sao boas para classfiicar funcionalidade das proteinas.
		blastpAnnotation(linhas, maxDespcritions, nomeOutput) # chamada da funcao para encontrar as descricoes de subject que mais aparecem para cada query. Mas excluindo descricoes promiscuas (aparecem muito no arquivo de saida do blast)