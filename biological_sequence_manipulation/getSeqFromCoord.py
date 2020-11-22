"""
==================
getSeqFromCoord.py
==================

PT-BR:
Script para extrair sequencias fasta a partir de um arquivo GFF3 de coordenadas. Além disso, é possível extrair de sequêncas de um arquivo de coordenadas genérico, bastando ajustar quais colunas estão as coordenadas ou outros atributos importantes como a fita (+/-).

EN:
Script to extract fasta sequences from a GFF3 file. In addition, it is possible to extract sequences from a generic coordinate file by simply adjusting which columns bear the coordinates or other important attributes such as the strand (+/-).

Type python3 getSeqFromCoord.py -help to see the instructions.

"""

print(__doc__)

usage = "\tgetSeqFromCoord.py\n\n\
		Script para extrair sequencias a partir de um arquivo de coordenadas. É possível extrair de sequêncas de qualquer tipo de arquivo de coordenadas, bastando ajustar quais colunas estão as coordenadas ou outros atributos importantes como a fita (+/-)\n\n\
		*** Argumento obrigatório\n\
		--input <arquivo> (arquivo com coordenadas)\n\
		--fasta <arquivo> (arquivo fasta para extrair as sequencias de acordo com as coordenadas)\n\
		*** Argumentos opicionais\n\
		--coordenadas <int,int> (numeros inteiros referentes as colunas das coordenadas inicial e final, nessa ordem). Default: 3,4 (arquivo GFF3, iniciando de 0)\n\
		--fita <int|none> (numero inteiro referente a coluna da fita OU 'none' caso todas as features estejam na fita positiva (5'3'), isto é, a fita não é relevante). Default: 6 (arquivo GFF3, iniciando de 0)\n\n"

import sys, re, subprocess, comparaCoordenadas, formataFasta

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

	if '--fasta' in sys.argv:
		argFASTA = sys.argv.index('--fasta')+1
		file = open(sys.argv[argFASTA],'r')
		fileFASTA = file.readlines()
		file.close()
	else:
		print(usage)
		exit()

	if '--coordenadas' in sys.argv:
		argCOORD = sys.argv.index('--coordenadas')+1
		coordenadas = str(sys.argv[argCOORD])
		col_coord_inicial = int(coordenadas.split(',')[0])-1
		col_coord_final = int(coordenadas.split(',')[1])-1
	else:
		col_coord_inicial = 3
		col_coord_final = 4

	if '--fita' in sys.argv:
		argFITA = sys.argv.index('--fita')+1
		col_fita = str(sys.argv[argFITA])
		if col_fita == 'none':
			pass
		else:
			col_fita = int(sys.argv[argFITA])
	else:
		col_fita = 6

	fasta_formatado = formataFasta.SeqOneLine(fileFASTA)

	arquivo = open(str(sys.argv[argINPUT])+'.coord.seqs', 'w')

	contig = ''
	for i,record in enumerate(fileINPUT):
		if record == '\n':
			continue
		record = record.rstrip()
		record_ID = record.split('\t')[0]
		coord_inicial = record.split('\t')[col_coord_inicial]
		coord_final = record.split('\t')[col_coord_final]
		if col_fita == 'none':
			fita = '+'
		else:
			fita = record.split('\t')[col_fita]
		
		if contig != record_ID:
			sequencia = ''
			for ID in fasta_formatado.keys():
				ID_splitted = ID.split()[0] # Os programas preditores de TE ou tandem repeats podem quebrar o ID no espaço em branco (fazem truncagem de ID). Como o modulo formatadaFasta_17Mar2020 pega o nome completo do contig, tenho que fazer essa busca tambem
				if ID == record_ID or ID == '>'+record_ID:
					sequencia = fasta_formatado[ID]
					break
				elif ID_splitted == record_ID or ID_splitted == '>'+record_ID:
					sequencia = fasta_formatado[ID]
					break
			contig = record_ID

		# em alguns casos, features na fita reversa (-) tem coordenadas invertidas. Basta mudar suas posições de lugar para extrair a sequencia.
		if int(coord_inicial) > int(coord_final):
			subsequencia = sequencia[int(coord_final)-1:int(coord_inicial)]
		else:
			subsequencia = sequencia[int(coord_inicial)-1:int(coord_final)]
		if fita == '-': # se a feature esta na fita antisenso devemos fazer a complementar reversa
			subsequencia = subsequencia[::-1] # ...fazemos a sequencia reversa...
			subsequencia_rev = ''
			for nuc in subsequencia: #... e a complementar.
				if nuc.lower() == 'a':
					subsequencia_rev += 't'
				elif nuc.lower() == 't':
					subsequencia_rev += 'a'
				elif nuc.lower() == 'c':
					subsequencia_rev += 'g'
				elif nuc.lower() == 'g':
					subsequencia_rev += 'c'
			subsequencia = subsequencia_rev.upper()

		arquivo.writelines(">"+record_ID+'_coordenadas:'+coord_inicial+'-'+coord_final+'[fita:'+fita+']\n'+subsequencia+'\n')
	
	arquivo.close()