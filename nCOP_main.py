# Incluindo bibliotecas
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import pandas as pd
import time
import random
import copy
import json
import collections
import sys

InformationsToFinalReport = {}

###########################################################################################################################

def IniciaAlgoritmo(networkPath, mutationalFilePath, weightFilePath, startAlpha, outputNameFile):

    # Gerando a rede biológica
    Graph = nx.read_edgelist(networkPath, delimiter = ' ', data = ('pacients',[]))
    for i in Graph.nodes:
        Graph.nodes[i]['pacients'] = []
    print("Tamanho da rede: " + str(len(Graph.nodes)))
    
    InformationsToFinalReport["lenghtOfNetwork"] = str(len(Graph.nodes))
    
    #Abrindo o arquivo de Genes mutados por pacientes
    with open(mutationalFilePath,'r') as f:
        lines = f.readlines()

    # Criando uma tabela de genes por pacientes
    InformationsToFinalReport["numberOfMutatedGenes"] = len(lines)
    InformationsToFinalReport["numberOfMutatedGenesNotInGraph"] = 0
    TabelaPacienteGene = []
    for line in lines:
        gene = line.split()[0]
        
        if gene not in Graph: InformationsToFinalReport["numberOfMutatedGenesNotInGraph"] += 1
        
        lista = line.split()[1:]
        for patient in lista:
            index = next((index for (index, d) in enumerate(TabelaPacienteGene) if d["patient"] == patient), None)
            if index == None:
                TabelaPacienteGene.append({'patient':patient, 'mutations': [gene]})
            else:
                TabelaPacienteGene[index]['mutations'] += [gene]
    
    InformationsToFinalReport["numberOfPatients"] = len(list(TabelaPacienteGene))
    # Criando uma tabela com o tamanho dos genes
    TabelaGeneTamanho = {}
    
    if(weightFilePath != None):
        print("Utilizando o arquivo de pesos...")
        arqWeights = open(weightFilePath,'r')
        lines = arqWeights.readlines()
        for line in lines:
            gene = line.split()[0]
            weight = float(line.split()[1])
            if gene in Graph:
                TabelaGeneTamanho[gene] = weight
        arqWeights.close()
    else:
        for gene in Graph:
            TabelaGeneTamanho[gene] = 1


    print("Comecando o processo de encontrar a constante de nomalização...")
    
    normalizationConstant = DescobrirMaiorSubGrafo(Graph, TabelaPacienteGene, TabelaGeneTamanho)

    print("Constante de normalização: " + str(normalizationConstant))

    InformationsToFinalReport["normalizationConstant"] = normalizationConstant

    print("Ajustando pesos de acordo com a constante de normalização...")
    for gene in TabelaGeneTamanho:
        TabelaGeneTamanho[gene] /= normalizationConstant

    # Calcula o melhor alpha
    
    bestAlpha = startAlpha
    InformationsToFinalReport["bestAlpha"] = "None"
    if(bestAlpha == None):
    
        print("Comecando o processo para descobrir o melhor alpha...")

        bestAlpha = DescobrirMelhorAlpha(Graph, TabelaPacienteGene, TabelaGeneTamanho)
        
        InformationsToFinalReport["bestAlpha"] = bestAlpha
    # Calcula 1000 interacoes
    print("Comecando o processo de priorização dos genes...")

    frequency = AlgoritmoFinal(Graph, TabelaPacienteGene, bestAlpha, TabelaGeneTamanho)

    if(outputNameFile == None): outputNameFile = "output"
    print("Escrevendo os genes no arquivo de saida...")
    output = open(outputNameFile + ".txt",'w')
    for item in frequency:
        key = list(item.keys())[0]
        output.write("{: <12} {: <12} \n".format(key, "{:.2f}".format(item[key]) + '%'))
    output.close()

###########################################################################################################################

def SorteiaNoInicial(GraphModify):
    # Gera subgrafo com o melhor de todos (aleatorizar para os 5 melhores)
    ordem = [0,0,0,0,0]
    bestNodes = ["","","","",""]
    AuxiliarGraph = copy.deepcopy(GraphModify)

    for i in range(0,5,1):
        bestNode = ""
        for node in AuxiliarGraph:
            num = len(AuxiliarGraph.nodes[node]['pacients'])
            if (num > ordem[i]):
                bestNodes[i] = node
                ordem[i] = num
        AuxiliarGraph.remove_node(bestNodes[i])
    
    randNumber = random.randint(0,np.sum(ordem))
    escolhido = 0
    soma = 0
    soma += ordem[escolhido]
    while(randNumber > soma):
        escolhido += 1
        soma += ordem[escolhido]
    return bestNodes[escolhido]


###########################################################################################################################

def DescobrirMaiorSubGrafo(Graph, TabelaPacienteGene, TabelaGeneTamanho):
    
    Conjunto = random.sample(TabelaPacienteGene,len(TabelaPacienteGene))

    GraphModify = copy.deepcopy(Graph)
    
    for dicts in Conjunto:
        patient = [dicts['patient']]
        for gene in dicts['mutations']:
            if gene in GraphModify:
                GraphModify.nodes[gene]['pacients'] += patient
    
    soma = 0
    nVezes = 10
    for i in range(0,nVezes,1):
        
        GraphModify = copy.deepcopy(Graph)
    
        for dicts in Conjunto:
            patient = [dicts['patient']]
            for gene in dicts['mutations']:
                if gene in GraphModify:
                    GraphModify.nodes[gene]['pacients'] += patient
                
                
        nodesQueApareceram = AlgoritmoV1_12(1, GraphModify, TabelaGeneTamanho)
        
        sizeOfGh = 0
        for gene in nodesQueApareceram:
            sizeOfGh += TabelaGeneTamanho[gene]
        soma += sizeOfGh
    soma /= nVezes
    #print("Constante de normalizacao: " + str(sizeOfGh))
    return soma


###########################################################################################################################


def AlgoritmoV1_12(alpha, GraphModify, TabelaGeneTamanho):
    ini = time.time()

    noInicial = SorteiaNoInicial(GraphModify)

    Gh = []
    Gh.append(noInicial)
    
    # Descobre quantos pacients existem no total no grafo
    listPatientsInGraph = []
    for node in GraphModify:
        listPatientsInGraph += list(set(GraphModify.nodes[node]['pacients'])-set(listPatientsInGraph))
    maxPatients = len(listPatientsInGraph)
    sizeOfAllGraph = len(GraphModify.nodes)

    antigaMelhorNota = np.inf

    while(True):

        # Criar lista com os pacientes cobertos por Gh na interação atual
        listPatientsInGh = []
        nodesGraphModify = GraphModify.nodes
        for node in Gh:
            listPatientsInGh += list(set(nodesGraphModify[node]['pacients']) - set(listPatientsInGh))
        patientsInGh = len(listPatientsInGh)
        
        # Define o tamanho de Gh na interação atual
        sizeOfGh = 0
        for gene in Gh:
            sizeOfGh += TabelaGeneTamanho[gene]
        
        if(antigaMelhorNota == np.inf):
            antigaMelhorNota = alpha*(1 - patientsInGh/maxPatients) + (1-alpha)*(sizeOfGh)
        
        menorNota = np.inf

        nodesNewSubgraph = []

        for gene in Gh:
            
            neighborsGene = GraphModify.neighbors(gene)
            
            for son in neighborsGene:
                
                if son in Gh:
                    continue
                
                coveredPatient = patientsInGh + len(set(nodesGraphModify[son]['pacients'])-set(listPatientsInGh))
                nota = alpha*(1 - coveredPatient/maxPatients) + (1-alpha)*(sizeOfGh + TabelaGeneTamanho[son])
            
                if(nota < menorNota):
                    menorNota = nota
                    nodesNewSubgraph = [son]
                else:
                    if(nota == menorNota):
                        nodesNewSubgraph += [son]

            
        if(menorNota != np.inf and menorNota < antigaMelhorNota):
            Gh = Gh + [nodesNewSubgraph[random.randint(0,len(nodesNewSubgraph)-1)]]
            #Gh = Gh + [nodesNewSubgraph[0]]
            antigaMelhorNota = menorNota;
            continue
        ##################### FIM DA ANALISE DE DISTANCIA 1 #####################
        
        nodesNewSubgraph = []
        for gene in Gh:
            
            neighborsGene = GraphModify.neighbors(gene)
            
            for son in neighborsGene:
                
                if son in Gh:
                    continue
                
                listPatientsInSon = nodesGraphModify[son]['pacients']

                neighborsOfSon = list(GraphModify.neighbors(son))

                listOfGrandSon = []

                listOfGrandSon = list(set(neighborsOfSon)-set(Gh)-set([son]))

                difference_son = len(set(listPatientsInSon)-set(listPatientsInGh))

                for grandSon in listOfGrandSon:

                    listPatientsInGrandSon = nodesGraphModify[grandSon]['pacients']

                    difference_grandSon = len(set(listPatientsInGrandSon)-set(listPatientsInSon)-set(listPatientsInGh))

                    nota = alpha*(1 - (patientsInGh + difference_son + difference_grandSon)/maxPatients) + (1-alpha)*(sizeOfGh + TabelaGeneTamanho[son] + TabelaGeneTamanho[grandSon])

                    if(nota < menorNota):
                        menorNota = nota
                        nodesNewSubgraph = [(son,grandSon)]
                    else:
                        if(nota == menorNota):
                            nodesNewSubgraph.append((son,grandSon))
            
            
            ######
        if(menorNota == np.inf or menorNota >= antigaMelhorNota):
            break

        randNumber = random.randint(0,len(nodesNewSubgraph)-1)
        Gh = Gh + [nodesNewSubgraph[randNumber][0]] + [nodesNewSubgraph[randNumber][1]]
        #Gh = Gh + [nodesNewSubgraph[0][0]] + [nodesNewSubgraph[0][1]]
        antigaMelhorNota = menorNota
        
    #print(antigaMelhorNota)
    return Gh



###########################################################################################################################

def InteracaoAlgoritmoFinal(Graph, TabelaPacienteGene, alpha, TabelaGeneTamanho):
    
    proportion = 0.85
    tamanhoLista = int(len(TabelaPacienteGene)*proportion)
    ConjuntoTreinamento = random.sample(TabelaPacienteGene,tamanhoLista)

    GraphModify = copy.deepcopy(Graph)
    
    for dicts in ConjuntoTreinamento:
        patient = [dicts['patient']]
        for gene in dicts['mutations']:
            if gene in GraphModify:
                GraphModify.nodes[gene]['pacients'] += patient
    
    nodesQueApareceram = []
    
    init = time.time()
    nodesQueApareceram = AlgoritmoV1_12(alpha, GraphModify, TabelaGeneTamanho)
    
    return nodesQueApareceram

###########################################################################################################################

def AlgoritmoFinal(Graph, TabelaPacienteGene, alpha, TabelaGeneTamanho):
    counterGenes = []
    nodesQueApareceram = []
    numberInterations = 1000
    for i in range(0,numberInterations,1):
        
        nodesQueApareceram += InteracaoAlgoritmoFinal(Graph, TabelaPacienteGene, alpha, TabelaGeneTamanho)
        
        if i%50 == 0:
            print("Concluido a iteração " + str(i))
    auxiliar = collections.Counter(nodesQueApareceram).most_common()
    frequency = []
    for item in auxiliar:
        frequency.append({item[0]:100*item[1]/numberInterations})
    return frequency

###########################################################################################################################

def InteracaoDescobrirAlpha(Graph, TabelaPacienteGene, alpha, TabelaGeneTamanho):

    ConjuntoTeste = random.sample(TabelaPacienteGene,int(len(TabelaPacienteGene)*0.10))
    resto = [elem for elem in TabelaPacienteGene if elem not in ConjuntoTeste]
    ConjuntoTreinamento = random.sample(resto,int(len(resto)*0.8))
    ConjuntoValidacao = [elem for elem in resto if elem not in ConjuntoTreinamento]

    GraphModify = copy.deepcopy(Graph)
    
    for dicts in ConjuntoTreinamento:
        for gene in dicts['mutations']:
            if gene in GraphModify:
                GraphModify.nodes[gene]['pacients'] += [dicts['patient']]
                
    nodesQueApareceram = AlgoritmoV1_12(alpha, GraphModify, TabelaGeneTamanho)

    coberturaConjuntoTreinamento = 0
    for dicts in ConjuntoTreinamento:
        if len(set(dicts['mutations']) & set(nodesQueApareceram)) != 0:
            coberturaConjuntoTreinamento += 1
    coberturaConjuntoTreinamento = 100*coberturaConjuntoTreinamento/len(ConjuntoTreinamento)
    
    coberturaConjuntoValidacao = 0
    for dicts in ConjuntoValidacao:
        if len(set(dicts['mutations']) & set(nodesQueApareceram)) != 0:
            coberturaConjuntoValidacao += 1
    coberturaConjuntoValidacao = 100*coberturaConjuntoValidacao/len(ConjuntoValidacao)
    
    coberturaConjuntoTeste = 0
    for dicts in ConjuntoTeste:
        if len(set(dicts['mutations']) & set(nodesQueApareceram)) != 0:
            coberturaConjuntoTeste += 1
    coberturaConjuntoTeste = 100*coberturaConjuntoTeste/len(ConjuntoTeste)
    
    return len(nodesQueApareceram), coberturaConjuntoTreinamento, coberturaConjuntoValidacao, coberturaConjuntoTeste

###########################################################################################################################

def DescobrirMelhorAlpha(Graph, TabelaPacienteGene, TabelaGeneTamanho):
    init = time.time()
    alphas = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]
    numInteracoes = 100
    NumNodes = []
    ConjuntoDeTreinamento = []
    ConjuntoDeValidacao = []
    ConjuntoDeTeste = []
    
    medianNumNodes = []
    medianConjuntoDeTreinamento = []
    medianConjuntoDeValidacao = []
    for i in alphas:
        medianNumNodes.append(0)
        medianConjuntoDeTreinamento.append(0)
        medianConjuntoDeValidacao.append(0)
    
    cursor = 0
    for alpha in alphas:

        for i in range(0,numInteracoes,1):

            a,b,c,d = InteracaoDescobrirAlpha(Graph, TabelaPacienteGene, alpha, TabelaGeneTamanho)
            if( i == 0 ):
                NumNodes.append([a])
                ConjuntoDeTreinamento.append([b])
                ConjuntoDeValidacao.append([c])
                ConjuntoDeTeste.append([d])
            else:
                NumNodes[cursor].append(a)
                ConjuntoDeTreinamento[cursor].append(b)
                ConjuntoDeValidacao[cursor].append(c)
                ConjuntoDeTeste[cursor].append(d)
                
            medianNumNodes[cursor] += a
            medianConjuntoDeTreinamento[cursor]  += b
            medianConjuntoDeValidacao[cursor]  += c
        print("Finalizado para o alpha = " + str(round(alpha,3)))
        cursor += 1
    
    medianNumNodes = [x/numInteracoes for x in medianNumNodes]
    medianConjuntoDeTreinamento = [x/numInteracoes for x in medianConjuntoDeTreinamento]
    medianConjuntoDeValidacao = [x/numInteracoes for x in medianConjuntoDeValidacao]
    
    bestAlpha = -1
    maxCoberturaValidacao = max(medianConjuntoDeValidacao)
    for i in range(0,len(alphas),1):
        if medianConjuntoDeValidacao[i]/maxCoberturaValidacao > 0.9:
            if abs(medianConjuntoDeTreinamento[i] - medianConjuntoDeValidacao[i]) > 0.05:
                bestAlpha = alphas[i]
                break
    print("Melhor alpha: " + str(bestAlpha))
    print("Tempo para achar o melhor alpha: " + str(time.time()-init))
    
    for i,alpha in enumerate(alphas):
        if alpha == bestAlpha:
            InformationsToFinalReport["medianConjuntoDeTreinamento"] = medianConjuntoDeTreinamento[i]
            InformationsToFinalReport["medianConjuntoDeValidacao"] = medianConjuntoDeValidacao[i]
            InformationsToFinalReport["medianNumNodes"] = medianNumNodes[i]
            InformationsToFinalReport["maxCoberturaValidacao"] = maxCoberturaValidacao
            break

    return bestAlpha

###########################################################################################################################

def Parser_MAF_File(Path, output_file='output'):
    
    '''
    The objective of this function is transformate one .maf file in a .txt file that serves as input to nCOP algorithm
    Path = the path of .maf file
    output_file = the name of the output file
    '''
    
    print ("Reading arquive...")
    dados = pd.read_csv(Path,comment='#', sep='\t',header=0,index_col=False, low_memory = False)
    
    
    print ("Creating new dataFrame with only essential columns")
    df = pd.concat([dados['Hugo_Symbol'],dados['Tumor_Sample_Barcode']],axis=1)
    
    print ("Formating output file...")
    
    output = open(output_file + '.txt','w')
    list_temp = []
    i = 0
    while i < len(df.index):

        temp = df['Hugo_Symbol'][i]

        while i < len(df.index) and temp == df['Hugo_Symbol'][i]:
            if(list_temp.count(df['Tumor_Sample_Barcode'][i]) == 0):
                list_temp.append(df['Tumor_Sample_Barcode'][i])
            i += 1
            
        output.write(temp + ' ' + ' '.join(list_temp) + '\n')
        list_temp = []
    
    print ("Done.")
    output.close()


###########################################################################################################################




if __name__ == "__main__":
    
    initTime = time.time()
    
    Arguments = sys.argv
    
    if(len(Arguments) <= 1):
        
        print("Error, no arguments")
        
    elif(len(Arguments) == 2):
        
        print("Error, only one argument passed")
        
    else:

        networkPath = str(Arguments[1])
        mutationalFilePath = str(Arguments[2])
        
        from datetime import datetime


        InformationsToFinalReport["initialDateAndHour"] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        
        
        InformationsToFinalReport["mutationFileWasConverted"] = "No"
        if(mutationalFilePath[-4:] == ".maf"):
            Parser_MAF_File(mutationalFilePath, mutationalFilePath[:-4])
            mutationalFilePath = mutationalFilePath[:-4] + ".txt"
            
            InformationsToFinalReport["mutationFileWasConverted"] = "Yes"
        
        
        weightFilePath = None
        
        if(len(Arguments) == 4 and Arguments[3] != "None"):
            weightFilePath = str(Arguments[3])
        
        startAlpha = None
        InformationsToFinalReport["startAlpha"] = None
        
        if(len(Arguments) == 5 and Arguments[4] != "None"):
            
            startAlpha = float(Arguments[4])
            
            if(startAlpha <= 0 or startAlpha > 1):
                print("Invalid alpha initial value")
                exit(-1)
            else:
                InformationsToFinalReport["startAlpha"] = startAlpha
        
        outputNameFile = None
        
        if(len(Arguments) == 6 and Arguments[5] != "None"):
            
            outputNameFile = str(Arguments[5])
        
        
        IniciaAlgoritmo(networkPath, mutationalFilePath, weightFilePath, startAlpha, outputNameFile)
        
        # Fazer em um write só
        report = open("Report.txt", 'w')
        report.write(
        f"Date and Hour of execution: {InformationsToFinalReport['initialDateAndHour']}\n" +
        f"Input mutation file are .maf? {InformationsToFinalReport['mutationFileWasConverted']}\n" +
        f"Lenght of network: {str(InformationsToFinalReport['lenghtOfNetwork'])}\n" +
        f"Number of mutated genes: {str(InformationsToFinalReport['numberOfMutatedGenes'])}\n" +
        f"Number of patients: {str(InformationsToFinalReport['numberOfPatients'])}\n" +
        f"Genes not found in the network and excluded: {str(InformationsToFinalReport['numberOfMutatedGenesNotInGraph'])} \n" +
        f"Constant of normalization: {str(InformationsToFinalReport['normalizationConstant'])}\n" +
        f"Start alpha: {str(InformationsToFinalReport['startAlpha'])}\n" +
        f"Best alpha: {str(InformationsToFinalReport['bestAlpha'])}\n"
        )
        print(InformationsToFinalReport)
        report.close()
        
        print(f"Tempo total: {time.time() - initTime}")

###########################################################################################################################