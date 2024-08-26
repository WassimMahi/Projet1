#Mahieddine Wassim 20142332

# Ce programme permet de partir d'un brin d'ADN pour arriver au protéines
# codées par les gènes contenu dans le brin d'ADN. Cela se fait en créant le 
#brin complémentaire du brin d'ADN fourni. Puis en transcrivant chacun des 
#gènes en ARN. Puis en traduisant l'ARN en proteine à l'aide du tableau 
#associatif et afficher les proteines sous forme de chaine de caractères et en 
#les dessiner  l'aide de la tortue

adn="TCGACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGACG\
GCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGAGTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGAATGCCAGCCAGC\
CAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGAACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGA\
ACTCGACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCATCCCAGCGATACCC\
AGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAG\
CCAGCGAACTCGTCTGCGTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGC\
GATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGTATGCCAGCC\
AGCATCCCAGCGA"

codons_aa = {
    "UUU": "Phénylalanine",
    "UUC": "Phénylalanine",
    "UUA": "Leucine",
    "UUG": "Leucine",
    "CUU": "Leucine",
    "CUC": "Leucine",
    "CUA": "Leucine",
    "CUG": "Leucine",
    "AUU": "Isoleucine",
    "AUC": "Isoleucine",
    "AUA": "Isoleucine",
    "AUG": "Méthionine (Start)",
    "GUU": "Valine",
    "GUC": "Valine",
    "GUA": "Valine",
    "GUG": "Valine",
    "UCU": "Sérine",
    "UCC": "Sérine",
    "UCA": "Sérine",
    "UCG": "Sérine",
    "CCU": "Proline",
    "CCC": "Proline",
    "CCA": "Proline",
    "CCG": "Proline",
    "ACU": "Thrénine",
    "ACC": "Thrénine",
    "ACA": "Thrénine",
    "ACG": "Thrénine",
    "GCU": "Alanine",
    "GCC": "Alanine",
    "GCA": "Alanine",
    "GCG": "Alanine",
    "UAU": "Tyrosine",
    "UAC": "Tyrosine",
    "UAA": "Stop",
    "UAG": "Stop",
    "CAU": "Histidine",
    "CAC": "Histidine",
    "CAA": "Glutamine",
    "CAG": "Glutamine",
    "AAU": "Asparagine",
    "AAC": "Asparagine",
    "AAA": "Lysine",
    "AAG": "Lysine",
    "GAU": "Aspartate",
    "GAC": "Aspartate",
    "GAA": "Glutamate",
    "GAG": "Glutamate",
    "UGU": "Cystéine",
    "UGC": "Cystéine",
    "UGA": "Stop",
    "UGG": "Tryptophane",
    "CGU": "Arginine",
    "CGC": "Arginine",
    "CGA": "Arginine",
    "CGG": "Arginine",
    "AGU": "Sérine",
    "AGC": "Sérine",
    "AGA": "Arginine",
    "AGG": "Arginine",
    "GGU": "Glycine",
    "GGC": "Glycine",
    "GGA": "Glycine",
    "GGG": "Glycine"}

lettreAa = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "UGU": "C",
    "UGC": "C",
    "UGA": "*",
    "UGG": "W",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"
}

def antisens(brinAdn): # Fonction pour obtenir le brin d'ADN complémentaire
    adnFourni = list(brinAdn)
    adnComp = []
    for letter in adnFourni: # On parcourt chaque lettre du brin d'ADN
        if letter == "A":
            adnComp.append("T")
        elif letter == "T":
            adnComp.append("A")
        elif letter == "C":
            adnComp.append("G")
        elif letter == "G":
            adnComp.append("C")
    adnComp.reverse() # On inverse la liste pour obtenir le brin complémentaire
    return adnComp

def trouveDebut(brinAdn): # Fonction pour trouver les codons de début
    adnFourni = list(brinAdn)
    codonStartFourni = []
    for i in range(len(adnFourni) - 2):
        if adnFourni[i] == "T" and adnFourni[i+1] == "A" and adnFourni[i+2] == "C":
            codonStartFourni.append(i)

    adnComp = antisens(brinAdn)
    codonStartComp = []
    for i in range(len(adnComp) - 2):
        if adnComp[i] == "T" and adnComp[i+1] == "A" and adnComp[i+2] == "C":
            codonStartComp.append(i)

    return codonStartFourni, codonStartComp

def trouveFin(brinAdn): # Fonction pour trouver les codons de fin
    adnFourni = list(brinAdn)
    codonStopFourni = []
    for i in range(len(adnFourni) - 2):
        if adnFourni[i:i+3] == ["A", "T", "T"] or adnFourni[i:i+3] == ["A", "T", "C"] or adnFourni[i:i+3] == ["A", "C", "T"]:
            codonStopFourni.append(i)

    adnComp = antisens(brinAdn) 
    codonStopComp = []
    for i in range(len(adnComp) - 2):
        if adnComp[i:i+3] == ["A", "T", "T"] or adnComp[i:i+3] == ["A", "T", "C"] or adnComp[i:i+3] == ["A", "C", "T"]:
            codonStopComp.append(i)

    return codonStopFourni, codonStopComp

def trouveGene(debut, fin): # Fonction pour trouver les gènes
    genes = []
    for debutGene in debut:
        for finGene in fin:
            if finGene > debutGene and (finGene - debutGene) % 3 == 0:
                genes.append((debutGene, finGene))
                break  # Sortir de la boucle intérieure dès qu'un gène est trouvé
    return genes

def transcrire(brinAdn): # Fonction pour transcrire le brin d'ADN en ARN
    return brinAdn.replace("T", "U")

global lettres
lettres=[]


def carre(longueur, nombre):#fct pour dessiner les petits carres
    axeY=0
    pu()
    goto(-150,0)
    pd()
    for i in range(1,nombre+1):
        for _ in range(4) :
            pd()
            fd(longueur)
            rt(90)
        if i%15==0 and i!=0 :#saut a la ligne apres 15 carres
            axeY+=1
            pu()
            goto(-150,-axeY*longueur)
            pd()
        else :   
            fd(longueur)
    #deuxieme boucle qui va mettre les lettres dans chaque carre
    axeY=1 
    pu()
    goto(-150,-10)
    fd(longueur/2)
    pd()
    for i in range(len(lettres)):#on parcours lettres[]  

        if i%15==0 and i!=0 :#idem
            write(lettres[i])
            axeY+=2
            pu()
            goto(-150,-10*axeY)
            fd(longueur/2)
            pd()
        else :  
            
            write(lettres[i])
            pu()
            fd(longueur)
            pd()

            
def traduire(brinArn):  
    # Traduit une séquence d'ARN en acides aminés
    amino_acids = []
    
    # Trouver le premier codon de début AUG
    start = brinArn.find("AUG")
    if start != -1:  # Si un codon de début est trouvé
        # Commencer la traduction à partir du codon de début
        for i in range(start, len(brinArn), 3):
            codon = brinArn[i:i+3]
            amino_acid = codons_aa.get(codon)
            lettres1=lettreAa.get(codon)
            if amino_acid == "Stop":  # Arrêtez la traduction au codon stop
                break
            elif amino_acid:  # Si le codon est connu et non un codon stop
                # Ajoute le nom complet de l'acide aminé
                amino_acids.append(amino_acid)
                lettres.append(lettres1)
    # Joindre les acides aminés avec des tirets et supprimer le "(Start)" pour l'affichage
    protein_sequence = "-".join(amino_acids).replace(" (Start)", "")
    print(lettres)
    print(protein_sequence)
    carre(20,len(lettres))
    return 

   

def fonction_mere(brinAdn): #fonction finale qui orchestre toutes les fcts.
    debutsFourni, debutsComp = trouveDebut(brinAdn)
    finsFourni, finsComp = trouveFin(brinAdn)
    genesFourni = trouveGene(debutsFourni, finsFourni)
    genesComp = trouveGene(debutsComp, finsComp)

    for debut, fin in genesFourni:
        arn = transcrire(brinAdn[debut:fin+3])
        traduire(arn)

    for debut, fin in genesComp:
        arn_comp = transcrire(brinAdn[debut:fin+3])
        traduire(arn_comp)


fonction_mere(adn)



# test unitaires
def test_antisens():
    # Test avec brin ADN simple
    assert antisens("CAGT") == ["A", "C", "T", "G"]

    #Test avec brin ADN qui contient des caractères spéciaux
    assert antisens("CAGT") == ["A", "C", "T", "G","!"]
    
    # Test avec brin ADN vide
    assert antisens("") == []

def test_trouveDebut():
    # Test avec un brin ADN qui contient un seul codon de départ
    assert trouveDebut("TACCATG") == ([1], [4])

    # Test avec un brin ADN qui contient plusieurs codons de départ
    assert trouveDebut("TACCATGACCATG") == ([1, 8], [4, 11])

    # Test avec brin ADN qui ne contient pas de codons de départ
    assert trouveDebut("ACCGTAA") == ([], [])

def test_trouveFin():
    # Test avec brin ADN qui contient un seul codon de terminaison
    assert trouveFin("TACCATG") == ([5], [2])

    
    # Test avec brin ADN qui contient plusieurs codons de terminaison
    assert trouveFin("TACCATGACCATG") == ([5, 12], [2, 9])


    # Test avec un brin ADN ne contient  pas de codons de terminaison
    assert trouveFin("ACCGTAA") == ([], [])

def test_trouveGene():
   # Test avec un brin ADN qui contient un seul gène
    assert trouveGene([1], [5]) == [(1, 5)]

    # Test avec brin ADN qui contient  plusieurs gènes
    assert trouveGene([1, 8], [5, 12]) == [(1, 5), (8, 12)]

    #Test d'un brin ADN sans genes
    assert trouveGene([], []) == []

def test_transcrire():
    # Test avec un brin ADN simple
    assert transcrire("CAGT") == "CAGU"

    # Test avec un brin ADN qui contient un caractère spécial
    assert transcrire("CAGT!") == "CAGU!"

    # Test avec un brin ADN vide
    assert transcrire("") == ""

