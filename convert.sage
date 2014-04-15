files = attached_files()

'''
Le but va être d'implémenter la "conversion" en base normale, i.e. 
le passage de la base x**i à la base v**p**i pour v un élément normale
'''

#v élément normale de F, z l'élément à exprimer 
#par rapport à la base normale
def convert(z, v, F):   
    x = F.gen()
    n = F.degree()
    p = F.characteristic()

    #L'anneau cyclotomique isomorphe aux matrices circulantes n x n
    temp.<U> = PolynomialRing(GF(p))
    Cycl.<u> = temp.quo(U**n - 1) 

    base_normale = [v]
    for i in range(n-1):
        base_normale.append(base_normale[-1]**p)


    B = []
    
#Calculer uniquement Tr(v*v^(p^(n-i))) puisque la matrice est circulante
    for i in range(n):    
        B.append((v*base_normale[-i]).trace())


#calcul de la base de l'anneau cyclotomique
#TODO : La calculer de façon plus optimisée
    base_cycl = [1]
    for i in range(n):
        base_cycl.append(u*base_cycl[-1])
        
#On va inverser alors l'élément de l'anneau correspondant à la matrice
    temp_elem = 0
    for i in range(n):
        temp_elem = temp_elem + B[i]*base_cycl[i]
        
#TODO : Pour que A*B = Id, on a via l'isomorphisme et les propriétés 
#sur les matrices circulantes que c'est équivalent à (pour une matrice 3x3):
#a_00.b_00 + a_02.b_01 + a_01.b_02 = 1
#On pourra y appliquer un calcule de pgcd
        
    temp = temp_elem**(-1)
    inv_list = temp.list()  

    val_trz = []        #Liste contenant les valeurs de Tr(v.z^(p^(n-1)))

    for i in range(n):
        temp = (v*(z**p**(n-i)))
        val_trz.append(temp.trace())
    
    c = []

    #On écrit l'élément en fonction de la base normale
    for j in range(n):                

        temp = 0

        #On calcul les c_i tels que z = sum_ i c_i*v**p**i:
        for i in range(n):
            temp =  temp + inv_list[i]*val_trz[i]

        c.append(temp)

#On passe à la ligne suivante en utilisant la structure des matrices circulante, 
#i.e. b_ij = b_i+1j+1. Ou pour résumer encore mieux, on applique une permutation 
#vers la droite
        temp_coeff = inv_list.pop(-1)       
        inv_list.insert(0, temp_coeff)


    return [c, base_normale, B, inv_list]
        
        
files = attached_files()

def isom_normal(v, w, F, G, base_normale_v = zero, base_normale_w = zero):
    p = F.characteristic()
    n = F.degree()

    #Ces deux boucles seraient éventuellement à améliorer et si possible
    #les récupérer de précédentes fonctions ou alors les renvoyer pour
    #éviter de les recalculer à chaque fois
    if base_normale_v == zero:
        base_normale_v = [v]
        for i in range(n):
            base_normale_v.append(base_normale_v[-1]**p)

    if base_normale_w == zero:
        base_normale_w = [w]
        for i in range(n):
            base_normale_w.append(base_normale_w[-1]**p)

    A = matrix(GF(p), n, n)
    B = matrix(GF(p), n, n)

    for i in range(n):
        A[i,:] = base_normale_v[i].vector() 

    #Il y a probablement mieux à faire qu'inverser la matrice de passage;
    #peut-être en s'inspirant de ce qui a été fait au-dessus
    try:
        Ainv = A.inverse()
    except ZeroDivisionError:
        print 'erreur'
        return A

    for i in range(n):
        B[i,:] = base_normale_w[i].vector()


    return Ainv*B

