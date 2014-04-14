files = attached_files()

#Le but va être d'implémenter la "conversion" en base normale, i.e. le passage de la base x**i à la base v**p**i pour v un élément normale

def convert(z, v, F): #v élément normale de F, z l'élément à exprimer par rapport à la base normale
    x = F.gen()
    n = F.degree()
    p = F.characteristic()

    temp.<U> = PolynomialRing(GF(p))
    Cycl.<u> = temp.quo(U**n - 1)           #L'anneau cyclotomique qui est isomorphe à l'ensemble des matrices circulantes de tailles nxn

    base_normale = [v]
    for i in range(n):
        base_normale.append(base_normale[-1]**p)


    b = []
    for i in range(n):             #On va calculer uniquement Tr(v*v^(p^(n-i))) puisque la matrice est circulante
        b.append((v*base_normale[n-i]).trace())


#calcul de la base de l'anneau cyclotomique, TODO : La calculer de façon plus optimisée (machine rapide avec la trace)
    base_cycl = [1]
    for i in range(n):
        base_cycl.append(u*base_cycl[-1])
        
#On va inverser alors l'élément  de l'anneau correspondant à la matrice
    temp_elem = 0
    for i in range(n):
        temp_elem = temp_elem + b[i]*base_cycl[i]
        
#TODO : Pour que A*B = Id, on a via l'isomorphisme et les propriétés sur les matrices circulantes que c'est équivalent à : a_00.b_00 + a_02.b_01 + a_01.b_02 = 1 ou encore trouver les b_0i qui ont un pgcd égal à 1 avec comme coefficient de bezout les a_0i; alors un algorithme d'euclide "inverse" est envisageable. Mais pour le moment on va pas se fatiguer :
        
    temp = temp_elem**(-1)
    inv_list = temp.list()  

    val_trz = []                    #Liste contenant les valeurs de Tr(v.z^(p^(n-1)))

    for i in range(n):
        temp = (v*(z**p**(n-i)))
        val_trz.append(temp.trace())
    
    c = []

    #On écrite l'élément en fonction de la base normale
    for j in range(n):                

        temp = 0

        #On calcul les c_i tels que z = sum_ i c_i*v**p**i:
        for i in range(n):
            temp =  temp + inv_list[i]*val_trz[i]

        c.append(temp)

#On passe à la ligne suivante en utilisant la structure des matrices circulante, i.e. b_ij = b_i+1j+1. Ou pour résumer encore mieux, on applique une permutation vers la droite
        temp_coeff = inv_list.pop(-1)       
        inv_list.insert(0, temp_coeff)


    return [c, base_normale]
        
