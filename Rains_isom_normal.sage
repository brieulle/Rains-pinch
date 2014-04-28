files = attached_files()



def isom_normal(v, w, F, G, base_normale_w = zero, base_normale_v = zero):
    '''
    On prend deux éléments normaux qui correspondent via un isomorphisme phi.
    Le but est de récupérer cet isomorphisme.
    
    Concrètement, l'isomorphisme étant entièrement déterminer par l'image de x,
    il suffit d'exprimer x en fonction de la base normale v^{p^i} puis en
    fonction de la base normale w^{p^i}, en tenant compte du fait que :
    
    x = sum_i c_i*v^{p^i} => phi(x) = sum_i c_i*w^{p^i}
    '''

    p = F.characteristic()
    n = F.degree()

    #Cette boucle serait éventuellement à améliorer et si possible
    #les récupérer de précédentes fonctions ou alors les renvoyer pour
    #éviter de les recalculer à chaque fois

    if base_normale_w == zero:
        base_normale_w = [w]
        for i in range(n):
            base_normale_w.append(base_normale_w[-1]**p)

    temp_normal = convert(F.gen(), v, F)[0] #On récupère les coeff de x
                                            #en fonction de la base normale
                                            #définie par v

    return sum([temp_normal[i]*base_normale_w[i]    #On renvoie l'image de x
                for i in range(len(temp_normal))])
                
#v élément normale de F, z l'élément à exprimer 
#par rapport à la base normale
def convert(z, v, F, base_normale = zero):
    '''
    Le but va être d'implémenter la "conversion" en base normale, i.e. 
    le passage de la base x**i à la base v**p**i pour v un élément normale
    '''
    n = F.degree()
    p = F.characteristic()

    #L'anneau cyclotomique isomorphe aux matrices circulantes n x n
    temp.<U> = PolynomialRing(GF(p))
    Cycl.<u> = temp.quo(U**n - 1) 

    if base_normale == zero:
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
#On pourra y appliquer un calcul de pgcd
        
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

        c.append(sum([inv_list[i]*var_trz[i] for i in range(n)]))

#On passe à la ligne suivante en utilisant la structure des matrices circulante, 
#i.e. b_ij = b_i+1j+1. Ou pour résumer encore mieux, on applique une permutation 
#vers la droite
        temp_coeff = inv_list.pop(-1)       
        inv_list.insert(0, temp_coeff)


    return [c, base_normale, B, inv_list]                
                
def calcul_isom_normal(elem, F, G, img_x):
    '''
    Fonction qui prend un élément elem de F et exprime son image dans G par 
    l'isomorphisme défini par l'image du générateur x, img_x := phi(x).
    '''
    
    n = F.degree()
    
    elem_vector = elem.vector()
    
    puis_img = [1]
    
    for i in range(n):
        puis_img.append(puis_img[-1]*img_x)
        
    return sum([elem_vector[i]*puis_img[i] for i in range(n)])
