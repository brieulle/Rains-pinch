files = attached_files()

#On va commencer par donner une valeur à p et n, respectivement la caractéristique du
# corps premier et le degré de l'extension qu'on va considérer, pour pouvoir 
#commencer à taper un peu de code

p = 11
n = 23		#On se base sur l'exemple de Pinch, je pense que c'est le plus simple
m = 829
R.<X> = PolynomialRing(GF(p))

def	pinch_method(p, n, m, f = zero, g = zero):		#le m est temporaire

	c, w = cputime(), walltime()
	R.<x> = PolynomialRing(GF(p))
	
	#Pour commencer on met deux polynômes au hasard, on pourra peut-être préciser 
	#deux arguments au cas où on veut donner des polynômes particuliers
	
	if f == zero:
		f = R.irreducible_element(n, algorithm='random')
	if g == zero:
		g = R.irreducible_element(n, algorithm='random')
#f = x**23 + 8*x**2 + x + 9
#g = y**23 + 3*y**2 + 4*y + 9

#On ne veut pas avoir la même extension (il suffirait de prendre l'identité comme isomorphisme):
	while f == g:
		g = R.irreductible_element(n, algorithm='random')


#On créé les deux extensions
	F.<a> = GF(p^n, modulus = f)
	G.<b> = GF(p^n, modulus = g)

#TODO : Il faut trouver un moyen de chercher un m "petit" qui diviser l'ordre de F* 
#(et donc de G*)  pour appliquer la méthode cyclotomique : Appliquer une 
#factorisation ? Chercher le modulo de l'ordre jusqu'à un certain nombre ?
# Ou bien carrément essayer de diviser l'ordre par tous les nombres en-dessous d'une 
#certaine borne ? 


	cofact = F.cardinality() // m	#On peut ne le faire qu'une fois puisque F & G 
                                    #ont le même cardinal

	
	#Une fois qu'on a m, on prend un élément au hasard dans F* et on l'élève 
	#à la puissance (F.order() - 1)/m en espérant tomber sur une racine primitive qui 
	#en plus engendre F (En gros, le m doit diviser qu'une seule fois l'ordre du groupe, 
	#si j'ai bien compris)

	fact = m.factor()
	rootmf = 1 			#Recherche d'une racine primitive de l'unité dans F
	while any(rootmf**(m // k[0]) == 1 for k in fact):
			rootmf = F.random_element()**cofact

	rootmg = 1 	    	#Recherche d'une racine primitive de l'unité dans G
	while any(rootmg**(m // k[0]) == 1 for k in fact):
			rootmg = G.random_element()**cofact


#On défini les matrices qui contiendront les coefficiens qui nous intéressent, 
#c'est le début de la partie "algèbre linéaire"
	A = matrix(GF(p), n, n)
	B = matrix(GF(p), n, n)

	for i in range(n):  #Calcul de la matrice A qui contient les coefficients 
                            #de chaque élément de la base par rapport à l'élement 
                            #qui engendre l'extension (la classe de X)

		temp = rootmf**(p**i)	#En fait, c'est une base si alpha est normal,
                                # ce dont on est pas assuré 
                                #(encore, il s'agira de l'implémentation de Rains)
		A[i,:] = temp.vector()

	try:
		Ainv = A.inverse()  	#Ça ne fonctionnera que si alpha est 
                                #effectivement normal
	except ZeroDivisionError:
		print 'erreur'

#La puissance est celle qui sert à calculer le beta pour trouver 
#effectivement l'isomorphisme
	puissance = 1	

#Dans la méthode de Pinch, il est précisé que si alpha est une puissance primitive 
#m-ième de l'unité alors son image par l'isomorphisme est égal à une puissance de 
#beta; si c'est effectivement un isomorphisme, c'est ce que cette boucle essaie de 
#vérifier

	while puissance <= m :	
		for i in range(n):
			temp = rootmg**(puissance*p**i)
			B[i,:] = temp.vector()


		C = Ainv*B	

		v = C[1,:]	#On prend la deuxième ligne qui correspond à l'image de x
		res = 0

		bs = [1]
		for i in range(n-1):
			bs.append(b*bs[-1])

		for k in range(n):				
			res = res + v[0,k]*bs[k]

		if f(res) == 0:
			print 'CPU %s, Wall %s' % (cputime(c), walltime(w))	
			return [res,C,rootmf,rootmg, puissance,f, F, G, bs]

		puissance = puissance + 1
	
	print 'Pas trouver isomorphisme'
	
#TODO : ptit rappel : trouver une façon plus pratiquer d'exprimer un vecteur 
#sous la forme une somme de puissance de a

def calcul_img(mat, elem, F, G, bs = zero):   #le but c'est de calculer phi(elem)
    c, w = cputime(), walltime()
    R.<x> = PolynomialRing(GF(p))

    n = F.degree()

    if bs == zero:
        bs = [1]
        for i in range(n-1):
            bs.append(G.gen()*bs[i])


#    bs = [1]
#    for i in range(n):
#		bs.append(b*bs[-1])


    tempvec = elem.vector()
    res = 0
    
    for i in range(n):
        if tempvec[i] != 0:
            v = mat[i,:]
            
            for j in range(n):
                res = res + v[0,j]*bs[j]

    print 'CPU %s, Wall %s' % (cputime(c), walltime(w))
    return res

    
    
    
    
    
