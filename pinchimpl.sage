files = attached_files()
#Il faudra pouvoir préciser qu'on utilise ce code si algorithm est égal à "rains" dans la fonction qui va bien (probableme Hom(E,F) ou une autre)

#On va commencer par donner une valeur à p et n, respectivement la caractéristique du corps premier et le degré de l'extension qu'on va considérer, pour pouvoir commencer à taper un peu de code

p = 11
n = 23		#On se base sur l'exemple de Pinch, je pense que c'est le plus simple
m = 829
R.<x> = PolynomialRing(GF(p))

def	pinch_method(p, n, m, f = zero, g = zero):		#le m est temporaire

	R.<x> = PolynomialRing(GF(p))
#Pour commencer on met deux polynômes au hasard, on pourra peut-être préciser deux arguments au cas où on veut donner des polynômes particuliers
	if f == zero:
		f = R.irreducible_element(n, algorithm='random')
	if g == zero:
		g = R.irreducible_element(n, algorithm='random')
#f = x**23 + 8*x**2 + x + 9
#g = y**23 + 3*y**2 + 4*y + 9

	while f == g:		#On ne veut pas avoir la même extension (il suffirait de prendre l'identité comme isomorphisme)
		g = R.irreducible_element(n, algorithm='random')


#On créé les deux extensions
	F.<a> = GF(p^n, modulus = f)
	G.<b> = GF(p^n, modulus = g)

#TODO : Il faut trouver un moyen de chercher un m "petit" qui diviser l'ordre de F* (et donc de G*)  pour appliquer la méthode cyclotomique : Appliquer une factorisation ? Chercher le modulo de l'ordre jusqu'à un certain nombre ? Ou bien carrément essayer de diviser l'ordre par tous les nombres en-dessous d'une certaine borne ? 

#Remarque : Si m est premier alors la chance d'avoir une racine primitive est quasi sûre (= phi(m)/m = (m-1)/m)


#Une fois qu'on a m, on prend un élément au hasard dans F* et on l'élève à la puissance (F.order() - 1)/m en espérant tomber sur une racine primitive qui en plus engendre F (En gros, le m doit diviser qu'une seule fois l'ordre du groupe, si j'ai bien compris)

	temp = 0
	while temp == 0:
		temp = F.random_element()

	rootmf = temp**((F.order() - 1)/m)

#rootmf = a**((F.order() - 1)/m)				

	while cyclotomic_polynomial(m)(rootmf) != 0:				#On vérifie qu'elle est primitive (faute de mieux je garde cette méthode)
		rootmf = temp**((F.order() - 1)/m)

#Pour le moment on applique la méthode Finch, alors on utilise exactement le même procédé dans G

	temp = 0
	while temp == 0:
		temp = G.random_element()

	rootmg = temp**((G.order() - 1)/m)

#rootmg = b**((G.order() - 1)/m)

	while cyclotomic_polynomial(m)(rootmg) != 0:				#On vérifie qu'elle est primitive (faute de mieux je garde cette méthode)
		rootmg = temp**((F.order() - 1)/m)

#On défini les matrices qui contiendront les coefficients qui nous intéressent, c'est le début de la partie "algèbre linéaire"
	A = matrix(GF(p), n, n)
	B = matrix(GF(p), n, n)

	for i in range(n):				#Calcul de la matrice A qui contient les coefficients de chaque élément de la base par rapport à l'élement qui engendre l'extension (la classe de X)
		temp = rootmf**(p**i)
		A[i,:] = temp.vector()


	success, puissance = 0, 1	#success sert à déterminer si on a effectivement trouver une racine m-ième correspondante

	while puissance <= m and success !=1:	#Dans la méthode de Pinch, il est précisé que si alpha est une puissance primitive m-ième de l'unité alors son image par l'isomorphisme est égal à une puissance de beta; si c'est effectivement un isomorphisme, c'est ce que cette boucle essaie de vérifier
		for i in range(n):
			temp = rootmg**(puissance*p**i)
			B[i,:] = temp.vector()


		C = A.inverse()*B	

		for j in range(n):	#Le but de la boucle est de tester chaque ligne de coefficients afin de voir si elle fournit le candidat désiré
			v = C[j,:]		#Vu que la première ligne fourni le coefficient dans F_p directement (pour b^0), je ne sais pas si ça vaut vraiment le coup de le tester 
			res = 0

			for k in range(n):				
				res = res + v[0,k]*b**k

			if f(res) == 0:
				success = 1
				break

		puissance = puissance + 1
	
	return [res, f, g, A, B, C]




