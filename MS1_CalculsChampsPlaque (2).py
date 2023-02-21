import numpy as np
import matplotlib.pyplot as plt

"""
    On va travailler en mètres, pour plus de facilité
"""
#########################################################
########### Définit l'émeteur et le récepteur ###########
#########################################################
E = np.array([0.05,0,0])
Rec = np.array([0,0,0])

#########################################################
################ On positionne la plaque ################
#########################################################
Pc = np.array([0.05, 10, 0]) # Centre de la plaque
theta = -(15/180)*np.pi # Angle de la plaque par rapport à l'axe x
a = 1 # Dimension de la plaque carrée [m]

#########################################################
########### On crée un maillage de la plaque ############
#########################################################
"""
    Rappel : on doit faire un maillage avec des points séparés entre eux d'environ du 10ème
        de la longueur d'onde (rappel, on travaille à 2.4GHz)
"""
#########################################################
############## On définit les constantes ################
#########################################################
h = 0.125
pas = 0.0125 # pas pour faire le maillage, donc en gros on aura 80 points
k = (2*np.pi*2.4*1e9)/(3*1e8)
I = np.array([0,0,1]) # courant dans l'émetteur
const1 = 1j * k * (376.730 / (4*np.pi))
n = np.array([np.sin(theta), np.cos(theta), 0]) # vecteur normal à la plaque
xshift = 0.05
yshift = 10

#########################################################
######## On définit les coordonnées de la plaque ########
#########################################################
x = np.linspace(-0.5+pas/2, 0.5-pas/2, 80)
z = np.linspace(-0.5+pas/2, 0.5-pas/2, 80)
y = np.linspace(0,0,80)

"""
Rappel de la matrice de rotation:
    R = np.matrix([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    Pour plus de facilité, on va d'abord centrer la plaque, lui faire subir une rotation, puis la faire subir une translation
    On stockera alors dans x et y les coordonnées des points sur lesquels nous calculerons des valeurs.
"""
x = np.cos(theta)*x[:] + np.linspace(xshift, xshift, x.size)
y = np.sin(theta)*x[:] + np.linspace(yshift,yshift,x.size)


#########################################################
########### On crée un maillage sur la plaque ###########
#########################################################

matrix = []
for i in range(x.size):
    a = []
    for j in range(x.size):
        a.append([x[i], y[i], z[j]])
    matrix.append(a)
#print(matrix)
# Dans matrix se trouve tous les points sur lesquels il faut intégrer.
#   Par exemple, matrix[0][0] est le point inférieur gauche de la plaque.

#########################################################
############ Représentation de la situation #############
#########################################################

plt.figure(figsize = (16,10))
plt.plot(x,y)
plt.scatter([0.05], [10], color = "magenta") # Centre de la plaque
plt.scatter([0], [0], color = "green") # Récepteur
plt.scatter([0.05], [0], color = "red") #Emetteur
plt.xlim((-1,1))
plt.ylim((0, 10.5))
plt.xlabel("Axe x - [m]")
plt.ylabel("Axe y - [m]")
plt.title("Situation avec l'émetteur en (0.05, 0) et le centre de la plaque en (0.05,10)")
plt.savefig("Situation.pdf")
plt.show()

#########################################################
######################### Code ##########################
#########################################################
E_p = []
H_p = []
J_p = []

for i in range(x.size):
    E_pi = []
    H_pi = []
    J_pi = []
    for j in range(x.size):
        R = np.linalg.norm((matrix[i][j] - E))
        u = (matrix[i][j] - E) / R
        # Rappel, impédance du vide ≃ 376.730
        a = const1 * np.cross(u, np.cross(u,I)) * np.exp(-1j*k*R) / R
        b = np.cross(u, a)/376.730
        
        E_pi.append(a)
        H_pi.append(b)
        J_pi.append(2*np.cross(n,b))
    E_p.append(E_pi)
    H_p.append(H_pi)
    J_p.append(J_pi)
#print(np.abs(J_p))
# Rappel: en chaque point de la plaque, E est un vecteur. Par conséquent, notre vecteur E a cette forme : 81x81x3 (oui c'est énorme)


# X, Y = np.meshgrid(x, y)
# plt.pcolor(X, Y, )
# plt.show()


def champ_récepteur(rec):
    E = np.array([0, 0, 0])
    for i in range(x.size):
        for j in range(x.size):
            R = np.linalg.norm((matrix[i][j] - rec))
            u = (matrix[i][j] - rec) / R
            a = const1 * np.cross(u, np.cross(u,J_p[i][j])) * np.exp(-1j*k*R) / R
            E = E + a*pas*pas
    return E
E_rec = champ_récepteur(Rec)
