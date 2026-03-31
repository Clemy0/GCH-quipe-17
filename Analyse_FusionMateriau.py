import numpy as np
import matplotlib.pyplot as plt
from fct_FusionMateriau import *

class parametres():
    L = 1                                           # Longueur
    T_l = 0.000000001                                      # Température du liquidus(état liquide)
    T_s = 0                                         # Température de solidus(état solide)
    T_c = 1                                         # Température imposée
    k = 10                                          # conductivité thermique
    delta_H = 10                                 # Différence d'enthalpie entre les états solide et liquide
    rho = 1000                                      # Masse volumique
    Beta = [0,1]                                    # Paramètre [beta_min, beta_max]
    C_pl = 1                                        # Capacité calorifique à phase liquide
    N = 100                                         # Nombre de noeuds
    t = 465                                        # Temps total de la simulation

class parametres1(parametres): #Ste = 0.1
    delta_H = 10
    
class parametres2(parametres): #Ste = 0.01
    delta_H = 100

class parametre3(parametres):  #Ste = 0.001
    delta_H = 1000



prm  = parametres()
prm1 = parametres1()
prm2 = parametres2()
prm3 = parametre3()

beta, f = funcbetaArray(prm)



t,x,T_Euler_explicite = Euler_explicite(prm)
T_euler_explicite_Ste01 = Euler_explicite(prm1)[2]
T_euler_explicite_Ste001 = Euler_explicite(prm2)[2]
T_euler_explicite_Ste0001 = Euler_explicite(prm3)[2]

T_theo = func_T(x,t,prm)

x_front_analytique, x_front_numérique = front_evolution(x, t, T_Euler_explicite, prm)

#Graphique pour déterminer l'intervalle de valeur initiale pour trouver beta avec la méthode bissection
plt.figure(1)
plt.plot(beta, f, label='f(beta)')
plt.xlabel('beta')
plt.ylabel('f(beta)')
plt.title('Fonction f(beta)')
plt.grid()
plt.legend()

#Comparaison analytique + Numérique, il manque Euler implicite
plt.figure(2)
plt.plot(x, T_Euler_explicite[len(t)-1,:], label='T(x) à t=50s')
plt.plot(x, T_theo[len(t)-1,:], label='T(x) à t=50s')
plt.xlabel('x')
plt.ylabel('T')
plt.title('Température en fonction de la position')
plt.grid()
plt.legend()

#Graph du résultat obtenu en fonction de Ste, à voir s'il faut mettre implicite aussi!!??! ­>_<
plt.figure(3)
plt.plot(x, T_euler_explicite_Ste01[len(t)-1,:], label='T(x) à t=50s et Ste = 0.1')
plt.plot(x, T_euler_explicite_Ste001[len(t)-1,:], label='T(x) à t=50s et Ste = 0.01')
plt.plot(x, T_euler_explicite_Ste0001[len(t)-1,:], label='T(x) à t=50s et Ste = 0.001')
plt.xlabel('x')
plt.ylabel('T')
plt.title('Température en fonction de la position')
plt.grid()
plt.legend()

#Comparaison du front de fusion trouver analytiquement et numériquement
plt.figure(4)
plt.plot(t, x_front_analytique, label='Front de fusion analytique en fonction du temps')
plt.plot(t, x_front_numérique, label='Front de fusion numérique en fonction du temps')
plt.xlabel('t')
plt.ylabel('x_front')
plt.title('Évolution du front de fusion')
plt.grid()
plt.legend()
plt.show()
