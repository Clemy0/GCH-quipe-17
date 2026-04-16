import numpy as np
import matplotlib.pyplot as plt
from fct_FusionMateriau import *

class parametres():
    L = 1                                           # Longueur
    T_l = 0.000000000001                                      # Température du liquidus(état liquide)
    T_s = 0                                         # Température de solidus(état solide)
    T_c = 1                                         # Température imposée
    k = 10                                          # conductivité thermique
    delta_H = 80                                 # Différence d'enthalpie entre les états solide et liquide
    rho = 1000                                      # Masse volumique
    Beta = [0,1]                                    # Paramètre [beta_min, beta_max]
    C_pl = 1                                        # Capacité calorifique à phase liquide
    N = 50                                         # Nombre de noeuds
    t = 50                                        # Temps total de la simulation

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
beta_0 = bissection(prm.Beta[0],prm.Beta[1],prm)
print(beta_0)
t,x,T_Euler_explicite = Euler_explicite(prm)
t_implicite,x_implicite,T_implicite = Euler_implicite(prm, fake_stability_value=0.5) 

t_instable, x_instable, T_Euler_explicite_instable =  Euler_explicite(prm, 0.6) 
T_euler_explicite_Ste01 = Euler_explicite(prm1)[2]
T_euler_explicite_Ste001 = Euler_explicite(prm2)[2]
T_euler_explicite_Ste0001 = Euler_explicite(prm3)[2]

T_theo = func_T(x,t,prm)

x_front_analytique, x_front_numérique = front_evolution(x, t, T_Euler_explicite, prm)

#Graphique pour déterminer l'intervalle de valeur initiale pour trouver beta avec la méthode bissection
plt.figure(1)
plt.plot(beta, f, label= r'f($\beta$)')
plt.plot(prm.Beta[0], funcbeta(prm.Beta[0],prm), linestyle = "None", marker = "o", color = "r", label = r"$\beta_1$")
plt.text(prm.Beta[0], funcbeta(prm.Beta[0],prm), r"$\beta_1$", va = "bottom")
plt.text(prm.Beta[0], funcbeta(prm.Beta[0],prm), rf"f($\beta_1$) = {funcbeta(prm.Beta[0],prm) :.2f}", ha = "right",  va = "top")

plt.plot(beta_0, funcbeta(beta_0,prm), linestyle = "None", marker = "o", color = "black", label = r"$\beta_0$")
plt.text(beta_0, funcbeta(beta_0,prm), r"$\beta_0$", va = "bottom")

plt.plot(prm.Beta[1], funcbeta(prm.Beta[1],prm), linestyle = "None", marker = "o", color = "b", label = r"$\beta_2$")
plt.text(prm.Beta[1], funcbeta(prm.Beta[1],prm), r"$\beta_2$", va = "bottom")
plt.text(prm.Beta[1], funcbeta(prm.Beta[1],prm), rf"f($\beta_2$) = {funcbeta(prm.Beta[1],prm) :.2f}", ha = "right",  va = "top")


plt.plot(np.linspace(0,1,100),np.zeros(100), color = "black", linestyle = "dotted")
plt.xlabel(r'$\beta$')
plt.ylabel(r'f($\beta$)')

plt.xlim(left = -0.25)
plt.grid()
plt.legend()

#Comparaison analytique + Numérique, il manque Euler implicite
plt.figure(2)
plt.plot(x, T_Euler_explicite[len(t)-1,:], label='T(x) explicite à t=50s')
plt.plot(x_implicite, T_implicite[len(t_implicite)-1,:], label='T(x) implicite à t=50s')
plt.plot(x, T_theo[len(t)-1,:], label='T(x) théorique à t=50s')
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
#Stabilité du schéma d'Euler explicite, on voit que pour S = 0.6 > 0.5, la température devient instable à t=50s alors que pour S = 0.5, la température est stable à t=50s
plt.figure(5)
plt.plot(x, T_theo[len(t)-1,:], label='T(x) théorique à t=50s')
plt.plot(x, T_Euler_explicite[len(t)-1,:], label='S = 0.5, T(x) stable à t=50s')
plt.plot(x_instable, T_Euler_explicite_instable[len(t_instable)-1,:], label='S = 0.6 > 0.5, T(x) à t=50s')
plt.xlabel('x')
plt.ylabel('T')
plt.title('Température en fonction de la position pour différente valeur de stabilité')
plt.grid()
plt.legend()
#plt.show()

''' Plot pour les résidus '''
# Solution explicite
t_exp, x_exp, T_exp = Euler_explicite(prm)

# Solution implicite
t_imp, x_imp, T_imp = Euler_implicite(prm, fake_stability_value=0.5)

# Résidus
R_exp, res_exp = residus_explicite(t_exp, x_exp, T_exp, prm)
R_imp, res_imp = residus_implicite(t_imp, x_imp, T_imp, prm)

#résidu explicite
plt.figure(6)
plt.semilogy(t_exp[1:], res_exp, label="Euler explicite")
plt.xlabel("temps")
plt.ylabel("résidu")
plt.title("Résidu du schéma explicite")
plt.grid(True)
plt.legend()

#résidu implicite
plt.figure(7)
plt.semilogy(t_imp[1:], res_imp, label="Euler implicite")
plt.xlabel("temps")
plt.ylabel("résidu")
plt.title("Résidu du schéma implicite")
plt.grid(True)
plt.legend()

# Comparaison
plt.figure(8)
plt.semilogy(t_exp[1:], res_exp, label="Explicite")
plt.semilogy(t_imp[1:], res_imp, label="Implicite")
plt.xlabel("temps")
plt.ylabel("résidu")
plt.title("Comparaison des résidus")
plt.grid(True)
plt.legend()

plt.show()
