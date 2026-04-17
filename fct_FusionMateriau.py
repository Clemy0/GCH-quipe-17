import numpy as np
from math import erf

def Stefunc(prm):
    """Fonction qui calcule le nombre de Stefan en fonction des paramètres d'entrée

    Entrées:
        - prm : Objet class parametres()
            - L         : Longueur
            - T_l       : Température du liquidus(état liquide)
            - T_s       : Température de fusion(état solide) 
            - T_c       : Température imposée 
            - k         : conductivité thermique 
            - delta_H   : Différence d'enthalpie entre les états solide et liquide
            - rho       : Masse volumique
            - beta      : Paramètre [beta_min, beta_max]
            - C_pl      : Capacité calorifique à phase liquide
            - N         : Nombre de noeuds
            - t         : Temps total de la simulation

    Sorties (dans l'ordre énuméré ci-bas):
        - Ste : Nombre de Stefan (sans unité)
     """
    # Calcul du nombre de Stefan
    Ste = (prm.C_pl *(prm.T_c-prm.T_s)) / prm.delta_H
    return Ste

def C_psl(prm):
    """Fonction qui calcule la capacité calorifique de changement de phase

    Entrées:
        - prm : Objet class parametres()
            - L         : Longueur
            - T_l       : Température du liquidus(état liquide)
            - T_s       : Température de fusion(état solide) 
            - T_c       : Température imposée 
            - k         : conductivité thermique 
            - delta_H   : Différence d'enthalpie entre les états solide et liquide
            - rho       : Masse volumique
            - beta      : Paramètre [beta_min, beta_max]
            - C_pl      : Capacité calorifique à phase liquide
            - N         : Nombre de noeuds
            - t         : Temps total de la simulation

    Sorties (dans l'ordre énuméré ci-bas):
        - C_psl : Capacité calorifique de changement de phase
    """
    # Calcul de la capacité calorifique de changement de phase
    C_psl = (prm.C_pl*(prm.T_l-prm.T_s) + prm.delta_H)/(prm.T_l - prm.T_s)
    return C_psl

def funcbeta(beta, prm):
    """ Fonction qui évalue la fonction étudiée f(beta) qui correspond à un nombre de Stefan donné )
    
    Entrées:
        - beta : Valeur du paramètre beta pour lequel on veut évaluer la fonction
        - prm : Objet class parametres()
            - L         : Longueur
            - T_l       : Température du liquidus(état liquide)
            - T_s       : Température de fusion(état solide)
            - T_c       : Température imposée
            - k         : conductivité thermique
            - delta_H   : Différence d'enthalpie entre les états solide et liquide
            - rho       : Masse volumique
            - Beta      : Paramètre [beta_min, beta_max]
            - C_pl      : Capacité calorifique à phase liquide
            - N         : Nombre de noeuds
            - t         : Temps total de la simulation
    Sorties:
        - f : Valeur de la fonction pour un beta donné
    """
    #calcul du nombre de Stefan
    Ste = Stefunc(prm)

    #calul de la fonction f(beta)
    f = Ste/(np.sqrt(np.pi)) - beta*np.exp(beta**2)*erf(beta)
    return f

def funcbetaArray(prm):
    """Fonction étudiée f(b) qui correspond à un nombre de Stefan donné, pour génération de graphiques
    Entrées:
        - prm : Objet class parametres()
                - L         : Longueur
                - T_l       : Température du liquidus(état liquide)
                - T_s       : Température de fusion(état solide)
                - T_c       : Température imposée
                - k         : conductivité thermique
                - delta_H   : Différence d'enthalpie entre les états solide et liquide
                - rho       : Masse volumique
                - Beta      : Paramètre [beta_min, beta_max]
                - C_pl      : Capacité calorifique à phase liquide
                - N         : Nombre de noeuds
                - t         : Temps total de la simulation
        Sorties:
        -beta   : Array(len = 100) 1D de valeurs de beta entre beta_min et beta_max
        -f      : Array(len = 100) 1D de valeurs de la fonction f(beta)
        """
    #On définit le vecteur beta
    beta = np.linspace(prm.Beta[0], prm.Beta[1], 100)

    #On évalue la fonction f(beta) pour chaque valeur de beta dans le vecteur
    f = np.array([funcbeta(b, prm) for b in beta])
    return beta, f

def bissection(beta1, beta2, prm):
    """Fonction qui utilise la méthode de la bissection pour trouver la valeur de beta qui correspond à un nombre de Stefan donné
    
    Entrées:
        - beta1 : Valeur initiale de beta [borne inférieure de l'intervalle de recherche]
        - beta2 : Valeur initiale de beta [borne supérieure de l'intervalle de recherche]
        - prm : Objet class parametres()
            - L         : Longueur
            - T_l       : Température du liquidus(état liquide)
            - T_s       : Température de fusion(état solide)
            - T_c       : Température imposée
            - k         : conductivité thermique
            - delta_H   : Différence d'enthalpie entre les états solide et liquide
            - rho       : Masse volumique
            - Beta      : Paramètre [beta_min, beta_max]
            - C_pl      : Capacité calorifique à phase liquide             
            - N         : Nombre de noeuds
            - t         : Temps total de la simulation
    Sorties:
        - beta_m : Valeur de beta qui correspond à un nombre de Stefan donné, trouvée par la méthode de bissection
    """
    f1 = funcbeta(beta1, prm)
    error = 1 

    while error > 1e-6:
        beta_m = (beta1 + beta2) / 2
        fm = funcbeta(beta_m, prm)
        if fm * f1 < 0:
            beta2 = beta_m
        else:
            beta1 = beta_m
        error = abs(fm)
    return beta_m

def func_T(t, x, prm):
    """Fonction qui donne la distribution analytique de la température en fonction du temps
    Entrées:
        - t : Array de temps de taille (Nt) 1D 
        - x : Array de position de taille (N) 1D
        - prm : Objet class parametres()
            - L         : Longueur
            - T_l       : Température du liquidus(état liquide)
            - T_s       : Température de fusion(état solide)
            - T_c       : Température imposée
            - k         : conductivité thermique
            - delta_H   : Différence d'enthalpie entre les états solide et liquide
            - rho       : Masse volumique
            - Beta      : Paramètre [beta_min, beta_max]
            - C_pl      : Capacité calorifique à phase liquide               
            - N         : Nombre de noeuds
            - t         : Temps total de la simulation
    Sorties:
        -T  : Array de taille (Nt, N) de la distribution analytique de la température en fonction du temps 
    """
    #Définition des paramètres importants
    Nx = len(x)
    Nt = len(t)

    alpha = prm.k/(prm.rho*prm.C_pl)
    
    T = np.empty((Nt,Nx))
    
    beta = bissection(prm.Beta[0], prm.Beta[1], prm)
    
    x_front_anlytique = func_Evolution_frontAnalytique(t, prm)
    for i in range(Nt):
        for j in range(Nx):
            if i == 0:
                T[i,j] = prm.T_s
            else: 
                if j == 0:
                    T[i, j] = prm.T_c
                elif j == Nx-1:
                    T[i,j] = prm.T_s
                elif x[j] < x_front_anlytique[i]:
                    T[i,j] = (prm.T_s - prm.T_c)*erf(x[j]*beta/(x_front_anlytique[i]))/erf(beta) + prm.T_c 
                else:
                    T[i,j] = prm.T_s
    return T

def func_Evolution_frontAnalytique(t, prm):
    """Fonction qui calcule l'évolution du front de fusion analytique en fonction du temps
    Entrées:
        - t : Array de temps de taille (Nt) 1D
        - prm : Objet class parametres()
                - L         : Longueur
                - T_l       : Température du liquidus(état liquide)
                - T_s       : Température de solidus(état solide)
                - T_c       : Température imposée
                - k         : conductivité thermique
                - delta_H   : Différence d'enthalpie entre les états solide et liquide
                - rho       : Masse volumique
                - Beta      : Paramètre [beta_min, beta_max]
                - C_pl      : Capacité calorifique à phase liquide
                - C_psl     : Capacité calorifique de changement de phase                
                - N         : Nombre de noeuds
                - t         : Temps total de la simulation
    Sorties:
        - x_front_Analytique:  Array de position du front de fusion en fonction du temps, de taille (Nt) 1D
    """    
    #Calcule des paramètres importants 
    Nt = len(t)
    t = np.linspace(0, prm.t, Nt)
    
    beta = bissection(prm.Beta[0], prm.Beta[1], prm)
    alpha = prm.k/(prm.rho*prm.C_pl)
    
    x_front_Analytique = np.empty(Nt)

    #Calcules du front
    for i in range(Nt):
        x_front_Analytique[i] = 2*beta*np.sqrt(alpha*t[i])
    
    return x_front_Analytique
 
def func_Evolution_frontNumerique(t, x, T, prm):
    """Fonction qui calcule l'évolution du front de fusion en fonction du temps
    Entrées:
        - x : Array de position de taille (N) 1D
        - t : Array de temps de taille (N) 1D
        - T : Array de temp de taille (Nt, N) 
        - prm : Objet class parametres()
                - L         : Longueur
                - T_l       : Température du liquidus(état liquide)
                - T_s       : Température de solidus(état solide)
                - T_c       : Température imposée
                - k         : conductivité thermique
                - delta_H   : Différence d'enthalpie entre les états solide et liquide
                - rho       : Masse volumique
                - Beta      : Paramètre [beta_min, beta_max]
                - C_pl      : Capacité calorifique à phase liquide
                - C_psl     : Capacité calorifique de changement de phase                
                - N         : Nombre de noeuds
                - t         : Temps total de la simulation
    Sorties:
        - x_front_Numérique: Array de position du front de fusion en fonction du temps, de taille (Nt) 1D
    """
    Nt = len(t)
    x_front_numerique = np.empty(Nt)

    for i in range(Nt):
        i_front = np.argmax(T[i,:] < prm.T_l) 
        x_front_numerique[i] = x[i_front]
    return x_front_numerique

def Euler_explicite(prm, S = 0.5): 
    """Fonction qui utilise la méthode d'Euler explicite pour résoudre l'équation de la fusion d'un matériau en fonction du temps et de x
    Entrées:
        - prm : Objet class parametres()
                - L         : Longueur
                - T_l       : Température du liquidus(état liquide)
                - T_s       : Température de solidus(état solide)
                - T_c       : Température imposée
                - k         : conductivité thermique
                - delta_H   : Différence d'enthalpie entre les états solide et liquide
                - rho       : Masse volumique
                - Beta      : Paramètre [beta_min, beta_max]
                - C_pl      : Capacité calorifique à phase liquide               
                - N         : Nombre de noeuds
                - t         : Temps total de la simulation
        - S : Coefficient de stabilité du Schéma d'Euler explicite qui sert également à ajuster le pas de temps
    Sorties:
        - t : Array de temps de taille (N) 1D 
        - x : Array de position de taille (N) 1D
        - T : Array de température de taille (N,N) 2D, où T[i,j] est la température au temps t[i] et à la position j
        """
    #Définition des paramètres et des variables
    x = np.linspace(0, prm.L, prm.N)
    dx = x[1]- x[0]
    dt = (S*prm.C_pl*prm.rho*dx**2)/prm.k
    
    Nt = int(prm.t/dt)
    t = np.linspace(0, prm.t, Nt)

    T = np.zeros((Nt, prm.N))
    #Boucle itérative en t et x pour créer la matrice de coefficients associées aux différences finies, i est l'indice de temps et j' est l'indice de position
    for i in range(Nt):
        if i == 0:
            T[i,:] = prm.T_s
        else:
            A = np.zeros((prm.N, prm.N))
            b = np.zeros(prm.N).T
            for j in range(prm.N):
                A[j,j] = 1 
                if j == 0:
                    b[j] = prm.T_c
                elif j == prm.N - 1:
                    b[j] = prm.T_s
                elif T[i-1, j] < prm.T_l: #Si la température à l'itération précédente est inférieure à la température de liquidus, on utilise la capacité calorifique de changement de phase
                    b[j] = (dt*(prm.k/(C_psl(prm)*prm.rho))*(T[i-1, j+1] - 2*T[i-1,j] + T[i-1, j-1])/dx**2) + T[i-1, j]
                else:                     #Sinon, on utilise la capacité calorifique à phase liquide     
                    b[j] = (dt*(prm.k/prm.rho)*(T[i-1, j+1] - 2*T[i-1,j] + T[i-1, j-1])/dx**2)/prm.C_pl + T[i-1, j]            
            T[i,:] = np.linalg.solve(A, b) #On résout le système d'équations pour trouver le profil de température à chaque pas de temps
    return t, x, T

def Euler_implicite(prm, S = 0.5):
    """Fonction qui utilise la méthode d'Euler implicite pour résoudre l'équation de la fusion d'un matériau en fonction du temps et de x.
    Entrées:
        - prm : Objet class parametres()
                - L         : Longueur
                - T_l       : Température du liquidus (état liquide)
                - T_s       : Température du solidus (état solide)
                - T_c       : Température imposée
                - k         : Conductivité thermique
                - delta_H   : Différence d'enthalpie entre solide et liquide
                - rho       : Masse volumique
                - Beta      : Paramètre [beta_min, beta_max]
                - C_pl      : Capacité calorifique à phase liquide
                - N         : Nombre de noeuds
                - t         : Temps total de la simulation
        - S : Coefficient de stabilité fictif pour ajuster le pas de temps et comparer les résultats avec la méthode d'Euler explicite (par défaut à 0.5)
    Sorties:
        - t   : Array de temps de taille (Nt) 1D
        - x   : Array de position de taille (N) 1D
        - T   : Array de température de taille (Nt, N) 2D
    """
    #Définition des paramètres et des variables
    x = np.linspace(0, prm.L, prm.N)
    dx = x[1]- x[0]
    dt = (S*prm.C_pl*prm.rho*dx**2)/prm.k

    Nt = int(prm.t/dt)
    t = np.linspace(0, prm.t, Nt)

    T = np.empty((Nt, prm.N))
    #Boucle iterative en t et x pour créer la matrice de coefficients associées aux différences finies, i est l'indice de temps et j' est l'indice de position
    for i in range(Nt):
        if i == 0:
            T[i,:] = prm.T_s
        else:
            A = np.zeros((prm.N, prm.N))
            b = np.zeros(prm.N).T
            for j in range(prm.N):
                A[j,j] = 1 
                if j == 0:
                    A[j,j] = 1 
                    b[j] = prm.T_c
                elif j == prm.N - 1:
                    A[j,j] = 1
                    b[j] = prm.T_s
                elif T[i-1, j] < prm.T_l:                               #Si la température à l'itération précédente est inférieure à la température de liquidus, on utilise la capacité calorifique de changement de phase        
                    A[j,j-1] = -dt*prm.k/(C_psl(prm)*prm.rho*dx**2) 
                    A[j,j] = 1+ 2*dt*prm.k/(C_psl(prm)*prm.rho*dx**2)
                    A[j,j+1] = -dt*prm.k/(C_psl(prm)*prm.rho*dx**2) 
                    b[j] = T[i-1, j]
                else:                                                   #Sinon, on utilise la capacité calorifique à phase liquide
                    A[j,j-1] = -dt*prm.k/(prm.C_pl*prm.rho*dx**2) 
                    A[j,j] = 1+ 2*dt*prm.k/(prm.C_pl*prm.rho*dx**2)
                    A[j,j+1] = -dt*prm.k/(prm.C_pl*prm.rho*dx**2) 
                    b[j] = T[i-1, j]            
            T[i,:] = np.linalg.solve(A, b)                              #On résout le système d'équations pour trouver le profil de température à chaque pas de temps
    return t, x, T

def erreur(t, T_analytique, T_numerique):
    """Fonction qui calcule l'erreur quadratique moyenne à un temps spécifique, on choisit arbitrairement le temps final. Ceci est suffisant pour pouvoir le représenter dans un graph après.
    Entrées:
        - t                 : Array de temps de taille (Nt) 1D
        - T_analytique      : Array de températures solution analytique (Nt, N) 2D
        - T_numerique       : Array de températures solution numérique  (Nt, N) 2D
    Sorties:
        - dt        : pas de temps
        - erreur    : Erreur quadratique moyenne de la température au temps final
    """
    Nt = len(t)
    N = T_analytique.shape[1]
    dt = t[1] - t[0]    
    
    """ erreurt = np.empty(Nt) """
    
    """ for i in range(Nt):
        erreurt[i] = abs(np.mean(T_analytique[Nt-1,:] - T_numerique[Nt-1,:]))/N """
    erreur = np.linalg.norm(T_analytique[Nt-1,:] - T_numerique[Nt-1,:])/Nt

    return dt, erreur

def Vecteur_erreur(prm, schema = "explicite"):

    S = np.linspace(0.1,0.5, 10)

    dt_array = np.empty(len(S))
    erreur_array = np.empty(len(S))
    
    for i in range(len(S)):
        if schema == "explicite":
            t, x, T_numerique = Euler_explicite(prm, S[i])

            T_analytique = func_T(t, x, prm)

            dt_array[i], erreur_array[i] = erreur(t, T_analytique, T_numerique)

        elif schema == "implicite":
            t, x, T_numerique = Euler_implicite(prm, S[i])

            T_analytique = func_T(t, x, prm)
            dt_array[i], erreur_array[i] = erreur(t, T_analytique, T_numerique)
        else:
            return "Schema non reconnu, veuillez choisir entre 'explicite' et 'implicite'"
    return dt_array, erreur_array


















""" def residus_explicite(t, x, T, prm):
   
    Nt, N = T.shape
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    R = np.zeros((Nt - 1, N - 2))

    
    for n in range(Nt - 1):
        Cp_n = cp_effectif_array(T[n, :], prm)
        
        '''Actual calculs avec formules'''
        for j in range(1, N - 1):
            diffusion_exp = (prm.k / (prm.rho * Cp_n[j])) * (T[n, j+1] - 2*T[n, j] + T[n, j-1]) / dx**2

            derivee_temp = (T[n+1, j] - T[n, j]) / dt

            R[n, j-1] = derivee_temp - diffusion_exp

    norme_L2 = np.linalg.norm(R, axis=1)
    return R, norme_L2 """


""" def residus_implicite(t, x, T, prm):

    Nt, N = T.shape
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    R = np.zeros((Nt - 1, N - 2))

    for n in range(Nt - 1):
        Cp_n = cp_effectif_array(T[n+1, :], prm)

        '''Formules vues en cours'''
        for j in range(1, N - 1):
            diffusion_imp = (prm.k / (prm.rho * Cp_n[j])) * (T[n+1, j+1] - 2*T[n+1, j] + T[n+1, j-1]) / dx**2

            derivee_temp = (T[n+1, j] - T[n, j]) / dt

            R[n, j-1] = derivee_temp - diffusion_imp

    norme_L2 = np.linalg.norm(R, axis=1)
    return R, norme_L2
         """
""" #fonction pour calculer le résidus 
#Setup du CP
def cp_effectif_array(T_line, prm):
    Cp = np.empty_like(T_line)
    
    for j in range(len(T_line)):
        if T_line[j] < prm.T_l:
            Cp[j] = C_psl(prm)
        else:
            Cp[j] = prm.C_pl
            
    return Cp """