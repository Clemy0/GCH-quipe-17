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
    
    
    C_psl = (prm.C_pl*(prm.T_l-prm.T_s) + prm.delta_H)/(prm.T_l - prm.T_s)
    return C_psl


def funcbeta(beta, prm):

    """ Foncction qui évalue la fonction étudiée pour trouver la valeur de beta qui correspond à un nombre de Stefan donné
    
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
        - f : Valeur de la fonction pour ce beta
    """
    Ste = Stefunc(prm)
    f = Ste/(np.sqrt(np.pi)) - beta*np.exp(beta**2)*erf(beta)
    return f

def funcbetaArray(prm):
    """Fonction étudiée pour trouver la valeur de beta qui correspond à un nombre de Stefan donné, pour génération de graphiques
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
        Sorties:
        -beta   : Array(len = 100) 1D de valeurs de beta entre beta_min et beta_max
        -f      : Array(len = 100) 1D de valeurs de la fonction pour un beta donné
        """
    
    beta = np.linspace(prm.Beta[0], prm.Beta[1], 100)
    f = np.array([funcbeta(b, prm) for b in beta])
    return beta, f


def bissection(beta0, beta1, prm):
    """Fonction qui utilise la méthode de bissection pour trouver la valeur de beta qui correspond à un nombre de Stefan donné
    
    Entrées:
        - beta0 : Valeur initiale de beta pour le début de l'intervalle de recherche
        - beta1 : Valeur initiale de beta pour la fin de l'intervalle de recherche
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
    Sorties:
        - beta_milieu : Valeur de beta qui correspond à un nombre de Stefan donné, trouvée par la méthode de bissection
    """
    f0 = funcbeta(beta0, prm)
    error = 1 

    while error > 1e-6:
        beta_m = (beta0 + beta1) / 2
        fm = funcbeta(beta_m, prm)
        if fm * f0 < 0:
            beta1 = beta_m
        else:
            beta0 = beta_m
        error = abs(fm)
    return beta_m

def func_T(x,t, prm):
    """Fonction qui donne la distribution de la température en fonction du temps
    Entrées:
        - x : Array de position de taille (N) 1D
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
            - N         : Nombre de noeuds
            - t         : Temps total de la simulation
    Sorties:
        -T  : Array de taille (Nt, N) de la distribution analytique de la température en fonction du temps 
    """
    Nx = len(x)
    Nt = len(t)
    alpha = prm.k/(prm.rho*prm.C_pl)
    T = np.empty((Nt,Nx))
    beta = bissection(prm.Beta[0], prm.Beta[1], prm)
    for i in range(Nt):
        for j in range(Nx):
            if j == 0:
                T[i, j] = prm.T_c
            elif i == Nx-1:
                T[i,j] = prm.T_s
            elif x[j] < 2*beta*np.sqrt(alpha*t[i]):
                T[i,j] = (prm.T_s - prm.T_c)*erf(x[j]/(2*np.sqrt(alpha*t[i])))/erf(beta) + prm.T_c 
            else:
                T[i,j] = prm.T_s
    return T





def Euler_explicite(prm, stability_value = 0.5): 
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
    Sorties:
        - t : Array de temps de taille (N) 1D 
        - x : Array de position de taille (N) 1D
        - T : Array de température de taille (N,N) 2D, où T[i,j] est la température au temps t[i] et à la position j
        
        """
    x = np.linspace(0, prm.L, prm.N)
    dx = x[1]- x[0]
    
    dt = (stability_value*prm.C_pl*prm.rho*dx**2)/prm.k
    Nt = int(prm.t/dt)
    t = np.linspace(0, prm.t, Nt)
    
    T = np.empty((Nt, prm.N))
    
    print(dx,dt)
        
    for i in range(Nt):
        for j in range(prm.N):
            if i == 0:
                T[i,j] = prm.T_s
            else:
                if j == 0:
                    T[i,j] = prm.T_c
                elif j == prm.N - 1:
                    T[i,j] = prm.T_s
                elif T[i-1, j] < prm.T_l:
                    T[i,j] = (dt*(prm.k/(C_psl(prm)*prm.rho))*(T[i-1, j+1] - 2*T[i-1,j] + T[i-1, j-1])/dx**2) + T[i-1, j]
                else:
                    T[i,j] = (dt*(prm.k/prm.rho)*(T[i-1, j+1] - 2*T[i-1,j] + T[i-1, j-1])/dx**2)/prm.C_pl + T[i-1, j]
    return t, x, T

#------------ Euler implicite essai Matthieu --------------------
def Euler_implicite(prm):
    """Fonction qui utilise la méthode d'Euler implicite pour résoudre 
    l'équation de la fusion d'un matériau en fonction du temps et de x.
    
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
    Sorties:
        - t   : Array de temps de taille (Nt) 1D
        - x   : Array de position de taille (N) 1D
        - T   : Array de température de taille (Nt, N) 2D
    """
    x  = np.linspace(0, prm.L, prm.N)
    dx = x[1] - x[0]

    dt = 0.99*(prm.rho*prm.C_pl / (2*prm.k))*dx**2
    Nt = int(prm.t / dt) + 1
    t = np.linspace(0, prm.t, Nt)

    T = np.empty((Nt, prm.N))

    for j in range(prm.N):
        T[0, j] = prm.T_s

    N_interieur = prm.N - 2

    for i in range(1, Nt):

        Cp = np.empty(N_interieur)
        for j in range(N_interieur):
            if T[i-1, j+1] < prm.T_l:
                Cp[j] = C_psl(prm)
            else:
                Cp[j] = prm.C_pl

        a = np.empty(N_interieur)
        for j in range(N_interieur):
            a[j] = prm.k*dt / (prm.rho*Cp[j]*dx**2)

        diag_sup = np.empty(N_interieur)
        diag_milieu = np.empty(N_interieur)        
        diag_inf = np.empty(N_interieur)
        B = np.empty(N_interieur)

        for j in range(N_interieur):
            diag_sup[j] = -a[j]
            diag_milieu[j] = 1 + 2*a[j]            
            diag_inf[j] = -a[j]
            B[j] = T[i-1, j+1]

        B[0] += a[0]*prm.T_c
        B[N_interieur-1] += a[N_interieur-1]*prm.T_s

        for j in range(1, N_interieur):
            facteur = diag_inf[j] / diag_milieu[j-1]
            diag_milieu[j] = diag_milieu[j] - facteur*diag_sup[j-1]
            B[j] = B[j] - facteur*B[j-1]

        T_interieur = np.empty(N_interieur)
        T_interieur[N_interieur-1] = B[N_interieur-1] / diag_milieu[N_interieur-1]
        for j in range(N_interieur-2, -1, -1):
            T_interieur[j] = (B[j] - diag_sup[j]*T_interieur[j+1]) / diag_milieu[j]

        T[i, 0]  = prm.T_c
        T[i, -1] = prm.T_s
        for j in range(N_interieur):
            T[i, j+1] = T_interieur[j]

    return t, x, T


def front_evolution(x, t, T, prm):
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
        - Array x_front_Analytique de position du front de fusion en fonction du temps, de taille (N) 1D
        - Array x_front_Numérique de position du front de fusion en fonction du temps, de taille (N) 1D
    """
    x_front_numerique = np.empty(len(t))
    x_front_analytique = np.empty(len(t))
    Nt = len(t)
    beta = bissection(prm.Beta[0], prm.Beta[1], prm)
    alpha = prm.k/(prm.rho*prm.C_pl)
    for i in range(Nt):
        i_front = np.argmax(T[i,:] < prm.T_l) 
        x_front_numerique[i] = x[i_front]
        x_front_analytique[i] = 2*beta*np.sqrt(alpha*t[i])
    return x_front_analytique, x_front_numerique

#------------ Euler implicite essai Aurélie --------------------
def euler_implicite(prm): 
    """Fonction qui utilise la méthode d'Euler implicite pour résoudre l'équation de la fusion d'un matériau en fonction du temps et de x
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
     Sorties:
        - t : Array de temps de taille (N) 1D 
        - x : Array de position de taille (N) 1D
        - T : Array de température de taille (N,N) 2D, où T[i,j] est la température au temps t[i] et à la position j
        
        """
# Discrétisation spatiale
    x = np.linspace(0, prm.L, prm.N)
    dx = x[1] - x[0]

    # Pas de temps (stable mais pas nécessaire en implicite)
    Stability_value = prm.rho * prm.C_pl / (2 * prm.k)
    dt = 0.99 * Stability_value * dx**2
    Nt = int(prm.t / dt)
    t = np.linspace(0, prm.t, Nt)

    # Initialisation température
    T = np.empty((Nt, prm.N))
    T[0, :] = prm.T_s
    T[0, 0] = prm.T_c
    T[0, -1] = prm.T_s

    Nint = prm.N - 2  # nombre de noeuds intérieurs

    # Boucle temporelle
    for n in range(1, Nt):

        T_old = T[n-1, :]

        # Cp pogné au temps précédent
        Cp = np.zeros(Nint)
        for j in range(Nint):
            if T_old[j+1] < prm.T_l:
                Cp[j] = C_psl(prm)
            else:
                Cp[j] = prm.C_pl

        # coefficient diffusion
        a = prm.k * dt / (prm.rho * Cp * dx**2)

        # Construction matrice A et vecteur b
        A = np.zeros((Nint, Nint))
        b = np.zeros(Nint)

        for j in range(Nint):

            A[j, j] = 1 + 2*a[j]

            if j > 0:
                A[j, j-1] = -a[j]

            if j < Nint - 1:
                A[j, j+1] = -a[j]

            b[j] = T_old[j+1]

        # Conditions limites
        b[0] += a[0] * prm.T_c
        b[-1] += a[-1] * prm.T_s

        T_new_int = np.linalg.solve(A, b)

        T[n, 0] = prm.T_c
        T[n, -1] = prm.T_s
        T[n, 1:-1] = T_new_int

    return t, x, T


        
#fonction pour calculer le résidus 
#Setup du CP
def cp_effectif_array(T_line, prm):
    Cp = np.empty_like(T_line)
    
    for j in range(len(T_line)):
        if T_line[j] < prm.T_l:
            Cp[j] = C_psl(prm)
        else:
            Cp[j] = prm.C_pl
            
    return Cp


def residus_explicite(t, x, T, prm):
   
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
    return R, norme_L2


def residus_implicite(t, x, T, prm):

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
        
