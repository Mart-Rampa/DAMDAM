# yolo
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math

# il faut vérifier toutes les valeurs

tau = 9  # @valeur taux compression@ #[-]
D = 0.076  # @valeur alesage@ #[m]
C = 0.08  # @valeur course@ #[m]
L = 0.15  # @valeur longueur bielle@ #[m]
mpiston = 0.5  # @valeur masse piston@ #[kg]
mbielle = 0.3  # @valeur masse bielle@ #[kg]
# @valeur chaleur emise par fuel par kg de melange admis@ #[J/kg_inlet gas]
Q = 2.8

R = C/2
beta = L/R
gamma = 1.3


def Vc():
    return ((np.pi/D ** 2)/4)*2*R


def W(rpm):
    return rpm*60


def Qtot():
    return Q


def V(theta):
    zp = R*(1 - np.cos(theta) + beta - (beta**2 - (np.sin(theta))**2)**1/2)
    return Vc()/2 * zp + Vc()/(tau - 1)


def Q_fct(theta, thetaC, deltaThetaC):
    resultat = Q*0.5*(1 - np.cos(np.pi*((theta - thetaC)/deltaThetaC)))
    return resultat


def dVdTheta(theta):
    return Vc()*0.5*(np.sin(theta)+ ( ( (beta ** 2) -(np.sin(theta))**2)**(-0.5) ) * np.cos(theta) * np.sin(theta))


def dQDTheta(theta, thetaC, deltaThetaC):
    resultat = (Q_fct(theta, thetaC, deltaThetaC)*np.pi/(2*deltaThetaC)) * np.sin(math.pi*((theta-thetaC)/deltaThetaC))
    return resultat


def p_fct(rpm, s, theta, thetaC, deltaThetaC):

    def model(p, theta):
        dpdTheta = -gamma * (p/V(theta))*dVdTheta(theta) + (gamma - 1) * 1/V(theta) * dQDTheta(theta, thetaC, deltaThetaC)
        return dpdTheta

    # initial condition
    # vérifier que ca vaut s

    result_solve_ivp = solve_ivp(model, (theta[0], theta[-1]), [s], t_eval=theta).y[0]

    return result_solve_ivp


"""
INPUT :

rpm : la vitesse du moteur
s : le taux de suralimentation  (le taux de suralimentation multiplie la pression atmosphérique, supposée égale à 1 bar, pour obtenir la pression d'admission)
theta : l'angle de rotation (theta, en degré, sous forme d'un vecteur de -180 degrés à 180 degrés)
thethaC : l'angle d'allumage (thetaC en degré d'angle de villebrequin avant point mort haut - thetaC = 30 veut dire 30 degrés avant PMH),
deltaThetaC : la durée de la combustion (deltaThetaC en degré d'angle de villebrequin).

OUTPUT :

V_output, en [m3] : l'évolution du volume du cylindre en fonction de l'angle de rotation (theta)
Q_output, en [J] : l'évolution le l'apport de chaleur dans le cylindre en fonction de l'angle de rotation (theta)
F_pied_output et F_tete_output, en [N] : l'évolution des differentes forces s'appliquant sur la bielle en fonction de l'angle de rotation (theta)
p_output, en [Pa] : l'évolution de la pression dans le cylindre en fonction de l'angle de rotation (theta)
t en [m] : l'épaisseur critique, (voir schéma de l'énoncé)
"""


def myfunc(rpm, s, theta, thetaC, deltaThetaC):
    #VOTRE CODE

    #V_output = Vc()/2 * (1 - cos(theta) + beta - sqrt(beta^2 - (sin(theta))^2) ) + Vc()/(tau - 1)

    # j'ai mis ca plus haut
    #zp = R*(1 - cos(theta) + beta - sqrt(beta^2 - (sin(theta))^2))
    #V_output = Vc()/2 * zp + Vc()/(tau - 1)

    V_output = V(theta)

    # Fp
    p_theta = p_fct(rpm, s, theta, thetaC, deltaThetaC)
    fraction = ((math.pi*(D ** 2))/4)*p_theta
    RW2Cos = R*(W(rpm)) ** 2*np.cos(theta)

    F_pied_output = fraction - mpiston*RW2Cos
    F_tete_output = -fraction+(mpiston+mbielle)*RW2Cos

    Q_output = Q_fct(theta, thetaC, deltaThetaC)

    p_output = p_fct(rpm, s, theta, thetaC, deltaThetaC)

    SommeF = mbielle*RW2Cos

    t = 1

    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t)



rpm = 4106
s = 2.5
theta = np.linspace(-180, 180, 1000)
thetaC = 21
deltaThetaC = 48

thetaR = np.deg2rad(theta)
theta = thetaR


v, q, fpied, ftete, p, dq = myfunc(rpm, s, theta, thetaC, deltaThetaC)


# plot results
plt.plot(theta, p)
plt.xlabel('theta')
plt.ylabel('p')
plt.show()

plt.plot(theta, v)
plt.xlabel('theta')
plt.ylabel('v')
plt.show





print('end')
