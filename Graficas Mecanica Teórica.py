# -*- coding: utf-8 -*-
"""

@author: USER
"""
#Librerias utilizadas
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Elige el modelo a graficar
print("Este programa esta hecho para graficar Pendulo simple, pendulo doble, maquina de Atwood y maquina de Atwood doble")
modelo = int(input("En ese orden elige de 1 al 4 para la grafica que deseas"))


####################Pendulo simple#############################
# Ecuaciones del sistema: [theta, omega]
if modelo == 1:

    # Parámetros
    g = 9.81
    l = float(input("longitud"))  
    
    # Ecuaciones del sistema: [theta, omega]
    def pendulo(t, y):
        theta, omega = y
        dtheta_dt = omega
        domega_dt = - (g / l) * np.sin(theta)
        return [dtheta_dt, domega_dt]
    
    # Condiciones iniciales: theta0 = pi/4, omega0 = 0
    y0 = [np.pi / 4, 0]
    t_span = (0, 10)
    t_eval = np.linspace(*t_span, 300)
    
    # Solución
    sol = solve_ivp(pendulo, t_span, y0, t_eval=t_eval)
    
    # Graficar
    plt.plot(sol.t, sol.y[0], label='θ(t)')
    plt.plot(sol.t, sol.y[1], label='ω(t)')
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Ángulo / Velocidad angular')
    plt.legend()
    plt.grid(True)
    plt.title('Péndulo simple - Método de Lagrange')
    plt.show()
###################Pendulo doble##############################
if modelo == 2:
    # Parámetros físicos
    g = 9.81
    l = float(input("longitud"))  
    m = float(input("masa"))  
    
    def deriv(t, y):
        θ1, ω1, θ2, ω2 = y
    
        Δ = θ2 - θ1
    
        den1 = (m * (2 - np.cos(2*Δ)))
        den2 = (m * (2 - np.cos(2*Δ)))
    
        dω1_dt = (
            -g*(2)*np.sin(θ1)
            - g*np.sin(θ1 - 2*θ2)
            - 2*np.sin(Δ)*(ω2**2*l + ω1**2*l*np.cos(Δ))
        ) / (l * den1)
    
        dω2_dt = (
            2*np.sin(Δ)*(ω1**2*l*m*(1) + g*m*np.cos(θ1) + ω2**2*l*m*np.cos(Δ))
        ) / (l * den2)
    
        return [ω1, dω1_dt, ω2, dω2_dt]

    # Condiciones iniciales: [θ1, ω1, θ2, ω2]
    y0 = [np.pi / 2, 0, np.pi / 2, 0]
    
    # Intervalo de tiempo
    t_span = (0, 20)
    t_eval = np.linspace(*t_span, 2000)
    
    # Resolver ODE
    sol = solve_ivp(deriv, t_span, y0, t_eval=t_eval)
    
    # Graficar
    plt.plot(sol.t, sol.y[0], label='θ₁(t)')
    plt.plot(sol.t, sol.y[2], label='θ₂(t)')
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Ángulo (rad)')
    plt.legend()
    plt.grid(True)
    plt.title('Péndulo doble - Lagrangiano')
    plt.show()

###################Maquina Atwood#############################

if modelo == 3:
    # Parámetros
    g = 9.81
    m1 = float(input("masa 1"))  
    m2 = float(input("masa 1"))  
    # Aceleración constante
    a = (m2 - m1) * g / (m1 + m2)
    
    def atwood(t, y):
        x, v = y
        return [v, a]  # dx/dt = v, dv/dt = a
    
    # Condiciones iniciales: x = 0 (posición), v = 0 (velocidad)
    y0 = [0.0, 0.0]
    t_span = (0, 5)
    t_eval = np.linspace(*t_span, 300)
    
    # Resolver
    sol = solve_ivp(atwood, t_span, y0, t_eval=t_eval)
    
    # Graficar posición
    plt.plot(sol.t, sol.y[0], label='Posición de m1 (x)')
    plt.plot(sol.t, sol.y[1], label='Velocidad de m1 (v)')
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Magnitud')
    plt.legend()
    plt.title('Máquina de Atwood - Lagrangiano')
    plt.grid(True)
    plt.show()
################Maquina doble atwood#########################################
if modelo == 4:
    # Parámetros
    g = 9.81
    m1 = float(input("masa 1"))  
    m2 = float(input("masa 2"))  
    m3 = float(input("masa 3"))  
    
    def atwood_doble(t, y):
        x, dx, y_, dy = y
    
        # Ecuaciones obtenidas del lagrangiano (simplificadas)
        M = m1 + m2 + m3
        A = m1 + m2 + m3
        B = m2 - m3
        C = m2 + m3
    
        ddx = (m1 - m2 - m3) * g / A
        ddy = - (m2 - m3) * g / C
    
        return [dx, ddx, dy, ddy]
    
    # Condiciones iniciales: [x, dx/dt, y, dy/dt]
    y0 = [0.0, 0.0, 0.1, 0.0]  # Posición inicial de polea y diferencia entre m2 y m3
    
    # Tiempo
    t_span = (0, 10)
    t_eval = np.linspace(*t_span, 1000)

    # Resolver
    sol = solve_ivp(atwood_doble, t_span, y0, t_eval=t_eval)
    
    # Graficar
    plt.plot(sol.t, sol.y[0], label='x(t) - Pos. polea móvil')
    plt.plot(sol.t, sol.y[2], label='y(t) - Dif. entre m2 y m3')
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Posición (m)')
    plt.title('Máquina de Atwood Doble - Lagrangiano')
    plt.legend()
    plt.grid(True)
    plt.show()