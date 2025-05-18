# -*- coding: utf-8 -*-
"""
Created on Thu May 15 22:30:56 2025

@author: USER
"""

# PIA
#Librerias
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
#Eligiendo modelo a graficar
print ("1 Pendulo simple, 2 maquina de atwood, 3 pendulo doble")
modelo = int(input("en ese mismo orden elegir del 1 al 3 para la grafica que desees"))

#Pendulo simple#
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

# Maquina Atwood
if modelo == 2:
    # Parametros
    g=9.8
    m1 = float(input("masa 1"))
    m2 = float(input("masa 2"))
    # Aceleracion Constante
    a = (m2 - m1) * g/ (m1 + m2)

    def atwood(t,y):
      x,v=y
      return [v,a] # dx /dt = v, dv/dt=a

    #Condiciones iniciales: x = 0 (posicion), v = 0 (velcidad)
    y0 = [0.0, 0.0]
    t_span = (0, 5)
    t_eval = np.linspace(*t_span, 300)
    # Resolver
    sol = solve_ivp(atwood, t_span, y0, t_eval=t_eval)
    #Graficar posicion
    plt.plot(sol.t, sol.y[0], label='posicion de m1 (x)')
    plt.plot(sol.t, sol.y[1], label='velocidad de m1 (v)')
    plt.xlabel('tiempo (s)')
    plt.ylabel('magnitud')
    plt.legend()
    plt.title('maquina de atwood')
    plt.grid(True)
    plt.show()
  
#pendulo_doble
if modelo == 3:
    # Parámetros físicos
    g = 9.81
    l = float(input("longitud"))  
    m = float(input("masa"))  
    
    def deriv(t, y):
        θ1, ω1, θ2, ω2 = y
    
        Δ = θ2 - θ1
    
        den1 = (m * (2 - np.cos(2*Δ)))
        den2 = (m * (2 - np.cos(2*Δ)))
    
        dω1_dt = (-g*(2)*np.sin(θ1)- g*np.sin(θ1 - 2*θ2) - 2*np.sin(Δ)*(ω2**2*l + ω1**2*l*np.cos(Δ))) / (l * den1)
        dω2_dt = (2*np.sin(Δ)*(ω1**2*l*m*(1) + g*m*np.cos(θ1) + ω2**2*l*m*np.cos(Δ))) / (l * den2)
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