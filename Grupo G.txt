import numpy as np
import matplotlib.pyplot as plt

def f(t, y):
  return np.cos(2*t) + np.sin(3*t)

def metodo_extrapolacion(a, b, alfa, TOL, hmax, hmin):
  NK = np.array([2, 4, 6, 8, 12, 16, 24, 32])  # Define NK array
  Qi = np.zeros((7, 7))
  for i in range(7):
    for j in range(i + 1):
      Qi[i, j] = (NK[i + 1] / NK[j])**2  # Calcula Qi

  TO = a
  WO = alfa
  h = hmax
  FLAG = 1

  # Corrección: Inicializar T como una lista vacía
  T = []
  # Corrección: Eliminar la línea W = [WO]
  # W = [WO]

  while FLAG == 1:  # Bucle principal
    k = 1
    NFLAG = 0
    while k <= 8 and NFLAG == 0: # Bucle para cada k
      HK = h / NK[k - 1]
      # Corrección: No reasignar T, solo actualizar TO
      TO = TO
      W2 = WO
      W3 = W2 + HK * f(TO, W2)  # Primer paso de Euler
      TO = TO + HK

      for j in range(1, NK[k - 1]): # Bucle para cada j
        W1 = W2
        W2 = W3
        W3 = W1 + 2 * HK * f(TO, W2)  # Método del punto medio
        TO = TO + (j + 1) * HK

      yk = (W3 + W2 + HK * f(TO, W3)) / 2  # Corrección de extremo

      if k >= 2: # Calcula extrapolaciones
        j = k
        # Corrección: Usar yk en lugar de y1
        v = yk
        # Inicializa yj_1 con yk
        yj_1 = yk
        while j >= 2:
          # Corrección: Usar yk en lugar de yj
          yj_1 = yk + (yk - yj_1) / (Qi[k - 2, j - 2] - 1) # Extrapolación
          j -= 1
        if abs(yk - v) <= TOL:
          NFLAG = 1 # Precisión deseada alcanzada

      k += 1

    k -= 1

    if NFLAG == 0: # Resultado rechazado
      h /= 2
      if h < hmin:
        print('hmin excedida') # **Indentación corregida**
        FLAG = 0
    else: # Resultado aceptado
      WO = yk
      TO += h
      # Corrección: Agregar TO a la lista T
      T.append(TO)
      W.append(WO)

      if TO >= b:
        FLAG = 0 # Procedimiento completado
      elif TO + h > b:
        h = b - TO # Termina en t = b
      elif k <= 3 and h < 0.5 * hmax:
        h *= 2 # Incrementa tamaño de paso

  return T, W, h

def solucion_real(t):
  return (1/2) * np.sin(2*t) - (1/3) * np.cos(3*t) + 4/3

# Parámetros del problema
a = 0
b = 1
alfa = 1
TOL = 1e-9
hmax = 0.05
hmin = 0.005

# Aplicar el método de extrapolación
T, W, h = metodo_extrapolacion(a, b, alfa, TOL, hmax, hmin)

# Calcular la solución real
t_real = np.linspace(a, b, 100)
y_real = solucion_real(t_real)

# Imprimir los resultados
print("T:", T)
print("W:", W)
print("h:", h)

# Graficar los resultados
plt.plot(T, W, label='Aproximación')
plt.plot(t_real, y_real, label='Solución Real')
plt.xlabel('t')
plt.ylabel('y(t)')
plt.title('Comparación de la Solución Aproximada y Real')
plt.legend()
plt.show()
