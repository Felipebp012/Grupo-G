import numpy as np
import matplotlib.pyplot as plt

def f(t,y):
    return np.cos(2*t) + np.sin(3*t)

def g(t,y):
    return -20*y + 20 *np.sin(t) +np.cos(t)

def extrapolacion(a,b,alpha,TOL,hmax,hmin):
    NK = np.array([2,4,6,8,12,16,24,32])
    Q = np.zeros((7,7))
    TO = a
    WO = alpha
    h = hmax
    y = np.zeros(7)
    T_k=[]
    W = []
    FLAG = 1
    for i in range(7):
        for j in range(i):
            Q[i,j] = (NK[i+1]/NK[j])**2
    while (FLAG == 1) :
        k=1
        NFLAG = 0
        while(k<=8 and NFLAG==0):
            HK = h/NK[k]
            T = TO
            W2 = WO
            W3 = W2 + HK * f(T,W2)#Primer paso de euler
            T = TO + HK
            for j in range(NK[k]-1):
                W1 = W2
                W2 = W3
                W3 = W1 + 2*HK * f(T,W2) # Metodo de punto medio 
                T = TO + (j+1) * HK
            y[k] = (W3 + W2 + HK * f(T,W3)) /2
            if (k>=2):
                j = k
                v = y[1]
                while(j>=2):
                    y[j-1] = y[j] + (y[j]-y[j-1])/(Q[k-1,j-1] -1)
                    j = j-1
                if abs(y[1]-v) <= TOL:
                    NFLAG = 1
            k+=1
        k-=1
        if (NFLAG == 0):
            
            h = h/2
            if h < hmin :
                print("hmin excedida")
                FLAG = 0 
        else:
            WO = y[1]
            TO += h
            T_k.append(TO)
            W.append(WO)
            if TO >= b:
                FLAG = 0
                print("Finalizado con exito")
            elif (TO + h > b):
                h = b - TO
            elif (k <= 3 and h < 0.5 * hmax):
                h = 2 * h
    
    return T_k, W, h

def soluciong(t):
    return np.sin(t) + 1/np.exp(20*t)

def solucionf(t):
    return (1/2) * np.sin(2*t) - (1/3) * np.cos(3*t) + 4/3

a = 0
b = 1
alpha = 1
TOL = 1e-9
hmax = 0.05
hmin = 0.005

T, W, h = extrapolacion(a, b, alpha, TOL, hmax, hmin)

t_real = np.linspace(a, b, 100)
y_real = solucionf(t_real)

print("T:", T)
print("W:", W)
print("h:", h)

plt.plot(T, W, label='Aproximación')
plt.plot(t_real, y_real, label='Solución Real')
plt.xlabel('t')
plt.ylabel('y(t)')
plt.title('Comparación de la Solución Aproximada y Real')
plt.legend()
plt.show()

tabla_data = []
yt = [solucionf(x) for x in T]
for i in range(len(T)):
    tabla_data.append([T[i],yt[i],W[i]])

fig, ax = plt.subplots(figsize=(10, 6))
ax.axis("off")

table = ax.table(
    cellText=tabla_data,
    colLabels=["t", "f(t)", "W = y(t)"],
    cellLoc="center",
    loc="center"
)

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.2)

plt.show()