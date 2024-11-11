def extrapolacion(a,b,alpha,TOL,hmax,hmin):
    NK = [2,4,6,8,12,16,24,32]
    TO = a
    WO = alpha
    h = hmax
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
            for j in range(HK-1):
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
            WO = y[1]
            TO += h
            if TO >= b:
                FLAG = 0
                print("Finalizado con exito")
            elif (TO + h > b):
                h = b - TO
            elif (k <= 3 and h < 0.5 * hmax):
                h = 2 * h
            return(TO,WO,h)
        else:
            h = h/2
            if h < hmin :
                print("hmin excedida")
                FLAG = 0 
    return