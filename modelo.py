import math
import numpy as np
import matplotlib.pyplot as plt


gbar_K=36;
gbar_Na=120;
g_L=.3;
E_K = -12;
E_Na=115;
E_L=10.6;
C=1;

def estimulos(tiempo1=500, valor1=50, tiempo2=2000, valor2=50, tiempoGeneral=10000):
    resultado = []
    contador = 0
    while(contador <= tiempo1):
        resultado += [valor1]
        contador += 1
    while(contador <= tiempo2):
        resultado += [0]
        contador += 1
    while(contador <= tiempoGeneral):
        resultado += [valor2]
        contador += 1
    return resultado

def inicio():
    an = 0;
    bn = 0;
    am = 0;
    bm = 0;
    ah = 0;
    bh = 0;

    resultadoPot = []

    estim = estimulos()
    tiempoF = 100; #el tiempo esta en milisegundos
    cambioTiempo = 0.01; #tiempo en milisegundos
    tiempo = valoresTiempo(0, tiempoF, cambioTiempo)
    
    voltajes=[0]
    contador = 0;
    an = alphaN(voltajes[contador])
    bn = betaN(voltajes[contador])
    am = alphaM(voltajes[contador])
    bm = betaM(voltajes[contador])
    ah = alphaH(voltajes[contador])
    bh = betaH(voltajes[contador])
    n = an/(an+bn);
    m = am/(am+bm);
    h = ah/(ah+bh);
    for t in range(len(tiempo)-1):
        an = alphaN(voltajes[contador])
        bn = betaN(voltajes[contador])
        am = alphaM(voltajes[contador])
        bm = betaM(voltajes[contador])
        ah = alphaH(voltajes[contador])
        bh = betaH(voltajes[contador])

        iNa = (m**3) * gbar_Na * h * (voltajes[contador]-E_Na)
        iK = (n**4) * gbar_K * (voltajes[contador]-E_K)
        iF = g_L *(voltajes[contador]-E_L)
        iIon = estim[contador] - iNa - iK -iF

        voltajes += [voltajes[contador] + cambioTiempo*iIon/C]
        n = n + cambioTiempo*(an *(1-n) - bn * n)
        m = m + cambioTiempo*(am *(1-m) - bm * m)
        h = h + cambioTiempo*(ah *(1-h) - bh * h)
        contador += 1

    for voltaje in voltajes:
        resultadoPot+=[voltaje-70]
    plt.plot(tiempo, resultadoPot, 'r', label='Potencial de Membrana')
    plt.xlabel("t (ms)")
    plt.ylabel("Potencial (mV)")
    plt.grid(True)
    plt.title("Modelo Hodkin-Huxley")
    plt.legend(loc='upper right')
    plt.show()
    



def valoresTiempo(inicio, fin, cambio):
    contador = inicio
    resultado = []
    while(contador <= fin):
        resultado += [contador]
        contador+=cambio
    return resultado

def alphaN(V):
    return  0.01 * ( (10-V) / (math.exp((10-V)/10)-1) );

def betaN(V):
    return 0.125*math.exp(-V/80); 

def alphaM(V):
    return 0.1*( (25-V) / (math.exp((25-V)/10)-1) );
    
def betaM(V):
    return 4*math.exp(-V/18);

def alphaH(V):
    return 0.07*math.exp(-V/20);

def betaH(V):
    return 1/(math.exp((30-V)/10)+1);



inicio()