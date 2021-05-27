import numpy as np
import matplotlib.pyplot as plt
#from scipy.fft import fft, fftfreq

# 1) Se inicializan las constantes tomadas de la tabla 3.
gbar_K = 36    # (1/mΩ*cm^2)
gbar_Na = 120  # (1/mΩ*cm^2)
g_L = 0.3      # (1/mΩ*cm^2)
E_K = -12      # (mV)
E_Na = 115     # (mV)
E_L = 10.6     # (mV)
C = 1          # (μF/cm^2)

# 2) Se define una función que generará los valores del tiempo.
def Genera_Tiempo(inicio, fin, cambio):
    '''
    Esta función crea el dominio del tiempo
    return: un arreglo t que contiene los tiempos a evaluar las ecuaciones.
    '''
    contador = inicio
    tiempo = []
    while (contador <= fin):
        tiempo += [contador]
        contador += cambio
    return tiempo

# 3) Se define una función que generaá los valores de la corriente de estímulo.
def J_estimulo(estimulo):
    '''
    Esta función crea la corriente de estímulo
    return: un arreglo corriente que coniene los estímulos
    '''
    lapso1 = 500    # (ms)
    lapso2 = 2000   # (ms)
    magnitud_j = estimulo # (mA/cm^2) del estímulo
    corriente = np.zeros(10000)
    for i in range(len(corriente)):
        if i <= lapso1:
            corriente[i] = magnitud_j
        elif lapso1 < i <= lapso2:
            corriente[i] = corriente[i]
        else:
            corriente[i] = magnitud_j/2
    return corriente

#-------------------4) Sección de funciones de las constantes dependidentes de voltaje--------------
def alfa_n(V):
    return 0.01 * ((10 - V) / (np.exp((10 - V) / 10) - 1))
def beta_n(V):
    return 0.125 * np.exp(-V / 80)
def alfa_m(V):
    return 0.1 * ((25 - V) / (np.exp((25 - V) / 10) - 1))
def beta_m(V):
    return 4 * np.exp(-V / 18)
def alfa_h(V):
    return 0.07 * np.exp(-V / 20)
def beta_h(V):
    return 1 / (np.exp((30 - V) / 10) + 1)
#---------------------------------------------------------------------------------------------------

# 5) Se define la función que ejecutará la aproximación del modelo
def Ejecuta_Modelo_Hodkin_Huxley(estimulo):

    # Se inicializan los valores de las constantes
    an = 0
    bn = 0
    am = 0
    bm = 0
    ah = 0
    bh = 0

    # Se inicia el arreglo que contendrá los valores para el potencial de acción en la membrana.
    potencial_accion = [] # (mV)
    conductancia_K = []
    conductancia_Na = []

    # Se llama a la densidad de corriente de estímulo.
    j_estim = J_estimulo(estimulo) # Recordar que es una variable tipo arreglo.

    # Se llama al los valores del timepo
    t_inicio = 0  # (ms)
    t_final = 100 # (ms)
    deltaT = 0.01 # (ms)
    tiempo = Genera_Tiempo(t_inicio, t_final, deltaT) # Recordar que es una variable tipo arrglo

    # Se inicializa una variable que construirá el potencial de acción
    voltajes = [0] # (mV)

    # Se definen las constantes con sus respectivas funciones dependientes de voltaje.
    contador = 0
    # Se definen las constantes con sus respectivas funciones dependientes de voltaje.
    an = alfa_n(voltajes[contador])
    bn = beta_n(voltajes[contador])
    am = alfa_m(voltajes[contador])
    bm = beta_m(voltajes[contador])
    ah = alfa_h(voltajes[contador])
    bh = beta_h(voltajes[contador])
    n = an/(an+bn)
    m = am/(am+bm)
    h = ah/(ah+bh)
    for t in range(len(tiempo)-1):
        an = alfa_n(voltajes[contador])
        bn = beta_n(voltajes[contador])
        am = alfa_m(voltajes[contador])
        bm = beta_m(voltajes[contador])
        ah = alfa_h(voltajes[contador])
        bh = beta_h(voltajes[contador])

        conductancia_Na += [(m**3) * gbar_Na * h * (voltajes[contador]-E_Na)]
        conductancia_K += [(n**4) * gbar_K * (voltajes[contador]-E_K)]
        jF = g_L *(voltajes[contador]-E_L)
        jIon = j_estim[contador] - conductancia_Na[contador] - conductancia_K[contador] -jF

    # Para resolver las ecuaciones diferenciales se utiliza el método
    # de diferencias hacia adelante, como se muestra a continuación.
        voltajes += [voltajes[contador] + deltaT*jIon/C]
        n = n + deltaT*(an *(1-n) - bn * n)
        m = m + deltaT*(am *(1-m) - bm * m)
        h = h + deltaT*(ah *(1-h) - bh * h)
        contador += 1

    conductancia_Na += [(m**3) * gbar_Na * h * (voltajes[contador]-E_Na)]
    conductancia_K += [(n**4) * gbar_K * (voltajes[contador]-E_K)]
    jF = g_L *(voltajes[contador]-E_L)
    jIon = j_estim[contador] - conductancia_Na[contador] - conductancia_K[contador] -jF

    # Finalmente se construye el potencial de acción
    for voltaje in voltajes:
        potencial_accion += [voltaje-70]    # El - 70 (mV) inidica el potencial de reposo o inicial.

    potencial_accion = voltajes
    # Graficación del potencial de acción:
    

    return([potencial_accion, conductancia_Na, conductancia_K, tiempo])

def calculo_de_frecuencias():
    estimulo = [10]
    frecuencias = [0]
    contador = 0 
    while(contador < 11):
        
        respuestas_modelo = Ejecuta_Modelo_Hodkin_Huxley(estimulo[contador])
        '''plt.plot(respuestas_modelo[3], respuestas_modelo[0], 'b', label='Potencial de acción')
        plt.xlabel("t (ms)")
        plt.ylabel("Potencial (mV)")
        plt.grid(True)
        plt.title("Modelo Hodkin-Huxley")
        plt.legend(loc='upper right')
        plt.show()

        plt.plot(respuestas_modelo[3], respuestas_modelo[1], 'r', label='Conductancia sodio')
        plt.plot(respuestas_modelo[3], respuestas_modelo[2], 'g', label='Conductancia potasio')
        plt.xlabel("t (ms)")
        plt.ylabel("Conductancia (1/mΩ*cm^2)")
        plt.grid(True)
        plt.title("Modelo Hodkin-Huxley")
        plt.legend(loc='upper right')
        plt.show()'''
        
        volt_max = max(respuestas_modelo[0])
        frecuencias += [busqueda_maximos(respuestas_modelo[0], volt_max)]
        estimulo += [estimulo[contador]+30]
        contador += 1
        

    plt.plot(estimulo, frecuencias, 'g', label='Respuesta de frecuencias')
    plt.xlabel("Densidad corriente (mA/cm^2)")
    plt.ylabel("Frecuencia")
    plt.grid(True)
    plt.title("Modelo Hodkin-Huxley")
    plt.legend(loc='upper right')
    plt.show()
    
def busqueda_maximos(voltajes, maximo):
    menor = voltajes[0]
    menor_siguiente = voltajes[0]
    presente = voltajes[0]
    resultado = 0
    for voltaje in voltajes:
        if(menor_siguiente > menor and menor_siguiente > presente):
            if(menor_siguiente >= maximo*0.4):
                resultado += 1
        menor = menor_siguiente
        menor_siguiente = presente
        presente = voltaje
    return resultado

    


#-----Ejecusión principal--------
calculo_de_frecuencias()