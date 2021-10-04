import numpy as np
import math
import matplotlib.pyplot as plt
from copy import deepcopy


def calcular_variables(perfil):
    # Asignamos los valores a las coordenadas de cada punto X e Y
    X = perfil[0]
    Y = perfil[1]

    # Ahora calculamos la longitud de onda de la señal
    lamda = Vh / frecuencia
    print('Longitud de Onda ', lamda, ' metros')

    # calculamos la posición (X) de inicio del perfil, donde (y)=Zs
    posini = np.where(Y < Ys)[0][0]
    RXini = X[posini]
    print('Posición inicial del receptor = ', RXini, 'metros')

    # Ahora hacemos lo mismo con la posición del receptor al final del perfil.
    posfin = np.where(Y < Ys)[0][-1]
    RXfin = X[posfin]
    print('Posición final del receptor = ', RXfin, ' metros')

    # Para calcular la distancia total del perfil solo tenemos que restar la posicion inicial a la posicion final
    D = RXfin - RXini
    print('Distancia total del perfil = ', D, ' metros')

    # Y si sabemos que haremos una traza cada 10 metros (Dx), entonces...
    Ntrazas = int(D / Dx)
    print('Distancia entre trazas = ', Dx, ' metros')
    print('Numero total de trazas = ', Ntrazas)

    """
    Para seguir tenemos que definir un tiempo razonable para capturar la traza, 
    de tal forma que la onda electromagnetica pueda llegar hasta los puntos mas 
    lejanos del perfil y retornar al sensor. Lo calculamos de esta manera:
    sabiendo que el perfil tiene 3190.1m de longitud (D), y asumiendo que las ondas 
    tienen que ir y volver, asumiremos que el tiempo mínimo tiene que ser suficiente 
    para recorrer esta distancia 2 veces:
    """
    Dtr = D * 2
    print('Distancia transmisor receptor = ', Dtr, ' metros')

    """
    Para ir sobre seguro vamos a utilizar una distancia un poco mayor a Dtr, la 
    distancia  Dtrc 8000m por ejemplo nos irá bién
    """
    Dtrc = 6500
    print('Distancia transmisor receptor corregida = ', Dtrc, ' metros')

    # Ahora si calculamos el tiempo mínimo:
    tiempoRadargrama = Dtrc / Vh
    print('Tiempo min de la Traza = ', tiempoRadargrama, ' milisegundos')

    # Procedemos a inicializar la variable que contendrá nuestro radargrama
    min_interval = tiempoRadargrama / Nmuestras
    return min_interval, Ntrazas, RXini, RXfin, X, Y, posini, posfin, tiempoRadargrama


def calcular_distancia(posX_suelo, posY_suelo, posX_perf, posY_perf):
    distancia = math.sqrt((posX_suelo - posX_perf) ** 2 + (posY_suelo - posY_perf) ** 2)
    return distancia


def calcular_tiempo(distancia, velocidad):
    tiempo = distancia / velocidad
    return tiempo


def calcular_tiempo_propagacion(distancia, posY_n):
    if posY_n < Ys:
        t_propag = calcular_tiempo(distancia, Vh)
    else:
        t_propag = calcular_tiempo(distancia, Va)
    return t_propag


def calcular_dist_total(posx_emisor, posx_perf, posy_perf):
    posx_receptor = posx_emisor + Dx

    # Se calcula la distancia total entre emisor-perfil (ida) y perfil-receptor (vuelta)
    d_emisor = calcular_distancia(posx_emisor, Ys, posx_perf, posy_perf)
    d_receptor = calcular_distancia(posx_receptor, Ys, posx_perf, posy_perf)
    return d_emisor + d_receptor


def calcular_intervalo(t, min_interval):
    if t < min_interval:
        pos = 0
    elif t > min_interval * 512:
        pos = None
    else:
        pos = int(t / min_interval)
    return pos


def calc_todos_tiempos(todas_posX, todas_posY, min_interval, Ntrazas, RXini):
    radargrama = np.zeros([Nmuestras, Ntrazas])
    radargrama1, radargrama2, radargrama3 = deepcopy(radargrama), deepcopy(radargrama), deepcopy(radargrama)
    for n_trazas in range(Ntrazas):
        pos_x = RXini + (n_trazas * Dx)
        for n_posc in range(len(todas_posX)):
            dist_tot = calcular_dist_total(pos_x, todas_posX[n_posc], todas_posY[n_posc])
            t_posc = calcular_tiempo_propagacion(dist_tot, todas_posY[n_posc])
            posicion = calcular_intervalo(t_posc, min_interval)
            radargrama1[posicion, n_trazas] = radargrama1[posicion, n_trazas] + 1
            radargrama2[posicion, n_trazas] = radargrama2[posicion, n_trazas] + (1 / dist_tot ** 2)
            radargrama3[posicion, n_trazas] = radargrama3[posicion, n_trazas] + math.log(1 / dist_tot ** 2, 10)
    return radargrama1, radargrama2, radargrama3


def graficar_radargrama(x, y, pos_ini, pos_final, radargrama_i, profundidad_muestras, details_plot):
    # Grafique los datos intensidad 1
    plt.figure()
    plt.figure(1)
    plt.imshow(radargrama_i, extent=[pos_ini, pos_final, profundidad_muestras[-1], profundidad_muestras[0]])
    plt.plot(x, y, "-b", label=details_plot['legend'])
    plt.legend(loc="lower left")
    plt.xlabel(details_plot['x_axis'])
    plt.ylabel(details_plot['y_axis'])
    plt.title('Intensity = {} , Perfil {}'.format(details_plot['intensidad'], details_plot['perfil']))
    plt.show()
    plt.savefig("mygraph{}{}.png".format(details_plot['intensidad'], details_plot['perfil']))
    plt.close()


def ejecutar_calculo(perfil, perf_ref):
    min_interval, Ntrazas, RXini, RXfin, X, Y, posini, posfin, tiempoRadargrama = calcular_variables(
        perfil)

    radargrama_i1, radargrama_i2, radargrama_i3 = calc_todos_tiempos(X, Y, min_interval, Ntrazas, RXini)
    # print("Todos los tiempos len : {}".format(len(todos_tiempos)))

    # Calculamos la profundida de cada muestra considerando la velocidad de la onda en el  hielo
    ProfundidadMuestras = Ys - np.linspace(0, Vh * tiempoRadargrama / 2, Nmuestras, endpoint=True)

    # Radargrama intensidad 1
    plot1 = {'x_axis': 'Valores X', 'y_axis': 'Valores Y', 'perfil': perf_ref, 'intensidad': '1',
             'legend': 'Leyenda N'}
    graficar_radargrama(X[posini:posfin], Y[posini:posfin], RXini, RXfin, radargrama_i1, ProfundidadMuestras, plot1)
    # Radargrama intensidad 1/d^2
    plot2 = {'x_axis': 'Valores X', 'y_axis': 'Valores Y', 'perfil': perf_ref, 'intensidad': '1d2',
             'legend': 'Leyenda N'}
    graficar_radargrama(X[posini:posfin], Y[posini:posfin], RXini, RXfin, radargrama_i2, ProfundidadMuestras, plot2)
    # Radargrama intensidad log(1/d^2)
    plot3 = {'x_axis': 'Valores X', 'y_axis': 'Valores Y', 'perfil': perf_ref, 'intensidad': 'log10(1d2)',
             'legend': 'Leyenda N'}
    graficar_radargrama(X[posini:posfin], Y[posini:posfin], RXini, RXfin, radargrama_i3, ProfundidadMuestras, plot3)

################ CONSTANTES #####################
c = 300e6
frecuencia = 7.5e6
Ys = 1200
Dx = 10
Vh = 160e6
Va = c
Nmuestras = 512

################ PERFIL 1 #####################
perfil1 = np.transpose(np.genfromtxt('perfil.csv', delimiter=',', skip_header=1))
ejecutar_calculo(perfil1, "1")

################ PERFIL 2 ######################
perfil2 = np.transpose(np.genfromtxt('perfil2.csv', delimiter=',', skip_header=1))
ejecutar_calculo(perfil2, "2")
