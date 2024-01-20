import math
import random

from fft_pot import fft,ifft,listprint

def bluestein(a):
    #Calculos basicos
    n = len(a)
    xi = pow(math.e,1j * math.tau/n)

    #Calculamos u(n) y v(n)
    u = [a[i] * pow(xi,-(i*i)/2) for i in range (n)]
    v = [pow(xi,(i *i)/2) for i in range (n)]
    v_estrella = [pow(xi,-(i * i)/2) for i in range (n)]

    #Extension a potencia de 2
    l = 2 ** math.ceil(math.log2(2 * n + 1))
    xi_l = math.cos(math.tau / l) + 1j * math.sin(math.tau / l)

    #Extension de u
    u_l = u + [0] * (l - n)

    #Extension de v
    aux = v[1:]
    aux.reverse()
    v_l = v + [0] * (l - 2*n + 1) + aux

    #Calculo de la convolucion

    dft_u_l = fft(u_l, xi_l) 

    dft_v_l = fft(v_l, xi_l)

    uv_l = [i*j for i,j in zip(dft_u_l,dft_v_l)]

    ift_l = ifft(uv_l , xi_l) 

    ift = ift_l[:n]

    # multiplicamos por el factor v*
    dft = [i*j for i,j in zip(v_estrella,ift)]

    return dft


# Comparar con la versión O(n²)
def check_dft(a, dft):
    dft_correcta = [
        sum([a[i] * pow(math.e, -math.tau * 1j * j * i / len(a)) for i in range(len(a))])
        for j in range(len(a))
    ]
    print(sum([abs(dft[i] - dft_correcta[i]) for i in range(len(a))]))

a = [i + 1 for i in range(11)]
import numpy as np
#print(fft([(1+0j), (1.0000000000000002-1.7320508075688772j), (-1.5000000000000009+2.5980762113533156j), 0, 0, 0, 0, 0], np.exp(math.tau*-1j/8)))
listprint(bluestein([i for i in range(1, 12)]))

# list(map(lambda x : str(round(x.real, ndigits=2)) + " + " + str(round(x.imag, ndigits=2)) + "i", bluestein(a)))