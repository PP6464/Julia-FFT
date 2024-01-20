import numpy as np
import math

def listprint(x):
    for i in x:
        print(i)

# De clase: FFT para potencias de 2
def fft(pol, xi):
    n = len(pol)
    if n == 1:
        return pol
    
    #print("FFT of:")
    #listprint(pol)
    #print(f"with root of unity: {xi}")

    pol_par = [pol[2 * i] for i in range(n // 2)]
    pol_impar = [pol[2 * i + 1] for i in range(n // 2)]
    res_par = fft(pol_par, xi ** 2)
    res_impar = fft(pol_impar, xi ** 2)
    res = [0j] * n

    #print("X evens: ")
    #listprint(res_par)
    #print("X odds: ")
    #listprint(res_impar)

    for i in range(n // 2):
        #print(f"res[{i}]: ")
        #print(res_par[i] + xi ** i * res_impar[i])
        #print(f"res[{i + n // 2}]")
        #print(res_par[i] - xi ** i * res_impar[i])
        res[i] = res_par[i] + xi ** i * res_impar[i]
        res[i + n // 2] = res_par[i] - xi ** i * res_impar[i]
    
    return res


# De clase: Inversa de FFT para potencias de 2
def ifft(vec, xi):
    n = len(vec)
    res = fft(vec, 1 / xi)
    for i in range(n):
        res[i] /= n
    return res

#listprint(ifft([(0.4999999999999987+0.8660254037844383j), (3.04088222241137-3.0355339059327404j), (11.232050807568882+4.133974596215561j), (5.408607520371811-4.0355339059327395j), (16.5+0.8660254037844384j), (-2.2370346451180003+4.035533905932738j), (7.767949192431123-5.866025403784442j), (5.787544902334822+3.0355339059327404j)], (0.7071067811865476+0.7071067811865476j)))