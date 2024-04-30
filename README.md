# SpinSim
Libreria de codigo abierto de python para construir y calcular predicciones termodinamicas de sistemas de espines.

## Casos de uso
### Ejemplo de construir hamiltoniano
Considerando el siguiente hamiltoniano de espines 0.5 - 1.0 - 0.5

$\mathcal{H} = -J_1(S_1 \cdot S_2) -J_2(S_2 \cdot S_3) -J_3(S_1 \cdot S_3)$

Codigo de python para construir el hamiltoniano

```python
from spinsim.operators import set_sij_vector
from spinsim.hamiltonian import construct_hamiltonian

# Valores bases que son la descripcion del hamiltoniano
indices = [ (0,1), (1,2), (0,2) ]
exchanges = [J1, J2, J3]
espines = [0.5, 1.0, 0.5]

#Construye los productos Si * Sj siguiendo los indices
SiSjvectores = set_sij_vector(indices, 3) 

#Agrupar los terminos con su exchange
terminos = []
for J, op in zip(exchanges, SiSjvectores):
    terminos += [ [J, op[0]], [J, op[1]], [J, op[2]] ] 

#Obtener la matriz (numpy array) asociada al hamiltoniano
H = construct_hamiltonian(terminos, espines) 
```
De esta forma obtenemos la representacion matricial del hamiltoniano usando los indices del operador $S_iS_j$ (siempre respetando que $i$ sea menor o igual que $j$ ).

Dependiendo de la naturaleza del hamiltoniano, este se puede volver mas o menos complejo, pero la idea general es la simplicidad para construirlo, sin tener la necesidad de construir las matrices de pauli y hacer los productos kronecker respectivos.

En caso de querer las matrices asociadas a cada termino, se puede usar la funcion (considerando querer variar exchanges y no querer estar reconstruyendo todo en cada iteracion):

```python
from spinsim.hamiltonian import construct_term

#Retorna las matrices en formato sparse de cada termino
H = [ construct_term(t, espines) for t in terminos ]
```

### Ejemplo de calcular un observable
Considerando el hamiltoniano del ejemplo anterior, aca se usan unidades de *eV/K* para la constante de Boltzmann.

```python
from spinsim.thermodynamic import hermitian_specific_heat

temperatura = np.linspace(0.00001, 100, 10000)
valores = hermitian_specific_heat(H, temperatura) 
```

De esta manera, podemos calcular el calor especifico del hamiltoniano entre $0.00001$ a $100$ (K).
