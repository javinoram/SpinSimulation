"""
input:
    - size (int): Numero de espines del sistema
output:
    - Vector de 3 posiciones, cada una almacena una lista con todos
    los operadores.
"""
def magnetic_vector(size: int) -> list:
    # Construir listas de string que representan los operadores
    x_pos = [ list("I"*size) for _ in range(size) ]
    y_pos = [ list("I"*size) for _ in range(size) ]
    z_pos = [ list("I"*size) for _ in range(size) ]

    # Reemplazar los operador de Identidad por el que corresponde
    # en la posicion que corresponde
    for i in range(size):
        x_pos[i][i] = "X"
        y_pos[i][i] = "Y"
        z_pos[i][i] = "Z"
    
    # Unir todos los elementos para que sigan la representacion
    x_pos =  [ "".join( elem ) for elem in x_pos ]
    y_pos =  [ "".join( elem ) for elem in y_pos ]
    z_pos =  [ "".join( elem ) for elem in z_pos ]
    return [x_pos, y_pos, z_pos]

"""
Dzyaloshinskii–Moriya effect, no esta considerado el efecto del signo, es decir,
hay que agregar que el segundo operador de cada lista tiene un signo negativos.
input:
    - index (tuple): Indices de los espines considerados en la interaccion.
    - size (int): Numero de espines del sistema.
output:
    - Vector de 3 posiciones con los operadores de Dzyaloshinskii–Moriya de cada eje.
"""
def antisymmetric_exchange_vector(index: tuple, size: int) -> list:
    # Construir listas de string que representan los operadores
    i_vector = [ list("I"*size) for _ in range(2) ]
    j_vector = [ list("I"*size) for _ in range(2) ]
    k_vector = [ list("I"*size) for _ in range(2) ]
    i,j=index

    i_vector[0][i] = "Y"
    i_vector[0][j] = "Z" 
    i_vector[1][i] = "Z"
    i_vector[1][j] = "Y" 

    j_vector[0][i] = "X"
    j_vector[0][j] = "Z" 
    j_vector[1][i] = "Z"
    j_vector[1][j] = "X" 

    k_vector[0][i] = "X"
    k_vector[0][j] = "Y" 
    k_vector[1][i] = "Y"
    k_vector[1][j] = "X" 

    # Unir todos los elementos para que sigan la representacion
    i_vector =  [ "".join( elem ) for elem in i_vector ]
    j_vector =  [ "".join( elem ) for elem in j_vector ]
    k_vector =  [ "".join( elem ) for elem in k_vector ]

    return [i_vector, j_vector, k_vector]


"""
Dzyaloshinskii–Moriya effect, no esta considerado el efecto del signo, es decir,
hay que agregar que el segundo operador de cada lista tiene un signo negativos
input:
    - indexes (list): Lista de tuplas con los indices de los pares de espines con la
    interaccion.
    - size (int): Numero de espines del sistema.
output:
    - Lista de los vectores de la funcion antisymmetric_exchange_vector.
"""
def set_antisymmetric_exchange_vector(indexes: list, size: int) -> list:
    operadores = []
    for index in indexes:
        operadores.append( antisymmetric_exchange_vector(index, size) )
    return operadores


"""
Calcular el producto punto del operador Si con Sj
input:
    - index (tuple): Indices de los espines considerados en la interaccion.
    - size (int): Numero de espines del sistema.
output:
    - Vector de 3 posiciones con los operadores de la interaccion Sij
"""
def sij_vector(index: tuple, size: int) -> list:
    x_pos = list("I"*size)
    y_pos = list("I"*size)
    z_pos = list("I"*size)
    i,j = index

    x_pos[i] = "X"
    x_pos[j] = "X"

    y_pos[i] = "Y"
    y_pos[j] = "Y"

    z_pos[i] = "Z"
    z_pos[j] = "Z"

    x_pos = "".join( x_pos )
    y_pos = "".join( y_pos )
    z_pos = "".join( z_pos )
    return [x_pos, y_pos, z_pos]


"""
Calcular el producto punto del operador Si con Sj de un conjunto de indices
input:
    - indexes (list): Lista de tuplas con los indices de los pares de espines con la
    interaccion.
    - size (int): Numero de espines del sistema.
output:
    - Lista de vector de 3 posiciones con los operadores de la interaccion Sij
"""
def set_sij_vector(indexes: list, size: int) -> list:
    operadores = []
    for index in indexes:
        operadores.append( sij_vector(index, size) )
    return operadores


"""
Toma el sistema ingresado y genera una copia no interactuante, es decir, el nuevo sistema esta compuesto de los 
dos pequeños y no interactuan entre si.
input:
    - operators (list): Lista de operadores con sus valores de constantes
    - copies (int): Numero de copias
output:
    - Lista de los nuevos operadores replicados en n copias
"""
def replicate_system(operators: list, copies: int) -> list:
    new_system = []
    size = len( operators[0][1] )
    for i in range(copies+1):
        for op in operators:
            new_system.append( [op[0], "I"*size*i + op[1] + "I"*size*(copies-i)] )
    return new_system

