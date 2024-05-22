#****************************************************************************#
#                        Máster en Bioinformática                            #
#                        DisulfideBridgeFinder.py                            #
#                                                                            #
#                       Trabajo 1 Script en Python                           #
#                                                                            #
#  Escribir un script en Python (o cualquier otro lenguaje de programación)  #
#   que lea una estructura 3D de una proteína en formato PDB y reporte       #
#    POTENCIALES puentes di-sulfuro que pudiesen darse en esa estructura.    #
#  En concreto, los criterios para encontrar esos potenciales                #
#    puentes diS serían:                                                     #
#                                                                            #
#    · Distancia entre los azufres de las dos cisteinas: 1.5-2.5 A.          #
#    · Ángulo diedro C1-S1-S2-C2: 84-96º                                     #
#                                                                            #
# Para evitar regiones potencialmente móviles, los residuos                  #
#     con B-factor >30 A2 (si la estructura es experimental)                 #
#     o pLDDT <50 (si es un modelo de AlphaFold) se descartarán.             #
#                                                                            #
#  Version 1.0                                                               #
#                                                                            #
#  Realizado por: Alejandro Virués Morales                                   #
#  Asignatura: Bioinformática Estructural                                    #
#  Curso: 2023/2024                                                          #
#                                                                            #
#****************************************************************************#

import math
from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral, Vector

def mostrar_bienvenida():
    '''Imprime un mensaje de bienvenida y una breve explicación 
    del programa.'''
    print("*********************************************************")
    print("* Bienvenido a DisulfideBridgeFinder                    *")
    print("* Este programa identifica posibles puentes disulfuro   *")
    print("* en estructuras proteicas obtenidas de archivos PDB.   *")
    print("* Soporta tanto datos experimentales como predicciones  *")
    print("* de AlphaFold.                                         *")
    print("*********************************************************\n")

# Función para calcular la distancia entre dos átomos
def calcular_distancia(atomo1, atomo2):
    '''
    Calcula la distancia entre dos átomos.
    Input: atomo1 (objeto Atom), atomo2 (objeto Atom)
    Output: distancia (float)
    '''
    return atomo1 - atomo2

# Función para calcular el ángulo diedro entre cuatro átomos
def calcular_diedro(atomo1, atomo2, atomo3, atomo4):
    '''
    Calcula el ángulo diedro entre cuatro átomos.
    Input: atomo1 (objeto Atom), atomo2 (objeto Atom), atomo3 (objeto Atom),
    atomo4 (objeto Atom)
    Output: angulo (float, en grados)
    '''
    vec1 = atomo1.get_vector()
    vec2 = atomo2.get_vector()
    vec3 = atomo3.get_vector()
    vec4 = atomo4.get_vector()
    angulo = math.degrees(calc_dihedral(vec1, vec2, vec3, vec4))
    return angulo

# Leer y parsear el archivo PDB
def leer_pdb(ruta_archivo):
    '''
    Lee y parsea el archivo PDB.
    Input: ruta_archivo (str, ruta al archivo PDB)
    Output: estructura (objeto Structure)
    '''
    parser = PDBParser(QUIET=True)
    try:
        estructura = parser.get_structure("proteina", ruta_archivo)
        return estructura
    except Exception as e:
        print(f"Error al leer el archivo PDB: {e}")
        return None

# Verificar si el archivo PDB es de AlphaFold
def es_alphafold(ruta_archivo):
    '''
    Verifica si el archivo PDB es de AlphaFold.
    Input: ruta_archivo (str, ruta al archivo PDB)
    Output: booleano (True si es de AlphaFold, False en caso contrario)
    '''
    try:
        with open(ruta_archivo, 'r') as archivo:
            for linea in archivo:
                if 'ALPHAFOLD' in linea.upper():
                    return True
        return False
    except FileNotFoundError:
        print(f"Archivo no encontrado: {ruta_archivo}")
        return False
    except Exception as e:
        print(f"Error al verificar el archivo PDB: {e}")
        return False

# Identificar cisteínas y sus átomos de azufre
def encontrar_cisteinas(estructura, es_alphafold):
    '''
    Identifica cisteínas y sus átomos de azufre en la estructura.
    Input: estructura (objeto Structure), es_alphafold (booleano)
    Output: lista de tuplas (atomo, residuo, cadena)
    '''
    atomos_azufre_cys = []
    for modelo in estructura:
        for cadena in modelo:
            for residuo in cadena:
                if residuo.get_resname() == 'CYS':
                    if 'SG' in residuo:
                        atomo = residuo['SG']
                        if es_alphafold:
                            # Para AlphaFold, usar pLDDT (suponiendo que se 
                            # almacena en el B-factor)
                            pLDDT = atomo.get_bfactor()
                            if pLDDT >= 50:
                                atomos_azufre_cys.append((atomo, residuo,
                                                           cadena))
                        else:
                            # Para datos experimentales, usar B-factor
                            b_factor = atomo.get_bfactor()
                            if b_factor <= 30:
                                atomos_azufre_cys.append((atomo, residuo,
                                                           cadena))
    return atomos_azufre_cys

# Encontrar posibles puentes disulfuro
def encontrar_puentes_disulfuro(atomos_azufre_cys):
    '''
    Encuentra posibles puentes disulfuro basados en la distancia 
    y el ángulo diedro.
    Input: atomos_azufre_cys (lista de tuplas (atomo, residuo, cadena))
    Output: lista de posibles puentes (tuplas con información de los residuos,
    cadenas, distancia y ángulo)
    '''
    posibles_puentes = []
    for i, (atomo1, residuo1, cadena1) in enumerate(atomos_azufre_cys):
        for j, (atomo2, residuo2, cadena2) in enumerate(atomos_azufre_cys):
            if i >= j:
                continue
            distancia = calcular_distancia(atomo1, atomo2)
            if 1.5 <= distancia <= 2.5:
                # Calcular ángulo diedro
                c1 = residuo1['CB']
                s1 = atomo1
                s2 = atomo2
                c2 = residuo2['CB']
                angulo = calcular_diedro(c1, s1, s2, c2)
                if 84 <= abs(angulo) <= 96:
                    posibles_puentes.append(((residuo1, residuo2, cadena1,
                                               cadena2), distancia, angulo))
    return posibles_puentes

# Reportar los posibles puentes disulfuro
def reportar_puentes(puentes):
    '''
    Imprime los posibles puentes disulfuro encontrados.
    Input: puentes (lista de posibles puentes)
    Output: None (imprime los resultados)
    '''
    if puentes:
        for (residuo1, residuo2, cadena1, cadena2), distancia, angulo \
            in puentes:
            id_res1 = residuo1.get_id()[1]
            id_res2 = residuo2.get_id()[1]
            id_cadena1 = cadena1.get_id()
            id_cadena2 = cadena2.get_id()
            print(f"Posible puente disulfuro entre CYS {id_res1} en cadena "
                  f"{id_cadena1} y CYS {id_res2} en cadena {id_cadena2}")
            print(f"Distancia: {distancia:.2f} Å, Ángulo "
                  f"diedro: {angulo:.2f}°")
    else:
        print("No se encontraron potenciales puentes disulfuro "
              "en este fichero.")

def ejecutar_programa():
    '''Función principal para ejecutar el programa de análisis 
    de puentes disulfuro.'''
    mostrar_bienvenida()
    while True:
        # Solicitar la ruta del archivo PDB al usuario
        ruta_archivo_pdb = input("Introduce la ruta completa"
                                 " al archivo PDB: ")

        # Verificar si el archivo PDB es de AlphaFold
        alphafold = es_alphafold(ruta_archivo_pdb)

        # Leer estructura PDB
        estructura = leer_pdb(ruta_archivo_pdb)
        if estructura is None:
            print("No se pudo procesar el archivo PDB.")
        else:
            # Identificar cisteínas
            cisteinas = encontrar_cisteinas(estructura, alphafold)

            # Encontrar posibles puentes disulfuro
            puentes_disulfuro = encontrar_puentes_disulfuro(cisteinas)

            # Reportar resultados
            reportar_puentes(puentes_disulfuro)

        # Preguntar si se desea analizar otro archivo
        while True:
            continuar = input("¿Deseas analizar otro archivo? "
                              "(S/N): ").strip().upper()
            if continuar in ['S', 'N']:
                break
            print("Por favor, introduce 'S' para sí o 'N' para no.")
        
        if continuar == 'N':
            print("Gracias por usar DisulfideBridgeFinder. ¡Hasta luego!")
            break

# Ejecutar el programa
if __name__ == "__main__":
    ejecutar_programa()
