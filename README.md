# DisulfideBridgeFinder

## Descripción
**DisulfideBridgeFinder** es un script en Python que identifica potenciales puentes disulfuro en estructuras proteicas obtenidas de archivos PDB. El programa es capaz de analizar tanto datos experimentales como predicciones de modelos generados por AlphaFold.

## Requisitos
- **Python 3.x**
- **Biblioteca BioPython** (instalable a través de pip: `pip install biopython`)

## Uso
1. **Ejecuta el script** `DisulfideBridgeFinder.py`.
2. **Ingresa la ruta completa** al archivo PDB que contiene la estructura proteica.
3. El programa analizará la estructura y mostrará los posibles puentes disulfuro encontrados, si los hay.
4. Si deseas analizar otro archivo, responde '**S**' o '**s**' cuando se te solicite.

## Funcionalidades
- **calcular_distancia(atomo1, atomo2)**: Calcula la distancia entre dos átomos.
- **calcular_diedro(atomo1, atomo2, atomo3, atomo4)**: Calcula el ángulo diedro entre cuatro átomos.
- **leer_pdb(ruta_archivo)**: Lee y parsea un archivo PDB.
- **es_alphafold(ruta_archivo)**: Verifica si el archivo PDB es de AlphaFold.
- **encontrar_cisteinas(estructura, es_alphafold)**: Identifica cisteínas y sus átomos de azufre en la estructura.
- **encontrar_puentes_disulfuro(atomos_azufre_cys)**: Encuentra posibles puentes disulfuro basados en la distancia y el ángulo diedro.
- **reportar_puentes(puentes)**: Imprime los posibles puentes disulfuro encontrados.

## Bibliotecas utilizadas
- **Bio.PDB**: Utilizada para leer archivos PDB y manipular estructuras proteicas.
  - **PDBParser**: Para leer y parsear archivos PDB.
  - **calc_dihedral**: Para calcular ángulos diedros entre átomos.
