# Trabajo Fin de Grado (TFG) - Víctor Alonso Rodríguez
## The Chromatic surface brightness MODulation (CMOD)
### Metodología para el análisis de la evolución de los parámetros estructurales de las galaxias en función del desplazamiento al rojo

En este repositorio se pueden ecnontrar los [programas](https://github.com/victoralonsorodriguez/TFG/tree/main/Programas_creados) creados durante la elaboración del TFG y que han sido creados desde cero con el fin de analizar cómo cambian los parámetros estructurales de las galaxias con el desplzamiento al rojo. Todos los programas están escritos en Python aunque hay documentos ShellScript para automatizar procesos por lo que es recomendable ejcutar estos programas en Linux o en macOS.

Para poder hacer una prueba se pueden descargar los cubos de datos de la galaxia [M84](https://github.com/victoralonsorodriguez/TFG/tree/main/M84_PSF_large_V5) y tras descomprimirlos se puede iniciar el análisis ejecutando uno de los archivos .sh con la versión que se quiera utilizar. 

Los programas crearán varios subdirectorios, cada uno con el nombre del filtro. Dentro de ellos se encuentran los datos de salida del análisis. Lo resultados se pueden visualizar en los directorios que comienzan con el nombre de la galaxia seguido de _plot_ y _galfit_ o _ratios_ dependiendo de los resultados que se quieran visualizar. 

El respositorio se irá completando y actualizando con los últimos cambios y versiones del código. 