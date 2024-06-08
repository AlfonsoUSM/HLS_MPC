## 2. Implementación de acelerador usando HLS
La síntesis de alto nivel permite obtenera descripción RTL que implemente un algoritmo escrito en un lenguaje de alto nivel (en este caso en C++).
Comenzamos esta etapa con el código C++ correspondiente al procesamiento on-line de MPC y los archivos de parámetros y referencias obtenidos de la etapa anterior.
<br><br>

La herramienta Vitis HLS permite al usuario afectar la arquitectura y el grado de paralelización del acelerador mediante <b>pragmas</b>. 
Los reportes de resultados muestran el impacto de estos pragmas en latencia del acelerador y uso de recutsos de la FPGA objetivo.
Mediante un proceso iterativo, se puede encontrar una combinación de pragmas tal que la implementación cumpla con las especificaciones.
<br><br>

Adicionalmente, es posible automatizar el proceso de exploración escribiendo scripts `.tcl` para probar mayor número de combinaciones de pragmas.
 

### Estructura

A continuación se describe el contenido de los directorios:
- source: archivos fuente de mpc para motor DC.
- source_der: archivos fuente de mpc para DER.
- source_mvm: ejemplo para analizar impacto de pragmas en multiplicación matriz-vector. 
- scripts: archivos `.tcl` utilizados para exploración rápida de soluciones de HLS.