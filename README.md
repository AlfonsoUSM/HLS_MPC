# Implementación de MPC usando HLS

Este repositorio forma parte de la tesis titulada <b>Implementación de Aceleradores de Cómputo en Lógica Reconfigurable para Aplicaciones de Control Predictivo por Modelo Utilizando Síntesis de Alto Nivel</b>,
para optar al grado de Magister en Ciencias de la Ingeniería Electrónica.
<br><br>

Este proyecto se desarrolla dentro del grupo de investiación de los profesores Gonzalo Carvajal, César Silva y Juan Agüero, de la Universidad Técnica Federico Santa María.


## Flujo

El flujo utilizado para el diseño de controladores con aceleración en hardware usando Síntesis de Alto Nivel se basa en las herramientas de AMD-Xilinx.
A partir de un diseño funcional del algoritmo en MatLab, se verifica el control deseado de la planta y se obtienen valores de parámetros y referencias funcionales.
<br><br>

Luego, se transcribe la sección on-line del algoritmo a C++ para realizar la síntesis de alto nivel usando Vitis HLS. Agregando diferentes <b>pragmas</b> al código se observan cambios en las estimaciones de latencia y uso de recursos.
Además, Vitis HLS permite simular la implementación en RTL usando Vivado, para verificar que se mantenga la funcionalidad. 
Los distintos pasos se pueden ejecutar manualmente desde un proyecto de Vitis HLS, pero también es posible usar scripts tcl para probar distintas combinaciones de pragmas y facilitar la exploración de soluciones.
El diseño en RTL elegido se exporta para realizar la implementación física en Vivado, obtienendo el bitstream que configurará la FPGA. 
<br><br>

Por último, se diseña la aplicación de software que utilizará el acelerador desde la CPU del controlador. En este trabajo se consideraron dos alternativas. 
Por una parte, el uso del framework Pynq facilita el diseño de la aplicación y el uso del acelerador con Python. Por otra parte, una aplicación baremetal diseñada con Vitis permite evitar el overhead de usar un sistema operativo.
En ambos casos, resulta necesario considerar la latencia de la comunicación entre la CPU y la FPGA. Para esto se diseña un método para medir los tiempos de comunicación. 


## Estructura

A continuación se describe el contenido de los directorios:
- matlab: implementación en matlab de mpc implícito usando admm
- hls: archivos fuente (C++) de mpc para motor DC, DER y ejemplo de multiplicación matriz-vector, junto con sctipts tcl
- pynq: jupyter notebooks para utilizar el acelerador usando Pynq