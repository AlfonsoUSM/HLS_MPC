## 1. Diseño funcional en MatLab

El presente trabajo se centra en el procesamiento requerido para usar MPC. Se cuenta inicialmente con un modelo de la planta a controlar y se considera que el controlador MPC
ya fue diseñado. Sólo se analizan aquellos decisiones de diseño que afectan el costo computacional del algoritmo. 
<br><br>

Los archivos en este directorio permiten simular el funcionamiento de dos plantas (un motor DC y un filtro de Distributed Energy Resource) junto al controlador MPC correspondiente 
en diferentes escenarios y producir los dos archivos necesarios para el resto del flujo:
- Se genera un archivo `.cpp` con los valores de parámetros constantes que pueden calcularse antes de la operación de la planta;
- Se genera un archivo `.bin` que contiene la referencia funcional de entradas y salidas (golden reference) para la simulación del acelerador.


### Archivos

A continuación se describen los archivos principales:
- `mpc_motor_dense.m` : Simulación de controlador MPC para motor DC con formulación densa.
- `mpc_der_dense.m` : Simulación de controlador MPC para DER con formulación densa.
- `fx_dense_matrices.m` : Función que arma matrices del funcional de costo de acuerdo a la formulación densa de MPC.
- `fx_stationary.m` : Función que calcula el cambio de variables necesario para realizar seguimiento de referencia.
- `fx_qp_admm.m` : Función que implementa solver QP basado en el Alternating Directions Method of Multipliers.
