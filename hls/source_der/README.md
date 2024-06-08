### 2.2. Archivos fuente de controlador MPC para DER

Los archivos en este directorio corresponden al código C++ de MPC, considerando la formulación del problema QP requerida para DER (específica para el tipo de restricciones).
Algunos pragmas utilizados sistemáticamente se incluyen directamente en el código C++ para afectar todas las pruebas, por ejemplo, aquellos que definen la interfaz del acelerador.

### Archivos

A continuación se describen los archivos principales:
- `mpc.cpp` y `mpc.hpp` : Código de MPC y ADMM.
- `system.hpp` : Declaración de constantes y parámetros del sistema (en particular dimensiones del problema QP).
- `system.cpp` : Definición de constantes, a reemplazar con los valores obtenidos de MatLab en el paso 1.
- `utils.hpp` : Código de funciones genéricas de operaciones entre vectores y matrices.