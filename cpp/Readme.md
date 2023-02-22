# Model Predictive Control C++

Se implementa MPC con C++ básico, es decir, sin utilizar threads ni alguna librería que permita paralelismo. Además se miden y guardan los tiempos de ejecución.

## Requerimientos

Si se quiere ejecutar en otro procesador se requieren los compiladores cruzados para C++ y fortran.

### Instalar compiladores cruzados

Si se requiere que este código se ejecute en un procesador ARM hay que instalar los siguientes compiladores cruzados:

  `sudo apt install gcc-aarch64-linux-gnu`

  `sudo apt install gfortran-aarch64-linux-gnu`

## Muestras (samples)
Se tienen que tener las muestras generadas por el código *motorMPCReferenceTracking.m* con el horizonte de predicción que se va a utilizar. Las instrucciones para generar las muestras se encuentran en *GenerateSamplesMATLAB/Readme.md*.

## Makefile

Se provee un makefile para compilar el código. Se le pueden dar argumentos para cambiar el horizonte de predicción, el tipo de dato, el linear solver a utilizar, o si se quiere hacer compilación cruzada para un procesador ARM cortex-a53. Si no se indica algún parámetro se utiliza el por defecto indicado en el header file *source/specs.h*.

Existen 3 Linear Solvers de los cuales se debe escoger uno al momento de compilar utilizando *-DLS=X*, donde *X* puede tomar uno de los siguientes valores:
  - MINRES = 1
  - CGRAD = 2
  - CHOL = 3


Para compilar para ejecutar localmente, con datos *float*, horizonte N=4 y linear solver "MINRES":

  `make DEFINES="-DELEM_FLOAT -DN_QP=4  -DLS=1"`

Para compilar para ejecutar en procesador ARM cortex-a53, con datos *double*, horizonte N=10 y linear solver "CHOL":

  `make arm DEFINES="-DELEM_DOUBLE -DN_QP=10  -DLS=3"`

Para compilar para ejecutar localmente y parámetros por defecto:

  `make`

## Ejecutar el Código

Se puede compilar utilizando el makefile provisto o el script *compile.sh* que compila para varios horizontes de predicción y tipos de datos (*float* y *double*). 

Para ejecutar el código se puede utilizar el siguiente comando (remplazando con el N correspondiente), donde los tiempos van a ser guardados en *t.csv*:

 `./build/bin/runner ../GenerateSamplesMATLAB/samples/samplesMPC_N[N].bin t.csv`

 El output en la consola debería ser parecido al siguiente:

 ```console
MPC testbench
Finished processing 10000 samples
Number of differences between expected and calculated:  0
Threshold: 0.0001

Mean time processing a sample:          67.0074878787879µs
Standard deviation processing a sample: 14.2091157762009µs
Max running processing a sample:        302.861µs
Min running processing a sample:        61.319µs
```

 O bien si se compiló con el script *compile.sh* se puede ejecutar el script *benchmark.sh*.

