## 4. Diseño de aplicación funcional con Pynq

El acelerador diseñado puede funcionar como un periférico disponible para ser usado desde la CPU. Para ello, es necesario programar una aplicación de software.
En el presente trabajo, se considera que los datos de entrada del controlador provienen de la CPU y la salida del controlador va debe ir a la CPU.
Una forma simple de acceder al acelerador como periférico es utilizar el framework Pynq, que permite hacerlo mediante Python.
<br><br>

Los Jupyter Notebooks en este directorio se utilizaron para verificar rápidamente la funcionalidad de los aceleradores MPC implementados. 
En general, es necesario ajustar el driver del acelerador de acuerdo a su interfaz. 
Adicionalmente, el Notebook `axi_time_mpc.ipynb` se utilizó para medir el overhead de la comunicación AXI entre la CPU y la FPGA, en implementaciones que incluyen el módulo sniffer.
