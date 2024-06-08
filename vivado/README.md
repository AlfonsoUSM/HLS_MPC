## 3. Implementación con Vivado

Para terminar de implementar el acelerador y obtener el archivo bitstream necesario para configurar la FPGA, hay que seguir los siguientes pasos (en Vivado 2022.2):

- Crear un proyecto en Vivado para la FPGA y tarjeta objetivo.
- Crear un `Block Design`
- En Tools>Settings>IP>Repository agregar dirección de la(s) IPs que se desean incluir.
- Agregar Processing System y ejecutar Block Automation
- Agregar la IP del acelerador (y del sniffer si se desea) y ejecutar Connection Automation
- Validate Design
- En Design Source>Create HDL Wrapper
- Generar Bitstream
- Exportar diseño con File>Export Block Design


En este directorio se incluye la IP sniffer desarrollada para registrar los intervalos de tiempo entre transacciones AXI.