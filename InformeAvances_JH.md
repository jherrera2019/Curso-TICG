# Informe de avances
---------------------------------------------
Se ha realizado el tutorial de *Expresión Diferencial* del Módulo 7 el cual tenía por objetivo introducir técnicas de el análisis de datos de microarreglos para detectar genes diferencialmente expresados producto de factores experimentales o de sus interacciones. El análisis ha sido ejecutado en el servidor remoto de la universidad bajo el lenguaje de programación UNIX y R. 
----------------------------------------------
En el estudio se perfilaron ocho ratones machos adultos de dos cepas, C57BL/6J y C57BL/6J-chrY<A/J/NaJ> , denominadas B y BY respectivamente. De cada cepa (genotipo), cuatro animales fueron castrados y cuatro animales quedaron intactos a modo control. El ARN se hibridizó a BeadChips Illumina MouseRef-8 v2.0 que contienen ocho microarreglos con 25,697 sondas cada uno. Solo se seleccionaron arbitrariamente 5000 sondas para este tutorial (Figura 1).

Objetivos del análisis de datos:

- Determinar si existe expresión diferencial entre genotipos.
- Determinar si existe expresión diferencial entre tratamientos.
- Evaluar las diferencias en la respuesta al tratamiento entre los dos genotipos (o sea si es que existe interacción)

En uno de los pasos que involucra la lectura de datos en bruto (no normalizados), se observa que no todas las sondas muestran la misma calidad al ser alineadas con el genoma de referencia utilizado.

```
annot <- read.delim("MouseRef-8_annot.txt") #se creó un objeto que contiene el genoma de referencia.
table(annot$ProbeQuality)
```

        Bad        Good    Good****    No match     Perfect  Perfect*** Perfect**** 
        289          60          15           5        4468          53         110 

Por lo tanto se agruparon las sondas "Bad" con "No match" y lo demás como "Good probes".

```probe_qc <- ifelse(annot$ProbeQuality %in% c("Bad", "No match"), "Bad probes","Good probes")```

## Control de calidad

Se crearon gráficos de caja coloreados por calidad de la sonda y calidad por tratamiento. 

![A](A.png)

**Figura A**. Diagramas de caja de datos sin procesar en escala log por microarreglo y calidad de sonda. 

En la figura A tenemos en el eje X las sondas divididas por calidad para cada array, en donde las de mala calidad tienen menor intensidad (rojo).

**FIgura B**. Diagramas de caja de datos en bruto por microarreglo.

En la figura B tenemos en el eje X, la posición de las matrices Illumina desde la A a H, mientras que en el eje Y el logaritmo en base 2 de la intensidad.  Las cajas están coloreadas según tratamiento. Se puede observar que hay una leve tendencia a mayores valores de intensidad en posiciones de castrados.

## Pruebas de expresión diferencial

Dado que el estudio utilizó un diseño factorial, esposible hacer variaas preguntas a partir de los datos y ser puestas a prueba mediante una matriz de contrastes entre grupos experimentales, permitiendo evaluar si esas comapraciones explican una proporción significativa de la varianza en los datos de expresión.

Los contrastes son vectores de coeficientes que cuando se multiplican con un vector de promedios por grupo experimental, crean contrastes que pueden evaluarse estadísticamente.

Se probó cada contraste utilizando 200 permutaciones, aplicando pruebas de F, utilizando una estimación de varianza residual por sonda (F1) y una estimación basada en contracción de varianza residual que utiliza información de múltiples sondas (Fs) 
