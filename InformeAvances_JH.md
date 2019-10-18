# Informe de avances
---------------------------------------------
Se ha realizado el tutorial de *Expresión Diferencial* del Módulo 7 el cual tenía por objetivo introducir técnicas de el análisis de datos de microarreglos para detectar genes diferencialmente expresados producto de factores experimentales o de sus interacciones. El análisis ha sido ejecutado en el servidor remoto de la universidad. 
----------------------------------------------
Se perfilaron ocho ratones machos adultos de dos cepas, C57BL/6J y C57BL/6J-chrY<A/J/NaJ> , denominadas B y BY respectivamente. De cada cepa (genotipo), cuatro animales fueron castrados y cuatro animales quedaron intactos a modo control. El ARN se hibridizó a BeadChips Illumina MouseRef-8 v2.0 que contienen ocho microarreglos con 25,697 sondas cada uno. Solo se seleccionaron arbitrariamente 5000 sondas para este tutorial (Figura 1).

Objetivos del análisis de datos:

- Determinar si existe expresión diferencial entre genotipos.
- Determinar si existe expresión diferencial entre tratamientos.
- Evaluar las diferencias en la respuesta al tratamiento entre los dos genotipos (es lo mismo que decir que existe interacción)

En uno de los pasos que involucra la lectura de datos en bruto, se observa que donde no todas las sondas muestran la misma calidad al ser alineadas con el genoma dado por el genoma de referencia utilizado.

```
annot <- read.delim("MouseRef-8_annot.txt") #se creó un objeto que contiene el genoma de referencia.
table(annot$ProbeQuality)
```

        Bad        Good    Good****    No match     Perfect  Perfect*** Perfect**** 
        289          60          15           5        4468          53         110 

Por lo tanto se agruparon las sondas "Bad" con "No match" y lo demás como "Good probes".

```probe_qc <- ifelse(annot$ProbeQuality %in% c("Bad", "No match"), "Bad probes","Good probes")```

Respecto al control de calidad, se crearon gráficos de caja coloreados por calidad de la sonda y calidad por tratamiento. En el eje X se observa la posición de las matrices desde la A a H
Y logaritmo en base 2 de la intensidad para cada muestra. estan divididos por anotación= malas y buenas dado por un paquete de biocontactor. 
mala calidad tienen menor intensidad, calidad lo det por secuencia. mapean secuencia contra el genoma de referencia --> malas mapean una vez o nunca o multipes veces en el genoma. *sugieren colorear por genotipo, plot of data raw by array colored by treatment*

AVG_Signal --> intensidad promedio entre beats
Detection_Pval --> importante para determinar si un gen está expresado o no, cmompara el valor de intesidad con la sonda negativa (no unión de moléculas de sample RNA) y calcula si lo registrado es background o señal

Matriz de diseño asocia información a las muestras, a los array, cada fila es un array, y tiene el nombre de la muestra el serial del bead chip.

-estableceremos contrastes entre grupos: estos son geno, trt, etc.
en la matriz de contrastes, 2 vectores se multiplicaron por 5; aquí voy a calcular constrastes entre el promedio de estos 2 grupos (B.C y B.I) (tienen en común el genotipo b ) y el promedio de los otros 2 (BY.C y BY.I). --> por lo que el contraste es la expresión promedio del grupo B menos la expresión promedio del grupo BY. Así, que si en un gen la expresión promedio del grupo B es > a la del BY va a tener un FC posiivo  

en análisis de expresión se usa en escala logaritmica, puesto que el enfoque es ver efectos a nivel aditivo y que puede ser modelado con un análsis anova, y hace se comporte con una distribución normal los datos. 

permutaciones para calcular la significancia de las diferencis de cada contraste, esto tiene una ventaja y es que no necesito asumir una función en particular. es otra forma de calcular la significancia de un valor de F.  En la misma línea calcula t, para calcular diferencias entre pares de grupos para cada fila de la matriz de contrastes, por lo que para cada contraste me dará un valor de t; significancia entre los tto, genotipos o la interacción. Fs es de las permutaciones



tenemos la distribución de los valores de p de F1 tabulares, falores de p de F1 permutado. Diferencias entre F1 y Fs es que el primero es el estándar que se hace gen por gen, osea toda la información para calcular el F de un gen x, en el segnudo se utiliza la información de error de otros GENES QUE se parecen al gen que esty midiendo en cuato al nivel de expresi´n que estoy midiendo (un shinkrage) --> forma artifical de aumentar el tamaño muestral, pq generalemne son bajos (podría tener un poco más de poder).

FDR --> ajuste con método adaptativo para el tutorial.

para los gráficos de venn
tenemos que cada fila en la matriz de datos (resultados) no son genes, son sondas, y un mismo gen puede estar representado por una misma sonda, hay repetición ,hay redundancia de los array pq quizás queremos detectar distintos trasncritos o el mismo mas de una vez. aquí resumiremos todas las sondas de un mismo gen en un resultado único por gen o sea, si cualquier sonda de un gen está detectando expresión diferencial, es que ese gen está diferencialmente expresado (caso particular para este tutorial).  
