# Informe módulo 2.- 

Joel Herrera Soto.- Curso Tópicos de investigación en 
Estudiante 1 año de Magister en Genética
Profesor: Ricardo Verdugo

---------------------------------------------
# Introducción

El concepto de expresión diferencial hace alusión al cambio en la expresión de genes que es explicado por la interación de agentes ambientales y los genes, o sea cuando hay expresión diferencial, es que hubo interacción genético ambiental, permitiendo entender la variación fenotípica que pueda tener el modelo de estudio. Gracias a las aplicaciones bioinformáticas, es que podemos estudiar los cambios de expresión, permitiéndonos clasificar los genes en vías de señalización particulares y por lo tanto entregarnos un panorama de los proocesos que se vean alterados bajo los distintos escenarios experimentales a los que estén inmersos nuestro modelo de estudio.

Para este módulo se ha realizado el tutorial de **[Expresión Diferencial](https://github.com/AliciaMstt/BioinfinvRepro/blob/master/Unidad7/Tutorial_de_expresion_diferencial_en_R.md)** del Módulo 7 el cual tenía por objetivo introducir técnicas de el análisis de datos de microarreglos para detectar genes diferencialmente expresados producto de factores experimentales o de sus interacciones. Se generaron figuras para describir control de calidad de las sondas y para describir el análsis de expresión diferencial mismo. El análisis ha sido ejecutado en el servidor remoto de la universidad bajo el lenguaje de programación UNIX y R.  

----------------------------------------------------
# Desarrollo

En el estudio se perfilaron ocho ratones machos adultos de dos cepas, C57BL/6J y C57BL/6J-chrY<A/J/NaJ> , denominadas B y BY respectivamente. De cada cepa (genotipo), cuatro animales fueron castrados y cuatro animales quedaron intactos a modo control. El ARN se hibridizó a BeadChips Illumina MouseRef-8 v2.0 que contienen ocho microarreglos con 25,697 sondas cada uno. Sólo se seleccionaron arbitrariamente 5000 sondas para este tutorial (Figura 1).

Objetivos del análisis de datos:

- Determinar si existe expresión diferencial entre genotipos.
- Determinar si existe expresión diferencial entre tratamientos.
- Evaluar las diferencias en la respuesta al tratamiento entre los dos genotipos (o sea si es que existe interacción)

## Lectura de datos en bruto

En la lectura de los datos en bruto (no normalizados) se observa que no todas las sondas muestran la misma calidad al ser alineadas con el genoma de referencia utilizado.

```
annot <- read.delim("MouseRef-8_annot.txt") #se creó un objeto que contiene el genoma de referencia.
table(annot$ProbeQuality)
```

        Bad        Good    Good****    No match     Perfect  Perfect*** Perfect**** 
        289          60          15           5        4468          53         110 

Por lo tanto se agruparon las sondas "Bad" con "No match" y lo demás como "Good probes". Son malas sondas porque hace match con secuencias repetidas, regiones intergénicas o intrónicas, o es poco problable que entregue una señal sensiblle y específica para cualquier trancrito al ser alineados contra el genoma de referencia (Barbosa-Morais et al. 2010). 

```probe_qc <- ifelse(annot$ProbeQuality %in% c("Bad", "No match"), "Bad probes","Good probes")```

## Control de calidad

Se crearon gráficos de caja coloreados por calidad de la sonda y calidad por tratamiento. 

![A](A.png)

**Figura 1**. Diagramas de caja de datos sin procesar en escala log por microarreglo y calidad de sonda. 

En la figura A tenemos en el eje X las sondas divididas por calidad para cada array, en donde las de mala calidad tienen menor intensidad (rojo).

![B](B.png)

**FIgura 2**. Diagramas de caja de datos en bruto por microarreglo.

En la figura B tenemos en el eje X, la posición de las matrices Illumina desde la A a H, mientras que en el eje Y el logaritmo en base 2 de la intensidad.  Las cajas están coloreadas según tratamiento. Se puede observar que hay una leve tendencia a mayores valores de intensidad en posiciones de castrados.

## Pruebas de expresión diferencial

Dado que el estudio utilizó un diseño factorial, es posible hacer varias preguntas a partir de los datos y ser puestas a prueba mediante una matriz de contrastes entre grupos experimentales, permitiendo evaluar si esas compraciones explican una proporción significativa de la varianza en los datos de expresión.

Los contrastes son vectores de coeficientes que cuando se multiplican con un vector de promedios por grupo experimental, crean contrastes que pueden evaluarse estadísticamente.

Se probó cada contraste utilizando 200 permutaciones, aplicando pruebas de F, utilizando una estimación de varianza residual por sonda (F1) y una estimación basada en contracción de varianza residual que utiliza información de múltiples sondas (Fs).

```
test.cmat <- matest(madata, fit.fix, term="Group", Contrast=cmat, n.perm=200, 
+                     test.type = "ttest", shuffle.method="sample", verbose=TRUE)
Doing F-test on observed data ...
Doing permutation. This may take a long time ... 
Finish permutation #  100 
Finish permutation #  200 

```
La función `matest` del paquete R/MAanova puede utilizar una matriz de contrastes y evaluar si esas comparaciones explican una proporción significativa de la varianza en los datos de expresión. 
El objeto de clase `madata` incluye las sondas que sólo detectaron transcritos (contiene la matriz de datos y la tabla de diseño experimental) las cuales son 4706. El objeto de clase `fitmaanova` llamado fit.fix contiene el modelo ajustado a ANOVA. El argumento `contrast` define la matriz de contraste a ser utilizada. El argumento `verbose`muestra mensajes del progreso de cálculo. 

## Contar genes expresados diferencialmente.

En este experimento se buscaba resolver las preguntas ¿existe un efecto de interacción entre el genotipo y el tratamiento sobre la expresión génica en los cardiomiocitos? además de evaluar la naturaleza de la interacción. ¿Ambos genotipos responden al tratamiento pero en direcciones opuestas? ¿O el tratamiento tiene un efecto en uno de los genotipos y no en el otro? 

Hay que tener en cuenta que cada sonda está detectando diferentes señales biológicas, por lo que se cuenta un gen como "seleccionado" si se selecciona alguno de los transcritos (sondas). Esto se logra con la base de datos de transcritos de referencia [RefSeq](http://www.ncbi.nlm.nih.gov/RefSeq/), la cual es curada, para genes conocidos y se diseño para evitar la redundancia.

Utilizamos la librería limma para crear diagramas de Venn. Luego contaremos los genes para cada combinación de efectos marginales y de interacción.

```
> library(limma)
> Counts.DE <- vennCounts(Genes.DE)
> print(Counts.DE)
  FDR.Geno FDR.Trt FDR.Int Counts
1        0       0       0    672
2        0       0       1     29
3        0       1       0    404
4        0       1       1     25
5        1       0       0    388
6        1       0       1     41
7        1       1       0    288
8        1       1       1     57
attr(,"class")
[1] "VennCounts"
```
Contar los genes DE entre niveles de un factor condicional en el otro factor.

```
> Counts.Int_Geno <- vennCounts(Genes.Int_Geno)
> print(Counts.Int_Geno)
  FDR.Geno_I FDR.Geno_C Counts
1          0          0      3
2          0          1    100
3          1          0     35
4          1          1     14
attr(,"class")
[1] "VennCounts"
> 
> Counts.Int_Trt  <- vennCounts(Genes.Int_Trt) 
> print(Counts.Int_Trt)
  FDR.Trt_B FDR.Trt_BY Counts
1         0          0      0
2         0          1    103
3         1          0     31
4         1          1     18
attr(,"class")
[1] "VennCounts"

```

Figura 3. Genes Diferencialemente expresados (DE) por efectos marginales y de interacción.

Figura 4. Genes DE por efectos de interacción, divididos por tratamiento (izquierda) y genotipo (derecha) .
