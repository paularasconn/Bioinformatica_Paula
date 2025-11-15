# PaulaRasconSousa_Trabajo2.R
# Trabajo final Bioinformática - Curso 25/26
# Análisis de parámetros biomédicos por tratamiento
getwd()

# 1. Cargar librerías (si necesarias) y datos del archivo "datos_biomed.csv". (0.5 pts)
install.packages("readr")
library(readr)

# Antes de empezar debemos asegurarnos de que el archivo que vamos a usar se ha guardado bien, para ello observamos la ruta de archivo:
ruta_archivo <- "C:/Users/Paula/OneDrive/Escritorio/datos_biomed.csv"

# Leemos los datos y los guardamos:
datos_biomed <- read_csv(ruta_archivo)


# 2. Exploración inicial con las funciones head(), summary(), dim() y str(). ¿Cuántas variables hay? ¿Cuántos tratamientos? (0.5 pts)
# Para asegurarnos de que todo esta bien, observamos las primeras filas del archivo:
head(datos_biomed)

# Summary lo usamos para obtener un resumen estadístico por columnas:
summary(datos_biomed)

# dim lo usamos para obtener dimensiones del data frame, es decir, las filas y columnas
dim(datos_biomed)  
# el resultado obtenido: [1] 100   5 (tratamientos o filas y variables o columnas)

# También podemos observar la estructura del archivo:
str(datos_biomed)

# Extraemos el número de variables o columnas, que en este caso vemos que son 5:
num_variables <- ncol(datos_biomed)
cat("Número de variables:", num_variables, "\n")

# Extraemos también el número de tratamientos o filas, que en este caso son 100:
num_tratamientos <- nrow(datos_biomed)
cat("Número de tratamientos:", num_tratamientos, "\n")


# 3. Una gráfica que incluya todos los boxplots por tratamiento. (1 pt)
# Cargamos la librería ggplot2 que nos permite usarla para gráficas avanzadas
install.packages("ggplot2")
library(ggplot2)

# A continuación, convertimos los datos a formato largo usando tidyr de la siguiente manera:
install.packages("tidyr")
library(tidyr)

# Seleccionamos las columnas numéricas ( aquellas que hacen referencia a la glucosa, presion y colesterol) y tratamiento:
datos_largos <- datos_biomed %>%
  pivot_longer(cols = c(Glucosa, Presion, Colesterol),
               names_to = "Variable",
               values_to = "Valor")

# Una vez realizado lo anterior, creamos el boxplot: Valor vs Tto, separado por Variable
ggplot(datos_largos, aes(x = Tratamiento, y = Valor, fill = Tratamiento)) +
  geom_boxplot() +
  facet_wrap(~Variable, scales = "free_y") +  
  theme_bw() +                           
  labs(title = "Boxplots por tratamiento para cada variable",
       x = "Tratamiento",
       y = "Valor")

# En cuanto a lo que podemos observar en el boxplot:
	# En relación al colesterol mencionamos que le fármaco A, tiene los valores más altos de colesterol de los tres grupos.
	# En cuanto a la glucosa, obseervamos que los tres grupos tienen valores similares (poca diferencia), por lo que los tratamientos no parecen producir diferencias importantes en los niveles de glucosa.
	# Por último, en relación a la presión, vemos que el fármaco B aumenta ligeramente la presión en comparación con el fármaco A, pero notoriamente frente al grupo placebo.


# 4. Realiza un violin plot (investiga qué es). (1 pt)
# Este tipo de gráficos son de gran utilidad para comparar la distribucion de distintas variables, en este caso, los tratamiento con el colesterol, glucosa y la presión.
# Usamos de nuevo datos en formato largo teniendo en cuenta lo siguiente:
	# 'Variable' = tipo de medición,
	# 'Valor' = valor numérico, 
	# 'Tratamiento' = grupo

# con los siguientes códigos definimos el tipo de gráfico, la mediana y el rango intercuartílico, y el tema que emplearemos.
ggplot(datos_largos, aes(x = Tratamiento, y = Valor, fill = Tratamiento)) +
  geom_violin(trim = FALSE) +                   
  geom_boxplot(width = 0.1, fill = "white") +   
  facet_wrap(~Variable, scales = "free_y") +    
  theme_bw() +                             
  labs(title = "Violin plots por tratamiento para cada variable",
       x = "Tratamiento",
       y = "Valor")

# Sobre el gráfico podemos comentar lo siguiente:
	# En cuanto al colesterol, el fármaco A se correlaciona con unos niveles mayores, en comparación con los otros dos grupos
	# En cuanto a la glucosa, tanto el fármaco A como el B demuestran una leve disminución en sus niveles en comparación con el grupo placebo.
	# Por último, vemos que la presión se ve aumentada de forma muy similar en los fármacos A y B.


# 5. Realiza un gráfico de dispersión "Glucosa vs Presión". Emplea legend() para incluir una leyenda en la parte superior izquierda. (1 pt)
# Antes de nada definimos los colores de los distintos tratamientos o grupos a comparar:
colores <- c("FarmacoA" = "red", "FarmacoB" = "green", "Placebo" = "blue")

# Una vez asignados los colores, creamos un gráfico base de la siguiente manera:
# Explicación de los códigos empleados a continuación: col para darle color  cada tratamiento, pch para definir el tipo de punto empleado en el gráfico, xlab e ylab para definir los eejr y main para añadirle un título.
plot(datos_biomed$Glucosa, datos_biomed$Presion,
     col = colores[datos_biomed$Tratamiento],  
     pch = 16,                                 
     xlab = "Glucosa",                        
     ylab = "Presión",                         
     main = "Dispersión Glucosa vs Presión")  

# Ahora, añadimos una leyenda en la esquina superior izquierda 
#con legend definimos la posicion de la leyenda y les ponemos nombre a los tratamientos. También hay que asignar loscolores dentro de la leyenda.
legend("topleft",                                         
       legend = levels(as.factor(datos_biomed$Tratamiento)),  
       col = c("red","green","blue"),                         
       pch = 16)       

# En cuanto a los resultados del gráfico podemos afirmar que debida a la gran dispersión y separacion de los puntos o valores, no podemos definir una tendencia con seguridad.


# 6. Realiza un facet Grid (investiga qué es): Colesterol vs Presión por tratamiento. (1 pt)
# Se trata de un tipo de gráfico comparativo que se usa en bibliotecas de visualización como ggplot2, en estos casos, se dividen los datos en subgraficos según las variables y se comparan de forma visual.
# Cargamos ggplot2 antes de nada:
library(ggplot2)

# Creamos un gráfico con facetas por tratamiento
# en este caso, los puntos de dispersión son definidos según los colores de los tratamientos y de la misma forma que hemos hecho anteriormente, definimos el tema.
ggplot(datos_biomed, aes(x = Colesterol, y = Presion)) +
  geom_point(aes(color = Tratamiento), size = 2) +   
  facet_grid(. ~ Tratamiento) +                      
  theme_minimal() +                                  
  labs(title = "Colesterol vs Presión por tratamiento",
       x = "Colesterol",
       y = "Presión") +
  theme(legend.position = "none")    

# El gráfico nos muestra la distribución y dispersión según tratamiento y grupo.


# 7. Realiza un histogramas para cada variable. (0.5 pts)
# Primero de todo, instalamos los paquetes necesarios ya explicados en anteriores apartados:
install.packages("ggplot2")
install.packages("tidyr")

library(ggplot2)
library(tidyr)

# definimos las columnas numéricas
vars_numericas <- c("Glucosa", "Presion", "Colesterol")

# Ahora, al igual que en otros puntos, reorganizamos los datos en formato largo con tidyr::pivot_longer
datos_largo <- pivot_longer(
  datos_biomed,
  cols = all_of(vars_numericas),
  names_to = "Variable",
  values_to = "Valor"
)

# A continuación, podemos crear los histogramas deseados:
ggplot(datos_largo, aes(x = Valor, fill = Variable)) +
  geom_histogram(bins = 30, alpha = 0.8, color = "black") +
  facet_wrap(~ Variable, scales = "free") +
  labs(
    title = "Distribución de valores - Histogramas por variable",
    x = "Valor medido",
    y = "Frecuencia"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# En los histogramas proporcionamos observamos que en el caso de la glucosa, se ve una distribucion mas simétrica que en el resto de grupos. Aparentemente no se obseervan valores anormales ni sesgados.
# este tipo de gráficos nos muestran la distribución de los datos por variables.


# 8. Crea un factor a partir del tratamiento. Investiga factor(). (1 pt)
# Convertimos la variable tratamiento en un factor con levels y el siguiente código:
datos_biomed$Tratamiento_factor <- factor(datos_biomed$Tratamiento,
                                          levels = c("Placebo", "FarmacoA", "FarmacoB"))

# Verificamos los datos y generamos una tabla para ver la cantidad de sujetos por grupo de tratamiento:
str(datos_biomed$Tratamiento_factor)
table(datos_biomed$Tratamiento_factor)  

# Observamos que la distribución de sujetos por grupo es correcta debido a que tienen un tamaño muestral similar:
	# 33 sujetos en Placebo
	# 36 sujetos en Fármaco A
	# 31 sujetos en Fármaco B


# 9. Obtén la media y desviación estándar de los niveles de glucosa por tratamiento. Emplea aggregate() o apply(). (0.5 pts)
# Media de Glucosa por tratamiento
media_glucosa <- aggregate(Glucosa ~ Tratamiento, 
                           data = datos_biomed, 
                           FUN = mean)
media_glucosa

# Desviación estándar de Glucosa por tratamiento
sd_glucosa <- aggregate(Glucosa ~ Tratamiento, 
                        data = datos_biomed, 
                        FUN = sd)
sd_glucosa

# Los resultados obtenidos son los siguientes:
	# En cuanto a la media:
		# FarmacoA 110.4750
		# FarmacoB 105.5839
		# Placebo 103.0455
	# En cuanto a la desviación estandar:
		# FarmacoA 13.58309
		# FarmacoB 12.15064
		# Placebo 17.18486

# En cuanto a la interpretación de los valores, podemos comentar que el fármaco A al tener unos valores de la media más elevados, interpretamos que este tratamiento podría asociarse con valores de glucosa superiores en comparación con los otros grupos.Por otro lado, vemos que el grupo placebo contiene una desviación estandar más alta, lo que indica que hay una mayor dispersión y que el grupo es más heterogéneo.


# 10. Extrae los datos para cada tratamiento y almacenalos en una variable. Ejemplo todos los datos de Placebo en una variable llamada placebo. (1 pt)
# Primero, usamos subset() para filtrar los datos y quedarnos solo con las filas que cumplen la condición deseada:
# Grupo placebo:
placebo <- subset(datos_biomed, Tratamiento == "Placebo")

# Grupo fármaco A
farmacoA <- subset(datos_biomed, Tratamiento == "FarmacoA")

# Grupo fármaco B
farmacoB <- subset(datos_biomed, Tratamiento == "FarmacoB")

# Para confirmar que está bien lo que hemos hecho, al igual que hemos hecho en puntos anteriores, observamos las primeras filas de cada uno:
head(placebo)
head(farmacoA)
head(farmacoB)

# En cuanto a los resultados obtenidos, mencionamos:
	# El fármaco A tiene los valores más altos de colesterol y mayor variabilidad general.
	# El fármaco B muestra valores intermedios con la glucosa algo baja.
	# El grupo placebo muestra una dispersión mayor a los otros dos grupos sin valores extremos o incoherentes.

# 11. Evalúa si los datos siguen una distribución normal y realiza una comparativa de medias acorde. (1 pt)
# Primero cargamos la librería para pruebas estadísticas:
install.packages("ggpubr")
library(ggpubr)

# Una vez instalados los paquetes, evaluamos la normalidad por grupo con shapiro:
shapiro.test(placebo$Glucosa)  
shapiro.test(farmacoA$Glucosa) 
shapiro.test(farmacoB$Glucosa) 
# Entre estos resultados:
	Grupo		W		p-value
	Placebo	0.967		0.402
	Fármaco A	0.965		0.301
	Fármaco B	0.966		0.411
	# Esto nos indica que al tener todos los grupos un p-value mayor a 0.5, los valores cumplen el requisito para aplicar pruebas paramétricas al no ser rechazada la hipótesis de normalidad.

# 12. Realiza un ANOVA sobre la glucosa para cada tratamiento. (1 pt)
# Comparativa de medias usando ANOVA
anova_glucosa <- aov(Glucosa ~ Tratamiento, data = datos_biomed)
summary(anova_glucosa)
# Entre estos resultados:
	Fuente		Df	SumSq	  Mean Sq	F value	Pr(>F)
	Tratamiento		2	989	  494.3	2.358 	0.1
	# Esto nos indica que el al ser el valor pr (avlor P asociado a F) mayor a O.5, no hay diferencias estadisticas notorias entre los grupos de tratamiento.


# Post-hoc: comparación por pares con TukeyHSD, este tipo de análisis lo realizamos para hacer comparaciones por pares de grupos:
tukey_glucosa <- TukeyHSD(anova_glucosa)
tukey_glucosa

# El resultado se basa en una tabla en el que se mencionan los siguientes parámetros:
	# diff: es la diferencia media entre dos grupos que comparamos
	# lwr / upr: se trata del límite inferior y superior del intervalo de confianza del 95%
	# p adj: es el valor p ajustado para comparaciones de varios grupos simultaneas
# Este análisis post-hoc con TukeyHSD nos sirve para comparar los niveles de glucosa entre tratamientos. 
# En este caso, ninguna de las comparaciones mostró significación estadística, lo que nos confirma que ninguno de los fármacos modifican notoriamente los niveles de glucosa en comparación al placebo.

