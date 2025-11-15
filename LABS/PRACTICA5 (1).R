#############################################################################
#
# PRACTICA R
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data)
head(data)
tail(data)

# Hacemos un primer histograma para explorar los datos
hist(data, col = "gray", main="GSE5583 - Histogram")

# Transformamos los datos con un logaritmo 
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
# Aunque para los datos brutos no se realiza ninguna transformación normalmente, al hacerla nos ayuda a visualizar mejor los datos.
data2 = log2(data)
hist(data2, col = "gray", main="GSE5583 (log2) - Histogram")


# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
# ¿Qué es un boxplot?
# Un boxplot se trata de una representación gráfica que resume la distribución de un conjunto de datos numéricos.
# La mediana es la línea que vemos dentro de la caja.
# Los cuartiles son los propios bordes de la caja(Q1=25%, Q3=75%).
# El rango intercuartílico es la altura de la caja(IQR = Q3-Q1).
# Los bigotes van desde los cuartiles hasta 1.5 × IQR o hasta el valor más extremo dentro de ese rango.
# Los outliers son los puntos que aparecen fuera de los bigotes.
# En nuestro caso vemos que nos crea un boxplot para cada columna, con el título GSE5583-boxplots. 
# Las=2 sirve para que aparezcan los nombres en vertical
boxplot(data2, col=c("blue", "blue", "blue",
	"orange", "orange", "orange"),
	main="GSE5583 - boxplots", las=2)
boxplot(data, col=c("blue", "blue", "blue",
                     "orange", "orange", "orange"),
        main="GSE5583 - boxplots", las=2)
#Lo tranformamos para que los datos puedan visualizarse bien. 	


# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación ç
# de los valores de expresión. ¿Es correcta la separación?
# El clustering se hace para comprobar que hemos conseguido las muestras bien con patrones simples. 
# Nos dice si hemos cogido muestras con genes similares, agrupa WT en un cluster y KO en otro. En este caso como lo divide en dos grupos, podríamos decir que es correcta la separación.
hc = hclust(as.dist(1-cor(data2)))
plot(hc, main="GSE5583 - Hierarchical Clustering")


#######################################
# Análisis de Expresión Diferencial 
#######################################
# Ahora tenemos un listado de genes con expresiones de KO y WT. Lo que queremos ver es qué genes son suficientemente distintos, para así, saber cuales tienen un papel importante en la expresión de ese gen.
# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
# Tenemos WT con sus tres variables y KO con otras tres.
wt <- data[,1:3]
ko <- data[,4:6]
class(wt)
head(wt)
head(ko)
# De esta forma, lo que conseguimos es crear una matriz de datos

# Calcula las medias de las muestras para cada condición. Usa apply
wt.mean = apply(wt, 1, mean)
ko.mean = apply(ko, 1, mean)
head(wt.mean)
head(ko.mean)
head(data)

# ¿Cuál es la media más alta?
limit = max(wt.mean, ko.mean)
limit
# La media más alta es 37460.5
# No la calculamos por ninguna razón en especial, lo hacemos con la función de max (o min) y seguido entre parentesis los objetos deseados.
# En el eje X hemos puesto los valores del WT y en el eje Y los de KO. 
# Vemos que aquí hay una asociación.

# Ahora hacemos un scatter plot (gráfico de dispersión)
plot(ko.mean ~ wt.mean, xlab = "WT", ylab = "KO",
	main = "GSE5583 - Scatter", xlim = c(0, limit), ylim = c(0, limit))
#el xlab y el ylab nos personaliza el texto de los ejes.

# Añadir una línea diagonal con abline
abline(0, 1, col = "red")
# También podríamos añadir nosotros más líneas. Tenemos que darle la pendiente y la ordenada.
# Si la correlación fuese perfecta todos los puntitos estarían pegados a la línea.

# ¿Eres capaz de añadirle un grid?
# Sí, nos añade un fondo de celdas.
grid()
#abline(a, b): línea de pendiente b y ordenada en el origen a
#abline(h=y): línea horizontal
#abline(v=x): línea vertical
abline(1, 2, col = "red")     # línea y = 2x + 1
abline(h = 2, col = "green")  # línea y = 2
abline(v = 3, col = "violet") # línea x = 3
# Hay que eliminar el # de delante de las 3 ultimas para que nos aparezcan las rectas.

# Calculamos la diferencia entre las medias de las condiciones
diff.mean = wt.mean - ko.mean
# Da igual lo que pongamos delante o detrás mientras lo sepamos nosotros.

# Hacemos un histograma de las diferencias de medias
hist(diff.mean, col = "gray")

# Calculamos la significancia estadística con un t-test.
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
# ¿Cuántas valores tiene cada muestra?
pvalue = NULL 
tstat = NULL 
for(i in 1 : nrow(data)) { #Para cada gen
	x = wt[i,] # gene wt número i
	y = ko[i,] # gene ko número i
	
	# Hacemos el test
	t = t.test(x, y)
	
	# Añadimos el p-value a la lista
	pvalue[i] = t$p.value
	# Añadimos las estadísticas a la lista
	tstat[i] = t$statistic
}

head(pvalue)
# Hemos creado un vector de p values vacío, otro de estadísticas vacío y con un bucle hemos hecho para cada gen de nuestra tabla calcúlame para cada valor, los tres valores del WT y los KO.

# Ahora comprobamos que hemos hecho TODOS los cálculos
length(pvalue)
# Como no hemos evaluado la normalidad no sabemos si tiene una distribución normal o no.
# Misma longitud que genes que tenga. En este caso a pesar de que la t student nos sirve de ejemplo, no estaría del todo bien.

# Hacemos un histograma de los p-values.
# ¿Qué pasa si le ponemos con una transformación de -log10?
# Transformamos los datos, de esta forma nos sirve para visualizar mejor los datos.
hist(pvalue,col="gray")
hist(-log10(pvalue), col = "gray")

# El primer histograma es de p values brutos y el segundo de transformados.

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano")

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
diff.mean_cutoff = 2
pvalue_cutoff = 0.01
abline(v = diff.mean_cutoff, col = "blue", lwd = 3)

#abline(v = -diff.mean_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)

# El volcano plot es muy usado en transcriptómica o al analizar expresión de genes. 
# En este caso, en el eje X tenemos diferencia de medias, en el eje de la Y los p values. Las bolitas negras son los genes.
# Los genes que tienen la misma expresión en WT y KO estarán cerca del 0.
# Tenemos el p value y una línea en el 2, que es l quivalente al 0,01 a la transformación. Esto significa que los genes que pasen el corte se situarán por encima de la línea.

# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean = abs(diff.mean) >= diff.mean_cutoff
dim(data[filter_by_diff.mean, ])

# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data[filter_by_pvalue, ])

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?
# Lo cumplen 426 genes.
filter_combined = filter_by_diff.mean & filter_by_pvalue
filtered = data[filter_combined,]
dim(filtered)
head(filtered)

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #2")
points (diff.mean[filter_combined], -log10(pvalue[filter_combined]),col = "red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés?
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #3")
points (diff.mean[filter_combined & diff.mean < 0],
	-log10(pvalue[filter_combined & diff.mean < 0]), col = "red")
points (diff.mean[filter_combined & diff.mean > 0],
	-log10(pvalue[filter_combined & diff.mean > 0]), col = "blue")

# El filter tiene todos los genes que estan expresados de forma diferencial. En la consola nos debería aparecer cuantos han pasado el corte (426).
# Los genes rojos estan sobreexpresados y los azules son reprimidos en WT.
# En este caso, los reprimidos aparecen como valores positivos.

# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
# El mapa de calor es un ejemplo de clustering. Nos clasifica los datos. 
# En este caso las tres replicas de WT en un blowue y las de KO en un bloque. Los hacemos unicamente con los genes diferencialmente expresados. Si un gen esta sobreexpresado en WT en KO estara reprimido y del reves.

# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,labRow=FALSE)

heatmap(filtered)
# No ordenamos por genes. Es un ejemplo de clustering. Le metes las columnas y filas y este, las clasifica en su sitio. 
# Los árboles de arriba son el clustering.

# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)
library(RColorBrewer)

# Hacemos nuestro heatmap
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row")

# Lo guardamos en un archivo PDF
pdf ("GSE5583_DE_Heatmap.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row",labRow=FALSE)
dev.off()
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,col = redgreen(75), scale = "row",labRow=FALSE)

# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table (filtered, "GSE5583_DE.txt", sep = "\t",quote = FALSE)
