###1. Cargar los datos y examinarlos
data <- read.table (file='datos-trabajoR.txt', header=TRUE)
head(data)

summary(data)

##Filas y columnas: 4 columnas por lo que hay 4 variables
dim(data)
str(data)
tail(data)

length(unique(data$Tratamiento))
##Hay 6 tratamientos


###2. Crear un boxplot
colorcuatro <- c("red","blue", "green", "purple")
colortres <- c("red","blue","green")
boxplot(data$Wildtype, data$Sequia, data$ExcesoRiego, col=colortres, main= "Tratamiento- Boxplot", las=2)


###3. Graficos de dispersion
colores <- c("black", "red", "green", "cyan", "blue", "pink")
plot(data$Wildtype, data$Sequia, col=colores, xlab= "Wildtype", ylab="Sequia", main="Dispersion- Sequia vs Wildtype")
plot(data$Wildtype, data$ExcesoRiego, col=colores, xlab= "Wildtype", ylab="ExcesoRiego", main="Dispersion- ExcesoRiego vs Wildtype")


###4. Leyenda
legend("bottomright", legend= paste("Tto", unique(data$Tratamiento)), col=colores, pch=15, cex=0.8, title= "Tratamientos")


###5. Histograma para cada variable
str(data)
hist (data$Tratamiento, col= colores , main="Histograma- Tratamiento", xlab= "Tratamiento")
hist (data$Wildtype, col= colores , main="Histograma- WildType", xlab= "Wildtype")
hist (data$Sequia, col= colores , main="Histograma- Sequia", xlab= "Sequia")
hist (data$ExcesoRiego, col= colores, main="Histograma- ExcesoRiego", xlab= "ExcesoRiego")


###6. Factor en la columna tratamiento
tratamiento_factor <- factor(data$Tratamiento)
print(tratamiento_factor)


###7. Media y desviación estándar
##Con aggregate()
names(data)

#Wildtype
media_wildtype <- aggregate(Wildtype ~ Tratamiento, data,  FUN = mean)
print(media_wildtype)
desviacion_estandar_wt <- aggregate(Wildtype ~ Tratamiento, data, FUN= sd)
print(desviacion_estandar_wt)

#Sequia
media_sequia <- aggregate(Sequia ~ Tratamiento, data,  FUN = mean)
print(media_sequia)
desviacion_estandar_sequia <- aggregate(Sequia ~ Tratamiento, data, FUN= sd)
print(desviacion_estandar_sequia)

#Exceso de Riego
media_excesoriego <- aggregate(ExcesoRiego ~ Tratamiento, data,  FUN = mean)
print(media_excesoriego)
desviacion_estandar_excesoriego <- aggregate(ExcesoRiego ~ Tratamiento, data, FUN= sd)
print(desviacion_estandar_excesoriego)


##Con tapply()
#Wildtype
media_tap_wt <- tapply(data$Wildtype, data$Tratamiento, FUN= mean)
print(media_tap_wt)
desviacion_tap_wt <- tapply(data$Wildtype, data$Tratamiento, FUN= sd)
print(desviacion_tap_wt)

#Sequia
media_tap_sequia <- tapply(data$Sequia, data$Tratamiento, FUN= mean)
print(media_tap_sequia)
desviacion_tap_sequia <- tapply(data$Sequia, data$Tratamiento, FUN= sd)
print(desviacion_tap_sequia)

#ExcesoRiego
media_tap_er <- tapply(data$ExcesoRiego, data$Tratamiento, FUN= mean)
print(media_tap_er)
desviacion_tap_er <- tapply(data$ExcesoRiego, data$Tratamiento, FUN= sd)
print(desviacion_tap_er)


###8.Table(). Cuantos elementos tiene cada tratamiento
elementos_tratamiento <- table(data$Tratamiento)
print(elementos_tratamiento)


###9. Sacar los datos del tratamiento 1 y 4 y guardarlos en una variable
tratamiento_uno <- subset(data, Tratamiento ==1)
print(tratamiento_uno)
tratamiento_cuatro <- subset (data, Tratamiento==4)
print(tratamiento_cuatro)

tratamiento_dos <- subset(data, Tratamiento ==2)
print(tratamiento_dos)
tratamiento_tres <- subset (data, Tratamiento==3)
print(tratamiento_tres)
tratamiento_cinco <- subset (data, Tratamiento==5)
print(tratamiento_cinco)


###10. Queremos comprobar que hay diferencias significativas para el tratamiento 1 y el tratamiento 5 entre Wildtype y
 Sequia, y entre Wildtype y ExcesoRiego. Primero, necesitaríamos comprobar si los datos se distribuyen de forma
 normal. En función de los resultados de la prueba de normalidad, ¿qué test usarías para cada comparativa? ¿Puedes
 comparar también Sequia con ExcesoRiego en ambos tratamientos? En general, asumimos que las muestras son
 independientes, pero ¿son sus varianzas iguales? Actúa de acuerdo con tus resultados. 
##Para comprobar la normalidad, Shapiro test
#Para tratamiento uno
shapiro.test (tratamiento_uno$Wildtype)
shapiro.test (tratamiento_uno$Sequia)
shapiro.test (tratamiento_uno$ExcesoRiego)

#Para tratamiento 5
shapiro.test (tratamiento_cinco$Wildtype)
shapiro.test (tratamiento_cinco$Sequia)
shapiro.test (tratamiento_cinco$ExcesoRiego)

##Según nuestros resultados, para considerar distribuición normal el p-value debe ser mayor a 0.05

##Utilizaria el t.test() para tratamiento uno (p-value > 0.05) y utilizaría Wilcox cuando p-value < 0.05)
#Tratamiento uno
t_test_tto1_wt_sq <- t.test (tratamiento_uno$Wildtype, tratamiento_uno$Sequia, var.equal=TRUE)
print(t_test_tto1_wt_sq)
t_test_tto1_wt_er <- t.test (tratamiento_uno$Wildtype, tratamiento_uno$ExcesoRiego, var.equal=TRUE)
print(t_test_tto1_wt_er)
t_test_tto1_sq_er <- t.test (tratamiento_uno$Sequia, tratamiento_uno$ExcesoRiego, var.equial=TRUE)
print(t_test_tto1_sq_er)

#Tratamiento cinco
wilcox_tto5_wt_sq <- wilcox.test (tratamiento_cinco$Wildtype, tratamiento_cinco$Sequia, alternative="g")
print(wilcox_tto5_wt_sq)
wilcox_tto5_wt_er <- wilcox.test (tratamiento_cinco$Wildtype, tratamiento_cinco$ExcesoRiego, alternative="g")
print(wilcox_tto5_wt_er)
wilcox_tto5_sq_er <- wilcox.test (tratamiento_cinco$Sequia, tratamiento_cinco$ExcesoRiego, alternative="g")
print(wilcox_tto5_sq_er)

##¿Son sus varianzas iguales?
varianza_wt <- var.test(tratamiento_uno$Wildtype, tratamiento_cinco$Wildtype)
print(varianza_wt)
#La varianza del Wildtype es muy diferente debido a que el valor F está muy alejado del 1

varianza_sequia <- var.test(tratamiento_uno$Sequia, tratamiento_cinco$Sequia)
print(varianza_sequia)
#La varianza en sequía es muy diferente debido a que el valor F está muy alejado del 1

varianza_excesoriego <- var.test(tratamiento_uno$ExcesoRiego, tratamiento_cinco$ExcesoRiego)
print(varianza_excesoriego)
#La varianza en Exceso de Riego es muy diferente ya que el valor F es muy alejado de 1



### 11. Realiza un ANOVA para comparar el tratamiento 1 en las tres condiciones. Pista: primero separa los valores de
 tratamiento1 en Wildtype, Sequia y ExcesoRiego en variables separadas. Luego fíjate en el archivo “datos
anova.txt” y trata de colocar los datos de esa forma en una tabla. Por último, ejecuta el test.
data
wildtype_tto_uno <- c(2.2, 2.3, 2.4, 4.5, 5.4, 5.4, 5.4, 4.0, 4.0, 4.4)
sequia_tto_uno <- c(0.20, 0.30, 0.40, 0.50, 0.20, 0.40, 0.20, 1.00, 0.90, 0.80)
excesoriego_tto_uno <- c(6.0, 6.2, 6.4, 6.5, 5.8, 6.1, 6.4, 5.6, 5.4, 5.3)

valores <- c(wildtype_tto_uno, sequia_tto_uno, excesoriego_tto_uno)
tratamiento <- as.factor(rep(c("Wildtype", "Sequia", "Exceso Riego"), each=length(wildtype_tto_uno)))

#boxplot(valores~tratamiento, col=c("blue","yellow","pink"), ylab= "Valores medios")

datos_anova <- read.table("datos-trabajoR.txt",sep=" ",header=TRUE)
anova <- aov(valores ~ tratamiento, data=datos_anova)
summary(anova)


###1. Cargar los datos y examinarlos
data <- read.table (file='datos-trabajoR.txt', header=TRUE)
head(data)

summary(data)

##Filas y columnas: 4 columnas por lo que hay 4 variables
dim(data)
str(data)
tail(data)

length(unique(data$Tratamiento))
##Hay 6 tratamientos


###2. Crear un boxplot
colorcuatro <- c("red","blue", "green", "purple")
colortres <- c("red","blue","green")
boxplot(data$Wildtype, data$Sequia, data$ExcesoRiego, col=colortres, main= "Tratamiento- Boxplot", las=2)


###3. Graficos de dispersion
colores <- c("black", "red", "green", "cyan", "blue", "pink")
plot(data$Wildtype, data$Sequia, col=colores, xlab= "Wildtype", ylab="Sequia", main="Dispersion- Sequia vs Wildtype")
plot(data$Wildtype, data$ExcesoRiego, col=colores, xlab= "Wildtype", ylab="ExcesoRiego", main="Dispersion- ExcesoRiego vs Wildtype")


###4. Leyenda
legend("bottomright", legend= paste("Tto", unique(data$Tratamiento)), col=colores, pch=15, cex=0.8, title= "Tratamientos")


###5. Histograma para cada variable
str(data)
hist (data$Tratamiento, col= colores , main="Histograma- Tratamiento", xlab= "Tratamiento")
hist (data$Wildtype, col= colores , main="Histograma- WildType", xlab= "Wildtype")
hist (data$Sequia, col= colores , main="Histograma- Sequia", xlab= "Sequia")
hist (data$ExcesoRiego, col= colores, main="Histograma- ExcesoRiego", xlab= "ExcesoRiego")


###6. Factor en la columna tratamiento
tratamiento_factor <- factor(data$Tratamiento)
print(tratamiento_factor)


###7. Media y desviación estándar
##Con aggregate()
names(data)

#Wildtype
media_wildtype <- aggregate(Wildtype ~ Tratamiento, data,  FUN = mean)
print(media_wildtype)
desviacion_estandar_wt <- aggregate(Wildtype ~ Tratamiento, data, FUN= sd)
print(desviacion_estandar_wt)

#Sequia
media_sequia <- aggregate(Sequia ~ Tratamiento, data,  FUN = mean)
print(media_sequia)
desviacion_estandar_sequia <- aggregate(Sequia ~ Tratamiento, data, FUN= sd)
print(desviacion_estandar_sequia)

#Exceso de Riego
media_excesoriego <- aggregate(ExcesoRiego ~ Tratamiento, data,  FUN = mean)
print(media_excesoriego)
desviacion_estandar_excesoriego <- aggregate(ExcesoRiego ~ Tratamiento, data, FUN= sd)
print(desviacion_estandar_excesoriego)


##Con tapply()
#Wildtype
media_tap_wt <- tapply(data$Wildtype, data$Tratamiento, FUN= mean)
print(media_tap_wt)
desviacion_tap_wt <- tapply(data$Wildtype, data$Tratamiento, FUN= sd)
print(desviacion_tap_wt)

#Sequia
media_tap_sequia <- tapply(data$Sequia, data$Tratamiento, FUN= mean)
print(media_tap_sequia)
desviacion_tap_sequia <- tapply(data$Sequia, data$Tratamiento, FUN= sd)
print(desviacion_tap_sequia)

#ExcesoRiego
media_tap_er <- tapply(data$ExcesoRiego, data$Tratamiento, FUN= mean)
print(media_tap_er)
desviacion_tap_er <- tapply(data$ExcesoRiego, data$Tratamiento, FUN= sd)
print(desviacion_tap_er)


###8.Table(). Cuantos elementos tiene cada tratamiento
elementos_tratamiento <- table(data$Tratamiento)
print(elementos_tratamiento)


###9. Sacar los datos del tratamiento 1 y 4 y guardarlos en una variable
tratamiento_uno <- subset(data, Tratamiento ==1)
print(tratamiento_uno)
tratamiento_cuatro <- subset (data, Tratamiento==4)
print(tratamiento_cuatro)

tratamiento_dos <- subset(data, Tratamiento ==2)
print(tratamiento_dos)
tratamiento_tres <- subset (data, Tratamiento==3)
print(tratamiento_tres)
tratamiento_cinco <- subset (data, Tratamiento==5)
print(tratamiento_cinco)


###10. Queremos comprobar que hay diferencias significativas para el tratamiento 1 y el tratamiento 5 entre Wildtype y
 Sequia, y entre Wildtype y ExcesoRiego. Primero, necesitaríamos comprobar si los datos se distribuyen de forma
 normal. En función de los resultados de la prueba de normalidad, ¿qué test usarías para cada comparativa? ¿Puedes
 comparar también Sequia con ExcesoRiego en ambos tratamientos? En general, asumimos que las muestras son
 independientes, pero ¿son sus varianzas iguales? Actúa de acuerdo con tus resultados. 
##Para comprobar la normalidad, Shapiro test
#Para tratamiento uno
shapiro.test (tratamiento_uno$Wildtype)
shapiro.test (tratamiento_uno$Sequia)
shapiro.test (tratamiento_uno$ExcesoRiego)

#Para tratamiento 5
shapiro.test (tratamiento_cinco$Wildtype)
shapiro.test (tratamiento_cinco$Sequia)
shapiro.test (tratamiento_cinco$ExcesoRiego)

##Según nuestros resultados, para considerar distribuición normal el p-value debe ser mayor a 0.05

##Utilizaria el t.test() para tratamiento uno (p-value > 0.05) y utilizaría Wilcox cuando p-value < 0.05)
#Tratamiento uno
t_test_tto1_wt_sq <- t.test (tratamiento_uno$Wildtype, tratamiento_uno$Sequia, var.equal=TRUE)
print(t_test_tto1_wt_sq)
t_test_tto1_wt_er <- t.test (tratamiento_uno$Wildtype, tratamiento_uno$ExcesoRiego, var.equal=TRUE)
print(t_test_tto1_wt_er)
t_test_tto1_sq_er <- t.test (tratamiento_uno$Sequia, tratamiento_uno$ExcesoRiego, var.equial=TRUE)
print(t_test_tto1_sq_er)

#Tratamiento cinco
wilcox_tto5_wt_sq <- wilcox.test (tratamiento_cinco$Wildtype, tratamiento_cinco$Sequia, alternative="g")
print(wilcox_tto5_wt_sq)
wilcox_tto5_wt_er <- wilcox.test (tratamiento_cinco$Wildtype, tratamiento_cinco$ExcesoRiego, alternative="g")
print(wilcox_tto5_wt_er)
wilcox_tto5_sq_er <- wilcox.test (tratamiento_cinco$Sequia, tratamiento_cinco$ExcesoRiego, alternative="g")
print(wilcox_tto5_sq_er)

##¿Son sus varianzas iguales?
varianza_wt <- var.test(tratamiento_uno$Wildtype, tratamiento_cinco$Wildtype)
print(varianza_wt)
#La varianza del Wildtype es muy diferente debido a que el valor F está muy alejado del 1

varianza_sequia <- var.test(tratamiento_uno$Sequia, tratamiento_cinco$Sequia)
print(varianza_sequia)
#La varianza en sequía es muy diferente debido a que el valor F está muy alejado del 1

varianza_excesoriego <- var.test(tratamiento_uno$ExcesoRiego, tratamiento_cinco$ExcesoRiego)
print(varianza_excesoriego)
#La varianza en Exceso de Riego es muy diferente ya que el valor F es muy alejado de 1



### 11. Realiza un ANOVA para comparar el tratamiento 1 en las tres condiciones. Pista: primero separa los valores de
 tratamiento1 en Wildtype, Sequia y ExcesoRiego en variables separadas. Luego fíjate en el archivo “datos
anova.txt” y trata de colocar los datos de esa forma en una tabla. Por último, ejecuta el test.
data

tto_uno_Wildtype <- tratamiento_uno$Wildtype
print(tto_uno_Wildtype)
tto_uno_Sequia <- tratamiento_uno$Sequia
print(tto_uno_Sequia)
tto_uno_ExcesoRiego <- tratamiento_uno$ExcesoRiego
print(tto_uno_ExcesoRiego)

datosAnova <- data.frame("Condicion" = c("W","W","W","W","W","W","W","W","W","W","S","S","S","S","S","S","S","S","S","S","E","E","E","E","E","E","E","E","E","E"),"Valor"=c(tto_uno_Wildtype, tto_uno_Sequia, tto_uno_ExcesoRiego))
anova <- aov(Valor ~ Condicion, data=datosAnova)
summary(anova)
