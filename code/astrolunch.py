#CODIGO DE ASTROLUNCH (Nicolas Rocha)

#Importa las librerias necesarias para analizar la imagen
import numpy as np
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

#Guarda los valores de la imagen en una matriz
hdulist  = fits.open("../data/bbso_tio_pcosr_20130902_162238.fts")
image_data = hdulist[0].data
print ""
print "La imagen ha sido cargada"

#Imprime informacion de la imagen que se cargo
print ""
print "Informacion de la imagen:"
print hdulist.info()
print(image_data.shape)

#Eliminamos los bordes de la imagen
#TODO Revisar valores en caso que el borde cambie de grosor
new_image_data = 1.0*image_data[ 100:-100 , 100:-100 ]


#Muestra la imagen
plt.imshow( new_image_data , cmap='gray')
plt.colorbar()
plt.show()

#Definimos la funcion que crea la matriz de matrices para ejecutar el algoritmo
def algoritmo( matriz ):
	derivada = np.zeros( ( 1843 , 1843 ) , dtype = object )
	for x in range( 1 , 1842 ):
		for y in range( 1 , 1842 ):
			pixel = np.zeros( (2 , 2) )
			pixel[0][0] = matriz[x+1][y] - 2.0*matriz[x][y] + matriz[x-1][y]
			pixel[1][0] = ( matriz[x+1][y+1] - matriz[x+1][y-1] - matriz[x-1][y+1] + matriz[x-1][y-1] ) / 4
			pixel[0][1] = ( matriz[x+1][y+1] - matriz[x+1][y-1] - matriz[x-1][y+1] + matriz[x-1][y-1] ) / 4
			pixel[1][1] = matriz[x][y+1] - 2.0*matriz[x][y] + matriz[x][y-1]
			derivada[x][y] = pixel
	return derivada

#Calculamos la matriz asociada al algoritmo
derivada = algoritmo(new_image_data)

#np.linalg.eig( matriz )
#devuelve 1x2 ==> [1][0] Autovalores y [1][1] Autovectores
print derivada[100][100]
print np.linalg.eig( derivada[100][100] )

#Definimos la funcion que calcula y devuelve los autovalores y autovectores
def auto( matriz ):
	autovalores = np.zeros( ( 1843 , 1843 ) , dtype = object )
	for x in range( 0 , 1843 ):
		for y in range( 0 , 1843 ):
			eigenstuff = np.zeros( ( 3 , 3 ) )
			#En esta linea esta el error
			eigenstuff = np.linalg.eig( matriz[x][y] )
			autovalor  = eigenstuff[1][0]
			autovalores[x][y] = autovalor
	return autovalores


print auto( derivada )

