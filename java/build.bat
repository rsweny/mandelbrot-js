del /Q /S classes\fractal
del /Q /S classes\math

javac -d classes -classpath classes -sourcepath src src/fractal/Mandelbrot3D.java
javac -d classes -classpath classes -sourcepath src src/fractal/Mandelbrot.java
javac -d classes -classpath classes -sourcepath src src/fractal/Buddhabrot.java
javac -d classes -classpath classes -sourcepath src src/fractal/Newton.java
javac -d classes -classpath classes -sourcepath src src/fractal/Newtonbrot.java
javac -d classes -classpath classes -sourcepath src src/fractal/Flame.java
javac -d classes -classpath classes -sourcepath src src/fractal/Taylor.java

jar cvf Mandelbrot3D.jar -C classes fractal -C classes math
jarsigner Mandelbrot3D.jar fractalkey
move Mandelbrot3D.jar classes

time /t
