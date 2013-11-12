CXX=g++

LIBS=

raytracer: raytracer.cpp mat4.o vec4.o material.o object.o mesh.o sphere.o pixelmap.o raytracer.h
	${CXX} -o raytracer raytracer.cpp mat4.o vec4.o material.o object.o mesh.o sphere.o pixelmap.o raytracer.h ${LIBS}

mat4.o: mat4.cpp mat4.h
	${CXX} -c mat4.cpp mat4.h ${LIBS}

vec4.o: vec4.cpp vec4.h mat4.h
	${CXX} -c vec4.cpp vec4.h ${LIBS}

material.o: material.cpp material.h vec4.h
	${CXX} -c material.cpp material.h ${LIBS}

object.o: object.cpp object.h vec4.h material.h
	${CXX} -c object.cpp object.h ${LIBS}

mesh.o: mesh.cpp mesh.h mat4.h vec4.h material.h object.h
	${CXX} -c mesh.cpp mesh.h object.h ${LIBS}

sphere.o: sphere.cpp sphere.h mat4.h vec4.h material.h object.h
	${CXX} -c sphere.cpp sphere.h object.h ${LIBS}

pixelmap.o: pixelmap.cpp pixelmap.h vec4.h
	${CXX} -c pixelmap.cpp pixelmap.h ${LIBS}

clean:
	rm -f *.o *.gch raytracer