#ifndef INC_MALLA_H
#define INC_MALLA_H

#include <vector>
#include <iostream>
#include <fstream>


using namespace std;

typedef vector<vector<double> > mdouble;
typedef vector<vector<int> > mint;

#include "Nodo.h"
#include "Arista.h"
#include "Elemento.h"
#include "vec.h"
#include "mat.h"

class Malla{
	
private:
	
	vector<Nodo> nodos;
	vector<Arista> aristas;
	vector<Elemento> elementos;
	double max,min,h, maxp,minp;
	mat K, U;
	vec F, B;
public:
	double error;
	
	Malla(){};
	
	double makeMalla(mdouble nods, mint e, mdouble cond);
	
	/// IMPLEMENTAR LO DE NS ACA
	void iterar();
	void corregir();
	void asignar();
	
	/// ESCRITURA DEL ARCHIVO DE SALIDA
	void write(FILE *fs);
	
	void findVel(double x, double y);
	
	void show(){
		int n = elementos.size();
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				printf("%.5f ",K(i,j));
			}
			printf("\n");
		}
		
		
	}
};


#endif
