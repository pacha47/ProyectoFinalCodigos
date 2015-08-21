#ifndef INC_MALLA_H
#define INC_MALLA_H

#include <vector>
#include <iostream>

using namespace std;

typedef vector<vector<double> > mdouble;
typedef vector<vector<int> > mint;


#include "vec.h"
#include "mat.h"

#include "Nodo.h"
#include "Arista.h"
#include "Elemento.h"

class Malla{
	
private:
	
	vector<Nodo> nodos;
	vector<Arista> aristas;
	vector<Elemento> elementos;
	
	mat M;
	vec F, P;
	
	double max,min,h;
	bool flag;
public:
	
	/// 	FUNCIONES N-S FRACTIONA STEP
	void operador_CD( double dt );
	void operador_CD_AB(double dt);
	void operador_P1( double dt );
	double operador_P2( double dt );
	void primeraIte();
	
	void addvel();
	void defVel();
	
	double iterar(double dt);
	
	/// CONSTRUCTORES
	Malla();
	double makeMalla(mdouble nods, mint e);
	double makeMalla(mdouble nods, mint e, mdouble ncond);
	
	void write(FILE *fs);

};


#endif
