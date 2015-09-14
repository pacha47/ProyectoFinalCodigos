#ifndef OPENGL_H
#define OPENGL_H



#include <GL/glut.h>

int w=640,h=480, el=1;
double escala = 400.0;
extern Malla m;

void regen() {
	
	// matriz de proyeccion
	glMatrixMode(GL_PROJECTION);  glLoadIdentity();
	double w0=(double)w/2.0/escala,h0=(double)h/2.0/escala; // semiancho y semialto en el target
	glOrtho(-w0,w0,-h0,h0,0,1);
	glTranslated(-w0+.1,-h0+.1,0);
	glMatrixMode(GL_MODELVIEW); glLoadIdentity(); // matriz del modelo->view
	glutPostRedisplay(); // avisa que hay que redibujar
}

void reshape_cb (int w, int h) {
	if (w==0||h==0) return;
	glViewport(0,0,w,h);
	regen();
	glLoadIdentity ();
}

void display_cb() {
//	glClear(GL_COLOR_BUFFER_BIT);
//	glColor3f(1,1,0); glLineWidth(3);
//	glBegin(GL_LINE_STRIP);
//		glVertex2i(130,060); glVertex2i( 50,060);
//		glVertex2i(130,150); glVertex2i( 50,150);
//	glEnd();
//	glBegin(GL_LINES);
//		glVertex2i( 70,100); glVertex2i(110,100);
//		glVertex2i(150,100); glVertex2i(230,100);
//		glVertex2i(190,140); glVertex2i(190,070);
//		glVertex2i(250,100); glVertex2i(330,100);
//		glVertex2i(290,140); glVertex2i(290,070);
//	glEnd();
	
	m.dibele();
//	m.dibquick();
	
	glutSwapBuffers();
}

void Keyboard_cb(unsigned char key,int x=0,int y=0) {
	glClear(GL_COLOR_BUFFER_BIT);
	switch(key){
	case '+':
		escala+=2;
		cout<<escala<<endl;
		break;
	case '-':
		escala-=2;
		break;
	case 'a':
		cout<<"dib "<<el<<endl;
		m.dibquick(el++);
		break;
	case 'A':
		cout<<"dib "<<el<<endl;
		m.dibquick(el--);
		break;
	case 27: // escape => exit2
		exit(EXIT_SUCCESS);
		break;
	}
	regen();
	glutPostRedisplay();
	
}


void initialize() {
	glutInitDisplayMode (GLUT_RGBA|GLUT_DOUBLE);
	glutInitWindowSize (w,h);
	glutInitWindowPosition (100,100);
	glutCreateWindow ("Ventana OpenGL");
	glutDisplayFunc (display_cb);
	glutReshapeFunc (reshape_cb);
	glutKeyboardFunc(Keyboard_cb);
	glClearColor(1.f,1.f,1.f,1.f);
}



#endif
