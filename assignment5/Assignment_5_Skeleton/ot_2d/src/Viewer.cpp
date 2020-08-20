/*! \file viewer.cpp
*   \brief 3D mesh viewer
*   \author David Gu
*   \date documented on 10/12/2010
*
*	3D triangle mesh viewer
*
*  Visualize a textured mesh. The input is a mesh with uv coordinates and texture bmp image. Command <P>
*  <B> ViewerConverter.exe sophie.uv.obj sophie.bmp <B><P>
*  <B> ViewerConverter.exe sophie.uv.m sophie.bmp <B><P>
*  <P>
* <table border="0" cellspacing="0">
* <tr>
* <td><IMG SRC="Viewer/Sophie_geom.png" alt="Input Mesh" width="320" height="320"></td>
* <td><IMG SRC="Viewer/Sophie.png" alt="Input Mesh" width="320" height="320"></td>
* <td><IMG SRC="Viewer/Sophie_texture.png" alt="Input Mesh" width="320" height="320"></td>
* </tr>
* <tr>
* <td><IMG SRC="Viewer/alex_geom.png" alt="Input Mesh" width="320" height="320"></td>
* <td><IMG SRC="Viewer/alex.png" alt="Input Mesh" width="320" height="320"></td>
* <td><IMG SRC="Viewer/alex_texture.png" alt="Input Mesh" width="320" height="320"></td>
* </tr>
* <tr>
* <td><IMG SRC="Viewer/susan_geom.png" alt="Input Mesh" width="320" height="320"></td>
* <td><IMG SRC="Viewer/susan.png" alt="Input Mesh" width="320" height="320"></td>
* <td><IMG SRC="Viewer/susan_texture.png" alt="Input Mesh" width="320" height="320"></td>
* </tr>
* </table>
*

*/

/*******************************************************************************
 *      Viewer - 3D triangle mesh viewer
 *
 *       Copyright (c) CCGL
 *
 *    Purpose:
 *       Display 3D triangle meshes
 *
 *       David Gu June 27, 2008
 *
 *******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef MAC_OS
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif // MAC_OS

#include "Operator/Operator.h"
#include "viewer/Arcball.h" /*  Arc Ball  Interface         */
#include "CDomainOptimalTransport.h"

using namespace MeshLib;

/* window width and height */
int win_width, win_height;
int gButton;
int startx, starty;
int texture_mode = 0;
int display_mode = 0;
int show_cell = 0;
bool light_edit = false;
// texture name
GLuint texture_name;

/* rotation quaternion and translation vector for the object */
CQrot ObjRot(0, 0, 1, 0);
CPoint ObjTrans(0, 0, 0);
CQrot LightRot(1, 0, 0, 0);

CPoint LightPosition;

/* global mesh */
COMTMesh mesh;
CDomainOptimalTransport* pOT = NULL;

/* arcball object */
CArcball arcball;

/*! setup the object, transform from the world to the object coordinate system */
void setupObject(void)
{
    double rot[16];

    glTranslated(ObjTrans[0], ObjTrans[1], ObjTrans[2]);
    ObjRot.convert(rot);
    glMultMatrixd((GLdouble*) rot);
}

/*! the eye is always fixed at world z = +5 */

void setupEye(void)
{
    glLoadIdentity();
    gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);
}

/*! setup light */
void setupLight()
{
    CQrot v(0, 0, 0, 1);
    CQrot CL(LightRot.m_w, -LightRot.m_x, -LightRot.m_y, -LightRot.m_z);
    CL = LightRot * v * CL;

    LightPosition = CPoint(CL.m_x, CL.m_y, CL.m_z);

    GLfloat lightOnePosition[4] = {(GLfloat) LightPosition[0], (GLfloat) LightPosition[1], (GLfloat) LightPosition[2], 0.0f};
    glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);
}

void draw_voronoi_cell(COMTMesh* pMesh)
{
    glDisable(GL_LIGHTING);

    if (show_cell)
    {
        glColor3f(1, 0, 0);
        glLineWidth(1.0);
        glBegin(GL_LINES);
        for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
        {
            COMTMesh::CVertex* pV = *viter;

            CPolygon3D& pg = pV->dual_cell();

            for (size_t i = 0; i < pg.edges().size(); i++)
            {
                CSegment3D& s = pg.edges()[i];
                CPoint p = s.start();
                CPoint q = s.end();
                glColor3f(1, 0, 0);
                glVertex3d(p[0], p[1], -2e-3);
                glVertex3d(q[0], q[1], -2e-3);
            }
        }
        glEnd();
        glLineWidth(1.0);
    }
    glFrontFace(GL_CW);
    glEnable(GL_LIGHTING);
    glBegin(GL_TRIANGLES);
    for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        COMTMesh::CVertex* pV = *viter;

        CPolygon3D& pg = pV->dual_cell();
        glNormal3f(pV->normal()[0], pV->normal()[1], pV->normal()[2]);
        glColor3f(pV->rgb()[0], pV->rgb()[1], pV->rgb()[2]);
        for (size_t i = 0; i < pg.edges().size(); i++)
        {
            CSegment3D& s = pg.edges()[i];
            CPoint p = s.start();
            CPoint q = s.end();
            CPoint r = pV->dual_center();

            glVertex3f(r[0], r[1], 0);
            glVertex3f(p[0], p[1], 0);
            glVertex3f(q[0], q[1], 0);
        }
    }
    glEnd();
}

void draw_weighted_delaunay(COMTMesh* pMesh)
{
    glDisable(GL_LIGHTING);
    glColor3f(0, 1, 0);
    glLineWidth(1.0);
    glBegin(GL_LINES);

    for (COMTMesh::MeshEdgeIterator eiter(pMesh); !eiter.end(); ++eiter)
    {
        COMTMesh::CEdge* pE = *eiter;

        COMTMesh::CVertex* pV0 = pMesh->edgeVertex1(pE);
        COMTMesh::CVertex* pV1 = pMesh->edgeVertex2(pE);

        CPoint2 q0 = CPoint2(pV0->dual_center()[0], pV0->dual_center()[1]);
        CPoint2 q1 = CPoint2(pV1->dual_center()[0], pV1->dual_center()[1]);

        glVertex3d(q0[0], q0[1], 0);
        glVertex3d(q1[0], q1[1], 0);
    }
    glEnd();
    glEnable(GL_LIGHTING);
}

/*! display call back function
 */
void display()
{
    /* clear frame buffer */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    setupLight();
    /* transform from the eye coordinate system to the world system */
    setupEye();
    glPushMatrix();
    /* transform from the world to the ojbect coordinate system */
    setupObject();
    /* draw the mesh */
    switch (display_mode)
    {
        case 0:
            draw_voronoi_cell(pOT->pWeightedDT());
            break;
        case 1:
            draw_weighted_delaunay(pOT->pWeightedDT());
            break;
    }

    glPopMatrix();
    glutSwapBuffers();
}

/*! Called when a "resize" event is received by the window. */

void reshape(int w, int h)
{
    float ar;

    win_width = w;
    win_height = h;

    ar = (float) (w) / h;
    glViewport(0, 0, w, h); /* Set Viewport */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // magic imageing commands
    gluPerspective(40.0,   /* field of view in degrees */
                   ar,     /* aspect ratio */
                   0.0001, /* Z near */
                   100.0 /* Z far */);

    glMatrixMode(GL_MODELVIEW);

    glutPostRedisplay();
}

/*! helper function to remind the user about commands, hot keys */
void help()
{
    printf("-----------------------------------------------------------------\n");
    printf(" Demo for Optimal Transport Map from the unit disk to the planar domain \n");
    printf(" Input  : A planar mesh\n");
    printf(" Output : OT Map \n");
    printf(" Command: *.exe input_mesh\n");
    printf(" Debug  : Press '&' for Solving Monge-Amper Equation\n");
    printf("------------------------------------------------------------------\n");
    
    printf(" w  -  Wireframe Display\n");
    printf(" f  -  Flat Shading \n");
    printf(" s  -  Smooth Shading\n");
    
    printf(" !  -  Gradient Descent Method\n");
    printf(" &  -  Newton's Method\n");
    printf(" e  -  Toggle showing cell\n");
    printf(" g  -  Switch display modes\n");

    printf(" ?  -  Help Information\n");
    printf(" esc - quit\n");
    printf("------------------------------------------------------------------\n");
}

void key_process(unsigned char key)
{
    switch (key)
    {
        case 'f':
            // Flat Shading
            glPolygonMode(GL_FRONT, GL_FILL);
            // configure.m_shade_mode = 0;
            break;
        case 's':
            // Smooth Shading
            glPolygonMode(GL_FRONT, GL_FILL);
            // configure.m_shade_mode = 1;
            break;
        case 'w':
            // Wireframe mode
            glPolygonMode(GL_FRONT, GL_LINE);
            break;
        case '?':
            help();
            break;
        case 27:
            exit(0);
            break;
        case 'g':
            display_mode = (display_mode + 1) % 2;
            break;
        case 'e':
            show_cell = (show_cell + 1) % 2;
            break;
        case '!':
        {
            COMTMesh* pM = NULL;
            pOT->__gradient_descend(pOT->pWeightedDT(), pM);
            delete pOT->pWeightedDT();
            pOT->pWeightedDT() = pM;
            pOT->_compute_error(pM);
        }
        break;
        case '&':
        {
            COMTMesh* pM = NULL;
            pOT->__newton(pOT->pWeightedDT(), pM);
            delete pOT->pWeightedDT();
            pOT->pWeightedDT() = pM;
            pOT->_compute_error(pM);
        }
        break;
    }
}
/*! Keyboard call back function */

void keyBoard(unsigned char key, int x, int y)
{
    key_process(key);
    glutPostRedisplay();
}

/*! setup GL states */
void setupGLstate()
{

    GLfloat lightOneColor[] = {1, 1, 1, 1};
    GLfloat globalAmb[] = {.1f, .1f, .1f, 1.0f};
    GLfloat lightOnePosition[] = {.0f, .0f, 1.0f, 0.0f};

    glEnable(GL_CULL_FACE);
    glFrontFace(GL_CCW);
    // glFrontFace(GL_CW);
    glEnable(GL_DEPTH_TEST);
    // glClearColor(0,0,0,0);
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glShadeModel(GL_SMOOTH);

    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);

    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightOneColor);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, globalAmb);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);

    const GLfloat specular[] = {1.0f, 1.0f, 1.0f, 1.0f};

    glLightfv(GL_LIGHT1, GL_SPECULAR, specular);

    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 64.0f);

    GLfloat mat_ambient[] = {0.0f, 0.0f, 0.0f, 1.0f};
    GLfloat mat_diffuse[] = {0.01f, 0.01f, 0.01f, 1.0f};
    GLfloat mat_specular[] = {0.5f, 0.5f, 0.5f, 1.0f};
    GLfloat mat_shininess[] = {32};

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
}

/*! mouse click call back function */
void mouseClick(int button, int state, int x, int y)
{

    /* set up an arcball around the Eye's center
            switch y coordinates to right handed system  */

    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        gButton = GLUT_LEFT_BUTTON;
        arcball = CArcball(win_width, win_height, x - win_width / 2, win_height - y - win_height / 2);
    }

    if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN)
    {
        startx = x;
        starty = y;
        gButton = GLUT_MIDDLE_BUTTON;
    }

    if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
    {
        startx = x;
        starty = y;
        gButton = GLUT_RIGHT_BUTTON;
    }
    return;
}

/*! mouse motion call back function */
void mouseMove(int x, int y)
{
    CPoint trans;
    CQrot rot;

    /* rotation, call arcball */
    if (gButton == GLUT_LEFT_BUTTON)
    {
        rot = arcball.update(x - win_width / 2, win_height - y - win_height / 2);
        if (light_edit)
            LightRot = rot * LightRot;
        else
            ObjRot = rot * ObjRot;
        glutPostRedisplay();
    }

    /*xy translation */
    if (gButton == GLUT_MIDDLE_BUTTON)
    {
        double scale = 10. / win_height;
        trans = CPoint(scale * (x - startx), scale * (starty - y), 0);
        startx = x;
        starty = y;
        ObjTrans = ObjTrans + trans;
        glutPostRedisplay();
    }

    /* zoom in and out */
    if (gButton == GLUT_RIGHT_BUTTON)
    {
        double scale = 10. / win_height;
        trans = CPoint(0, 0, scale * (starty - y));
        startx = x;
        starty = y;
        ObjTrans = ObjTrans + trans;
        glutPostRedisplay();
    }
}

void process_key(char* keys)
{
    for (size_t i = 0; i < strlen(keys); i++)
    {
        unsigned char ch = keys[i];
        key_process(ch);
    }
}

/*! main function for viewer
 */
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        printf("Usage: %s mesh_name\n", argv[0]);
        return -1;
    }

    std::string mesh_name(argv[1]);
    {
        if (strutil::endsWith(mesh_name, ".m"))
        {
            mesh.read_m(mesh_name.c_str());
        }
        else
        {
            printf("Only file format .m supported.\n");
            return EXIT_FAILURE;
        }

        COperator<COMTMesh> pS(&mesh);
        pS._normalize();
        pS._calculate_face_vertex_normal();
        pS._calculate_face_vertex_area();
        pS._set_vertex_default_color();
    }

    {
        pOT = new CDomainOptimalTransport(&mesh);
        pOT->_initialize();
    }

    /* glut stuff */
    glutInit(&argc, argv); /* Initialize GLUT */
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(600, 600);
    std::string title("OT Mesh Viewer - ");
    title = title + mesh_name;
    glutCreateWindow(title.c_str()); /* Create window with given title */
    glViewport(0, 0, 800, 800);

    glutDisplayFunc(display); /* Set-up callback functions */
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseClick);
    glutMotionFunc(mouseMove);
    glutKeyboardFunc(keyBoard);
    setupGLstate();

    glutMainLoop(); /* Start GLUT event-processing loop */

    delete pOT;
    return 0;
}