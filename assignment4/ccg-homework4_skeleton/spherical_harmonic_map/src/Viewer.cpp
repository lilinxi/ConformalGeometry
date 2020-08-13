#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef MAC_OS

#include <GLUT/glut.h>

#else
#include <GL/glut.h>
#endif // MAC_OS

#include "viewer/Arcball.h" /*  Arc Ball  Interface         */
#include "SphericalHarmonicMapMesh.h"
#include "SphericalHarmonicMap.h"

using namespace MeshLib;

/* window width and height */
int g_win_width, g_win_height;
int g_button;
int g_startx, g_starty;
int g_shade_flag = 0;
bool g_show_mesh = true;
bool g_show_uv = false;

/* rotation quaternion and translation vector for the object */
CQrot g_obj_rot(0, 0, 1, 0);
CPoint g_obj_trans(0, 0, 0);

/* arcball object */
CArcball g_arcball;

/* global g_mesh */
CSHMMesh g_mesh;
CSphericalHarmonicMap g_mapper;

/*! setup the object, transform from the world to the object coordinate system */
void setupObject(void) {
    double rot[16];

    glTranslated(g_obj_trans[0], g_obj_trans[1], g_obj_trans[2]);
    g_obj_rot.convert(rot);
    glMultMatrixd((GLdouble *) rot);
}

/*! the eye is always fixed at world z = +5 */
void setupEye(void) {
    glLoadIdentity();
    gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);
}

/*! setup light */
void setupLight() {
    GLfloat lightOnePosition[4] = {0, 0, 1, 0};
    GLfloat lightTwoPosition[4] = {0, 0, -1, 0};
    glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);
    glLightfv(GL_LIGHT2, GL_POSITION, lightTwoPosition);
}

/*! draw g_mesh */
void drawMesh() {
    glEnable(GL_LIGHTING);

    glLineWidth(1.0);
    glColor3f(229.0 / 255.0, 162.0 / 255.0, 141.0 / 255.0);
    for (CSHMMesh::MeshFaceIterator fiter(&g_mesh); !fiter.end(); ++fiter) {
        glBegin(GL_POLYGON);
        CSHMMesh::CFace *pF = *fiter;
        for (CSHMMesh::FaceVertexIterator fviter(pF); !fviter.end(); ++fviter) {
            CSHMMesh::CVertex *pV = *fviter;
            CPoint &p = pV->point();
            CPoint n;
            switch (g_shade_flag) {
                case 0:
                    n = pF->normal();
                    break;
                case 1:
                    n = pV->normal();
                    break;
            }
            glNormal3d(n[0], n[1], n[2]);
            glVertex3d(p[0], p[1], p[2]);
        }
        glEnd();
    }
}

/*! draw g_mesh */
void drawUv() {
    glEnable(GL_LIGHTING);

    glLineWidth(1.0);
    glColor3f(229.0 / 255.0, 162.0 / 255.0, 141.0 / 255.0);
    for (CSHMMesh::MeshFaceIterator fiter(&g_mesh); !fiter.end(); ++fiter) {
        glBegin(GL_POLYGON);
        CSHMMesh::CFace *pF = *fiter;
        for (CSHMMesh::FaceVertexIterator fviter(pF); !fviter.end(); ++fviter) {
            CSHMMesh::CVertex *pV = *fviter;
            CPoint &u = pV->u();
            CPoint n;
            switch (g_shade_flag) {
                case 0:
                    n = pF->normal();
                    break;
                case 1:
                    n = pV->normal();
                    break;
            }
            glNormal3d(n[0], n[1], n[2]);
            glVertex3d(u[0], u[1], u[2]);
        }
        glEnd();
    }
}

/*! display call back function
 */
void display() {
    /* clear frame buffer */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    setupLight();
    /* transform from the eye coordinate system to the world system */
    setupEye();
    glPushMatrix();
    /* transform from the world to the ojbect coordinate system */
    setupObject();

    /* draw the mesh */
    if (g_show_mesh)
        drawMesh();
    if (g_show_uv)
        drawUv();

    glPopMatrix();
    glutSwapBuffers();
}

/*! Called when a "resize" event is received by the window. */
void reshape(int w, int h) {
    float ar;

    g_win_width = w;
    g_win_height = h;

    ar = (float) (w) / h;
    glViewport(0, 0, w, h); /* Set Viewport */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(40.0, /* field of view in degrees */
                   ar,   /* aspect ratio */
                   0.1,  /* Z near */
                   100.0 /* Z far */);

    glMatrixMode(GL_MODELVIEW);

    glutPostRedisplay();
}

/*! helper function to remind the user about commands, hot keys */
void help() {
    printf("1  -  Show or hide original mesh\n");
    printf("2  -  Show or hide harmonic map\n");
    printf("n  -  Take next one step\n");
    printf("N  -  Take next twenty step\n");
    printf("h  -  Compute the spherical harmonic map\n");
    printf("w  -  Wireframe Display\n");
    printf("f  -  Flat Shading \n");
    printf("s  -  Smooth Shading\n");
    printf("?  -  Help Information\n");
    printf("esc - Quit\n");
}

/*! Keyboard call back function */
void keyBoard(unsigned char key, int x, int y) {
    switch (key) {
        case '1':
            // Show or hide mesh
            g_show_mesh = !g_show_mesh;
            break;
        case '2':
            // Show or hide uv
            g_show_uv = !g_show_uv;
            break;
        case 'n':
            // take next one step
            g_mapper.step_one();
            break;
        case 'N':
            // take next twenty steps
            g_mapper.step_one(20);
            break;
        case 'h':
            // spherical harmonic map
            g_mapper.map();
            break;
        case 'f':
            // Flat Shading
            glPolygonMode(GL_FRONT, GL_FILL);
            g_shade_flag = 0;
            break;
        case 's':
            // Smooth Shading
            glPolygonMode(GL_FRONT, GL_FILL);
            g_shade_flag = 1;
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
    }
    glutPostRedisplay();
}

/*! setup GL states */
void setupGLstate() {
    GLfloat lightOneColor[] = {1, 1, 1, 1.0};
    GLfloat globalAmb[] = {.1, .1, .1, 1};
    GLfloat lightOnePosition[] = {.0, 0.0, 1.0, 1.0};
    GLfloat lightTwoPosition[] = {.0, 0.0, -1.0, 1.0};

    glEnable(GL_CULL_FACE);
    glFrontFace(GL_CCW);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.35, 0.53, 0.70, 0);
    glShadeModel(GL_SMOOTH);

    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);

    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightOneColor);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, lightOneColor);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, globalAmb);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);
    glLightfv(GL_LIGHT2, GL_POSITION, lightTwoPosition);

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
void mouseClick(int button, int state, int x, int y) {
    /* set up an arcball around the Eye's center
    switch y coordinates to right handed system  */

    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        g_button = GLUT_LEFT_BUTTON;
        g_arcball = CArcball(g_win_width,
                             g_win_height,
                             x - g_win_width / 2,
                             g_win_height - y - g_win_height / 2);
    }

    if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
        g_startx = x;
        g_starty = y;
        g_button = GLUT_MIDDLE_BUTTON;
    }

    if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
        g_startx = x;
        g_starty = y;
        g_button = GLUT_RIGHT_BUTTON;
    }
    return;
}

/*! mouse motion call back function */
void mouseMove(int x, int y) {
    CPoint trans;
    CQrot rot;

    /* rotation, call g_arcball */
    if (g_button == GLUT_LEFT_BUTTON) {
        rot = g_arcball.update(x - g_win_width / 2, g_win_height - y - g_win_height / 2);
        g_obj_rot = rot * g_obj_rot;
        glutPostRedisplay();
    }

    /*xy translation */
    if (g_button == GLUT_MIDDLE_BUTTON) {
        double scale = 10. / g_win_height;
        trans = CPoint(scale * (x - g_startx), scale * (g_starty - y), 0);
        g_startx = x;
        g_starty = y;
        g_obj_trans = g_obj_trans + trans;
        glutPostRedisplay();
    }

    /* zoom in and out */
    if (g_button == GLUT_RIGHT_BUTTON) {
        double scale = 10. / g_win_height;
        trans = CPoint(0, 0, scale * (g_starty - y));
        g_startx = x;
        g_starty = y;
        g_obj_trans = g_obj_trans + trans;
        glutPostRedisplay();
    }
}

/*! Normalize g_mesh
 * \param pMesh the input g_mesh
 */
void normalizeMesh(CSHMMesh *pMesh) {
    CPoint s(0, 0, 0);
    for (CSHMMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter) {
        CSHMMesh::CVertex *v = *viter;
        s = s + v->point();
    }
    s = s / pMesh->numVertices();

    for (CSHMMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter) {
        CSHMMesh::CVertex *v = *viter;
        CPoint p = v->point();
        p = p - s;
        v->point() = p;
    }

    double d = 0;
    for (CSHMMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter) {
        CSHMMesh::CVertex *v = *viter;
        CPoint p = v->point();
        for (int k = 0; k < 3; k++) {
            d = (d > fabs(p[k])) ? d : fabs(p[k]);
        }
    }

    for (CSHMMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter) {
        CSHMMesh::CVertex *v = *viter;
        CPoint p = v->point();
        p = p / d;
        v->point() = p;
    }
};

/*! Compute the face normal and vertex normal
 * \param pMesh the input g_mesh
 */
void computeNormal(CSHMMesh *pMesh) {
    for (CSHMMesh::MeshFaceIterator fiter(pMesh); !fiter.end(); ++fiter) {
        CSHMMesh::CFace *pF = *fiter;
        CPoint p[3];
        CSHMMesh::CHalfEdge *he = pMesh->faceHalfedge(pF);

        for (int k = 0; k < 3; k++) {
            CSHMMesh::CVertex *pv = pMesh->halfedgeTarget(he);
            p[k] = pv->point();
            he = pMesh->halfedgeNext(he);
        }

        CPoint fn = (p[1] - p[0]) ^(p[2] - p[0]);
        pF->area() = fn.norm() / 2.0; // 计算面片的面积
        pF->normal() = fn / fn.norm(); // 计算面片的法线方向
    }

    for (CSHMMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter) {
        CSHMMesh::CVertex *v = *viter;
        CPoint n(0, 0, 0);
        double area = 0;
        for (CSHMMesh::VertexFaceIterator vfiter(v); !vfiter.end(); ++vfiter) {
            CSHMMesh::CFace *pF = *vfiter;
            n += pF->normal() * pF->area(); // 点的法线是周围相连接面的法线的面积加权平均
            area += pF->area();// 点的面积是周围相连接面的面积的平均
        }
        v->area() = area / 3.0;
        n = n / n.norm();
        v->normal() = n;
    }
};

void initOpenGL(int argc, char *argv[]) {
    /* glut stuff */
    glutInit(&argc, argv); /* Initialize GLUT */
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(600, 600);
    glutCreateWindow("Mesh Viewer"); /* Create window with given title */
    glViewport(0, 0, 600, 600);

    glutDisplayFunc(display); /* Set-up callback functions */
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseClick);
    glutMotionFunc(mouseMove);
    glutKeyboardFunc(keyBoard);
    setupGLstate();

    glutMainLoop(); /* Start GLUT event-processing loop */
}

/*! main function for viewer
 */
int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage: %s input.m\n", argv[0]);
        return EXIT_FAILURE;
    }

    std::string mesh_name(argv[1]);
    if (strutil::endsWith(mesh_name, ".m")) {
        g_mesh.read_m(mesh_name.c_str());
    } else {
        printf("Only file format .m supported.\n");
        return EXIT_FAILURE;
    }

    normalizeMesh(&g_mesh);
    computeNormal(&g_mesh);

    g_mapper.set_mesh(&g_mesh);

    help();

    initOpenGL(argc, argv);
    return EXIT_SUCCESS;
}
