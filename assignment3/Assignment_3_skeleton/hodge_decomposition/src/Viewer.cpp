#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef MAC_OS
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif // MAC_OS

#include "HodgeDecomposition.h"
#include "HodgeDecompositionMesh.h"
#include "WedgeProduct.h"
#include "bmp/RgbImage.h"
#include "viewer/Arcball.h" /*  Arc Ball  Interface         */

using namespace MeshLib;

/* window width and height */
int g_win_width, g_win_height;
int g_button;
int g_startx, g_starty;
int g_shade_flag = 0;
bool g_show_mesh = true;
bool g_show_uv = false;
bool g_show_boundary = false;

/* rotation quaternion and translation vector for the object */
CQrot g_obj_rot(0, 0, 1, 0);
CPoint g_obj_trans(0, 0, 0);

/* arcball object */
CArcball g_arcball;

/* global g_mesh */
CHodgeDecompositionMesh* g_domain_mesh = NULL;
std::vector<CHodgeDecompositionMesh*> g_meshes;
CHodgeDecomposition g_mapper;
int g_show_index = 0;

int g_texture_flag = 2;
/* texture id and g_image */
GLuint g_texture_name;
RgbImage g_image;

/*! initialize bitmap g_image texture */
void initializeBmpTexture()
{
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &g_texture_name);
    glBindTexture(GL_TEXTURE_2D, g_texture_name);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,   GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    int ImageWidth = g_image.GetNumCols();
    int ImageHeight = g_image.GetNumRows();
    GLubyte* ptr = (GLubyte*) g_image.ImageData();

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, ImageWidth, ImageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, ptr);

    if (g_texture_flag == 1)
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    else if (g_texture_flag == 2)
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glEnable(GL_TEXTURE_2D);
}

/*! setup the object, transform from the world to the object coordinate system */
void setupObject(void)
{
    double rot[16];

    glTranslated(g_obj_trans[0], g_obj_trans[1], g_obj_trans[2]);
    g_obj_rot.convert(rot);
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
    GLfloat lightOnePosition[4] = {0, 0, 1, 0};
    GLfloat lightTwoPosition[4] = {0, 0, -1, 0};
    glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);
    glLightfv(GL_LIGHT2, GL_POSITION, lightTwoPosition);
}

/*! draw g_mesh */
void drawMesh()
{
    glEnable(GL_LIGHTING);
    glBindTexture(GL_TEXTURE_2D, g_texture_name);
    if (g_texture_flag > 0)
        glEnable(GL_TEXTURE_2D);

    glLineWidth(1.0);
    for (CHodgeDecompositionMesh::MeshFaceIterator fiter(g_domain_mesh); !fiter.end(); ++fiter)
    {
        glBegin(GL_POLYGON);
        CHodgeDecompositionFace* pF = *fiter;
        for (CHodgeDecompositionMesh::FaceVertexIterator fviter(pF); !fviter.end(); ++fviter)
        {
            CHodgeDecompositionVertex* pV = *fviter;
            CPoint& p = pV->point();
            CPoint2& uv = pV->uv();
            CPoint& rgb = pV->rgb();
            CPoint n;
            switch (g_shade_flag)
            {
                case 0:
                    n = pF->normal();
                    break;
                case 1:
                    n = pV->normal();
                    break;
            }
            glNormal3d(n[0], n[1], n[2]);
            glTexCoord2d(uv[0], uv[1]);
            glColor3f(rgb[0], rgb[1], rgb[2]);
            glVertex3d(p[0], p[1], p[2]);
        }
        glEnd();
    }
}

/*! draw uv mesh */
void drawUv()
{
    glEnable(GL_LIGHTING);
    glBindTexture(GL_TEXTURE_2D, g_texture_name);
    if (g_texture_flag > 0)
        glEnable(GL_TEXTURE_2D);

    glLineWidth(1.0);
    glColor3f(229.0 / 255.0, 162.0 / 255.0, 141.0 / 255.0);
    for (CHodgeDecompositionMesh::MeshFaceIterator fiter(g_domain_mesh); !fiter.end(); ++fiter)
    {
        glBegin(GL_POLYGON);
        CHodgeDecompositionFace* pF = *fiter;
        for (CHodgeDecompositionMesh::FaceVertexIterator fviter(pF); !fviter.end(); ++fviter)
        {
            CHodgeDecompositionVertex* pV = *fviter;
            CPoint2& uv = pV->uv();
            CPoint& rgb = pV->rgb();
            CPoint n;
            switch (g_shade_flag)
            {
                case 0:
                    n = pF->normal();
                    break;
                case 1:
                    n = pV->normal();
                    break;
            }
            glNormal3d(n[0], n[1], n[2]);
            glTexCoord2d(uv[0], uv[1]);
            glColor3f(rgb[0], rgb[1], rgb[2]);
            glVertex3d(uv[0], uv[1], 0);
        }
        glEnd();
    }
}

/*! draw boundary 
 *  mode: 1, mesh; 2, uv
 */
void drawBoundary(int mode) 
{
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    glLineWidth(3.0);
    glColor3f(0, 0, 1);
    glBegin(GL_LINES);
    for (CHodgeDecompositionMesh::MeshEdgeIterator eit(g_domain_mesh); !eit.end(); ++eit)
    {
        CHodgeDecompositionMesh::CEdge* pE = *eit;
        if (!pE->boundary()) continue;

        CHodgeDecompositionMesh::CVertex* pA = g_domain_mesh->edgeVertex1(pE);
        CHodgeDecompositionMesh::CVertex* pB = g_domain_mesh->edgeVertex2(pE);

        if (mode == 1) // draw mesh
        {
            CPoint& a = pA->point();
            CPoint& b = pB->point();
            glVertex3d(a[0], a[1], a[2]);
            glVertex3d(b[0], b[1], b[2]);
        }
        else if (mode == 2) // draw uv
        {
            CPoint2& a = pA->uv();
            CPoint2& b = pB->uv();
            glVertex3d(a[0], a[1], 0);
            glVertex3d(b[0], b[1], 0);
        }
    }
    glEnd();
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
    switch (g_texture_flag)
    {
        case 0:
            glDisable(GL_TEXTURE_2D);
            break;
        case 1:
            glEnable(GL_TEXTURE_2D);
            glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
            break;
        case 2:
            glEnable(GL_TEXTURE_2D);
            glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
            break;
    }

    /* draw the mesh */
    if (g_show_mesh)
    {
        if (g_show_boundary)
            drawBoundary(1);
        drawMesh();
    }
    if (g_show_uv)
    {
        if (g_show_boundary)
            drawBoundary(2);
        drawUv();
    }

    glPopMatrix();
    glutSwapBuffers();
}

/*! Called when a "resize" event is received by the window. */
void reshape(int w, int h)
{
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
void help()
{
    printf("\n");
    printf("1  -  Show or hide mesh\n");
    printf("2  -  Show or hide uv\n");
    printf("n  -  Integrate the next 1-form\n");
    printf("w  -  Wireframe Display\n");
    printf("f  -  Flat Shading \n");
    printf("s  -  Smooth Shading\n");
    printf("t  -  Texture rendering mode\n");
    printf("b  -  Show or hide boundary\n");
    printf("?  -  Help Information\n");
    printf("esc - Quit\n");
}

/*! Keyboard call back function */
void keyBoard(unsigned char key, int x, int y)
{
    switch (key)
    {
        case '1':
            // Show or hide mesh
            g_show_mesh = !g_show_mesh;
            break;
        case '2':
            // Show or hide uv
            g_show_uv = !g_show_uv;
            break;
        case 'n':
            // integrate the next one
            g_show_index = (g_show_index + 1) % g_meshes.size();
            g_mapper.integration(g_meshes[g_show_index], g_domain_mesh);
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
        case 't':
            // Texture rendering mode
            g_texture_flag = (g_texture_flag + 1) % 3;
            break;
        case 'b':
            // Show or hide boundary
            g_show_boundary = !g_show_boundary;
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
void setupGLstate()
{
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
void mouseClick(int button, int state, int x, int y)
{
    /* set up an arcball around the Eye's center
    switch y coordinates to right handed system  */

    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        g_button = GLUT_LEFT_BUTTON;
        g_arcball = CArcball(g_win_width, 
                             g_win_height, 
                             x - g_win_width / 2, 
                             g_win_height - y - g_win_height / 2);
    }

    if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN)
    {
        g_startx = x;
        g_starty = y;
        g_button = GLUT_MIDDLE_BUTTON;
    }

    if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
    {
        g_startx = x;
        g_starty = y;
        g_button = GLUT_RIGHT_BUTTON;
    }
    return;
}

/*! mouse motion call back function */
void mouseMove(int x, int y)
{
    CPoint trans;
    CQrot rot;

    /* rotation, call g_arcball */
    if (g_button == GLUT_LEFT_BUTTON)
    {
        rot = g_arcball.update(x - g_win_width / 2, g_win_height - y - g_win_height / 2);
        g_obj_rot = rot * g_obj_rot;
        glutPostRedisplay();
    }

    /*xy translation */
    if (g_button == GLUT_MIDDLE_BUTTON)
    {
        double scale = 10. / g_win_height;
        trans = CPoint(scale * (x - g_startx), scale * (g_starty - y), 0);
        g_startx = x;
        g_starty = y;
        g_obj_trans = g_obj_trans + trans;
        glutPostRedisplay();
    }

    /* zoom in and out */
    if (g_button == GLUT_RIGHT_BUTTON)
    {
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
void normalizeMesh(CHodgeDecompositionMesh* pMesh)
{
    CPoint s(0, 0, 0);
    for (CHodgeDecompositionMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CHodgeDecompositionVertex* v = *viter;
        s = s + v->point();
    }
    s = s / pMesh->numVertices();

    for (CHodgeDecompositionMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CHodgeDecompositionVertex* v = *viter;
        CPoint p = v->point();
        p = p - s;
        v->point() = p;
    }

    double d = 0;
    for (CHodgeDecompositionMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CHodgeDecompositionVertex* v = *viter;
        CPoint p = v->point();
        for (int k = 0; k < 3; k++)
        {
            d = (d > fabs(p[k])) ? d : fabs(p[k]);
        }
    }

    for (CHodgeDecompositionMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CHodgeDecompositionVertex* v = *viter;
        CPoint p = v->point();
        p = p / d;
        v->point() = p;
    }
};

/*! Compute the face normal and vertex normal
 * \param pMesh the input g_mesh
 */
void computeNormal(CHodgeDecompositionMesh* pMesh)
{
    for (CHodgeDecompositionMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CHodgeDecompositionVertex* v = *viter;
        CPoint n(0, 0, 0);
        for (CHodgeDecompositionMesh::VertexFaceIterator vfiter(v); !vfiter.end(); ++vfiter)
        {
            CHodgeDecompositionFace* pF = *vfiter;

            CPoint p[3];
            CHalfEdge* he = pF->halfedge();
            for (int k = 0; k < 3; k++)
            {
                p[k] = he->target()->point();
                he = he->he_next();
            }

            CPoint fn = (p[1] - p[0]) ^ (p[2] - p[0]);
            pF->normal() = fn / fn.norm();
            n += fn;
        }

        n = n / n.norm();
        v->normal() = n;
    }
};

void initOpenGL(int argc, char* argv[])
{
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
    initializeBmpTexture();
    
    glutMainLoop(); /* Start GLUT event-processing loop */
}

/*!
 *   verify if the mesh is closed or with boundary
 *   return true if closed; false, if open;
 */
bool closed_mesh(CHodgeDecompositionMesh* pM)
{
    using M = CHodgeDecompositionMesh;
    for (M::MeshVertexIterator viter(pM); !viter.end(); viter++)
    {
        M::CVertex* pv = *viter;
        if (pv->boundary())
            return false;
    }
    return true;
}

/*!
 *   calculate the holomorphic 1-form   
 */
void calc_holo_1_form(const std::string& input_mesh)
{
    using M = CHodgeDecompositionMesh;
    M* pM = new M;
    pM->read_m(input_mesh.c_str());
    g_meshes.push_back(pM);

    std::srand((unsigned) time(NULL));
    if (!closed_mesh(pM))
    {
        M::CBoundary bnd(pM);
        size_t n = bnd.loops().size() - 1;

        for (size_t i = 1; i < 2 * n; i++)
        {
            pM = new M;
            pM->read_m(input_mesh.c_str());
            g_meshes.push_back(pM);
        }

        size_t i = 0;
        for (std::vector<M*>::iterator miter = g_meshes.begin(); miter != g_meshes.end(); miter++, i++)
        {
            M* pM = *miter;
            g_mapper.set_mesh(pM);
            
            if (i < n)
            {
                g_mapper.exact_harmonic_form(i + 1);
            }
            else
            {
                g_mapper.random_harmonic_form();
            }
        }

        CBaseHolomorphicForm<M> holo(g_meshes);
        holo.conjugate();

        // normalize boundary integration
        i = 0;
        for (std::vector<M*>::iterator iter = g_meshes.begin(); iter != g_meshes.end(); iter++, i++)
        {
            if (i >= n)
                continue;

            M* pM = *iter;
            M::CBoundary bnd(pM);
            M::CLoop* pL = bnd.loops()[0];
            std::list<M::CHalfEdge*>& hs = pL->halfedges();

            double s = 0;
            for (std::list<M::CHalfEdge*>::iterator hiter = hs.begin(); hiter != hs.end(); hiter++)
            {
                M::CHalfEdge* ph = *hiter;
                M::CEdge* pe = pM->halfedgeEdge(ph);
                s += pe->duv()[1];
            }
            std::cout << "integration along the boundary is " << s << std::endl;

            for (M::MeshEdgeIterator eiter(pM); !eiter.end(); eiter++)
            {
                M::CEdge* pe = *eiter;
                pe->duv() /= s;
            }

            s = 0;
            for (std::list<M::CHalfEdge*>::iterator hiter = hs.begin(); hiter != hs.end(); hiter++)
            {
                M::CHalfEdge* ph = *hiter;
                M::CEdge* pe = pM->halfedgeEdge(ph);
                s += pe->duv()[1];
            }
            std::cout << "integration along the boundary is " << s << std::endl;
        }
    }
    else
    {
        int Euler = pM->numVertices() + pM->numFaces() - pM->numEdges();
        int genus = (2 - Euler) / 2;

        for (int i = 1; i < 2 * genus; i++)
        {
            M* pM = new M;
            pM->read_m(input_mesh.c_str());
            g_meshes.push_back(pM);
        }

        for (std::vector<M*>::iterator iter = g_meshes.begin(); iter != g_meshes.end(); iter++)
        {
            M* pM = *iter;
            normalizeMesh(pM);
            computeNormal(pM);
            g_mapper.set_mesh(pM);
            g_mapper.random_harmonic_form();
        }

        CBaseHolomorphicForm<M> holo(g_meshes);
        holo.conjugate();
    }

    // integrate the first one
    if (g_meshes.size() > 0)
    {
        g_mapper.integration(g_meshes[0], g_domain_mesh);
    }
}

/*! main function for viewer
 */
int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        printf("Usage: %s input.m input_open.m output_file_prefix\n", argv[0]);
        return EXIT_FAILURE;
    }

    std::string input_mesh_name(argv[1]);
    if (!strutil::endsWith(input_mesh_name, ".m"))
    {
        printf("Usage: %s input.m input_open.m output_file_prefix\n", argv[0]);
        return EXIT_FAILURE;
    }

    g_domain_mesh = new CHodgeDecompositionMesh;
    g_domain_mesh->read_m(argv[2]);
    normalizeMesh(g_domain_mesh);
    computeNormal(g_domain_mesh);

    calc_holo_1_form(input_mesh_name);

    // load texture
    g_image.LoadBmpFile(argv[3]);
    
    help();
    
    initOpenGL(argc, argv);

    return EXIT_SUCCESS;
}
