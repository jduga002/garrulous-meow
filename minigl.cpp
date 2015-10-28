#include <iostream>
/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 * Do not use any additional files
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include "minigl.h"

using namespace std;

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

/**
 * Error messages
 */
const char* GL_BEGIN_ERROR = "ERROR: mglBegin already called. Must use function after call to mglEnd";
const char* GL_END_ERROR = "ERROR: mglBegin not called";
const MGLfloat identity[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,1,0,
                               0,0,0,1};
class vertex {
  public:
    MGLfloat x;
    MGLfloat y;
    MGLfloat z;
    MGLfloat w;
    MGLpixel color;

    vertex();
    vertex(MGLfloat x, MGLfloat y, MGLfloat z, 
           MGLbyte red, MGLbyte green, MGLbyte blue) {
        setValues(x,y,z,1);
        setColor(red,green,blue);
    }
    vertex(const vertex& other)
     : color(other.color) {
        setValues(other.x, other.y, other.z, 1);
    }
    void setValues(int x, int y, int z) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = 1;
    }
    void setValues(int x, int y, int z, int w) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = w;
    }
    void setColor(MGLbyte red, MGLbyte green, MGLbyte blue) {
        MGLpixel color=0;
        MGL_SET_RED(color, red);
        MGL_SET_GREEN(color, green);
        MGL_SET_BLUE(color, blue);
        this->color = color;
    }
};

#define TRIANGLE 3
#define QUADRILATERAL 4
#define UNKNOWN_POLY -1
class MGL_Polygon {
    vector<vertex> vertices;
  public:
    MGL_Polygon(vertex A, vertex B, vertex C) {
        vertices.push_back(A);
        vertices.push_back(B);
        vertices.push_back(C);
    }
    MGL_Polygon(vertex A, vertex B, vertex C, vertex D) {
        vertices.push_back(A);
        vertices.push_back(B);
        vertices.push_back(C);
        vertices.push_back(D);
    }
    vector<vertex> &getVertices() {
        return vertices;
    }
    int polyType() {
        int type = vertices.size();
        if (type != TRIANGLE && type != QUADRILATERAL) {
            type = UNKNOWN_POLY;
        }
        return type;
    }
};

class MGLObject {
    bool mglBeginCalled;
    MGLpoly_mode mglPoly_Mode;
    MGLmatrix_mode mglMatrix_Mode;

    vector<MGL_Polygon> polygonList;
    vector<vertex> vertexList;
    vector< vector<MGLfloat> > projection_matrixStack;
    vector< vector<MGLfloat> > modelview_matrixStack;

    MGLbyte currentColor[3];
  public:
    MGLObject()
     : mglBeginCalled(false), projection_matrixStack(1), modelview_matrixStack(1) {
        vector<MGLfloat> id(identity, identity + 16);
        projection_matrixStack.at(0) = id;
        modelview_matrixStack.at(0) = id;

        currentColor[0] = 255;
        currentColor[1] = 255;
        currentColor[2] = 255;
    }
    
    void readPixels(MGLsize width,
                       MGLsize height,
                       MGLpixel *data) {
        cout << "Reading pixels" << endl;
        //TODO: also do quadrilaterals
        for (unsigned i = 0; i < polygonList.size(); i += 9) {
            MGL_Polygon &polygon = polygonList.at(i);
            if (polygon.polyType() == TRIANGLE) {
                draw_triangle(polygon.getVertices().at(0),
                              polygon.getVertices().at(1),
                              polygon.getVertices().at(2),
                              width, height, data);
            }
            else if (polygon.polyType() == QUADRILATERAL) {
                draw_triangle(polygon.getVertices().at(0),
                              polygon.getVertices().at(1),
                              polygon.getVertices().at(2),
                              width, height, data);
                draw_triangle(polygon.getVertices().at(1),
                              polygon.getVertices().at(2),
                              polygon.getVertices().at(3),
                              width, height, data);
            }
        } 
        for (unsigned i = 0; i < width*height; i++ ) {
            if (data[i] != 0)
                cout << "Pixel " << i << ":\t" << data[i] << endl;
        }
    }
   
    void begin(MGLpoly_mode mode) {
        if (mglBeginCalled) {
            MGL_ERROR(GL_BEGIN_ERROR);
        }
        mglBeginCalled = true;
        mglPoly_Mode = mode;
    }
    
    void end() {
        if (!mglBeginCalled) {
            MGL_ERROR(GL_END_ERROR);
        }
        mglBeginCalled = false;
        if (mglPoly_Mode == MGL_TRIANGLES) {
            for (unsigned i = 0; i < vertexList.size() / 3; i++) {
                polygonList.push_back(MGL_Polygon(vertexList.at(i),
                                                  vertexList.at(i+1),
                                                  vertexList.at(i+2)));
            }
        }
        else if (mglPoly_Mode == MGL_QUADS) {
            for (unsigned i = 0; i < vertexList.size() / 4; i++) {
                polygonList.push_back(MGL_Polygon(vertexList.at(i),
                                                  vertexList.at(i+1),
                                                  vertexList.at(i+2),
                                                  vertexList.at(i+3)));
            }
        }
        vertexList.erase(vertexList.begin(), vertexList.end());
    }

    void vertex3(MGLfloat x,
                    MGLfloat y,
                    MGLfloat z) {
        if (mglBeginCalled) {
            MGLfloat vertice[4] = {x,y,z,1};
            //multiply by camera matrix, which in this case is just identity
            vector<MGLfloat> &proj_matrix = projection_matrixStack.back();
            //multiply by projection matrix,
            vector<MGLfloat> &modview_matrix = modelview_matrixStack.back();
            for (unsigned i = 0; i < 4; i++) {
            }
            vertexList.push_back(vertex(x,y,z,currentColor[0],
                                              currentColor[1],
                                              currentColor[2]));
        }
    }
    void setMatrixMode(MGLmatrix_mode mode) {
        mglMatrix_Mode = mode;
    }
    void pushMatrix() {
        if (mglMatrix_Mode == MGL_MODELVIEW) {
            vector<MGLfloat> back = modelview_matrixStack.back();
            modelview_matrixStack.push_back(back);
        }
        else if (mglMatrix_Mode == MGL_PROJECTION) {
            vector<MGLfloat> back = projection_matrixStack.back();
            projection_matrixStack.push_back(back);
        }
    }
    void popMatrix() {
        if (mglMatrix_Mode == MGL_MODELVIEW) {
            if (modelview_matrixStack.size() > 1) {
                modelview_matrixStack.pop_back();
            }
        }
        else if (mglMatrix_Mode == MGL_PROJECTION) {
            projection_matrixStack.pop_back();
        }
    }
    void loadIdentity() {
        loadMatrix(identity);
    }
    void loadMatrix(const MGLfloat *matrix) {
        if (mglBeginCalled) {
            MGL_ERROR(GL_BEGIN_ERROR);
        }
        vector<MGLfloat> matrix_vec(matrix, matrix + 16);
        if (mglMatrix_Mode == MGL_MODELVIEW) {
            modelview_matrixStack.back() = matrix_vec;
        }
        else if (mglMatrix_Mode == MGL_PROJECTION) {
            projection_matrixStack.back() = matrix_vec;
        }
    }
    void multMatrix(const MGLfloat *matrix) {
        if (mglBeginCalled) {
            MGL_ERROR(GL_BEGIN_ERROR);
        }
        vector<MGLfloat> mult_matrix(16);
        vector<MGLfloat> *back = NULL;
        if (mglMatrix_Mode == MGL_MODELVIEW)
            back = &(modelview_matrixStack.back());
        else if (mglMatrix_Mode == MGL_PROJECTION)
            back = &(projection_matrixStack.back());
        for (unsigned i = 0; i < 4; i++) {
            for (unsigned j = 0; i < 4; i++) {
                MGLfloat sum = 0.0;
                for (unsigned k = 0; k < 4; k++) {
                    sum += (back->at(i+4*k))*(matrix[k+4*j]);
                }
                mult_matrix.at(i+4*j) = sum;
            }
        }
        *(back) = mult_matrix;
    }

    void translate(MGLfloat x,
                   MGLfloat y,
                   MGLfloat z) {
        //TODO: Implement me!!
        MGLfloat trans_matrix[16] = {1,0,0,0,
                                     0,1,0,0,
                                     0,0,1,0,
                                     x,y,z,1};
        multMatrix(trans_matrix);
    }
    void rotate(MGLfloat angle,
                   MGLfloat x,
                   MGLfloat y,
                   MGLfloat z) {
        //TODO: Implement me!!
    }
    void scale(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z) {
        //TODO: Implement me!!
    }
    void frustum(MGLfloat left,
                    MGLfloat right,
                    MGLfloat bottom,
                    MGLfloat top,
                    MGLfloat near,
                    MGLfloat far) {
        //TODO: Implement me!!
    }
    void ortho(MGLfloat left,
                  MGLfloat right,
                  MGLfloat bottom,
                  MGLfloat top,
                  MGLfloat near,
                  MGLfloat far) {
        //TODO: Implement me!!
    }
    void color(MGLbyte red,
                  MGLbyte green,
                  MGLbyte blue) {
        //TODO: Implement me!!
        //MGLbyte *currentColor = mgl.getCurrentColor();
        //currentColor[0] = red;
        //currentColor[1] = green;
        //currentColor[2] = blue;
    }
    void set_pixel(int x, int y, unsigned width, unsigned height, MGLpixel *data) {
    unsigned index = y*width+ x;
    //unsigned index = x*height + y;
    if (index < width*height) {
        MGLpixel pixel = 0;
        MGL_SET_RED(pixel, 255);
        MGL_SET_GREEN(pixel, 255);
        MGL_SET_BLUE(pixel, 255);
        data[index] = pixel;
    }
}

void draw_line(vertex& v0, vertex& v1, unsigned width, unsigned height, MGLpixel* data)
{
    float dx = v1.x - v0.x;
    float dy = v1.y - v0.y;
    
    if (dx == 0) {
        if (dy == 0) {
            set_pixel(v0.x, v0.y, width, height, data);
            return;
        }
        int x = v0.x;
        if (dy > 0) {
            for (int y = v0.y; y < v1.y; y++)
                set_pixel(x, y, width, height, data);
        }
        else if (dy < 0) {
            for (int y = v0.y; y > v1.y; y--)
                set_pixel(x, y, width, height, data);
        }
        return;
    }

    float m = dy/dx;
    
    if (m <= 1 && m >= -1) {
        if (dx > 0) {
            float y = v0.y;
            for(int x = v0.x; x < v1.x; ++x) {
                set_pixel(x, static_cast<int>(y + 0.5), width, height, data);
                y += m;
            }
         }

        else if (dx < 0) {
            float y = v0.y;
            for(int x = v0.x; x > v1.x; --x) {
                set_pixel(x, static_cast<int>(y + 0.5), width, height, data);
                y -= m;
            }
        }
    }
    else { // absolute value of slope is greater than 1
        // dy != 0 so we do not need to check it
        float m_inverse = dx/dy; // 1/m
        float x = v0.x;
        if (dy > 0) {
            for (int y = v0.y; y < v1.y; ++y) {
                 set_pixel(static_cast<int>(x + 0.5), y, width, height, data);
                 x += m_inverse;
            }
        }
        else { // dy < 0
            for (int y = v0.y; y > v1.y; --y) {
                 set_pixel(static_cast<int>(x + 0.5), y, width, height, data);
                 x -= m_inverse;
            }
        }
    }

    return;
}

/**
 * Helper function for drawing triangles
 */
void draw_triangle(vertex &v1, vertex &v2, vertex &v3,
                   const int width, const int height, MGLpixel* data) {
    draw_line(v1, v2, width, height, data);
    draw_line(v1, v3, width, height, data);
    draw_line(v2, v3, width, height, data);
}
};

static MGLObject mgl;

/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
    mgl.readPixels(width, height, data);
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
    mgl.begin(mode);
}

/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
    mgl.end();
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
    mglVertex3(x, y, 0.0);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
    mgl.vertex3(x,y,z);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
    mgl.setMatrixMode(mode);
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
    mgl.pushMatrix();
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
    mgl.popMatrix();
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
    mgl.loadIdentity();
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
    mgl.loadMatrix(matrix);
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
    mgl.multMatrix(matrix);
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
    mgl.translate(x,y,z);
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
    mgl.rotate(angle, x, y, z);
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
    mgl.scale(x,y,z);
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
    mgl.frustum(left,right,bottom,top,near,far);
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
    mgl.ortho(left,right,bottom,top,near,far);
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLbyte red,
              MGLbyte green,
              MGLbyte blue)
{
    mgl.color(red,green,blue);
}
