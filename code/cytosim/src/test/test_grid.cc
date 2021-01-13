// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on January 2008

#define DIM 2

#include "vector.h"
#include "random.h"
#include "glapp.h"
#include "glut.h"
#include "real.h"
#include "tictoc.h"
#include "gle.h"
using namespace gle;

#include "grid.h"
#include "grid_display.h"

extern Random RNG;

//area of the grid
const int range = 5;
real  left[] = {  -range,  -range,  0 };
real right[] = {   range,   range,  1 };
int   size[] = { 2*range,   range,  1 };


typedef Grid<real, DIM, unsigned> grid_type;
grid_type myGrid;

grid_type::index_type cell_indx;
int coord[DIM];
Vector3 pos(0,0,0);
Vector3 nod(0,0,0);


#define TEST_REGIONS
real  regionRadius = 1.5;

//------------------------------------------------------------------------------
void throwMarbles(int cnt)
{
    myGrid.setValues(0);
    real w[3] = { 0, 0, 0 };
    for( int n = 0; n < cnt; ++n )
    {
        w[0] = range * RNG.sreal();
        w[1] = range * RNG.sreal();
        w[2] = range * RNG.sreal();
        ++myGrid(w);
    }
}


void processNormalKey(unsigned char c, int x=0, int y=0)
{
    switch (c)
    {
        case 'p':
            myGrid.periodic(!myGrid.periodic());
            break;
            
        case 'i':
            //decrease region-radius
            if ( regionRadius > 1 )
                regionRadius -= 0.25;
            myGrid.createRoundRegions(regionRadius);
            //myGrid.createSquareRegions(regionRadius);
            glApp::flashText("radius = %f", regionRadius);
            break;

        case 'o':
            // increase region-radius
            regionRadius += 0.25;
            myGrid.createRoundRegions(regionRadius);
            //myGrid.createSquareRegions(regionRadius);
            glApp::flashText("radius = %f", regionRadius);
            break;

        case 'r':
            myGrid.createRoundRegions(regionRadius);
            glApp::flashText("radius = %f", regionRadius);
            break;

        case 's':
            myGrid.createSideRegions(regionRadius);
            break;
            
        case 'h':
            printf("Shift-click to position test-point\n");
            return;

        case ' ': 
            throwMarbles(32);
            break;

        default:
            glApp::processNormalKey(c,x,y);
            return;
    }
    glutPostRedisplay();
}


//------------------------------------------------------------------------------

///set callback for shift-click, with unprojected click position
void  processMouseClick(int, int, const Vector3 & a, int)
{
    pos = a;
    cell_indx = myGrid.index(pos);
    myGrid.setPositionFromIndex(nod, cell_indx, 0);
    myGrid.setCoordinatesFromIndex(coord, cell_indx);
    
    char str[32];

    if ( myGrid.hasRegions() )
    {
        real num = myGrid.sumValuesInRegion(cell_indx);
        snprintf(str, sizeof(str), "cell %u : coord %i %i : %.0f marbles",
                 cell_indx, coord[0], coord[1], num);
    } 
    else
    {
        snprintf(str, sizeof(str), "cell %u : coord %i %i", cell_indx, coord[0], coord[1]);
    }
    
    glApp::setMessage(str);
    glutPostRedisplay();
}

///set callback for shift-drag, with unprojected mouse and click positions
void  processMouseDrag(int mx, int my, Vector3 & a, const Vector3 & b, int m)
{
    processMouseClick(mx, my, b, m);
}

//------------------------------------------------------------------------------

static bool field_color2(const real& val, Vector2 const& pos)
{
    glColor3f(val/5.0, 0, 0);
    return true;
}

static bool field_color3(const real& val, Vector3 const& pos)
{
    glColor3f(val/5.0, 0, 0);
    return true;
}


void display()
{
    char str[16];
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

#if ( DIM == 3 )
    Vector3 dir(0,0,1);
    drawValues(myGrid, field_color3, dir);
#else
    drawValues(myGrid, field_color2);
#endif
    
    //--------------draw a grid in gray:
    glColor4f(1,1,1,0.6);
    glLineWidth(1.0);
    drawEdges(myGrid);

    //--------------draw content of cells
    glColor3f(0.5,0.5,1);
    for ( grid_type::index_type c = 0 ; c < myGrid.nbCells(); ++c )
    {
        int val = myGrid.icell(c);
        if ( val > 0 )
        {
            Vector x;
            myGrid.setPositionFromIndex(x, c, 0.5);
            glPointSize(sqrt(16*val));
            glBegin(GL_POINTS);
            gleVertex(x);
            glEnd();
        }
    }

    //-------------draw selected-cell
    glPointSize(8.0);
    glBegin(GL_POINTS);
    glColor4f(1,1,1,1);
    gleVertex(pos);
    glColor4f(1,1,0,1);
    gleVertex(nod);
    glEnd();

    //-------------draw region
    if ( myGrid.hasRegions() )
    {
        int * offset = 0;
        int nb = myGrid.getRegion(offset, cell_indx);

        glColor4f(1,1,1,0.7);
        for ( int ii = 0; ii < nb; ++ii )
        {
            Vector x;
            myGrid.setPositionFromIndex(x, cell_indx+offset[ii], 0.4);
            snprintf(str, sizeof(str), "%i", ii);
            gleDrawText(x, str, GLUT_BITMAP_HELVETICA_10);
        }
    }
    else
    {
        real vi = myGrid.interpolate(pos);
        snprintf(str, sizeof(str), "cell %u %f", cell_indx, vi);
        glApp::setMessage(str);
        
#if ( DIM == 2 )
        vi = myGrid.interpolate2D(pos);
#else
        vi = myGrid.interpolate3D(pos);
#endif
        snprintf(str, sizeof(str), "cell %u %f", cell_indx, vi);
        glApp::setMessageT(str);
    }
}


void speedTest()
{
    printf("Real test...");

    real  left[] = { 0, 0, 0};
    real right[] = { 1, 1, 1};
    int   size[] = { 10, 10, 10};

    Grid<float, 3, unsigned int> map;
    map.setDimensions(left, right, size);
    map.createCells();
    map.setValues(0);

    real w[3];
    for ( int cc=0; cc<10000; ++cc )
    {
        w[0] = RNG.preal();
        w[1] = RNG.preal();
        w[2] = RNG.preal();
        for ( int x = 0; x < 1000; ++x )
        {
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
        }
    }

    FILE* test = fopen("testgrid.out","w");
    map.printValues(test, 0);
    fclose(test);
    printf("wrote file testgrid.out\n");
}


void testInterpolate()
{
    real  left[] = { 0, 0, 0 };
    real right[] = { 1, 1, 1 };
    int   size[] = { 100, 100, 100 };
    
    const int MAX = 1 << 14;
    real  rand[MAX+3] = { 0 };
    for ( int i = 0; i < MAX+3; ++i )
        rand[i] = RNG.preal();
    
    Grid<double, 3, unsigned int> map;
    map.setDimensions(left, right, size);
    map.createCells();
    map.setValues(0);
    
    const int CNT = 1000000;
    for ( int cc = 0; cc < CNT; ++cc )
    {
        real w[] = { RNG.preal(), RNG.preal(), RNG.preal() };
        ++map( w );
    }

    
    real * vec[CNT];
    for ( int i = 0; i < CNT; ++i )
        vec[i] = rand + ( RNG.pint() & MAX );

    real sum = 0;
    TicToc::tic();
    for ( int r = 0; r < 100; ++r )
        for ( int cc = 0; cc < CNT; ++cc )
            sum += map.interpolate3D(vec[cc]) + map.interpolate3D(vec[cc]);
    printf("sum = %f\n", sum);
    TicToc::toc("interpolate3D", "\n");
    
    real som = 0;
    TicToc::tic();
    for ( int r = 0; r < 100; ++r )
        for ( int cc = 0; cc < CNT; ++cc )
            som += map.interpolate(vec[cc]) + map.interpolate(vec[cc]);
    printf("som = %f\n", sum);
    TicToc::toc("interpolate  ", "\n");
}


int main(int argc, char* argv[])
{
    RNG.seedTimer();

    if ( argc > 1 )
    {
        testInterpolate();
        return 0;
    }

    //initialize the grid:
    myGrid.setDimensions(left, right, size);
    myGrid.createCells();
    //myGrid.periodic(1);
    throwMarbles(8*(1<<DIM));

    glutInit(&argc, argv);
    glApp::setDimensionality(DIM);
    glApp::attachMenu(GLUT_RIGHT_BUTTON);
    glApp::actionFunc(processMouseClick);
    glApp::actionFunc(processMouseDrag);
    glApp::normalKeyFunc(processNormalKey);
    glApp::createWindow(display);
    glApp::setScale(2*range+1);

    glutMainLoop();

    return 0;
}
