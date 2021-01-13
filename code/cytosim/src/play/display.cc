// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "display.h"
#include "hand_prop.h"
#include "sphere_prop.h"
#include "fiber_prop.h"
#include "point_disp.h"
#include "fiber_disp.h"
#include "line_disp.h"
#include "opengl.h"
#include "modulo.h"
#include "simul.h"
#include "gle.h"
#include "gle_color.h"
#include "gle_color_list.h"
#include "glapp.h"
#include "glut.h"

using namespace gle;
extern Modulo const* modulo;


//------------------------------------------------------------------------------
#pragma mark -


Display::Display(DisplayProp const* dp)
: pixelSize(1), uFactor(1), sFactor(1), prop(dp)
{
    assert_true(dp);
    
    fiber_prep_time = -1;
}

void Display::setPixelFactors(real ps, real u)
{
    pixelSize = ps;
    uFactor   = u;
    /*
     the 0.5 comes from the fact that glPointSize uses diameter
     while most gle::primitives use radius as arguments
     */
    sFactor    = u * 0.5 * pixelSize;
}


void Display::displayObjects(Simul const& sim)
{
    gleDrawText("Empty Display::display", GLUT_BITMAP_8_BY_13);
}


void Display::display(Simul const& sim)
{
    // clear list of transparent objects
    zObjects.clear();

#if ( DIM == 3 )
    glEnable(GL_LIGHTING);
#else
    glDisable(GL_LIGHTING);
#endif
    
    displaySpaces(sim.spaces);
    
    /*
     Draw opaque objects:
     - depth buffer is writable
     - glColor specifies the Front material color
     - back material color is set as 'inner_color'
     */

#if ( DIM == 3 )
    
    glEnable(GL_LIGHTING);
    glDepthMask(GL_TRUE);
    
#endif
    
    displayObjects(sim);
    
    /*
     Draw translucent objects:
     - make depth buffer readible only
     - objects are depth-sorted, from far to near
     - Dual pass is used to display back before front
     */

#if ( DIM == 3 )
    
    glEnable(GL_LIGHTING);
    glDepthMask(GL_FALSE);
    
    displayTransparentObjects(zObjects);
    
    displaySpacesF(sim.spaces);

#endif
    
    glEnable(GL_LIGHTING);
    glDepthMask(GL_TRUE);
    
#ifndef NDEBUG
    gleReportErrors(stderr, "in Display::display()");
#endif
}


/**
 To get correct display, it would be necessary to display all opaque objects first,
 and then all transparent objects for all tiles. Here, we calls Display::display()
 a number of times, and objects are sorted within each tile. The result is not perfect.
 */
void Display::displayTiled(Simul const& sim, int nine)
{
    assert_true(modulo);
    
    int l[3] = { 0 };
    int u[3] = { 0 };
    
    for ( int d = 0; d < DIM; ++d )
    {
        if ( modulo->isPeriodic(d) )
        {
            l[d] = (nine>1) ? -1 : 0;
            u[d] = +1;
        }
    }
    
    glMatrixMode(GL_MODELVIEW);
    
    Vector px = modulo->periodicity(0);
    Vector py = modulo->periodicity(1);
    Vector pz = modulo->periodicity(2);

    for ( int dx = l[0]; dx <= u[0]; ++dx )
    for ( int dy = l[1]; dy <= u[1]; ++dy )
    for ( int dz = l[2]; dz <= u[2]; ++dz )
    {
        Vector T = dx * px + dy * py + dz * pz;
        gleTranslate( T);
        display(sim);
        gleTranslate(-T);
    }
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 Create a FiberDisp for this Property if necessary
 */
void Display::prepareFiberDisp(FiberProp * p, PropertyList& alldisp, gle_color col)
{
    assert_true(p);
    
    FiberDisp *& disp = p->disp;
    if ( disp == 0 )
    {
        disp = new FiberDisp(p->name());
        alldisp.push_back(disp);
        // set default:
        disp->color       = col;
        disp->inner_color = col.darken(0.5);
        disp->point_size  = prop->point_size;
        disp->line_width  = prop->line_width;
    }
    
    if ( p->display_fresh )
    {
        try {
            disp->read_string(p->display);
        } catch(Exception & e) {
            std::clog << "fiber:display: " << e.what() << std::endl;
        }
        p->display_fresh = false;
    }
    
    if ( disp->coloring == FiberDisp::COLORING_FLECK )
        fiber_prep |= 1;
    
    if ( disp->line_style == 2 )
        fiber_prep |= 2;

    if ( disp->coloring == FiberDisp::COLORING_AGE )
        fiber_prep |= 4;
}


/**
 set LineDisp for given Fiber
 */
void Display::prepareLineDisp(const Fiber * fib)
{
    assert_true(fib->prop);
    FiberDisp const* disp = fib->prop->disp;
    
    if ( fib->disp == 0 )
        fib->disp = new LineDisp();
    LineDisp * self = fib->disp;

    // default values set from class:
    self->visible = disp->visible;
    self->color   = disp->color;
    
    // change body color depending on coloring mode:
    switch ( disp->coloring )
    {
        case FiberDisp::COLORING_RANDOM:
            self->color = gle::bright_color(fib->signature());
            break;
        case FiberDisp::COLORING_DIRECTION:
            // use the scalar product with direction defined by FiberDisp:right
            if ( fib->diffPoints(0) * disp->right_direction > 0 )
                self->color = disp->color_right;
            else
                self->color = disp->color_left;
            break;
        case FiberDisp::COLORING_MARK:
            self->color = gle::nice_color(fib->mark());
            break;
        case FiberDisp::COLORING_FLECK:
            self->color = gle::std_color(fib->flag());
            break;
        case FiberDisp::COLORING_AGE:
            self->color = gle::jet_color((fib->age()-age_min)*age_range, 1);
            break;
    }
    
#if ( 0 )
    // colors of ends set to match body color:
    self->end_color[0] = self->color;
    self->end_color[1] = self->color;
#else
    // colors of ends for non-dynamic filaments:
    self->end_color[0] = disp->end_color[0];
    self->end_color[1] = disp->end_color[0];
#endif
    
    // For dynamic Fibers, change colors of tips according to state:
    if ( fib->dynamicStateP() > 0 )
        self->end_color[0] = disp->end_color[fib->dynamicStateP()%5];
    
    if ( fib->dynamicStateM() > 0 )
        self->end_color[1] = disp->end_color[fib->dynamicStateM()%5];

    // hide right or left-pointing fibers:
    if ( ( disp->exclude & 1 )  &&  fib->diffPoints(0)*disp->right_direction < 0 )
        self->color = disp->colorT;
    if ( ( disp->exclude & 2 )  &&  fib->diffPoints(0)*disp->right_direction > 0 )
        self->color = disp->colorT;
    
#if ( DIM == 2 )
    // hide clockwise or counter-clockwise orientated fibers:
    if ( ( disp->exclude & 4 )  &&  cross(fib->posP(0), fib->diffPoints(0)) < 0 )
        self->color = disp->colorT;
    if ( ( disp->exclude & 8 )  &&  cross(fib->posP(0), fib->diffPoints(0)) > 0 )
        self->color = disp->colorT;
#elif ( DIM == 3 )
    // hide clockwise or counter-clockwise orientated fibers in the XY plane
    if ( ( disp->exclude & 4 )  &&  cross(fib->posP(0), fib->diffPoints(0)).ZZ < 0 )
        self->color = disp->colorT;
    if ( ( disp->exclude & 8 )  &&  cross(fib->posP(0), fib->diffPoints(0)).ZZ > 0 )
        self->color = disp->colorT;
#endif
    
#if ( 1 )
    // hide fibers depending on mask
    if ( fib->signature() & disp->mask_bitfield )
        self->color = disp->colorT;
#else
    if ( fib->mark() & disp->mask_bitfield )
        self->color = disp->colorT;
#endif

#if ( 0 )
    // hide fibers which are not growing
    if ( fib->dynamicStateP() == STATE_WHITE )
    {
        PRINT_ONCE("non-growing fibers made invisible\n");
        self->visible = -1;
        self->color = disp->colorT;
        self->end_color[0] = self->color.alpha_one();
    }
#endif

    // set visible flag according to body color:
    if ( !self->color.visible() )
        self->visible = 0;
    else if ( self->color.transparent() )
        self->visible = -1;

    // set parameters for exploded display
    if ( disp->explode )
        self->explode_shift = ( lcrng3(fib->signature(),3) * 0x1p-32 - 0.5 ) * disp->explode_range;
    else
        self->explode_shift = 0;
}


/**
 Create a PointDisp for this Property if necessary
 */
template < typename T >
void Display::preparePointDisp(T * p, PropertyList& alldisp, gle_color col)
{
    assert_true(p);
        
    PointDisp *& disp = p->disp;
    if ( disp == 0 )
    {
        disp = new PointDisp(p->category()+":display", p->name());
        alldisp.push_back(disp);
        // set default:
        disp->color  = col;
        disp->color2 = col.darken(0.5);
        disp->size   = prop->point_size;
        if ( p->category() == "hand" )
            disp->width = prop->link_size;
        else
            disp->width = prop->line_width;
    }
    
    // parse display string once:
    if ( p->display_fresh )
    {
        p->display_fresh = false;
        //std::clog << "reading display=(" << p->display << ") for Property " << p->name() << std::endl;
        try {
            disp->read_string(p->display);
        } catch(Exception & e) {
            std::cerr << "Error while reading " << disp->category() << ": " << e.what() << std::endl;
        }
    }
    
    disp->prepare(uFactor, sFactor);
}

/**
 Perform the operations that are necessary to display the simulation:
 - create FiberDisp, HandDisp, SphereDisp, etc. (one per Property)
 - create LineDisp (one per Fiber)
 - set default values,
 - parse display strings
 .
*/
void Display::prepareForDisplay(Simul const& sim, PropertyList& alldisp)
{
    //std::clog << "prepareForDisplay" << std::endl;

    if ( prop->fold )
        sim.foldPosition();
    
    // counter to give different colors to the objects
    unsigned int idx = 0;
    
    fiber_prep = 0;
    PropertyList plist = sim.properties.find_all("fiber");
    
    // create a FiberDisp for each FiberProp:
    for ( PropertyList::iterator pi = plist.begin(); pi < plist.end(); ++pi )
        prepareFiberDisp(static_cast<FiberProp*>(*pi), alldisp, gle::nice_color(idx++));

    // the cluster analysis only needs to be done once per state:
    if ( fiber_prep & 1  &&  fiber_prep_time != sim.time() )
    {
        sim.flagClusters();
        fiber_prep_time = sim.time();
    }
    
    // if fiber tensions are used for display, recompute them now:
    if ( fiber_prep & 2 )
        sim.computeForces();
    
    // calculate age() range and set color scaling factor:
    if ( fiber_prep & 4 )
    {
        unsigned cnt;
        real avg, dev, mn, mx;
        FiberSet::infoBirthtime(sim.fibers.collect(), cnt, avg, dev, mn, mx);
        if ( mx > mn )
        {
            age_min = mn;
            age_range = 5.0 / ( mx - mn );
            if ( mn != age_min )
                std::clog << "age min " << mn << " max " << mx << " range =" << age_range << std::endl;
        }
        else
        {
            age_min = 0;
            age_range = 1;
        }
    }
    
    // attribute LineDisp, and set individual display values for all fibers
    for ( Fiber * fib = sim.fibers.first(); fib; fib = fib->next() )
        prepareLineDisp(fib);
    
    //create a PointDisp for each HandProp:
    plist = sim.properties.find_all("hand");
    for ( PropertyList::iterator pi = plist.begin(); pi < plist.end(); ++pi )
        preparePointDisp(static_cast<HandProp*>(*pi), alldisp, gle::nice_color(idx++));
    
    //create a PointDisp for each SphereProp:
    plist = sim.properties.find_all("sphere");
    for ( PropertyList::iterator pi = plist.begin(); pi < plist.end(); ++pi )
        preparePointDisp(static_cast<SphereProp*>(*pi), alldisp, gle::bright_color(idx++));
    
    //create a PointDisp for each BeadProp:
    plist = sim.properties.find_all("bead", "solid");
    for ( PropertyList::iterator pi = plist.begin(); pi < plist.end(); ++pi )
        preparePointDisp(static_cast<BeadProp*>(*pi), alldisp, gle::bright_color(idx++));
    
    //create a PointDisp for each SpaceProp:
    gle_color col(DIM==3?0x000000AA:0xAAAAAAFF);
    plist = sim.properties.find_all("space");
    for ( PropertyList::iterator pi = plist.begin(); pi < plist.end(); ++pi )
        preparePointDisp(static_cast<SpaceProp*>(*pi), alldisp, col);
}


/**
 if `coloring` is enabled, this loads the N-th bright color,
 otherwise load the object' display color
 */
void Display::bodyColor(PointDisp const* disp, ObjectID nb) const
{
    if ( disp->coloring )
    {
        bright_color(nb).load_both();
        bright_color(nb).load();
    }
    else
    {
        disp->color2.load_back();
        disp->color.load_load(1.0);
    }
}

/**
 if `coloring` is enabled, this loads the N-th bright color,
 with an alpha value matched to the one of the object's display color.
 */
void Display::bodyColorT(PointDisp const* disp, ObjectID nb) const
{
    if ( disp->coloring )
    {
        gle_color col = bright_color(nb);
        col.load_load(disp->color.alpha());
        col.load_back(disp->color.alpha());
    }
    else
    {
        disp->color2.load();
        //disp->color.load_emission();
        if ( disp->color.transparent() )
            disp->color.load_both();
        else
        {
            disp->color2.load_back();
            disp->color.load_front();
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark -


void Display::displaySpaces(SpaceSet const& set)
{
#if ( DIM == 3 )
    
    // draw all non-transparent Spaces first:
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
    glDepthMask(GL_TRUE);
    for ( Space * obj = set.first(); obj; obj=obj->next() )
    {
        const PointDisp * disp = obj->prop->disp;
        
        if ( obj->prop->disp->visible & 1 && !disp->color2.transparent() )
        {
            lineWidth(disp->width);
            disp->color2.load_back();
            obj->display();
        }
    }
    
    // draw all transparent Spaces first:
    glDepthMask(GL_FALSE);
    for ( Space * obj = set.first(); obj; obj=obj->next() )
    {
        const PointDisp * disp = obj->prop->disp;
        if ( obj->prop->disp->visible & 1 && disp->color2.transparent() )
        {
            lineWidth(disp->width);
            disp->color2.load_back();
            obj->display();
        }
    }
    glDepthMask(GL_TRUE);
    glDisable(GL_CULL_FACE);

#else
    
    for ( Space * obj = set.first(); obj; obj=obj->next() )
    {
        const PointDisp * disp = obj->prop->disp;
        if ( disp->visible )
        {
            lineWidth(disp->width);
            disp->color.load_load();
            obj->display();
        }
    }
    
#endif
}


/**
 This is called only in 3D
 Note that if the color is not transparent, this will hide everything inside.
 */
void Display::displaySpacesF(SpaceSet const& set)
{
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glDepthMask(GL_TRUE);
    for ( Space * obj = set.first(); obj; obj=obj->next() )
    {
        const PointDisp * disp = obj->prop->disp;
        if ( obj->prop->disp->visible & 2 && !disp->color.transparent() )
        {
            disp->color.load_load();
            obj->display();
        }
    }
    
    glDepthMask(GL_FALSE);
    for ( Space * obj = set.first(); obj; obj=obj->next() )
    {
        const PointDisp * disp = obj->prop->disp;
        if ( obj->prop->disp->visible & 2 && disp->color.transparent() )
        {
            disp->color.load_load();
            obj->display();
        }
    }
    glDepthMask(GL_TRUE);
    glDisable(GL_CULL_FACE);
}


/**
 This displays only one Field, specified by DisplayProp:field_number
 
 GL_CULL_FACE and GL_LIGHTING should be disabled
 */
void Display::displayFields(FieldSet const& set)
{
#if ( DIM == 3 )
    // get current modelview transformation:
    GLfloat mat[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, mat);
    
    // extract axis corresponding to vertical direction:
    Vector3 dir(mat[2], mat[6], mat[10]);
    dir.normalize();
#else
    Vector3 dir(0,0,1);
#endif
    
#if ( 0 )
    
    for ( Field * obj = set.first(); obj; obj=obj->next() )
        if ( obj->hasField() && obj->prop->visible )
            obj->display(false, dir, 0);
    
#else
    
    Field * obj = set.first();
    
    if ( obj && obj->hasField() )
    {
        if ( obj->prop->visible == 1 )
            obj->display(false, dir, 0);
        else if ( obj->prop->visible == 2 )
            obj->display();
    }
#endif
}


//------------------------------------------------------------------------------
#pragma mark -

void Display::displayAverageFiber(ObjectList const& objs)
{
    Vector G, D, M, P;
    real S = FiberSet::infoPosition(objs, M, G, P);
    
    if ( S > REAL_EPSILON )
    {
        Vector MP = ( P - M ).normalized();
        gleCylinder(M, MP, 10*pixelSize);
        gleCone(P, MP, 10*pixelSize);
        gleObject(G, 10*pixelSize, gleSphere2B);
    }
}


bool selectR(Object const* obj, void const* arg)
{
    Fiber const* fib = static_cast<Fiber const*>(obj);
    return ( fib->prop == arg ) && ( fib->diffPoints(0)*fib->prop->disp->right_direction > 0 );
}

bool selectL(Object const* obj, void const* arg)
{
    Fiber const* fib = static_cast<Fiber const*>(obj);
    return ( fib->prop == arg ) && ( fib->diffPoints(0)*fib->prop->disp->right_direction < 0 );
}

void Display::displayAverageFiber1(FiberSet const& fibers, void const* arg)
{
    ObjectList objs = fibers.collect(match_property, arg);

#if ( 1 )
    // highlight with a black outline
    glLineWidth(3);
    glColor3f(0,0,0);
    glDepthMask(GL_FALSE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    displayAverageFiber(objs);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);        
    glDepthMask(GL_TRUE);
#endif
    
    glColor3f(1,1,1);
    displayAverageFiber(objs);
}


void Display::displayAverageFiber2(FiberSet const& fibers, void const* arg)
{
    ObjectList objsR = fibers.collect(selectR, arg);
    ObjectList objsL = fibers.collect(selectL, arg);

#if ( 1 )
    // highlight with a black outline
    glLineWidth(3);
    glColor3f(0,0,0);
    glDepthMask(GL_FALSE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    displayAverageFiber(objsR);
    displayAverageFiber(objsL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);        
    glDepthMask(GL_TRUE);
#endif
    
    // display right-pointing fibers in Red
    glColor3f(1,0,0);
    displayAverageFiber(objsR);
    
    // display left-pointing fibers in Green
    glColor3f(0,1,0);
    displayAverageFiber(objsL);
}    


void Display::displayMisc(Simul const& sim)
{
#if ( 0 )
    // display Steric Grid for visual debugging:
    glLineWidth(0.5);
    glColor3f(0, 0, 1);
    sim.stericGrid.display();
#endif
    
    PropertyList plist = sim.properties.find_all("fiber");
    
    for ( PropertyList::iterator pi = plist.begin(); pi < plist.end(); ++pi )
    {
        FiberProp* fp = static_cast<FiberProp*>(*pi);
        if ( fp->disp->show_average == 1 )
            displayAverageFiber1(sim.fibers, fp);
        else if ( fp->disp->show_average == 2 )
            displayAverageFiber2(sim.fibers, fp);
    }
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Display the MINUS_END of a Fiber, according to `style`:
 - 1: draw a sphere
 - 2: draw a cone
 - 3: draw a flat cylinder
 - 4: draw an arrow-head
 - 5: arrow-head in reverse direction
 .
 */
void Display::displayFiberEndMinus(Fiber const& fib) const
{
    const int IM = 1;
    const FiberDisp * disp = fib.prop->disp;
    
    real width = disp->end_size[IM]*sFactor;
    if ( disp->end_style[IM]  &&  width > 0 )
    {
        lineWidth(disp->line_width);
        switch(disp->end_style[IM])
        {
            case 1:
                gleObject(fib.posP(0), width, gleSphere2B);
                break;
            case 2:
                gleCone(fib.posP(0), -fib.dirPoint(0), width);
                break;
            case 3:
                gleCylinder(fib.posP(0), -fib.dirPoint(0), width);
                break;
            case 4:
                gleArrowTail(fib.posP(0), fib.dirPoint(0), width);
                break;
            case 5:
                gleArrowTail(fib.posP(0), -fib.dirPoint(0), width);
                break;
        }
    }
}


/**
 Display the PLUS_END of a Fiber, according to `style`:
 - 1: draw a sphere
 - 2: draw a cone
 - 3: draw a flat cylinder
 - 4: draw an arrow-head
 - 5: arrow-head in reverse direction
 .
 */
void Display::displayFiberEndPlus(Fiber const& fib) const
{
    const int IP = 0;
    const FiberDisp * disp = fib.prop->disp;
    
    real width = disp->end_size[IP]*sFactor;
    if ( disp->end_style[IP]  &&  width > 0 )
    {
        lineWidth(disp->line_width);
        switch(disp->end_style[IP])
        {
            case 1:
                gleObject(fib.posEndP(), width, gleSphere2B);
                break;
            case 2:
                gleCone(fib.posEndP(), fib.dirEndP(), width);
                break;
            case 3:
                gleCylinder(fib.posEndP(), fib.dirEndP(), width);
                break;
            case 4:
                gleArrowTail(fib.posEndP(), fib.dirEndP(), width);
                break;
            case 5:
                gleArrowTail(fib.posEndP(), -fib.dirEndP(), width);
                break;
        }
    }
}


void Display::displayFiberLinesMinus(Fiber const& fib, real abs, real width) const
{
    if ( abs < fib.abscissaM() )
        return;
    
    lineWidth(width);
    glBegin(GL_LINE_STRIP);
    unsigned ii = 0;
    real a = fib.abscissaM();
    while ( a < abs  &&  ii < fib.nbPoints() )
    {
        gleVertex(fib.posP(ii));
        a += fib.segmentation();
        ++ii;
    }
    if ( ii < fib.nbPoints() )
        gleVertex(fib.pos(abs));
    glEnd();
}


void Display::displayFiberLinesPlus(Fiber const& fib, real abs, real width) const
{
    if ( abs > fib.abscissaP() )
        return;
    
    lineWidth(width);
    glBegin(GL_LINE_STRIP);
    int ii = fib.lastPoint();
    real a = fib.abscissaP();
    while ( abs < a  &&  ii >= 0 )
    {
        gleVertex(fib.posP(ii));
        a -= fib.segmentation();
        --ii;
    }
    if ( ii >= 0 )
        gleVertex(fib.pos(abs));
    glEnd();
}


void Display::displayFiberLines(Fiber const& fib) const
{
    const FiberDisp * disp = fib.prop->disp;
    
    if ( disp->line_style == 1 )
    {
        // display plain lines:
        lineWidth(disp->line_width);
        glBegin(GL_LINE_STRIP);
        for ( unsigned ii = 0; ii < fib.nbPoints(); ++ii )
            gleVertex( fib.posP(ii) );
        glEnd();
    }
    else if ( disp->line_style == 2 )
    {
        // display segments with color indicating internal tension
        lineWidth(disp->line_width);
        const GLfloat alpha=disp->color.alpha();
        glBegin(GL_LINES);
        for ( unsigned ii = 0; ii < fib.lastPoint(); ++ii )
        {
            // the Lagrange multipliers are negative under compression
            jet_color(fib.tension(ii)/disp->tension, alpha).load();
            gleVertex( fib.posP(ii) );
            gleVertex( fib.posP(ii+1) );
        }
        glEnd();
    }
    else if ( disp->line_style == 3 )
    {
        // display segments with color indicating the curvature
        lineWidth(disp->line_width);
        const GLfloat alpha=disp->color.alpha();
        glBegin(GL_LINE_STRIP);
        if ( fib.nbPoints() > 2 )
            jet_color(fib.curvature(1), alpha).load();
        else
            jet_color(0, alpha).load();
        gleVertex(fib.posP(0));
        for ( unsigned ii = 1; ii < fib.lastPoint(); ++ii )
        {
            jet_color(fib.curvature(ii), alpha).load();
            gleVertex(fib.posP(ii));
        }
        gleVertex(fib.posP(fib.lastPoint()));
        glEnd();
    }
    else if ( disp->line_style == 4 )
    {
        // color according to the angle with respect to the XY-plane:
        lineWidth(disp->line_width);
        glBegin(GL_LINES);
        for ( unsigned ii = 0; ii < fib.lastPoint(); ++ii )
        {
            hue_color(fib.angleXY(ii), 1.0).load();
            gleVertex(fib.posP(ii));
            gleVertex(fib.posP(ii+1));
        }
        glEnd();
    }
}


void Display::displayFiberSpeckles(Fiber const& fib) const
{
    const FiberDisp * disp = fib.prop->disp;
    
    // display random speckles:
    if ( disp->speckle_style == 1 )
    {
        /*
         A simple random number generator seeded by fib.signature()
         is used to distribute points always at the same position
         with respect to the lattice of each fiber.
         */
        pointSize(disp->point_size);
        glBegin(GL_POINTS);
        
        const real spread = disp->speckle_interval;
        const real S = 0x1p-32;
        // draw speckles below the origin of abscissa
        if ( fib.abscissaM() < 0 )
        {
            uint32_t z = fib.signature();
            real a = spread * log(z*S);
            while ( a > fib.abscissaP() )
            {
                z = lcrng2(z);
                a += spread * log(z*S);
            }
            while ( a >= fib.abscissaM() )
            {
                gleVertex(fib.pos(a));
                z = lcrng2(z);
                a += spread * log(z*S);
            }
        }
        // draw speckles above the origin of abscissa
        if ( fib.abscissaP() > 0 )
        {
            uint32_t z = ~fib.signature();
            real a = -spread * log(z*S);
            while ( a < fib.abscissaM() )
            {
                z = lcrng1(z);
                a -= spread * log(z*S);
            }
            while ( a <= fib.abscissaP() )
            {
                gleVertex(fib.pos(a));
                z = lcrng1(z);
                a -= spread * log(z*S);
            }
        }
        glEnd();
    }
    else if ( disp->speckle_style == 2 )
    {
        // display regular speckles
        pointSize(disp->point_size);
        glBegin(GL_POINTS);
        //we distribute points regularly along the tube
        const real grad = disp->speckle_interval;
        real ab = grad * ceil( fib.abscissaM() / grad );
        while ( ab <= fib.abscissaP() ) {
            gleVertex( fib.pos(ab) );
            ab += grad;
        }
        glEnd();
    }
}


void Display::displayFiberPoints(Fiber const& fib) const
{
    const FiberDisp * disp = fib.prop->disp;

    if ( disp->point_style == 1 )
    {
        // display model-points:
        pointSize(disp->point_size);
        glBegin(GL_POINTS);
        for ( unsigned ii = 0; ii < fib.nbPoints(); ++ii )
            gleVertex( fib.posP(ii) );
        glEnd();
    }
    else if ( disp->point_style == 2 )
    {
        // display arrowheads along the fiber:
        real ab = ceil( fib.abscissaM() );
        const real grad = disp->point_interval;
        while ( ab <= fib.abscissaP() )
        {
            gleCone(fib.pos(ab), fib.dir(ab), 1.5*disp->point_size*sFactor);
            ab += grad;
        }
    }
    else if ( disp->point_style == 3 )
    {
        // display middle of fiber:
        pointSize(2*disp->point_size);
        glBegin(GL_POINTS);
            gleVertex( fib.posMiddle() );
        glEnd();
    }
}


void Display::displayFiberLattice(Fiber const& fib) const
{
    const FiberDisp * disp = fib.prop->disp;
    const gle_color col = fib.disp->color;
    
    FiberLattice const* lat = fib.lattice();
    assert_true(lat->ready());
    
    const real uni = lat->unit();
    const int inf = lat->index(fib.abscissaM());
    const int sup = lat->index(fib.abscissaP());
    assert_true( inf <= sup );
    const real fac = 1.0 / ( uni * disp->lattice_scale );

    /*
     This style uses one vertex for each site, positionned at the center of the range
     OpenGL will interpolate the colors, and each site will be covered by a gradient.
     */
    if ( disp->lattice_style & 1 )
    {
        lineWidth(disp->line_width);
        glBegin(GL_LINE_STRIP);
        if ( inf == sup )
        {
            //the Fiber is entirely covered by one site!
            real s = fib.abscissaP() - fib.abscissaM();
            col.load_darken(fac * lat->value(inf) * uni / s);
            gleVertex(fib.posEndM());
            gleVertex(fib.posEndP());
        }
        else
        {
            // the terminal site may be truncated
            real s = lat->abscissa(inf+1) - fib.abscissaM();
            if ( s > 0 )
                col.load_darken(fac * lat->value(inf) * uni / s);
            gleVertex(fib.posEndM());
            if ( uni*(inf+0.5) > fib.abscissaM() )
                gleVertex(fib.pos(uni*(inf+0.5)));

            for ( int h = inf+1; h < sup; ++h )
            {
                col.load_darken(fac * lat->value(h));
                gleVertex(fib.pos(uni*(h+0.5)));
            }

            // the terminal site may be truncated
            s = fib.abscissaP() - lat->abscissa(sup);
            if ( s > 0 )
                col.load_darken(fac * lat->value(sup) * uni / s);
            if ( uni*(sup+0.5) < fib.abscissaP() )
                gleVertex(fib.pos(uni*(sup+0.5)));
            gleVertex(fib.posEndP());
        }
        glEnd();
    }
    /*
     This style, uses two vertices for each site, positionned at the extremity of the range,
     and each site is entirely covered by the color corresponding to the value.
     */
    else if ( disp->lattice_style & 2 )
    {
        lineWidth(disp->line_width);
        glBegin(GL_LINE_STRIP);
        if ( inf == sup )
        {
            //the Fiber is entirely covered by one site!
            real s = fib.abscissaP() - fib.abscissaM();
            col.load_darken(fac * lat->value(inf) * uni / s);
            gleVertex(fib.posEndM());
            gleVertex(fib.posEndP());
        }
        else
        {
            // the terminal site may be truncated
            real s = lat->abscissa(inf+1) - fib.abscissaM();
            if ( s > 0 )
            col.load_darken(fac * lat->value(inf) * uni / s);
            gleVertex(fib.posEndM());
            
            for ( int h = inf+1; h < sup; ++h )
            {
                gleVertex(fib.pos(uni*h));
                col.load_darken(fac * lat->value(h));
                gleVertex(fib.pos(uni*h));
            }
            
            // the terminal site may be truncated
            s = fib.abscissaP() - lat->abscissa(sup);
            if ( s > 0 )
            col.load_darken(fac * lat->value(sup) * uni / s);
            gleVertex(fib.pos(uni*sup));
            gleVertex(fib.posEndP());
        }
        glEnd();
    }
    
    /*
     Indicate the edges between sites with small white dots
     */
    if ( disp->lattice_style & 4 )
    {
        disp->inner_color.load();
        glPointSize(2);
        glBegin(GL_LINES);
        for ( int h = inf+1; h <= sup; h+=8 )
        {
            gleVertex(fib.pos(uni*h));
            gleVertex(fib.pos(uni*h)+ Vector(uni,0,0));
        }
        glEnd();
    }
}


void Display::displayFiberLabels(Fiber const& fib) const
{
    const FiberDisp * disp = fib.prop->disp;

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    if ( disp->label_style == 1 )
    {
        // indicate fiber reference at every model point
        std::string str = "f"+sMath::repr(fib.identity()); //fib.reference();
        for ( unsigned ii = 0; ii < fib.nbPoints(); ++ii )
            gleDrawText(fib.posP(ii), str.c_str(), GLUT_BITMAP_8_BY_13);
    }
    else if ( disp->label_style == 2 )
    {
        // indicate model-points with numbers
        char tmp[16];
        for ( unsigned ii = 0; ii < fib.nbPoints(); ++ii )
        {
            snprintf(tmp, sizeof(tmp), "%i", ii);
            gleDrawText(fib.posP(ii), tmp, GLUT_BITMAP_8_BY_13);
        }
    }
    else if ( disp->label_style == 3 )
    {
        // display integral abscissa along the fiber
        char tmp[16];
        snprintf(tmp, sizeof(tmp), "%.2f", fib.abscissaM());
        gleDrawText(fib.posEndM(), tmp, GLUT_BITMAP_HELVETICA_10);
        
        int s = ceil( fib.abscissaM() );
        int e = floor( fib.abscissaP() );
        for ( int a = s; a <= e; ++a )
        {
            snprintf(tmp, sizeof(tmp), "%i", a);
            gleDrawText(fib.pos(a), tmp, GLUT_BITMAP_HELVETICA_10);
        }
        
        snprintf(tmp, sizeof(tmp), "%.2f", fib.abscissaP());
        gleDrawText(fib.posEndP(), tmp, GLUT_BITMAP_HELVETICA_10);
    }
    glPopAttrib();
}


// display forces acting on the points
void Display::displayFiberForces(Fiber const& fib) const
{
    const FiberDisp * disp = fib.prop->disp;
    
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    lineWidth(disp->point_size);
    glBegin(GL_LINES);
    for ( unsigned ii = 0; ii < fib.nbPoints(); ++ii )
    {
        gleVertex( fib.posP(ii) );
        gleVertex( fib.posP(ii) + disp->forces * fib.netForce(ii) );
    }
    glEnd();
    glPopAttrib();
}

#define DRAW_SPHERE gleSphere2B

/**
 This renders 26 spheres positionned on a right-handed helix,
 making one turn every 74nm, with a max width of ~ 9nm.
 This is roughly Ken Holmes' model of F-actin:
 Nature 347, 44 - 49 (06 September 1990); doi:10.1038/347044a0
 which shows half a turn in 37nm containing 13 monomers.
 */
void Display::displayActin(Fiber const& fib,
                           gle_color const& color1,
                           gle_color const& color2,
                           gle_color const& colorE) const
{    
    // axial translation between two sucessive monomers:
    const real dab = 0.00275;
    // enlarge radius of monomers to make them overlap
    const real rad = 1.3 * dab;
    // distance from central axis to center of monomers
    real off = 0.0045 - dab;
    
    // rotation angle between consecutive monomers
    const real dan = 166 * M_PI / 180;
    const real cs = cos(dan);
    const real sn = sin(dan);

    real ab = 0;
    Vector3 d(fib.dirEndM(), DIM);   // unit tangent to centerline
    Vector3 n(fib.normal, DIM);      // normal direction
    if ( n.normSqr() < 0.5 )
        n = d.orthogonal(1.0);
    n = d.orthogonal(n, 1.0);
    fib.normal.get(n, DIM);
    //std::clog << fib.reference() << " " << n << "    " << n.normSqr() << '\n';
    
    Vector3 p, q;
    
    glEnable(GL_CLIP_PLANE4);
    
    int cnt = 0;
    // rotate until we reach the MINUS_END
    while ( ab <= fib.abscissaM() )
    {
        ++cnt;
        ab += dab;
        n = d.rotate(n, cs, sn);
    }
    p = Vector3(fib.pos(ab), DIM) + off * n;
    // draw the monomers until the PLUS_END:
    while ( ab < fib.abscissaP() )
    {
        q = p;
        ab += dab;
        d = Vector3(fib.dir(ab), DIM);
        n = d.rotate(n, cs, sn);
        p = Vector3(fib.pos(ab), DIM) + off * n;
        
        // use different tones to individualize the two strands:
        if ( ++cnt & 1 )
            color1.load_load();
        else
            color2.load_load();

        // change color for the last monomer:
        if ( ab + dab > fib.abscissaP() )
        {
            colorE.load_load();
            glDisable(GL_CLIP_PLANE4);
        }
        
        // set clipping plane with the next monomer
        setClipPlane(GL_CLIP_PLANE4, (q-p).normalized(), (p+q)*0.5);
        
        gleObject(q, rad, DRAW_SPHERE);
        
        // set cliping plane with the previous:
        setClipPlane(GL_CLIP_PLANE5, (p-q).normalized(), (p+q)*0.5);
        
        glEnable(GL_CLIP_PLANE5);
    }
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
}


/**
 This renders a Microtubule using sphere of alternating colors
 */
void Display::displayMicrotubule(Fiber const& fib,
                                 gle_color const& color1,
                                 gle_color const& color2,
                                 gle_color const& colorE) const
{
    real da[] = {0,0.000923,0.001846,0.002769,0.003692,0.004615,0.005538,0.006461,0.007384,0.008308,0.009231,0.010154,0.011077};
    real dx[] = {0.8855,0.5681,0.1205,-0.3546,-0.7485,-0.9709,-0.9709,-0.7485,-0.3546,0.1205,0.5681,0.8855,1.0000};
    real dy[] = {-0.4647,-0.8230,-0.9927,-0.9350,-0.6631,-0.2393,0.2393,0.6631,0.9350,0.9927,0.8230,0.4647,0};

    // axial translation between two sucessive monomers:
    const real dab = 0.004;
    // enlarged radius of monomers makes them overlap slighlty
    const real rad = 0.7 * dab;
    // distance from central axis to center of monomers
    real off = 0.025 / 2 - rad;
    
    real ab = dab * ceil( fib.abscissaM() / dab );
    Vector3 d(fib.dir(ab), DIM);   // unit tangent vector
    Vector3 n(fib.normal, DIM);    // normal direction

    // adjust normal direction:
    if ( n.normSqr() < 0.5 )
        n = d.orthogonal(1.0);
    n = d.orthogonal(n, 1.0);
    fib.normal.get(n, DIM);
    
    const real abmax = fib.abscissaP();

    int cnt = 0;
    while ( ab < abmax )
    {
        d = Vector3(fib.dir(ab), DIM);
        Vector3 p(fib.pos(ab), DIM);
        
        // adjust 'n' to keep it orthogonal to 'd':
        n = d.orthogonal(n, 1.0);
        
        // set two vectors orthogonal to 'd' of norm 'off':
        Vector3 e = n * off;
        Vector3 f = cross(d, e);

        // color of alpha subunit:
        if ( 0 == ( cnt & 1 ) )
            color1.load_load();
       
        real up = abmax - ( cnt & 1 ? 0 : dab );
        
        for ( int i = 0; i < 13; ++i )
        {
            if ( ab + da[i] < up )
            {
                if ( cnt & 1 )
                {
                    // change color for beta subunit:
                    if ( ab + da[i] + 4 * dab > abmax )
                        colorE.load_load();
                    else
                        color2.load_load();
                }
                gleObject(p+dx[i]*e+dy[i]*f+da[i]*d, rad, DRAW_SPHERE);
            }
        }

        ab += dab;
        ++cnt;
    }
}


void Display::displayFiber(Fiber const& fib)
{
    const FiberDisp * disp = fib.prop->disp;
    disp->inner_color.load_back();
    
#if ( DIM == 3 )
    if ( fib.disp->color.transparent() )
    {
        zObjects.push_back(zObject(&fib, fib.nbPoints()/2));
    }
    else
#endif
    {
        if ( fib.lattice()  &&  disp->lattice_style )
        {
            fib.disp->color.load_load();
            displayFiberLines(fib);
            displayFiberLattice(fib);
        }
        else if ( disp->line_style )
        {
            fib.disp->color.load_load();

            if ( disp->line_style == 1 && disp->style == 1 )
                displayActin(fib, fib.disp->color, fib.disp->color.darken(0.625), fib.disp->end_color[0]);
            else if ( disp->line_style == 1 && disp->style == 2 )
                displayMicrotubule(fib, fib.disp->color, fib.disp->color.darken(0.625), fib.disp->end_color[0]);
            else
                displayFiberLines(fib);
        }
    }

    if ( disp->end_length[1] > 0 )
    {
        fib.disp->end_color[1].load_load();
        displayFiberLinesMinus(fib, fib.abscissaM()+disp->end_length[1], disp->end_size[1]);
    }
    
    if ( disp->end_length[0] > 0 )
    {
        fib.disp->end_color[0].load_load();
        displayFiberLinesPlus(fib, fib.abscissaP()-disp->end_length[0], disp->end_size[0]);
    }

    if ( disp->point_style > 0 )
    {
        fib.disp->color.load_load();
        displayFiberPoints(fib);
    }
    
    if ( disp->speckle_style > 0 )
    {
        fib.disp->color.load_load();
        displayFiberSpeckles(fib);
    }

    // draw other fiber elements only if fiber is fully visible:
    if ( fib.disp->visible > 0 )
    {
        if ( disp->label_style )
        {
            fib.disp->color.alpha_set(0.5).load();
            displayFiberLabels(fib);
        }
        
        if ( disp->end_style[1] )
        {
            fib.disp->end_color[1].load_load();
            displayFiberEndMinus(fib);
        }
        
        if ( disp->end_style[0] )
        {
            fib.disp->end_color[0].load_load();
            displayFiberEndPlus(fib);
        }
        
        if ( disp->forces )
        {
            disp->forces_color.load();
            displayFiberForces(fib);
        }
    }
}


void Display::displayFibers(FiberSet const& set)
{
#if ( 0 )
    // display Fibers in a random (ever changing) order:
    for ( Fiber * fib = set.first(); fib ; fib=fib->next() )
    {
        if ( fib->disp->visible )
            displayFiber(*fib);
    }
#else
    // display the Fiber always in the same order:
    Fiber * fib = set.firstID();
    while ( fib )
    {
        if ( fib->disp->visible )
            displayFiber(*fib);
        fib = set.nextID(fib);
    }
#endif
}



//------------------------------------------------------------------------------
#pragma mark -

void Display::displayCouplesF(CoupleSet const& set)
{
    if ( prop->couple_flip )
        displayCouplesF2(set);
    else
        displayCouplesF1(set);
}


void Display::displaySolids(SolidSet const& set)
{
    for ( Solid * obj = set.first(); obj; obj=obj->next() )
    {
        if ( obj->prop->disp->visible )
        {
            displaySolid(*obj);
#if ( DIM == 3 )
            if ( obj->prop->disp->color.transparent() )
            {
                for ( unsigned ii = 0; ii < obj->nbPoints(); ++ii )
                    zObjects.push_back(zObject(obj, ii));
            }
            else
#endif
            {
                for ( unsigned ii = 0; ii < obj->nbPoints(); ++ii )
                    displaySolidT(*obj, ii);
            }
        }
    }
}


void Display::displayBeads(BeadSet const& set)
{
    for ( Bead * obj = set.first(); obj; obj=obj->next() )
    {
        if ( obj->prop->disp->visible )
        {
            displayBead(*obj);
#if ( DIM == 3 )
            if ( obj->prop->disp->color.transparent() )
                zObjects.push_back(zObject(obj));
            else
#endif
                displayBeadT(*obj);
        }
    }
}


void Display::displaySpheres(SphereSet const& set)
{
    for ( Sphere * obj=set.first(); obj ; obj=obj->next() )
    {
        if ( obj->prop->disp->visible )
        {
            displaySphere(*obj);
#if ( DIM == 3 )
            if ( obj->prop->disp->color.transparent() )
                zObjects.push_back(zObject(obj));
            else
#endif
                displaySphereT(*obj);
        }
    }
}


void Display::displayOrganizers(OrganizerSet const& set)
{
    for ( Organizer * obj=set.first(); obj ; obj=obj->next() )
        displayOrganizer(*obj);
}


//------------------------------------------------------------------------------
#pragma mark -


/// display sub-part `inx` of object `obj`
void Display::zObject::display(Display * disp)
{
    switch( mObj->tag() )
    {
        case Fiber::TAG: {
            //\todo we should depth-sort segments of the fibers independently
            Fiber const* fib = static_cast<const Fiber*>(mObj);
            fib->disp->color.load_both();
            disp->displayFiberLines(*fib);
        } break;
            
        case Solid::TAG:
            disp->displaySolidT(*static_cast<const Solid*>(mObj), mIdx);
            break;
            
        case Bead::TAG:
            disp->displayBeadT(*static_cast<const Bead*>(mObj));
            break;
            
        case Sphere::TAG:
            disp->displaySphereT(*static_cast<const Sphere*>(mObj));
            break;
            
        default:
            std::cerr << "Internal error: unknown zObject pointer" << std::endl;
    }
}


/**
 This display objects in `zObjects` from back to front

 Depth-sorting is used in 3D to display transparent objects
 from the furthest to the nearest.
*/
void Display::displayTransparentObjects(Array<zObject>& list)
{
#if ( DIM == 3 )
    
    // get current modelview transformation:
    GLfloat mat[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, mat);
    
    // extract axis corresponding to vertical direction:
    Vector3 vertical(mat[2], mat[6], mat[10]);
    
    for ( zObject * i = list.begin(); i < list.end(); ++i )
        i->depth( i->position() * vertical );
    
    // depth-sort objects:
    qsort(list.data(), list.size(), sizeof(zObject), &closer);

    /*
     Enable polygon offset to avoid artifacts with objects of same size,
     particularly the ends of filaments with their shaft.
     */
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);

    for ( zObject * i = list.begin(); i < list.end(); ++i )
        i->display(this);
    
    glDisable(GL_POLYGON_OFFSET_FILL);

#endif
}


