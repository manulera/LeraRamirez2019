// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "simul.h"
#include "display3.h"
#include "modulo.h"

#include "fake.h"
#include "fiber_disp.h"
#include "line_disp.h"
#include "point_disp.h"

#include "opengl.h"
#include "gle.h"
#include "gle_color_list.h"
#include "glapp.h"
#include "glut.h"

using namespace gle;
extern Modulo const* modulo;


Display3::Display3(DisplayProp const* dp) : Display(dp)
{
}

//------------------------------------------------------------------------------
#pragma mark -

void Display3::displayObjects(Simul const& sim)
{
    glDepthMask(GL_FALSE);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    displayFields(sim.fields);
    
    glEnable(GL_LIGHTING);
    glDepthMask(GL_TRUE);
    
    View& view = glApp::currentView();

    if ( DIM == 3 && view.stencil )
    {
        /*
         We use here the stencil test to make sure that nothing else is drawn
         where the inner side of the fibers is visible. This improves the
         display with clipping planes, as fibers appear as cut solid objects
         */
        glClearStencil(0);
        glClear(GL_STENCIL_BUFFER_BIT);
        glEnable(GL_STENCIL_TEST);
        glStencilFunc(GL_ALWAYS, 1, ~0);
        
        // draw inner surfaces of fibers:
        glEnable(GL_CULL_FACE);
        glCullFace(GL_FRONT);
        
        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
        displayFibers(sim.fibers);
        
        // draw front surfaces of fibers:
        glCullFace(GL_BACK);
        glStencilOp(GL_KEEP, GL_KEEP, GL_ZERO);
        displayFibers(sim.fibers);
        
        glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
        glStencilFunc(GL_EQUAL, 0, ~0);
    }
    else
    {
        glDisable(GL_CULL_FACE);
        displayFibers(sim.fibers);
    }
    
    glEnable(GL_LIGHTING);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    view.inner_color.load_back();

    displayBeads(sim.beads);
    displaySolids(sim.solids);
    displaySpheres(sim.spheres);
    
    if ( prop->single_select & 1 )
        displaySinglesF(sim.singles);
    
    if ( prop->couple_select & 1 )
        displayCouplesF(sim.couples);

    if ( prop->couple_select & 2 )
        displayCouplesA(sim.couples);

    if ( prop->couple_select & 4 )
        displayCouplesB(sim.couples);

    if ( prop->single_select & 2 )
        displaySinglesA(sim.singles);
    
    if ( DIM == 3 && view.stencil )
    {
        glClearStencil(0);
        glDisable(GL_STENCIL_TEST);
    }

    displayOrganizers(sim.organizers);
    displayMisc(sim);
}


//------------------------------------------------------------------------------
#pragma mark -

void Display3::displayBall(Vector const& pos, real radius)
{
    assert_true(glIsEnabled(GL_CULL_FACE));
    assert_true(glIsEnabled(GL_LIGHTING));
    glPushMatrix();
    gleTranslate(pos);
    gleScale(radius);
    glCullFace(GL_FRONT);
    gleSphere4B();
    glCullFace(GL_BACK);
    gleSphere8B();
    glPopMatrix();
}


void Display3::displayPoint(Vector const& pos, PointDisp const* disp)
{
    if ( disp->perceptible )
    {
        glPushMatrix();
        gleTranslate(pos);
        gleScale(disp->size*sFactor);
        gleSphere2B();
        
#if ( 0 )
        if ( disp->symbol )
        {
            glDisable(GL_LIGHTING);
            disp->symbol_color.load();
            glRasterPos2f(0,0);
            glBitmap(0,0,0,0,-5,-4,0);
            glutBitmapCharacter(GLUT_BITMAP_9_BY_15, disp->symbol);
            glEnable(GL_LIGHTING);
        }
#endif
        glPopMatrix();
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void Display3::displayFiberLinesMinus(Fiber const& fib, real abs, real width) const
{
    if ( abs < fib.abscissaM() )
        return;
    
    const real rad = width * sFactor;
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(-1.0, 1.0);
    
    //close the MINUS_END with a sphere or a disc
    if ( fib.prop->disp->line_caps == 2 )
    {
        glEnable(GL_CLIP_PLANE4);
        setClipPlane(GL_CLIP_PLANE4, -fib.dirEndM(), fib.posEndM());
        gleObject(fib.posEndM(), rad, gleSphere2);
        glDisable(GL_CLIP_PLANE4);
    }
    else
        gleObject(fib.posEndM(), -fib.dirEndM(), rad, gleCircleSB);

    unsigned ii = 0;
    real a = fib.abscissaM() + fib.segmentation();
    while ( a < abs  &&  ii < fib.nbPoints() )
    {
        gleTube(fib.posP(ii), fib.posP(ii+1), rad, gleTube2B);
        a += fib.segmentation();
        ++ii;
    }
    
    if ( ii < fib.nbPoints() )
    {
        gleTube(fib.posP(ii), fib.pos(abs), rad, gleTube2B);
        //close with a disc
        glDisable(GL_POLYGON_OFFSET_FILL);
        gleObject(fib.pos(abs),  fib.dir(abs), rad, gleCircleSB);
    }
    else
    {
        gleObject(fib.posEndP(),  fib.dirEndP(), rad, gleCircleSB);
        glDisable(GL_POLYGON_OFFSET_FILL);
    }
}


void Display3::displayFiberLinesPlus(Fiber const& fib, real abs, real width) const
{
    if ( abs > fib.abscissaP() )
        return;

    const real rad = width * sFactor;
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(-1.0, 1.0);
    
    //close the PLUS_END with a sphere or a disc
    if ( fib.prop->disp->line_caps == 2 )
    {
        glEnable(GL_CLIP_PLANE4);
        setClipPlane(GL_CLIP_PLANE4, fib.dirEndP(), fib.posEndP());
        gleObject(fib.posEndP(), rad, gleSphere2);
        glDisable(GL_CLIP_PLANE4);
    }
    else
        gleObject(fib.posEndP(), fib.dirEndP(), rad, gleCircleSB);

    int ii = fib.lastSegment();
    real a = fib.abscissaP() - fib.segmentation();
    while ( abs < a  &&  0 < ii )
    {
        gleTube(fib.posP(ii), fib.posP(ii+1), rad, gleTube2B);
        a -= fib.segmentation();
        --ii;
    }

    if ( 0 <= ii )
    {
        gleTube(fib.pos(abs), fib.posP(ii+1), rad, gleTube2B);
        //close with a disc:
        glDisable(GL_POLYGON_OFFSET_FILL);
        gleObject(fib.pos(abs), -fib.dir(abs), rad, gleCircleSB);
    }
    else
    {
        gleObject(fib.posEndM(), -fib.dirEndM(), rad, gleCircleSB);
        glDisable(GL_POLYGON_OFFSET_FILL);
    }
}


void Display3::displayFiberLines(Fiber const& fib) const
{
    const FiberDisp * disp = fib.prop->disp;
    // radius of lines and points in space units:
    real rad = disp->line_width * sFactor;
    
    if ( disp->line_style == 1 )
    {
#if ( DIM > 1 )
        // this code makes nice joints for bent tubes, by using CLIP_PLANES
        
        // draw MINUS_END
        Vector dir = fib.dirEndM();
        Vector pos = fib.posEndM();
        
        if ( disp->line_caps == 1 )
            gleObject(pos, -dir, rad, gleCircleSB);

        glEnable(GL_CLIP_PLANE4);
        if ( disp->line_caps == 2 )
        {
            setClipPlane(GL_CLIP_PLANE4, -dir, pos);
            gleObject(pos, rad, gleSphere2);
        }

        setClipPlane(GL_CLIP_PLANE4, dir, pos);
        glEnable(GL_CLIP_PLANE5);
        
        // draw inner segments
        for ( unsigned ii = 0; ii < fib.lastSegment(); ++ii )
        {
            pos = fib.posP(ii+1);
            dir = ( fib.diffPoints(ii) + fib.diffPoints(ii+1) ).normalized();
            setClipPlane(GL_CLIP_PLANE5, -dir, pos);
            gleTube(fib.posP(ii), pos, rad, gleLongTube2B);
            setClipPlane(GL_CLIP_PLANE4,  dir, pos);
        }
        
        // draw last segment:
        dir = fib.dirEndP();
        pos = fib.posEndP();
        setClipPlane(GL_CLIP_PLANE5, -dir, pos);
        gleTube(fib.posP(fib.lastSegment()), pos, rad, gleLongTube2B);
        
        glDisable(GL_CLIP_PLANE4);
        // draw PLUS_END:
        if ( disp->line_caps == 2 )
        {
            setClipPlane(GL_CLIP_PLANE5, dir, pos);
            gleObject(pos, rad, gleSphere2);
        }
        
        glDisable(GL_CLIP_PLANE5);
        if ( disp->line_caps == 1 )
            gleObject(pos, dir, rad, gleCircleSB);
#else
        for ( unsigned ii = 0; ii < fib.nbSegments(); ++ii )
            gleTube(fib.posP(ii), fib.posP(ii+1), rad, gleTube2B);
#endif
    }
    else if ( disp->line_style == 2 )
    {
        const GLfloat alpha=disp->color.alpha();
#if ( DIM > 1 )
        // display internal tensions
        jet_color(fib.tension(0)/disp->tension, alpha).load_front();
        // draw MINUS_END
        Vector dir = fib.dirEndM();
        Vector pos = fib.posEndM();
        
        if ( disp->line_caps == 1 )
            gleObject(pos, -dir, rad, gleCircleSB);
        
        glEnable(GL_CLIP_PLANE4);
        if ( disp->line_caps == 2 )
        {
            setClipPlane(GL_CLIP_PLANE4, -dir, pos);
            gleObject(pos, rad, gleSphere2);
        }
        
        setClipPlane(GL_CLIP_PLANE4, dir, pos);
        glEnable(GL_CLIP_PLANE5);
        
        for ( unsigned ii = 0; ii < fib.lastSegment(); ++ii )
        {
            pos = fib.posP(ii+1);
            dir = ( fib.diffPoints(ii) + fib.diffPoints(ii+1) ).normalized();
            setClipPlane(GL_CLIP_PLANE5, -dir, pos);
            jet_color(fib.tension(ii)/disp->tension, alpha).load_front();
            gleTube(fib.posP(ii), pos, rad, gleLongTube2B);
            setClipPlane(GL_CLIP_PLANE4,  dir, pos);
        }
        // draw last segment:
        dir = fib.dirEndP();
        pos = fib.posEndP();
        setClipPlane(GL_CLIP_PLANE5, -dir, pos);
        jet_color(fib.tension(fib.lastSegment())/disp->tension, alpha).load_front();
        gleTube(fib.posP(fib.lastSegment()), pos, rad, gleLongTube2B);
        
        glDisable(GL_CLIP_PLANE4);
        // draw PLUS_END:
        if ( disp->line_caps == 2 )
        {
            setClipPlane(GL_CLIP_PLANE5, dir, pos);
            gleObject(pos, rad, gleSphere2);
        }
        
        glDisable(GL_CLIP_PLANE5);
        if ( disp->line_caps == 1 )
            gleObject(pos, dir, rad, gleCircleSB);
#else
        for ( unsigned ii = 0; ii < fib.nbSegments(); ++ii )
        {
            // the Lagrange multipliers are negative under compression
            jet_color(fib.tension(ii)/disp->tension, alpha).load_front();
            gleTube(fib.posP(ii), fib.posP(ii+1), rad, gleTube2B);
        }
#endif
    }
    else if ( disp->line_style == 3 )
    {
        // display segments with color indicating the curvature
        gle_color col1, col2;
        if ( fib.nbPoints() > 2 )
            col2 = jet_color(fib.curvature(1), 1.0);
        else
            col2 = jet_color(0, 1);
        for ( unsigned ii = 0; ii < fib.lastSegment(); ++ii )
        {
            col1 = col2;
            col2 = jet_color(fib.curvature(ii+1), 1.0);
            gleTube(fib.posP(ii), rad, col1, fib.posP(ii+1), rad, col2);
        }
        unsigned ii = fib.lastSegment();
        gleTube(fib.posP(ii), rad, col1, fib.posP(ii+1), rad, col2);
    }
    else if ( disp->line_style == 4 )
    {
        // color according to the angle with respect to the XY-plane:
        for ( unsigned ii = 0; ii < fib.nbSegments(); ++ii )
        {
            hue_color(fib.angleXY(ii), 1.0).load_front();
            gleTube(fib.posP(ii), fib.posP(ii+1), rad, gleTube2B);
        }
    }
}


/**
 Display the MINUS_END of a Fiber, according to `style`:
 - 1: draw a sphere
 - 2: draw a cone
 - 3: draw a flat cylinder
 - 4: draw an arrow-head
 - 5: arrow-head in reverse direction
 .
 with 3D objects
 */
void Display3::displayFiberEndMinus(Fiber const& fib) const
{
    const int IM = 1;
    const FiberDisp * disp = fib.prop->disp;
    
    real width = disp->end_size[IM]*sFactor;
    if ( width > 0 )
    {
        switch (disp->end_style[IM])
        {
            case 1:
                gleObject(fib.posP(0), width, gleSphere2B);
                break;
            case 2:
                gleObject(fib.posP(0), -fib.dirPoint(0), width, gleCone1B);
                break;
            case 3:
                gleObject(fib.posP(0), -fib.dirPoint(0), width, gleCylinderB);
                break;
            case 4:
                gleObject(fib.posP(0),  fib.dirPoint(0), width, gleArrowTail2B);
                break;
            case 5:
                gleObject(fib.posP(0), -fib.dirPoint(0), width, gleArrowTail2B);
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
 with 3D objects
 */
void Display3::displayFiberEndPlus(Fiber const& fib) const
{
    const int IP = 0;
    const FiberDisp * disp = fib.prop->disp;
    
    real width = disp->end_size[IP]*sFactor;
    if ( width > 0 )
    {
        switch(disp->end_style[IP])
        {
            case 1:
                gleObject(fib.posEndP(), width, gleSphere2B);
                break;
            case 2:
                gleObject(fib.posEndP(), fib.dirEndP(), width, gleCone1B);
                break;
            case 3:
                gleObject(fib.posEndP(), fib.dirEndP(), width, gleCylinderB);
                break;
            case 4:
                gleObject(fib.posEndP(), fib.dirEndP(), width, gleArrowTail2B);
                break;
            case 5:
                gleObject(fib.posEndP(), -fib.dirEndP(), width, gleArrowTail2B);
                break;
        }
    }
}


void Display3::displayFiberLattice(Fiber const& fib) const
{
    const FiberDisp * disp = fib.prop->disp;
    // radius of lines and points in space units:
    const real rad = disp->line_width * sFactor;
    const gle_color col = fib.disp->color;

    FiberLattice const* lat = fib.lattice();
    const real fac = 1.0 / ( lat->unit() *  disp->lattice_scale );

    glPushAttrib(GL_LIGHTING_BIT|GL_ENABLE_BIT);
    GLfloat mat[] = { 0.0, 0.0, 0.0, 1.0 };
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION,  mat);
    glMateriali (GL_FRONT_AND_BACK, GL_SHININESS, 32);
    
    const real uni = lat->unit();
    int   ix = lat->index(fib.abscissaM());
    real abs = lat->abscissa(ix+1);
    real sup = fib.abscissaP() - uni;
    
    // correct for the length of the truncated section:
    real s = abs - fib.abscissaM();
    if ( s > 0 )
    {
        col.darken(fac * lat->value(ix) * uni / s).load_both();
        gleTube(fib.posEndM(), fib.pos(abs), rad, gleTube2B);
    }
    while ( abs < sup )
    {
        ++ix;
        col.darken(fac * lat->value(ix)).load_both();
        gleTube(fib.pos(abs), fib.pos(abs+uni), rad, gleTube2B);
        abs += uni;
    }
    // correct for the length of the truncated section:
    s = fib.abscissaP() - abs;
    if ( s > 0 )
    {
        col.darken(fac * lat->value(ix+1) * uni / s).load_both();
        gleTube(fib.pos(abs), fib.posEndP(), rad, gleTube2B);
    }
    glPopAttrib();
}



void Display3::displayFiberSpeckles(Fiber const& fib) const
{
    const FiberDisp * disp = fib.prop->disp;
    // diameter of lines and points in space units:
    real rad = disp->point_size * sFactor;
    
    if ( disp->point_size * uFactor < 2 )
        return;

    // display random speckles:
    if ( disp->speckle_style == 1 )
    {
        /*
         A simple random number generator seeded by fib.signature()
         is used to distribute points always at the same position
         with respect to the lattice of each fiber.
         */
        
        const real spread = disp->speckle_interval;
        const real S = 0x1p-32;
        // draw speckles below the origin of abscissa:
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
                gleObject(fib.pos(a), rad, gleSphere1B);
                z = lcrng2(z);
                a += spread * log(z*S);
            }
        }
        // draw speckles above the origin of abscissa:
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
                gleObject(fib.pos(a), rad, gleSphere1B);
                z = lcrng1(z);
                a -= spread * log(z*S);
            }
        }
    }
    else if ( disp->speckle_style == 2 )
    {
        //we distribute points regularly along the tube
        const real grad = disp->speckle_interval;
        real ab = grad * ceil( fib.abscissaM() / grad );
        while ( ab <= fib.abscissaP() )
        {
            gleObject(fib.pos(ab), rad, gleSphere1B);
            ab += grad;
        }
    }
}


void Display3::displayFiberPoints(Fiber const& fib) const
{
    const FiberDisp * disp = fib.prop->disp;
    // diameter of lines and points in space units:
    real rad = disp->point_size * sFactor;
    
    if ( disp->point_size * uFactor < 2 )
        return;

    if ( disp->point_style == 1 )
    {
        // display model-points:
        for ( unsigned ii = 0; ii < fib.nbPoints(); ++ii )
            gleObject(fib.posP(ii), rad, gleSphere2B);
    }
    else if ( disp->point_style == 2 )
    {
        // display an arrow-head every micro-meter:
        real ab = ceil( fib.abscissaM() );
        const real grad = disp->point_interval;
        while ( ab <= fib.abscissaP() )
        {
            gleCone(fib.pos(ab), fib.dir(ab), 1.5*rad);
            ab += grad;
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark -

void Display3::displayBead(Bead const& obj)
{
    const PointDisp * disp = obj.prop->disp;

    // display center:
    if ( disp->style & 2 )
    {
        bodyColor(disp, obj.signature());
        displayPoint(obj.position(), disp);
    }
    
#if ( DIM == 2 )
    // display outline:
    if ( disp->style & 4 )
    {
        bodyColor(disp, obj.signature());
        lineWidth(disp->width);
        gleObject(obj.position(), obj.radius(), gleCircleLB);
    }
#endif
}

/**
 Display a bead as a sphere
 */
void Display3::displayBeadT(Bead const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    if ( disp->style & 1 )
    {
        bodyColorT(disp, obj.signature());
        displayBall(obj.position(), obj.radius());
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void Display3::displaySolid(Solid const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    //display points:
    if ( disp->style & 2  &&  disp->size > 0 )
    {
        bodyColor(disp, obj.signature());
        for ( unsigned ii = 0; ii < obj.nbPoints(); ++ii )
            displayPoint(obj.posP(ii), disp);
    }
    
#if ( DIM == 3 )
    //special display for ParM simulations (DYCHE)
    if ( obj.mark()  &&  disp->style & 4  &&  obj.nbPoints() >= 3 )
    {
        bodyColor(disp, obj.signature());
        gleObject(obj.posP(0), obj.diffPoints(1, 0), obj.radius(0), gleCircleLB);
    }
#endif
    
    //display a signature for each Solid
    if ( disp->style & 8 )
    {
        char tmp[8];
        snprintf(tmp, sizeof(tmp), "%lu", obj.identity());
        bodyColor(disp, obj.signature());
        gleDrawText(obj.posP(0), tmp, GLUT_BITMAP_HELVETICA_10);
    }
}


/**
 Display a semi-transparent disc / sphere
 */
void Display3::displaySolidT(Solid const& obj, unsigned int ii)
{
    const PointDisp * disp = obj.prop->disp;

    if ( disp->style & 1  &&  obj.radius(ii) > 0 )
    {
        bodyColorT(disp, obj.signature());
        displayBall(obj.posP(ii), obj.radius(ii));
    }
}


//------------------------------------------------------------------------------
#pragma mark -

void Display3::displaySphere(Sphere const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    //display center and surface points
    if ( disp->size > 0  &&  disp->style & 2 )
    {
        bodyColor(disp, obj.signature());
        displayPoint(obj.posP(0), disp);
        for ( unsigned ii = obj.nbRefPoints(); ii < obj.nbPoints(); ++ii )
            displayPoint(obj.posP(ii), disp);
    }

    //display reference points
    if ( disp->size > 0  &&  disp->style & 8 )
    {
        bodyColor(disp, obj.signature());
        for ( unsigned ii = 1; ii < obj.nbRefPoints(); ii++ )
            displayPoint(obj.posP(ii), disp);
    }
}

void Display3::displaySphereT(Sphere const& obj)
{
    const PointDisp * disp = obj.prop->disp;

    //display the envelope
    if ( disp->style & 5 )
    {
        bodyColorT(disp, obj.signature());
        lineWidth(disp->width);
        
        glPushMatrix();
        
#if ( DIM < 3 )
        
        gleTranslate(obj.posP(0));
        glutSolidTorus(disp->size*sFactor, obj.radius(), 2*gle::finesse, 10*gle::finesse);
        
#else
        
        /* Note: The rotation matrix for the sphere calculated below from the
            reference points, includes scaling by the radius of the sphere.
            We then use a primitive for a sphere of radius 1.
            */
        const Vector C = obj.posP(0);
        gleTransRotate(obj.posP(1)-C, obj.posP(2)-C, obj.posP(3)-C, C);
        
        if ( disp->style & 1 )
        {
            glCullFace(GL_FRONT);
            gleSphere4B();
            glCullFace(GL_BACK);
            gleSphere8B();
        }
        if ( disp->style & 4 )
        {
            disp->color2.load_front();
            gleThreeBands(64);
        }
        
#endif
        glPopMatrix();
    }
}

//------------------------------------------------------------------------------

void Display3::displayOrganizer(Organizer const& obj)
{
    const PointDisp * disp = obj.disp();
    
    if ( !disp )
        return;

    const real w = disp->width*sFactor;

    if ( disp->style & 2 )
    {
        Vector P, Q;
        unsigned inx = 0;
        bodyColor(disp, obj.signature());
        
        while ( obj.getLink(inx++, P, Q) )
            displayPoint(P, disp);
        inx = 0;

        while ( obj.getLink(inx++, P, Q) )
        {
            if (modulo) modulo->fold(Q, P);
            gleTube(P, Q, w, gleTube1B);
        }
    }
    /**
     This displays the Solid connecting two Aster as a spindle.
     Used for Cleo Kozlowski simulation of C. elegans
     */
    if ( disp->style & 1 && obj.tag() == Fake::TAG )
    {
        Solid const* so = static_cast<const Fake&>(obj).solid();
        if ( so && so->nbPoints() >= 4 )
        {
            bodyColor(so->prop->disp, so->signature());
#if ( DIM == 3 )
            glPushMatrix();
            Vector3 a = 0.5 * (so->posP(0) + so->posP(2));
            Vector3 b = 0.5 * (so->posP(1) + so->posP(3));
            const real diam = 1;
            Vector3 dir = b-a;
            Vector3 P1  = dir.orthogonal(diam);
            Vector3 P2  = cross(dir, P1).normalized(diam);
            gleTransRotate(P1, P2, dir, a);
            glColor3f(0.6,0.6,0.6);
            gleDualPass(gleBarrel1);
            glPopMatrix();
#else
            for ( unsigned ii = 0; ii < so->nbPoints(); ii+=2 )
                gleTube(so->posPoint(ii), so->posPoint(ii+1), w, gleHexTube1);
#endif
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void Display3::displaySinglesF(SingleSet const& set)
{
    for ( Single * obj=set.firstF(); obj ; obj=obj->next() )
    {
        obj->disp()->color2.load_both();
        displayPoint(obj->posFoot(), obj->disp());
    }
}


void Display3::displaySinglesA(SingleSet const& set)
{
    // display the Hands
    for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
    {
        if ( !obj->fiber()->disp->visible )
            continue;

        const PointDisp * disp = obj->disp();
        Vector ph = obj->posHand();
        
        disp->color.load_both();

        if ( obj->hasForce() && disp->width > 0 )
        {
            Vector pf = obj->posFoot();
            if (modulo)
                modulo->fold(pf, ph);
#if ( DIM >= 3 )
            gleTube(ph, pf, disp->width*sFactor);
#else
            gleBand(ph, disp->width*sFactor, disp->color, pf, disp->width*sFactor, disp->color.alpha_scaled(0.5));
#endif
        }
        displayPoint(ph, disp);
    }
}


//------------------------------------------------------------------------------
#pragma mark -
/**
 Display always Hand1 of Couple
 */
void Display3::displayCouplesF1(CoupleSet const& set)
{
    for ( Couple * cx = set.firstFF(); cx ; cx=cx->next() )
    {
        cx->disp1()->color2.load_both();
        displayPoint(cx->posFree(), cx->disp1());
    }
}


/**
 Display either Hand1 or Hand2, exposing both sides with equal chances.
 This gives the impression that Couple flicker randomly between frames,
 as if they were two-sided balls 'rotating' very fast.
 */
void Display3::displayCouplesF2(CoupleSet const& set)
{
    Couple * nxt;
    Couple * obj = set.firstFF();
    
    if ( set.sizeFF() % 2 )
    {
        nxt = obj->next();
        obj->disp1()->color2.load_both();
        displayPoint(obj->posFree(), obj->disp1());
        obj = nxt;
    }
    while ( obj )
    {
        nxt = obj->next();
        obj->disp2()->color2.load_both();
        displayPoint(obj->posFree(), obj->disp2());
        obj = nxt->next();
        nxt->disp1()->color2.load_both();
        displayPoint(nxt->posFree(), nxt->disp1());
    }
}


void Display3::displayCouplesA(CoupleSet const& set)
{
    for ( Couple * cx=set.firstAF(); cx ; cx=cx->next() )
    {
        if ( cx->fiber1()->disp->visible )
        {
#if ( 0 )
            // ENDOCYTOSIS 2015
            if ( cx->fiber1()->disp->color.transparent() )
                cx->disp1()->color.load_both(cx->fiber1()->disp->color.alpha());
            else
#endif
            cx->disp1()->color.load_both();
            displayPoint(cx->posHand1(), cx->disp1());
        }
    }
    
    for ( Couple * cx=set.firstFA(); cx ; cx=cx->next() )
    {
        if ( cx->fiber2()->disp->visible )
        {
#if ( 0 )
            // ENDOCYTOSIS 2015
            if ( cx->fiber2()->disp->color.transparent() )
                cx->disp1()->color.load_both(cx->fiber2()->disp->color.alpha());
            else
#endif
            cx->disp2()->color.load_both();
            displayPoint(cx->posHand2(), cx->disp2());
        }
    }
}


void Display3::displayCouplesB(CoupleSet const& set)
{
    for ( Couple * cx=set.firstAA(); cx ; cx=cx->next() )
    {
        // do not display Couple if the associated Fibers are both hidden
        if ( !cx->fiber1()->disp->visible && !cx->fiber2()->disp->visible )
            continue;
        // only display if bridging two anti-parallel filaments
        if ( prop->couple_select & 8  && cx->cosAngle() > 0 )
            continue;
        // only display if bridging two parallel filaments
        if ( prop->couple_select & 16 && cx->cosAngle() < 0 )
            continue;
        
        const PointDisp * pd1 = cx->disp1();
        const PointDisp * pd2 = cx->disp2();
        
        Vector p1 = cx->posHand1();
        Vector p2 = cx->posHand2();
        if (modulo)
            modulo->fold(p2, p1);
        
        Vector dir = p2 - p1;
        real dn = dir.norm();
        
#if ( 1 )
        if ( dn > 1e-5 )
        {
            // position the heads at the surface of the filaments:
            const real rad1 = cx->fiber1()->prop->disp->line_width * sFactor;
            const real rad2 = cx->fiber2()->prop->disp->line_width * sFactor;
            p1 += dir * ( rad1 / dn );
            p2 -= dir * ( rad2 / dn );
        }
#endif

        if ( pd1 == pd2 )
        {
#if ( 0 )
            // ENDOCYTOSIS 2015
            if ( cx->fiber1()->disp->color.transparent() )
            {
                pd1->color.load_both(cx->fiber1()->disp->color.alpha());
                glDepthMask(GL_FALSE);
                gleTube(p1, p2, pd2->width*sFactor, gleHexTube1);
                glDepthMask(GL_TRUE);
                continue;
            }
#endif
            pd1->color.load_both();
            gleTube(p1, p2, pd2->width*sFactor, gleHexTube1);
            displayPoint(p1, pd1);
            displayPoint(p2, pd2);
        }
        else if ( dn > 1e-5 )
        {
            Vector mid = 0.5 * ( p1 + p2 );
            gleTube(p1, pd1->width*sFactor, pd1->color, p2, pd2->width*sFactor, pd2->color);

            glEnable(GL_CLIP_PLANE5);
            if ( pd1->visible )
            {
                setClipPlane(GL_CLIP_PLANE5, -dir, mid);
                pd1->color.load_both();
                displayPoint(p1, pd1);
            }

            if ( pd2->visible )
            {
                setClipPlane(GL_CLIP_PLANE5,  dir, mid);
                pd2->color.load_both();
                displayPoint(p2, pd2);
            }
            glDisable(GL_CLIP_PLANE5);
        }
        else
        {
            pd1->color.load_both();
            displayPoint(p1, pd1);
            pd2->color.load_both();
            displayPoint(p2, pd2);
        }
    }
}

