// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "sim.h"
#include "simul.h"
#include "display1.h"
#include "modulo.h"

#include "fake.h"
#include "line_disp.h"
#include "point_disp.h"

#include "opengl.h"
#include "gle.h"
#include "gle_color_list.h"
#include "glapp.h"
#include "glut.h"

using namespace gle;
extern Modulo const* modulo;

//------------------------------------------------------------------------------

Display1::Display1(DisplayProp const* dp) : Display(dp)
{
}


void Display1::displayObjects(Simul const& sim)
{
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    displayFields(sim.fields);
    
    if ( prop->couple_select & 1 )
        displayCouplesF(sim.couples);
    
    if ( prop->couple_select & 2 )
        displayCouplesA(sim.couples);
    
    if ( prop->single_select & 1 )
        displaySinglesF(sim.singles);
    
    displayFibers(sim.fibers);
    
#if ( DIM == 3 )
    
    glEnable(GL_LIGHTING);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    View& view = glApp::currentView();
    view.inner_color.load_back();

#endif
    
    displayBeads(sim.beads);
    displaySolids(sim.solids);
    displaySpheres(sim.spheres);
    
#if ( DIM == 3 )
    
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);

#endif

    if ( prop->couple_select & 4 )
        displayCouplesB(sim.couples);
    
    if ( prop->single_select & 2 )
        displaySinglesA(sim.singles);
    
#if ( DIM == 3 )
    
    glEnable(GL_LIGHTING);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

#endif

    displayOrganizers(sim.organizers);
    displayMisc(sim);
}

//------------------------------------------------------------------------------
#pragma mark -

void Display1::displayBall(Vector const& pos, real radius)
{
    glPushMatrix();
    gleTranslate(pos);
    gleScale(radius);
    if ( DIM == 3 )
    {
        glCullFace(GL_FRONT);
        gleSphere2B();
        glCullFace(GL_BACK);
        gleSphere4B();
    }
    else
        gleCircleSB();
    glPopMatrix();
}


/// draw a little sphere
inline void Display1::displayPoint(Vector const& pos, PointDisp const* disp)
{
    glPushMatrix();
    gleTranslate(pos);
    gleScale(disp->size*sFactor);
    gleSphere1B();
    glPopMatrix();
}


//------------------------------------------------------------------------------
#pragma mark -

void Display1::displayBead(Bead const& obj)
{
    const PointDisp * disp = obj.prop->disp;

    // display center:
    if ( disp->style & 2  && disp->perceptible )
    {
        bodyColor(disp, obj.signature());
        displayPoint(obj.position(), disp);
    }
    
#if ( DIM == 2 )
    // display outline:
    if ( disp->style & 4 )
    {
        lineWidth(disp->width);
        bodyColor(disp, obj.signature());
        gleObject(obj.position(), obj.radius(), gleCircleLB);
    }
#endif
}


/**
 Display a semi-transparent disc / sphere
 */
void Display1::displayBeadT(Bead const& obj)
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

void Display1::displaySolid(Solid const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    //display points of the Solids
    if ( disp->style & 2  &&  disp->size > 0  && disp->perceptible )
    {
        bodyColor(disp, obj.signature());
        for ( unsigned ii = 0; ii < obj.nbPoints(); ++ii )
            displayPoint(obj.posP(ii), disp);
    }
    
    //display outline of spheres
    if ( disp->style & 4 )
    {
        lineWidth(disp->width);
        bodyColor(disp, obj.signature());
#if ( DIM == 2 )
        for ( unsigned ii = 0; ii < obj.nbPoints(); ++ii )
        {
            if ( obj.radius(ii) > 0 )
                gleObject(obj.posP(ii), obj.radius(ii), gleCircleLB);
        }
#elif ( DIM == 3 )
        //special display for ParM simulations (DYCHE)
        if ( obj.mark()  &&  obj.nbPoints() >= 3 )
            gleObject(obj.posP(0), obj.diffPoints(1, 0), obj.radius(0), gleCircleLB);
#endif
    }
    
    //print the number for each solid
    if ( disp->style & 8 )
    {
        char tmp[8];
        bodyColor(disp, obj.signature());
        snprintf(tmp, sizeof(tmp), "%lu", obj.identity());
        gleDrawText(obj.posP(0), tmp, GLUT_BITMAP_HELVETICA_10);
    }
    
    //draw polygon around model points of Solid
    if ( disp->style & 16 )
    {
        lineWidth(disp->width);
        bodyColor(disp, obj.signature());
        glBegin(GL_LINE_LOOP);
        for ( unsigned ii = 0; ii < obj.nbPoints(); ++ii )
            gleVertex(obj.posPoint(ii));
        glEnd();
    }
}

/**
 Display a semi-transparent disc / sphere
 */
void Display1::displaySolidT(Solid const& obj, unsigned int ii)
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

void Display1::displaySphere(Sphere const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    //display center and surface points
    if ( disp->style & 2  &&  disp->perceptible )
    {
        bodyColor(disp, obj.signature());
        displayPoint(obj.posP(0), disp);
        for ( unsigned ii = obj.nbRefPoints(); ii < obj.nbPoints(); ii++ )
            displayPoint(obj.posP(ii), disp);
    }
    
    //display reference points
    if ( disp->style & 8  &&  disp->perceptible  )
    {
        bodyColor(disp, obj.signature());
        for ( unsigned ii = 1; ii < obj.nbRefPoints(); ii++ )
            displayPoint(obj.posP(ii), disp);
    }
}


void Display1::displaySphereT(Sphere const& obj)
{
    const PointDisp * disp = obj.prop->disp;

    if ( disp->style & 5 )
    {
        bodyColorT(disp, obj.signature());
        lineWidth(disp->width);
        
#if ( DIM < 3 )
        
        if ( disp->style & 1 )
            gleObject(obj.posP(0), obj.radius(), gleCircleSB);
        else
            gleObject(obj.posP(0), obj.radius(), gleCircleLB);
        
#else
        
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        /* Note: The rotation matrix for the sphere calculated below from the
         reference points, includes scaling by the radius of the sphere.
         We then use a primitive for a sphere of radius 1.
         */
        const Vector C = obj.posP(0);
        gleTransRotate(obj.posP(1)-C, obj.posP(2)-C, obj.posP(3)-C, C);
        
        //draw transparent envelope
        if ( disp->style & 1 )
            gleDualPass(gleSphere8B);
        if ( disp->style & 4 )
        {
            disp->color2.load_front();
            gleThreeBands(64);
        }
        glPopMatrix();
        
#endif
    }
}

//------------------------------------------------------------------------------
void Display1::displayOrganizer(Organizer const& obj)
{
    const PointDisp * disp = obj.disp();
    
    if ( !disp )
        return;

    if ( disp->style & 2 )
    {
        Vector P, Q;
        unsigned inx = 0;
        glDisable(GL_LIGHTING);
        
        while ( obj.getLink(inx++, P, Q) )
            disp->drawA(P);
        inx = 0;

        lineWidth(disp->width);
        bodyColor(disp, obj.signature());
        glBegin(GL_LINES);
        while ( obj.getLink(inx++, P, Q) )
        {
            if (modulo) modulo->fold(Q, P);
            gleVertex(P);
            gleVertex(Q);
        }
        glEnd();
    }

    /**
     This display the Solid connecting two Aster as a spindle.
     Used for Cleo Kozlowski simulation of C. elegans
     */
    if ( disp->style & 1 && obj.tag() == Fake::TAG )
    {
        Solid const* so = static_cast<const Fake&>(obj).solid();
        if ( so && so->nbPoints() >= 4 )
        {
            bodyColor(disp, obj.signature());
#if ( DIM == 3 )
            glEnable(GL_LIGHTING);
            glPushMatrix();
            Vector3 a = 0.5*(so->posP(0) + so->posP(2));
            Vector3 b = 0.5*(so->posP(1) + so->posP(3));
            const real diam = 1;
            Vector3 dir = b-a;
            Vector3 P1  = dir.orthogonal(diam);
            Vector3 P2  = cross(dir, P1).normalized(diam);
            gleTransRotate(P1, P2, dir, a);
            gleDualPass(gleBarrel1);
            glPopMatrix();
            glDisable(GL_LIGHTING);
#else
            glBegin(GL_LINES);
            for ( unsigned ii = 0; ii < so->nbPoints(); ++ii )
                gleVertex(so->posPoint(ii));
            glEnd();
#endif
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 display the attached position of free singles
 */
void Display1::displaySinglesF(const SingleSet & set)
{
    for ( Single * obj=set.firstF(); obj ; obj=obj->next() )
        obj->disp()->drawF(obj->posFoot());
}


void Display1::displaySinglesA(const SingleSet & set)
{
    // display the Hands
    for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
    {
        const PointDisp * disp = obj->disp();
        if ( disp->perceptible  &&  obj->fiber()->disp->visible )
        {
            Vector ph = obj->posHand();
            
            disp->drawA(ph);
            
            if ( obj->hasForce() && disp->width > 0 )
            {
                Vector ps = obj->posSide();
                Vector pf = obj->posFoot();
                if ( modulo )
                {
                    modulo->fold(pf, ph);
                    modulo->fold(ps, ph);
                }
                
#if ( DIM >= 3 )
                gleTube(ph, disp->width*sFactor, disp->color, pf, disp->width*sFactor, disp->color.alpha_scaled(0.5));
#else
                disp->color.load();
                gleBand(ph, disp->width*sFactor, ps, disp->width*sFactor);
                gleBand(ps, disp->width*sFactor, disp->color, pf, disp->width*sFactor, disp->color.alpha_scaled(0.5));
#endif
            }
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Always display Hand1 of the Couple.
 */
void Display1::displayCouplesF1(CoupleSet const& set)
{
    for ( Couple * cx = set.firstFF() ; cx ; cx=cx->next() )
    {
        if ( cx->active() )
            cx->disp1()->drawF(cx->posFree());
        else
            cx->disp1()->drawI(cx->posFree());
    }
}

/**
 We display either Hand1 or Hand2, exposing both sides with equal chances.
 This gives the impression that Couple flicker randomly between frames,
 as if they were two-sided balls 'rotating' very fast.
 */
void Display1::displayCouplesF2(CoupleSet const& set)
{
    Couple * nxt;
    Couple * obj = set.firstFF();
    // this loop is unrolled, processing objects 2 by 2:
    if ( set.sizeFF() % 2 )
    {
        nxt = obj->next();
        if ( obj->active() )
            obj->disp1()->drawF(obj->posFree());
        else
            obj->disp1()->drawI(obj->posFree());
        obj = nxt;
    }
    while ( obj )
    {
        nxt = obj->next();
        if ( obj->active() )
            obj->disp2()->drawF(obj->posFree());
        else
            obj->disp2()->drawI(obj->posFree());
        obj = nxt->next();
        if ( nxt->active() )
            nxt->disp1()->drawF(nxt->posFree());
        else
            nxt->disp1()->drawI(nxt->posFree());
    }
}


void Display1::displayCouplesA(CoupleSet const& set)
{
    // display bound couples
    for ( Couple * cx=set.firstAF(); cx ; cx=cx->next() )
    {
        if ( cx->fiber1()->disp->visible )
        {
            if ( cx->active() )
                cx->disp1()->drawF(cx->posHand1());
            else
                cx->disp1()->drawI(cx->posHand1());
        }
    }
    
    for ( Couple * cx=set.firstFA(); cx ; cx=cx->next() )
    {
        if ( cx->fiber2()->disp->visible )
        {
            if ( cx->active() )
                cx->disp2()->drawF(cx->posHand2());
            else
                cx->disp2()->drawI(cx->posHand2());
            
        }
    }
}


void Display1::displayCouplesB(CoupleSet const& set)
{
    for ( Couple * cx=set.firstAA(); cx ; cx=cx->next() )
    {
        // do not display Couple if the associated Fibers are both hidden
        if ( cx->fiber1()->disp->visible==0  &&  cx->fiber2()->disp->visible==0 )
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
        
        if ( pd1 == pd2 )
        {
            if ( pd1->perceptible )
            {
                pd1->color.load();
#if ( DIM == 2 )
                //gleBand(p1, pd1->width*sFactor, p2, pd2->width*sFactor);
                real rad = pd1->width*sFactor;
                gleMan(p1, rad * cx->dirFiber1(), p2, rad * cx->dirFiber2());
#else
                lineWidth(pd1->width);
                glBegin(GL_LINES);
                gleVertex(p1);
                gleVertex(p2);
                glEnd();
#endif
            }
        }
        else
        {
            if ( pd1->perceptible || pd2->perceptible )
            {
#if ( DIM == 2 )
                //gleBand(p1, pd1->width*sFactor, pd1->color, p2, pd2->width*sFactor, pd2->color);
                gleMan(p1, ( pd1->width * sFactor ) * cx->dirFiber1(), pd1->color,
                       p2, ( pd2->width * sFactor ) * cx->dirFiber2(), pd2->color);
#else
                lineWidth(pd1->width);
                glBegin(GL_LINES);
                gleVertex(p1);
                gleVertex(p2);
                glEnd();
#endif
            }
        }
        
        if ( cx->active() )
        {
            pd1->drawA(p1);
            pd2->drawA(p2);
        }
        else
        {
            pd1->drawI(p1);
            pd2->drawI(p2);
        }
    }
}

