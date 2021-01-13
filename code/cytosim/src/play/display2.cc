// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "sim.h"
#include "simul.h"
#include "display2.h"
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



#define EXPLODE_DISPLAY


Display2::Display2(DisplayProp const* dp) : Display(dp)
{
}


void Display2::displayObjects(Simul const& sim)
{
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    displayFields(sim.fields);
    
    if ( prop->couple_select & 1 )
        displayCouplesF(sim.couples);

    if ( prop->single_select & 1 )
        displaySinglesF(sim.singles);
    
    displayFibers(sim.fibers);
    
    if ( prop->couple_select & 2 )
        displayCouplesA(sim.couples);
    
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
    {
        displaySinglesA(sim.singles);
#ifdef TRAP_SINGLES

        displaySinglesTrapped(sim.singles);

#endif
    }
    
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

void Display2::displayBall(Vector const& pos, real radius)
{
    glPushMatrix();
    gleTranslate(pos);
    gleScale(radius);
    if ( DIM == 3 )
    {
        assert_true(glIsEnabled(GL_CULL_FACE));
        glCullFace(GL_FRONT);
        gleSphere2B();
        glCullFace(GL_BACK);
        gleSphere4B();
    }
    else
        gleCircleSB();
    glPopMatrix();
}


/// this version usually draws a little sphere
inline void Display2::displayPoint(Vector const& pos, PointDisp const* disp)
{
    if ( disp->perceptible )
    {
#if ( 0 )
        // draw a OpenGL point
        pointSize(disp->size);
        glBegin(GL_POINTS);
        gleVertex(pos);
        glEnd();
#else
        /// draw a little sphere
        glPushMatrix();
        gleTranslate(pos);
        gleScale(disp->size*sFactor);
        gleSphere1B();
        glPopMatrix();
#endif
    }
}

//------------------------------------------------------------------------------
#pragma mark -


void Display2::displayFiber(Fiber const& fib)
{
#ifdef EXPLODE_DISPLAY
    //translate whole display to display the Fiber
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    gleTranslate(0, fib.disp->explode_shift, 0);
#endif
    Display::displayFiber(fib);
#ifdef EXPLODE_DISPLAY
    glPopMatrix();
#endif
}

//------------------------------------------------------------------------------
#pragma mark -

void Display2::displayBead(Bead const& obj)
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
 Display a semi-transparent disc / sphere
 */
void Display2::displayBeadT(Bead const& obj)
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

void Display2::displaySolid(Solid const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    //display points:
    if ( disp->style & 2  &&  disp->size > 0 )
    {
        bodyColor(disp, obj.signature());
        for ( unsigned ii = 0; ii < obj.nbPoints(); ++ii )
            displayPoint(obj.posP(ii), disp);
    }
    
    //display outline of spheres
    if ( disp->style & 4 )
    {
        bodyColor(disp, obj.signature());
        lineWidth(disp->width);
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
    
    //print the number for each Solid
    if ( disp->style & 8 )
    {
        char tmp[8];
        snprintf(tmp, sizeof(tmp), "%lu", obj.identity());
        bodyColor(disp, obj.signature());
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
void Display2::displaySolidT(Solid const& obj, unsigned int ii)
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

void Display2::displaySphere(Sphere const& obj)
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
    if ( disp->style & 8  &&  disp->perceptible )
    {
        bodyColor(disp, obj.signature());
        for ( unsigned ii = 1; ii < obj.nbRefPoints(); ii++ )
            displayPoint(obj.posP(ii), disp);
    }
}


void Display2::displaySphereT(Sphere const& obj)
{
    const PointDisp * disp = obj.prop->disp;

    if ( disp->style & 5 )
    {
        bodyColorT(disp, obj.signature());
        lineWidth(disp->width);
        
#if (DIM <= 2)
        
        if ( disp->style & 1 )
            gleObject(obj.posP(0), obj.radius(), gleCircleSB);
        else
            gleObject(obj.posP(0), obj.radius(), gleCircleLB);
        
#elif (DIM == 3)
        
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
void Display2::displayOrganizer(Organizer const& obj)
{
    const PointDisp * disp = obj.disp();
    
    if ( !disp )
        return;

    if ( disp->style & 2 )
    {
        Vector P, Q;
        unsigned inx = 0;
        glDisable(GL_LIGHTING);

        bodyColor(disp, obj.signature());
        pointSize(disp->size);
        glBegin(GL_POINTS);
        while ( obj.getLink(inx++, P, Q) )
            gleVertex(P);
        inx = 0;
        glEnd();

        bodyColor(disp, obj.signature());
        lineWidth(disp->width);
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

#ifdef EXPLODE_DISPLAY

inline void shiftedVertex(Vector const& pos, real shift)
{
#if ( DIM == 3 )
#ifdef REAL_IS_FLOAT
    glVertex3f(pos.XX, pos.YY+shift, pos.ZZ);
#else
    glVertex3d(pos.XX, pos.YY+shift, pos.ZZ);
#endif
#elif ( DIM == 2 )
#ifdef REAL_IS_FLOAT
    glVertex2f(pos.XX, pos.YY+shift);
#else
    glVertex2d(pos.XX, pos.YY+shift);
#endif
#else
#ifdef REAL_IS_FLOAT
    glVertex2f(pos.XX, shift);
#else
    glVertex2d(pos.XX, shift);
#endif
#endif
}

inline void drawVertex(Vector const& pos, const Fiber * fib, const PointDisp* disp)
{
    if ( disp->perceptible && fib->disp->visible )
    {
        disp->color.load();
        shiftedVertex(pos, fib->disp->explode_shift);
    }
}


inline void drawVertex2(Vector const& pos, const Fiber * fib, const PointDisp* disp)
{
    if ( disp->perceptible && fib->disp->visible )
    {
        disp->color2.load();
        shiftedVertex(pos, fib->disp->explode_shift);
    }
}


inline void drawLink(const Vector & a, const Fiber * fib, const PointDisp* disp, const Vector & b)
{
    
    if ( disp->perceptible && fib->disp->visible )
    {
        disp->color.load();
        shiftedVertex(a, fib->disp->explode_shift);
        disp->color2.load();
        shiftedVertex(b, fib->disp->explode_shift);
    }
}


/**
 Draw two segments if `explode' is enabled
 */
inline void drawLink(const Vector & a, const Fiber * fibA, const PointDisp* dispA,
                     const Vector & b, const Fiber * fibB, const PointDisp* dispB)
{
    
    if ( dispA->perceptible && fibA->disp->visible )
    {
        dispA->color.load();
        shiftedVertex(a, fibA->disp->explode_shift);
        dispB->color.load();
        shiftedVertex(b, fibB->disp->explode_shift);
    }
    if ( dispB->perceptible && fibB->disp->visible && fibB->prop->disp->explode )
    {
        dispA->color.load();
        shiftedVertex(a, fibA->disp->explode_shift);
        dispB->color.load();
        shiftedVertex(b, fibB->disp->explode_shift);
    }
}

inline void drawLinkTrap(const Vector & a, const Fiber * fib, const PointDisp* disp, const Vector3 & b)
{
    
    if ( disp->perceptible && fib->disp->visible )
    {
        disp->color.load();
        shiftedVertex(a, fib->disp->explode_shift);
        disp->color2.load();
#if (DIM==1)
        shiftedVertex(Vector(b), b.YY);
#else
        shiftedVertex(Vector(b),0);
#endif
    }
}

#else

// define macros without spatial shift:

inline void drawVertex(Vector const& pos, const Fiber * fib, const PointDisp* disp)
{
    if ( disp->perceptible && fib->disp->visible )
    {
        disp->color.load();
        gleVertex(pos);
    }
}


inline void drawVertex2(Vector const& pos, const Fiber * fib, const PointDisp* disp)
{
    if ( disp->perceptible && fib->disp->visible )
    {
        disp->color2.load();
        gleVertex(pos);
    }
}


inline void drawLink(const Vector & a, const Fiber * fib, const PointDisp* disp, const Vector & b)
{
    if ( disp->perceptible && fib->disp->visible )
    {
        disp->color.load();
        gleVertex(a);
        disp->color2.load();
        gleVertex(b);
    }
}

inline void drawLink(const Vector & a, const Fiber * fibA, const PointDisp* dispA,
                     const Vector & b, const Fiber * fibB, const PointDisp* dispB)
{
    if (   dispA->perceptible && fibA->disp->visible
        && dispB->perceptible && fibB->disp->visible )
    {
        dispA->color.load();
        gleVertex(a);
        dispB->color.load();
        gleVertex(b);
    }
}

#endif


/// call glVertex using PointDisp::color
inline void drawVertex(const Vector & pos, const PointDisp* disp)
{
    if ( disp->perceptible )
    {
        disp->color.load();
        gleVertex(pos);
    }
}

/// call glVertex using PointDisp::color2
inline void drawVertex2(const Vector & pos, const PointDisp* disp)
{
    if ( disp->perceptible )
    {
        disp->color2.load();
        gleVertex(pos);
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void Display2::displaySinglesF(const SingleSet & set)
{
    //display the attached position of free singles:
    if ( prop->point_size > 0 )
    {
        pointSize(prop->point_size);
        glBegin(GL_POINTS);
        for ( Single * obj=set.firstF(); obj ; obj=obj->next() )
        {
#if ( DIM == 1 ) && defined EXPLODE_DISPLAY
            obj->disp()->color2.load();
            gleVertex(obj->posFoot().XX, ( obj->signature() * 0x1p-31 - 1 ) * 10);
#else
            drawVertex2(obj->posFoot(), obj->disp());
#endif
        }
        glEnd();
    }
}


void Display2::displaySinglesA(const SingleSet & set)
{
    // display the positions
    if ( prop->point_size > 0 )
    {
        pointSize(prop->point_size);
        glBegin(GL_POINTS);
        for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
            drawVertex(obj->posHand(), obj->fiber(), obj->disp());
        glEnd();
    }
    
    // display the links
    if ( prop->link_size > 0 )
    {
        lineWidth(prop->link_size);
        glBegin(GL_LINES);
        for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
            if ( obj->hasForce() )
            {
                Vector ph = obj->posHand();
                Vector pf = obj->posFoot();
                if (modulo) modulo->fold(pf, ph);
                drawLink(ph, obj->fiber(), obj->disp(), pf);
            }
        glEnd();
    }
}
#ifdef TRAP_SINGLES

void Display2::displaySinglesTrapped(const SingleSet & set)
{
    // display the positions
    if ( prop->point_size > 0 )
    {
        pointSize(prop->point_size);
        glBegin(GL_POINTS);
        for ( Single * obj=set.firstTrapped(); obj ; obj=obj->next() )
        {
            
            
            if (obj->attached())
                drawVertex2(obj->posHand(), obj->fiber(), obj->disp());
            else
            {

                Couple * cx = static_cast<Couple *>(obj->trapped_haMon);
                if (cx->linking())
                {
                    Vector3 center = cx->trapCenter(cx->fiber1()->disp->explode_shift,
                                                cx->fiber2()->disp->explode_shift,
                                                0);
                    obj->disp()->color2.load();
#if ( DIM != 3 )
                    gleVertex(center.XX, center.YY);
#else
                    gleVertex(center.XX, center.YY, center.ZZ);
#endif
                }
            }
        }
        glEnd();
    }
#if (TRAP_SINGLES ==1)
    // display the links
    if ( prop->link_size > 0 )
    {
        lineWidth(prop->link_size);
        glBegin(GL_LINES);
        for ( Single * obj=set.firstTrapped(); obj ; obj=obj->next() )
        {
            Vector ph = obj->posHand();
            Vector pf = obj->trappedOtherHandPos();
            if (modulo) modulo->fold(pf, ph);
            gle_color c = obj->disp()->color;
            obj->disp()->color = obj->disp()->color2;
            Hand * other_h = obj->trappedOtherHand();
            drawLink(ph, obj->fiber(), obj->disp(), pf,other_h->fiber(),obj->disp());
            obj->disp()->color = c;
        }
        glEnd();
    }
#endif
}
#endif
//------------------------------------------------------------------------------
#pragma mark -
/**
Always display Hand1 of Couple
 */
void Display2::displayCouplesF1(CoupleSet const& set)
{
    if ( prop->point_size > 0 )
    {
        pointSize(prop->point_size);
        
        glBegin(GL_POINTS);
        for ( Couple * cx = set.firstFF(); cx ; cx=cx->next() )
        {
#if ( DIM == 1 ) && defined EXPLODE_DISPLAY
            PointDisp * disp = cx->disp1();
            
            if ( disp->perceptible )
            {
                disp->color2.load();
                gleVertex(cx->posFree().XX, ( cx->signature() * 0x1p-31 - 1 ) * 10);
            }
#else
            drawVertex2(cx->posFree(), cx->disp1());
#endif
        }
        glEnd();
        
#if ( DIM > 1 )
        // display inactive Couples with bitmap:
        for ( Couple * cx = set.firstFF(); cx ; cx=cx->next() )
            if ( !cx->active() && cx->disp1()->perceptible )
                cx->disp1()->drawI(cx->posFree());
#endif
    }
}


/**
 Display either Hand1 or Hand2, exposing both sides with equal chances.
 This gives the impression that Couple flicker randomly between frames,
 as if they were two-sided balls 'rotating' very fast.
 */
void Display2::displayCouplesF2(CoupleSet const& set)
{
    if ( prop->point_size > 0 )
    {
        Couple * nxt;
        Couple * obj = set.firstFF();

        pointSize(prop->point_size);
        glBegin(GL_POINTS);
        if ( set.sizeFF() % 2 )
        {
            nxt = obj->next();
            drawVertex2(obj->posFree(), obj->disp1());
            obj = nxt;
        }
        while ( obj )
        {
            nxt = obj->next();
            drawVertex2(obj->posFree(), obj->disp2());
            obj = nxt->next();
            drawVertex2(nxt->posFree(), nxt->disp1());
        }
        glEnd();
    }
}


void Display2::displayCouplesA(CoupleSet const& set)
{
    if ( prop->point_size > 0 )
    {
        // display bound couples
        pointSize(prop->point_size);
        glBegin(GL_POINTS);

        for ( Couple * cx=set.firstAF(); cx ; cx=cx->next() )
            drawVertex2(cx->posHand1(), cx->fiber1(), cx->disp1());
        
        for ( Couple * cx=set.firstFA(); cx ; cx=cx->next() )
            drawVertex2(cx->posHand2(), cx->fiber2(), cx->disp2());
        glEnd();
#ifdef TRAP_SINGLES
        glBegin(GL_LINES);
        for ( Couple * cx=set.firstAF(); cx ; cx=cx->next() )
        {
            
            if (cx->trapped())
            
            {
                Hand * h_trap = 0, * dummy = 0;
                cx->trapped_haMon->getHands(h_trap, dummy);
                if(h_trap->attached())
                    drawLink(cx->posHand1(), cx->fiber1(), cx->disp1(), h_trap->pos(), h_trap->fiber(), h_trap->prop->disp);
            }
        }
        
        for ( Couple * cx=set.firstFA(); cx ; cx=cx->next() )
        {
            
            if (cx->trapped()) {
                Hand * h_trap = 0, * dummy = 0;
                cx->trapped_haMon->getHands(h_trap, dummy);
                if(h_trap->attached())
                    drawLink(cx->posHand2(), cx->fiber2(), cx->disp2(), h_trap->pos(), h_trap->fiber(), h_trap->prop->disp);
            }
        }
        glEnd();
#endif

    }
}


void Display2::displayCouplesB(CoupleSet const& set)
{
    // display bridging couples
    if ( prop->point_size > 0 )
    {
        pointSize(prop->point_size);
        glBegin(GL_POINTS);
        for ( Couple * cx=set.firstAA(); cx ; cx=cx->next() )
        {
            // only display if bridging two anti-parallel filaments
            if ( prop->couple_select & 8  && cx->cosAngle() > 0 )
                continue;
            // only display if bridging two parallel filaments
            if ( prop->couple_select & 16 && cx->cosAngle() < 0 )
                continue;
            
            drawVertex(cx->posHand1(), cx->fiber1(), cx->disp1());
            drawVertex(cx->posHand2(), cx->fiber2(), cx->disp2());
        }
        glEnd();
    }
    
    // display the link for bridging couples
    if ( prop->link_size > 0 )
    {
        lineWidth(prop->link_size);
        glBegin(GL_LINES);
        for ( Couple * cx=set.firstAA(); cx ; cx=cx->next() )
        {
            // only display if bridging two anti-parallel filaments
            if ( prop->couple_select & 8  && cx->cosAngle() > 0 )
                continue;
            // only display if bridging two parallel filaments
            if ( prop->couple_select & 16 && cx->cosAngle() < 0 )
                continue;
            
            Vector P = cx->posHand1();
            Vector Q = cx->posHand2();
            
#ifndef TRAP_SINGLES
            if (modulo) modulo->fold(Q, P);
            drawLink(P, cx->fiber1(), cx->disp1(), Q, cx->fiber2(), cx->disp2());
#else
            // A bit ugly, but avoids having even more virtual functions in couple
            
            bool normal_link = true;
            Hand * h_trap = 0, * dummy = 0;
            if (cx->trapped())
            {
                cx->trapped_haMon->getHands(h_trap, dummy);
                normal_link = !h_trap->attached();
            }
            
            if (normal_link)
            {
                if (modulo) modulo->fold(Q, P);
                drawLink(P, cx->fiber1(), cx->disp1(), Q, cx->fiber2(), cx->disp2());
            }
            else
            {
                //@todo Here there should be a propper implementation when boundary conditions are on   
                Vector R = h_trap->pos();
                Vector3 center = cx->trapCenter(cx->fiber1()->disp->explode_shift,
                                                cx->fiber2()->disp->explode_shift,
                                                h_trap->fiber()->disp->explode_shift);
                drawLinkTrap(P, cx->fiber1(), cx->disp1(), center);
                drawLinkTrap(Q, cx->fiber2(), cx->disp2(), center);
                drawLinkTrap(R, h_trap->fiber(), h_trap->prop->disp, center);
            }
            
#endif
//            displayCoupleTension(P, cx->fiber1(), cx->disp1(), Q, cx->fiber2(), cx->disp2(), cx->force().norm());
        }
        glEnd();
    }
}

void Display2::displayCoupleTension(const Vector & Pa, Fiber * fba, PointDisp * dispa, const Vector&  Pb, Fiber * fbb, PointDisp * dispb , real force)
{
    gle_color old_ca = dispa->color;
    real old_sa = dispa->size;
    gle_color old_cb = dispb->color;
    real old_sb = dispb->size;
    
    // Display tension in the link
    gle_color col1 =jet_color(force, 0.3);
    dispa->color = col1;
    dispb->color = col1;
    // Make point size zero, so we dont get the jet dots on top of the other dots (pre_set in color)
    dispa->size = 0.;
    dispb->size = 0.;
    
    drawLink(Pa, fba, dispa, Pb, fbb, dispb);
    
    dispa->color = old_ca;
    dispb->color = old_cb;
    dispa->size = old_sa;
    dispb->size = old_sb;
}
