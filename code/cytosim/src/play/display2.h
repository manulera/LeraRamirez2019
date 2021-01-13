// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef DISPLAY2_H
#define DISPLAY2_H

#include "display.h"
class PointDisp;

///Cytosim display class for style=2
/**
 This style produces a fast 2D display.
 Some of the parameters in PointDisp are ignored.

 Point-like objects are rendered using OpenGL::Points and display-lists.
 All points are therefore displayed with the same size.
 */
class Display2 : public Display
{
    ///display a ball
    void displayBall(Vector const&, real radius);
    
    ///display a point
    void displayPoint(Vector const&, PointDisp const*);
    
public:
    
    ///constructor
    Display2(DisplayProp const*);
    
    ///destructor
    ~Display2() {}
    
     
    ///display the given simulation state using OpenGL commands
    void displayObjects(Simul const&);
   
    ///display Fibers with offset
    void displayFiber(Fiber const&);
    
    ///display the Solids
    void displaySolid(Solid const&);
 
    ///display the transparent part for the Solids
    void displaySolidT(Solid const&, unsigned int);
    
    ///display a Bead
    void displayBead(Bead const&);
    
    ///display transparent membrane of Bead
    void displayBeadT(Bead const&);
    
    ///display a Sphere
    void displaySphere(Sphere const&);
    
    ///display transparent membrane of Sphere
    void displaySphereT(Sphere const&);
    
    ///display the free Singles
    void displaySinglesF(SingleSet const&);

    ///display the attached Singles
    void displaySinglesA(SingleSet const&);
#ifdef TRAP_SINGLES
    ///display the trapped Singles
    void displaySinglesTrapped(SingleSet const&);

#endif
    ///display the free Couples
    void displayCouplesF1(CoupleSet const&);

    ///display the free Couples, randomizing which Hand is drawn
    void displayCouplesF2(CoupleSet const&);

    ///display the attached Couples
    void displayCouplesA(CoupleSet const&);

    ///display the bridging Couples
    void displayCouplesB(CoupleSet const&);
    
    ///display an Organizer
    void displayOrganizer(Organizer const&);
    
    ///display the tension of a couple;
    void displayCoupleTension(const Vector & Pa, Fiber * fba, PointDisp * dispa, const Vector&  Pb, Fiber * fbb, PointDisp * dispb , real force);
};

#endif

