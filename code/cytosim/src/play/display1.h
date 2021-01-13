// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef DISPLAY1_H
#define DISPLAY1_H

#include "display.h"
class PointDisp;

///Cytosim display class for style=1
/**
 This is the standard 2D display.
 It implements most of the characteristics in PointDisp and FiberDisp
 */
class Display1 : public Display
{
    ///display a ball
    void displayBall(Vector const&, real radius);
    
    ///display a point
    void displayPoint(Vector const&, PointDisp const*);

public:
    
    ///constructor
    Display1(DisplayProp const*);
    
    ///destructor
    ~Display1() {}
    
    
    ///display the given simulation state using OpenGL commands
    void displayObjects(Simul const&);
    
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

    ///display the free Couples
    void displayCouplesF1(CoupleSet const&);
    
    ///display free Couple, randomizing which Hand is drawn
    void displayCouplesF2(CoupleSet const&);
    
    ///display the attached Couples
    void displayCouplesA(CoupleSet const&);
    
    ///display the bridging Couples
    void displayCouplesB(CoupleSet const&);
    
    ///display an Organizer
    void displayOrganizer(Organizer const&);
};

#endif

