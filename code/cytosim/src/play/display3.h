// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef DISPLAY3_H
#define DISPLAY3_H

#include "display.h"
#include "real.h"
#include "opengl.h"
#include "vector.h"
class PointDisp;

///Cytosim display class for style=3
/**
 This style is for rendering in 3D.
 It uses Lighting for better volume rendering
 */
class Display3 : public Display
{
private:
    
    ///display a spherical object, finer than a point
    void displayBall(Vector const&, real radius);

    ///display a point with a small sphere
    void displayPoint(Vector const&, PointDisp const*);

public:
        
    ///constructor
    Display3(DisplayProp const*);
    
    ///destructor
    ~Display3() {}
    
    ///display the given simulation state using OpenGL commands
    void displayObjects(Simul const&);
    
    ///display Fiber PLUS_END
    void displayFiberEndMinus(Fiber const&) const;
    
    ///display Fiber MINUS_END
    void displayFiberEndPlus(Fiber const&) const;
    
    ///display Fiber linear features
    void displayFiberLines(Fiber const&) const;
    
    ///display Fiber linear features located between [MINUS_END, abs]
    void displayFiberLinesMinus(Fiber const&, real abs, real width) const;
    
    ///display Fiber linear features located between [abs, PLUS_END]
    void displayFiberLinesPlus(Fiber const&, real abs, real width) const;

    ///display Fiber Lattice
    void displayFiberLattice(Fiber const&) const;
    
    ///display Fiber point-like features
    void displayFiberPoints(Fiber const&) const;
    
    ///display Fiber speckles
    void displayFiberSpeckles(Fiber const&) const;
    
    ///display the Solids
    void displaySolid(Solid const&);
 
    ///display the transparent parts of Solid
    void displaySolidT(Solid const&, unsigned int);

    ///display a Bead
    void displayBead(Bead const&);
    
    ///display transparent membrane of Bead
    void displayBeadT(Bead const&);
    
    ///display a Sphere
    void displaySphere(Sphere const&);
    
    ///display transparent membrane of Sphere
    void displaySphereT(Sphere const&);
    
    ///display the free Single
    void displaySinglesF(SingleSet const&);

    ///display the attached Single
    void displaySinglesA(SingleSet const&);

    ///display free Couple
    void displayCouplesF1(CoupleSet const&);

    ///display free Couple, randomizing which Hand is drawn
    void displayCouplesF2(CoupleSet const&);

    ///display attached Couple
    void displayCouplesA(CoupleSet const&);

    ///display bridging Couple
    void displayCouplesB(CoupleSet const&);
 
    //display an Organizer
    void displayOrganizer(Organizer const&);
};

#endif

