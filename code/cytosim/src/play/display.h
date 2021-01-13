// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef DISPLAY_H
#define DISPLAY_H

#include "real.h"
#include "array.h"
#include "mecable.h"
#include "display_prop.h"
#include "property_list.h"

class Simul;
class Mecable;
class SingleSet;
class CoupleSet;
class Fiber;
class FiberSet;
class Solid;
class SolidSet;
class Organizer;
class OrganizerSet;
class Space;
class SpaceSet;
class Sphere;
class SphereSet;
class Bead;
class BeadSet;
class FieldSet;
class FiberProp;
class FiberDisp;
class PointDisp;
class LineDisp;

/// defining the DISPLAY keyword enables display code in included files
#define DISPLAY


///Base class to display Cytosim state using OpenGL
class Display
{
protected:
    
    /**
     @brief A display element with a depth coordinate
     
     zObject is used to depth-sort and display transparent objects
     */
    class zObject
    {
        /// position of the center of gravity
        GLfloat        zPos;
        
        /// pointer to object
        Mecable const* mObj;
        
        /// index of specific point
        unsigned       mIdx;
    public:
        
        zObject()
        {
            mObj = 0;
            mIdx = 0;
        }
        
        zObject(Mecable const* o, unsigned i = 0)
        {
            mObj = o;
            mIdx = i;
        }
        
        /// position
        Vector position() const { return mObj->posPoint(mIdx); }
        
        /// height
        GLfloat depth()   const { return zPos; }
        
        /// set height
        void depth(real z)      { zPos = z; }
        
        ///display
        void display(Display*);

    };
    
    /// function to sort zObjects according to their position 'depth'
    static int closer(const void * ap, const void * bp)
    {
        zObject const* a = static_cast<const zObject*>(ap);
        zObject const* b = static_cast<const zObject*>(bp);
        
        if ( a->depth() > b->depth() ) return  1;
        if ( a->depth() < b->depth() ) return -1;
        return 0;
    }

    /// array of transparent objects to be displayed last
    Array<zObject> zObjects;

private:
    
    /// set default value of FiberProp
    void prepareFiberDisp(FiberProp*, PropertyList&, gle_color);

    /// set values of fiber's LineDisp
    void prepareLineDisp(Fiber const*);

    template < typename T >
    void preparePointDisp(T * prop, PropertyList&, gle_color);

protected:
    
    /// the pixel size for this particular display
    real           pixelSize;
    
    /// scaling factors to convert 'size' parameter into pixels used in glPointSize() and glLineWidth()
    real           uFactor;
    
    /// scaling factors to convert 'size' parameter into real dimensions used in glScale()
    real           sFactor;
    
    
    /// flag used to calculate clusterAnalysis only once
    unsigned       fiber_prep;
    
    /// used to calculate clusterAnalysis only once
    real           fiber_prep_time;
    
    /// min and max age used to adjust color range with COLORING_AGE
    real           age_range, age_min;

public:
    
    /// associated parameters
    DisplayProp const* prop;
    
    ///constructor
    Display(DisplayProp const*);
    
    ///destructor
    virtual ~Display() {}
    
    
    ///display opaque internal objects using OpenGL commands
    virtual void displayObjects(Simul const&);

    /// set current pixel-size and the value of the point in pixels
    void setPixelFactors(real pixel_size, real uFactor);
    
    ///prepare to display
    void prepareForDisplay(Simul const&, PropertyList&);
    
    ///display the whole simulation
    void display(Simul const&);
    
    ///display for periodic systems
    void displayTiled(Simul const&, int nine);
    
        
    /// find an individual color
    void bodyColor(PointDisp const*, ObjectID) const;

    /// find an individual color that may be transparent
    void bodyColorT(PointDisp const*, ObjectID) const;

    /// set OpenGL line width
    void lineWidth(real w) const { glLineWidth(std::max(w*uFactor, 0.25)); }

    /// set OpenGL point size
    void pointSize(real w) const { glPointSize(std::max(w*uFactor, 0.25)); }
    
    
    ///display a scalar field
    void displayFields(FieldSet const&);
    
    ///display all Spaces (in 3D the back side)
    void displaySpaces(SpaceSet const&);
    
    ///display the front-side of Spaces in 3D
    void displaySpacesF(SpaceSet const&);

    
    ///display Fiber MINUS_END
    virtual void displayFiberEndMinus(Fiber const&) const;
    
    ///display Fiber PLUS_END
    virtual void displayFiberEndPlus(Fiber const&) const;

    
    ///display Fiber linear features
    virtual void displayFiberLines(Fiber const&) const;
    
    ///display Fiber linear features located between [MINUS_END, abs]
    virtual void displayFiberLinesMinus(Fiber const&, real abs, real width) const;
    
    ///display Fiber linear features located between [abs, PLUS_END]
    virtual void displayFiberLinesPlus(Fiber const&, real abs, real width) const;

    /// actin-like rendering using a sphere to represent each monomer
    void         displayActin(Fiber const& fib, gle_color const&, gle_color const&, gle_color const&) const;
    
    /// microtubule-like rendering using a sphere to represent each monomer
    void         displayMicrotubule(Fiber const& fib, gle_color const&, gle_color const&, gle_color const&) const;
    
    ///display Fiber point-like features
    virtual void displayFiberPoints(Fiber const&) const;
    
    ///display Fiber Speckles
    virtual void displayFiberSpeckles(Fiber const&) const;
   
    /// display lattice subtance using color
    virtual void displayFiberLattice(Fiber const&) const;
    
    /// display Labels for a Fiber
    void         displayFiberLabels(Fiber const&) const;
    
    /// display forces acting on each point
    void         displayFiberForces(Fiber const&) const;
    
    ///display all features of Fiber
    virtual void displayFiber(Fiber const&);
    
    ///display Fibers
    virtual void displayFibers(FiberSet const&);
    
    ///display the average fiber for the pool defined by func(obj, val) == true
    void displayAverageFiber(ObjectList const&);
    
    ///display the averaged fiber
    void displayAverageFiber1(FiberSet const&, void const*);
    
    ///display the average for left-pointing and right-pointing fibers
    void displayAverageFiber2(FiberSet const&, void const*);

    
    ///display a Bead
    virtual void displayBead(Bead const&) = 0;

    ///display translucent elements of a Bead
    virtual void displayBeadT(Bead const&) = 0;
    
    ///display the Beads
    void displayBeads(BeadSet const&);

    
    ///display opaque elements of a Solid
    virtual void displaySolid(Solid const&) = 0;
    
    ///display translucent elements of a Solid
    virtual void displaySolidT(Solid const&, unsigned int) = 0;
    
    ///display the Solids
    void displaySolids(SolidSet const&);
    
    
    ///display the Sphere
    virtual void displaySphere(Sphere const&) = 0;

    ///display translucent elements of a Sphere
    virtual void displaySphereT(Sphere const&) = 0;
    
    ///display the Spheres
    void displaySpheres(SphereSet const&);
    
    
    ///display the free Singles
    virtual void displaySinglesF(SingleSet const&) = 0;
    
    ///display the attached Singles
    virtual void displaySinglesA(SingleSet const&) = 0;

    ///display the free Couples, showing Hand1
    virtual void displayCouplesF1(CoupleSet const&) = 0;
    
    ///display the free Couples, showing Hand1 or Hand2
    virtual void displayCouplesF2(CoupleSet const&) = 0;
    
    ///display the free Couples
    void displayCouplesF(CoupleSet const&);
    
    ///display the attached Couples
    virtual void displayCouplesA(CoupleSet const&) = 0;
    
    ///display the bridging Couples
    virtual void displayCouplesB(CoupleSet const&) = 0;

    
    ///display Organizer
    virtual void displayOrganizer(Organizer const&) = 0;
    
    ///display the Organizers
    void displayOrganizers(OrganizerSet const&);


    ///display translucent objects after depth-sorting
    void displayTransparentObjects(Array<zObject>&);

    
    ///display additional items
    void displayMisc(Simul const&);
    
};


#endif

