// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_prop.h"
#include "filepath.h"
#include "glossary.h"
#include "property_list.h"
#include "simul_prop.h"
#include "simul.h"
#include "sim.h"

#include "space.h"
#include "space_force.h"
#include "space_square.h"
#include "space_sphere.h"
#include "space_twinspheres.h"
#include "space_polygon.h"
#include "space_polygonZ.h"
#include "space_capsule.h"
#include "space_banana.h"
#include "space_torus.h"
#include "space_dice.h"
#include "space_strip.h"
#include "space_periodic.h"
#include "space_ellipse.h"
#include "space_cylinder.h"
#include "space_cylinderZ.h"
#include "space_cylinderP.h"
#include "space_ring.h"
#include "space_beads.h"
#include "space_tee.h"

#ifdef NEW_SPACES
#include "space_inflate.h"
#include "space_combine.h"
#include "space_rotate.h"
#endif

#ifdef NEW_DYNAMIC_SPACES
#include "space_disc.h"
#include "space_dynamic_sphere.h"
#include "space_dynamic_ellipse.h"
#endif

/**
 @defgroup SpaceGroup Space and Geometry
 @ingroup ObjectGroup
 @ingroup NewObject
 @brief A Space defines a confined region
 
 A Space is created by specifying a geometry:
 @code
 set space NAME
 {
    geometry = GEOMETRY DIMENSIONS
 }
 @endcode
 
 DIMENSIONS is usually a list of numbers.
 
 List of known `geometry`:
 
 GEOMETRY      |   Class              | DIMENSIONS
 --------------|----------------------|-----------------------------------
 `rectangle`   | SpaceSquare          | sizeX sizeY sizeZ
 `sphere`      | SpaceSphere          | radius
 `twin_spheres`| SpaceTwinSpheres     | radius, radius, overlap
 `polygon`     | SpacePolygon         | file_name height
 `polygonZ`    | SpacePolygonZ        | file_name
 `capsule`     | SpaceCapsule         | half_length radius
 `torus`       | SpaceTorus           | radius thickness
 `banana`      | SpaceBanana          | total_length width radius_of_curvature
 `dice`        | SpaceDice            | sizeX sizeY sizeZ radius
 `strip`       | SpaceStrip           | sizeX sizeY sizeZ
 `periodic`    | SpacePeriodic        | sizeX sizeY sizeZ
 `ellipse`     | SpaceEllipse         | sizeX sizeY sizeZ
 `cylinder`    | SpaceCylinder        | half_length radius
 `cylinderZ`   | SpaceCylinderZ       | radius bottom top
 `cylinderP`   | SpaceCylinderP       | half_length radius
 `ring`        | SpaceRing            | half_length radius
 `tee`         | SpaceTee             | half_length radius arm_position arm_length
 `beads`       | SpaceBeads           | custom class

 Dynamic Space with variable geometry:
 
 GEOMETRY           |   Class              | DIMENSIONS
 -------------------|----------------------|-------------------------
 `disc`             | SpaceDisc            | radius
 `dynamic_sphere`   | SpaceDynamicSphere   | radius
 `dynamic_ellipse`  | SpaceDynamicEllipse  | sizeX sizeY sizeZ
 
 Example:
 @code
 set space cell
 {
   geometry = sphere 5
 }
 @endcode
 */
Space * SpaceProp::newSpace(SpaceProp const* sp, Glossary& opt)
{
    const std::string s = sp->shape;
    
    if ( s=="rectangle" || s=="square" )       return new SpaceSquare(sp);
    if ( s=="circle" || s=="sphere" )          return new SpaceSphere(sp);
    if ( s=="twin_spheres" )                   return new SpaceTwinSpheres(sp);
    if ( s=="polygon" )                        return new SpacePolygon(sp, opt);
    if ( s=="polygonZ" )                       return new SpacePolygonZ(sp, opt);
    if ( s=="capsule" || s=="spherocylinder" ) return new SpaceCapsule(sp);
    if ( s=="banana" )                         return new SpaceBanana(sp);
    if ( s=="torus" )                          return new SpaceTorus(sp);
    if ( s=="dice" )                           return new SpaceDice(sp);
    if ( s=="strip" )                          return new SpaceStrip(sp);
    if ( s=="periodic" )                       return new SpacePeriodic(sp);
    if ( s=="ellipse" || s=="ellipsoid" )      return new SpaceEllipse(sp);
#if ( DIM == 3 )
    if ( s=="cubic" )                          return new SpaceSquare(sp);
    if ( s=="cylinder" )                       return new SpaceCylinder(sp);
    if ( s=="cylinderZ" )                      return new SpaceCylinderZ(sp);
    if ( s=="cylinderP" )                      return new SpaceCylinderP(sp);
#elif ( DIM == 2 )
    if ( s=="cylinder" )                       return new SpaceSquare(sp);
    if ( s=="cylinderP" )                      return new SpaceStrip(sp);
#else
    if ( s=="cylinder" )                       return new SpaceSquare(sp);
    if ( s=="cylinderP" )                      return new SpacePeriodic(sp);
#endif
    if ( s=="ring" )                           return new SpaceRing(sp);
    if ( s=="tee" )                            return new SpaceTee(sp);
#ifdef NEW_SPACES
    if ( s=="force" )                          return new SpaceForce(sp, opt);
    if ( s=="beads" )                          return new SpaceBeads(sp);
#endif
#ifdef NEW_DYNAMIC_SPACES
    if ( s=="disc" )                           return new SpaceDisc(sp);
    if ( s=="dynamic_sphere" )                 return new SpaceDynamicSphere(sp);
    if ( s=="dynamic_ellipse" )                return new SpaceDynamicEllipse(sp);
    // backward compatibility:
    if ( s=="contractile" )                    return new SpaceDynamicEllipse(sp);
#endif
    return 0;
}


Space * SpaceProp::newSpace(Glossary& opt) const
{
    Space * spc = newSpace(this, opt);
    
    if ( spc == 0 )
        throw InvalidParameter("unknown space:shape `"+shape+"'");
    
    // set dimensions:
    if ( dimensions.length() )
        spc->readLengths(dimensions);
    
    std::string dim;
    if ( opt.set(dim, "dimensions") )
        spc->readLengths(dim);

#if ( 0 )
    real len;
    if ( opt.set(len, "inflate") && len > 0 )
    {
        spc = new SpaceInflate(this, spc, len);
    }
#endif
    
    return spc;
}


//------------------------------------------------------------------------------

void SpaceProp::clear()
{
    geometry   = "";
    shape      = "none";
    dimensions = "";
    shape_spec = "";
    display    = "";
    
#ifdef NEW_DYNAMIC_SPACES
    viscosity  = INFINITY;
    viscosity_rot = INFINITY;
    tension    = 0;
    volume     = 0;
#endif
    
    dimensions_old = "";
}

void SpaceProp::read(Glossary& glos)
{    
    glos.set(shape,        "shape");
    glos.set(dimensions,   "dimension") || glos.set(dimensions, "dimensions");
#ifdef BACKWARD_COMPATIBILITY
    glos.set(dimensions,   "spec");  // format 36
#endif
    glos.set(geometry,     "geometry"); //deprecated as of 10.2015
    
#ifdef NEW_DYNAMIC_SPACES
    glos.set(tension,       "tension");
    glos.set(volume,        "volume");
    glos.set(viscosity,     "viscosity");
    glos.set(viscosity_rot, "viscosity", 1);
#endif
    
    if ( glos.set(display, "display") )
        display_fresh = true;
}

//------------------------------------------------------------------------------

void SpaceProp::complete(Simul const* sim)
{
    if ( !geometry.empty() )
    {
        std::istringstream iss(geometry);
        iss >> shape;
        
        if ( iss.fail() )
            throw InvalidParameter("invalid geometry `"+geometry+"' for Space");
        
        char c;
        while ( isspace(iss.peek()) )
            iss.get(c);
        
        c = iss.peek();
        if ( !isdigit(c) && c != '+' && c != '-' )
        {
            iss >> shape_spec;
        }
        
        while ( isspace(iss.peek()) )
            iss.get(c);

        // get remaining characters as a whole:
        if ( iss.good() )
            dimensions = geometry.substr(iss.tellg());
    }

#ifdef NEW_DYNAMIC_SPACES
    if ( viscosity > 0 )
        mobility_dt = sim->prop->time_step / viscosity;
    else
        throw InvalidParameter("space:viscosity must be > 0");
    
    if ( viscosity_rot > 0 )
        mobility_rot_dt = sim->prop->time_step / viscosity_rot;
    else
        throw InvalidParameter("space:viscosity[1] (rotational viscosity) must be > 0");
#endif

    /*
     If the dimensions have changed, update any Space with this property.
     This is necessary to make 'change space:dimension' work.
     */
    if ( sim && dimensions != dimensions_old )
    {
        for ( Space * spc = sim->spaces.first(); spc; spc=spc->next() )
        {
            if ( spc->prop == this )
            {
                spc->readLengths(dimensions);
                // allow Simul to update:
                if ( spc == sim->space() )
                    const_cast<Simul*>(sim)->changeSpace(spc);
            }
        }
        dimensions_old = dimensions;
    }
}

//------------------------------------------------------------------------------

void SpaceProp::write_values(std::ostream & os) const
{
    //write_value(os, "geometry",   geometry);
    write_value(os, "shape",      shape);
    write_value(os, "dimensions", dimensions);
#ifdef NEW_DYNAMIC_SPACES
    write_value(os, "tension",    tension);
    write_value(os, "volume",     volume);
    write_value(os, "viscosity",  viscosity, viscosity_rot);
#endif
    write_value(os, "display",    "("+display+")");
}



