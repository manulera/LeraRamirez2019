// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef POINT_DISP_H
#define POINT_DISP_H

#include "gle_color.h"
#include "property.h"
#include "vector.h"
#include "gle.h"

class Glossary;

#define POINTDISP_USES_PIXELMAPS


/// the parameters necessary to display a point-like object
/**
 Note that these parameters may be interpreted differently when displaying different classes. 
 For example `coloring` is only implemented for Sphere, Beads but not for Hands,
 while `shape` and `symbol` are only implemented for Hands.
 */
class PointDisp : public Property
{    
    /// used to differentiate between different uses of the class
    std::string mKind;
    
    /// size of feature in pixels
    unsigned    pixSize;
    
    /// draw basic shape
    void strokeShape() const;

    /// draw active state with OpenGL vector primitives
    void strokeA() const;
    
    /// draw inactive state with OpenGL vector primitives
    void strokeI() const;

#ifdef POINTDISP_USES_PIXELMAPS
    
    /// pointer to 3 square bitmaps with 4*nPix*nPix pixels each
    GLubyte   *bmp[3];

    /// index of the Pixel Buffer Objects on GPU
    GLuint     pbo[3];
    
    /// center of bitmap
    GLfloat    mOffs;
    
    /// allocated size of bitmap
    unsigned   nPix;

    /// allocate pixelmap memory
    void allocatePixelmap();
    
    /// release pixelmap memory
    void releasePixelmap();
    
    /// scale down pixelmap by factor 'bin'
    void downsampleRGBA(GLubyte*, unsigned, unsigned, GLubyte const*, unsigned bin);

    /// draw pixelmap
    void savePixelmap(GLubyte*, unsigned dim, GLuint);
    
    /// create the pixelmaps
    void makePixelmaps(real, unsigned supersampling);
    
    /// draw pixel map
    void drawPixelmap(Vector const& pos, unsigned ii) const;

#endif
    
    /// clear pointers
    void clearPixelmaps();

public:
    
    /**
     @defgroup PointDispPar Display parameters: Points
     @ingroup DisplayParameters
     @{
     */
    
    
    /// visibility flag : 0=hidden, 1=opaque
    /**
     For a Space, bit 1 controls display of the back side, and bit 2 the front side:
     and you can use visible=1 to display only the back, 2 to display only the front.
     */
    int        visible;
    
    /// color of object (in 3D display, the color of outer surfaces)
    gle_color  color;
    
    /// second color (set as color[1])
    /**
     This is used to display unattached Single and unbridging Couple, 
     and the inner surfaces of objects such as Sphere, Solid, Bead and Space.
     Unless it is set directly as color[1], `color2` is set to be a darker tone of `color`.
     */
    gle_color  color2;
    
    /// if true, use various colors to display different objects
    int        coloring;
    
    /// display diameter of points
    real       size;
    
    /// display width of lines
    real       width;
    
    /// 'c' for circle, 'h' for hexagon, 's' for star, etc.
    char       shape;
    
    /// a bitfield to set different display options
    int        style;
    
    /// character displayed (do not set, or set as 0 to disable this feature)
    char       symbol;
    
    /// color of symbol (set as symbol[1])
    gle_color  symbol_color;
    
    /// @}
    
    /// visible and big enough to be seen
    bool       perceptible;
    
    /// this is size in real unit
    real       realSize;
    
public:
    
    /// constructor
    PointDisp(const std::string& k, const std::string& n);
    
    /// copy constructor
    PointDisp(PointDisp const&);
    
    /// copy assignment
    PointDisp& operator =(PointDisp const&);

    /// destructor
    ~PointDisp();
    
    /// identifies the property
    std::string category() const { return mKind; }

    /// clear to default values
    void clear();
    
    /// set from glossary
    void read(Glossary&);
    
    /// return a carbon copy of object
    Property* clone() const { return new PointDisp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
    /// recalculate bitmaps
    void      prepare(real uf, real sf);

#ifdef POINTDISP_USES_PIXELMAPS
    
    /// draw inactive state
    void      drawI(Vector const& pos) const { if ( perceptible ) drawPixelmap(pos, 0); }

    /// draw active unbound state
    void      drawF(Vector const& pos) const { if ( perceptible ) drawPixelmap(pos, 1); }
    
    /// draw active bound state
    void      drawA(Vector const& pos) const { if ( perceptible ) drawPixelmap(pos, 2); }

#else
    
    /// draw inactive state
    void      drawI(Vector const& pos) const;

    /// draw active state
    void      drawF(Vector const& pos) const;
    
    /// draw active state
    void      drawA(Vector const& pos) const;

#endif

};


#endif
