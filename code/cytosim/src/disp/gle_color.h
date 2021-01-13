// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef GLE_COLOR_H
#define GLE_COLOR_H

#include "opengl.h"
#include <iostream>

/**
 gle_color implements colors with 4-components:
 - Red
 - Green
 - Blue
 - Alpha = transparency
 .
 
 This class implements the `RGBA` format using an 'unsigned integer'
 and an array of 4 GLfloat.
 
 F. Nedelec -- Merged two older color classes on 23 August 2015
 */
/// Color with 4 components: red, green, blue, alpha (RGBA)
class gle_color
{
    static uint32_t combine(GLubyte r, GLubyte g, GLubyte b, GLubyte a)
    {
        return ( uint32_t(r) << 24 ) | ( uint32_t(g) << 16 ) | ( uint32_t(b) << 8 ) | uint32_t(a);
    }
    
private:
        
    /// 32-bits integer containing 4 one-byte components
    uint32_t rgba;
    
    /// array of 4 float components
    GLfloat col[4];
    
public:
    
    /// set to white
    void set_white()
    {
        rgba=0xFFFFFFFF;
        col[0]=1;
        col[1]=1;
        col[2]=1;
        col[3]=1;
    }
    
    /// set to black
    void set_black()
    {
        rgba=0x000000FF;
        col[0]=0;
        col[1]=0;
        col[2]=0;
        col[3]=1;
    }

    /// specify floating point components
    void set_floats(GLfloat r, GLfloat g, GLfloat b, GLfloat a)
    {
        col[0] = r < 1.0 ? r : 1.0;
        col[1] = g < 1.0 ? g : 1.0;
        col[2] = b < 1.0 ? b : 1.0;
        col[3] = a < 1.0 ? a : 1.0;
        rgba = combine(col[0]*255, col[1]*255, col[2]*255, col[3]*255);
    }
    
    /// export floating point components
    void put_floats(GLfloat& r, GLfloat& g, GLfloat& b, GLfloat& a)
    {
        r = col[0];
        g = col[1];
        b = col[2];
        a = col[3];
    }
    
    /// export floating point components
    void put_floats(GLfloat* c)
    {
        c[0] = col[0];
        c[1] = col[1];
        c[2] = col[2];
        c[3] = col[3];
    }

    /// specify components with bytes
    void set_bytes(GLubyte r, GLubyte g, GLubyte b, GLubyte a)
    {
        rgba = combine(r, g, b, a);
        col[0] = r / 255.f;
        col[1] = g / 255.f;
        col[2] = b / 255.f;
        col[3] = a / 255.f;
    }

    /// export components as bytes
    void put_bytes(GLubyte& r, GLubyte& g, GLubyte& b, GLubyte& a)
    {
        r = 0xFF & ( rgba >> 24 );
        g = 0xFF & ( rgba >> 16 );
        b = 0xFF & ( rgba >> 8 );
        a = 0xFF & ( rgba );
    }
    
    /// set from 4-bytes integer
    void set_int(uint32_t u)
    {
        rgba = u;
        col[0] = ( 0xFF & ( u >> 24 ) ) / 255.f;
        col[1] = ( 0xFF & ( u >> 16 ) ) / 255.f;
        col[2] = ( 0xFF & ( u >>  8 ) ) / 255.f;
        col[3] = ( 0xFF & u ) / 255.f;
   }

    /// default constructor
    gle_color() : rgba(0)
    {
        col[0] = 0;
        col[1] = 0;
        col[2] = 0;
        col[3] = 0;
    }
    
    /// constructor
    gle_color(const uint32_t& u)
    {
        set_int(u);
    }
    
    /// constructor with Alpha component = 1.0
    gle_color(const GLfloat& r, const GLfloat& g, const GLfloat& b)
    {
        set_floats(r,g,b,1.0f);
    }

    /// constructor
    gle_color(const GLfloat& r, const GLfloat& g, const GLfloat& b, const GLfloat& a)
    {
        set_floats(r,g,b,a);
    }

    void operator = (const uint32_t& col)
    {
        set_int(col);
    }
    
    bool operator ==(const gle_color col) const { return rgba == col.rgba; }
    bool operator !=(const gle_color col) const { return rgba != col.rgba; }
    
    GLfloat const* data() const { return col; }
    
    GLfloat red()         const { return col[0]; }
    GLfloat green()       const { return col[1]; }
    GLfloat blue()        const { return col[2]; }
    GLfloat alpha()       const { return col[3]; }

    GLfloat brightness()  const { return ( col[0] + col[1] + col[2] ) * col[3]; }
    bool    opaque()      const { return ( (rgba & 0xFF) == 0xFF ); }
    bool    transparent() const { return ( (rgba & 0xFF) != 0xFF ); }
    bool    visible()     const { return ( rgba & 0xFF ); }
    bool    invisible()   const { return ( rgba & 0xFF ) == 0; }
    
    gle_color darken(GLfloat s) const
    {
        if ( s < 1 )
            return gle_color(s*col[0], s*col[1], s*col[2], col[3]);
        else
            return gle_color(*this);
    }
    
    gle_color lighten(GLfloat s) const
    {
        GLfloat r = s * col[0];
        GLfloat g = s * col[1];
        GLfloat b = s * col[2];
        return gle_color(r<1.0?r:1.0, g<1.0?g:1.0, b<1.0?b:1.0, col[3]);
    }
    
    gle_color alpha_scaled(GLfloat s) const
    {
        if ( s < 1 )
            return gle_color(col[0], col[1], col[2], s*col[3]);
        else
            return gle_color(*this);
    }
    
    gle_color alpha_set(GLfloat s) const
    {
        if ( s < 1 )
            return gle_color(col[0], col[1], col[2], s);
        else
            return gle_color(col[0], col[1], col[2], 1.0);
    }
    
    gle_color alpha_one() const
    {
        return gle_color(col[0], col[1], col[2], 1.0);
    }
    
    gle_color blend(gle_color c) const
    {
        GLfloat s = alpha() + c.alpha();
        GLfloat h = alpha()   / s;
        GLfloat g = c.alpha() / s;
        return gle_color(h*col[0]+g*c.col[0], h*col[1]+g*c.col[1], h*col[2]+g*c.col[2], h+g);
    }
    
    gle_color inverted() const
    {
        return gle_color(1.f-col[0], 1.f-col[1], 1.f-col[2], col[3]);
    }
    
    
    /// set current OpenGL color by calling glColor
    void load() const
    {
        //std::clog << "OpenGL load " << col[0] << " " << col[1] << " " << col[2] << " " << col[3] << "\n";
        glColor4fv(col);
        //glColor4f(col[0], col[1], col[2], col[3]);
    }
    
    /// set current OpenGL color, but with `s` as alpha component
    void load(GLfloat s) const
    {
        if ( s < 1.0 )
            glColor4f(col[0], col[1], col[2], s);
        else
            glColor4fv(col);
    }
    
    /// set current OpenGL color, to *this darkened by factor `s`
    void load_darken(GLfloat s) const
    {
        if ( s < 1.0 )
            glColor4f(s*col[0], s*col[1], s*col[2], col[3]);
        else
            glColor4fv(col);
    }
    
    
    void load_clear() const
    {
        glClearColor(col[0], col[1], col[2], col[3]);
    }
    
    /// set FRONT material property for lighting
    void load_front() const
    {
        GLfloat blk[4] = { 0, 0, 0, 1 };
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, col);
        glMaterialfv(GL_FRONT, GL_EMISSION, blk);
    }
    
    /// set FRONT material property for lighting, and current color
    void load_load() const
    {
        GLfloat blk[4] = { 0, 0, 0, 1 };
        glColor4fv(col);
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, col);
        glMaterialfv(GL_FRONT, GL_EMISSION, blk);
    }
    
    /// set FRONT material property for lighting, and current color
    void load_load(GLfloat s) const
    {
        GLfloat blk[4] = { 0, 0, 0, 1 };
        GLfloat mat[4] = { col[0], col[1], col[2], s < 1.0 ? s : 1.0 };
        glColor4fv(mat);
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat);
        glMaterialfv(GL_FRONT, GL_EMISSION, blk);
    }

    /// set front OpenGL color, but with `s` as alpha component
    void load_front(GLfloat s) const
    {
        GLfloat blk[4] = { 0, 0, 0, 1 };
        GLfloat mat[4] = { col[0], col[1], col[2], s < 1.0 ? s : 1.0 };
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat);
        glMaterialfv(GL_FRONT, GL_EMISSION, blk);
    }

    /// set BACK material property for lighting
    void load_back() const
    {
        GLfloat blk[4] = { 0, 0, 0, 1 };
        glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, blk);
        glMaterialfv(GL_BACK, GL_EMISSION, col);
    }
    
    /// set BACK material property for lighting
    void load_back(GLfloat s) const
    {
        GLfloat blk[4] = { 0, 0, 0, 1 };
        GLfloat mat[4] = { col[0], col[1], col[2], s < 1.0 ? s : 1.0 };
        glMaterialfv(GL_BACK, GL_EMISSION, blk);
        glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, mat);
    }
    
    /// set FRONT and BACK material property for lighting
    void load_both() const
    {
        GLfloat blk[4] = { 0, 0, 0, 1 };
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, col);
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, blk);
        //glMaterialfv(GL_BACK, GL_EMISSION, col);
    }
    
    
    /// set FRONT and BACK material property for lighting
    void load_emission() const
    {
        GLfloat blk[4] = { 0, 0, 0, 1 };
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, blk);
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, col);
        //glMaterialfv(GL_BACK, GL_EMISSION, col);
    }

    /// set FRONT and BACK material property for lighting
    void load_both(GLfloat s) const
    {
        GLfloat blk[4] = { 0, 0, 0, 1 };
        GLfloat mat[4] = { col[0], col[1], col[2], s < 1.0 ? s : 1.0 };
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat);
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, blk);
        //glMaterialfv(GL_BACK, GL_EMISSION, mat);
    }


    /// conversion to string
    std::string to_string() const;
    
    /// 
    static void print_materials(std::ostream & os)
    {
        GLfloat mat[4] = { 0 };
        glGetMaterialfv(GL_FRONT, GL_AMBIENT, mat);
        os << "front  ambient "<<mat[0]<<" "<<mat[1]<<" "<<mat[2]<<" "<<mat[3]<<std::endl;
        glGetMaterialfv(GL_FRONT, GL_DIFFUSE, mat);
        os << "front  diffuse "<<mat[0]<<" "<<mat[1]<<" "<<mat[2]<<" "<<mat[3]<<std::endl;
        glGetMaterialfv(GL_FRONT, GL_EMISSION, mat);
        os << "front emission "<<mat[0]<<" "<<mat[1]<<" "<<mat[2]<<" "<<mat[3]<<std::endl;

        glGetMaterialfv(GL_BACK, GL_AMBIENT, mat);
        os << "back   ambient "<<mat[0]<<" "<<mat[1]<<" "<<mat[2]<<" "<<mat[3]<<std::endl;
        glGetMaterialfv(GL_BACK, GL_DIFFUSE, mat);
        os << "back   diffuse "<<mat[0]<<" "<<mat[1]<<" "<<mat[2]<<" "<<mat[3]<<std::endl;
        glGetMaterialfv(GL_BACK, GL_EMISSION, mat);
        os << "back  emission "<<mat[0]<<" "<<mat[1]<<" "<<mat[2]<<" "<<mat[3]<<std::endl;
    }
};



/// input operator:
std::istream & operator >> (std::istream&, gle_color&);

/// output operator:
std::ostream & operator << (std::ostream&, const gle_color&);



namespace gle
{
    /// conversion function from RGB to HSV color space
    void RGB2HSV(GLfloat r, GLfloat g, GLfloat b, GLfloat* h, GLfloat* s, GLfloat* v);
    
    /// conversion functions from HSV to RGB color space
    void HSV2RGB(GLfloat h, GLfloat s, GLfloat v, GLfloat* r, GLfloat* g, GLfloat* b);
    
    
    /// set a RGB color from a factor in [-PI, PI], continuously varying through all colors
    void set_hue_color(GLfloat& r, GLfloat& g, GLfloat& b, GLfloat h);
    
    /// return saturated color with given Hue in [-PI, PI]
    gle_color hue_color(GLfloat h, GLfloat alpha);
    
    /// set a RGB color from a factor in [0, 1], continuously varying between blue, green, red
    void set_jet_color(GLfloat& r, GLfloat& g, GLfloat& b, GLfloat h, GLfloat hm);

    /// return jet color
    gle_color jet_color(GLfloat h, GLfloat alpha);

    /// set a RGB color from a factor in [0, 1], continuously varying through blue, green, red, yellow, white
    void set_jet_colorE(GLfloat& r, GLfloat& g, GLfloat& b, GLfloat h, GLfloat hm);
    
    /// return jet color2
    gle_color jet_colorE(GLfloat h, GLfloat alpha);
}

#endif
