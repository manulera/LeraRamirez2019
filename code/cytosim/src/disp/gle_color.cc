// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "gle_color.h"
#include "gle_color_list.h"
#include "smath.h"
#include "stream_func.h"
#include "exceptions.h"
#include <iomanip>
#include <cctype>

#pragma mark Input/Output

char int2hex(int d)
{
    if ( d > 9 )
        return 'A' + (d-10);
    else
        return '0' + d;
}


std::string gle_color::to_string() const
{
#if ( 0 )
    std::ostringstream oss;
    oss << "0x" << std::setw(8) << std::setfill('0') << std::hex << std::uppercase << col.rgba;
    return oss.str();
#else
    char str[16] = { 0 };
    str[0] = '0';
    str[1] = 'x';
    uint32_t n = rgba;
    for ( int c = 0; c < 8; ++c )
    {
        uint32_t d = ( n >> 28 );
        n <<= 4;
        str[c+2] = int2hex(d);
    }
    return str;
#endif
}


int hex2int(char c)
{
    if ( isalpha(c) )
        return tolower(c) - 'a' + 10;
    else
        return c - '0';
}

int hex2int(char a, char b)
{
    return hex2int(a)*16 + hex2int(b);
}

/**
 A color is composed of 4 components (Red, Green, Blue, Alpha),
 and can be specified in different ways:
 -# with a hexadecimal integer: 0xFF0000FF  or  0xff0000ff
 -# with 3 or 4 floats: (1 0 0) or (1 0 0 1)
 -# with a name:  red
 -# with a number: #1
 .
*/
std::istream & operator >> (std::istream& is, gle_color& col)
{
    col.set_white();

    std::istream::sentry s(is);
    if ( s )
    {
        char c = is.get();
        char d = is.peek();
        
        if ( isalpha(c) )
        {
            is.unget();
            std::string name;
            is >> name;
            try {
                col = gle::std_color(name);
            }
            catch ( InvalidSyntax & e )
            {
                std::cerr << e.what() << " : using white" << std::endl;
                gle::print_std_colors(std::cerr);
                col = 0xFFFFFFFF;
            }
        }
        else if ( '#'==c  &&  isdigit(d) )
        {
            unsigned int nb;
            is >> nb;
            col = gle::alt_color(nb);
        }
        else if ( '0'==c  &&  'x'==d )
        {
            is.get();
            GLubyte u[4] = { 255 };
            int i = 0;
            while ( i < 4 )
            {
                is >> c >> d;
                if ( is.fail() )
                {
                    if ( i < 3 )
                        throw InvalidSyntax("incomplete hexadecimal color specification");
                    is.clear();
                    break;
                }
                if ( !isxdigit(c) || !isxdigit(d) )
                    throw InvalidSyntax("invalid hexadecimal digit");
                u[i++] = hex2int(c, d);
            }
            col.set_bytes(u[0], u[1], u[2], u[3]);
        }
        else if ( isdigit(c) )
        {
            is.unget();
            GLfloat r, g, b, a=1;
            is >> r >> g >> b;
            if ( ! is.fail() )
            {
                is >> a;
                is.clear();
                col.set_floats(r,g,b,a);
            }
        }
    }
    return is;
}


std::ostream& operator << (std::ostream& os, gle_color const& x)
{
    os << x.to_string();
    return os;
}

//-----------------------------------------------------------------------
//=========================  RAINBOW EFFECTS  ===========================
//-----------------------------------------------------------------------
#pragma mark -


/**
 r,g,b values are from 0 to 1
 h = [0, 360], s = [0,1], v = [0,1]
 if s == 0, then h = -1 (undefined)
*/

void gle::RGB2HSV(const GLfloat r, const GLfloat g, const GLfloat b, GLfloat* h, GLfloat* s, GLfloat* v)
{
    GLfloat mn, mx, delta;
    mn = std::min(r, std::min(g, b));
    mx = std::max(r, std::max(g, b));
    *v = mx;
    delta = mx - mn;
    if ( mx != 0 )
        *s = delta / mx;
    else {
        *s = 0;
        *h = -1;
        return;
    }
    if ( r == mx )
        *h = ( g - b ) / delta;       // between yellow & magenta
    else if ( g == mx )
        *h = 2 + ( b - r ) / delta;   // between cyan & yellow
    else
        *h = 4 + ( r - g ) / delta;   // between magenta & cyan
    *h *= 60;                         // degrees
    if ( *h < 0 )
        *h += 360;
}



void gle::HSV2RGB( const GLfloat h, const GLfloat s, const GLfloat v, GLfloat* r, GLfloat* g, GLfloat* b )
{
    int i;
    GLfloat f, p, q, t;
    if ( s == 0 ) {
        // achromatic (gray)
        *r = *g = *b = v;
        return;
    }
    GLfloat hc = h/60;               // sector 0 to 5
    i = (int)floor( hc );
    f = hc - i;                      // factorial part of h
    p = v * ( 1 - s );
    q = v * ( 1 - s * f );
    t = v * ( 1 - s * ( 1 - f ) );
    switch( i )
    {
        case 0:  *r = v; *g = t; *b = p; break;
        case 1:  *r = q; *g = v; *b = p; break;
        case 2:  *r = p; *g = v; *b = t; break;
        case 3:  *r = p; *g = q; *b = v; break;
        case 4:  *r = t; *g = p; *b = v; break;
        case 5:  *r = v; *g = p; *b = q; break;
        default: *r = 1; *g = 1; *b = 1; break;
    }
}



/**
 set a RGB color as a function of a Hue value `a` in [-PI, PI].
 The colors follow in this order: red, green, blue, red ...
*/
void gle::set_hue_color(GLfloat& r, GLfloat& g, GLfloat& b, GLfloat h)
{
    GLfloat x = 3 * ( h * M_1_PI + 1 );
    int i = floor(x);
    GLfloat f = x-i, u = 1-f;
    switch( i )
    {
        case 0:  r = 1;  g = f;  b = 0;  break;
        case 1:  r = u;  g = 1;  b = 0;  break;
        case 2:  r = 0;  g = 1;  b = f;  break;
        case 3:  r = 0;  g = u;  b = 1;  break;
        case 4:  r = f;  g = 0;  b = 1;  break;
        case 5:  r = 1;  g = 0;  b = u;  break;
        default: r = 1;  g = 0;  b = 0;  break;
    }
}


/**
 set a RGB color as a function of a Hue value `h` in [-PI, PI],
 with a alpha-component equal to `a`.
 The colors follow in this order: red, green, blue, red ...
 */
gle_color gle::hue_color(GLfloat h, GLfloat a)
{
    GLfloat r, g, b;
    set_hue_color(r, g, b, h);
    return gle_color(r, g, b, a);
}



/**
 set a RGB color as a function of a value h in [0, 4].
 The result vary from dark-blue, blue, cyan, yellow, orange to red:
 - 0 : dark-blue
 - 1 : blue
 - 2 : green
 - 3 : red
 - 4 : dark red
 */
void gle::set_jet_color(GLfloat& r, GLfloat& g, GLfloat& b, const GLfloat h, const GLfloat hm)
{
    if ( h <= hm )
    {
        r = 0;
        g = 0;
        b = hm;
    }
    else
    {
        int i = (int)floor(h);
        GLfloat f = h-i;
        switch( i )
        {
            case 0:  r = 0;   g = 0;   b = f;   break;
            case 1:  r = 0;   g = f;   b = 1;   break;
            case 2:  r = f;   g = 1;   b = 1-f; break;
            case 3:  r = 1;   g = 1-f; b = 0;   break;
            case 4:  r = 1-f; g = 0;   b = 0;  if ( r < hm ) r = hm; break;
            default: r = hm;  g = 0;   b = 0;   break;
        }
    }
}


gle_color gle::jet_color(GLfloat h, GLfloat a)
{
    GLfloat r, g, b;
    set_jet_color(r, g, b, h, 0.3);
    return gle_color(r, g, b, a);
}


/**
 set a RGB color as a function of a value h in [0, 4].
 The result vary from black, blue, cyan, yellow, orange to red:
 - 0 : black
 - 1 : blue
 - 2 : green
 - 3 : red
 - 4 : yellow
 - 5 : white
 */
void gle::set_jet_colorE(GLfloat& r, GLfloat& g, GLfloat& b, const GLfloat h, const GLfloat hm)
{
    if ( h <= hm )
    {
        r = 0;
        g = 0;
        b = hm;
    }
    else
    {
        int i = (int)floor(h);
        GLfloat f = h-i;
        switch( i )
        {
            case 0:  r = 0;   g = 0;   b = f;   break;
            case 1:  r = 0;   g = f;   b = 1-f; break;
            case 2:  r = f;   g = 1-f; b = 0;   break;
            case 3:  r = 1;   g = f;   b = 0;   break;
            case 4:  r = 1;   g = 1;   b = f;   break;
            default: r = 1;   g = 1;   b = 1;   break;
        }
    }
}


gle_color gle::jet_colorE(GLfloat h, GLfloat a)
{
    GLfloat r, g, b;
    set_jet_colorE(r, g, b, h, 0.1);
    return gle_color(r, g, b, a);
}

