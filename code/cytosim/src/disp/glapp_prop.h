// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef GLAPP_PROP_H
#define GLAPP_PROP_H

#include "real.h"
#include "gle_color.h"
#include "vector3.h"
#include "quaternion.h"
#include "property.h"

/// Parameter set for glApp
class glAppProp : public Property
{    
public:
    
    /**
     @defgroup glAppPar Display parameters: Graphics
     @ingroup DisplayParameters
     @{
     */
    
    
    /// flag to use a double buffer for smoother rendering (default=1)
    /**
    http://en.wikipedia.org/wiki/Multiple_buffering#Double_buffering_in_computer_graphics
    */
    int          buffered;
    
    /// flag to switch to full-screen mode
    int          full_screen;
    
    /// flag to display information on screen
    int          show_message;
    
    /// text displayed in center of window
    std::string  message;
    
    /// @}
    
public:
    
    /// constructor
    glAppProp(const std::string& n) : Property(n)  { clear(); }
    
    /// destructor
    ~glAppProp()  { }
    
    /// identifies the property
    std::string category() const { return "simul:display"; }
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// return a carbon copy of object
    Property* clone() const { return new glAppProp(*this); }

    /// write all values
    void write_values(std::ostream &) const;
    
};

#endif
