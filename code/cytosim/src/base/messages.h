// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef  MESSAGES_H
#define  MESSAGES_H

#include <cstdio>
#include <cstdarg>
#include <string>
#include "filewrapper.h"


/// Control output file with different levels of verbosity
/**
 @todo build Messages from a OutputStream
 */
class Messages : public FileWrapper
{
    
    ///max. number of output-warnings
    static const int max_warnings = 50;
      
private:
    
    ///verbose level, PRINT has a option specifying a level which is compared to mVerbose
    int verbose_;
    
    /// number of warnings already given
    int num_warnings_;

public:
        
    ///Constructor
    Messages();
    
    ///Destructor 
    virtual ~Messages() {}
    
    /// suppress all output by setting Verbose to -1
    void    silent()               { verbose_ = -1; }

    ///suppresses most output by setting Verbose to 0
    void    quiet()                { verbose_ = 0; }

    /// return the verbose level
    int     verboseLevel()         { return verbose_; }
    
    /// set verbose to level m
    void    setVerbose(int m)      { verbose_ = m; }
        
    ///convenient access to print() with the () operator
    void    operator()(const char* fmt, ...);

    ///convenient access to print() with the () operator, with a verbose level
    void    operator()(int, const char* fmt, ...);
    
    ///warning() is equivalent to print() with "Warning:" in front
    void    warning(const char* fmt, ...);

    ///warning() is equivalent to print() with "Warning:" in front
    void    warning(std::string const&);
};

/// global instantiation used for Cytosim output
/** C++ makes no guaranty about the order in which global variables are initialized.
Therefore an error may occur if another constructor uses MSG before it is constructed */
extern Messages MSG;

#endif
