// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "play_prop.h"
#include "glossary.h"
#include "saveimage.h"
#include "smath.h"

//------------------------------------------------------------------------------
void PlayProp::clear()
{
    play         = 0;
    loop         = 0;
    exit_at_eof  = false;
    period       = 1;
    delay        = 32;
    live         = 0;
    
    report_index = 0;
    report       = "";
    report0      = "";
    report1      = "fiber:lengths";
    report2      = "fiber:dynamics";
    report3      = "single";
    report4      = "couple";
    report5      = "inventory";
    report6      = "space";
    report7      = "fiber:age";
    report7      = "fiber:segment";
    report9      = "fiber:distribution";

    for ( int k = 0; k < NB_MAGIC_KEYS; ++k )
    {
        magic_key[k]  = 0;
        magic_code[k] = "";
    }
    
    save_images = false;
    if ( SaveImage::supported("png") )
        image_format = "png";
    else
        image_format = "ppm";

    image_dir    = "";
    downsample   = 1;
    image_index  = 0;
    poster_index = 0;
}

//------------------------------------------------------------------------------
void PlayProp::read(Glossary& glos)
{
    glos.set(play,         "play");
    glos.set(loop,         "loop");
    glos.set(period,       "period");
    if ( period < 1 ) period = 1;
    glos.set(delay,        "delay");
    glos.set(save_images,  "save_images");
    glos.set(image_format, "image_format");
    glos.set(image_dir,    "image_dir");
    glos.set(downsample,   "downsample");
    glos.set(downsample,   "downsampling");
    
    if ( ! SaveImage::supported(image_format.c_str()) )
        throw InvalidParameter("unsupported image format");
    
    std::string var = "magic_key";
    for ( int k = 0; k < NB_MAGIC_KEYS; ++k )
    {
        glos.set(magic_key[k],  var);
        glos.set(magic_code[k], var, 1);
        var = "magic_key" + sMath::repr(k);
    }
    
    glos.set(report,      "report");
    glos.set(report0,     "report0");
    glos.set(report1,     "report1");
    glos.set(report2,     "report2");
    glos.set(report3,     "report3");
    glos.set(report4,     "report4");
    glos.set(report5,     "report5");
    glos.set(report6,     "report6");
    glos.set(report7,     "report7");
    glos.set(report8,     "report8");
    glos.set(report9,     "report9");
}


//------------------------------------------------------------------------------

void PlayProp::write_values(std::ostream & os) const
{
    write_value(os, "play",   play);
    write_value(os, "loop",   loop);
    write_value(os, "period", period);
    write_value(os, "delay",  delay);
    write_value(os, "report", report);
    write_value(os, "save_images", save_images);
    write_value(os, "image_format", image_format);
    write_value(os, "image_dir", image_dir);
    write_value(os, "downsample", downsample);

    for ( int k = 0; k < NB_MAGIC_KEYS; ++k )
    {
        std::string var = "magic_key" + sMath::repr(k);
        write_value(os, var, magic_key[k], "("+magic_code[k]+")");
    }
}

//------------------------------------------------------------------------------

void PlayProp::toggleReport(bool alt)
{
    if ( alt )
    {
        if ( report_index < 5 )
            report_index = 5;
        else
            report_index = 5 + ( report_index - 4 ) % 5;
    }
    else
        report_index = ( report_index + 1 ) % 5;

    switch( report_index )
    {
        case 0: report = report0; break;
        case 1: report = report1; break;
        case 2: report = report2; break;
        case 3: report = report3; break;
        case 4: report = report4; break;
        case 5: report = report5; break;
        case 6: report = report6; break;
        case 7: report = report7; break;
        case 8: report = report8; break;
        case 9: report = report9; break;
    }
}

