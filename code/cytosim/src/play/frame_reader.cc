// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "frame_reader.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"


// Use the second definition to get some verbose reports:
#define VLOG(ARG) ((void) 0)
//#define VLOG(ARG) std::clog << ARG;

//------------------------------------------------------------------------------

FrameReader::FrameReader() : input(0)
{
    clearPositions();
}

void FrameReader::rewind()
{
    input.rewind();
    clearPositions();
}

void FrameReader::openFile(std::string& file)
{
    clearPositions();
    
    int error = input.open(file.c_str(), "rb");
    
    if ( error )
    {
        //file was not found, we try 'gunzip'
        std::string tmp = file + ".gz";
        FILE* fp = fopen(tmp.c_str(), "r");
        if ( fp )
        {
            fclose(fp);
            tmp = "gunzip " + tmp;
            std::clog << tmp << std::endl;
            
            if ( 0 == system(tmp.c_str()) )
                input.open(file.c_str(), "rb");
        }
    }
    
    if ( !input.file() )
        throw InvalidIO("file `"+file+"' not found");

    if ( input.error() )
        throw InvalidIO("file `"+file+"' is invalid");
 
    input.vectorSize(DIM);
    //std::clog << "FrameReader: has openned " << obj_file << std::endl;
}


int FrameReader::badFile()
{
    if ( 0 == input.file() )
        return 8;
    
    if ( input.eof() )
        input.clearerr();
    
    if ( ! input.good() )
        return 7;
    
    return 0;
}


void FrameReader::checkFile()
{
    if ( 0 == input.file() )
        throw InvalidIO("No open file");
    
    if ( input.eof() )
        input.clearerr();
    
    if ( ! input.good() )
        throw InvalidIO("File has errors");
}

//------------------------------------------------------------------------------
#pragma mark -

void FrameReader::clearPositions()
{
    VLOG("FrameReader: clear\n");
    
    frameIndex = -1;
    framePos.clear();
    framePos.reserve(1024);
}


void FrameReader::savePos(int frm, const fpos_t& pos, int s)
{
    if ( frm < 0 )
        return;
    
    unsigned inx = frm;
    
    if ( inx >= framePos.capacity() )
    {
        const unsigned chunk = 1024;
        unsigned sz = ( inx + chunk - 1 ) & ~( chunk -1 );
        framePos.reserve(sz);
    }
    
    if ( inx >= framePos.size() )
    {
        unsigned int i = framePos.size();
        framePos.resize(inx+1);
        while ( i <= inx )
            framePos[i++].status = 0;
    }
    
    if ( framePos[inx].status < s )
    {
        framePos[inx].status = s;
        framePos[inx].value = pos;
    
        VLOG("FrameReader: learned position of frame " << frm << '\n');
    }
}


/**
 This uses the current knowledge to move to a position
 in the file where we should find frame `frm`.
*/
int FrameReader::seekPos(int frm)
{
    if ( input.eof() )
        input.clearerr();
    
    if ( frm < 1 || framePos.size() < 1 )
    {
        VLOG("FrameReader: seekPos rewind\n");
        input.rewind();
        return 0;
    }
    
    unsigned int inx = frm;
    
    if ( inx >= framePos.size() ) 
        inx = framePos.size() - 1;

    while ( inx > 0  &&  framePos[inx].status == 0 )
        --inx;
    
    //check if we know already were the frame starts:
    if ( 0 < inx )
    {
        VLOG("FrameReader: using known position of frame " << inx << '\n');
        input.set_pos(framePos[inx].value);
        return inx;
    }
    else {
        VLOG("FrameReader: rewind\n");
        input.rewind();
        return 0;
    }
}


int FrameReader::lastFrame() const
{
    int res = framePos.size() - 1;
    while ( 0 < res  &&  framePos[res].status < 2 )
        --res;
    return res;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 scan file forward from current position to find the next occurence of FRAME_TAG
 @return 0 if the frame was found
*/
int FrameReader::seekFrame(const int frm)
{        
    VLOG("FrameReader::seekFrame("<< frm <<")\n");
    
    int inx = seekPos(frm);
    
    if ( inx == frm )
        return 0;
    
    fpos_t pos;
    std::string line;

    while ( ! input.eof() )
    {
        do {
            
            input.get_pos(pos);
            input.get_line(line);
            
            if ( input.eof() )
                return 1;
            
#ifdef BACKWARD_COMPATIBILITY
            if ( 0 == line.compare(0, 7, "#frame ") )
                break;
#endif
            
        } while ( line[0] != '#'
                 || line.compare(1, strlen(FRAME_TAG), FRAME_TAG) );

        VLOG("FrameReader: " << line << '\n');

        if ( ! input.eof() )
        {
            savePos(inx, pos, 2);
            if ( inx == frm )
            {
                input.set_pos(pos);
                return 0;
            }
            ++inx;
        }
    }
    
    VLOG("FrameReader::seekFrame("<< frm <<") encountered EOF\n");
    return 1;
}

//------------------------------------------------------------------------------
/** 
 returns 0 for success, an error code, or throws an exception
 */
int FrameReader::loadFrame(Simul& sim, int frm, const bool reload)
{
    if ( badFile() )
        return 7;

    VLOG("FrameReader::loadFrame("<<frm<<", " << reload <<")\n");
    
    // a negative index is counted from the last known frame
    if ( frm < 0 )
    {
        VLOG("FrameReader: going down from frame " << lastFrame() << '\n');
        frm += 1 + lastFrame();
        if ( frm < 0 )
            frm = 0;
    }
    
    // what we are looking for might already be in the buffer:
    if ( frm == frameIndex && ! reload )
        return 0;
    
    // it might be the next one in the buffer:
    if ( frm == 1+frameIndex )
        return loadNextFrame(sim);

    // try to find the start tag from there:
    
    if ( 0 != seekFrame(frm) )
        return 1;
    
    VLOG("FrameReader: reading frame "<< frm << '\n');
    
    fpos_t pos;
    input.get_pos(pos);

    // switch to cytosim at this point:
    if ( 0 == sim.reloadObjects(input) )
    {
        VLOG("FrameReader::loadFrame("<< frm <<") successful\n");
        frameIndex = frm;
        savePos(frameIndex, pos, 3);

        // the next frame should start at the current position:
        input.get_pos(pos);
        savePos(frameIndex+1, pos, 1);
        return 0;
    }
    else
    {
        VLOG("FrameReader::loadFrame("<< frm <<") EOF at frame " << frm << '\n');
        return 1;
    }
}


/** 
 returns 0 for success, an error code, or throws an exception
 */
int FrameReader::loadNextFrame(Simul& sim)
{
    if ( badFile() )
        return 7;
    
    fpos_t pos;
    input.get_pos(pos);

    if ( 0 == sim.reloadObjects(input) )
    {
        ++frameIndex;
        
        // the position we used was good, to read this frame
        savePos(frameIndex, pos, 3);

        VLOG("FrameReader::loadNextFrame() read frame " << currFrame() << '\n');
        
        // the next frame should start from the current position:
        input.get_pos(pos);
        savePos(frameIndex+1, pos, 1);
        return 0;
    } 
    else
    {
        VLOG("FrameReader::loadNextFrame() EOF after frame " << currFrame() << '\n');
        return 1;
    }
}

