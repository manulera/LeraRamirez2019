// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "saveimage.h"
#include <cstdlib>
#include <cstring>

#ifndef NO_OPENGL
#  include "opengl.h"
#endif

/// output stream
FILE * SaveImage::err = stderr;


bool SaveImage::supported(const char format[])
{
#ifdef HAS_PNG
    if ( 0 == strcasecmp(format, "png") )
        return true;
#endif
    if ( 0 == strcasecmp(format, "ppm") )
        return true;
    
    return false;
}


//------------------------------------------------------------------------------
#pragma mark Images

/**
 saveImage(...) will read pixels from the current OpenGL read buffer,
 and save them in a file with the requested format
 */
int SaveImage::saveImage(const char * filename,
                         const char format[],
                         int downsample)
{
#ifndef NO_OPENGL
    GLint vp[4];
    glGetIntegerv(GL_VIEWPORT, vp);
    //printf("saveImage viewport %i %i %i %i\n", vp[0], vp[1], vp[2], vp[3]);
    
    return saveImage(filename, format, vp[0], vp[1], vp[2], vp[3], downsample);
#endif
    return 1;
}


/**
 saveImage(...) will read pixels from the current OpenGL read buffer,
 and save them in a file with the requested format
 */
int SaveImage::saveImage(const char * filename,
                         const char format[],
                         const int xpos, const int ypos,
                         const int width, const int height,
                         int downsample)
{
    int res = 1;

#ifndef NO_OPENGL
    //allocate memory to hold image:
    GLubyte* pixels = new GLubyte[width*height*3];

    //set the alignment to double-words
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    
#if ( 0 )
    GLint readbuf = 0, drawbuf = 0;
    glGetIntegerv(GL_READ_BUFFER, &readbuf);
    glGetIntegerv(GL_DRAW_BUFFER, &drawbuf);
    printf("framebuffers: read %i draw %i\n", readbuf, drawbuf);
#endif
    
    //read the pixel values, from top-left corner:
    glReadPixels(xpos, ypos, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
    
    GLenum glError = glGetError();
    
    if ( glError != GL_NO_ERROR )
    {
        if ( err )
            fprintf(err, "Error: SaveImage could not read pixels (error %i)\n", glError);
        delete[] pixels;
        return glError;
    }
    
    res = savePixels(filename, format, pixels, width, height, downsample);
    
    delete[] pixels;

#else
    if ( err )
        fprintf(err, "Error: SaveImage was built without OpenGL!\n");
#endif

    return res;
}



/**
 After setting a higher resolution, this will translate the ModelView to produce several
 images that will be stiched together in memory, into an image with higher resolution.
 This works even if the image is larger than the maximum OpenGL viewPort,
 but there can be artifacts due to stitching imperfections.
 */
int SaveImage::saveCompositeImage(const int mag,
                                  const char * filename,
                                  const char format[],
                                  const int width, const int height,
                                  const double pixel_size,
                                  void (*display)(int, void *), void * arg,
                                  int downsample)
{
    if ( ! supported(format) )
        return -1;
    
    int res = 1;

#ifndef NO_OPENGL

    int mw = mag * width;
    int mh = mag * height;
    
    const int pix = 3;  //number of bytes for each pixel
    GLubyte* pixels = new GLubyte[mw*mh*pix];
    GLubyte* sub    = new GLubyte[width*height*pix];
    
    const double cc = ( mag - 1 ) * 0.5;
    const double dx = width*pixel_size / mag;
    const double dy = height*pixel_size / mag;
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glScalef(mag, mag, mag);

    for ( int iy = 0; iy < mag ; ++iy )
    {
        for ( int ix = 0; ix < mag; ++ix )
        {
            glPushMatrix();
            glTranslated((cc-ix)*dx, (cc-iy)*dy, 0);
            display(mag, arg);
            glPopMatrix();
            
            glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, sub);

            GLubyte * dst = &pixels[width*pix*(ix+mh*iy)];
            for ( int ii=0; ii<height; ++ii )
                memcpy(&dst[ii*mw*pix], &sub[ii*width*pix], width*pix);
        }
    }
    glPopMatrix();
    delete[] sub;

    res = savePixels(filename, format, pixels, mw, mh, downsample);
    
    delete[] pixels;
    
#endif
    return res;
}


/**
 This adjusts the Viewport to produce an image with higher resolution.
 The result should be better than saveCompositeImage, but uses more
 memory on the graphic card.
 
 */
int SaveImage::saveMagnifiedImage(const int mag,
                                  const char * filename,
                                  const char format[],
                                  const int width, const int height,
                                  void (*display)(int, void *), void * arg,
                                  int downsample)
{
    if ( ! supported(format) )
        return -1;
    
    int res = 1;
    
#ifndef NO_OPENGL
    GLint vp[4];
    glGetIntegerv(GL_VIEWPORT, vp);

    int W = mag * width;
    int H = mag * height;

    GLint maxdim[2];
    glGetIntegerv(GL_MAX_VIEWPORT_DIMS, maxdim);
    if ( W > maxdim[0]  ||  H > maxdim[1] )
    {
        if ( err )
            fprintf(err, "SaveImage:: exceeding maximum supported dimensions (%ix%i)\n", (int)maxdim[0], (int)maxdim[1]);
        return 1;
    }
    
    const int pix = 3;  //number of bytes for each pixel
    //allocate memory to hold the full image:
    GLubyte* pixels = new GLubyte[W*H*pix];
    GLubyte* sub    = new GLubyte[width*height*pix];

    //set the alignment ?
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    
    for ( int iy = 0; iy < mag ; ++iy )
    {
        for ( int ix = 0; ix < mag; ++ix )
        {
            glViewport(-ix*width, -iy*height, W, H);
            display(mag, arg);
            glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, sub);
            
            GLubyte * dst = &pixels[width*pix*(ix+H*iy)];
            for ( int ii=0; ii<height; ++ii )
                memcpy(&dst[ii*W*pix], &sub[ii*width*pix], width*pix);
        }
    }
    delete[] sub;
    
    res = savePixels(filename, format, pixels, W, H, downsample);

    delete[] pixels;
    //restore initial viewport:
    glViewport(vp[0], vp[1], vp[2], vp[3]);
#endif
    
    return res;
}


//------------------------------------------------------------------------------
#pragma mark Pixels

FILE * SaveImage::openFile(const char * filename)
{
    if ( filename[0] == '\0' )
        return 0;
    FILE * file = fopen(filename, "wb");
    if ( file && ferror(file) )
    {
        fclose(file);
        file = 0;
    }
    return file;
}


int SaveImage::savePixels(const char format[],
                          FILE * file,
                          GLubyte* pixels,
                          int width, int height)
{
    if ( 0 == strcasecmp(format, "ppm") )
        return saveColorPPM(file, pixels, width, height);
    
    if ( 0 == strcasecmp(format, "png") )
        return saveColorPNG(file, pixels, width, height);

    return -1;
}


int SaveImage::savePixels(const char * filename,
                          const char format[],
                          GLubyte* pixels,
                          int width, int height)
{
    FILE * file = openFile(filename);
    
    if ( file )
    {
        int res = savePixels(format, file, pixels, width, height);
        
        if ( err )
        {
            if ( res == 0 )
                fprintf(err, " saved %ix%i snapshot %s\n", width, height, filename);
            else
                fprintf(err, " error %i occured while saving %s\n", res, filename);
        }
        
        fclose(file);
        
        // remove file if an error occured:
        if ( res )
            remove(filename);
        
        return res;
    }
    
    return 3;
}


//------------------------------------------------------------------------------
#pragma mark - Downsample

/**
 This will downsample pixelmap `src` and set destination `dst`. The pixel
 array `src` should be of size `3*W*H` with 3 bytes per pixels: R, G, B,
 while `src` will be `bin*bin` times smaller. Pixels are stored in row order
 from the lowest to the highest row, left to right in each row (as in OpenGL).
 The pixels components of `src` are averaged to produce `dst`.
 */
void SaveImage::downsampleRGB(GLubyte* dst, GLubyte const* src,
                              unsigned W, unsigned H, unsigned bin)
{
    const unsigned s = bin * bin;
    const unsigned sx = W / bin;
    const unsigned sy = H / bin;
    
#if ( 0 )
    //reset destination:
    for ( unsigned u = 0; u < sx*sy; ++u )
    {
        dst[3*u  ] = 0xFF;
        dst[3*u+1] = 0xFF;
        dst[3*u+2] = 0xFF;
    }
#endif
    
    for ( unsigned x = 0; x < sx; ++x )
        for ( unsigned y = 0; y < sy; ++y )
        {
            unsigned r = 0, g = 0, b = 0;
            for ( unsigned dx = 0; dx < bin; ++dx )
                for ( unsigned dy = 0; dy < bin; ++dy )
                {
                    GLubyte const* p = src + 3 * ( dx + bin*x + W*(dy+bin*y) );
                    r += p[0];
                    g += p[1];
                    b += p[2];
                }
            
            dst[3*(x+sx*y)  ] = r / s;
            dst[3*(x+sx*y)+1] = g / s;
            dst[3*(x+sx*y)+2] = b / s;
        }
}


int SaveImage::savePixels(const char * filename,
                           const char format[],
                           GLubyte* pixels,
                           int width, int height,
                           int downsample)
{
    int res = 0;
    if ( downsample > 1 )
    {
        //printf("downsampling %i to : %i %i\n", downsample, mw, mh);
        int W = width / downsample;
        int H = height / downsample;
        GLubyte* img = new GLubyte[3*W*H];
        downsampleRGB(img, pixels, width, height, downsample);
        res = savePixels(filename, format, img, W, H);
        delete[] img;
    }
    else
    {
        res = savePixels(filename, format, pixels, width, height);
    }
    return res;
}


//------------------------------------------------------------------------------
//------------------------------- PPM FORMAT -----------------------------------
//------------------------------------------------------------------------------
#pragma mark - PPM

/**
 Write the image in the Portable Pixmap format, also Netpbm format (man -k ppm).
 We use here the 'raw' binary format starting with P6 
 */
int SaveImage::saveColorPPM(FILE* file,
                            const GLubyte pixels[],
                            const int width, const int height)
{
    if ( file==0 || ferror(file) )
        return 1;
    
    fprintf(file, "P6\n");
    fprintf(file, "%i %i\n", width, height);
    fprintf(file, "255\n");
    
    //write the pixels binary, line by line:
    for ( int ii = height-1; ii >= 0; --ii )
        fwrite(&pixels[3*ii*width], 1, 3*width, file);
    return 0;
}



//------------------------------------------------------------------------------
//-------------------------------- PNG FORMAT ----------------------------------
//------------------------------------------------------------------------------
#pragma mark - PNG

#ifndef HAS_PNG

int SaveImage::savePNG(FILE*, void* pixels,
                       const int, const int,
                       const int, const int)
{
    if ( err )
        fprintf(err, "PNG format not supported (recompilation needed)\n");
    return -1;
}

#else

#include <png.h>

/*
int SaveImage::readPNG(FILE* fp, png_bytep *& row_pointers, int& bit_depth, int& color_mode, int& width, int& height)
{
    if (!fp)
        return 9;
    
    char header[8];    // 8 is the maximum size that can be checked
    fread(header, 1, 8, fp);
    
    if (png_sig_cmp(header, 0, 8))
        return 8;
    
    // initialize stuff 
    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    
    if (!png_ptr)
        return 7;
    
    png_infop info_ptr = png_create_info_struct(png_ptr);
    
    if (!info_ptr)
        return 6;
    
    if (setjmp(png_jmpbuf(png_ptr)))
        return 5;
    
    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 8);
    
    png_read_info(png_ptr, info_ptr);
    
    width = png_get_image_width(png_ptr, info_ptr);
    height = png_get_image_height(png_ptr, info_ptr);
    
    color_type = png_get_color_type(png_ptr, info_ptr);
    bit_depth = png_get_bit_depth(png_ptr, info_ptr);
    
    int number_of_passes = png_set_interlace_handling(png_ptr);
    png_read_update_info(png_ptr, info_ptr);
    
    
    // read file
    if (setjmp(png_jmpbuf(png_ptr)))
        return 4;
    
    if ( row_pointers )
        free(row_pointers);
    
    row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    
    for (y=0; y<height; y++)
        row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr,info_ptr));
    
    png_read_image(png_ptr, row_pointers);
}
*/


int savePNG(FILE* file,
            png_bytep row_pointers[],
            const int bit_depth, const int nb_colors,
            const int width, const int height)
{
    if ( !file )
        return 9;
    
    if ( bit_depth != 8 && bit_depth != 16 )
        return 9;
    
    int color_type = -1;
    
    if ( nb_colors == 1 )
        color_type = PNG_COLOR_TYPE_GRAY;
    
    if ( nb_colors == 3 )
        color_type = PNG_COLOR_TYPE_RGB;
    
    if ( nb_colors == 4 )
        color_type = PNG_COLOR_TYPE_RGBA;
    
    if ( color_type < 0 )
        return 8;
    
    /* initialize stuff */
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    
    if (!png_ptr)
        return 7;
    
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
        return 6;
    
    if (setjmp(png_jmpbuf(png_ptr)))
        return 5;
    
    png_init_io(png_ptr, file);
    
    /* write header */
    if (setjmp(png_jmpbuf(png_ptr)))
        return 4;
    
    png_set_IHDR(png_ptr, info_ptr, width, height,
                 bit_depth, color_type,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    
    png_write_info(png_ptr, info_ptr);
    
    /* write bytes */
    if (setjmp(png_jmpbuf(png_ptr)))
        return 3;
    
    if ( nb_colors == 1 && bit_depth == 16 )
        png_set_swap(png_ptr);
    
    png_write_image(png_ptr, row_pointers);
    
    
    /* end write */
    if (setjmp(png_jmpbuf(png_ptr)))
        return 2;
    
    png_write_end(png_ptr, NULL);
 
    return 0;
}


int SaveImage::savePNG(FILE* file, void* pixels,
                       const int bit_depth, const int nb_colors,
                       const int width, const int height)
{
    png_bytep * rows = new png_bytep[height];
    
    int bytes_per_row = ( bit_depth / 8 ) * nb_colors * width;
    
    png_byte * start = (png_byte*)pixels;
    
    for ( int y = 0; y < height; ++y )
        rows[y] = start + bytes_per_row * ( height-y-1 );
    
    int res = savePNG(file, rows, bit_depth, nb_colors, width, height);
    
    delete[] rows;

    return res;
}

#endif


/**
 Save RGBA image, 4 x 8-bits per pixel
 */
int SaveImage::saveAlphaPNG(FILE* file,  void* pixels,
                           const int width, const int height)
{
    return savePNG(file, pixels, 8, 4, width, height);
}

/**
 Save RGB image, 3 x 8-bits per pixel
 */
int SaveImage::saveColorPNG(FILE* file,  void* pixels,
                            const int width, const int height)
{
    return savePNG(file, pixels, 8, 3, width, height);
}

/**
 Save 16-bits gray-level image
 */
int SaveImage::saveGrayPNG(FILE* file,  void* pixels,
                           const int width, const int height)
{    
    return savePNG(file, pixels, 16, 1, width, height);
}

