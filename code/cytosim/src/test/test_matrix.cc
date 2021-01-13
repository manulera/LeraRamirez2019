// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <sys/time.h>
#include "dim.h"
#include "exceptions.h"
#include "random.h"
#include "tictoc.h"
#include "vecprint.h"

#include "matsparsesym.h"
#include "matsparsesym1.h"
#include "matsparsesym2.h"
#include "matsparsesymblk.h"
#include "matsparsesym4.h"
#include "matsparsesym5.h"

typedef MatrixSparseSymmetricBlock MatrixSparseSymmetric3;

extern Random RNG;

real alpha = 1.0;

void zero(unsigned s, real * vec)
{
    for ( int n=0; n<s; ++n )
        vec[n] = 0;
}

//------------------------------------------------------------------------------

void compare(const int size,  Matrix & mat1, Matrix& mat2, unsigned fill)
{
    real * tmp1 = new real[size*size];
    real * tmp2 = new real[size*size];
    
    mat1.makeZero();
    mat2.makeZero();
    
    mat1.resize(size);
    mat2.resize(size);
    
    for ( int n=0; n<fill; ++n )
    {
        real a = RNG.preal();
        unsigned ii = RNG.pint(size);
        unsigned jj = RNG.pint(size);
        mat1(ii, jj) += a;
        mat2(ii, jj) += a;
    }
    
    unsigned inx = DIM * RNG.pint(size/DIM-1);
    unsigned nbc = DIM * RNG.pint((size-inx)/DIM) + DIM;
    std::clog << "size " << size << " inx " << inx << " nbc " << nbc << "\n";
    
    zero(size*size, tmp1);
    
    mat1.addDiagonalBlock(tmp1, size, inx, nbc);
    
    //std::clog<<"mat1:\n";
    //VecPrint::matPrint(std::clog, nbc, nbc, tmp1, size);

    zero(size*size, tmp2);
    mat2.addDiagonalBlock(tmp2, size, inx, nbc);
    
    //std::clog<<"mat2:\n";
    //VecPrint::matPrint(std::clog, nbc, nbc, tmp2, size);
    
    for ( int n=0; n<size*size; ++n )
        tmp1[n] -= tmp2[n];

    std::clog<<"diff:\n";
    VecPrint::matPrint(std::clog, nbc, nbc, tmp1, size);
    
    delete[] tmp1;
        delete[] tmp2;
}




void testMatrix(const int RUNCNT, Matrix & mat,
                const int size, real * x, real * y,
                const int fill, int inx[], int iny[])
{
    zero(DIM*size, y);
    mat.resize(DIM*size);

    TicToc::tic();
    for ( int ii=0; ii<RUNCNT; ++ii )
    {
        mat.makeZero();
        for ( int kk=0; kk<fill; ++kk ) 
        {
#if ( DIM == 3 )
            mat(DIM*inx[kk]  , DIM*inx[kk]  ) += alpha;
            mat(DIM*inx[kk]+1, DIM*inx[kk]  ) += alpha;
            mat(DIM*inx[kk]+2, DIM*inx[kk]  ) += alpha;
            mat(DIM*iny[kk]  , DIM*inx[kk]  ) += alpha;
            mat(DIM*iny[kk]+1, DIM*inx[kk]  ) += alpha;
            mat(DIM*iny[kk]+2, DIM*inx[kk]  ) += alpha;

            mat(DIM*inx[kk]+1, DIM*inx[kk]+1) -= alpha;
            mat(DIM*inx[kk]+2, DIM*inx[kk]+1) -= alpha;
            mat(DIM*iny[kk]  , DIM*inx[kk]+1) -= alpha;
            mat(DIM*iny[kk]+1, DIM*inx[kk]+1) -= alpha;
            mat(DIM*iny[kk]+2, DIM*inx[kk]+1) -= alpha;
 
            mat(DIM*inx[kk]+2, DIM*inx[kk]+2) -= alpha;
            mat(DIM*iny[kk]  , DIM*inx[kk]+2) -= alpha;
            mat(DIM*iny[kk]+1, DIM*inx[kk]+2) -= alpha;
            mat(DIM*iny[kk]+2, DIM*inx[kk]+2) -= alpha;
 
            mat(DIM*iny[kk]  , DIM*iny[kk]  ) -= alpha;
            mat(DIM*iny[kk]+1, DIM*iny[kk]  ) -= alpha;
            mat(DIM*iny[kk]+2, DIM*iny[kk]  ) -= alpha;
 
            mat(DIM*iny[kk]+1, DIM*iny[kk]+1) -= alpha;
            mat(DIM*iny[kk]+2, DIM*iny[kk]+1) -= alpha;
 
            mat(DIM*iny[kk]+2, DIM*iny[kk]+2) -= alpha;
#elif ( DIM == 2 )
            mat(DIM*inx[kk]  , DIM*inx[kk]  ) += alpha;
            mat(DIM*inx[kk]+1, DIM*inx[kk]  ) += alpha;
            mat(DIM*iny[kk]  , DIM*inx[kk]  ) += alpha;
            mat(DIM*iny[kk]+1, DIM*inx[kk]  ) += alpha;
            
            mat(DIM*inx[kk]+1, DIM*inx[kk]+1) -= alpha;
            mat(DIM*iny[kk]  , DIM*inx[kk]+1) -= alpha;
            mat(DIM*iny[kk]+1, DIM*inx[kk]+1) -= alpha;
            
            mat(DIM*iny[kk]  , DIM*iny[kk]  ) -= alpha;
            mat(DIM*iny[kk]+1, DIM*iny[kk]  ) -= alpha;
            
            mat(DIM*iny[kk]+1, DIM*iny[kk]+1) -= alpha;
#else
            mat(inx[kk], inx[kk]) += alpha;
            mat(iny[kk], inx[kk]) += alpha;
            mat(iny[kk], iny[kk]) -= alpha;
#endif
        }
    }
    double ts = TicToc::toc();
    
    if ( mat.size() < 1 )
    {
        std::cout << "mat:" << mat.what() << "\n";
        mat.printSparse(std::cout);
    }

    
    TicToc::tic();
    for ( int ii=0; ii<RUNCNT; ++ii )
    {
        mat.prepareForMultiply();
        for ( int ii=0; ii<16; ++ii )
            mat.vecMulAdd(x, y);
    }
    double t1 = TicToc::toc();
    
    real sum = 0;
    for ( int kk=0; kk<size; ++kk )
        sum += y[kk];

    printf("Matrix %20s : ", mat.what().c_str());
    printf("set %6.0f  mul %6.0f", ts, t1);

    mat.resize(size);
    
    TicToc::tic();
    for ( int ii=0; ii<RUNCNT; ++ii )
    {
        mat.makeZero();
        for ( int kk=0; kk<fill; ++kk )
        {
            mat(inx[kk], inx[kk]) += alpha;
            mat(iny[kk], iny[kk]) += alpha;
            mat(iny[kk], inx[kk]) -= alpha;
        }
    }
    double t2 = TicToc::toc();
    
    TicToc::tic();
    mat.prepareForMultiply();
#if ( DIM == 3 )
    for ( int ii=0; ii<32*RUNCNT; ++ii )
       mat.vecMulAddIso3D(x, y);
#elif ( DIM == 2 )
    for ( int ii=0; ii<32*RUNCNT; ++ii )
        mat.vecMulAddIso2D(x, y);
#endif
    double t3 = TicToc::toc();
    
    printf(" isoset %6.0f isomul %6.0f  tot %6.0f", t2, t3, ts+t1+t2+t3);
    printf("  res %.20f\n", sum);
}



void testMatrix(const int size, const int fill)
{
    printf("\n **** Matrix size %i  filled %.1f %% :\n", size, 100*fill/size/double(size));
    MatrixSparseSymmetric  mat0;
    MatrixSparseSymmetric1 mat1;
    MatrixSparseSymmetric2 mat2;
    MatrixSparseSymmetric3 mat3;
    MatrixSparseSymmetric4 mat4;
    MatrixSparseSymmetric5 mat5;
    
    int * inx = new int[fill];
    int * iny = new int[fill];

    for ( int kk=0; kk<fill; ++kk )
    {
        inx[kk] = RNG.pint(size);
        iny[kk] = RNG.pint(size);
    }
    
    real * x = new real[DIM*size];
    real * y = new real[DIM*size];
    
    for ( int ii=0; ii<DIM*size; ++ii )
        x[ii] = RNG.sreal();
    alpha = RNG.sreal();
   
    unsigned run = 1;
    testMatrix(run, mat0, size, x, y, fill, inx, iny);
    testMatrix(run, mat1, size, x, y, fill, inx, iny);
    testMatrix(run, mat2, size, x, y, fill, inx, iny);
    testMatrix(run, mat3, size, x, y, fill, inx, iny);
    testMatrix(run, mat4, size, x, y, fill, inx, iny);
    //testMatrix(run, mat5, size, x, y, fill, inx, iny);
    
    delete[] x;
    delete[] y;
    delete[] inx;
    delete[] iny;
}



 int main( int argc, char* argv[] )
{
    RNG.seedTimer();

    if ( 0 )
    {
        MatrixSparseSymmetric1 mat1;
        MatrixSparseSymmetric3 mat3;
        compare(12, mat1, mat3, 1);
        testMatrix(6, 1);
        testMatrix(DIM*7, 1);
        testMatrix(DIM*17, 2);
    }
    else
    {
        testMatrix(DIM*13, 234);
        testMatrix(DIM*33, 1111);
        testMatrix(DIM*197, 1<<17);
        testMatrix(DIM*437, 1<<19);
        testMatrix(DIM*1111, 1<<21);
    }
    
    return EXIT_SUCCESS;
}


