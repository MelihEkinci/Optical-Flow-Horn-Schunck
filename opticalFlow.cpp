# include "iostream"
#include <memory>
#include <stdexcept>
#include <algorithm>
#include "opticalFlow.hpp"
#include <vector>
#include <omp.h>
#include <thread>
#include <cstdint>
#include <iomanip>
#include <string.h>

# define cimg_display 0 
# define THREADS 8
# include "CImg-v.3.1.6/CImg.h"
using namespace cimg_library;

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int32_t  
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Core"  
#include "eigen-3.4.0/Eigen/Sparse" 
#include "eigen-3.4.0/Eigen/IterativeLinearSolvers"
#include "eigen-3.4.0/Eigen/SparseCholesky"
#include "eigen-3.4.0/Eigen/SparseLU"
#include "eigen-3.4.0/Eigen/SparseQR"

//Global variable for threading
std::atomic<int> count_a{0};

//preparing linear system of equations A,b
std::tuple<Eigen::SparseMatrix<float>,Eigen::SparseMatrix<float>> fill_A_b(CImg <float> Ix, CImg <float> Iy, CImg <float> It, int a) {
    double height=Ix.height();
    double width=Ix.width();
    double num_pixel=height*width;
	int count_a = 0;

    std::cout<<"initializing A and b"<<std::endl;
    Eigen::SparseMatrix<float> matA(2*num_pixel , 2*num_pixel);
    Eigen::SparseMatrix<float> matb(2*num_pixel ,1);
    matA.setZero();
    matb.setZero();

    std::cout<<"Filling A and b with coefficients"<<std::endl;
    /*
    #pragma omp parallel for simd  \
                collapse(2)\
                schedule(simd:static,THREADS)
    */



    for (int i=0; i!=Ix.height() ;++i){
        for (int j=0; j!=Ix.width() ;++j){
            //TODO CHECK WHETHER MATRIX IS FILLED CORRECTLY
            //first equation
            matA.coeffRef(count_a,i*width+j)=Ix(j,i)*Ix(j,i)+4*a;

            if (i*width+j+1>=0 && i*width+j+1<2*num_pixel)
                matA.coeffRef(count_a,i*width+j+1)=-a;

            if (i*width+j-1>=0 && i*width+j-1<2*num_pixel)
                matA.coeffRef(count_a,i*width+j-1)=-a;

            if ((i+1)*width+j>=0 && (i+1)*width+j<2*num_pixel)
                matA.coeffRef(count_a,(i+1)*width+j)=-a;

            if ((i-1)*width+j>=0 && (i-1)*width+j<2*num_pixel)
                matA.coeffRef(count_a,(i-1)*width+j)=-a;
            
            matA.coeffRef(count_a,num_pixel+i*width+j)=Ix(j,i)*Iy(j,i);

            //second equation
            matA.coeffRef(num_pixel+count_a,num_pixel+i*width+j)=Iy(j,i)*Iy(j,i)+4*a;

            if (num_pixel+i*width+j+1>=0 && num_pixel+i*width+j+1<2*num_pixel)
                matA.coeffRef(num_pixel+count_a,num_pixel+i*width+j+1)=-a;

            if (num_pixel+i*width+j-1>=0 && num_pixel+i*width+j-1<2*num_pixel)
                matA.coeffRef(num_pixel+count_a,num_pixel+i*width+j-1)=-a;

            if (num_pixel+(i+1)*width+j>=0 && num_pixel+(i+1)*width+j<2*num_pixel)
                matA.coeffRef(num_pixel+count_a,num_pixel+(i+1)*width+j)=-a;

            if (num_pixel+(i-1)*width+j>=0 && num_pixel+(i-1)*width+j<2*num_pixel)
                matA.coeffRef(num_pixel+count_a,num_pixel+(i-1)*width+j)=-a;

            matA.coeffRef(num_pixel+count_a,(i)*width+j)=Ix(j,i)*Iy(j,i);

            //filling_b
            matb.coeffRef(count_a,0)=-Ix(j,i)*It(j,i);
            matb.coeffRef(num_pixel+count_a,0)=-Iy(j,i)*It(j,i);

            ++count_a;
            //std::cout<<j<<'-'<<i<<std::endl;
            
            }
    }
    std::cout<<"A and b is ready"<<std::endl;

    return std::make_tuple(matA,matb);

}


Eigen::SparseMatrix<float> linear_solve(Eigen::SparseMatrix<float> matA, Eigen::SparseMatrix<float> matb) {
    std::cout << "solving started" << std::endl;
    Eigen::SparseLU <Eigen::SparseMatrix<float>> solver;
    solver.compute(matA);
    std::cout << "compute is done" << std::endl;
    Eigen::SparseMatrix<float> x = solver.solve(matb);
    std::cout << "solved" << std::endl;
    return x;
}

std::tuple<CImg < float >,CImg < float >> calculate_u_v (Eigen::SparseMatrix<float> x, int height, int width){

    CImg <float> u(width,height);
    CImg <float> v(width,height);

    #pragma omp parallel for simd  \
            collapse(2)\
            schedule(simd:static)
    for (int i=0 ; i != u.height(); ++i){
        for(int j=0 ; j != u.width(); ++j){
            u(j,i)=x.coeff(i*u.width()+j,0);
            v(j,i)=x.coeff(u.width()*u.height()+i*u.width()+j,0);
            
        }
    }
    return std::make_tuple(u,v);
}

int main(int argc, char* argv[]){
    std::cout<<argv[1]<<std::endl;
    std::cout<<argv[2]<<std::endl;
    std::cout<<argv[3]<<std::endl;
    std::cout<<argv[4]<<std::endl;

    // create image objects
    CImg < float > in_image1(argv[1]);
    CImg < float > in_image2(argv[2]);
    std::cout << "Step 1 : Read Images" << std::endl;

    /*
    for (int i=0 ; i != in_image1.height(); ++i){
        for(int j=0 ; j != in_image1.width(); ++j){
            std::cout << "pixel value: "<<in_image1(j,i) << std::endl;
        }
    }
    */

    //stacking two images to a single 3d object
    CImg <float> img_3d(in_image1.width(),in_image1.height(),2);
    int count=0;
      #pragma omp parallel for simd  \
                collapse(2)\
                schedule(simd:static)
      for (int i=0 ; i != in_image1.height(); ++i){
        for(int j=0 ; j != in_image1.width(); ++j){
            count+=1;
            img_3d(j,i,0)=in_image1(j,i);
            img_3d(j,i,1)=in_image2(j,i);
            
        }
    }
    //rescaling pixels between 0 and 1
    in_image1.normalize(0,1);
    in_image2.normalize(0,1);
    
    // get gradient
    const CImgList<> grad=img_3d.get_gradient("xyz");
    std::cout << "Step 2 : Get Gradients" << std::endl;

    std::cout << "Step 3 : Get A and b" << std::endl;
    Eigen::initParallel();
    auto [A,b]=fill_A_b(grad[0], grad[1], grad[2],1);

    std::cout << "Step 4 : Get u and v" << std::endl;
    Eigen::SparseMatrix<float> x=linear_solve(A,b);
    auto [u,v]=calculate_u_v(x,in_image1.height(),in_image1.width());

    // Normalize the flow field so that the largest vector (ui, vi) has magnitude 1
    // version 1
    double max_magnitude = 0.0;
    for (int y = 0; y < u.height(); ++y){
        for (int x = 0; x < u.width(); ++x){
            double magnitude = std::sqrt(u(x, y) * u(x, y) + v(x, y) * v(x, y));
            max_magnitude = std::max(max_magnitude, magnitude);
        }
    }
    u = u / max_magnitude;
    v = v / max_magnitude;
    

    //clamp and quantize
    //original version
    /*
    for (int i=0 ; i != u.height(); ++i){
        for(int j=0 ; j != u.width(); ++j){
            u(j,i)=std::max<unsigned char>(std::min<unsigned char>(255*(0.5*(u(j,i)+1)), 255), 0);
            v(j,i)=std::max<unsigned char>(std::min<unsigned char>(255*(0.5*(v(j,i)+1)), 255), 0);
        }
    }
    */
    
    
    // clamp and quantize version 2 
    for (int i=0 ; i != u.height(); ++i){
        for(int j=0 ; j != u.width(); ++j){
            u(j,i)=std::clamp<unsigned char>((255 * (0.5*(u(j,i)  + 1))), 0, 255);
            v(j,i)=std::clamp<unsigned char>((255 * (0.5*(v(j,i)  + 1))), 0, 255);
        }
    }
    
    
    std::cout << "Step 5 : Save u and v" << std::endl;
    u.save(argv[3]);
    v.save(argv[4]);

    return 0;

}
