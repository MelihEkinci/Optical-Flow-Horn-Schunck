# include "iostream"
# include "opticalFlow.hpp"

int main(){
    // create image objects
    //call caculate funtion to up
    CImg < double > in_image1("test_images/rects_640x480_ref_v.bmp");
    std::cout << "hallo";
    std::cout << static_cast<int>(in_image1(60,200));
    //std::cout << in_image1.width() << std::endl;

    Tensor *pic1 = new Tensor(in_image1);
    //std::cout << pic1->width << std::endl;
    //std::cout << pic1->height << std::endl;

    //double x =  pic1.data[30];
    //std::cout << static_cast<int>(x);  
    //vec = std::vector<double>(5) ;
    //vec.push_back(52.5);
    //std::cout << vec[0];
}
