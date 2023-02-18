#include <iostream>
#include <vector>

# define cimg_display 0 
# include "CImg-v.3.1.6/CImg.h"
using namespace cimg_library;


class Tensor{
    public:
        Tensor();
        Tensor(CImg <float> input_image);
        ~Tensor();
    
        std::vector< float > data;
        int width = 0;
        int height = 0;
        long unsigned int num_pixel = 0;

    
};

// TODO check dimension
bool operator==(const Tensor &a, const Tensor &b);
double calculate_x(const Tensor &a, int i, int j);
double calculate_y(const Tensor &a, int i, int j);

Tensor diff_x(const Tensor &a, const Tensor &b);
Tensor diff_y(const Tensor &a, const Tensor &b);
Tensor diff_t(const Tensor &a, const Tensor &b);



