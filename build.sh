echo "The project should be built now..."
g++ -Wall -Werror -fsanitize=address -g -O3 -fopenmp -std=c++20 -L/usr/X11R6/lib -lm -lpthread -o executable opticalFlow.cpp

