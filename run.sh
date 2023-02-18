#! /bin/bash
echo "This should run your application."

image1=$1;
image2=$2;
v=$3;
u=$4;

./executable $image1 $image2 $v $u
