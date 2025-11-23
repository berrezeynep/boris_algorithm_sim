#include <stdio.h>

float cube(float x);

int main(void) {
    float x = 4.0f;
    printf("The cube of x is: %f\n", cube(x));
    return 0;
}

float cube(float x) {
    return x * x * x;
}
