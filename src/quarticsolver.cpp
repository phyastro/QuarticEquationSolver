#include <iostream>
#include <math.h>

#define PI 3.14159265358979323846

double sign(double x) {
    if (x != 0.0) {
        return x / abs(x);
    } else {
        return 1.0;
    }
}

struct complex{
    double real;
    double imaginary;
};

complex add(complex z, double c) {
    return complex(z.real + c, z.imaginary);
}

complex add(complex z1, complex z2) {
    return complex(z1.real + z2.real, z1.imaginary + z2.imaginary);
}

complex invert(complex z) {
    return complex(-z.real, -z.imaginary);
}

complex multiply(complex z, double c) {
    return complex(z.real * c, z.imaginary * c);
}

complex multiply(complex z1, complex z2) {
    return complex(z1.real * z2.real - z1.imaginary * z2.imaginary, z1.real * z2.imaginary + z1.imaginary * z2.real);
}

complex divide(complex z, double c) {
    return complex(z.real / c, z.imaginary / c);
}

complex divide(complex z1, complex z2) {
    return complex((z1.real * z2.real + z1.imaginary * z2.imaginary) / (z2.real * z2.real + z2.imaginary * z2.imaginary), (z1.imaginary * z2.real - z1.real * z2.imaginary) / (z2.real * z2.real + z2.imaginary * z2.imaginary));
}

double modulus(complex z) {
    return sqrt(z.real * z.real + z.imaginary * z.imaginary);
}

complex pow(complex z, double n) {
    // Convert Cartesian Form Into Polar Form
    double r = modulus(z);
    double theta_real = acos(z.real / r);
    double theta_imaginary = theta_real;
    if ((sign(z.imaginary)) < 0.0) {
        theta_imaginary = asin(z.imaginary / r);
    }

    // De Moivre's Theorem
    complex new_z = complex(cos(n * theta_real) * pow(r, n), sin(n * theta_imaginary) * pow(r, n));
    if ((new_z.real < 1e-8) and (new_z.real > -1e-8)) {
        new_z.real = 0.0;
    }
    if ((new_z.imaginary < 1e-8) and (new_z.imaginary > -1e-8)) {
        new_z.imaginary = 0.0;
    }

    return new_z;
}

void solveQuadratic(double a, complex b, complex c, complex &r1, complex &r2) {
    complex sqrtdiscriminant = pow(add(multiply(b, b),invert(multiply(c, 4.0 * a))), 0.5);
    r1 = divide(add(invert(sqrtdiscriminant), invert(b)), 2.0 * a);
    r2 = divide(add(sqrtdiscriminant, invert(b)), 2.0 * a);
}

void solveCubic(double a, double b, double c, double d, complex &r1, complex &r2, complex &r3) {
    double D0 = b * b - 3.0 * a * c;
    double D1 = 2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d;
    complex C = pow(divide(add(pow(complex(D1 * D1 - 4.0 * D0 * D0 * D0, 0.0), 0.5), D1), 2.0), 1.0 / 3.0);

    complex e1 = complex(-0.5, 0.5 * sqrt(3.0));
    complex e2 = complex(-0.5, -0.5 * sqrt(3.0));
    complex e1C = multiply(e1, C);
    complex e2C = multiply(e2, C);
    r1 = invert(divide(add(add(divide(complex(D0, 0.0), C), C), b), 3.0 * a));
    r2 = invert(divide(add(add(divide(complex(D0, 0.0), e1C), e1C), b), 3.0 * a));
    r3 = invert(divide(add(add(divide(complex(D0, 0.0), e2C), e2C), b), 3.0 * a));

    if ((r1.imaginary < 1e-12) and (r1.imaginary > -1e-12)) {
        r1.imaginary = 0.0;
    }
    if ((r2.imaginary < 1e-12) and (r2.imaginary > -1e-12)) {
        r2.imaginary = 0.0;
    }
    if ((r3.imaginary < 1e-12) and (r3.imaginary > -1e-12)) {
        r3.imaginary = 0.0;
    }
}

void solveQuartic(double alpha, double beta, double gamma, double delta, complex &a, complex &b, complex &c, complex &d) {
    complex r1 = complex(0.0, 0.0);
    complex r2 = complex(0.0, 0.0);
    complex r3 = complex(0.0, 0.0);
    solveCubic(1.0, -beta, alpha * gamma - 4.0 * delta, -(alpha * alpha * delta - 4.0 * beta * delta + gamma * gamma), r1, r2, r3);

    complex r11 = complex(0.0, 0.0);
    complex r12 = complex(0.0, 0.0);
    complex r21 = complex(0.0, 0.0);
    complex r22 = complex(0.0, 0.0);
    complex r31 = complex(0.0, 0.0);
    complex r32 = complex(0.0, 0.0);
    solveQuadratic(1.0, invert(r1), complex(delta, 0.0), r31, r32);
    solveQuadratic(1.0, invert(r2), complex(delta, 0.0), r11, r12);
    solveQuadratic(1.0, invert(r3), complex(delta, 0.0), r21, r22);
    a = divide(multiply(add(r21, r31), -gamma), add(multiply(add(r32, r22), add(add(r21, r31), r12)), multiply(r12, add(r21, r31))));
    b = divide(r11, a);
    c = divide(r21, a);
    d = divide(r31, a);
}

int main() {
    complex r1 = complex(0.0, 0.0);
    complex r2 = complex(0.0, 0.0);
    complex r3 = complex(0.0, 0.0);
    complex r4 = complex(0.0, 0.0);
    solveQuartic(-2.0, 0.0, 2.0, 0.0, r1, r2, r3, r4);
    std::cout << r1.real << " + " << r1.imaginary << "i" << std::endl;
    std::cout << r2.real << " + " << r2.imaginary << "i" << std::endl;
    std::cout << r3.real << " + " << r3.imaginary << "i" << std::endl;
    std::cout << r4.real << " + " << r4.imaginary << "i" << std::endl;

    return 0;
}