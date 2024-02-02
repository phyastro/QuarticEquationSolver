#include <iostream>
#include <math.h>

static const bool SOLVERTYPEREAL = true;

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

struct vec3{
    double x;
    double y;
    double z;
};

struct vec4{
    double x;
    double y;
    double z;
    double w;
};

struct bvec3{
    bool x;
    bool y;
    bool z;
};

struct bvec4{
    bool x;
    bool y;
    bool z;
    bool w;
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
    double theta = atan2(z.imaginary, z.real);

    // De Moivre's Theorem
    complex new_z = complex(cos(n * theta) * pow(r, n), sin(n * theta) * pow(r, n));
    if ((new_z.real < 1e-12) and (new_z.real > -1e-12)) {
        new_z.real = 0.0;
    }
    if ((new_z.imaginary < 1e-12) and (new_z.imaginary > -1e-12)) {
        new_z.imaginary = 0.0;
    }

    return new_z;
}

void solveQuadratic(double a, complex b, complex c, complex &r1, complex &r2) {
    complex sqrtdiscriminant = pow(add(multiply(b, b), invert(multiply(c, 4.0 * a))), 0.5);
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

void solveQuartic(double a4, double a3, double a2, double a1, double a0, complex &a, complex &b, complex &c, complex &d) {
    double alpha = a3 / a4;
	double beta = a2 / a4;
	double gamma = a1 / a4;
	double delta = a0 / a4;

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

bvec3 solveCubic(double b, double c, double d, vec3 &roots) {
    // https://arxiv.org/abs/1903.10041
    // Solves Cubic Equation By Combining Two Different Methods To Find Real Roots
    double bdiv3 = b / 3.0;
    double Q = c / 3.0 - bdiv3 * bdiv3;
    double R = 0.5 * bdiv3 * c - bdiv3 * bdiv3 * bdiv3 - 0.5 * d;
    if ((Q == 0.0) && (R == 0.0)) {
        // If All The Roots Of The Cubic Equation Are Equal
        roots = {-bdiv3, -bdiv3, -bdiv3};
        return {true, true, true};
    }
    double D = Q * Q * Q + R * R;
    if (D > 0.0) {
        double S = cbrt(R + sqrt(D));
        double T = cbrt(R - sqrt(D));
        roots.x = S + T - bdiv3;
        return {true, false, false};
    }
    double sqrtnegQ = sqrt(-Q);
    double thetadiv3 = acos(R / (sqrtnegQ * sqrtnegQ * sqrtnegQ)) / 3.0;
    roots.x = 2.0 * sqrtnegQ * cos(thetadiv3) - bdiv3;
    roots.y = 2.0 * sqrtnegQ * cos(thetadiv3 + 2.0 * PI / 3.0) - bdiv3;
    roots.z = 2.0 * sqrtnegQ * cos(thetadiv3 + 4.0 * PI / 3.0) - bdiv3;
    return {true, true, true};
}

bvec4 solveQuartic(double a, double b, double c, double d, double e, vec4 &roots) {
    // https://www.mdpi.com/2227-7390/10/14/2377
    // Solves Quartic Equation By Using The Method Given In The Above Paper
    double inva = 1.0 / a;
    double inva2 = inva * 0.5;
    double inva2a2 = inva2 * inva2;
    double bb = b * b;
    double p = -1.5 * bb * inva2a2 + c * inva;
    double q = bb * b * inva2a2 * inva2 - b * c * inva * inva2 + d * inva;
    double r = -0.1875 * bb * bb * inva2a2 * inva2a2 + 0.5 * c * bb * inva2a2 * inva2 - b * d * inva2a2 + e * inva;
    vec3 s = {0.0, 0.0, 0.0};
    solveCubic(0.5 * -p, -r, 0.5 * p * r - 0.125 * q * q, s);
    double s2subp = 2.0 * s.x - p;
    if (s2subp < 0.0) {
        return {false, false, false, false};
    }
    double invs2subp= -2.0 * s.x - p;
    double sqrts2subp = sqrt(s2subp);
    double q2divsqrt = 2.0 * q / sqrts2subp;
    double invaddq2div = invs2subp + q2divsqrt;
    double invsubq2div = invs2subp - q2divsqrt;
    double bdiv4a = 0.25 * inva * b;
    bvec4 isReal = {false, false, false, false};
    if (invaddq2div >= 0.0) {
        double sqrtinvadd = sqrt(invaddq2div);
        roots.x = 0.5 * (-sqrts2subp + sqrtinvadd) - bdiv4a;
        roots.y = 0.5 * (-sqrts2subp - sqrtinvadd) - bdiv4a;
        isReal.x = true;
        isReal.y = true;
    }
    if (invsubq2div >= 0.0) {
        double sqrtinvsub = sqrt(invsubq2div);
        roots.z = 0.5 * (sqrts2subp + sqrtinvsub) - bdiv4a;
        roots.w = 0.5 * (sqrts2subp - sqrtinvsub) - bdiv4a;
        isReal.z = true;
        isReal.w = true;
    }
    return isReal;
}

int main() {
    if (!SOLVERTYPEREAL) {
        // Complex Roots Solver For Quartic Equation With Real Coefficients
        complex r1 = complex(0.0, 0.0);
        complex r2 = complex(0.0, 0.0);
        complex r3 = complex(0.0, 0.0);
        complex r4 = complex(0.0, 0.0);
        solveQuartic(0.01, 1.0, -1.0, 2.0, 1.0, r1, r2, r3, r4);
        std::cout << r1.real << " + " << r1.imaginary << "i" << std::endl;
        std::cout << r2.real << " + " << r2.imaginary << "i" << std::endl;
        std::cout << r3.real << " + " << r3.imaginary << "i" << std::endl;
        std::cout << r4.real << " + " << r4.imaginary << "i" << std::endl;
    }

    if (SOLVERTYPEREAL) {
        // Real Roots Solver For Quartic Equation With Real Coefficients
        vec4 roots = {0.0, 0.0, 0.0, 0.0};
        bvec4 isReal = solveQuartic(1.5, 1.5, -1.9, -0.5, 0.3, roots);
        if (isReal.x)
            std::cout << roots.x << std::endl;
        if (isReal.y)
            std::cout << roots.y << std::endl;
        if (isReal.z)
            std::cout << roots.z << std::endl;
        if (isReal.w)
            std::cout << roots.w << std::endl;
    }

    return 0;
}