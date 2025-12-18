#include <iostream>
#include <array>
#include <autodiff/forward/dual.hpp>

using autodiff::dual2nd;

/*
// --------------------------------------------------------
// Displacement functor u(x,y) → {u1, u2}
// Replace this with your own function
// --------------------------------------------------------
struct U {
    template<typename T>
    std::array<T,2> operator()(const T& x, const T& y) const {
        // Example: u1 = x^2 y,  u2 = x y^2
        return { x*x*y, x*y*y };
    }
};

// --------------------------------------------------------
// 3×3 elasticity matrix (Voigt ordering)
// --------------------------------------------------------
struct DMatrix {
    double D[3][3];

    template<typename T>
    std::array<T,3> operator*(const std::array<T,3>& e) const {
        std::array<T,3> s{};
        for(int i=0;i<3;i++){
            s[i] = T(0);
            for(int j=0;j<3;j++)
                s[i] += D[i][j] * e[j];
        }
        return s;
    }
};

// --------------------------------------------------------
// Compute stress and div(sigma) with autodiff
// --------------------------------------------------------
void computeStressAndDiv(const U& ufun, const DMatrix& D, double x0, double y0)
{
    using T = dual2nd;   // 2nd order automatic differentiation

    // promote x,y to dual numbers to compute derivatives
    T x = x0;
    T y = y0;

    // displacement
    auto u = ufun(x, y);
    auto u1 = u[0];
    auto u2 = u[1];

    // First derivatives (strain components)
    T ux_x = derivative([&](auto xx){ return ufun(xx, y)[0]; }, wrt(x), x);
    T ux_y = derivative([&](auto yy){ return ufun(x, yy)[0]; }, wrt(y), y);
    T uy_x = derivative([&](auto xx){ return ufun(xx, y)[1]; }, wrt(x), x);
    T uy_y = derivative([&](auto yy){ return ufun(x, yy)[1]; }, wrt(y), y);

    // Strain in Voigt notation: [εxx, εyy, γxy]
    std::array<T,3> eps {
        ux_x,
        uy_y,
        ux_y + uy_x
    };

    // Stress = D * strain
    auto sigma = D * eps;
    T sxx = sigma[0];
    T syy = sigma[1];
    T sxy = sigma[2];

    // Divergence of stress
    T div1 = derivative([&](auto xx){ 
        return (D * std::array<T,3>{
            derivative([&](auto x2){ return ufun(x2, y)[0]; }, wrt(xx), xx),
            derivative([&](auto y2){ return ufun(xx, y2)[1]; }, wrt(y), y),
            derivative([&](auto y2){ return ufun(xx, y2)[0]; }, wrt(y), y)
              + derivative([&](auto x2){ return ufun(x2, y)[1]; }, wrt(xx), xx)
        })[0]; 
    }, wrt(x), x)
    +
    derivative([&](auto yy){
        return (D * std::array<T,3>{
            derivative([&](auto x2){ return ufun(x2, yy)[0]; }, wrt(x), x),
            derivative([&](auto y2){ return ufun(x, y2)[1]; }, wrt(yy), yy),
            derivative([&](auto y2){ return ufun(x, y2)[0]; }, wrt(yy), yy)
              + derivative([&](auto x2){ return ufun(x2, yy)[1]; }, wrt(x), x)
        })[2];
    }, wrt(y), y);

    T div2 = derivative([&](auto xx){
        return (D * std::array<T,3>{
            derivative([&](auto x2){ return ufun(x2, y)[0]; }, wrt(xx), xx),
            derivative([&](auto y2){ return ufun(xx, y2)[1]; }, wrt(y), y),
            derivative([&](auto y2){ return ufun(xx, y2)[0]; }, wrt(y), y)
              + derivative([&](auto x2){ return ufun(x2, y)[1]; }, wrt(xx), xx)
        })[2];
    }, wrt(x), x)
    +
    derivative([&](auto yy){
        return (D * std::array<T,3>{
            derivative([&](auto x2){ return ufun(x2, yy)[0]; }, wrt(x), x),
            derivative([&](auto y2){ return ufun(x, y2)[1]; }, wrt(yy), yy),
            derivative([&](auto y2){ return ufun(x, y2)[0]; }, wrt(yy), yy)
              + derivative([&](auto x2){ return ufun(x2, yy)[1]; }, wrt(x), x)
        })[1];
    }, wrt(y), y);

    std::cout << "strain = [" 
              << val(eps[0]) << ", " 
              << val(eps[1]) << ", " 
              << val(eps[2]) << "]\n";

    std::cout << "stress = [" 
              << val(sxx) << ", " 
              << val(syy) << ", " 
              << val(sxy) << "]\n";

    std::cout << "div sigma = ("
              << val(div1) << ", "
              << val(div2) << ")\n";
}

// --------------------------------------------------------
int main()
{
    // Example: plane strain isotropic D-matrix
    double E = 200e9;
    double nu = 0.3;
    double c = E/((1+nu)*(1-2*nu));

    DMatrix D {{
        { c*(1-nu), c*nu,        0 },
        { c*nu,     c*(1-nu),    0 },
        { 0,        0,     c*(1-2*nu)/2 }
    }};

    U u;   // displacement field

    computeStressAndDiv(u, D, 1.0, 2.0);

    return 0;
}
    */
   /*
// C++ includes
#include <iostream>
using namespace std;

// autodiff include
#include <iostream>
#include <autodiff/reverse/var.hpp>
// #include <autodiff/reverse/var/eigen.hpp> // Includes support for variables and derivatives

using namespace autodiff;

// A sample function of two variables
var function_f(var x, var y) {
    return x * x + y * y + exp(x * y);
}

int main() {
    // 1. Define input variables with their values
    var x = 1.0;
    var y = 2.0;

    // 2. Mark variables as independent and evaluate the function
    var u = function_f(x, y);

    // 3. Compute all derivatives (up to the order supported by the 'var' implementation)
    // The 'derivatives' function traverses the graph once to get all first derivatives
    auto dud = derivativesx(u, wrt(x, y));

    // Extracting the first derivatives
    auto dudx = dud[0]; // ∂u/∂x
    auto dudy = dud[1]; // ∂u/∂y

    std::cout << "Value of u: " << u << std::endl;
    std::cout << "du/dx at (1, 2): " << dudx << std::endl;
    std::cout << "du/dy at (1, 2): " << dudy << std::endl;

    return 0;
}
*/
/*
#include <autodiff/reverse/var.hpp>
using namespace autodiff;

std::vector<var> f(std::vector<var>& x)
{
    var u = sin(x[0]) * cos(x[1]);
    var v = u * u;
    return {u, v};
}
int main()
{
    var x = 0.5;  // the input variable x
    var y = 0.1;  // the input variable y
    std::vector<var> vars = {x, y};
    auto u = f(vars);

    std::vector<autodiff::reverse::detail::Variable<double> > vars2(2);
    auto [ux, uy] = derivativesx(u[0], wrt(x, y));  // evaluate the first order derivatives of u
    vars2[0] = ux;
    auto [uxx, uxy] = derivativesx(vars2[0], wrt(x, y));  // evaluate the second order derivatives of ux

    cout << "u = " << u[0] << " " << u[1] << endl;  // print the evaluated output variable u
    cout << "ux(autodiff) = " << ux << endl;  // print the evaluated first order derivative ux
    cout << "ux(exact) = " << 1 - 2*sin(x)*sin(x) << endl;  // print the exact first order derivative ux
    cout << "uxx(autodiff) = " << uxx << endl;  // print the evaluated second order derivative uxx
    cout << "uxx(exact) = " << -4*cos(x)*sin(x) << endl;  // print the exact second order derivative uxx
}

*/

// C++ includes
#include <iostream>

// autodiff include
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
using namespace autodiff;

// The scalar function for which the gradient is needed
template<int val>
dual2nd f(const Vector3dual2nd& x)
{
    Vector3dual2nd disp;
    for (int i = 0; i < x.size(); ++i) {
        disp[0] = cos(M_PI * x[0]) * sin(2 * M_PI * x[1]);
        disp[1] = cos(M_PI * x[1])* sin(M_PI * x[0]);
        disp[2] = cos(M_PI * x[2])* sin(M_PI * x[1]);
    }
    return disp[val];
}

int main()
{
    using Eigen::Matrix3d;

    Vector3dual2nd x(3); // the input vector x with 3 variables
    x << 0.2, 0.3, 0.6; // x = [0.2, 0.3, 0.6]

    std::vector<dual2nd> u(3); // the output scalar u = f(x) evaluated together with Hessian below
    std::vector<Vector3dual> g(3);
    std::vector<Matrix3d> H(3);

    #define Hess(dof) H[dof] = hessian(f<dof>, wrt(x), at(x), u[dof], g[dof]); // evaluate the function value u and its Hessian matrix H
    Hess(0)
    Hess(1)
    Hess(2)
    #undef Hess

    for(int i=0; i<3; i++) {
        std::cout << "----------------------------------" << std::endl;
        std::cout << "Evaluating at dof = " << i << std::endl;
        std::cout << "u = "   << u[i] << std::endl; // print the evaluated output u
        std::cout << "g = \n" << g[i] << std::endl; // print the evaluated gradient vector g = du/dx
        std::cout << "H = \n" << H[i] << std::endl; // print the evaluated Hessian matrix H = d²u/dx²
    }
}
