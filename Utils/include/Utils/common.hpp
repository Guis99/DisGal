#ifndef DISGAL_COMMON_INCLUDE
#define DISGAL_COMMON_INCLUDE

#include <QTM/include/QTM.hpp>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

#ifdef USE_SINGLE_PRECISION
    using Real_b = float;
    static_assert(sizeof(Real_b) == 4);
#else
    using Real_b = double;
    static_assert(sizeof(Real_b) == 8);
#endif

typedef std::vector<Eigen::Triplet<Real_b>> TrDV; // Vector of Real_b triplets for initializing sparse matrices
typedef Eigen::SparseMatrix<Real_b> SpD; // Dynamically-sized matrix of Real_b
typedef Eigen::Matrix<Real_b, Eigen::Dynamic, Eigen::Dynamic> DD; // Dynamically-sized vector of Real_b
typedef Eigen::Matrix<Real_b, Eigen::Dynamic, 1> DvD;

#ifndef VARIADIC_PRINT
#define VARIADIC_PRINT
template<typename... Args>
void var_print(Args&&... args) {
    (std::cout << ... << args) << std::endl; // Fold expression to print all arguments
}
#endif

#ifndef DEBUG_PRINT
    #ifdef VERBOSE
        template<typename... Args>
        void var_print(Args&&... args) {
            (std::cout << ... << args) << std::endl;
        }

        #define DEBUG_PRINT(...) var_print(__VA_ARGS__)
    #else
        #define DEBUG_PRINT(...)
    #endif
#endif

#ifdef MULTITHREAD_ENABLED
    #define CONFIRM_MT(numThreads) var_print("Multithreading active. Using ", numThreads, " threads")
    #define MT_ACTIVE 1
#else
    #define CONFIRM_MT(numThreads)
    #define MT_ACTIVE 0
#endif

#endif // includeguard