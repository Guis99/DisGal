#ifndef DISGAL_COMMON_INCLUDE
#define DISGAL_COMMON_INCLUDE

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