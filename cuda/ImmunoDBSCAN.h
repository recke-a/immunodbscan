// the configured options and settings for Tutorial
#define ImmunoDBSCAN_VERSION_MAJOR 1
#define ImmunoDBSCAN_VERSION_MINOR 0 


#define GPUERRCHK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

