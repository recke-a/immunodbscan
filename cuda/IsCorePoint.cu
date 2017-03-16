// Dieses ist eine Testimplementation der Kernsuche
// Es sollen verschiedene Methoden untersucht werden.

#include <iostream>
#include <math.h>
#include <thrust/device_vector.h>
#include <thrust/count.h>
#include <thrust/execution_policy.h>
#include <thrust/reduce.h>
#include <curand.h>
#include <curand_kernel.h>

#include <chrono>


__global__ void setup_kernel(thrust::device_vector<curandState>::iterator state, int maxN, unsigned long seed)
{
    int id = threadIdx.x + blockIdx.x * 64;
    if (id < maxN)
    {
	    /* Each thread gets same seed, a different sequence 
	       number, no offset */
	    curand_init(seed, id, 0, thrust::raw_pointer_cast(&(*(state+id))));
    }
}


__global__ void fill_states(thrust::device_vector<curandState>::iterator state, thrust::device_vector<bool>::iterator result, double threshold, int maxN)
{
	int id = threadIdx.x + blockIdx.x * 64;
	if (id < maxN)
    {
		*(result+id) = threshold > curand_uniform(thrust::raw_pointer_cast(&(*(state+id))));
	}
}

__global__ void search_kernel_1(thrust::device_vector<bool>::iterator input_list, int maxN, int required, int *found, int *ueberfluessig)
{
	if (*found < required)
	{
		int id = threadIdx.x + blockIdx.x * 64;
		if (id < maxN)
	    {
			unsigned int x = 0; 
			while (x < 100000) ++x;
			if (*(input_list + id)) atomicAdd(found, 1);
		}
	} else
	{
		atomicAdd(ueberfluessig, 1);
	}
}

__global__ void search_kernel_2(thrust::device_vector<bool>::iterator input_list, int maxN, int required, int *found, int repeated)
{
	int basis = 0;
	while (*found < required && basis < repeated)
	{
		int id = basis + threadIdx.x * repeated + blockIdx.x * 64 * repeated;
		if (id < maxN)
	    {
			unsigned int x = 0; 
			while (x < 100000) ++x;
			if (*(input_list + id)) atomicAdd(found, 1);
		}
		++basis;
	} 
}


int main()
{
	const int N=100000;
	thrust::device_vector<curandState> devStates(N);
	thrust::device_vector<bool> dv(N);
	
	double threshold = 0.1;
	setup_kernel<<<(N + 64 - 1)/64, 64>>>(devStates.begin(), N, unsigned(time(NULL)));
	fill_states<<<(N + 64 - 1)/64, 64>>>(devStates.begin(), dv.begin(), threshold, N);
	
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	// calculations
	unsigned summe = thrust::count(thrust::device, dv.begin(), dv.end(), true);
	// finito
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "Zeit für Suche: " << time_span.count() << " seconds." << std::endl;
    std::cout << "Summe bei Schwelle " << threshold << " ist " << summe << "\n";	
	
	
	t1 = std::chrono::high_resolution_clock::now();
	// calculations
	int ergebnis = 0;
	thrust::device_ptr<int> dp_result = thrust::device_malloc<int>(1);
	thrust::device_ptr<int> dp_ueberfluessig = thrust::device_malloc<int>(1);
	// cudaMemcpy(thrust::raw_pointer_cast(dp_result), &ergebnis, sizeof(int), cudaMemcpyHostToDevice);
	*dp_result = ergebnis;
	*dp_ueberfluessig = 0;
	search_kernel_1<<<N,1>>>(dv.begin(), N, 10, thrust::raw_pointer_cast(dp_result), thrust::raw_pointer_cast(dp_ueberfluessig));
	//search_kernel_2<<<(N+1000-1)/1000, 1>>>(dv.begin(), N, 10, thrust::raw_pointer_cast(dp_result), 1000);
	// cudaMemcpy(&ergebnis, thrust::raw_pointer_cast(dp_result),  sizeof(int), cudaMemcpyDeviceToHost);
	ergebnis = *dp_result;
	int ueberfluessig = *dp_ueberfluessig;
	//thrust::device_free(dp_result);
	//thrust::device_free(dp_ueberfluessig);
	// finito
	t2 = std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "Zeit für Suche: " << time_span.count() << " seconds." << std::endl;
    std::cout << "Ergebnis = " << ergebnis << "\n";	
    std::cout << "Insgesamt ueberfluessig = " << ueberfluessig << "\n";	
	
	   
	
	return 0;
}
