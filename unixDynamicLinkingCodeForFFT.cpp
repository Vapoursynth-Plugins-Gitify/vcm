//Mac OS GNUC UNIX dynamic linking constructor code
//#include <dlfcn.h>
//#include <string>
bool okay = true;
std::string sDLLName = "libfftw3f3.3.so";
/*
dlopen(sDLLName.c_str(), 2);
// the above returns a void pointer to the library
return dlsym(Lib, FnName);
// close free library
dlclose(void* hdll)
*/
// mac shared library constructor code
void* hDLL;
hDLL = dlopen(sDLLName.c_str(), 2);
if (hDll != 0)
{
	d->hinstLib = hDll;
	d->fftwf_free = dlsym(d->hinstLib, "fftwf_free");
	d->fftwf_malloc = dlsym(d->hinstLib, "fftwf_malloc");
	// for 1 d fft
	d->fftwf_plan_dft_r2c_1d = dlsym(d->hinstLib, "fftwf_plan_dft_r2c_1d");
	d->fftwf_plan_dft_c2r_1d = dlsym(d->hinstLib, "fftwf_plan_dft_c2r_1d");
	// for 2d fft
	d->fftwf_plan_dft_r2c_2d = dlsym(d->hinstLib, "fftwf_plan_dft_r2c_2d");
	d->fftwf_plan_dft_c2r_2d = dlsym(d->hinstLib, "fftwf_plan_dft_c2r_2d");

	d->fftwf_destroy_plan = dlsym(d->hinstLib, "fftwf_destroy_plan");
	d->fftwf_execute = dlsym(d->hinstLib, "fftwf_execute");

	d->fftwf_execute_dft_r2c = dlsym(d->hinstLib, "fftwf_execute_dft_r2c");
	d->fftwf_execute_dft_c2r = dlsym(d->hinstLib, "fftwf_execute_dft_c2r");
	// required for internal threading
//	fftwf_init_threads = dlsym(hinstLib, "fftwf_init_threads");
//	fftwf_plan_with_nthreads = dlsym(hinstLib, "fftwf_plan_with_nthreads");
//	fftwf_cleanup_threads = dlsym(hinstLib, "fftwf_cleanup_threads");
//	fftwf_cleanup = dlsym(hinstLib, "fftwf_cleanup");

}
if (d->hinstLib == NULL || d->fftwf_free == NULL || d->fftwf_malloc == NULL || d->fftwf_plan_dft_r2c_2d == NULL ||
	d->fftwf_plan_dft_c2r_2d == NULL || d->fftwf_destroy_plan == NULL || d->fftwf_execute == NULL
	|| d->fftwf_execute_dft_r2c == NULL || d->fftwf_execute_dft_c2r == NULL
	// required for internal threading
//	|| d->fftwf_init_threads == NULL || d->fftwf_plan_with_nthreads == NULL
//	|| d->fftwf_cleanup_threads == NULL || d->fftwf_cleanup == NULL
)
{
	okay = false;

}

