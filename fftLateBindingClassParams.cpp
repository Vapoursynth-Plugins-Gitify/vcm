

// late binding of fft dll
HINSTANCE hinstLib;

fftwf_malloc_proc fftwf_malloc;
fftwf_free_proc fftwf_free;
// for id fft
fftwf_plan_dft_r2c_1d_proc fftwf_plan_dft_r2c_1d;
fftwf_plan_dft_c2r_1d_proc fftwf_plan_dft_c2r_1d;
// for 2d fft
fftwf_plan_dft_r2c_2d_proc fftwf_plan_dft_r2c_2d;
fftwf_plan_dft_c2r_2d_proc fftwf_plan_dft_c2r_2d;


fftwf_destroy_plan_proc fftwf_destroy_plan;
fftwf_execute_proc fftwf_execute;

fftwf_execute_dft_r2c_proc  fftwf_execute_dft_r2c;
fftwf_execute_dft_c2r_proc  fftwf_execute_dft_c2r;

fftwf_init_threads_proc fftwf_init_threads;
fftwf_plan_with_nthreads_proc fftwf_plan_with_nthreads;
fftwf_cleanup_threads_proc fftwf_cleanup_threads;
fftwf_cleanup_proc fftwf_cleanup;
