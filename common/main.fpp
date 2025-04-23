program nonbond_kernel
#ifdef HAVE_PERF
    use perf_helper_mod
#endif
    use gparameter
    use module_pointers
    implicit none
    integer  :: step
    type(s_genesis_kernel_param):: gparam

    call read_data_file (gparam)

#ifdef HAVE_PERF
    call perf_initialize()
#endif

    do step=1,1000

      call kernel(gparam)
      
    end do

#ifdef HAVE_PERF
    call perf_finalize()
#endif
    write(6,*) 'time=',gparam%time
    call check_validation(gparam%MaxAtom,gparam%ncell,gparam%nthread,gparam%force)

end program nonbond_kernel


