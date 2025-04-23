
subroutine kernel(gparam)
#ifdef HAVE_PERF
    use perf_helper_mod
#endif
    use gparameter
    use  module_pointers
    implicit none
    ! formal arguments
    type(s_genesis_kernel_param),   target, intent(inout)    :: gparam

    ! local variables
    real(wp)                  :: dij(1:nxyz)
    real(wp)                  :: dij1,dij2,dij3, rij2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  ::  grad_coef
    real(wp)                  :: work(1:nxyz)
    real(wp)                  :: rtmp(1:nxyz), qtmp, jqtmp
    real(wp)                  :: force_local(1:nxyz)
    real(wp)                  :: viri(1:nxyz)
    real(wp)                  :: ieps, jeps, eps, irmin, jrmin, rmin
    integer                  :: i, ix, iy, j, k, ij, L, L1, ii, ki, kki
#ifdef OMP_PARALLEL
    integer                  :: omp_get_thread_num
#endif
    integer                  :: id, ik
    integer                  :: iatmcls,jatmcls
    integer                  :: maxcell,ncell,ncell_local
    integer                  :: num_atom_domain
    integer                  :: nthread
    real(wp)                  :: density, cutoff, cutoff2
    real(wp)                  :: inv_MaxAtom
    integer,   save           :: init_flg = 0
    integer*8                :: CountPerSec, CountMax
    real(dp), save           :: ratio_usec
    real(dp), save           :: time = 0
    real(dp)                 :: time_st, time_en
    integer*8 time_c

!    real(8)                  :: Val, expected_result, check_result

    integer, pointer,contiguous   :: natom(:)
    integer, pointer,contiguous   :: start_atom(:)
    integer, pointer,contiguous   :: atmcls(:,:)
    integer, pointer,contiguous   :: num_nb15_calc(:,:)
    integer, pointer,contiguous   :: nb15_calc_list(:,:,:)
    integer, pointer,contiguous   :: atmcls_pbc(:)
    real(wp), pointer,contiguous  :: charge(:,:)
    real(wp), pointer,contiguous  :: coord(:,:,:)
    real(wp), pointer,contiguous  :: coord_pbc(:,:,:)
    real(wp), pointer,contiguous  :: trans1(:,:,:)
    real(wp), pointer,contiguous  :: table_grad(:)
    real(wp), pointer,contiguous  :: force(:,:,:,:)
    real(wp), pointer,contiguous  :: nonb_lj12(:,:)
    real(wp), pointer,contiguous  :: nonb_lj6(:,:)
    real(wp), pointer,contiguous  :: virial(:,:,:)

    maxcell         = gparam%maxcell
    ncell           = gparam%ncell
    num_atom_domain = gparam%num_atom_domain
    nthread         = gparam%nthread
    cutoff          = gparam%cutoff 
    cutoff2         = gparam%cutoff2
    density         = gparam%density

    coord              => gparam%coord
    coord_pbc          => gparam%coord_pbc
    trans1             => gparam%trans1
    natom              => gparam%natom
    start_atom         => gparam%start_atom
    atmcls             => gparam%atmcls
    num_nb15_calc      => gparam%num_nb15_calc
    nb15_calc_list     => gparam%nb15_calc_list
    atmcls_pbc         => gparam%atmcls_pbc
    charge             => gparam%charge     
    table_grad         => gparam%table_grad
    nonb_lj12          => gparam%nonb_lj12
    nonb_lj6           => gparam%nonb_lj6 
    virial             => gparam%virial   
    force              => gparam%force    

    force=0.0
    coord_pbc=0.0
    virial=0.0

    if( init_flg == 0 ) then
            init_flg = 1
            call system_clock(time_c, CountPerSec, CountMax)
            ratio_usec = 1.0d+6 / dble(CountPerSec)
    else
            call system_clock(time_c)
    endif
    time_st = dble(time_c)*ratio_usec

#ifdef OMP_PARALLEL
    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, k, ki, kki, iy, ij, j, rij2, L, L1, R, rtmp,  &
    !$omp         qtmp, jqtmp, iatmcls, jatmcls, term_lj12, term_lj6,      &
    !$omp         term_elec, grad_coef, work, force_local, lj6, lj12,      &
    !$omp         dij, viri)
    id = omp_get_thread_num()
#else
    id = 0
#endif

#ifdef HAVE_PERF
    call perf_start_section(0)
#endif

    do i = id+1, ncell, nthread
      ki = start_atom(i)
#ifdef DIR
!ocl norecurrence
!!ocl simd
!ocl swp
!GCC$ ivdep
!GCC$ vector 
#endif
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(kki,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(kki,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(kki,3,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(kki,4,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

#ifdef OMP_PARALLEL
   !$omp barrier
#endif

    do kki = id+1, num_atom_domain, nthread

      rtmp(1) = coord_pbc(kki,1,1)
      rtmp(2) = coord_pbc(kki,2,1)
      rtmp(3) = coord_pbc(kki,3,1)
      qtmp    = coord_pbc(kki,4,1)
      iatmcls = atmcls_pbc(kki)

      force_local(1) = 0.0
      force_local(2) = 0.0
      force_local(3) = 0.0
      viri(1)        = 0.0
      viri(2)        = 0.0
      viri(3)        = 0.0

#ifdef DIR
!ocl norecurrence
!!ocl simd
!ocl swp
!GCC$ ivdep
!GCC$ vector 
#endif
      do k = 1, num_nb15_calc(kki,1)

        ij = nb15_calc_list(k,kki,1)

        dij(1) = rtmp(1) - coord_pbc(ij,1,1)
        dij(2) = rtmp(2) - coord_pbc(ij,2,1)
        dij(3) = rtmp(3) - coord_pbc(ij,3,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2    = cutoff2*density / rij2

        jatmcls = atmcls_pbc(ij)
        lj12    = nonb_lj12(jatmcls,iatmcls)
        lj6     = nonb_lj6 (jatmcls,iatmcls)
        jqtmp   = coord_pbc(ij,4,1)

        L  = int(rij2)
        R  = rij2 - L
        L1 = 3*L - 2

        term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1))
        term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
        term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
        grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(ij,1,1,id+1) = force(ij,1,1,id+1) + work(1)
        force(ij,2,1,id+1) = force(ij,2,1,id+1) + work(2)
        force(ij,3,1,id+1) = force(ij,3,1,id+1) + work(3)

      end do

      force(kki,1,1,id+1) = force(kki,1,1,id+1) + force_local(1)
      force(kki,2,1,id+1) = force(kki,2,1,id+1) + force_local(2)
      force(kki,3,1,id+1) = force(kki,3,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    end do

#ifdef HAVE_PERF
    call perf_stop_section(0)
#endif
#ifdef OMP_PARALLEL
    !$omp end parallel
#endif

    call system_clock(time_c)
    time_en = dble(time_c)*ratio_usec
    time = time + (time_en-time_st)*1.0d-6
    gparam%time = time


!    PROF_STOP("Nonb15F")
!    PROF_STOP_ALL
!    PROF_FINALIZE

    return

end subroutine kernel

