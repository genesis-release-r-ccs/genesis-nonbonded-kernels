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
    real(wp)                  :: dij(1:nxyz), rij2
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:nxyz)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:nxyz), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(nxyz)
    real(wp)                  :: save_table(6)
    real(wp)                  :: dij_list(3,262)
    real(wp)                  :: rij2_list(262), force_localj(3,262)
    real(wp)                  :: rij2_inv, within_cutoff
    integer                   :: i, ix, iy, j, k, ij, iix, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id
#ifdef OMP_PARALLEL
    integer                  :: omp_get_thread_num
#endif
    integer                   :: ncell, ncell_local
    integer                   :: j_list(262)

    integer,   save           :: init_flg = 0
    integer*8                :: CountPerSec, CountMax
    real(dp), save           :: ratio_usec
    real(dp), save           :: time = 0
    real(dp)                 :: time_st, time_en
    integer*8 time_c

    integer, pointer,contiguous   :: natom(:)
    integer, pointer,contiguous   :: atmcls(:,:)
    integer, pointer,contiguous   :: num_nb15_calc1(:,:)
    integer, pointer,contiguous   :: num_nb15_calc(:,:)
    integer, pointer,contiguous   :: nb15_calc_list1(:,:)
    integer, pointer,contiguous   :: nb15_calc_list(:,:)
    integer, pointer,contiguous   :: nb15_list(:,:)
    integer, pointer,contiguous   :: nb15_cell(:)
    integer, pointer,contiguous   :: cell_pairlist(:,:)
    integer, pointer,contiguous   :: cell_move(:,:,:)
    real(wp), pointer,contiguous  :: charge(:,:)
    real(wp), pointer,contiguous  :: coord(:,:,:)
    real(wp), pointer,contiguous  :: coord_pbc(:,:,:)
    real(wp), pointer,contiguous  :: system_size(:)
    real(wp), pointer,contiguous  :: trans1(:,:,:)
    real(wp), pointer,contiguous  :: table_grad(:)
    real(wp), pointer,contiguous  :: force(:,:,:,:)
    real(wp), pointer,contiguous  :: nonb_lj12(:,:)
    real(wp), pointer,contiguous  :: nonb_lj6(:,:)
    real(wp), pointer,contiguous  :: virial(:,:)

    real(wp) :: density, cutoff
    integer :: maxcell, num_atom_domain, nthread

    maxcell         = gparam%maxcell
    ncell           = gparam%ncell
    ncell_local     = gparam%ncell_local
    num_atom_domain = gparam%num_atom_domain
    nthread         = gparam%nthread
    cutoff          = gparam%cutoff 
    cutoff2         = gparam%cutoff2
    density         = gparam%density

    coord              => gparam%coord
    coord_pbc          => gparam%coord_pbc
    trans1             => gparam%trans1
    cell_move          => gparam%cell_move
    system_size        => gparam%system_size
    natom              => gparam%natom
    num_nb15_calc      => gparam%num_nb15_calc
    num_nb15_calc1     => gparam%num_nb15_calc1
    nb15_calc_list     => gparam%nb15_calc_list
    nb15_calc_list1    => gparam%nb15_calc_list1
    nb15_list          => gparam%nb15_list      
    nb15_cell          => gparam%nb15_cell      
    cell_pairlist      => gparam%cell_pairlist
    atmcls             => gparam%atmcls
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
    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef, work, &
    !$omp         ij, j, trans_x, trans_y, trans_z, iix, force_local, lj6,     &
    !$omp         lj12, L1, elec_temp, evdw_temp, dij, &
    !$omp          dij_list, rij2_list, force_localj, j_list)
    id = omp_get_thread_num()
#else
    id = 0
#endif

#ifdef HAVE_PERF
    call perf_start_section(0)
#endif

    do i = id+1, ncell, nthread
#ifdef DIR
!ocl norecurrence
!ocl simd
!ocl swp
!GCC$ ivdep
!GCC$ vector 
#endif
      do ix = 1, natom(i)
        coord_pbc(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

#ifdef OMP_PARALLEL
   !$omp barrier
#endif
    ! energy within a cell
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, natom(i) - 1

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

#ifdef DIR
!ocl norecurrence
!ocl simd
!ocl swp
!GCC$ ivdep
!GCC$ vector 
#endif
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)


          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          
          rij2  = cutoff2*density/rij2
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))

          L    = int(rij2)
          R    = rij2 - L
          L1   = 3*L - 2


          term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)

          force(iy,1,i,id+1) = force(iy,1,i,id+1) + work(1)
          force(iy,2,i,id+1) = force(iy,2,i,id+1) + work(2)
          force(iy,3,i,id+1) = force(iy,3,i,id+1) + work(3)
        end do

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)

      end do

    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)
!      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix   = nb15_list(iix,ij)

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15 = fin_nb15

        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

#ifdef DIR
!ocl norecurrence
!ocl simd
!ocl swp
!GCC$ ivdep
!GCC$ vector 
#endif
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          rij2  = cutoff2*density/rij2
          lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))

          L    = int(rij2)
          R    = rij2 - L

          L1   = 3*L - 2
          term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))

          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)

          force(iy,1,j,id+1) = force(iy,1,j,id+1) + work(1)
          force(iy,2,j,id+1) = force(iy,2,j,id+1) + work(2)
          force(iy,3,j,id+1) = force(iy,3,j,id+1) + work(3)
!
        end do


!       if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
!        eelec(id+1) = eelec(id+1) + elec_temp
!        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

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


    return

end subroutine kernel
