module gparameter

  implicit none
  private

  integer,  public, parameter :: dp = selected_real_kind(15, 307)
  integer,  public, parameter :: sp = selected_real_kind(6, 37)
  integer,  public, parameter :: wp = sp
  integer, public, parameter  :: nxyz= 3
  integer, public, parameter  :: ncq= 4
! apoa1
  real(wp), public, parameter :: expected_result(1:3)=[-157.941620,337.906647,-217.176392]



  type,public  :: s_genesis_kernel_param

    integer                     :: maxcell
    integer                     :: ncell
    integer                     :: ncell_local
    integer                     :: num_atom_domain
    integer                     :: MaxAtom
    integer                     :: MaxNb15
    integer                     :: MaxNb15_Fugaku
    integer                     :: nthread
    integer                     :: cutoff_int2
    integer                     :: max_class
    real(wp)                    :: density
    real(wp)                    :: cutoff
    real(wp)                    :: cutoff2
    real(dp)                    :: time
    integer,allocatable         :: natom(:)
    integer,allocatable         :: num_nb15_calc1(:,:)
    integer,allocatable         :: num_nb15_calc(:,:)
    integer,allocatable         :: nb15_calc_list1(:,:)
    integer,allocatable         :: nb15_calc_list(:,:)
    integer,allocatable         :: nb15_list(:,:)
    integer,allocatable         :: nb15_cell(:)
    integer,allocatable         :: atmcls(:,:)
    integer,allocatable         :: cell_pairlist(:,:)
    integer,allocatable         :: cell_move(:,:,:)
    real(wp),allocatable        :: charge(:,:)
    real(wp),allocatable        :: coord(:,:,:)
    real(wp),allocatable        :: coord_pbc(:,:,:)
    real(wp),allocatable        :: system_size(:)
    real(wp),allocatable        :: force(:,:,:,:)
    real(wp),allocatable        :: trans1(:,:,:)
    real(wp),allocatable        :: nonb_lj6(:,:)
    real(wp),allocatable        :: nonb_lj12(:,:)
    real(wp),allocatable        :: table_grad(:)
    real(wp),allocatable        :: virial(:,:)

  end type s_genesis_kernel_param

end module gparameter
module module_pointers
use gparameter
contains


subroutine read_data_file (gparam)
    type(s_genesis_kernel_param), intent(inout) :: gparam
    integer   :: alloc_stat
    integer   :: maxcell, ncell_local, ncell,                        &
                 num_atom_domain, MaxNb15_Fugaku,MaxNb15,            &
                 MaxAtom, nthread, cutoff_int2, &
                 max_class

#ifdef OMP_PARALLEL
    integer   :: omp_get_num_threads 
#endif

    ! local variables
    real(dp)                  :: Val

    character*30 :: filename
    integer :: iunit, is_ok
    iunit=77
    write(filename,'(a)') "../data/data_kernel_generic"

    open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
    if (is_ok.ne.0) then
        write(*,'(a,a)') "*** Error. failed to open file: ", filename
    endif
#ifdef OMP_PARALLEL
    !$omp parallel shared(nthread)
    gparam%nthread = omp_get_num_threads()
    !$omp end parallel
#else
    gparam%nthread = 1
#endif
    call sub_i4_data_read(iunit, gparam%maxcell,          1 )
    call sub_i4_data_read(iunit, gparam%ncell,            1 )
    call sub_i4_data_read(iunit, gparam%ncell_local,      1 )
    call sub_i4_data_read(iunit, gparam%MaxAtom,          1 )
    call sub_i4_data_read(iunit, gparam%MaxNb15,          1 )
    call sub_i4_data_read(iunit, gparam%max_class,        1 )
    call sub_i4_data_read(iunit, gparam%cutoff_int2,      1 )
    call sub_i4_data_read(iunit, gparam%num_atom_domain,  1 )

    alloc_stat = 0
    maxcell         = gparam%maxcell
    ncell           = gparam%ncell
    ncell_local     = gparam%ncell_local
    MaxAtom         = gparam%MaxAtom
    MaxNb15         = gparam%MaxNb15
    cutoff_int2     = gparam%cutoff_int2
    nthread         = gparam%nthread
    num_atom_domain = gparam%num_atom_domain
    max_class       = gparam%max_class

    allocate(gparam%natom(1:ncell),                           &
             gparam%charge(1:MaxAtom,1:ncell),                &
             gparam%atmcls(1:MaxAtom,1:ncell),                &
             gparam%nonb_lj6(1:max_class,1:max_class),        &
             gparam%nonb_lj12(1:max_class,1:max_class),       &
             gparam%num_nb15_calc1(1:MaxAtom,1:ncell_local),         &
             gparam%num_nb15_calc(1:MaxAtom,1:maxcell),         &
             gparam%nb15_calc_list1(1:MaxNb15,1:ncell_local),  &
             gparam%nb15_calc_list(1:MaxNb15,1:maxcell),      &
             gparam%nb15_cell(1:maxcell),                     &
             gparam%nb15_list(1:MaxAtom,1:maxcell),           &
             gparam%cell_pairlist(1:2,1:maxcell),             &
             gparam%coord(1:nxyz,1:MaxAtom,1:ncell),          &
             gparam%coord_pbc(1:MaxAtom,1:ncq,1:ncell),       &
             gparam%force(1:MaxAtom,1:nxyz,1:ncell,1:nthread),     &
             gparam%trans1(1:nxyz,1:MaxAtom,1:ncell),         &
             gparam%cell_move(1:nxyz,1:ncell,1:ncell),         &
             gparam%table_grad(1:cutoff_int2*6),              &
             gparam%virial(1:nxyz,1:maxcell),          &
             gparam%system_size(1:nxyz),                      &
             stat = alloc_stat)


    call sub_i4_data_read(iunit, gparam%natom        , ncell )
    call sub_R4_data_read(iunit, gparam%charge       , (MaxAtom*ncell) )
    call sub_i4_data_read(iunit, gparam%atmcls       , (MaxAtom*ncell) )
    call sub_R4_data_read(iunit, gparam%nonb_lj6     , (max_class*max_class) )
    call sub_R4_data_read(iunit, gparam%nonb_lj12    , (max_class*max_class) )
    call sub_i4_data_read(iunit, gparam%num_nb15_calc1  , (MaxAtom*ncell_local) )
    call sub_i4_data_read(iunit, gparam%num_nb15_calc  , (MaxAtom*maxcell) )
    call sub_i4_data_read(iunit, gparam%nb15_calc_list1 , (MaxNb15*ncell_local))
    call sub_i4_data_read(iunit, gparam%nb15_calc_list , (MaxNb15*maxcell) )
    call sub_i4_data_read(iunit, gparam%nb15_cell , (maxcell) )
    call sub_i4_data_read(iunit, gparam%nb15_list , (MaxAtom*maxcell) )
    call sub_i4_data_read(iunit, gparam%cell_pairlist, (2*maxcell) )
    call sub_R4_data_read(iunit, gparam%coord        , (nxyz*MaxAtom*ncell) )
    call sub_R4_data_read(iunit, gparam%trans1       , (nxyz*MaxAtom*ncell) )
    call sub_i4_data_read(iunit, gparam%cell_move    , (nxyz*ncell*ncell) )
    call sub_R4_data_read(iunit, gparam%system_size    , (nxyz) )
    call sub_r4_data_read(iunit, gparam%table_grad   , (cutoff_int2*6) )
    call sub_r4_data_read(iunit, gparam%density      , 1)
    call sub_r4_data_read(iunit, gparam%cutoff       , 1)
    call sub_r4_data_read(iunit, gparam%cutoff2      , 1)

    close(iunit)

    write(*,*) "maxcell=", maxcell, "   loc(maxcell)=", loc(maxcell)

end subroutine


subroutine check_validation (MaxAtom, ncell, nthread, force)
    use gparameter
    real(wp) :: force(1:MaxAtom,1:nxyz,1:ncell,1:nthread)
    real(wp) :: val(1:3), check_result
    integer  :: j, k, l,ll

    val(1:3)=0.0d0
    do l=1,nthread
      do i=1,3
        do ll = 1, ncell
          do k=1,MaxAtom,5   ! sum of all forces should be zero, then "5" is required
            val(i)=val(i)+real(force(k,i,ll,l),wp)
          enddo
        enddo
      enddo
    enddo

!    expected_result(:) = 0.0
    do i=1,3
      check_result = val(i) / expected_result(i)

      if ( 0.999 < check_result .and. check_result < 1.001 ) then
          write(*,*) "val=", val(i), "    expected_result=", expected_result(i)
          write(*,*) "the computed result seems to be OK."
      else
          write(*,*) "val=", val(i), "    expected_result=", expected_result(i)
          write(*,*) "the computed result is not close enough to the expected value."
      endif
    end do

end subroutine

end module

