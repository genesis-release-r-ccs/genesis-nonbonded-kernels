module gparameter

  implicit none
  private

  integer,  public, parameter :: dp = selected_real_kind(15, 307)
  integer,  public, parameter :: sp = selected_real_kind(6, 37)
  integer,  public, parameter :: wp = sp
  integer, public, parameter  :: nxyz= 3
  integer, public, parameter  :: ncq= 4
! lysozyme
!  real(wp), public, parameter :: expected_result(1:3)=[-699.488770, -318.171570,-273.786652]
 ! apoa1
   real(wp), public, parameter :: expected_result(1:3)=[-283.682770,315.114319,-51.6235008]

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
    integer,allocatable         :: start_atom(:)
    integer,allocatable         :: num_nb15_calc(:,:)
    integer,allocatable         :: nb15_calc_list(:,:,:)
    integer,allocatable         :: atmcls(:,:)
    integer,allocatable         :: atmcls_pbc(:)
    real(wp),allocatable        :: charge(:,:)
    real(wp),allocatable        :: coord(:,:,:)
    real(wp),allocatable        :: coord_pbc(:,:,:)
    real(wp),allocatable        :: force(:,:,:,:)
    real(wp),allocatable        :: trans1(:,:,:)
    real(wp),allocatable        :: nonb_lj6(:,:)
    real(wp),allocatable        :: nonb_lj12(:,:)
    real(wp),allocatable        :: table_grad(:)
    real(wp),allocatable        :: virial(:,:,:)

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

    character*26 :: filename
    integer :: iunit, is_ok
    iunit=77
    write(filename,'(a)') "../data/data_kernel_Oct"

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
    call sub_i4_data_read(iunit, gparam%MaxNb15_Fugaku,   1 )
    call sub_i4_data_read(iunit, gparam%max_class,        1 )
    call sub_i4_data_read(iunit, gparam%cutoff_int2,      1 )
    call sub_i4_data_read(iunit, gparam%num_atom_domain,  1 )

    alloc_stat = 0
    maxcell         = gparam%maxcell
    ncell           = gparam%ncell
    ncell_local     = gparam%ncell_local
    MaxAtom         = gparam%MaxAtom
    MaxNb15         = gparam%MaxNb15
    MaxNb15_Fugaku  = gparam%MaxNb15_Fugaku
    cutoff_int2     = gparam%cutoff_int2
    nthread         = gparam%nthread
    num_atom_domain = gparam%num_atom_domain
    max_class       = gparam%max_class

    allocate(gparam%natom(1:ncell),                           &
             gparam%start_atom(1:ncell),                      &
             gparam%charge(1:MaxAtom,1:ncell),                &
             gparam%atmcls(1:MaxAtom,1:ncell),                &
             gparam%atmcls_pbc(1:MaxAtom*ncell),              &
             gparam%nonb_lj6(1:max_class,1:max_class),        &
             gparam%nonb_lj12(1:max_class,1:max_class),       &
             gparam%num_nb15_calc(1:MaxAtom,1:ncell),         &
             gparam%nb15_calc_list(1:MaxNb15_Fugaku,1:MaxAtom*ncell,1:1),  &
             gparam%coord(1:nxyz,1:MaxAtom,1:ncell),          &
             gparam%coord_pbc(1:MaxAtom*ncell,1:ncq,1:1),     &
             gparam%force(1:MaxAtom*ncell,1:nxyz,1:1,1:nthread),     &
             gparam%trans1(1:nxyz,1:MaxAtom,1:ncell),         &
             gparam%table_grad(1:cutoff_int2*6),              &
             gparam%virial(1:nxyz,1:nxyz,1:nthread),          &
             stat = alloc_stat)


    call sub_i4_data_read(iunit, gparam%natom        , ncell )
    call sub_R4_data_read(iunit, gparam%charge       , (MaxAtom*ncell) )
    call sub_i4_data_read(iunit, gparam%atmcls       , (MaxAtom*ncell) )
    call sub_R4_data_read(iunit, gparam%nonb_lj6     , (max_class*max_class) )
    call sub_R4_data_read(iunit, gparam%nonb_lj12    , (max_class*max_class) )
    call sub_i4_data_read(iunit, gparam%num_nb15_calc  , (MaxAtom*ncell) )
    call sub_i4_data_read(iunit, gparam%nb15_calc_list , (MaxNb15_Fugaku*MaxAtom*ncell) )
    call sub_i4_data_read(iunit, gparam%start_atom   , ncell )
    call sub_R4_data_read(iunit, gparam%coord        , (nxyz*MaxAtom*ncell) )
    call sub_R4_data_read(iunit, gparam%trans1       , (nxyz*MaxAtom*ncell) )
    call sub_r4_data_read(iunit, gparam%table_grad   , (cutoff_int2*6) )
    call sub_r4_data_read(iunit, gparam%density      , 1)
    call sub_r4_data_read(iunit, gparam%cutoff       , 1)
    call sub_r4_data_read(iunit, gparam%cutoff2      , 1)

    close(iunit)

    write(*,*) "maxcell=", maxcell, "   loc(maxcell)=", loc(maxcell)

end subroutine


subroutine check_validation (MaxAtom, ncell, nthread, force)
    use gparameter
    real(wp) :: force(1:MaxAtom*ncell,1:nxyz,1:1,1:nthread)
    real(wp) :: val(1:3), check_result
    integer  :: j, k, l

    val(1:3)=0.0d0
    do l=1,nthread
      do i=1,3
        do k=1,ncell*MaxAtom,5   ! sum of all forces should be zero, then "5" is required
          val(i)=val(i)+real(force(k,i,1,l),wp)
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

