
subroutine sub_i1_data_write(iunit, iarray, n)
        integer :: iunit, n
        integer(1) :: iarray(n)
        integer :: i
        write(iunit,'(i10)') n
        write(iunit,'(i10)') (iarray(i),i=1,n)
        return
end subroutine
subroutine sub_i4_data_write(iunit, iarray, n)
        integer :: iunit, n
        integer(4) :: iarray(n)
        integer :: i
        write(iunit,'(i10)') n
        write(iunit,'(i10)') (iarray(i),i=1,n)
        return
end subroutine
subroutine sub_r4_data_write(iunit, array, n)
        integer :: iunit, n
        real(4) :: array(n)
        integer :: i
        write(iunit,'(i10)') n
        write(iunit,'(1pe25.15)') (array(i),i=1,n)
        return
end subroutine


subroutine sub_i1_data_read(iunit, iarray, n)
        integer :: iunit, n
        integer(1) :: iarray(n)
        integer :: i, nr
        character(20) :: a
        read(iunit,'(a20,i10)') a,nr
        if (n.ne.nr) then
            write(*,*) "** Warning. The records do not match array size. n, nr=", n, nr
        endif
        read(iunit,'(i1)') (iarray(i),i=1,nr)
        return
end subroutine
subroutine sub_i4_data_read(iunit, iarray, n)
        integer :: iunit, n
        integer(4) :: iarray(n)
        integer :: i, nr
        character(20) :: a
        read(iunit,'(a20,i10)') a,nr
        if (n.ne.nr) then
            write(*,*) "** Warning. The records do not match array size. n, nr=", n, nr
        endif
        read(iunit,'(i10)') (iarray(i),i=1,nr)
        return
end subroutine
subroutine sub_r4_data_read(iunit, array, n)
        integer :: iunit, n
        real(4) :: array(n)
        character(20) :: a
        integer :: i, nr
        read(iunit,'(a20,i10)') a,nr
        if (n.ne.nr) then
            write(*,*) "** Warning. The records do not match array size. n, nr=", n, nr
        endif
        read(iunit,'(1pe25.15)') (array(i),i=1,nr)
        return
end subroutine

subroutine sub_r8_data_read(iunit, array, n)
        integer :: iunit, n
        real(8) :: array(n)
        character(20) :: a
        integer :: i, nr
        read(iunit,'(a20,i10)') a,nr
        if (n.ne.nr) then
            write(*,*) "** Warning. The records do not match array size. n, nr=", n, nr
        endif
        read(iunit,'(1pe25.15)') (array(i),i=1,nr)
        return
end subroutine

