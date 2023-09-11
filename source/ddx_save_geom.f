      subroutine ddx_save_geom
      use atomid
      use atoms
      use ddx_interface
      use units
      implicit none
      integer i, j, k, iunit
   10 format(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)
      iunit = 200
      open(unit=iunit,file=ddx_cavity_file)
      j = 1
      do i = 1, n
         write(iunit,10) 'ATOM  ',j,name(i),' ',' MM','A',1,' ',x(i),
     &     y(i),z(i),1.0d0,0.0d0,' ','  '
         j = j + 1
      end do

      do i = 1, ddx_data%constants%ncav
         write(iunit,10) 'ATOM  ',j,'HE',' ','CAV','A',2,' ',
     &     bohr*ddx_data%constants%ccav(1,i),
     &     bohr*ddx_data%constants%ccav(2,i),
     &     bohr*ddx_data%constants%ccav(3,i),1.0d0,0.0d0,' ','  '
         j = j + 1
      end do
      close(iunit)
      end subroutine ddx_save_geom
