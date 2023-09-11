c
c     ###############################################################
c     ##                                                           ##
c     ## subroutine ddx_build_multipoles -- ddX initial multipoles ##
c     ##                                                           ##
c     ###############################################################
c
c     "ddx_build_multipoles" build the initial arrays of multipoles
c     and of their coordinates based on the forcefield and on the
c     geometry of the system
c
      subroutine ddx_build_multipoles
      use atoms
      use charge
      use ddx_interface
      use iounit
      use mpole
      use units
      use math
      implicit none
      integer :: i, ii, info
      real*8 sqrt3, fac0, fac1, fac2
      sqrt3 = sqrt(3.0d0)
      fac0 = 1.0d0/sqrt(4.0d0*pi)
      fac1 = 1.0d0/(sqrt(4.0d0*pi/3.0)*bohr)
      fac2 = 1.0d0/(sqrt(4.0d0*pi/5.0)*bohr*bohr)
c
c     allocate the multipole arrays
c
      nm = n
      if (nion.ne.0) then
         mmax = 0
      end if
      if (npole.ne.0) then
         mmax = 2
      end if
      allocate(ddx_multipoles((mmax+1)**2,nm),ddx_mcoords(3,nm),
     &   stat=info)
      if (info.ne.0) then
         write(iout,*) 'Allocation failed in ddx_build_multipoles'
         call fatal
      end if
      ddx_multipoles = 0.0d0
c
c     fill the array of multipoles
c
      if (nion.ne.0) then
         do i = 1, nion
            ii = iion(i)
            if (ii .eq. 0) cycle
            ddx_multipoles(1,ii) = ddx_multipoles(1,ii) + fac0*pchg(i)
         end do
      else if (npole.ne.0) then
         do i = 1, npole
            ii = pollist(i)
            if (ii .eq. 0) cycle
            ddx_multipoles(1,ii) = ddx_multipoles(1,ii)
     &         + fac0*rpole(1,i)
            ddx_multipoles(2,ii) = ddx_multipoles(2,ii)
     &         + fac1*rpole(3,i)
            ddx_multipoles(3,ii) = ddx_multipoles(3,ii)
     &         + fac1*rpole(4,i)
            ddx_multipoles(4,ii) = ddx_multipoles(4,ii)
     &         + fac1*rpole(2,i)
            ddx_multipoles(5,ii) = ddx_multipoles(5,ii)
     &         + fac2*rpole(6,i)*sqrt3
            ddx_multipoles(6,ii) = ddx_multipoles(6,ii)
     &         + fac2*rpole(10,i)*sqrt3
            ddx_multipoles(7,ii) = ddx_multipoles(7,ii)
     &         + fac2*rpole(13,i)*3.0d0/2.0d0
            ddx_multipoles(8,ii) = ddx_multipoles(8,ii)
     &         + fac2*rpole(7,i)*sqrt3
            ddx_multipoles(9,ii) = ddx_multipoles(9,ii)
     &         +fac2*(rpole(5,i)-rpole(9,i))*sqrt3/2.0d0
         end do
      end if
c
c     copy and convert the coordinates
c
      ddx_mcoords(1,:) = x(:n)/bohr
      ddx_mcoords(2,:) = y(:n)/bohr
      ddx_mcoords(3,:) = z(:n)/bohr
c
      end subroutine ddx_build_multipoles
