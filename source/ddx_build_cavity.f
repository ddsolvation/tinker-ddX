c
c     #######################################################
c     ##                                                   ##
c     ## subroutine ddx_build_cavity -- ddX initial cavity ##
c     ##                                                   ##
c     #######################################################
c
c     "ddx_build_cavity" build the initial arrays of centers
c     and radii of the spheres composing the cavity
c
      subroutine ddx_build_cavity
      use atomid
      use atoms
      use ddx_interface
      use units
      implicit none
      integer i, info
c
c     allocate the cavity (array of radii and of coordinates)
c
      allocate(ddx_radii(n),ddx_coords(3,n),stat=info)
      if (info.ne.0) then
         write(6,*) 'Allocation in ddx_cavity failed'
         call fatal
      end if
c
c     note the table contains the diameters in A, we are converting to
c     radii in bohr
c
      do i = 1, n
         ddx_radii(i) = 0.5d0*rii(atomic(i))/bohr
      end do
c
c     modify according to cavity options
c
      if (cavity_typ(1:3) .eq. 'SAS') then
         ddx_radii = ddx_radii + sas_r/bohr
      else if (cavity_typ .eq. 'VDW') then
         ddx_radii = 1.1d0*ddx_radii
      else
         write(6,*) 'Cavity type not supported'
         call fatal
      end if
c
c     copy and convert to Bohr the coordinates
c
      nsph = n
      ddx_coords(1,:) = x(:n)/bohr
      ddx_coords(2,:) = y(:n)/bohr
      ddx_coords(3,:) = z(:n)/bohr
      centered_on_m = .true.
c
      end subroutine ddx_build_cavity
