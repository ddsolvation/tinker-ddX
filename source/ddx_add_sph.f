c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine ddx_add_sph -- include additional spheres ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "ddx_add_sph" copies the additional spheres read from the input
c     file to the appropriate arrays
c
      subroutine ddx_add_sph
      use ddx_interface
      implicit none
      real*8, allocatable :: tmp_coords(:,:), tmp_radii(:)
      integer isph, info
c
c     check if there are additional spheres
c
      if (nsph_add.le.0) return
c
c     backup the existing coordinates
c
      allocate(tmp_coords(3,nsph),tmp_radii(nsph),stat=info)
      if (info.ne.0) then
         write(6,*) 'Reallocation failed in ddx_add_sph'
         call fatal
      end if
      tmp_coords = ddx_coords
      tmp_radii = ddx_radii
      deallocate(ddx_coords,ddx_radii,stat=info)
      if (info.ne.0) then
         write(6,*) 'Reallocation failed in ddx_add_sph'
         call fatal
      end if
c
c     allocate the new arrays
c
      allocate(ddx_coords(3,nsph+nsph_add),ddx_radii(nsph+nsph_add),
     &   stat=info)
      if (info.ne.0) then
         write(6,*) 'Reallocation failed in ddx_add_sph'
         call fatal
      end if
c
c     copy the data
c
      ddx_coords(:,:nsph) = tmp_coords
      ddx_radii(:nsph) = tmp_radii
      do isph = 1, nsph_add
         ddx_coords(:,nsph+isph) = ddx_add_coords(:,isph)
         ddx_radii(nsph+isph) = ddx_add_radii(isph)
      end do
      nsph = nsph + nsph_add
c
c     deallocate the temporaries
c
      deallocate(tmp_coords,tmp_radii,stat=info)
      if (info.ne.0) then
         write(6,*) 'Reallocation failed in ddx_add_sph'
         call fatal
      end if
c
      end subroutine ddx_add_sph
