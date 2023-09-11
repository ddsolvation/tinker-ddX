      subroutine ddx_free
      use ddx_interface
      use iounit
      implicit none
      integer info
c
c     deallocate solutions and temporary arrays
c
      call deallocate_state(ddx_state, ddx_error)
c
c     deallocate ddx_data object
c
      call deallocate_model(ddx_data, ddx_error)
c
c     check if everything is deallocated
c
      call check_error(ddx_error)
c
c     dellocate psi
c
      deallocate(psi,stat=info)
      if (info .ne. 0) then
         write(iout,*) 'Deallocation in ddx_free failed'
         call fatal
      end if
c
c
      end subroutine ddx_free
