c
c     #################################################
c     ##                                             ##
c     ##  subroutine ddx_init -- ddX initialization  ##
c     ##                                             ##
c     #################################################
c
c     "ddx_init" initialize the ddx data structures and precomputes
c     all the necessary data. Almost everything depends on the geometry
c     of the system, so this routine must be called at every step.
c
      subroutine ddx_init
      use atoms
      use ddx_interface
      use iounit
      use math
      use units
      implicit none
      integer info
      character (len=2047) :: banner
c
c     call the ddx initialization routine
c
      call ddinit(model,nsph,ddx_coords,ddx_radii,eps,ddx_data,
     &   ddx_error,force=force,kappa=kappa,eta=eta,shift=se,lmax=lmax,
     &   ngrid=ngrid,incore=matvecmem,maxiter=maxiter,
     &   jacobi_ndiis=jacobi_ndiis,enable_fmm=fmm,pm=pm,pl=pl,
     &   nproc=ddx_nproc,logfile=ddx_output_file)
c
c     check if initialization was successful
c
      call check_error(ddx_error)
c
c     print the ddx header
c
      if (.not.printed_header) then
         call get_banner(banner)
         write(iout,*) trim(banner)
         printed_header = .true.
      end if
c
c     if required save the cavity
c
      if (save_cavity) call ddx_save_geom
c
c     init also the object that contains the solutions
c
      call init_state(ddx_data%params,ddx_data%constants,ddx_state,
     &   ddx_error)
c
c     check if initialization was successful
c
      call check_error(ddx_error)
c
      end subroutine ddx_init
