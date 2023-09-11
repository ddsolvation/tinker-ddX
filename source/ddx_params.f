c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ddx_params -- ddX parameters initialization  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ddx_params" assigns default values to the model parameters or
c     reads them from the input. Most of the required initialization
c     depends on the coordinates and hence it is done together with
c     the energy and force computation
c
      subroutine ddx_params
      use ddx_interface
      use iounit
      use keys
      use openmp
      use units
      implicit none
      integer info, i, k, next
      real*8 sph_x, sph_y, sph_z, sph_r
      character*20 keyword
      character*20 value
      character*240 record
      character*240 string
c
c     set the model
c
      if (ddx_typ(1:5) .eq. 'COSMO') then
         model = 1
         se = -1.0d0
      else if (ddx_typ(1:3) .eq. 'PCM') then
         model = 2
         se = 0.0d0
      else if (ddx_typ(1:3) .eq. 'LPB') then
         model = 3
         se = 0.0d0
      else
         write(iout,*) 'ddX model not supported'
         call fatal
      end if
c
c     set some defaults
c
      ddx_output_file = ''
      lmax = 6
      ngrid = 302
      force = 0
      fmm = 0
      pm = -1
      pl = -1
      eta = 0.01d0
      eps = 78.0d0
      kappa = 0.104d0
      matvecmem = 0
      maxiter = 200
      jacobi_ndiis = 20
      cavity_typ = 'VDW'
      sas_r = 1.4
      ddx_conv = 1.0d-8
c
c     treat properly ddx_nproc and force as they must be set by Tinker
c
      ddx_nproc = nthread
c
c     parse the keys from the input file and possibly override
c     the default values
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:3) .eq. 'DDX') then
            call gettext (record,value,next)
            call upcase (value)
            if (value(1:4) .eq. 'LMAX') then
               call getnumb (record,k,next)
               lmax = k
            else if (value(1:5) .eq. 'NGRID') then
               call getnumb (record,k,next)
               ngrid = k
            else if (value(1:3) .eq. 'FMM') then
               call getnumb (record,k,next)
               fmm = k
            else if (value(1:2) .eq. 'PM') then
               call getnumb (record,k,next)
               pm = k
            else if (value(1:2) .eq. 'PL') then
               call getnumb (record,k,next)
               pl = k
            else if (value(1:3) .eq. 'ETA') then
               string = record(next:)
               read (string,*,err=10,end=10) eta
            else if (value(1:2) .eq. 'SE') then
               string = record(next:)
               write(iout,*) 'Use option `ddx se` only for debug'
               read (string,*,err=10,end=10) se
            else if (value(1:3) .eq. 'EPS') then
               string = record(next:)
               read (string,*,err=10,end=10) eps
            else if (value(1:5) .eq. 'KAPPA') then
               string = record(next:)
               read (string,*,err=10,end=10) kappa
            else if (value(1:9) .eq. 'MATVECMEM') then
               call getnumb (record,k,next)
               matvecmem = k
            else if (value(1:7) .eq. 'MAXITER') then
               call getnumb (record,k,next)
               maxiter = k
            else if (value(1:12) .eq. 'JACOBI_NDIIS') then
               call getnumb (record,k,next)
               jacobi_ndiis = k
            else if (value(1:5) .eq. 'NPROC') then
               call getnumb (record,k,next)
               write(iout,*) 'Use option `ddx nproc` only for debug'
               ddx_nproc = k
            else if (value(1:7) .eq. 'SAVECAV') then
               save_cavity = .true.
               call gettext (record,value,next)
               call lowcase(value)
               ddx_cavity_file = value
            else if (value(1:6) .eq. 'OUTPUT') then
                call gettext (record,value,next)
                call lowcase(value)
                ddx_output_file = value
            else if (value(1:10) .eq. 'ADDITIONAL') then
               string = record(next:)
               read (string,*,err=10,end=10) sph_x, sph_y, sph_z, sph_r
               nsph_add = nsph_add + 1
               if (nsph_add .gt. max_nsph_add) then
                  write(iout,*) 'Increase the maximum amount of addition
     &l spheres'
                  call fatal
               end if
               ddx_add_coords(1,nsph_add) = sph_x/bohr
               ddx_add_coords(2,nsph_add) = sph_y/bohr
               ddx_add_coords(3,nsph_add) = sph_z/bohr
               ddx_add_radii(nsph_add) = sph_r
            else if (value(1:6) .eq. 'CAVITY') then
               call gettext (record,value,next)
               call upcase (value)
               cavity_typ = value
               if (value(1:3) .eq. 'SAS') then
                  string = record(next:)
                  read (string,*,err=10,end=10) sas_r
               end if
            else if (value(1:4) .eq. 'CONV') then
               string = record(next:)
               read (string,*,err=10,end=10) ddx_conv
            end if
   10       continue
         end if
      end do
c
c     set the defaults for pl and pm depending on lmax
c
      if (pl .eq. -1) pl = lmax
      if (pm .eq. -1) pm = lmax
c
c     find the closest lebedev grid
c
      call closest_supported_lebedev_grid(ngrid)
c
      printed_header = .false.
c
      end subroutine ddx_params
