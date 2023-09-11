c     #####################################################
c     ##                                                 ##
c     ##  module ddx -- workspaces, parameters and RHSs  ##
c     ##                                                 ##
c     #####################################################
c
      module ddx_interface
c     external module from libddx.so
      use ddx_core , only : ddx_type, ddx_state_type, deallocate_model,
     &   deallocate_state, allocate_state, get_banner,
     &   closest_supported_lebedev_grid, ddx_electrostatics_type
      use ddx_multipolar_solutes , only: multipole_electrostatics,
     &   multipole_force_terms, multipole_psi
      use ddx_errors , only : ddx_error_type, check_error
      use ddx , only: ddinit, ddrun
      implicit none
c
c     flags used by Tinker
c
      logical use_ddx
      character*8 ddx_typ
      logical save_cavity
      character*120 ddx_cavity_file
      logical printed_header
      character*8 cavity_typ
      real*8 sas_r, ddx_conv
c
c     data structures for ddx
c
      type(ddx_type) ddx_data
      type(ddx_state_type) ddx_state
      type(ddx_error_type) ddx_error
      type(ddx_electrostatics_type) ddx_electrostatics
c
c     control parameters for ddx
c
      real*8 eta, eps, se, kappa
      character*255 ddx_output_file
      integer lmax, ngrid, force, fmm, pm, pl, matvecmem,
     &   maxiter, jacobi_ndiis, ddx_nproc, model
c
c     information about the multipoles:
c     mmax            maximum angular momentum of the multipoles
c     nm              number of multipoles
c     ddx_mcoords     coordinates of the multipoles, size (3,nm)
c     ddx_multipoles  multipoles in real spherical harmonics,
c                     size ((mmax+1)**2,nm)
c
      integer mmax, nm
      real*8, allocatable :: ddx_mcoords(:,:), ddx_multipoles(:,:)
c
c     information about the cavity (union of spheres):
c     nsph            number of spheres
c     ddx_coords      coordinates of the centers (3,nsph)
c     ddx_radii       radii of the spheres (nsph)
c
      integer nsph
      real*8, allocatable :: ddx_radii(:), ddx_coords(:,:)
c
c     RHSs for ddX: Phi (potential) and Psi
c
      real*8, allocatable :: phi_cav(:), psi(:,:), e_cav(:,:),
     $   g_cav(:,:,:)
c
c     UFF diameters in A for atoms from H to Kr.
c
      real*8 rii(36)
      save rii
      data rii/   2.886d+00,2.362d+00,2.451d+00,2.745d+00,4.083d+00,
     $  3.851d+00,3.660d+00,3.500d+00,3.364d+00,3.243d+00,2.983d+00,
     $  3.021d+00,4.499d+00,4.295d+00,4.147d+00,4.035d+00,3.947d+00,
     $  3.868d+00,3.812d+00,3.399d+00,3.295d+00,3.175d+00,3.144d+00,
     $  3.023d+00,2.961d+00,2.912d+00,2.872d+00,2.834d+00,3.495d+00,
     $  2.763d+00,4.383d+00,4.280d+00,4.230d+00,4.205d+00,4.189d+00,
     $  4.141d+00/
c
c     additional spheres
c
      integer nsph_add
      integer, parameter :: max_nsph_add = 100
      real*8 ddx_add_coords(3,max_nsph_add), ddx_add_radii(max_nsph_add)
      end module ddx_interface
