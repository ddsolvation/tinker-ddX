c
c     #####################################################
c     ##                                                 ##
c     ## subroutine ddx_run -- compute energy and forces ##
c     ##                                                 ##
c     #####################################################
c
c     "ddx_run" compute the solvation energy and, if requested, the
C     correspondent forces using the specified model
c
      subroutine ddx_run(do_forces)
      use atoms
      use charge
      use ddx_interface
      use deriv
      use energi
      use iounit
      use inform
      use units
      implicit none
      logical :: do_forces
      real*8 fac, wall, cpu, total_force(3)
      integer i, info
      real*8, allocatable :: forces(:,:)
   10 format(X,A34,F20.10)
   20 format(X,A34,I20)
c
      if (do_forces) then
         force = 1
         allocate(forces(3, n),stat=info)
      else
         force = 0
      end if
c
c     assemble the multipoles for ddx depending on the force field
c
      call ddx_build_multipoles
c
c     assemble the cavity
c
      call ddx_build_cavity
c
c     include the additional spheres if present
c
      call ddx_add_sph
c
c     Init the ddx data structures
c
      call settime
      call ddx_init
      call gettime(wall, cpu)
      if (verbose) write(iout,10) 'ddx_init time(s):', wall
c
c     build the RHS and the representation of the solute density
c
      call settime
      call multipole_electrostatics(ddx_data%params,ddx_data%constants,
     &   ddx_data%workspace,ddx_multipoles,mmax,ddx_electrostatics,
     &   ddx_error)
      call check_error(ddx_error)
      call gettime(wall, cpu)
      if (verbose) write(iout,10) 'build_phi time(s):', wall
c
      call settime
      call multipole_psi(ddx_data%params,ddx_multipoles,mmax,psi)
      call gettime(wall, cpu)
      if (verbose) write(iout,10) 'build_psi time(s):', wall
c
c     run ddX
c
      call ddrun(ddx_data,ddx_state,ddx_electrostatics,psi,ddx_conv,es,
     &   ddx_error,forces,.false.)
      call check_error(ddx_error)
c
c     convert the energy to kcal/mol, apply the cosmo scaling factor
c
      es = es*hartree
      if (model .eq. 1) es = es*(eps - 1.0d0)/eps
      if (verbose) write(iout,10) 'Solvation energy (Hartree):',
     &   es/hartree
c
c     contributions coming from the derivatives of the electrostatic
c     potential phi, this is a helper routine from the module ddx_rhs
c     which is specific of a multipolar distribution of any order
c
      if (do_forces) then
         call settime
         call multipole_force_terms(ddx_data%params,ddx_data%constants,
     &      ddx_data%workspace,ddx_state,mmax,ddx_multipoles,forces,
     &      ddx_error)
         call check_error(ddx_error)
         call gettime(wall, cpu)
         if (verbose) write(iout,10) 'forces time(s):', wall
c
c        convert the forces to kcal/mol/A, apply the ddcosmo scaling
c        factor
c
         if (model .eq. 1) forces = forces*(eps - 1.0d0)/eps
         des = des + forces*hartree/bohr
c
c        sanity check on the forces
c
         total_force = 0.0d0
         do i = 1, ddx_data%params%nsph
            total_force(:) = total_force(:) + forces(:,i)
         end do
         if ((abs(total_force(1)) + abs(total_force(2)) +
     &         abs(total_force(3))) .gt. 1d-5) then
            write(iout,*) "Warning: the total ddX force is summing to",
     &         total_force
         end if
      end if
c
c     print some performance information
c
      if (verbose) then
         if (model .eq. 1) then
            write(iout,10) 'ddcosmo_solve time (s):', ddx_state%xs_time
            write(iout,20) 'ddcosmo_solve iterations:',
     &         ddx_state%xs_niter
         else if (model .eq. 2) then
            write(iout,10) 'ddpcm_solve phieps time (s):',
     &         ddx_state%phieps_time
            write(iout,20) 'ddpcm_solve phieps iterations:',
     &         ddx_state%phieps_niter
            write(iout,10) 'ddpcm_solve x time (s):', ddx_state%xs_time
            write(iout,20) 'ddpcm_solve x iterations:',
     &         ddx_state%xs_niter
         else if (model .eq. 3) then
            write(iout,10) 'ddlpb_solve time (s):', ddx_state%x_lpb_time
            write(iout,20) 'ddlpb_solve iterations:',
     &         ddx_state%x_lpb_niter
         end if
      end if
c
      if (verbose .and. do_forces) then
         if (model .eq. 1) then
            write(iout,10) 'ddcosmo_adjoint time (s):',
     &         ddx_state%s_time
            write(iout,20) 'ddcosmo_adjoint iterations:',
     &         ddx_state%s_niter
         else if (model .eq. 2) then
            write(iout,10) 'ddpcm_solve_adjoint s time (s):',
     &         ddx_state%s_time
            write(iout,20) 'ddpcm_solve_adjoint s iterations:',
     &         ddx_state%s_niter
            write(iout,10) 'ddpcm_solve_adjoint y time (s):',
     &         ddx_state%y_time
            write(iout,20) 'ddpcm_solve_adjoint y iterations:',
     &         ddx_state%y_niter
         else if (model .eq. 3) then
            write(iout,10) 'ddlpb_solve_adjoint time (s):',
     &         ddx_state%x_adj_lpb_time
            write(iout,20) 'ddlpb_solve_adjoint iterations:',
     &         ddx_state%x_adj_lpb_niter
         end if
      end if
c
c     free the data structures
c
      if (allocated(forces)) then
         deallocate(forces,stat=info)
      end if
      call ddx_free
c
      end subroutine ddx_run
