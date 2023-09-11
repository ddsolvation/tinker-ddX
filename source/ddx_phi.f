      subroutine ddx_phi
      use atoms
      use ddx_interface
      use sizes
      use charge
      use mpole
      implicit none
      integer its, i, ii, info
      real*8  zero, pt5, one, three
      real*8  dx, dy, dz, dd, tmp, r, r2, r3, r5, ci, dix, diy, diz,
     $   qixx, qixy, qixz, qiyy, qiyz, qizz, qix, qiy, qiz, scd, scq,
     $   fc, fd, fq
      save zero, pt5, one, three
      data zero/0.0d0/, pt5/0.5d0/, one/1.0d0/, three/3.0d0/
c
c     clean RHS
c
      phi_cav = 0.0d0
      psi = 0.0d0
      if (nion.ne.0) then
c
c        charges only.
c
      do its = 1, ddx_data%constants%ncav
            do i = 1, nion
               ii = iion(i)
               dx =  ddx_data%constants%ccav(1,its) - x(ii)
               dy =  ddx_data%constants%ccav(2,its) - y(ii)
               dz =  ddx_data%constants%ccav(3,its) - z(ii)
               r = sqrt(dx*dx + dy*dy + dz*dz)
               phi_cav(its) = phi_cav(its) + pchg(i)/r
            end do
         end do
      else if (npole.ne.0) then
c
c        multipoles up to quadrupole
c
         do its = 1, ddx_data%constants%ncav
            tmp = zero
            do i = 1, npole
               ii = ipole(i)
               dx = ddx_data%constants%ccav(1,its) - x(ii)
               dy = ddx_data%constants%ccav(2,its) - y(ii)
               dz = ddx_data%constants%ccav(3,its) - z(ii)
               r2 = dx*dx + dy*dy + dz*dz
               r  = sqrt(r2)
               r3 = r*r2
               r5 = r3*r2
c
c              get the permanent multipoles
c
               ci   = rpole(1,i)
               dix  = rpole(2,i)
               diy  = rpole(3,i)
               diz  = rpole(4,i)
               qixx = rpole(5,i)
               qixy = rpole(6,i)
               qixz = rpole(7,i)
               qiyy = rpole(9,i)
               qiyz = rpole(10,i)
               qizz = rpole(13,i)
c
c              intermediates:
c
               qix = qixx*dx + qixy*dy + qixz*dz
               qiy = qixy*dx + qiyy*dy + qiyz*dz
               qiz = qixz*dx + qiyz*dy + qizz*dz
               scd = dix*dx + diy*dy + diz*dz
               scq = qix*dx + qiy*dy + qiz*dz
               fc  = one/r
               fd  = one/r3
               fq  = pt5*three/r5
               tmp = tmp + fc*ci + fd*scd + fq*scq
            end do
            phi_cav(its) = tmp
         end do
      end if
      end subroutine ddx_phi
