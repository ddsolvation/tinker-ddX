      subroutine ddx_compute_phi
      use sizes
      use atoms
      use charge
      use ddx_interface
      use iounit
      use mpole
      use units
      implicit none
      integer isph
      integer its, i, ii
      real*8  zero, pt5, one, three
      real*8  dx, dy, dz, dd, tmp, r, r2, r3, r5, ci, dix, diy, diz, 
     $  qixx, qixy, qixz, qiyy, qiyz, qizz, qix, qiy, qiz, scd, scq, fc,
     $  fd, fq
      save zero, pt5, one, three
      data zero/0.0d0/, pt5/0.5d0/, one/1.0d0/, three/3.0d0/
c
c     clean up:
c
      do its = 1, ddx_data%constants%ncav
        phi_cav(its) = zero
      end do
c
      if (nion.ne.0) then 
c
c        charges only.
c
         do its = 1, ddx_data%constants%ncav
            tmp = zero
            do i = 1, nion
               ii = iion(i)
               dx = ddx_data%constants%ccav(1,its) - ddx_mcoords(1,ii)
               dy = ddx_data%constants%ccav(2,its) - ddx_mcoords(2,ii)
               dz = ddx_data%constants%ccav(3,its) - ddx_mcoords(3,ii)
               dd = sqrt(dx*dx + dy*dy + dz*dz)
               tmp = tmp + pchg(i)/dd
            end do
            phi_cav(its) = tmp
         end do
      else if (npole.ne.0) then
c
c        multipoles up to quadrupole
c
         do its = 1, ddx_data%constants%ncav
            tmp = zero
            do i = 1, npole
               ii = ipole(i)
               dx = ddx_data%constants%ccav(1,its) - ddx_mcoords(1,ii)
               dy = ddx_data%constants%ccav(2,its) - ddx_mcoords(2,ii)
               dz = ddx_data%constants%ccav(3,its) - ddx_mcoords(3,ii)
               r2 = one/(dx*dx + dy*dy + dz*dz)
               r  = sqrt(r2)
               r3 = r*r2
               r5 = r3*r2
c
c              get the permanent multipoles
c
               ci   = rpole(1,i)
               dix  = rpole(2,i)/bohr
               diy  = rpole(3,i)/bohr
               diz  = rpole(4,i)/bohr
               qixx = rpole(5,i)/(bohr*bohr)
               qixy = rpole(6,i)/(bohr*bohr)
               qixz = rpole(7,i)/(bohr*bohr)
               qiyy = rpole(9,i)/(bohr*bohr)
               qiyz = rpole(10,i)/(bohr*bohr)
               qizz = rpole(13,i)/(bohr*bohr)
               write(6,*) qixx, qixy, qixz, qiyy, qiyz, qizz
c
c              intermediates:
c
               qix = qixx*dx + qixy*dy + qixz*dz
               qiy = qixy*dx + qiyy*dy + qiyz*dz
               qiz = qixz*dx + qiyz*dy + qizz*dz
               scd = dix*dx + diy*dy + diz*dz
               scq = qix*dx + qiy*dy + qiz*dz
               fc  = one*r
               fd  = one*r3
               fq  = three*r5
               tmp = tmp + fc*ci + fd*scd + fq*scq
            end do
            phi_cav(its) = tmp
         end do
      end if
c

      end subroutine ddx_compute_phi
