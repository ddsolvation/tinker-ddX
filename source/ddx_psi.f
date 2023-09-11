      subroutine ddx_psi
      use atoms
      use sizes
      use charge
      use ddx_interface
      use mpole
      use polar
      use math
      use units
      implicit none
      integer nrhs
      integer isph, irhs, l, ind, mm, m
      real*8  zero, one, two, three, four, fac, sphpol(9), sq3, ri,
     &   sqrtpi
      save zero, one, two, three, four
      data zero/0.0d0/, one/1.0d0/, two/2.0d0/, three/3.0d0/,
     &   four/4.0d0/
c
      psi = zero
      sqrtpi = sqrt(four*atan(one))
c
      if (nion.ne.0 .and. npole.eq.0) then
c
c     charges
c
         do isph = 1, n
            fac = two*sqrtpi
            psi(1,isph) = fac*pchg(isph)
         end do
      else if (npole.ne.0) then
c         do isph = 1, n
c            sphpol = zero
c            sq3  = sqrt(three)
cc
cc           transform the multipoles into the spherical representation
cc
c            i = pollist(isph)
c            sphpol(1) = rpole(1,i)
c            sphpol(2) = rpole(3,i) + mu(2,irhs,isph)
c            sphpol(3) = rpole(4,i) + mu(3,irhs,isph)
c            sphpol(4) = rpole(2,i) + mu(1,irhs,isph)
c            sphpol(5) = rpole(6,i)*sq3
c            sphpol(6) = rpole(10,i)*sq3
c            sphpol(7) = rpole(13,i)*three/two
c            sphpol(8) = rpole(7,i)*sq3
c            sphpol(9) = (rpole(5,i)-rpole(9,i))*sq3/two
c            ri = one
c            do l = 0, 2
c               fac = two*ri*sqrtpi/(sqrt(two*dble(l)+one))
c               ind = l*l + l + 1
c               do mm = 0, 2*l
c                  m = mm - l
c                  do irhs = 1, nrhs
c                     vlm(irhs,ind+m) = fac*sphpol(irhs,ind+m)
c                  end do
c               end do
c               ri = ri/rsph(isph)
c            end do
c         end do
      end if
      end subroutine ddx_psi
