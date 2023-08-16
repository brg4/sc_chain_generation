module quaternions

  use params

  implicit none

contains

  pure function cross_prod(va,vb) result(v_c)

    implicit none

    real(rp), intent(in) :: va(1:3), vb(1:3)

    real(rp) :: v_c(1:3)

    v_c(1) = va(2)*vb(3)-va(3)*vb(2)
    v_c(2) = va(3)*vb(1)-va(1)*vb(3)
    v_c(3) = va(1)*vb(2)-va(2)*vb(1)

  end function cross_prod

  pure function q_conj(q) result(q_c)

    implicit none

    real(rp), intent(in) :: q(1:4)
    real(rp) :: q_c(1:4)

    q_c(1) = q(1)
    q_c(2:4) = -q(2:4)
    
  end function q_conj

  pure function q_norm(q) result(m)

    implicit none

    real(rp), intent(in) :: q(1:4)
    real(rp) :: m

    m = norm2(q)
    
  end function q_norm

  pure function q_inv(q) result(q_i)

    implicit none

    real(rp), intent(in) :: q(1:4)
    real(rp) :: q_i(1:4)

    q_i = q_conj(q)/dot_product(q,q)
    
  end function q_inv

  pure function q_mult(qa,qb) result(q_m)

    implicit none

    real(rp), intent(in) :: qa(1:4), qb(1:4)

    real(rp) :: q_m(1:4)

    q_m(1) = qa(1)*qb(1) - dot_product(qa(2:4),qb(2:4))

    q_m(2:4) = qa(1)*qb(2:4) + qb(1)*qa(2:4) + cross_prod(qa(2:4),qb(2:4))

  end function q_mult

  pure function q_unit(q) result(q_u)

    implicit none

    real(rp), intent(in) :: q(1:4)

    real(rp) :: q_u(1:4)

    q_u = q/q_norm(q)

  end function q_unit

  pure function Rmat_to_q(R) result(q)

    implicit none

    real(rp), intent(in) :: R(1:3,1:3)

    real(rp) :: q(1:4), tr, S

    q = 0.0
    
    tr = R(1,1) + R(2,2) + R(3,3)

    if (tr .gt. 0.0d0) then
       
       S = 2.0*sqrt(1.0+tr)
       
       q(1) = 0.25*S
       q(2) = (R(3,2)-R(2,3))/S
       q(3) = (R(1,3)-R(3,1))/S
       q(4) = (R(2,1)-R(1,2))/S
       
    elseif ((R(1,1) .gt. R(2,2)) .and. (R(1,1) .gt. R(3,3))) then

       S = 2.0*sqrt(1.0+R(1,1)-R(2,2)-R(3,3))

       q(1) = (R(3,2)-R(2,3))/S
       q(2) = 0.25*S
       q(3) = (R(1,2)+R(2,1))/S
       q(4) = (R(1,3)+R(3,1))/S
       
    elseif (R(2,2) .gt. R(3,3)) then
       
       S = 2.0*sqrt(1.0+R(2,2)-R(1,1)-R(3,3))

       q(1) = (R(1,3)-R(3,1))/S
       q(2) = (R(1,2)+R(2,1))/S
       q(3) = 0.25*S
       q(4) = (R(2,3)+R(3,2))/S

    else
       
       S = 2.0*sqrt(1.0+R(3,3)-R(1,1)-R(2,2))

       q(1) = (R(2,1)-R(1,2))/S
       q(2) = (R(1,3)+R(3,1))/S
       q(3) = (R(2,3)+R(3,2))/S
       q(4) = 0.25*S
       
    end if

  end function Rmat_to_q

end module quaternions
