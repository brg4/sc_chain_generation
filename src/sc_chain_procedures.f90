module sc_chain_procedures

  use params
  use sc_chain_object
  use quaternions

  implicit none

contains

  !function to measure closest distance between two lines
  pure function dist_line_to_line(v1,v2,w1,w2,&
       clamp_v1,clamp_v2,clamp_w1,clamp_w2) result(d)

    implicit none

    real(rp), intent(in) :: v1(1:3), v2(1:3), w1(1:3), w2(1:3)
    logical, intent(in) :: clamp_v1, clamp_v2, clamp_w1, clamp_w2

    logical :: c_v1, c_v2, c_w1, c_w2
    
    real(rp) :: m_dv, m_dw, dv(1:3), dw(1:3), de(1:3), udv(1:3), udw(1:3)
    real(rp) :: A(1:2,1:2), b(1:2), mu(1:2)
    real(rp) :: ev, ew, vw, vv, ww
    real(rp) :: r_v(1:3), r_w(1:3)
    
    real(rp) :: d, dot_temp

    ! calculate the difference vectors
    dv = v2 - v1
    dw = w2 - w1
    de = v1 - w1

    !calculate the magnitudes and unit vectors
    m_dv = norm2(dv)
    udv = dv/m_dv
    m_dw = norm2(dw)
    udw = dw/m_dw

    !calculate the dot products
    ev = dot_product(de,dv)
    ew = dot_product(de,dw)
    vw = dot_product(dv,dw)
    vv = dot_product(dv,dv)
    ww = dot_product(dw,dw)

    !create the linear system
    A(1,1) = vv
    A(1,2) = -vw
    A(2,1) = -vw
    A(2,2) = ww
    b(1) = -ev
    b(2) = ew

    !solve the linear system
    mu(1) = (A(2,2)*b(1)-A(1,2)*b(2))/(A(1,1)*A(2,2)-A(1,2)*A(2,1))
    mu(2) = (b(2)-A(2,1)*mu(1))/A(2,2)

    !apply the clamp conditions to line v
    c_v1 = .false.
    c_v2 = .false.
    if (clamp_v1 .and. (mu(1) .lt. 0.0d0)) then
       mu(1) = 0.0d0
       c_v1 = .true.
    end if
    if (clamp_v2 .and. (mu(1) .gt. 1.0d0)) then
       mu(1) = 1.0d0
       c_v2 = .true.
    end if

    !apply the clamp conditions to line w
    c_w1 = .false.
    c_w2 = .false.
    if (clamp_w1 .and. (mu(2) .lt. 0.0d0)) then
       mu(2) = 0.0d0
       c_w1 = .true.
    end if
    if (clamp_w2 .and. (mu(2) .gt. 1.0d0)) then
       mu(2) = 1.0d0
       c_w2 = .true.
    end if

    r_v = v1 + mu(1)*dv
    r_w = w1 + mu(2)*dw

    if (c_v1 .or. c_v2) then
       dot_temp = dot_product(udw, r_v-w1)/m_dw
       if (c_w1 .and. (dot_temp .lt. 0.0d0)) then
          dot_temp = 0.0d0
       elseif (c_w2 .and. (dot_temp .gt. 1.0d0)) then
          dot_temp = 1.0d0
       end if
       r_w = w1 + dot_temp*dw
    end if

    if (c_w1 .or. c_w2) then
       dot_temp = dot_product(udv, r_w-v1)/m_dv
       if (c_v1 .and. (dot_temp .lt. 0.0d0)) then
          dot_temp = 0.0d0
       elseif (c_v2 .and. (dot_temp .gt. 1.0d0)) then
          dot_temp = 1.0d0
       end if
       r_v = v1 + dot_temp*dv
    end if

    d = norm2(r_v - r_w)
    
  end function dist_line_to_line

  
  !function to measure closest distance between line and point
  pure function dist_line_to_point(v1, v2, w, clamp) result(d)

    implicit none

    logical, intent(in) :: clamp
    
    real(rp), intent(in) :: v1(1:3), v2(1:3), w(1:3)
    
    real(rp) :: r(1:3), dv(1:3), mu

    real(rp) :: d

    dv = v2 - v1

    mu = dot_product(w-v1,dv)/dot_product(dv,dv)

    if (clamp) then
       if (mu .gt. 1.0d0) then
          mu = 1.0d0
       elseif (mu .lt. 0.0d0) then
          mu = 0.0d0
       end if
    end if

    r = v1 + mu*dv - w

    d = norm2(r)
    
  end function dist_line_to_point

  !function to test if a line segment intersects a plane
  pure function line_intersect_plane(w1,w2,p1,p2,p3) result(intersect)

    implicit none

    real(rp), intent(in), dimension(1:3) :: w1, w2, p1, p2, p3

    logical :: intersect

    real(rp), dimension(1:3) :: w21, p21, p31, wp

    real(rp) :: t, u, v, det, tol

    intersect = .false.

    tol = 1.0d-5

    w21 = w2 - w1
    p21 = p2 - p1
    p31 = p3 - p1
    wp = w1 - p1

    det = dot_product(-w21,cross_prod(p21,p31))

    if (abs(det) .lt. tol) return

    t = dot_product(cross_prod(p21,p31),wp)/det
    if (.not. ((t .ge. 0.0) .and. (t .le. 1.0))) return
    
    u = dot_product(cross_prod(p31,-w21),wp)/det
    if (.not. ((u .ge. 0.0) .and. (u .le. 1.0))) return
    
    v = dot_product(cross_prod(-w21,p21),wp)/det
    if (.not. ((v .ge. 0.0) .and. (v .le. 1.0))) return

    if (u + v .lt. 1.0) then
       intersect = .true.
       return
    end if

  end function line_intersect_plane

  !function to test sc-sc clashes
  pure function sc_sc_clash(x_sc1, x_sc2, r_sc1, r_sc2) result(clash)

    implicit none

    real(rp), intent(in) :: x_sc1(1:3,1:2), x_sc2(1:3,1:2), r_sc1, r_sc2

    logical :: clash

    real(rp) :: d, test_d

    !assume no clashes occur
    clash = .false.

    !calculate minimum distance
    d = r_sc1 + r_sc2

    !calculate closest distance
    test_d = dist_line_to_line(x_sc1(:,1),x_sc1(:,2),&
         x_sc2(:,1),x_sc2(:,2),&
         .true.,.true.,&
         .true.,.true.)

    if (test_d .le. d) clash = .true.
    
  end function sc_sc_clash

  
  !function to test sc-s clashes
  pure function sc_s_clash(x_sc, x_s, r_sc, r_s) result(clash)

    implicit none

    real(rp), intent(in) :: x_sc(1:3,1:2), x_s(1:3), r_sc, r_s

    logical :: clash

    real(rp) :: d, test_d

    !assume no clashes occur
    clash = .false.

    !calculate minimum distance
    d = r_sc + r_s

    !calculate closest distance
    test_d = dist_line_to_point(x_sc(:,1),x_sc(:,2),&
         x_s,.true.)

    if (test_d .le. d) clash = .true.
    
  end function sc_s_clash

  
  !function to test s-s clashes
  pure function s_s_clash(x_s1, x_s2, r_s1, r_s2) result(clash)

    implicit none

    real(rp), intent(in) :: x_s1(1:3), x_s2(1:3), r_s1, r_s2

    logical :: clash

    real(rp) :: d

    !assume no clashes occur
    clash = .false.

    !calculate minimum distance
    d = r_s1 + r_s2

    if (norm2(x_s2-x_s1) .le. d) clash = .true.
    
  end function s_s_clash

  
  !function to test angles of segments
  pure function test_angle(x_sc1, x_sc2, cos_ang_thresh) result(clash)

    implicit none

    real(rp), intent(in) :: x_sc1(1:3,1:2), x_sc2(1:3,1:2), cos_ang_thresh
    real(rp) :: v(1:3), w(1:3)
    
    logical :: clash

    !assume no clashes occur
    clash = .false.

    v = x_sc1(:,2) - x_sc1(:,1)
    w = x_sc2(:,2) - x_sc2(:,1)

    !test if angle is insufficient
    if (dot_product(v,w)/(norm2(v)*norm2(w)) .lt. cos_ang_thresh) clash = .true.

  end function test_angle

  
  !function to initialize a randomly-rotated triangle in the xy-plane
  subroutine initialize_triangle(x,x0,L)

    implicit none

    integer :: i

    real(rp), intent(inout) :: x(1:3,1:3)
    real(rp), intent(in) :: x0(1:3), L
    real(rp) :: w0, w(1:3)

    x = 0.0d0

    w0 = r_rand(0.0d0,two_pi)
    w = (two_pi/3.0)*(/(i, i=0,2)/)

    x(1,:) = cos(w+w0)
    x(2,:) = sin(w+w0)

    x = (L/sqrt(3.0d0))*x

    do i=1,3
       x(:,i) = x(:,i) + x0
    end do
    
  end subroutine initialize_triangle

  
  !subroutine to initialize spherocylinder chain
  subroutine initialize_chain(sc_s,sc_p)

    implicit none

    type(sc_chain_state), intent(inout) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    logical :: initial_point, clash

    integer(ip) :: i, j
    integer(ip) :: i_chain, j_chain

    real(rp) :: L, R_b_sqrd, x0(1:3), x_sc_temp(1:3,1:2), sc2(1:3,1:2)

    !initialize the chain state
    call new_sc_chain_state(sc_s,sc_p)

    L = sc_p%g_p%stages(1,1)*(2*sc_p%r)
    R_b_sqrd = (sc_p%R_b-(L/sqrt(3.0d0))+sc_p%r_sc)**2.0d0

    do i_chain=1,sc_p%N_chains

       initial_point = .false.

       do while (initial_point .eqv. .false.)

          !sample an initial point
          do i = 1,3
             x0(i) = r_rand(-sc_p%R_b,sc_p%R_b)
          end do

          !test if initial point is within boundary
          if (dot_product(x0,x0) .lt. R_b_sqrd) then

             call initialize_triangle(sc_s%chain(i_chain)%x(1:3,1:3),x0,L)

             clash = .false.

             !test the new spherocylinder segments for clashes
             do i = 1,3

                x_sc_temp(:,1) = sc_s%chain(i_chain)%x(:,i)
                x_sc_temp(:,2) = sc_s%chain(i_chain)%x(:,mod(i,3)+1)

                !test for clashes with obstacles
                if ((clash .eqv. .false.) .and. (sc_p%obst%N .gt. 0)) then
                   do j = 1,sc_p%obst%N
                      if (sc_s_clash(x_sc_temp,sc_p%obst%x(:,j),sc_p%r_sc,sc_p%R_o)) then
                         clash = .true.
                         exit
                      end if
                   end do !end loop over obstacles
                end if

                !test for clashes with previously placed chains
                if ((clash .eqv. .false.) .and. (i_chain .gt. 1)) then
                   do j_chain = 1,(i_chain-1)

                      do j = 1,3
                         
                         sc2(:,1) = sc_s%chain(j_chain)%x(:,j)
                         sc2(:,2) = sc_s%chain(j_chain)%x(:,mod(j,3)+1)

                         if (sc_s_clash(x_sc_temp,sc_p%obst%x(:,j),sc_p%r_sc,sc_p%R_o)) then
                            clash = .true.
                            exit
                         end if

                      end do

                      if (clash .eqv. .true.) exit
                      
                   end do !end loop over previous chains
                end if

             end do

             if (clash .eqv. .false.) then
                initial_point = .true.
                sc_s%chain%N = 3
             end if

          end if

       end do !end while loop attempting initial triangle placements

    end do !end loop over chains

  end subroutine initialize_chain

  
  !subroutine to initialize random obstacles
  subroutine initialize_obstacles(sc_p)

    implicit none

    type(sc_system_params), intent(inout) :: sc_p

    logical :: valid
    
    integer(ip) :: i, j

    real(rp) :: R_b_sqrd

    R_b_sqrd = (sc_p%R_b - sc_p%R_o)**2.0d0

    sc_p%obst%N = sc_p%obst%N_max

    do i = 1,sc_p%obst%N

       valid = .false.
       
       do while (valid .eqv. .false.)
          
          do j = 1,3
             sc_p%obst%x(j,i) = r_rand(-sc_p%R_b,sc_p%R_b)
          end do

          if (dot_product(sc_p%obst%x(:,i),sc_p%obst%x(:,i)) .lt. R_b_sqrd) valid=.true.

       end do
       
    end do
    
  end subroutine initialize_obstacles

  
  !function to test clashes during spherocylinder segment insertion
  pure function test_inserted_segment_clash(sc_s,sc_p,i_chain,&
       i,ipp,x_new) result(clash)

    type(sc_chain_state), intent(in) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    integer(ip), intent(in) :: i_chain, i, ipp
    real(rp), intent(in) :: x_new(1:3)

    logical :: clash

    logical :: intersect_seg1, intersect_seg2, ang_1, ang_2, knot

    integer(ip) :: j_chain
    integer(ip) :: c, max_c
    integer(ip) :: j, k, kp

    real(rp) :: R_b_sqrd, cos_ang_thresh
    real(rp) :: sc_seg1(1:3,1:2), sc_seg2(1:3,1:2), sc2(1:3,1:2)

    clash = .false.

    R_b_sqrd = (sc_p%R_b - sc_p%r_sc)**2.0d0
    cos_ang_thresh = cos(2.0d0*third_pi)
 
    !test if new point is within boundary
    if (dot_product(x_new,x_new) .ge. R_b_sqrd) then
       clash = .true.
       return
    end if

    sc_seg1(:,1) = sc_s%chain(i_chain)%x(:,i)
    sc_seg1(:,2) = x_new

    sc_seg2(:,1) = x_new
    sc_seg2(:,2) = sc_s%chain(i_chain)%x(:,ipp)

    !test if spherocylinder segments intersect with the chain undergoing addition
    do j = 2,(sc_s%chain(i_chain)%N-2)

       k = circ_trans(i+j,sc_s%chain(i_chain)%N)
       kp = circ_trans(i+j+1,sc_s%chain(i_chain)%N)

       sc2(:,1) = sc_s%chain(i_chain)%x(:,k)
       sc2(:,2) = sc_s%chain(i_chain)%x(:,kp)

       !test for intersections with the first segment
       intersect_seg1 = sc_sc_clash(sc_seg1,sc2,&
            sc_p%r_sc,sc_p%r_sc)

       !test for intersections with the second segment
       intersect_seg2 = sc_sc_clash(sc_seg2,sc2,&
            sc_p%r_sc,sc_p%r_sc)

       if (intersect_seg1 .or. intersect_seg2) then
          clash = .true.
          return
       end if

    end do

    !test for bad angles and nearby intersections
    if (sc_s%chain(i_chain)%N > 4) then

       !test the angle between the first segment and the previous segment
       !test the clash between the second segment and the segment prior to the first

       k = circ_trans(i-1,sc_s%chain(i_chain)%N)
       kp = circ_trans(i,sc_s%chain(i_chain)%N)
       sc2(:,1) = sc_s%chain(i_chain)%x(:,k)
       sc2(:,2) = sc_s%chain(i_chain)%x(:,kp)

       intersect_seg2 = sc_sc_clash(sc_seg2,sc2,&
            sc_p%r_sc,sc_p%r_sc)

       ! ang_1 = test_angle(sc_seg1,cshift(sc2,1,dim=2),&
       !      cos_ang_thresh)

       ang_1 = test_angle(sc_seg1,sc2,&
            cos_ang_thresh)

       if (ang_1 .or. intersect_seg2) then
          clash = .true.
          return
       end if

       !test the angle between the second segment and the subsequent segment
       !test the clash between the first segment and the segment subsequent to the second

       k = circ_trans(ipp,sc_s%chain(i_chain)%N)
       kp = circ_trans(ipp+1,sc_s%chain(i_chain)%N)
       sc2(:,1) = sc_s%chain(i_chain)%x(:,k)
       sc2(:,2) = sc_s%chain(i_chain)%x(:,kp)

       intersect_seg1 = sc_sc_clash(sc_seg1,sc2,&
            sc_p%r_sc,sc_p%r_sc)

       ! ang_2 = test_angle(cshift(sc_seg2,1,dim=2),sc2,&
       !      cos_ang_thresh)

       ang_2 = test_angle(sc_seg2,sc2,&
            cos_ang_thresh)

       if (ang_2 .or. intersect_seg1) then
          clash = .true.
          return
       end if

    end if

    !test if spherocylinder segments intersect other chains in system
    if (sc_p%N_chains .gt. 1) then

       do j_chain = 1,sc_p%N_chains

          if (j_chain .ne. i_chain) then

             do k = 1,sc_s%chain(j_chain)%N

                kp = circ_trans(k+1,sc_s%chain(j_chain)%N)

                sc2(:,1) = sc_s%chain(j_chain)%x(:,k)
                sc2(:,2) = sc_s%chain(j_chain)%x(:,kp)

                !test for intersections with the first segment
                intersect_seg1 = sc_sc_clash(sc_seg1,sc2,&
                     sc_p%r_sc,sc_p%r_sc)

                !test for intersections with the second segment
                intersect_seg2 = sc_sc_clash(sc_seg2,sc2,&
                     sc_p%r_sc,sc_p%r_sc)

                if (intersect_seg1 .or. intersect_seg2) then
                   clash = .true.
                   return
                end if

             end do

          end if

       end do

    end if

    !test for knots
    if (sc_p%knot_prevent .eq. 1) then

       !test for self-knotting of chain undergoing addition
       if ((sc_s%chain(i_chain)%N > 5) .and.&
            (sc_p%g_p%stages(1,sc_s%stage) .ge. 3)) then

          do j = 3,(sc_s%chain(i_chain)%N-3)

             k = circ_trans(i+j,sc_s%chain(i_chain)%N)
             kp = circ_trans(i+j+1,sc_s%chain(i_chain)%N)

             sc2(:,1) = sc_s%chain(i_chain)%x(:,k)
             sc2(:,2) = sc_s%chain(i_chain)%x(:,kp)

             knot = line_intersect_plane(sc2(:,1),sc2(:,2),&
                  sc_seg1(:,1),sc_seg2(:,1),sc_seg2(:,2))

             if (knot) then
                clash = .true.
                return
             end if

          end do

       end if !end conditional for knots of single chain

       !test for knots with other chains
       if (sc_p%N_chains .gt. 1) then

          do j_chain = 1,sc_p%N_chains

             if (j_chain .ne. i_chain) then

                do k = 1,sc_s%chain(j_chain)%N

                   kp = circ_trans(k+1,sc_s%chain(j_chain)%N)

                   sc2(:,1) = sc_s%chain(j_chain)%x(:,k)
                   sc2(:,2) = sc_s%chain(j_chain)%x(:,kp)

                   knot = line_intersect_plane(sc2(:,1),sc2(:,2),&
                        sc_seg1(:,1),sc_seg2(:,1),sc_seg2(:,2))

                   if (knot) then
                      clash = .true.
                      return
                   end if

                end do

             end if

          end do

       end if !end tests for multiple chains

    end if

    !test for intersections with spherical obstacles
    if (sc_p%obst%N .gt. 0) then

       do j = 1,sc_p%obst%N

          intersect_seg1 = sc_s_clash(sc_seg1,sc_p%obst%x(:,j),&
               sc_p%r_sc,sc_p%R_o)

          intersect_seg2 = sc_s_clash(sc_seg2,sc_p%obst%x(:,j),&
               sc_p%r_sc,sc_p%R_o)

          if (intersect_seg1 .or. intersect_seg2) then
             clash = .true.
             return
          end if

       end do

    end if

    return

  end function test_inserted_segment_clash

  
  !subroutine to insert a spherocylinder segment at a random location
  subroutine insert_random_segment(sc_s,sc_p,i_chain,insert_fail)

    implicit none

    type(sc_chain_state), intent(inout) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    integer(ip), intent(in) :: i_chain

    logical, intent(out) :: insert_fail

    logical :: accepted, clash

    integer(ip) :: j_chain
    integer(ip) :: c, max_c, i, ipp
    integer(ip) :: j, k, kp

    real(rp) :: mid_dist, R_b_sqrd, cos_ang_thresh
    real(rp) :: t, p, c_t, s_t, c_p, s_p
    real(rp) :: dxu
    real(rp) :: dx(1:3), xm(1:3), u(1:3), x_new(1:3)
    real(rp) :: sc_seg1(1:3,1:2), sc_seg2(1:3,1:2), sc2(1:3,1:2)
    
    accepted = .false.

    max_c = 100
    c = 0

    mid_dist = (sqrt(3.0d0)/2.0d0)*(sc_p%g_p%stages(1,sc_s%stage)*(2*sc_p%r))

    do while ((accepted .eqv. .false.) .and. (c .lt. max_c))

       i = int_rand(1,sc_s%chain(i_chain)%N)
       ipp = circ_trans(i+1,sc_s%chain(i_chain)%N)

       dx = sc_s%chain(i_chain)%x(:,ipp)-sc_s%chain(i_chain)%x(:,i)
       xm = (sc_s%chain(i_chain)%x(:,ipp)+sc_s%chain(i_chain)%x(:,i))/2.0d0

       u = dx/norm2(dx)

       t = r_rand(0.0d0,two_pi)
       p = r_rand(0.0d0,pi)

       c_t = cos(t)
       s_t = sin(t)
       c_p = cos(p)
       s_p = sin(p)

       dx = (/c_p*s_t,s_p*s_t,c_t/)

       dxu = dot_product(dx,u)

       if (abs(dxu) .gt. 0.9) cycle

       dx = dx - dxu*u

       dx = dx/norm2(dx)

       x_new = xm + mid_dist*dx

       clash = test_inserted_segment_clash(sc_s,sc_p,i_chain,i,ipp,x_new)

       !cycle the main loop if there was a clash
       if (clash) then
          c = c + 1
       else
          accepted = .true.
       end if
       
    end do

    if (accepted) then

       if (i .eq. (sc_s%chain(i_chain)%N)) then
          sc_s%chain(i_chain)%x(:,sc_s%chain(i_chain)%N+1) = x_new
       else
          sc_s%chain(i_chain)%x(:,ipp+1:sc_s%chain(i_chain)%N+1) =&
               sc_s%chain(i_chain)%x(:,ipp:sc_s%chain(i_chain)%N)
          sc_s%chain(i_chain)%x(:,ipp) = x_new
       end if

       sc_s%chain(i_chain)%N = sc_s%chain(i_chain)%N + 1

       insert_fail = .false.
       
    else
       insert_fail = .true.
    end if

  end subroutine insert_random_segment

  
  !function to test clashes during spherocylinder segment insertion
  pure function test_inserted_sphere_clash(sc_s,sc_p,i_chain,&
       i,ipp,x_new) result(clash)

    type(sc_chain_state), intent(in) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    integer(ip), intent(in) :: i_chain, i, ipp
    real(rp), intent(in) :: x_new(1:3)

    logical :: clash

    logical :: intersect

    integer(ip) :: j_chain
    integer(ip) :: j, k

    real(rp) :: R_b_sqrd

    clash = .false.

    R_b_sqrd = (sc_p%R_b - sc_p%r_sc)**2.0d0
 
    !test if new point is within boundary
    if (dot_product(x_new,x_new) .ge. R_b_sqrd) then
       clash = .true.
       return
    end if

    !test if spheres intersect
    do j = 2,(sc_s%chain(i_chain)%N-1)

       k = circ_trans(i+j,sc_s%chain(i_chain)%N)

       intersect = s_s_clash(x_new,sc_s%chain(i_chain)%x(:,k),&
            sc_p%r_sc,sc_p%r_sc)

       if (intersect) then
          clash = .true.
          return
       end if

    end do

    !test if the new sphere intersects other chains in system
    if (sc_p%N_chains .gt. 1) then

       do j_chain = 1,sc_p%N_chains

          if (j_chain .ne. i_chain) then

             do k = 1,sc_s%chain(j_chain)%N

                !test for intersections with the second segment
                intersect = s_s_clash(x_new,sc_s%chain(j_chain)%x(:,k),&
                     sc_p%r_sc,sc_p%r_sc)

                if (intersect) then
                   clash = .true.
                   return
                end if

             end do

          end if

       end do

    end if

    !test for intersections with spherical obstacles
    if (sc_p%obst%N .gt. 0) then

       do j = 1,sc_p%obst%N

          intersect = s_s_clash(x_new,sc_p%obst%x(:,j),&
               sc_p%r_sc,sc_p%R_o)

          if (intersect) then
             clash = .true.
             return
          end if

       end do

    end if

    return

  end function test_inserted_sphere_clash

  
  !subroutine to insert a sphere at a random location
  subroutine insert_random_sphere(sc_s,sc_p,i_chain,insert_fail)

    implicit none

    type(sc_chain_state), intent(inout) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    integer(ip), intent(in) :: i_chain

    logical, intent(out) :: insert_fail

    logical :: accepted, clash

    integer(ip) :: c, max_c, i, ipp

    real(rp) :: mid_dist
    real(rp) :: t, p, c_t, s_t, c_p, s_p
    real(rp) :: dxu
    real(rp) :: dx(1:3), xm(1:3), u(1:3), x_new(1:3)
    
    accepted = .false.

    max_c = 20
    c = 0

    mid_dist = (sqrt(3.0d0)/2.0d0)*(2*sc_p%r)

    do while ((accepted .eqv. .false.) .and. (c .lt. max_c))

       i = int_rand(1,sc_s%chain(i_chain)%N)
       ipp = circ_trans(i+1,sc_s%chain(i_chain)%N)

       dx = sc_s%chain(i_chain)%x(:,ipp)-sc_s%chain(i_chain)%x(:,i)
       xm = (sc_s%chain(i_chain)%x(:,ipp)+sc_s%chain(i_chain)%x(:,i))/2.0d0

       u = dx/norm2(dx)

       t = r_rand(0.0d0,two_pi)
       p = r_rand(0.0d0,pi)

       c_t = cos(t)
       s_t = sin(t)
       c_p = cos(p)
       s_p = sin(p)

       dx = (/c_p*s_t,s_p*s_t,c_t/)

       dxu = dot_product(dx,u)

       if (abs(dxu) .gt. 0.9) cycle

       dx = dx - dxu*u

       dx = dx/norm2(dx)

       x_new = xm + mid_dist*dx

       clash = test_inserted_sphere_clash(sc_s,sc_p,i_chain,i,ipp,x_new)

       !cycle the main loop if there was a clash
       if (clash) then
          c = c + 1
       else
          accepted = .true.
       end if
       
    end do

    if (accepted) then

       if (i .eq. (sc_s%chain(i_chain)%N)) then
          sc_s%chain(i_chain)%x(:,sc_s%chain(i_chain)%N+1) = x_new
       else
          sc_s%chain(i_chain)%x(:,ipp+1:sc_s%chain(i_chain)%N+1) =&
               sc_s%chain(i_chain)%x(:,ipp:sc_s%chain(i_chain)%N)
          sc_s%chain(i_chain)%x(:,ipp) = x_new
       end if

       sc_s%chain(i_chain)%N = sc_s%chain(i_chain)%N + 1

       insert_fail = .false.
       
    else
       insert_fail = .true.
    end if

  end subroutine insert_random_sphere
  

  !subroutine to grow the chain by inserting spherocylinder segments
  subroutine grow_chain_segments(sc_s,sc_p,rep_id,grow_fail)

    implicit none

    type(sc_chain_state), intent(inout) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    logical, intent(out) :: grow_fail
    
    character(len=120), intent(in) :: rep_id

    logical :: insert_fail

    integer(ip) :: N_target, p, progress, dN, N_init, N_min
    integer(ip) :: fail_count, max_fail
    integer(ip) :: queue(1:sc_p%N_chains), i_chain

    grow_fail = .false.
    
    N_target = sc_p%g_p%stages(2,sc_s%stage)

    max_fail = 10000*sc_p%N_chains
    fail_count = 0
    
    N_init = sc_s%chain(sc_p%N_chains)%N
    dN = N_target - N_init
    progress = 0

    queue = (/(i_chain,i_chain=1,sc_p%N_chains)/)

    if (dN .gt. 0) then
       
       do while (sc_s%chain(queue(1))%N .lt. N_init + 1)
          
          call insert_random_segment(sc_s,sc_p,queue(1),insert_fail)
          if (insert_fail) then
             fail_count = fail_count + 1
          else
             queue = cshift(queue,shift=1)
             fail_count = 0
             exit
          end if

          if (fail_count .gt. max_fail) then
             grow_fail = .true.
             return
          end if
          
       end do
       
    else
       return
    end if
    
    !loop until target number of segments are inserted
    write(unit_log,"(A,A,A)")'   [grow_chain_segments] ',&
         trim(rep_id),&
         '--- begin loop over chain additions ---'
    
    do while (sc_s%chain(queue(1))%N .lt. N_target)

       call insert_random_segment(sc_s,sc_p,queue(1),insert_fail)
       
       if (insert_fail) then
          fail_count = fail_count + 1
       else !advance to next chain if insertion was successful
          queue = cshift(queue,shift=1)
       end if

       p = int(floor(1.0d1*(sc_s%chain(queue(1))%N-N_init)/dN))
       if (p > progress) then
          progress = p
          fail_count = 0
          write(unit_log,"(A,A,A,I3,A)")'   [grow_chain_segments] ',&
               trim(rep_id),&
               '--- progress: ',&
               10*p,&
               '/100 ---'
       end if

       if (fail_count .gt. max_fail) then
          grow_fail = .true.
          exit
       end if
       
    end do

  end subroutine grow_chain_segments

  
  !subroutine to grow the chain by inserting spheres
  subroutine grow_chain_spheres(sc_s,sc_p,rep_id,grow_fail)

    implicit none

    type(sc_chain_state), intent(inout) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    logical, intent(out) :: grow_fail
    
    character(len=120), intent(in) :: rep_id

    logical :: insert_fail

    integer(ip) :: N_target, p, progress, dN, N_init, N_min
    integer(ip) :: fail_count, max_fail
    integer(ip) :: queue(1:sc_p%N_chains), i_chain

    grow_fail = .false.
    
    N_target = sc_p%g_p%stages(2,sc_s%stage)

    max_fail = 1000*sc_p%N_chains
    fail_count = 0
    
    N_init = sc_s%chain(sc_p%N_chains)%N
    dN = N_target - N_init
    progress = 0

    queue = (/(i_chain,i_chain=1,sc_p%N_chains)/)

    if (dN .gt. 0) then
       
       do while (sc_s%chain(queue(1))%N .lt. N_init + 1)
          
          call insert_random_sphere(sc_s,sc_p,queue(1),insert_fail)
          if (insert_fail) then
             fail_count = fail_count + 1
          else
             queue = cshift(queue,shift=1)
             fail_count = 0
          end if

          if (fail_count .gt. max_fail) then
             grow_fail = .true.
             return
          end if
          
       end do
       
    else
       return
    end if
    
    !loop until target number of segments are inserted
    do while (sc_s%chain(queue(1))%N .lt. N_target)

       call insert_random_sphere(sc_s,sc_p,queue(1),insert_fail)
       
       if (insert_fail) then
          fail_count = fail_count + 1
       else !advance to next chain if insertion was successful
          queue = cshift(queue,shift=1)
       end if

       p = int(floor(1.0d1*(sc_s%chain(queue(1))%N-N_init)/dN))
       if (p > progress) then
          progress = p
          fail_count = 0
          write(unit_log,"(A,A,A,I3,A)")'   [grow_chain_spheres] ',&
               trim(rep_id),&
               '--- progress: ',&
               10*p,&
               '/100 ---'
       end if

       if (fail_count .gt. max_fail) then
          grow_fail = .true.
          exit
       end if
       
    end do

  end subroutine grow_chain_spheres

  
  !subroutine to interpolate to smaller segments
  subroutine interpolate_segments(sc_s,sc_p,rep_id)

    implicit none

    type(sc_chain_state), intent(inout) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    character(len=120), intent(in) :: rep_id

    integer(ip) :: s_i, N_new, i, i_dx, i_chain

    real(rp) :: s_r
    
    real(rp), allocatable :: dx(:,:), x0(:,:)

    if (sc_s%stage .lt. sc_p%g_p%N_stages) then

       s_i = sc_p%g_p%stages(1,sc_s%stage)/sc_p%g_p%stages(1,sc_s%stage+1)
       s_r = 1.0d0*s_i
       N_new = s_i*sc_s%chain(1)%N

       write(unit_log,"(A,A,A,I3,A,I3,A,I2)")'   [interpolate_segments] ',&
            trim(rep_id),&
            ' interpolation - ',&
            sc_p%g_p%stages(1,sc_s%stage),&
            "  >>>>",&
            sc_p%g_p%stages(1,sc_s%stage+1),&
            ", s=",s_i
       
       if (allocated(dx)) deallocate(dx)
       if (allocated(x0)) deallocate(x0)

       allocate(dx(1:3,sc_s%chain(1)%N))
       allocate(x0(1:3,sc_s%chain(1)%N))

       !loop over the total set of chains
       do i_chain = 1,sc_p%N_chains

          dx = cshift(sc_s%chain(i_chain)%x(:,1:sc_s%chain(i_chain)%N),shift=1,dim=2)&
               - sc_s%chain(i_chain)%x(:,1:sc_s%chain(i_chain)%N)
          x0 = sc_s%chain(i_chain)%x(:,1:sc_s%chain(i_chain)%N)

          !loop over chain segments to perform interpolation
          do i = 1,sc_s%chain(i_chain)%N

             sc_s%chain(i_chain)%x(:,(i-1)*s_i+1) = x0(:,i)

             do i_dx = 1,(s_i-1)
                sc_s%chain(i_chain)%x(:,(i-1)*s_i+1+i_dx) = x0(:,i) + (i_dx/s_r)*dx(:,i)
             end do

          end do !end loop over chain segments

          sc_s%chain(i_chain)%N = N_new

       end do !end loop over chains
    
       sc_s%stage = sc_s%stage + 1

    end if

  end subroutine interpolate_segments


  !function to test clashes during spherocylinder segment movement
  pure function test_move_segment_clash(sc_s,sc_p,i_chain,&
       i,imm,ipp,x_new) result(clash)

    type(sc_chain_state), intent(in) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    integer(ip), intent(in) :: i_chain, i, ipp, imm
    real(rp), intent(in) :: x_new(1:3)

    logical :: clash

    logical :: ang_1, ang_2, intersect_seg1, intersect_seg2, knot

    integer(ip) :: j_chain
    integer(ip) :: j, k, kp

    real(rp) :: R_b_sqrd, cos_ang_thresh
    real(rp) :: sc_seg1(1:3,1:2), sc_seg2(1:3,1:2), sc2(1:3,1:2)

    clash = .false.

    R_b_sqrd = (sc_p%R_b - sc_p%r_sc)**2.0d0
    cos_ang_thresh = cos(2.0d0*third_pi)
 
    !test if new point is within boundary
    if (dot_product(x_new,x_new) .ge. R_b_sqrd) then
       clash = .true.
       ! write(unit_log,*)"CLASH - outside boundary"
       return
    end if

    sc_seg1(:,1) = sc_s%chain(i_chain)%x(:,imm)
    sc_seg1(:,2) = x_new

    sc_seg2(:,1) = x_new
    sc_seg2(:,2) = sc_s%chain(i_chain)%x(:,ipp)

    !test if spherocylinder segments intersect with the chain undergoing addition
    do j = 3,(sc_s%chain(i_chain)%N-3)

       k = circ_trans(i+j,sc_s%chain(i_chain)%N)
       kp = circ_trans(i+j+1,sc_s%chain(i_chain)%N)

       sc2(:,1) = sc_s%chain(i_chain)%x(:,k)
       sc2(:,2) = sc_s%chain(i_chain)%x(:,kp)

       !test for intersections with the first segment
       intersect_seg1 = sc_sc_clash(sc_seg1,sc2,&
            sc_p%r_sc,sc_p%r_sc)

       !test for intersections with the second segment
       intersect_seg2 = sc_sc_clash(sc_seg2,sc2,&
            sc_p%r_sc,sc_p%r_sc)

       if (intersect_seg1 .or. intersect_seg2) then
          ! write(unit_log,*)"CLASH - segment (same chain)"
          clash = .true.
          return
       end if

    end do

    !test for bad angles and nearby intersections
    if (sc_s%chain(i_chain)%N > 4) then

       !test the angle between the first segment and the previous segment
       !test the clash between the second segment and the segment prior to the first

       k = circ_trans(imm-1,sc_s%chain(i_chain)%N)
       kp = circ_trans(imm,sc_s%chain(i_chain)%N)
       sc2(:,1) = sc_s%chain(i_chain)%x(:,k)
       sc2(:,2) = sc_s%chain(i_chain)%x(:,kp)

       intersect_seg2 = sc_sc_clash(sc_seg2,sc2,&
            sc_p%r_sc,sc_p%r_sc)

       ! ang_1 = test_angle(sc_seg1,cshift(sc2,1,dim=2),&
       !      cos_ang_thresh)

       ang_1 = test_angle(sc_seg1,sc2,&
            cos_ang_thresh)

       if (ang_1 .or. intersect_seg2) then
          ! write(unit_log,*)"CLASH - ang_1"
          clash = .true.
          return
       end if

       !test the angle between the second segment and the subsequent segment
       !test the clash between the first segment and the segment subsequent to the second

       k = circ_trans(ipp,sc_s%chain(i_chain)%N)
       kp = circ_trans(ipp+1,sc_s%chain(i_chain)%N)
       sc2(:,1) = sc_s%chain(i_chain)%x(:,k)
       sc2(:,2) = sc_s%chain(i_chain)%x(:,kp)

       intersect_seg1 = sc_sc_clash(sc_seg1,sc2,&
            sc_p%r_sc,sc_p%r_sc)

       ! ang_2 = test_angle(cshift(sc_seg2,1,dim=2),sc2,&
       !      cos_ang_thresh)

       ang_2 = test_angle(sc_seg2,sc2,&
            cos_ang_thresh)

       if (ang_2 .or. intersect_seg1) then
          ! write(unit_log,*)"CLASH - ang_2"
          clash = .true.
          return
       end if

    end if

    !test if spherocylinder segments intersect other chains in system
    if (sc_p%N_chains .gt. 1) then

       do j_chain = 1,sc_p%N_chains

          if (j_chain .ne. i_chain) then

             do k = 1,sc_s%chain(j_chain)%N

                kp = circ_trans(k+1,sc_s%chain(j_chain)%N)

                sc2(:,1) = sc_s%chain(j_chain)%x(:,k)
                sc2(:,2) = sc_s%chain(j_chain)%x(:,kp)

                !test for intersections with the first segment
                intersect_seg1 = sc_sc_clash(sc_seg1,sc2,&
                     sc_p%r_sc,sc_p%r_sc)

                !test for intersections with the second segment
                intersect_seg2 = sc_sc_clash(sc_seg2,sc2,&
                     sc_p%r_sc,sc_p%r_sc)

                if (intersect_seg1 .or. intersect_seg2) then
                   clash = .true.
                   return
                end if

             end do

          end if

       end do

    end if

    ! !test for knots
    ! if (sc_p%knot_prevent .eq. 1) then

    !    !test for self-knotting of chain undergoing addition
    !    if ((sc_s%chain(i_chain)%N > 5) .and.&
    !         (sc_p%g_p%stages(1,sc_s%stage) .ge. 3)) then

    !       do j = 3,(sc_s%chain(i_chain)%N-3)

    !          k = circ_trans(i+j,sc_s%chain(i_chain)%N)
    !          kp = circ_trans(i+j+1,sc_s%chain(i_chain)%N)

    !          sc2(:,1) = sc_s%chain(i_chain)%x(:,k)
    !          sc2(:,2) = sc_s%chain(i_chain)%x(:,kp)

    !          knot = line_intersect_plane(sc2(:,1),sc2(:,2),&
    !               sc_seg1(:,1),sc_seg2(:,1),sc_seg2(:,2))

    !          if (knot) then
    !             clash = .true.
    !             return
    !          end if

    !       end do

    !    end if !end conditional for knots of single chain

    !    !test for knots with other chains
    !    if (sc_p%N_chains .gt. 1) then

    !       do j_chain = 1,sc_p%N_chains

    !          if (j_chain .ne. i_chain) then

    !             do k = 1,sc_s%chain(j_chain)%N

    !                kp = circ_trans(k+1,sc_s%chain(j_chain)%N)

    !                sc2(:,1) = sc_s%chain(j_chain)%x(:,k)
    !                sc2(:,2) = sc_s%chain(j_chain)%x(:,kp)

    !                knot = line_intersect_plane(sc2(:,1),sc2(:,2),&
    !                     sc_seg1(:,1),sc_seg2(:,1),sc_seg2(:,2))

    !                if (knot) then
    !                   clash = .true.
    !                   return
    !                end if

    !             end do

    !          end if

    !       end do

    !    end if !end tests for multiple chains

    ! end if

    !test for intersections with spherical obstacles
    if (sc_p%obst%N .gt. 0) then

       do j = 1,sc_p%obst%N

          intersect_seg1 = sc_s_clash(sc_seg1,sc_p%obst%x(:,j),&
               sc_p%r_sc,sc_p%R_o)

          intersect_seg2 = sc_s_clash(sc_seg2,sc_p%obst%x(:,j),&
               sc_p%r_sc,sc_p%R_o)

          if (intersect_seg1 .or. intersect_seg2) then
             ! write(unit_log,*)"CLASH - obstacle"
             clash = .true.
             return
          end if

       end do

    end if

    return

  end function test_move_segment_clash


  !attempt a segment move in a specified chain
  subroutine attempt_move_in_chain(sc_s,sc_p,i_chain,i)

    implicit none

    type(sc_chain_state), intent(inout) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    integer(ip), intent(in) :: i_chain, i

    logical :: clash

    integer(ip) :: imm, ipp, immm, ippp

    real(rp), dimension(1:3) :: x_new, x, xmm, xmmm, xpp, xppp, dx, u
    real(rp), dimension(1:3) :: dx_m, dx_p
    real(rp), dimension(1:4) :: p, q

    real(rp) :: t
    real(rp) :: E_old, E_new, dE, L
    real(rp) :: r, a

    !propose move of point i in chain i_chain
    imm = circ_trans(i-1,sc_s%chain(i_chain)%N)
    ipp = circ_trans(i+1,sc_s%chain(i_chain)%N)

    x = sc_s%chain(i_chain)%x(:,i)
    xmm = sc_s%chain(i_chain)%x(:,imm)
    xpp = sc_s%chain(i_chain)%x(:,ipp)

    !determine rotation axis
    dx = xpp - xmm
    u = dx/norm2(dx)

    !determine rotation quaternion
    t = r_rand(0.0d0,two_pi)
    t = 0.5d0*t
    q(1) = cos(t)
    q(2:4) = u*sin(t)

    !build a quaternion to rotate from the displacement
    dx = x - xmm
    p(1) = 0.0
    p(2:4) = dx

    !rotate the quaternion and get the new displacement
    p = q_mult(q_mult(q,p),q_inv(q))
    dx = p(2:4)

    !calculate the new point
    x_new = xmm + dx

    !test move for clashes
    clash = test_move_segment_clash(sc_s,sc_p,i_chain,&
         i,imm,ipp,x_new)

    if (clash) then
       ! write(unit_log,*)'clash rejection!'
       return
    end if

    !calculate energy difference
    immm = circ_trans(i-2,sc_s%chain(i_chain)%N)
    ippp = circ_trans(i+2,sc_s%chain(i_chain)%N)
    
    xmmm = sc_s%chain(i_chain)%x(:,immm)
    xppp = sc_s%chain(i_chain)%x(:,ippp)

    E_old = 0.0d0
    !first angle
    dx_m = xmm - xmmm
    dx_p = x - xmm
    E_old = E_old + dot_product(dx_m,dx_p)/(norm2(dx_m)*norm2(dx_p))
    !second angle
    dx_m = x - xmm
    dx_p = xpp - x
    E_old = E_old + dot_product(dx_m,dx_p)/(norm2(dx_m)*norm2(dx_p))
    !third angle
    dx_m = xpp - x
    dx_p = xppp - xpp
    E_old = E_old + dot_product(dx_m,dx_p)/(norm2(dx_m)*norm2(dx_p))

    E_new = 0.0d0
    !first angle
    dx_m = xmm - xmmm
    dx_p = x_new - xmm
    E_new = E_new + dot_product(dx_m,dx_p)/(norm2(dx_m)*norm2(dx_p))
    !second angle
    dx_m = x_new - xmm
    dx_p = xpp - x_new
    E_new = E_new + dot_product(dx_m,dx_p)/(norm2(dx_m)*norm2(dx_p))
    !third angle
    dx_m = xpp - x_new
    dx_p = xppp - xpp
    E_new = E_new + dot_product(dx_m,dx_p)/(norm2(dx_m)*norm2(dx_p))
    
    L = (sc_p%g_p%stages(1,sc_s%stage)*(2*sc_p%r))
    dE = (E_new - E_old)*(sc_p%L_p/L)

    !accept move according to metropolis criterion
    a = min(1.0,exp(-sc_p%beta*dE))
    call random_number(r)

    !replace the coordinate
    if (r .lt. a) then

       sc_s%chain(i_chain)%x(:,i) = x_new

       !write(unit_log,*)'move accepted!'

    end if

  end subroutine attempt_move_in_chain


  !attempt a segment move across the set of chains
  subroutine attempt_move(sc_s,sc_p,N_total)

    type(sc_chain_state), intent(inout) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    integer(ip), intent(in) :: N_total
    
    integer(ip) :: M, i_chain, j_chain
    integer(ip) :: i_move

    i_chain = 1
    i_move = int_rand(1,N_total)

    if (sc_p%N_chains .gt. 1) then

       i_chain = 0
       M = 0
       do while (i_move .gt. M)
          i_chain = i_chain + 1
          M = M + sc_s%chain(i_chain)%N
       end do

       M = 0
       if (i_chain .gt. 1) then
          do j_chain = 1,(i_chain-1)
             M = M + sc_s%chain(j_chain)%N
          end do
       end if

       i_move = i_move - M

    end if

    call attempt_move_in_chain(sc_s,sc_p,i_chain,i_move)

  end subroutine attempt_move


  !move the chain
  subroutine move_segments(sc_s,sc_p,rep_id)

    implicit none

    type(sc_chain_state), intent(inout) :: sc_s
    type(sc_system_params), intent(in) :: sc_p
    
    character(len=120), intent(in) :: rep_id

    integer(ip) :: s, p, progress
    integer(ip) :: i_chain, N_total

    write(unit_log,"(A,A,A,I9)")'   [move_segments] ',&
         trim(rep_id),&
         ' total MC steps - ',&
         sc_p%MC_steps

    !count the total number of segments
    N_total = 0
    do i_chain = 1,sc_p%N_chains
       N_total = N_total + sc_s%chain(i_chain)%N
    end do
    
    s = 0
    progress = 0

    do while (s .lt. sc_p%MC_steps)

       call attempt_move(sc_s,sc_p,N_total)
       
       s = s + 1

       p = int(floor((1.0d1*s)/sc_p%MC_steps))
       if (p > progress) then
          progress = p
          write(unit_log,"(A,A,A,I3,A)")'   [move_segments] ',&
               trim(rep_id),&
               '--- progress: ',&
               10*p,&
               '/100 ---'
       end if

    end do


  end subroutine move_segments

  
  !subroutine to perform sequential chain growths
  subroutine sequential_chain_growth(sc_s,sc_p,rep_id,grow_fail)

    implicit none

    type(sc_chain_state), intent(inout) :: sc_s
    type(sc_system_params), intent(in) :: sc_p
    
    character(len=120), intent(in) :: rep_id

    logical, intent(out) :: grow_fail

    grow_fail = .false.

    do while ((sc_s%stage .lt. sc_p%g_p%N_stages) .and.&
         (grow_fail .eqv. .false.))

       if (sc_p%g_p%stages(1,sc_s%stage) .ne. 1) then

          call grow_chain_segments(sc_s,sc_p,rep_id,grow_fail)

          if (grow_fail .eqv. .false.) then

             call move_segments(sc_s,sc_p,rep_id)

             call interpolate_segments(sc_s,sc_p,rep_id)
             
          end if

       end if

    end do

    if (grow_fail .eqv. .false.) then
       call grow_chain_spheres(sc_s,sc_p,rep_id,grow_fail)
    end if

  end subroutine sequential_chain_growth


  !function calculate local rotation systems along chain
  !based on double reflection method of Wang et al., 2008
  !DOI: 10.1145/1330511.1330513
  function local_rot_systems(x,N) result(R)

    implicit none

    integer(ip), intent(in) :: N
    real(rp), intent(in) :: x(1:3,1:N)

    integer(ip) :: i

    real(rp) :: c1, c2, t1, t2, dot

    real(rp), dimension(1:3,1:N) :: xm2, xm1, xp1, xp2, t
    real(rp), dimension(1:3) :: v1, v2, pL, tL, p, s
    
    real(rp) :: R(1:3,1:3,1:N)

    xm2 = cshift(x,shift=-2,dim=2)
    xm1 = cshift(x,shift=-1,dim=2)
    xp1 = cshift(x,shift=1,dim=2)
    xp2 = cshift(x,shift=2,dim=2)

    !estimate tangent vectors
    t = xm2 - 8.0*xm1 + 8.0*xp1 - xp2

    !normalize tangent vectors
    do i = 1,N
       t(:,i) = t(:,i)/norm2(t(:,i))
    end do

    !prepare local rotation systems
    do i=1,N

       !initialize the first orientiation randomly
       if (i .eq. 1) then

          dot = 1.0

          do while (dot .gt. 0.95d0)
          
             t1 = r_rand(0.0d0,two_pi)
             t2 = r_rand(0.0d0,pi)

             p = (/cos(t1)*sin(t2),sin(t1)*sin(t2),cos(t2)/)

             dot = dot_product(p,t(:,i))

          end do
          
          p = p - dot*t(:,i)

          p = p/norm2(p)

       !otherwise generate using RMF
       else

          p = R(:,1,i-1)

          v1 = xp1(:,i) - x(:,i)

          c1 = dot_product(v1,v1)

          pL = p - (2.0/c1)*dot_product(v1,p)*v1

          tL = t(:,i-1) - (2.0/c1)*dot_product(v1,t(:,i-1))*v1

          v2 = t(:,i) - tL

          c2 = dot_product(v2,v2)

          p = pL - (2.0/c2)*dot_product(v2,pL)*v2
          
       end if

       s = cross_prod(t(:,i),p)
       s = s/norm2(s)

       R(:,1,i) = p
       R(:,2,i) = s
       R(:,3,i) = t(:,i)
       
    end do
    
  end function local_rot_systems

  
  !subroutine to generate quaternion system from coordinates
  subroutine generate_quats(sc_s,sc_p)

    implicit none

    type(sc_chain_state), intent(inout) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    integer(ip) :: i, i_chain

    real(rp), allocatable :: R(:,:,:)
    
    if (allocated(R)) deallocate(R)
    allocate(R(1:3,1:3,1:sc_s%chain(1)%N))

    do i_chain=1,sc_p%N_chains
    
       sc_s%quat(i_chain)%N = sc_s%chain(i_chain)%N

       R = local_rot_systems(sc_s%chain(i_chain)%x,sc_s%chain(i_chain)%N)

       do i = 1,sc_s%quat(i_chain)%N

          sc_s%quat(i_chain)%x(:,i) = q_unit(Rmat_to_q(transpose(R(:,:,i))))

       end do
    
    end do

  end subroutine
  
  !function to generate an integer on the interval [lb,ub]
  function int_rand(lb,ub) result(u_int)
    
    implicit none
    
    integer, intent(in) :: lb, ub
    real(rp) :: u
    integer :: u_int

    call random_number(u)
    u_int=lb+floor((ub+1-lb)*u)

  end function int_rand

  
  !function to generate an real number on the interval [lb,ub]
  function r_rand(lb,ub) result(u_r)
    
    implicit none
    
    real(rp), intent(in) :: lb, ub
    real(rp) :: u
    real(rp) :: u_r

    call random_number(u)
    u_r=lb+(ub-lb)*u

  end function r_rand

  
  !function to calculate circular coordinate along chain
  pure function circ_trans(i,N) result(ic)

    implicit none

    integer(ip), intent(in) :: i, N

    integer(ip) :: ic

    if (i .gt. N) then
       ic = i - N
    elseif (i .le. 0) then
       ic = i + N
    else
       ic = i
    end if

  end function circ_trans

  
  !function to calculate center of mass of chain
  pure function calc_com(sc_s,i_chain) result(r_com)

    implicit none

    type(sc_chain_state), intent(in) :: sc_s
    integer(ip), intent(in) :: i_chain

    real(rp) :: r_com(1:3)

    r_com = sum(sc_s%chain(i_chain)%x,dim=2)/sc_s%chain(i_chain)%N

  end function calc_com

  
  !subroutine to perform cyclic shifts
  subroutine cyclic_shift(sc_s,sc_p)

    implicit none

    type(sc_chain_state), intent(inout) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    integer :: t, mid, i_chain

    real(rp) :: dist(1:sc_s%chain(1)%N)

    mid = sc_s%chain(1)%N/2

    select case (sc_p%cyclic_shift)

    case (1) !ter closest to boundary

       do i_chain=1,sc_p%N_chains

          dist = sum(sc_s%chain(i_chain)%x**2.0d0,dim=1)

          t = maxloc(dist,dim=1)

          sc_s%chain(i_chain)%x = cshift(sc_s%chain(i_chain)%x,t,dim=2)

       end do

    case (2) !ori closest to boundary

       do i_chain=1,sc_p%N_chains

          dist = sum(sc_s%chain(i_chain)%x**2.0d0,dim=1)

          t = maxloc(dist,dim=1)

          sc_s%chain(i_chain)%x = cshift(sc_s%chain(i_chain)%x,t+mid,dim=2)

       end do

    case (3) !ter closest to center

       do i_chain=1,sc_p%N_chains

          dist = sum(sc_s%chain(i_chain)%x**2.0d0,dim=1)

          t = minloc(dist,dim=1)

          sc_s%chain(i_chain)%x = cshift(sc_s%chain(i_chain)%x,t,dim=2)

       end do

    case (4) !ori closest to center

       do i_chain=1,sc_p%N_chains

          dist = sum(sc_s%chain(i_chain)%x**2.0d0,dim=1)

          t = minloc(dist,dim=1)

          sc_s%chain(i_chain)%x = cshift(sc_s%chain(i_chain)%x,t+mid,dim=2)

       end do

    case (5) !ter closest to center of mass

       do i_chain=1,sc_p%N_chains

          dist = sum((sc_s%chain(i_chain)%x-&
               spread(calc_com(sc_s,i_chain),2,sc_s%chain(i_chain)%N))**2.0d0,&
               dim=1)

          t = minloc(dist,dim=1)

          sc_s%chain(i_chain)%x = cshift(sc_s%chain(i_chain)%x,t,dim=2)

       end do

    case (6) !ori closest to center of mass

       do i_chain=1,sc_p%N_chains

          dist = sum((sc_s%chain(i_chain)%x-&
               spread(calc_com(sc_s,i_chain),2,sc_s%chain(i_chain)%N))**2.0d0,&
               dim=1)

          t = minloc(dist,dim=1)

          sc_s%chain(i_chain)%x = cshift(sc_s%chain(i_chain)%x,t+mid,dim=2)

       end do

    case (7) !ter furthest from center of mass

       do i_chain=1,sc_p%N_chains

          dist = sum((sc_s%chain(i_chain)%x-&
               spread(calc_com(sc_s,i_chain),2,sc_s%chain(i_chain)%N))**2.0d0,&
               dim=1)

          t = maxloc(dist,dim=1)

          sc_s%chain(i_chain)%x = cshift(sc_s%chain(i_chain)%x,t,dim=2)

       end do

    case (8) !ori furthest from center of mass

       do i_chain=1,sc_p%N_chains

          dist = sum((sc_s%chain(i_chain)%x-&
               spread(calc_com(sc_s,i_chain),2,sc_s%chain(i_chain)%N))**2.0d0,&
               dim=1)

          t = maxloc(dist,dim=1)

          sc_s%chain(i_chain)%x = cshift(sc_s%chain(i_chain)%x,t+mid,dim=2)

       end do

    end select

  end subroutine cyclic_shift

  
  !function to calculate the chain with minimal length
  pure function calc_N_min(sc_s,sc_p) result(N_min)

    implicit none

    type(sc_chain_state), intent(in) :: sc_s
    type(sc_system_params), intent(in) :: sc_p
    integer(ip) :: i_chain

    integer(ip) :: N_min

    do i_chain = 1,sc_p%N_chains
       if (i_chain .eq. 1) then
          N_min = sc_s%chain(i_chain)%N
       elseif (sc_s%chain(i_chain)%N .lt. N_min) then
          N_min = sc_s%chain(i_chain)%N
       end if
    end do

  end function calc_N_min

end module sc_chain_procedures
