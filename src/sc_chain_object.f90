module sc_chain_object
  
  use params
  
  implicit none

  !derived type used to hold coordinate information
  type coords
     integer(ip) :: dim
     integer(ip) :: N_max, N
     real(rp), allocatable :: x(:,:)
  end type coords

  !derived type used to hold the growth information
  type growth_params
     integer(ip) :: N_stages
     integer(ip), allocatable :: stages(:,:)
  end type growth_params

  !derived type used to hold the chain state
  type sc_chain_state
     integer(ip) :: stage
     type(coords), allocatable :: chain(:), quat(:)
  end type sc_chain_state

  !derived type used to store the system description
  type sc_system_params
     real(rp) :: r, r_sc, R_o, R_b, L_p, beta
     integer(ip) :: run_flag, min_rep, max_rep, N_reps
     integer (ip) :: N_chains
     integer(ip) :: seed
     integer(ip) :: cyclic_shift, knot_prevent, MC_steps
     type(coords) :: obst
     type(growth_params) :: g_p
  end type sc_system_params
  
contains

  !subroutine to create new coords
  subroutine new_coords(c,&
       dim, N_max)
    
    implicit none
    integer(ip), intent(in) :: dim
    integer(ip) :: N_max
    type(coords), intent(inout) :: c
    
    c%dim = dim
    c%N_max = N_max
    c%N = 0
    
    if (allocated(c%x)) deallocate(c%x)
    if (N_max .gt. 0) then
       allocate(c%x(1:dim,1:N_max))
       c%x = 0.0d0
    end if
    
  end subroutine new_coords

  !subroutine to create new growth_params
  subroutine new_growth_params(g_p,&
       N_stages)

    implicit none
    integer(ip), intent(in) :: N_stages
    type(growth_params), intent(inout) :: g_p

    g_p%N_stages = N_stages
    
    if (allocated(g_p%stages)) deallocate(g_p%stages)
    if (N_stages .gt. 0) then
       allocate(g_p%stages(1:2,1:N_stages))
    end if
    
  end subroutine new_growth_params

  !subroutine to create new sc_chain_state
  subroutine new_sc_chain_state(sc_s,sc_p)

    implicit none
    
    type(sc_chain_state), intent(inout) :: sc_s
    type(sc_system_params), intent(in) :: sc_p

    integer(ip) :: i_chain


    ! if (allocated(sc_s%stage)) deallocate(sc_s%stage)
    ! allocate(sc_s%stage(1:sc_p%N_chains))
    if (allocated(sc_s%chain)) deallocate(sc_s%chain)
    allocate(sc_s%chain(1:sc_p%N_chains))
    if (allocated(sc_s%quat)) deallocate(sc_s%quat)
    allocate(sc_s%quat(1:sc_p%N_chains))

    sc_s%stage = 1

    do i_chain=1,sc_p%N_chains
       !create the coordinates object
       call new_coords(sc_s%chain(i_chain), 3, sc_p%g_p%stages(2,sc_p%g_p%N_stages))

       !create the coordinates object for the quaternions
       call new_coords(sc_s%quat(i_chain), 4, sc_s%chain(i_chain)%N_max)
    end do
    
  end subroutine new_sc_chain_state

  !subroutine to create new sc_system_params
  subroutine new_sc_system_params(sc_p,&
       run_flag,&
       seed,&
       min_rep,max_rep,&
       N_chains,&
       cyclic_shift,knot_prevent,MC_steps,&
       N_stages, stages, N_o,&
       r, r_sc, R_o, R_b, L_p, beta)

    implicit none
    integer(ip) :: run_flag, min_rep, max_rep
    integer(ip), intent(in) :: seed, cyclic_shift, knot_prevent, MC_steps
    integer(ip), intent(in) :: N_stages, N_o, N_chains
    integer(ip), allocatable, intent(in) :: stages(:,:)
    real(rp), intent(in) :: r, r_sc, R_o, R_b, L_p, beta
    
    type(sc_system_params), intent(inout) :: sc_p

    !set the run_flag
    sc_p%run_flag = run_flag

    !set the replicate numbers
    sc_p%min_rep = min_rep
    sc_p%max_rep = max_rep
    sc_p%N_reps = max_rep - min_rep + 1

    !set the seed
    sc_p%seed = seed

    !set the number of chains
    sc_p%N_chains = N_chains

    !set the cyclic shift and knot prevention
    sc_p%cyclic_shift = cyclic_shift
    sc_p%knot_prevent = knot_prevent
    sc_p%MC_steps = MC_steps
    
    !set the dimensions of the sphero-cylinders, monomers, obstacles, and boundary
    sc_p%r = r
    sc_p%r_sc = r_sc
    sc_p%R_o = R_o
    sc_p%R_b = R_b
    
    sc_p%L_p = L_p
    sc_p%beta = beta

    !create the growth_params object
    if (stages(1,N_stages) .ne. 1) then
       call new_growth_params(sc_p%g_p, N_stages+1)
       sc_p%g_p%stages(:,1:N_stages) = stages
       sc_p%g_p%stages(:,N_stages+1) = (/1,&
            sc_p%g_p%stages(1,N_stages)*sc_p%g_p%stages(2,N_stages)/)
    else
       call new_growth_params(sc_p%g_p, N_stages)
       sc_p%g_p%stages = stages
    end if

    !create the obstacles object
    call new_coords(sc_p%obst, 3, N_o)
    
  end subroutine new_sc_system_params
  
end module sc_chain_object
