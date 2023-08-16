module read_write

  !module used to handle reading the input file and
  !writing the restart and output files
  
  use params
  use sc_chain_object
  use sc_chain_procedures

  implicit none

contains

  !routine to load the input file
  subroutine load_input_file(sc_p,input_file)

    implicit none

    type(sc_system_params), intent(inout) :: sc_p
    character(len=120), intent(in) :: input_file

    integer(ip) :: unit_input, stat
    integer(ip) :: run_flag, min_rep, max_rep, N_reps
    integer(ip) :: seed, N_o, N_stages, N_chains
    integer(ip) :: cyclic_shift, knot_prevent, MC_steps
    integer(ip) :: i_stage, i_obst
    integer(ip), allocatable :: stages(:,:)

    real(rp) :: r, r_sc, R_o, R_b, L_p, beta
    
    character(len=120) :: temp_line, arg_specifier, arg, obstacle_file

    type(coords) :: obst_temp

    unit_input = offset_units

    run_flag = 0
    min_rep = 0
    max_rep = 0
    N_chains = 1
    cyclic_shift = 0
    MC_steps = 0
    beta = 1.0

    write(unit_log,*)input_file
    
    open(UNIT=unit_input,FILE=trim(input_file),ACTION='read',STATUS='old')

    do

       read(unit_input,"(A)",IOSTAT=stat)temp_line
       !write(unit_log,*)temp_line

       if (IS_IOSTAT_END(stat)) then
          exit
       else

          !take argument from terms with "="
          if (index(temp_line,"=") .ne. 0) then
             arg_specifier=trim(adjustl(temp_line(:index(temp_line,"=")-1)))
             arg=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
          else
             arg_specifier=trim(adjustl(temp_line))
          end if

          select case(arg_specifier)

          case ('run')

             read(arg,"(I5)")run_flag

          case ('min_rep')

             read(arg,"(I5)")min_rep

          case ('max_rep')

             read(arg,"(I5)")max_rep

          case ('seed')

             read(arg,"(I5)")seed
             
          case ('N_chains')

             read(arg,"(I5)")N_chains

          case ('cyclic_shift')

             read(arg,"(I5)")cyclic_shift

          case ('knot_prevent')

             read(arg,"(I5)")knot_prevent

          case ('MC_steps')

             read(arg,"(I9)")MC_steps

          case ('beta')

             read(arg,*)beta

          case ('N_o')

             read(arg,"(I5)")N_o

          case ('obstacle_file')

             read(arg,"(A)")obstacle_file

             if (N_o .eq. -1) then
                call read_obstacles(obst_temp,obstacle_file)
                N_o = obst_temp%N
             end if

          case ('r')

             read(arg,*)r

          case ('r_sc')

             read(arg,*)r_sc

          case ('R_o')

             read(arg,*)R_o

          case ('R_b')

             read(arg,*)R_b

          case ('L_p')

             read(arg,*)L_p

          case ('N_stages')

             read(arg,"(I5)")N_stages

             if (allocated(stages)) deallocate(stages)
             allocate(stages(1:2,1:N_stages))

             !read the stages list
             do i_stage=1,N_stages
                read(unit_input,*)stages(:,i_stage)
             end do

             !write(unit_log,*)stages
             
          end select

       end if
       
    end do

    write(unit_log,*)stages


    call new_sc_system_params(sc_p,&
         run_flag,&
         seed,&
         min_rep,max_rep,&
         N_chains,&
         cyclic_shift,knot_prevent,MC_steps,&
         N_stages,stages,&
         N_o,&
         r,r_sc,R_o,R_b,L_p,beta)

    deallocate(stages)

    close(unit_input)

    if (N_o .eq. -1) then
       do i_obst = 1,sc_p%obst%N
          sc_p%obst%x(:,i_obst) = obst_temp%x(:,i_obst)
       end do
       sc_p%obst%x = sc_p%r*sc_p%obst%x
    else
       call initialize_obstacles(sc_p)
    end if

    !write(unit_log,*)'finished loading input file'
    
  end subroutine load_input_file

  subroutine read_obstacles(obst,obstacle_file)

    implicit none

    type(coords), intent(inout) :: obst
    character(len=120), intent(in) :: obstacle_file

    integer :: unit_input
    integer(ip) :: i_obst, N, l

    character(len=120) :: temp_line

    unit_input = offset_units + 1

    open(UNIT=unit_input,FILE=trim(obstacle_file),ACTION='read',STATUS='old')

    !read the number of obstacles
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I6)")N

    !create a new set of coordinates
    call new_coords(obst,3,N)
    obst%N = N

    !read the obstacle coordinates
    do i_obst=1,N
       read(unit_input,*)l,obst%x(:,i_obst)
    end do

    close(unit_input)

  end subroutine read_obstacles

  subroutine write_coords_plain(c,unit,file)

    implicit none

    type(coords), intent(in) :: c
    integer(ip), intent(in) :: unit
    character(len=*), intent(in) :: file

    integer(ip) :: i, j
    character(len=120) :: num_string, format_string, x_string

    !open the file
    open(UNIT=unit,&
         FILE=trim(adjustl(file)),&
         STATUS='replace',&
         ACTION='write')


    num_string = "ES20.10"
    format_string="(I1,'(A1,("//trim(num_string)//"))')"
    write(x_string,format_string)(c%dim-1)
    x_string="("//trim(num_string)//","//trim(x_string)//")"

    !write the coordinates as plain-text
    do i=1,c%N
       write(unit,x_string)c%x(1,i),(",",c%x(j,i),j=2,c%dim)
    end do

    !close the file
    close(unit)
    
  end subroutine write_coords_plain

  subroutine write_catcoords_plain(c,Nc,unit,file)

    implicit none

    type(coords), intent(in) :: c(1:Nc)
    integer(ip), intent(in) :: Nc
    integer(ip), intent(in) :: unit
    character(len=*), intent(in) :: file

    type(coords) :: c_cat

    integer(ip) :: ic, N_cat
    
    character(len=120) :: num_string, format_string, x_string

    N_cat = 0
    do ic=1,Nc
       
       N_cat = N_cat + c(ic)%N
       
    end do

    call new_coords(c_cat,c(1)%dim,N_cat)
    c_cat%N = N_cat

    N_cat = 0
    do ic=1,Nc
       
       c_cat%x(:,N_cat+1:N_cat+c(ic)%N) = c(ic)%x
       N_cat = N_cat + c(ic)%N
       
    end do

    call write_coords_plain(c_cat,unit,file)
    
  end subroutine write_catcoords_plain

  subroutine write_coords_binary(c,unit,file)

    implicit none

    type(coords), intent(in) :: c
    integer(ip), intent(in) :: unit
    character(len=*), intent(in) :: file

    integer :: i, j
    character(len=120) :: num_string, format_string, x_string

    !open the file
    open(UNIT=unit,&
         FILE=trim(adjustl(file)),&
         STATUS='replace',&
         FORM='unformatted',&
         RECL=c%dim*c%N*rp,&
         ACCESS='direct')


    write(unit,rec=1)c%x(:,1:c%N)

    !close the file
    close(unit)
    
  end subroutine write_coords_binary

  subroutine write_catcoords_binary(c,Nc,unit,file)

    implicit none

    type(coords), intent(in) :: c(1:Nc)
    integer(ip), intent(in) :: Nc
    integer(ip), intent(in) :: unit
    character(len=*), intent(in) :: file

    type(coords) :: c_cat

    integer(ip) :: ic, N_cat
    
    character(len=120) :: num_string, format_string, x_string

    N_cat = 0
    do ic=1,Nc
       
       N_cat = N_cat + c(ic)%N
       
    end do

    call new_coords(c_cat,c(1)%dim,N_cat)
    c_cat%N = N_cat

    N_cat = 0
    do ic=1,Nc
       
       c_cat%x(:,N_cat+1:N_cat+c(ic)%N) = c(ic)%x
       N_cat = N_cat + c(ic)%N
       
    end do

    call write_coords_binary(c_cat,unit,file)
    
  end subroutine write_catcoords_binary

  subroutine write_coords_xyz(c,unit,file)

    implicit none

    type(coords), intent(in) :: c
    integer(ip), intent(in) :: unit
    character(len=*), intent(in) :: file

    integer(ip) :: i, j
    character(len=120) :: num_string, format_string, x_string

    !open the file
    open(UNIT=unit,&
         FILE=trim(adjustl(file)),&
         STATUS='replace',&
         ACTION='write')


    num_string = "ES20.10"
    format_string="(I1,'(A1,("//trim(num_string)//"))')"
    write(x_string,format_string)(c%dim)
    x_string="(A,"//trim(x_string)//")"

    !write the coordinates as an xyz file

    write(unit,*)c%N
    write(unit,*)
    
    do i=1,c%N
       write(unit,x_string)"C",(" ",c%x(j,i),j=1,c%dim)
    end do

    !close the file
    close(unit)
    
  end subroutine write_coords_xyz

  subroutine write_catcoords_xyz(c,Nc,unit,file)

    implicit none

    type(coords), intent(in) :: c(1:Nc)
    integer(ip), intent(in) :: Nc
    integer(ip), intent(in) :: unit
    character(len=*), intent(in) :: file

    type(coords) :: c_cat

    integer(ip) :: ic, N_cat
    
    character(len=120) :: num_string, format_string, x_string

    N_cat = 0
    do ic=1,Nc
       
       N_cat = N_cat + c(ic)%N
       
    end do

    call new_coords(c_cat,c(1)%dim,N_cat)
    c_cat%N = N_cat

    N_cat = 0
    do ic=1,Nc
       
       c_cat%x(:,N_cat+1:N_cat+c(ic)%N) = c(ic)%x
       N_cat = N_cat + c(ic)%N
       
    end do

    call write_coords_xyz(c_cat,unit,file)
    
  end subroutine write_catcoords_xyz

end module read_write
