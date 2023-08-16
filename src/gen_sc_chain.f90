program gen_sc_chain

  use params
  use sc_chain_object
  use sc_chain_procedures
  use quaternions
  use read_write
  use omp_lib

  implicit none

  logical :: log_specified, input_specified, seed_specified
  logical :: output_dir_specified, output_label_specified
  logical :: write_bin, write_xyz, write_dat
  logical :: grow_fail

  integer(ip) :: num_threads
  integer(ip) :: i_arg, num_cl_args
  integer(ip) :: N_stages, N_o
  integer(ip) :: n_seed, seed_in
  integer(ip) :: unit_rep
  integer(ip) :: i_rep
  integer(ip), allocatable :: seed(:)
  integer(ip), allocatable :: stages(:,:)

  real(rp) :: r, r_sc, R_o, R_b
  
  character(len=120) :: input_file, log_file, output_dir
  character(len=120) :: output_label, output_label_rep, rep_id
  character(len=120) :: arg_specifier, arg

  type(sc_chain_state) :: sc_s
  type(sc_system_params) :: sc_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! MAIN PROGRAM !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  input_specified = .false.
  log_specified = .false.
  output_dir_specified = .false.
  output_label_specified = .false.
  seed_specified = .false.
  write_bin = .false.
  write_xyz = .false.
  write_dat = .false.
    
  num_cl_args=command_argument_count()

  if (num_cl_args .gt. 0) then
     
     do i_arg=1,num_cl_args
        
        call getarg(i_arg,arg)
        
        !take argument from terms with "="
        if (index(arg,"=") .ne. 0) then
           arg_specifier=trim(adjustl(arg(:index(arg,"=")-1)))
           arg=trim(adjustl(arg(index(arg,"=")+1:)))
        else
           arg_specifier=trim(adjustl(arg))
        end if

        select case (arg_specifier)
           !test if input is specified
        case ("--input_file", "--i_f")

           input_file=trim(arg)
           write(6,*)"input_file = ",trim(input_file)
           input_specified=.true.

           !test if output_dir is specified
        case ("--output_dir", "--o_d")

           output_dir=trim(arg)
           write(6,*)"output_dir = ",trim(output_dir)
           output_dir_specified=.true.

           !test if output_label is specified
        case ("--output_label", "--o_l")

           output_label=trim(arg)
           write(6,*)"output_label = ",trim(output_label)
           output_label_specified=.true.

           !test if log file was specified
        case ("--log", "--l")

           log_file=trim(arg)
           write(6,*)"log_file = ",trim(log_file)
           log_specified=.true.

           !test if log file was specified
        case ("--seed", "--s")

           read(arg,"(I5)")seed_in
           write(6,"(A,I5)")"seed = ",seed_in
           seed_specified=.true.

           !test if the number of threads were specified
        case ("--num_threads", "--n_t")

           read(arg,"(I3)")num_threads
           num_threads=min(num_threads,max_num_threads)

        case ("--xyz")

           write_xyz = .true.

        case ("--bin")

           write_bin = .true.

        case ("--dat")

           write_dat = .true.
           
        end select
        
     end do
     
  end if

  !if the log unit is not stdout, open a log file at that unit number
  if (unit_log .ne. 6) then
     if (log_specified .eqv. .false.) then
        log_file="./run.log"
     end if
     open(UNIT=unit_log,FILE=trim(log_file),ACTION='write',STATUS='replace')
  end if

  write(unit_log,*)'[gen_sc_chain] BEGIN'

  if (output_dir_specified .eqv. .false.) output_dir = "../output/"

  if (output_label_specified) then
     output_label = "_"//trim(output_label)
  else
     output_label = ""
  end if
  

  write(unit_log,*)' [gen_sc_chain] create new chain from input file'

  call load_input_file(sc_p,input_file)

  if (seed_specified) then
     sc_p%seed = seed_in
  end if

  !initialize the PRNG
  call random_seed(size=n_seed)
  allocate(seed(n_seed))
  seed = sc_p%seed
  call random_seed(put=seed)
  deallocate(seed)


  write(unit_log,*)' [gen_sc_chain] SYSTEM PARAMETERS'
  write(unit_log,*)'   run_flag = ',sc_p%run_flag
  write(unit_log,*)'   min_rep = ',sc_p%min_rep
  write(unit_log,*)'   max_rep = ',sc_p%max_rep
  write(unit_log,*)'   N_reps = ',sc_p%N_reps
  write(unit_log,*)'   seed = ',sc_p%seed
  write(unit_log,*)'   cyclic_shift = ',sc_p%cyclic_shift
  write(unit_log,*)'   knot_prevent = ',sc_p%knot_prevent
  write(unit_log,*)'   MC_steps = ',sc_p%MC_steps
  write(unit_log,*)'   beta = ',sc_p%beta
  write(unit_log,*)'   N_chains = ',sc_p%N_chains
  write(unit_log,*)'   N_stages = ',sc_p%g_p%N_stages
  write(unit_log,*)'   N_final = ',sc_p%g_p%stages(2,sc_p%g_p%N_stages)
  write(unit_log,*)'   N_o = ',sc_p%obst%N
  write(unit_log,*)'   r = ',sc_p%r
  write(unit_log,*)'   r_sc = ',sc_p%r_sc
  write(unit_log,*)'   R_o = ',sc_p%R_o
  write(unit_log,*)'   R_b = ',sc_p%R_b
  write(unit_log,*)'   L_p = ',sc_p%L_p
  write(unit_log,*)'   write_bin = ',write_bin
  write(unit_log,*)'   write_dat = ',write_dat
  write(unit_log,*)'   write_xyz = ',write_xyz


  if (sc_p%run_flag .eq. 1) then

     do i_rep = 1,sc_p%N_reps

        write(rep_id,"(A,I5.5,A)")'(',(sc_p%min_rep+i_rep-1),')'
        rep_id=trim(rep_id)

        write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' BEGIN'

        unit_rep = offset_units

        output_label_rep = trim(output_label)//"_rep"
        write(output_label_rep,"(A,I5.5)")trim(output_label_rep),(sc_p%min_rep+i_rep-1)

        grow_fail = .true.
        
        do while (grow_fail .eqv. .true.)
        
           !initialize the chain
           write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' initialize chain'
           call initialize_chain(sc_s,sc_p)

           !perform sequential chain growths
           write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' sequential growth'
           call sequential_chain_growth(sc_s,sc_p,rep_id,grow_fail)

        end do

        !perform cyclic shift if needed
        write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' cyclic shifts'
        if (sc_p%cyclic_shift .ne. 0) call cyclic_shift(sc_s,sc_p)

        !calculate quaternions
        write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' generate quaternions'
        call generate_quats(sc_s,sc_p)

        !write plain-text files for python input
        if (write_dat .eqv. .true.) then
           write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' writing dat files'
           call write_catcoords_plain(sc_s%chain,sc_p%N_chains,unit_rep,&
                trim(output_dir)//'x_chain'//trim(output_label_rep)//'.dat')
           call write_coords_plain(sc_p%obst,unit_rep,&
                trim(output_dir)//'x_obst'//trim(output_label_rep)//'.dat')
           call write_catcoords_plain(sc_s%quat,sc_p%N_chains,unit_rep,&
                trim(output_dir)//'x_quat'//trim(output_label_rep)//'.dat')
        end if

        !write binary files
        if (write_bin .eqv. .true.) then
           write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' writing bin files'
           call write_catcoords_binary(sc_s%chain,sc_p%N_chains,unit_rep,&
                trim(output_dir)//'x_chain'//trim(output_label_rep)//'.bin')
           call write_coords_binary(sc_p%obst,unit_rep,&
                trim(output_dir)//'x_obst'//trim(output_label_rep)//'.bin')
           call write_catcoords_binary(sc_s%quat,sc_p%N_chains,unit_rep,&
                trim(output_dir)//'x_quat'//trim(output_label_rep)//'.bin')
        end if

        !write xyz files for VMD visualization
        if (write_xyz .eqv. .true.) then
           write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' writing xyz files'
           call write_catcoords_xyz(sc_s%chain,sc_p%N_chains,unit_rep,&
                trim(output_dir)//'x_chain'//trim(output_label_rep)//'.xyz')
           call write_coords_xyz(sc_p%obst,unit_rep,&
                trim(output_dir)//'x_obst'//trim(output_label_rep)//'.xyz')
        end if

        write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' END'

     end do

  elseif (sc_p%run_flag .eq. 2) then

     call omp_set_num_threads(num_threads)

     !$OMP PARALLEL DO &
     !$OMP DEFAULT(PRIVATE) &
     !$OMP SHARED(sc_p,&
     !$OMP output_dir, output_label,&
     !$OMP write_xyz,write_bin,write_dat)

     do i_rep = 1,sc_p%N_reps

        write(rep_id,"(A,I5.5,A)")'(',(sc_p%min_rep+i_rep-1),')'
        rep_id=trim(rep_id)

        write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' BEGIN'

        unit_rep = omp_get_thread_num() + offset_units

        output_label_rep = trim(output_label)//"_rep"
        write(output_label_rep,"(A,I5.5)")trim(output_label_rep),(sc_p%min_rep+i_rep-1)

        grow_fail = .true.
        
        do while (grow_fail .eqv. .true.)
        
           !initialize the chain
           write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' initialize chain'
           call initialize_chain(sc_s,sc_p)

           !perform sequential chain growths
           write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' sequential growth'
           call sequential_chain_growth(sc_s,sc_p,rep_id,grow_fail)

        end do

        !perform cyclic shift if needed
        write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' cyclic shifts'
        if (sc_p%cyclic_shift .ne. 0) call cyclic_shift(sc_s,sc_p)

        !calculate quaternions
        write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' generate quaternions'
        call generate_quats(sc_s,sc_p)

        !write plain-text files for python input
        if (write_dat .eqv. .true.) then
           write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' writing dat files'
           call write_catcoords_plain(sc_s%chain,sc_p%N_chains,unit_rep,&
                trim(output_dir)//'x_chain'//trim(output_label_rep)//'.dat')
           if (sc_p%obst%N .gt. 0) then
              call write_coords_plain(sc_p%obst,unit_rep,&
                   trim(output_dir)//'x_obst'//trim(output_label_rep)//'.dat')
           end if
           call write_catcoords_plain(sc_s%quat,sc_p%N_chains,unit_rep,&
                trim(output_dir)//'x_quat'//trim(output_label_rep)//'.dat')
        end if

        !write binary files
        if (write_bin .eqv. .true.) then
           write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' writing bin files'
           call write_catcoords_binary(sc_s%chain,sc_p%N_chains,unit_rep,&
                trim(output_dir)//'x_chain'//trim(output_label_rep)//'.bin')
           if (sc_p%obst%N .gt. 0) then
              call write_coords_binary(sc_p%obst,unit_rep,&
                   trim(output_dir)//'x_obst'//trim(output_label_rep)//'.bin')
           end if
           call write_catcoords_binary(sc_s%quat,sc_p%N_chains,unit_rep,&
                trim(output_dir)//'x_quat'//trim(output_label_rep)//'.bin')
        end if

        !write xyz files for VMD visualization
        if (write_xyz .eqv. .true.) then
           write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' writing xyz files'
           call write_catcoords_xyz(sc_s%chain,sc_p%N_chains,unit_rep,&
                trim(output_dir)//'x_chain'//trim(output_label_rep)//'.xyz')
           if (sc_p%obst%N .gt. 0) then
              call write_coords_xyz(sc_p%obst,unit_rep,&
                   trim(output_dir)//'x_obst'//trim(output_label_rep)//'.xyz')
           end if
        end if

        write(unit_log,*)' [gen_sc_chain] '//trim(rep_id)//' END'

     end do

  end if
  
  write(unit_log,*)'[gen_sc_chain] END'

  !if not writing to stdout, then close the log file
  if (unit_log .ne. 6) then
     close(unit_log)
  end if

end program gen_sc_chain
