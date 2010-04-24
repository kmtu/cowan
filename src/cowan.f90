!Coordinated Water Analysis (cowan)
!This program is used to analyze the coordinated water molecules around an ion

PROGRAM cowan
  IMPLICIT NONE
  INTEGER, PARAMETER :: REQUIRED_NUM_PAR = 7
  INTEGER, PARAMETER :: input_fileid = 10
  CHARACTER(LEN=128) :: input_filename
  INTEGER, PARAMETER :: output_fileid = 11
  CHARACTER(LEN=128) :: output_filename
  INTEGER :: count, i, j, k, stat
  CHARACTER(LEN=128) :: command, usage, arg
  CHARACTER(LEN=2) :: out_fmt
  CHARACTER(LEN=8) :: str_timestep, atom_name
  CHARACTER, PARAMETER :: OW_IN = '*', OW_OUT = '_'
  CHARACTER :: OW_state
  CHARACTER(LEN=2), PARAMETER :: OW_name = 'OW'
  INTEGER :: timestep, num_water, num_atoms, num_frames, init_timestep, delta_timestep, atom_index
  INTEGER :: init_OW_index, num_digits_atom_index
  REAL(KIND=8) :: hydration_radius, hydration_radius_square, dist_square
  REAL(KIND=8), DIMENSION(3) :: pos_ion, pos_water
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: hydration_table
  
  !-- read parameters from command line arguments --!
  count = COMMAND_ARGUMENT_COUNT()
  call GET_COMMAND_ARGUMENT(NUMBER=0, VALUE=command)
  usage = "Usage: " // TRIM(ADJUSTL(command)) // " input_filename &
          &-f <frame#> -n <H2O#> -r <hydration radius>"
  
  if (count < REQUIRED_NUM_PAR) then
     write(*,"(A,I1,A)") "At least ", (REQUIRED_NUM_PAR+1)/2, " arguments are needed."
     write(*,*) usage
     call EXIT(1)
  end if
  call GET_COMMAND_ARGUMENT(NUMBER=1, VALUE=input_filename, STATUS=stat)
  if (stat /= 0) then
     write(*,*) "Unable to get the input_filename."
     write(*,*) usage        
     call EXIT(1)
  end if
  
  i = 2
  do while (i <= count)
     call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=stat)
     select case (arg)
     case ('-n')
        i = i + 1
        call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=stat)
        if (stat /= 0) then
           write(*,*) "Unable to read the value of argument -n"
           write(*,*) usage
           call EXIT(1)
        end if
        read(arg,*, IOSTAT=stat) num_water
        if (stat/=0) then
           write(*,*) "Unable to read the value of argument -n, only INTEGER is allowed!"
           call EXIT(1)
        end if
     case ('-r')
        i = i + 1
        call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=stat)
        if (stat /= 0) then
           write(*,*) "Unable to read the value of argument -r"
           write(*,*) usage
           call EXIT(1)
        end if
        read(arg,*, IOSTAT=stat) hydration_radius
        if (stat/=0) then
           write(*,*) "Unable to read the value of argument -r, a real number is needed!"
           call EXIT(1)
        end if
     case ('-f')
        i = i + 1
        call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=stat)
        if (stat /= 0) then
           write(*,*) "Unable to read the value of argument -f"
           write(*,*) usage
           call EXIT(1)
        end if
        read(arg,*, IOSTAT=stat) num_frames
        if (stat/=0) then
           write(*,*) "Unable to read the value of argument -f, only INTEGER is allowed!"
           call EXIT(1)
        end if
     case default
        write(*,*) "Undefined argument: ", arg
        write(*,*) usage
        call EXIT(1)
     end select
     i = i + 1
  end do
  !-- end of reading parameters from command line arguments --!

  !output parameters being read for confirmation
  write(*,*) "input filename: ", TRIM(input_filename)
  write(*,*) "number of frames:", num_frames  
  write(*,"(' number of water molecules: ',I5)") num_water
  write(*,*) "hydration radius:", hydration_radius

  !square the hydration_radius for furthur comparison with the dist_square between water and ion
  hydration_radius_square = hydration_radius * hydration_radius

  allocate(hydration_table(num_water,num_frames), STAT=stat)
  if (stat /=0) then
     write(*,*) "Error: failed to allocate hydration_table(num_water, num_frames)!"
     call EXIT(1)
  end if
  
  open(UNIT=input_fileid, FILE=input_filename, STATUS='OLD', IOSTAT=stat, ACTION='READ')
  if (stat /= 0) then
     write(*,*) "Error: unable to open file: ", input_filename
     call EXIT(1)
  end if

  read(input_fileid, *) !title line
  read(input_fileid, *) i, j, num_atoms

  init_OW_index = num_atoms - num_water*3 + 1

  ! start reading records of each timestep
  do i = 1, num_frames
     read(input_fileid, *, IOSTAT=stat) str_timestep, timestep
     if (i == 1) then
        if (stat /= 0) then
           write(*,*) "Error: problem occurs while reading the 1st timestep"
           call EXIT(1)
        end if
        init_timestep = timestep
     else if (i == 2) then
        if (stat /= 0) then
           write(*,*) "Error: problem occurs while reading the 2nd timestep"
           call EXIT(1)
        end if
        delta_timestep = timestep - init_timestep
     else
        if (stat /= 0) then
           write(*,*) "Error: problem occurs while reading timestep:", &
                &init_timestep+(i-1)*delta_timestep
           call EXIT(1)
        end if
     end if

     !ignore cell vectors
     read(input_fileid, *) !cell vector
     read(input_fileid, *) !cell vector
     read(input_fileid, *) !cell vector

     !read ion position at current timestep
     read(input_fileid, *, IOSTAT=stat) atom_name, atom_index
     if (stat /= 0) then
        write(*,*) "Error: problem occurs while reading ion records at timestep:", timestep
        call EXIT(1)
     end if
     read(input_fileid, *, IOSTAT=stat) pos_ion
     if (stat /= 0) then
        write(*,*) "Error: problem occurs while reading ion position at timestep:", timestep
        call EXIT(1)
     end if

     !skip non-water atoms ()
     do j = 1, num_atoms - num_water*3 - 1
        read(input_fileid, *) !skip atom_name, atom_index
        read(input_fileid, *) !skip atom position
     end do

     !read water molecule position and decide if it is inside the first hydration layer
     do j = 1, num_water
        read(input_fileid, *, IOSTAT=stat) atom_name, atom_index
        if (stat /= 0) then
           write(*,*) "Error: problem occurs while reading water record at timestep:", timestep
           call EXIT(1)
        else if (atom_name /= OW_name .or. atom_index /= init_OW_index + (j-1)*3) then
           write(*,*) "Error: problem occurs while reading water record at timestep:", timestep
           write(*,*) "       the OW index should be ", init_OW_index + (j-1)*3, &
                & "but ", atom_index, " was read."
           call EXIT(1)
        else
           read(input_fileid, *, IOSTAT=stat) pos_water
           if (stat /= 0) then
              write(*,*) "Error: problem occurs while reading water position at timestep:", &
                   &timestep
              call EXIT(1)
           end if
           do k = 1, 4
              read(input_fileid, *) !omit the HW records
           end do
        end if

        dist_square = 0.0
        do k = 1, 3
           dist_square = dist_square + (pos_water(k) - pos_ion(k))*(pos_water(k) - pos_ion(k))
        end do
        
        if (dist_square <= hydration_radius_square) then
           hydration_table(j, i) = .true.
        else
           hydration_table(j, i) = .false.
        end if
     end do
  end do

  close(input_fileid)

  ! output results
  output_filename = TRIM(ADJUSTL(input_filename)) // ".cowan"
  open(UNIT=output_fileid, FILE=output_filename, STATUS='REPLACE', IOSTAT=stat, ACTION='WRITE')
  if (stat /= 0) then
     write(*,*) "Error: unable to create file: ", output_filename
     call EXIT(1)
  end if

  ! decide the max number of digits of timestep for output
  ! plus 2 is the safest way (because 10000 may be 9999.999 in REAL format)
!  max_num_timestep_digits = INT(LOG10(REAL(init_timestep + (num_frames-1)*delta_timestep))) + 2

  ! decide the number of digits of atom index for output
  ! plus 2 is the safest way (because 10000 may be 9999.999 in REAL format)
  num_digits_atom_index = INT(LOG10(REAL(num_atoms))) + 2
  
  do i = 1, num_water
     if (ANY(hydration_table(i,:))) then
        do j = 0, num_frames
           if (j == 0) then
              write(out_fmt, "(I2)") num_digits_atom_index
              write(output_fileid, "(I"//TRIM(ADJUSTL(out_fmt))//")", IOSTAT=stat, ADVANCE='NO') &
                   & init_OW_index + (i-1)*3
           else
              if (hydration_table(i,j)) then
                 OW_state = OW_IN
              else
                 OW_state = OW_OUT
              end if

              write(output_fileid, "(A1)", ADVANCE='NO') OW_state

              if (j == num_frames) then
                 write(output_fileid, *) ! newline
              end if
           end if
        end do
     end if
  end do

  deallocate(hydration_table)
  close(output_fileid)
END PROGRAM cowan
