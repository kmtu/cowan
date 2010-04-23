!Coordinated Water Analysis (cowan)
!This program is used to analyze the coordinated water molecules around an ion

PROGRAM cowan
  IMPLICIT NONE
  INTEGER, PARAMETER :: REQUIRED_NUM_PAR = 7
  INTEGER, PARAMETER :: input_fileid = 10
  CHARACTER(LEN=128) :: input_filename
  INTEGER :: count, i, j, stat
  CHARACTER(LEN=128) :: command, usage, arg
  CHARACTER(LEN=8) :: str_timestep, atom_name
  INTEGER :: timestep, num_water, num_atoms, num_frames, init_timestep, delta_timestep, atom_index
  INTEGER :: init_OW_index
  REAL(KIND=8) :: hydration_radius
  REAL(KIND=8), DIMENSION(3) :: pos
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
           write(*,*) "Error: problem occurs while reading timestep:", init_timestep+(i-1)*delta_timestep
           call EXIT(1)
        end if
     end if

     !ignore cell vectors
     read(input_fileid, *) !cell vector
     read(input_fileid, *) !cell vector
     read(input_fileid, *) !cell vector

     !read NA+ position at current timestep
     read(input_fileid, *, IOSTAT=stat) atom_name, atom_index
     if (stat /= 0) then
        write(*,*) "Error: problem occurs while reading NA+ records at timestep:", timestep
        call EXIT(1)
     end if
     read(input_fileid, *, IOSTAT=stat) pos
     if (stat /= 0) then
        write(*,*) "Error: problem occurs while reading NA+ position at timestep:", timestep
        call EXIT(1)
     end if

     !skip non-water atoms
     do j = 1, num_atoms - num_water*3 - 1
        read(input_fileid, *) !skip atom_name, atom_index
        read(input_fileid, *) !skip atom position
     end do

     do j = 1, num_water
     end do
  end do
END PROGRAM cowan
