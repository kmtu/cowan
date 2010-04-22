!Coordinated Water Analysis (cowan)
!This program is used to analyze the coordinated water molecules around an ion

PROGRAM cowan
  IMPLICIT NONE
  INTEGER, PARAMETER :: REQUIRED_NUM_PAR = 5
  INTEGER, PARAMETER :: input_fileid = 10
  CHARACTER(LEN=128) :: input_filename
  INTEGER :: count, i, j, stat
  CHARACTER(LEN=128) :: command, usage, arg
  CHARACTER(LEN=8) :: str_timestep
  INTEGER :: timestep, num_water, num_atoms
  REAL(KIND=8) :: hydration_radius
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: hydration_table
  
  !read the parameters from command line arguments
  count = COMMAND_ARGUMENT_COUNT()
  call GET_COMMAND_ARGUMENT(NUMBER=0, VALUE=command)
  usage = "Usage: " // TRIM(ADJUSTL(command)) // " input_filename &
          &-n <# of H2O> -r <hydration radius>"
  
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
     if (arg == '-n') then
        i = i + 1
        call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=stat)
        if (stat /= 0) then
           write(*,*) "Unable to read the value of argument -n"
           write(*,*) usage
           call EXIT(1)
        end if
        read(arg,*) num_water
     else if (arg == '-r') then
        i = i + 1
        call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=sntat)
        if (stat /= 0) then
           write(*,*) "Unable to read the value of argument -r"
           write(*,*) usage
           call EXIT(1)
        end if
        read(arg,*) hydration_radius
     else
        write(*,*) "Undefined argument: ", arg
        write(*,*) usage
        call EXIT(1)
     end if
     i = i + 1
  end do

  write(*,*) "input filename: ", input_filename
  write(*,"('number of water molecules: ',I5)") num_water
  
  open(UNIT=input_fileid, FILE=input_filename, STATUS='OLD', IOSTAT=stat, ACTION='READ')
  if (stat /= 0) then
     write(*,*) "Error: unable to open file: ", input_filename
     call EXIT(1)
  end if

  read(UNIT=input_fileid, *) !title line
  read(UNIT=input_fileid, *) i, j, num_atoms

  allocate()
  
  read(UNIT=input_fileid, *) str_timestep, timestep
  read(UNIT=input_fileid, *) !cell vector
  read(UNIT=input_fileid, *) !cell vector
  read(UNIT=input_fileid, *) !cell vector  

  read(
  
END PROGRAM cowan
