!Coordinated Water Analysis (cowan)
!This program is used to analyze the coordinated water molecules around an ion

PROGRAM cowan
  IMPLICIT NONE
  INTEGER, PARAMETER :: REQUIRED_NUM_PAR = 3
  INTEGER, PARAMETER :: input_fileid = 10
  CHARACTER(LEN=128) :: input_filename
  INTEGER :: count, i, stat, num_water
  CHARACTER(LEN=128) :: command, usage, arg
  !read the parameters from command line arguments
  count = COMMAND_ARGUMENT_COUNT()
  call GET_COMMAND_ARGUMENT(NUMBER=0, VALUE=command)
  usage = "Usage: " // TRIM(ADJUSTL(command)) // " input_filename &
          &-n <# of H2O>"
  
  if (count < REQUIRED_NUM_PAR) then
     write(*,"(A,I1,A)") "At least ", REQUIRED_NUM_PAR, " arguments are needed."
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
     if (stat /= 0) then
        write(*,*) "Something wrong while reading arguments!"
        write(*,*) usage
        call EXIT(1)
     end if
     if (arg == '-n') then
        i = i + 1
        call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=stat)
        read(arg,*) num_water
     else
        write(*,*) "Undefined argument: ", arg
        write(*,*) usage
        call EXIT(1)
     end if
     i = i + 1
  end do
  
END PROGRAM cowan
