MODULE filehandle
  CHARACTER(LEN=5), PARAMETER :: BEGIN_CHAR = "BEGIN"
  CHARACTER(LEN=3), PARAMETER :: END_CHAR = "END"

  INTEGER, PARAMETER :: ONE_LINE_LEN = 128
  CHARACTER(LEN=6), PARAMETER :: ONE_LINE = "(A128)"

  LOGICAL, SAVE, PRIVATE :: idsetted = .false.
  INTEGER, SAVE, PRIVATE :: setid = 0

  CHARACTER(LEN=50), SAVE, PRIVATE :: comment_symbol = ""

CONTAINS
!--- SUBROUTINE List ---
! fopen(unit,file,[status],[position],[sfxnum])
!    [status]: the same usage as STATUS in open()
!    [position]: the same usage as POSITION in open()
!    [sfxnum]: suffix number to the filename
! fsetid(unit)
! fread_begin()
! fread_end()
! fclose([unit])
!
! fread(string)
!
!--- FUNCTION List
! int2str(num)
!-----------------------

  SUBROUTINE fopen(unit,file,status,position,sfxnum)  !Open a file and get id
    IMPLICIT NONE
    LOGICAL alive, opened_f, opened_u
    INTEGER, INTENT(OUT) :: unit
    CHARACTER(LEN=*), INTENT(IN) :: file  !filename string
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: status
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: position
    INTEGER, INTENT(IN), OPTIONAL :: sfxnum !suffix number of filename
    CHARACTER(LEN=100):: filename
    !len(filename) = len(file) + len(string(sfxnum))
    CHARACTER(LEN=7) :: stat
    CHARACTER(LEN=6) :: pos
    CHARACTER(LEN=7), PARAMETER :: DEFFAULT_STAT = "UNKNOWN"
    CHARACTER(LEN=6), PARAMETER :: DEFFAULT_POS = "ASIS"
    INTEGER :: iostat
    unit = 100
    if( present(status) ) then
      stat = status
    else
      stat = DEFFAULT_STAT
    end if

    if( present(position) ) then
      pos = position
    else
      pos = DEFFAULT_POS
    end if

    if( present(sfxnum) ) then
      filename = trim(file)//int2str(sfxnum)
    else
      filename = file
    end if

    inquire( FILE = filename, EXIST = alive, OPENED = opened_f )
    if( ( stat == 'OLD' .or. stat == 'old' ) .and. (.not. alive) ) then
      write(*,*) trim(filename), " doesn't exist."
      call EXIT(1)
    else if( opened_f ) then
      write(*,*) trim(filename), " already opened."
      call EXIT(1)
    end if

    inquire( UNIT = unit, OPENED = opened_u )
    if( opened_u ) then
      do
        unit = unit + 1
        inquire( UNIT = unit, OPENED = opened_u )
        if( .not. opened_u )  exit
      end do
    end if

    open( UNIT=unit, FILE=filename, STATUS=stat, IOSTAT=iostat, POSITION=pos )
    if( iostat /= 0 ) then
      write(*,*) "Failed to open '",trim(filename),"' file."
      call EXIT(1)
    end if
  END SUBROUTINE fopen

  SUBROUTINE fsetid(unit)  !Set currently id in use
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    LOGICAL opened_u
    inquire( UNIT = unit, OPENED = opened_u )
    if( .not. opened_u ) then
      write(*,*) "unit=",unit," not opened yet, please call fopen(name,id) first."
      call EXIT(1)
    else
      setid = unit
      idsetted = .true.
    end if
  END SUBROUTINE fsetid

  SUBROUTINE fread_begin()
    IMPLICIT NONE
    CHARACTER(LEN=ONE_LINE_LEN) :: line
    INTEGER :: iostat

    if( .not. idsetted ) then
      write(*,*) "fread_begin: No unit is in use. Please call fsetid(id) first."
      call EXIT(1)
    end if

    do  !loop until it reads begin_tag
      read( setid, ONE_LINE, IOSTAT=iostat ) line
      if( iostat /= 0 ) then
        write(*,*) 'Reading file error, lacking of "',BEGIN_CHAR,'".'
        call EXIT(1)
      else if( trim( adjustl(line) ) == BEGIN_CHAR ) then
             !get rid of the space characters in front and end
        exit
      end if
    end do
  END SUBROUTINE fread_begin

  SUBROUTINE fread_end()
    IMPLICIT NONE
    CHARACTER(LEN=ONE_LINE_LEN) :: line
    INTEGER :: iostat
    if( .not. idsetted ) then
      write(*,*) "fread_end: No id is in use. Please call fsetid(id) first."
      stop
    end if

    read( setid, ONE_LINE, IOSTAT=iostat ) line
    if( iostat /= 0 .or. trim( adjustl(line) ) /= END_CHAR ) then
      write(*,*) 'Reading file error, lacking of "',END_CHAR,'".'
      stop
    end if
  END SUBROUTINE fread_end

  SUBROUTINE fclose(unit)
    IMPLICIT NONE
    INTEGER, OPTIONAL, INTENT(IN) :: unit
    if( present(unit) ) then
      if( unit == setid ) then
        idsetted = .false.
      end if
      close(unit)
    else if( idsetted ) then
      idsetted = .false.
      close(setid)
    end if
  END SUBROUTINE fclose

  PURE FUNCTION int2str(num)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: num
    CHARACTER(LEN = floor( log10( real(num) ) ) + 1 ) :: int2str
    INTEGER :: i,dum,r
    dum = num
    do i = len(int2str),1,-1
      r = mod(dum,10)
      int2str(i:i) = ACHAR(r+48)
      dum = dum/10
    end do
  END FUNCTION int2str

  SUBROUTINE fset_comment_symbol( sym_str )
    IMPLICIT NONE
    CHARACTER(LEN=*) :: sym_str
    comment_symbol = TRIM(ADJUSTL(sym_str))
  END SUBROUTINE fset_comment_symbol

  FUNCTION ffind(unit, target_str, mode)
    IMPLICIT NONE
    LOGICAL ffind
    INTEGER, OPTIONAL :: unit
    CHARACTER(LEN=*) :: target_str
    CHARACTER(LEN=*), OPTIONAL :: mode ![WORD][LINE]
    CHARACTER(LEN=*), PARAMETER :: MODE_S='STRING', MODE_L='LINE'
    !MODE_S: Search the first line CONTAINS the target string
    !MODE_L: Search the first line IS the target string
    INTEGER :: runit ! real unit in use
    CHARACTER(LEN=10) :: rmode ! real mode in use
    CHARACTER(LEN=ONE_LINE_LEN) :: line
    INTEGER :: com_pos ! comment position in a line
    INTEGER :: iostat

    if( present(unit) ) then
      runit = unit
    else if( idsetted ) then
      runit = setid
    else
      write(*,*) "ffind: No unit is in use. &
              &Please call fsetid(id) first or specify a unit."
      call EXIT(1)
    end if

    if( present(mode) ) then
      rmode = mode
    else
      rmode = MODE_S
    end if

    do
      read(unit, ONE_LINE, IOSTAT=iostat) line
      if( iostat > 0 ) then
        write(*,*) 'ffind: Reading file error'
        call EXIT(1)
      else if( iostat < 0 ) then
        ffind = .false.
        exit
      end if
      line = ADJUSTL(line)
      com_pos = SCAN( line, comment_symbol )
      if( com_pos == 1 ) then ! The whole line is a comment
        cycle
      else if( com_pos /= 0 ) then
        line = line(1:com_pos-1) ! Ignore those after comment symbol
      end if

      select case(rmode)
        case(MODE_S)
          if( INDEX(line, target_str) /= 0 ) then !Target discovered
            exit
          end if
        case(MODE_L)
          if( line == target_str ) then !Target discovered
            exit
          end if
      end select
    end do
    ffind = .true.
    RETURN
  END FUNCTION ffind

END MODULE filehandle
