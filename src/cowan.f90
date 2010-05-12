!Coordinated Water Analysis (cowan)
!This program is used to analyze the coordinated water molecules around an ion

PROGRAM cowan
  IMPLICIT NONE
  INTEGER, PARAMETER :: REQUIRED_NUM_PAR = 3
  INTEGER, PARAMETER :: input_fileid = 10
  CHARACTER(LEN=128) :: input_filename
  INTEGER, PARAMETER :: output_fileid = 11
  CHARACTER(LEN=128) :: output_filename
  INTEGER :: n, i, j, k, stat
  CHARACTER(LEN=128) :: command, usage, arg, line
  CHARACTER(LEN=2) :: out_fmt
  CHARACTER(LEN=8) :: str_timestep, atom_name
  CHARACTER, PARAMETER :: OW_IN = '*', OW_OUT = '_'
  CHARACTER :: OW_state
  CHARACTER(LEN=2), PARAMETER :: OW_name = 'OW'
  INTEGER :: timestep, num_water, num_atoms, num_frames, init_timestep
  INTEGER :: delta_timestep, atom_index, num_non_ion_water_atoms
  INTEGER :: init_OW_index, num_records
  REAL(KIND=8) :: hydration_radius, hydration_radius_square, dist_square
  REAL(KIND=8) :: tube_half_length, steptime
  REAL(KIND=8), DIMENSION(3) :: pos_ion, pos_tun1, pos_water, cell_vector
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: pos_ion_table
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: hydration_table
  INTEGER, DIMENSION(:), ALLOCATABLE :: hydration_sum_table
  
  !-- read parameters from command line arguments --!
  n = COMMAND_ARGUMENT_COUNT()
  call GET_COMMAND_ARGUMENT(NUMBER=0, VALUE=command)
  usage = "Usage: " // TRIM(ADJUSTL(command)) // " input_filename &
          &-r <hydration radius>"
  
  if (n < REQUIRED_NUM_PAR) then
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
  do while (i <= n)
     call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=stat)
     select case (arg)
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
     case default
        write(*,*) "Undefined argument: ", arg
        write(*,*) usage
        call EXIT(1)
     end select
     i = i + 1
  end do
  !-- end of reading parameters from command line arguments --!

  open(UNIT=input_fileid, FILE=input_filename, STATUS='OLD', IOSTAT=stat, ACTION='READ')
  if (stat /= 0) then
     write(*,*) "Error: unable to open file: ", input_filename
     call EXIT(1)
  end if  
  
  call count_num_frames(num_frames, steptime)
  call count_num_water(num_water, init_OW_index)
  
  !output parameters being read for confirmation
  write(*,*) "input filename: ", TRIM(input_filename)
  write(*,*) "number of frames:", num_frames
  write(*,*) "timestep::", steptime
  write(*,*) "number of water molecules:", num_water
  write(*,*) "initial OW index:", init_OW_index
  write(*,*) "hydration radius:", hydration_radius

  !square the hydration_radius for furthur comparison with the dist_square between water and ion
  hydration_radius_square = hydration_radius * hydration_radius

  allocate(hydration_table(num_water,num_frames), STAT=stat)
  if (stat /=0) then
     write(*,*) "Error: failed to allocate hydration_table(num_water, num_frames)!"
     call EXIT(1)
  end if

  allocate(pos_ion_table(num_frames), STAT=stat)
  if (stat /=0) then
     write(*,*) "Error: failed to allocate pos_ion_table(num_frames)!"
     call EXIT(1)
  end if

  do while(.true.)
     read(input_fileid, "(A)", IOSTAT=stat) line
     if (stat < 0) then
        write(*,*) "Error: problem occurs while reading the 1st timestep"
        write(*,*) "       End of file encounters!"
        call EXIT(1)
     end if
     read(line, *, IOSTAT=stat) str_timestep, timestep, num_atoms
     if (stat == 0) then
        if (str_timestep == "timestep") then
           init_timestep = timestep        
           exit
        end if
     end if
  end do
  num_non_ion_water_atoms = num_atoms - num_water*3 - 1

  ! start reading records of each timestep
  do i = 1, num_frames
     if (i > 1) then
        read(input_fileid, *, IOSTAT=stat) str_timestep, timestep
        if (i == 2) then
           if (stat /= 0 .OR. str_timestep /= "timestep") then
              write(*,*) "Error: problem occurs while reading the 2nd timestep"
              call EXIT(1)
           end if
           delta_timestep = timestep - init_timestep
        else
           if (stat /= 0 .OR. str_timestep /= "timestep") then
              write(*,*) "Error: problem occurs while reading timestep:", &
                   &init_timestep+(i-1)*delta_timestep
              call EXIT(1)
           end if
        end if
     end if

     !ignore cell vectors
     read(input_fileid, *) !cell vector
     read(input_fileid, *) !cell vector
     !read z vector for making postscript graph
     read(input_fileid, *) cell_vector !cell vector

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
     pos_ion_table(i) = pos_ion(3)

     if (num_non_ion_water_atoms /= 0) then
        !read the 1st TUN1 atom z position to know the length of the tube
        read(input_fileid, *, IOSTAT=stat) atom_name, atom_index
        if (stat /= 0) then
           write(*,*) "Error: problem occurs while reading TUN1 records at timestep:", timestep
           call EXIT(1)
        end if
        if (atom_name /= "TUN1") then
           write(*,*) "Error: problem occurs while reading TUN1 records at timestep:", timestep
           write(*,*) "       The atom_name should be TUN1, instead of ", atom_name
           call EXIT(1)
        end if
        read(input_fileid, *, IOSTAT=stat) pos_tun1
        if (stat /= 0) then
           write(*,*) "Error: problem occurs while reading TUN1 position at timestep:", timestep
           call EXIT(1)
        end if
        tube_half_length = ABS(pos_tun1(3))

        !skip non-water atoms ()
        do j = 1, num_non_ion_water_atoms-1
           read(input_fileid, *) !skip atom_name, atom_index
           read(input_fileid, *) !skip atom position
        end do
     else
        tube_half_length = -1
     end if

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

  ! make hydration_sum_table
  allocate(hydration_sum_table(num_frames), STAT=stat)
  if (stat /=0) then
     write(*,*) "Error: failed to allocate hydration_sum_table(num_frames)!"
     call EXIT(1)
  end if
  hydration_sum_table = COUNT(hydration_table, 1)

  ! output results
  output_filename = TRIM(ADJUSTL(input_filename)) // ".cowan"
  open(UNIT=output_fileid, FILE=output_filename, STATUS='REPLACE', IOSTAT=stat, ACTION='WRITE')
  if (stat /= 0) then
     write(*,*) "Error: unable to create file: ", output_filename
     call EXIT(1)
  end if

  num_records = 0
  do i = 1, num_water
     if (ANY(hydration_table(i,:))) then
        num_records = num_records + 1
        do j = 0, num_frames
           if (j == 0) then
              write(out_fmt, "(I2)") get_num_digits(num_atoms)
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

  call output_ps()

  deallocate(hydration_table)
  close(output_fileid)
CONTAINS
  SUBROUTINE output_ps()
    IMPLICIT NONE
    INTEGER, PARAMETER :: ps_output_fileid = 12
    CHARACTER(LEN=128) :: ps_output_filename
    CHARACTER(LEN=2) :: out_fmt2, out_fmt3
    INTEGER :: order
    
    ! output results
    ps_output_filename = TRIM(ADJUSTL(input_filename)) // ".ps"
    open(UNIT=ps_output_fileid, FILE=ps_output_filename, STATUS='REPLACE', IOSTAT=stat, ACTION='WRITE')
    if (stat /= 0) then
       write(*,*) "Error: unable to create file: ", ps_output_filename
       call EXIT(1)
    end if

    write(ps_output_fileid, "(&
         &'%!PS-Adobe-3.0'/&
         &'<< /PageSize [ 595 842 ] /ImagingBBox null >> setpagedevice'/&
         &'72 25.4 div dup scale %use mm-units'/&
         &''/&
         &'/pageHeight 297 def'/&
         &'/pageWidth 210 def'/&
         &''/&
         &'/marginLeft pageWidth 15 div def'/&
         &'/marginRight marginLeft def'/&
         &'/marginTop pageHeight 30 div def'/&
         &'/marginBot pageHeight 20 div def'/&
         &''/&
         &'/boundTop pageHeight marginTop sub def'/&
         &'/boundLeft marginLeft def'/&
         &'/boundBot marginBot def'/&
         &''/&
         &)")
    write(ps_output_fileid, "(&
         &'/numRecords ', I10, ' def'/&
         &'/numFrames ', I10, ' def'/&
         &)") num_records, num_frames
    write(ps_output_fileid, "(&
         &'/recordHeight pageHeight 2 div numRecords div def'/&
         &''/&
         &'/vertAxis boundLeft marginLeft add def'/&
         &'/vertAxis2 pageWidth marginRight 2 mul sub def'/&
         &''/&
         &'/waterLabelAlign {vertAxis (00000) stringwidth pop sub} bind def'/&
         &'recordHeight 1.2 div 4 lt {/waterLabelFontSize recordHeight 1.2 div def} '/&
         &'{/waterLabelFontSize 4 def} ifelse '/&
         &'/generalFontSize pageWidth 50 div def'/&
         &''/&
         &'/Helvetica findfont'/&
         &'generalFontSize scalefont'/&
         &'setfont'/&
         &''/&
         &'/minVertTick vertAxis pageWidth 40 div add def'/&
         &'/maxVertTick vertAxis2 pageWidth 40 div sub def'/&
         &'/numVertTicks 11 def'/&
         &'/vertTickSpace maxVertTick minVertTick sub numVertTicks 1 sub div def'/&
         &''/&
         &'/horiAxis1 boundBot def'/&
         &'/horiAxis2 horiAxis1 numRecords 2 add recordHeight mul add pageHeight 60 div add def'/&
         &'/horiAxis3 boundTop horiAxis2 add 2 div def'/&
         &''/&
         &'/minHoriTick horiAxis2 pageHeight 40 div add def'/&
         &'/maxHoriTick horiAxis3 pageHeight 40 div pageHeight 60 div add sub def'/&
         &'/minHoriTickNum ', I2, ' def'/&
         &'/maxHoriTickNum ', I2, ' def'/&
         &)") MINVAL(hydration_sum_table), MAXVAL(hydration_sum_table)
    write(ps_output_fileid, "(&
         &'/numHoriTicks maxHoriTickNum minHoriTickNum sub 1 add def'/&
         &'/horiTickSpace maxHoriTick minHoriTick sub numHoriTicks 1 sub div def'/&
         &''/&
         &'/dataPointSpace maxVertTick minVertTick sub numFrames div def'/&
         &''/&
         &'/HSBoxHeight recordHeight 2 div def'/&
         &'/HSBoxWidth dataPointSpace def'/&
         &'/HSBox {'/&
         &'newpath'/&
         &'moveto'/&
         &'HSBoxWidth 2 div 0 rlineto'/&
         &'0 HSBoxHeight rlineto'/&
         &'HSBoxWidth neg 0 rlineto'/&
         &'0 HSBoxHeight neg rlineto'/&
         &'closepath'/&
         &'fill'/&
         &'} bind def'/&
         &''/&
         &'/vertTickLength pageHeight 200 div def'/&
         &'/vertTick {'/&
         &'  newpath'/&
         &'  moveto'/&
         &'  0 vertTickLength rlineto'/&
         &'  stroke'/&
         &'} bind def'/&
         &''/&
         &'/horiTickLength pageHeight 200 div def'/&
         &'/horiTick {'/&
         &'  newpath'/&
         &'  moveto'/&
         &'  horiTickLength 0 rlineto'/&
         &'  stroke'/&
         &'} bind def'/&
         &''/&
         &'/arrowhead {% stack: s x1 y1, current point: x0 y0 '/&
         &'  gsave '/&
         &'    currentpoint % s x1 y1 x0 y0 '/&
         &'    4 2 roll exch % s x0 y0 y1 x1 '/&
         &'    4 -1 roll exch % s y0 y1 x0 x1 '/&
         &'    sub 3 1 roll % s x1-x2 y0 y1 '/&
         &'    sub exch % s y0-y1 x1-x2 '/&
         &'    atan rotate % rotate over arctan((y0-y1)/(x1-x2)) '/&
         &'    dup scale % scale by factor s'/&
         &'    currentpoint exch 2 add exch moveto'/&
         &'    -5 2 rlineto 1 -2 rlineto -1 -2 rlineto '/&
         &'    closepath fill % fill arrowhead '/&
         &'  grestore '/&
         &'} def'/&
         &'%%%%%% def end %%%%%%'/&
         &''/&
         &'%set suitable linewidth'/&
         &'pageHeight 600 div pageWidth 400 div lt'/&
         &'{pageHeight 600 div} {pageWidth 400 div} ifelse setlinewidth'/&
         &''/&
         &'%set line cap to round cap'/&
         &'1 setlinecap'/&
         &''/&
         &'%draw vertAxis'/&
         &'newpath'/&
         &'vertAxis horiAxis1 moveto'/&
         &'vertAxis horiAxis2 pageHeight 60 div sub lineto'/&
         &'stroke'/&
         &''/&
         &'newpath'/&
         &'vertAxis horiAxis2 moveto'/&
         &'vertAxis horiAxis3 pageHeight 60 div sub lineto'/&
         &'0.6 vertAxis horiAxis2 arrowhead'/&
         &'stroke'/&
         &'%draw horiTick with labels'/&
         &'0 1 numHoriTicks 1 sub {'/&
         &'  dup horiTickSpace mul minHoriTick add vertAxis exch horiTick'/&
         &'  newpath'/&
         &'  dup horiTickSpace mul minHoriTick add pageHeight 200 div sub vertAxis (00) stringwidth pop sub exch moveto'/&
         &'  minHoriTickNum add 1 string cvs show'/&
         &'} for'/&
         &'gsave'/&
         &'  /Helvetica findfont'/&
         &'  generalFontSize 1.3 mul scalefont'/&
         &'  setfont'/&
         &'  newpath'/&
         &'  boundLeft (0) stringwidth pop add horiAxis2 horiAxis3 add 2 div moveto'/&
         &'  (n) show'/&
         &'grestore'/&
         &''/&
         &'newpath'/&
         &'vertAxis horiAxis3 moveto'/&
         &'vertAxis boundTop lineto'/&
         &'0.6 vertAxis horiAxis3 arrowhead'/&
         &'stroke'/&
         &'%draw vertical axis label'/&
         &'gsave'/&
         &'  /Helvetica findfont'/&
         &'  generalFontSize 1.3 mul scalefont'/&
         &'  setfont'/&
         &'  newpath'/&
         &'  boundLeft (0) stringwidth pop add boundTop horiAxis3 add 2 div moveto'/&
         &'  (z) show'/&
         &'grestore'/&
         &''/&
         &'%draw horiAxis1 and arrowhead'/&
         &'newpath'/&
         &'boundLeft horiAxis1 moveto'/&
         &'vertAxis2 horiAxis1 lineto'/&
         &'0.6 vertAxis2 1 sub horiAxis1 arrowhead'/&
         &'stroke'/&
         &)")
    write(ps_output_fileid, "(&    
         &'%draw vertTick with labels'/&
         &'0 1 numVertTicks 1 sub {'/&
         &'  dup vertTickSpace mul minVertTick add horiAxis1 vertTick'/&
         &'  newpath'/&
         &'  dup ',G7.2,' numFrames mul numVertTicks 1 sub div mul'/&
         &'  dup exch 10 lt {pop /temp (0) def} {100 lt {/temp (00) def} {/temp (000) def} ifelse} ifelse'/&
         &'  dup vertTickSpace mul minVertTick add temp stringwidth pop 2 div sub'/&
         &'  horiAxis1 pageHeight 60 div sub moveto'/&
         &'  ',G7.2,' numFrames mul numVertTicks 1 sub div mul cvi 6 string cvs show'/&
         &'%use below if the timestep labels have decimal part.'/&
         &'% ',G7.2,' numFrames mul numVertTicks 1 sub div mul cvr 6 string cvs show'/&
         &'} for'/&
         &)") delta_timestep*steptime, delta_timestep*steptime,&
         &delta_timestep*steptime
    write(ps_output_fileid, "(&
         &'%draw horiAxis2 and arrowhead'/&
         &'newpath'/&
         &'boundLeft horiAxis2 moveto'/&
         &'vertAxis2 horiAxis2 lineto'/&
         &'0.6 vertAxis2 1 sub horiAxis2 arrowhead'/&
         &'stroke'/&
         &'%draw vertTick with labels'/&
         &'0 1 numVertTicks 1 sub {'/&
         &'  dup vertTickSpace mul minVertTick add horiAxis2 vertTick'/&
         &'  newpath'/&
         &'  dup ',G7.2,' numFrames mul numVertTicks 1 sub div mul'/&
         &'  dup exch 10 lt {pop /temp (0) def} {100 lt {/temp (00) def} {/temp (000) def} ifelse} ifelse'/&
         &'  dup vertTickSpace mul minVertTick add temp stringwidth pop 2 div sub'/&
         &'  horiAxis2 pageHeight 60 div sub moveto'/&
         &'  ',G7.2,' numFrames mul numVertTicks 1 sub div mul cvi 6 string cvs show'/&
         &'%use below if the timestep labels have decimal part.'/&
         &'% ',G7.2,' numFrames mul numVertTicks 1 sub div mul cvr 6 string cvs show'/&
         &'} for'/&
         &)") delta_timestep*steptime, delta_timestep*steptime,&
         &delta_timestep*steptime
    write(ps_output_fileid, "(&
         &'%draw horiAxis3 and arrowhead'/&
         &'newpath'/&
         &'boundLeft horiAxis3 moveto'/&
         &'vertAxis2 horiAxis3 lineto'/&
         &'0.6 vertAxis2 1 sub horiAxis3 arrowhead'/&
         &'stroke'/&
         &'%draw vertTick with labels'/&
         &'0 1 numVertTicks 1 sub {'/&
         &'  dup vertTickSpace mul minVertTick add horiAxis3 vertTick'/&
         &'  newpath'/&
         &'  dup ', G7.2,' numFrames mul numVertTicks 1 sub div mul'/&
         &'  dup exch 10 lt {pop /temp (0) def} {100 lt {/temp (00) def} {/temp (000) def} ifelse} ifelse'/&
         &'  dup vertTickSpace mul minVertTick add temp stringwidth pop 2 div sub'/&
         &'  horiAxis3 pageHeight 60 div sub moveto'/&
         &'  ', G7.2, ' numFrames mul numVertTicks 1 sub div mul cvi 6 string cvs show'/&
         &'%use below if the timestep labels have decimal parts.'/&
         &'% ', G7.2, ' numFrames mul numVertTicks 1 sub div mul cvr 6 string cvs show'/&
         &)") delta_timestep*steptime, delta_timestep*steptime,&
         &delta_timestep*steptime
    write(ps_output_fileid, "(&
         &'} for'/&
         &''/&
         &''/&
         &'newpath'/&
         &'vertAxis2 (0) stringwidth pop add horiAxis1 pageHeight 290 div sub moveto'/&
         &'(t/ps) show'/&
         &''/&
         &'newpath'/&
         &'vertAxis2 (0) stringwidth pop add horiAxis2 pageHeight 290 div sub moveto'/&
         &'(t/ps) show'/&
         &''/&
         &'newpath'/&
         &'vertAxis2 (0) stringwidth pop add horiAxis3 pageHeight 290 div sub moveto'/&
         &'(t/ps) show'/&
         &''/&
         &''/&
         &'/Helvetica findfont'/&
         &'waterLabelFontSize scalefont'/&
         &'setfont'/&
         &''/&
         &'%draw water records'/&
         &)")
    write(ps_output_fileid, "('%%%% Fortran loop starts')")
    order = 0
    write(out_fmt3, "(I2)") get_num_digits(num_frames)  
    do i = num_water, 1, -1
       if (ANY(hydration_table(i,:))) then
          order = order + 1
          write(out_fmt2, "(I2)") get_num_digits(init_OW_index + (i-1)*3)  
          write(ps_output_fileid, "(&
               &'newpath'/&
               &'/order ', I"//TRIM(ADJUSTL(out_fmt))//", ' def'/&
               &'waterLabelAlign horiAxis1 recordHeight order mul add moveto'/&
               &'(',I"//TRIM(ADJUSTL(out_fmt2))//",') show'/&
               &)") order, init_OW_index + (i-1)*3
          do j = 1, num_frames
             if (hydration_table(i,j)) then
                write(ps_output_fileid, "('minVertTick HSBoxWidth ', I"//TRIM(ADJUSTL(out_fmt3))//&
                     &",' mul add horiAxis1 recordHeight order mul add HSBox')") j
             end if
          end do
          write(ps_output_fileid,*)
       end if
    end do
    write(ps_output_fileid, "('%%%% Fortran loop ends'/)")
    write(ps_output_fileid, "(&
         &'%draw n'/&
         &'gsave'/&
         &'1 setlinejoin'/&
         &'%0.2 setlinewidth %use it if the lines are too dense'/&
         &'newpath'/&
         &)")
    write(ps_output_fileid, "('%%%% Fortran loop starts')")
    do i = 1, num_frames
       write(ps_output_fileid, "(&
            &'/data ', I2, ' def'/&
            &'/point data minHoriTickNum sub def')") hydration_sum_table(i)
       if (i == 1) then
          write(ps_output_fileid, "(&
               &'minVertTick dataPointSpace ', I5, ' mul add minHoriTick horiTickSpace point mul add moveto'/&
               &)") i
       else
          write(ps_output_fileid, "(&
               &'minVertTick dataPointSpace ', I5, ' mul add minHoriTick horiTickSpace point mul add lineto'/&
               &)") i
       end if
    end do  
    write(ps_output_fileid, "('%%%% Fortran loop ends'/)")
    
    write(ps_output_fileid, "(&
         &'stroke'/&
         &'grestore'//&
         &'%draw tube top and bottom line'/&
         &'/topLimit boundTop def'/&
         &'/botLimit horiAxis3 vertTickLength add def'/&
         &'/graphRange topLimit botLimit sub def'/&
         &)")
    write(ps_output_fileid, "(&
         &'/dataTopLimit ', F6.2, ' def'/&
         &'/dataBotLimit ', F6.2, ' def'/&
         &'/dataRange dataTopLimit dataBotLimit sub def'/&
         &''/&
         &'/changeToGraph {'/&
         &'  dataBotLimit sub dataRange div graphRange mul botLimit add'/&
         &'} def'/&
         &)") cell_vector(3)/2., -cell_vector(3)/2.
    if (tube_half_length > 0) then
       write(ps_output_fileid, "(&
            &'/dataTubeTop ', F6.2, ' def'/&
            &'/dataTubeBot ', F6.2, ' def'/&
            &'/tubeTopAxis dataTubeTop changeToGraph def'/&
            &'/tubeBotAxis dataTubeBot changeToGraph def'/&
            &)") tube_half_length, -tube_half_length
       write(ps_output_fileid, "(&
            &'gsave'/&
            &'  0 setlinecap'/&
            &'  [pageWidth 200 div dup] 0 setdash'/&
            &'  newpath'/&
            &'  vertAxis tubeTopAxis moveto'/&
            &'  vertAxis2 tubeTopAxis lineto'/&
            &'  stroke'/&
            &'  newpath'/&
            &'  vertAxis tubeBotAxis moveto'/&
            &'  vertAxis2 tubeBotAxis lineto'/&
            &'  stroke'/&
            &'grestore'/&
            &''/&
            &'%draw labels'/&
            &'/Helvetica findfont'/&
            &'generalFontSize scalefont'/&
            &'setfont'/&
            &'newpath'/&
            &'/label dataTubeTop 6 string cvs def'/&
            &'vertAxis label stringwidth pop (0) stringwidth pop 2 div add sub tubeTopAxis pageHeight 200 div sub moveto'/&
            &'label show'/&
            &'/label dataTubeBot 6 string cvs def'/&
            &'vertAxis label stringwidth pop (0) stringwidth pop 2 div add sub tubeBotAxis pageHeight 200 div sub moveto'/&
            &'label show'/&
            &'/label (AA) def'/&
            &'vertAxis label stringwidth pop (0) stringwidth pop 2 div add sub topLimit pageHeight 60 div sub moveto'/&
            &'label show'/&
            &)")
    else
        write(ps_output_fileid, "(&
             &'%draw labels'/&
             &'/Helvetica findfont'/&
             &'generalFontSize scalefont'/&
             &'setfont'/&
             &'newpath'/&
             &'/label dataTopLimit 6 string cvs def'/&
             &'vertAxis label stringwidth pop (0) stringwidth pop 2 div add sub &
             &dataTopLimit changeToGraph pageHeight 80 div sub moveto'/&
             &'label show'/&
             &'/label dataBotLimit 6 string cvs def'/&
             &'vertAxis label stringwidth pop (0) stringwidth pop 2 div add sub &
             &dataBotLimit changeToGraph pageHeight 200 div add moveto'/&
             &'label show'/&
             &'/label (AA) def'/&
             &'vertAxis label stringwidth pop (0) stringwidth pop 2 div add sub &
             &topLimit pageHeight 60 div sub moveto'/&
             &'label show'/&
             &)")
    end if
    write(ps_output_fileid, "(&
         &'gsave'/&
         &'1 setlinejoin'/&
         &'newpath'/&
         &)")
    write(ps_output_fileid, "('%%%% Fortran loop starts')")
    do i = 1, num_frames
       if (i == 1) then
          write(ps_output_fileid, "(&
               &'/data ', e20.12, ' def'/&
               &'minVertTick dataPointSpace ', I5, ' mul add data changeToGraph moveto'/&
               &)") pos_ion_table(i), i
       else
          write(ps_output_fileid, "(&
               &'/data ', e20.12, ' def'/&
               &'minVertTick dataPointSpace ', I5, ' mul add data changeToGraph lineto'/&
               &)") pos_ion_table(i), i
       end if
    end do
    write(ps_output_fileid, "('%%%% Fortran loop ends'/)")
    write(ps_output_fileid, "(&
         &'stroke'/&
         &'grestore'/&
         &''/&
         &'showpage'/&
         &)")
  END SUBROUTINE output_ps

  SUBROUTINE count_num_frames(total_num, dt)
    IMPLICIT NONE
    INTEGER :: total_num
    REAL(KIND=8) :: dt
    CHARACTER(LEN=128) :: line
    INTEGER :: dummy
    total_num = 0
    dt = -1.
    do while(.true.)
       read(input_fileid, "(A)", IOSTAT=stat) line
       if (stat < 0) then       !End of file
          exit
       end if
       read(line, *, IOSTAT=stat) str_timestep
       if (stat /= 0) then
          write(*,*) "Error: problem occurs while counting the num_frames"
          call EXIT(1)          
       else if (str_timestep == "timestep") then
          if (dt < 0) then
             read(line, *, IOSTAT=stat) str_timestep, dummy, dummy ,dummy,&
                  & dummy, dt
          end if
          total_num = total_num + 1
       end if
    end do
    REWIND(input_fileid)
  END SUBROUTINE count_num_frames

  SUBROUTINE count_num_water(total_num, init_num)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: total_num, init_num
    CHARACTER(LEN=128) :: line    
    do while(.true.)
       read(input_fileid, "(A)", IOSTAT=stat) line
       if (stat < 0) then       !End of file
          write(*,*) "Error: problem occurs while reading ion records for counting num_water."
          write(*,*) "       There is no water at all!"
          call EXIT(1)
       end if       
       read(line, *, IOSTAT=stat) atom_name
       if (stat /= 0) then
          write(*,*) "Error: problem occurs while reading ion records for counting num_water."
          call EXIT(1)
       else if (atom_name == "OW") then
          read(line, *, IOSTAT=stat) atom_name, atom_index
          init_num = atom_index
          total_num = 1
          read(input_fileid, *, IOSTAT=stat) pos_ion
          exit
       end if
    end do
    do while(.true.)
       do i = 1,2               !ignore HW records
          read(input_fileid, *, IOSTAT=stat) atom_name, atom_index
          if (stat /= 0) then
             write(*,*) "Error: problem occurs while reading ion records for counting num_water"
             call EXIT(1)
          end if
          if (atom_name /= "HW") then
             write(*,*) "Error: problem occurs while reading ion records for counting num_water"
             write(*,*) "atom_name should be HW, instead of "//atom_name
             call EXIT(1)
          end if
          read(input_fileid, *, IOSTAT=stat) pos_ion
       end do
       ! read OW record
       read(input_fileid, *, IOSTAT=stat) atom_name, atom_index
       if (stat /= 0) then
          write(*,*) "Error: problem occurs while reading ion records for counting num_water"
          call EXIT(1)
       end if
       if (atom_name == "OW") then
          total_num = total_num + 1
       else if (atom_name == "timestep") then
          exit
       else
          write(*,*) "Error: problem occurs while reading ion records for counting num_water:"
          write(*,*) "atom_name should be OW, instead of "//atom_name
          call EXIT(1)          
       end if
       read(input_fileid, *, IOSTAT=stat) pos_ion       
    end do
    REWIND(input_fileid)
  END SUBROUTINE count_num_water

  FUNCTION get_num_digits(ii)
    IMPLICIT NONE
    INTEGER :: get_num_digits
    INTEGER, INTENT(IN) :: ii
    get_num_digits = INT(LOG10(REAL(ii))) + 1
    if ( ii > 10.**get_num_digits - 1) then
       get_num_digits = get_num_digits + 1
    end if
  END FUNCTION get_num_digits
END PROGRAM cowan
