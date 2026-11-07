module umwm_cli
  ! Command-line interface of the umwm executable,
  !   umwm CONFIG [-o OUTPUT] [-i INPUT]
  ! which runs the configuration of CONFIG, a directory with the main.nml
  ! and spectrum.nml namelists, or a main namelist file whose directory holds
  ! the others; reads the input files from INPUT (CONFIG/input by default),
  ! where it prepares the missing ones of a configuration with a data source,
  ! such as ERA5; and writes the output and restart files, and copies of the
  ! namelists, to OUTPUT (the current directory by default).

  use datetime_module, only: datetime, strptime, timedelta
  use iso_c_binding, only: c_associated, c_char, c_int, c_null_char, c_ptr
  use iso_fortran_env, only: output_unit
  use umwm_config, only: config_type
  use umwm_env, only: mpisize, nproc
  use umwm_grid, only: input_grid_problem
#ifdef MPI
  use mpi
#endif

  implicit none

  private
  public :: parse_command_line, prepare_inputs, first_pending_input, write_namelist_copies

  character(*), parameter :: main_namelist = 'main.nml'
  character(*), parameter :: copy_prefix = 'umwmout_' ! of the namelist copies in OUTPUT
  character(*), parameter :: cases_index = 'cases.txt'
  integer, parameter :: max_cases = 100

  type :: case_type
    character(64) :: name = ''
    character(256) :: description = ''
  end type case_type

  interface
    type(c_ptr) function c_realpath(path, resolved) bind(c, name='realpath')
      import :: c_char, c_ptr
      character(kind=c_char), intent(in) :: path(*)
      character(kind=c_char), intent(out) :: resolved(*)
    end function c_realpath
    integer(c_int) function c_access(path, mode) bind(c, name='access')
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value :: mode
    end function c_access
    integer(c_int) function c_mkdir(path, mode) bind(c, name='mkdir')
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value :: mode
    end function c_mkdir
  end interface

contains

  subroutine parse_command_line(namelist, input_dir, output_dir, run)
    ! Sets namelist to the main namelist of the configuration given as the
    ! argument, CONFIG/main.nml of a directory or the file itself; input_dir
    ! to the directory of the input files, given with -i or else input/ next
    ! to the main namelist; and output_dir to the directory of the output and
    ! restart files, given with -o or else the current directory, which is
    ! created if it does not exist. run is .false. when the program should
    ! stop without running, after showing the help or the list of cases.
    ! Invalid arguments, a missing configuration, or an output directory that
    ! cannot be written to stop the program with an error.
    character(:), allocatable, intent(out) :: namelist, input_dir, output_dir
    logical, intent(out) :: run

    character(:), allocatable :: argument, config
    integer :: n

    run = .true.
    n = 1
    do while (n <= command_argument_count())
      argument = get_argument(n)
      select case (argument)
      case ('-h', '--help')
        call print_help()
        run = .false.
        return
      case ('-l', '--list-cases')
        call print_cases()
        run = .false.
        return
      case ('-o', '--output')
        if (n == command_argument_count()) call fail('option ' // argument // ' requires a directory')
        n = n + 1
        output_dir = get_argument(n)
      case ('-i', '--input')
        if (n == command_argument_count()) call fail('option ' // argument // ' requires a directory')
        n = n + 1
        input_dir = get_argument(n)
      case default
        if (index(argument, '--output=') == 1) then
          output_dir = argument(len('--output=')+1:)
        else if (index(argument, '--input=') == 1) then
          input_dir = argument(len('--input=')+1:)
        else if (index(argument, '-') == 1) then
          call fail('unknown option ' // argument // '; see umwm --help')
        else if (allocated(config)) then
          call fail('more than one configuration given: ' // config // ' and ' // argument &
                    // '; see umwm --help')
        else
          config = argument
        end if
      end select
      n = n + 1
    end do

    if (.not. allocated(config)) call fail('no configuration given; usage: umwm CONFIG ' &
                                           // '[-o OUTPUT] [-i INPUT], see umwm --help')
    if (len(config) == 0) call fail('the configuration is an empty path')
    namelist = config_namelist(config)
    if (.not. allocated(input_dir)) input_dir = namelist(1:index(namelist, '/', back=.true.)) // 'input'
    if (.not. allocated(output_dir)) output_dir = '.'
    if (len(input_dir) == 0) call fail('the input directory is an empty path')
    if (len(output_dir) == 0) call fail('the output directory is an empty path')
    call prepare_output_directory(output_dir)

  end subroutine parse_command_line


  function config_namelist(config) result(namelist)
    ! The main namelist of the configuration config: config/main.nml if it
    ! is a directory, or the file itself.
    character(*), intent(in) :: config
    character(:), allocatable :: namelist
    character(:), allocatable :: dir
    integer :: n

    n = len(config)
    do while (n > 1 .and. config(n:n) == '/')
      n = n - 1
    end do
    dir = config(1:n)
    if (path_exists(dir // '/.')) then
      namelist = dir // '/' // main_namelist
      if (.not. file_exists(namelist)) call fail('no ' // main_namelist // ' in the configuration ' &
                                                 // 'directory ' // dir)
    else if (path_exists(config)) then
      namelist = config
    else
      call fail('configuration not found: ' // config)
    end if

  end function config_namelist


  function examples_dir() result(dir)
    ! The directory of the example cases: $UMWM_EXAMPLES if set, otherwise
    ! examples/ in the current directory or next to the executable, with
    ! symbolic links resolved.
    character(:), allocatable :: dir
    integer :: length, stat

    call get_environment_variable('UMWM_EXAMPLES', length=length, status=stat)
    if (stat == 0 .and. length > 0) then
      allocate(character(length) :: dir)
      call get_environment_variable('UMWM_EXAMPLES', value=dir)
      if (.not. file_exists(dir // '/' // cases_index)) &
        call fail('no ' // cases_index // ' in UMWM_EXAMPLES=' // dir)
      return
    end if

    dir = 'examples'
    if (file_exists(dir // '/' // cases_index)) return

    dir = executable_dir()
    if (len(dir) > 0) then
      dir = dir // 'examples'
      if (file_exists(dir // '/' // cases_index)) return
    end if

    call fail('example cases not found; run from the UMWM directory, ' &
              // 'or set UMWM_EXAMPLES to the examples directory')

  end function examples_dir


  function executable_dir() result(dir)
    ! The directory of the executable, ending with /, with symbolic links
    ! resolved: that of the program path it was started with, or of the
    ! first directory in PATH that holds a program started by its name.
    ! Empty if the program is not found.
    character(:), allocatable :: dir
    character(:), allocatable :: path, search, entry
    integer :: length, stat, start, colon

    dir = ''
    path = get_argument(0)
    if (index(path, '/') == 0) then
      call get_environment_variable('PATH', length=length, status=stat)
      if (stat /= 0 .or. length == 0) return
      allocate(character(length) :: search)
      call get_environment_variable('PATH', value=search)
      start = 1
      do
        colon = index(search(start:), ':')
        if (colon == 0) then
          entry = search(start:)
        else
          entry = search(start:start+colon-2)
        end if
        if (len(entry) == 0) entry = '.'
        if (file_exists(entry // '/' // path)) exit
        if (colon == 0) return
        start = start + colon
      end do
      path = entry // '/' // path
    end if
    path = real_path(path)
    dir = path(1:index(path, '/', back=.true.))

  end function executable_dir


  function real_path(path) result(res)
    ! The absolute path with symbolic links resolved, by the POSIX realpath,
    ! or path itself if it cannot be resolved.
    character(*), intent(in) :: path
    character(:), allocatable :: res
    character(kind=c_char) :: buffer(4096)
    integer :: n

    if (.not. c_associated(c_realpath(path // c_null_char, buffer))) then
      res = path
      return
    end if
    n = 0
    do while (buffer(n+1) /= c_null_char)
      n = n + 1
    end do
    allocate(character(n) :: res)
    do n = 1, len(res)
      res(n:n) = buffer(n)
    end do

  end function real_path


  subroutine read_cases(dir, cases, ncases)
    ! Reads the names and descriptions of the example cases from the index
    ! of the examples directory, one case per line, the name first.
    character(*), intent(in) :: dir
    type(case_type), intent(out) :: cases(:)
    integer, intent(out) :: ncases
    character(512) :: line
    integer :: unit, stat, space

    ncases = 0
    open(newunit=unit, file=dir // '/' // cases_index, status='old', &
         action='read', iostat=stat)
    if (stat /= 0) call fail('cannot open ' // dir // '/' // cases_index)
    do
      read(unit, '(a)', iostat=stat) line
      if (stat /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      if (ncases == size(cases)) exit
      ncases = ncases + 1
      space = index(line, ' ')
      cases(ncases) % name = line(1:space-1)
      cases(ncases) % description = adjustl(line(space:))
    end do
    close(unit)

  end subroutine read_cases


  subroutine print_cases()
    character(:), allocatable :: dir
    type(case_type) :: cases(max_cases)
    integer :: ncases, n

    dir = examples_dir()
    call read_cases(dir, cases, ncases)
    if (nproc /= 0) return
    write(*, '(a)') 'Example cases in ' // dir // ':'
    do n = 1, ncases
      write(*, '(2x, a24, 1x, a)') cases(n) % name, trim(cases(n) % description)
    end do
    write(*, '(a)') 'Run one with: umwm ' // dir // '/NAME -o OUTPUT'

  end subroutine print_cases


  subroutine print_help()
    if (nproc /= 0) return
    write(*, '(a)') &
      'Usage: umwm CONFIG [-o OUTPUT] [-i INPUT]', &
      '       umwm --list-cases', &
      '', &
      'Runs the University of Miami Wave Model with the configuration CONFIG: a', &
      'directory with the main.nml namelist and the spectrum.nml points of the', &
      'spectrum output, or a main namelist file whose directory holds the others.', &
      '', &
      '  CONFIG              the configuration directory, or its main namelist file', &
      '  -o, --output DIR    write the output and restart files to DIR, created if', &
      '                      needed (default: the current directory), with copies', &
      '                      of the namelists, umwmout_main.nml and the others;', &
      '                      a restart reads its restart file from DIR too', &
      '  -i, --input DIR     read the input files from DIR (default: CONFIG/input)', &
      '  -l, --list-cases    list the example cases', &
      '  -h, --help          show this help', &
      '', &
      'The output directory may be the configuration or input directory: the files', &
      'the model writes there are all named umwmout_*, umwmspc_* or umwmrst_*.', &
      '', &
      'With a data source in the FORCING namelist, such as source = ''era5'', the', &
      'missing input files, and those it prepared for another grid, are prepared in', &
      'the input directory before the run, with the umwm Python package:', &
      '$UMWM_PYTHON -m umwm.cli forcing SOURCE CONFIG -i INPUT (python3 by default).', &
      '', &
      'The example cases are in $UMWM_EXAMPLES if set, otherwise in examples/ in the', &
      'current directory or next to the executable; run one with, for example,', &
      'umwm examples/jonswap -o runs/jonswap.'

  end subroutine print_help


  subroutine prepare_output_directory(dir)
    ! Creates the output directory dir and its parents if they do not exist,
    ! and stops with an error if it cannot be written to. The root process
    ! does it and tells the others, so that all of them stop.
    character(*), intent(in) :: dir
    integer :: unit, stat, n
#ifdef MPI
    integer :: ierr
#endif

    stat = 0
    if (nproc == 0) then
      if (.not. path_exists(dir // '/.')) then
        do n = 2, len(dir)
          if (dir(n:n) == '/' .and. dir(n-1:n-1) /= '/') stat = c_mkdir(dir(1:n-1) // c_null_char, &
                                                                       511_c_int) ! 0777, less the umask
        end do
        stat = c_mkdir(dir // c_null_char, 511_c_int) ! 0777, less the umask
      end if
      open(newunit=unit, file=dir // '/.umwm_write_test', status='replace', iostat=stat)
      if (stat == 0) close(unit, status='delete')
    end if
#ifdef MPI
    call mpi_bcast(stat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
    if (stat /= 0) call fail('cannot create or write to the output directory ' // dir)

  end subroutine prepare_output_directory


  subroutine write_namelist_copies(namelist, output_dir)
    ! Copies the main namelist, and the spectrum.nml next to it, into
    ! output_dir as umwmout_main.nml and umwmout_spectrum.nml, to record the
    ! configuration of the output, and removes a copy there of a namelist
    ! that the configuration does not have. Their names differ from those of a configuration, so that the
    ! output directory may be the configuration directory; a namelist that
    ! is already the copy is left as it is. Only the root process writes.
    use umwm_config, only: sibling_namelist
    character(*), intent(in) :: namelist, output_dir
    character(*), parameter :: names(2) = [character(12) :: 'main.nml', 'spectrum.nml']
    character(:), allocatable :: source, target
    integer :: n

    if (nproc /= 0) return
    do n = 1, size(names)
      if (n == 1) then
        source = namelist
      else
        source = sibling_namelist(namelist, trim(names(n)))
      end if
      target = output_dir // '/' // copy_prefix // trim(names(n))
      if (file_exists(source)) then
        if (real_path(source) /= real_path(target)) call copy_file(source, target)
      else if (file_exists(target)) then
        call delete_file(target)
      end if
    end do

  end subroutine write_namelist_copies


  subroutine copy_file(source, target)
    ! Copies the file source to target, replacing it.
    character(*), intent(in) :: source, target
    character(:), allocatable :: content
    integer :: unit, stat, size_bytes

    inquire(file=source, size=size_bytes)
    if (size_bytes < 0) call fail('cannot read ' // source)
    allocate(character(size_bytes) :: content)
    open(newunit=unit, file=source, access='stream', form='unformatted', status='old', &
         action='read', iostat=stat)
    if (stat /= 0) call fail('cannot read ' // source)
    if (size_bytes > 0) read(unit, iostat=stat) content
    close(unit)
    if (stat /= 0) call fail('cannot read ' // source)
    open(newunit=unit, file=target, access='stream', form='unformatted', status='replace', &
         action='write', iostat=stat)
    if (stat /= 0) call fail('cannot write ' // target)
    write(unit, iostat=stat) content
    close(unit)
    if (stat /= 0) call fail('cannot write ' // target)

  end subroutine copy_file


  subroutine delete_file(path)
    character(*), intent(in) :: path
    integer :: unit, stat
    open(newunit=unit, file=path, status='old', iostat=stat)
    if (stat == 0) close(unit, status='delete')
  end subroutine delete_file


  function get_argument(n) result(argument)
    integer, intent(in) :: n
    character(:), allocatable :: argument
    integer :: length

    call get_command_argument(n, length=length)
    allocate(character(length) :: argument)
    call get_command_argument(n, value=argument)

  end function get_argument


  logical function file_exists(path)
    character(*), intent(in) :: path
    inquire(file=path, exist=file_exists)
  end function file_exists


  logical function path_exists(path)
    ! Whether a file or directory exists at path, by the POSIX access.
    character(*), intent(in) :: path
    path_exists = c_access(path // c_null_char, 0_c_int) == 0
  end function path_exists


  subroutine prepare_inputs(namelist, input_dir)
    ! Prepares the input files that the configuration of namelist reads and
    ! that are missing from input_dir, or that are for another grid, from the
    ! data source of its forcing (source in the FORCING namelist, such as
    ! 'era5'), with the umwm Python package: the root process runs
    !   $UMWM_PYTHON -m umwm.cli forcing SOURCE NAMELIST -i INPUT_DIR
    ! with python3 by default, which prepares again the files for another
    ! grid that it prepared, and stops at those it did not. Does nothing
    ! without a data source or files to prepare, and stops with an error if
    ! they cannot be prepared, or in a parallel run, which does not start
    ! Python: its input files are prepared beforehand. An invalid
    ! configuration is left to the domain to report.
    character(*), intent(in) :: namelist, input_dir
    type(config_type) :: config
    character(:), allocatable :: pending, reason, source, command
    character(256) :: message
    integer :: stat, exitstat, cmdstat
#ifdef MPI
    integer :: ierr
#endif

    config = config_type(namelist)
    if (.not. config % read_ok) return
    if (.not. config % has_source()) return
    if (.not. config % validate(stop_on_error=.false., rank=-1)) return
    config % input_dir = input_dir
    pending = first_pending_input(config, reason)
    if (len(pending) == 0) return

    source = trim(adjustl(config % source))
    command = 'umwm forcing ' // shell_quote(source) // ' ' // shell_quote(namelist) &
      // ' -i ' // shell_quote(input_dir)
    if (mpisize > 1) call fail('input file ' // pending // ' ' // reason // '. A parallel run ' &
      // 'does not prepare its input files; prepare them beforehand with: ' // command)

    stat = 0
    if (nproc == 0) then
      command = environment('UMWM_PYTHON', 'python3') // ' -m umwm.cli forcing ' &
        // shell_quote(source) // ' ' // shell_quote(namelist) // ' -i ' &
        // shell_quote(input_dir)
      write(*, '(a)') 'umwm: input file ' // pending // ' ' // reason // '; preparing the ' &
        // 'input files from ' // source // ' with:', '  ' // command
      flush(output_unit)
      message = ''
      call execute_command_line(command, exitstat=exitstat, cmdstat=cmdstat, cmdmsg=message)
      if (cmdstat /= 0) then
        stat = 1
        write(*, '(a)') 'umwm: ' // trim(message)
      else if (exitstat /= 0) then
        stat = 1
      end if
    end if
#ifdef MPI
    call mpi_bcast(stat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
    if (stat /= 0) call fail('the input files could not be prepared from ' // source // '; see ' &
      // 'the error above. Their preparation needs the umwm Python package, with the ' &
      // 'dependencies of the data source (pip install ".[' // source // ']" in the UMWM ' &
      // 'repository): set UMWM_PYTHON to the Python interpreter that has it, or prepare the ' &
      // 'input files beforehand with: umwm forcing ' // source // ' ' // namelist // ' -i ' &
      // input_dir)

    pending = first_pending_input(config, reason)
    if (len(pending) > 0) call fail('input file ' // pending // ' ' // reason // &
                                    ', after the preparation of the input files from ' // source)

  end subroutine prepare_inputs


  function first_pending_input(config, reason) result(path)
    ! The first input file that a valid configuration reads that does not
    ! exist, or that is for another grid (see input_grid_problem of
    ! umwm_grid): the grid and bathymetry, and the forcing files of the start
    ! time and of the end of every global time step; empty if there is none.
    ! reason tells why it is pending. Restart files are not input files.
    type(config_type), intent(in) :: config
    character(:), allocatable, intent(out) :: reason
    character(:), allocatable :: path
    type(datetime) :: time, stop_time

    path = ''
    reason = ''
    if (config % topofromfile) then
      if (pending(config % input_path('umwm.gridtopo'))) return
    else if (config % gridfromfile) then
      if (file_exists(config % input_path('umwm.gridtopo'))) then
        if (pending(config % input_path('umwm.gridtopo'))) return
      else
        if (pending(config % input_path('umwm.grid'))) return
      end if
    end if
    if (config % forcing_from_file()) then
      if (config % steady) then
        if (pending(config % input_path('umwmin.nc'))) return
      else
        time = strptime(config % starttimestr, '%Y-%m-%d %H:%M:%S')
        stop_time = strptime(config % stoptimestr, '%Y-%m-%d %H:%M:%S')
        do
          if (pending(config % input_path('umwmin_' // time % strftime('%Y-%m-%d_%H:%M:%S') &
                                          // '.nc'))) return
          if (.not. time < stop_time) exit
          time = time + timedelta(seconds=nint(config % dtg))
        end do
      end if
    end if
    path = ''
    reason = ''

  contains

    logical function pending(file)
      ! Whether file is missing or for another grid, setting path and reason.
      character(*), intent(in) :: file
      path = file
      if (.not. file_exists(file)) then
        reason = 'is missing'
      else
        reason = input_grid_problem(config, file)
        if (len(reason) > 0) reason = 'is for another grid: ' // reason
      end if
      pending = len(reason) > 0
    end function pending

  end function first_pending_input


  function environment(name, default) result(value)
    ! The value of the environment variable name, or default if it is not
    ! set or empty.
    character(*), intent(in) :: name, default
    character(:), allocatable :: value
    integer :: length, stat
    call get_environment_variable(name, length=length, status=stat)
    if (stat /= 0 .or. length == 0) then
      value = default
      return
    end if
    allocate(character(length) :: value)
    call get_environment_variable(name, value=value)
  end function environment


  function shell_quote(text) result(quoted)
    ! text quoted for the POSIX shell, between single quotes.
    character(*), intent(in) :: text
    character(:), allocatable :: quoted
    integer :: n
    quoted = "'"
    do n = 1, len(text)
      if (text(n:n) == "'") then
        quoted = quoted // "'\''"
      else
        quoted = quoted // text(n:n)
      end if
    end do
    quoted = quoted // "'"
  end function shell_quote


  subroutine fail(message)
    character(*), intent(in) :: message
    if (nproc == 0) then
      write(*, '(a)') 'umwm: error: ' // message
      flush(output_unit)
    end if
    stop 1
  end subroutine fail

end module umwm_cli
