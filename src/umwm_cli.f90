module umwm_cli
  ! Command-line interface of the umwm executable: selects the main namelist
  ! to run, given as a path or as the name of an example case of the
  ! examples gallery, and checks the run directory.

  use iso_c_binding, only: c_associated, c_char, c_null_char, c_ptr
  use iso_fortran_env, only: output_unit
  use umwm_env, only: nproc
#ifdef MPI
  use mpi
#endif

  implicit none

  private
  public :: parse_command_line

  character(*), parameter :: default_namelist = 'namelists/main.nml'
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
  end interface

contains

  subroutine parse_command_line(namelist, input_dir, run)
    ! Sets namelist to the main namelist to run: the path given as the
    ! argument, the main.nml of the example case given with --case, or
    ! namelists/main.nml by default; and input_dir to the directory of the
    ! input files: input/ of the case, or input/ of the current directory.
    ! run is .false. when the program should stop without running, after
    ! showing the help or the list of cases. Invalid arguments, an unknown
    ! case, or a missing output directory stop the program with an error.
    character(:), allocatable, intent(out) :: namelist, input_dir
    logical, intent(out) :: run

    character(:), allocatable :: argument, case_name
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
      case ('-c', '--case')
        if (n == command_argument_count()) call fail('option ' // argument // ' requires a case name')
        n = n + 1
        case_name = get_argument(n)
      case default
        if (index(argument, '--case=') == 1) then
          case_name = argument(len('--case=')+1:)
        else if (index(argument, '-') == 1) then
          call fail('unknown option ' // argument // '; see umwm --help')
        else if (allocated(namelist)) then
          call fail('more than one namelist given; see umwm --help')
        else
          namelist = argument
        end if
      end select
      n = n + 1
    end do

    input_dir = 'input'
    if (allocated(case_name)) then
      if (allocated(namelist)) call fail('give either a namelist or --case, not both')
      namelist = case_namelist(case_name)
      input_dir = namelist(1:index(namelist, '/', back=.true.)) // 'input'
    else if (.not. allocated(namelist)) then
      namelist = default_namelist
    end if

    if (.not. file_exists(namelist)) call fail('namelist file not found: ' // namelist)
    call check_output_directory()

  end subroutine parse_command_line


  function case_namelist(name) result(namelist)
    ! The main namelist of the example case of the given name, or of the
    ! only case whose name starts with it.
    character(*), intent(in) :: name
    character(:), allocatable :: namelist
    character(:), allocatable :: dir
    type(case_type) :: cases(max_cases)
    integer :: ncases, n, match, matches

    dir = examples_dir()
    call read_cases(dir, cases, ncases)

    match = 0
    matches = 0
    do n = 1, ncases
      if (trim(cases(n) % name) == name) then
        match = n
        matches = 1
        exit
      end if
      if (index(cases(n) % name, name) == 1) then
        match = n
        matches = matches + 1
      end if
    end do

    if (matches /= 1) then
      if (nproc == 0) then
        if (matches == 0) then
          write(*, '(a)') 'umwm: error: unknown case ' // name // '; the example cases are:'
        else
          write(*, '(a)') 'umwm: error: ' // name // ' matches more than one case:'
        end if
        do n = 1, ncases
          if (matches == 0 .or. index(cases(n) % name, name) == 1) &
            write(*, '(2x, a)') trim(cases(n) % name)
        end do
        flush(output_unit)
      end if
      stop 1
    end if

    namelist = dir // '/' // trim(cases(match) % name) // '/main.nml'
    if (nproc == 0) write(*, '(a)') 'umwm: running example case ' &
      // trim(cases(match) % name) // ' from ' // namelist

  end function case_namelist


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
    write(*, '(a)') 'Run one in the current directory with: umwm --case NAME'

  end subroutine print_cases


  subroutine print_help()
    if (nproc /= 0) return
    write(*, '(a)') &
      'Usage: umwm [NAMELIST]', &
      '       umwm --case NAME', &
      '       umwm --list-cases', &
      '', &
      'Runs the University of Miami Wave Model in the current directory, which holds', &
      'the input/, output/ and restart/ directories of the run; an example case', &
      'reads its input files, if it has any, from its own input/ directory.', &
      '', &
      '  NAMELIST          the main namelist to run (default: ' // default_namelist // ')', &
      '  -c, --case NAME   run the example case NAME; a unique prefix of the name', &
      '                    is enough, such as --case titan', &
      '  -l, --list-cases  list the example cases', &
      '  -h, --help        show this help', &
      '', &
      'The example cases are read from $UMWM_EXAMPLES if set, otherwise from', &
      'examples/ in the current directory or next to the executable.'

  end subroutine print_help


  subroutine check_output_directory()
    ! Stops with an error if the output directory cannot be written to. The
    ! root process tests it and tells the others, so that all of them stop.
    integer :: unit, stat
#ifdef MPI
    integer :: ierr
#endif

    stat = 0
    if (nproc == 0) then
      open(newunit=unit, file='output/.umwm_write_test', status='replace', iostat=stat)
      if (stat == 0) close(unit, status='delete')
    end if
#ifdef MPI
    call mpi_bcast(stat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
    if (stat /= 0) call fail('cannot write to output/ in the current directory; ' &
                             // 'create the run directory with mkdir -p output restart')

  end subroutine check_output_directory


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


  subroutine fail(message)
    character(*), intent(in) :: message
    if (nproc == 0) then
      write(*, '(a)') 'umwm: error: ' // message
      flush(output_unit)
    end if
    stop 1
  end subroutine fail

end module umwm_cli
