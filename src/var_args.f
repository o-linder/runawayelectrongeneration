C######################################################################|
      module var_args
C----------------------------------------------------------------------|
C     This module is used to pass a variable number of arrays to an
C     arbitrary function accepting only one 1-dimensional array.
C----------------------------------------------------------------------|
      use double
      use math, only :
     >  max_a

      implicit none

      public ::
     >  var_args_r,
     >  get_dependent_indices_and_size

      private ::
     >  dp, max_a

      contains

C======================================================================|
      recursive function var_args_r(pntr_s, args, n_elems) result(res)
C----------------------------------------------------------------------|
C     This function can be used in a workflow to pass an array `args` 
C     of arbitrary values to a function accepting only scalar values:
C       var_args_r(args(:,:)) 
C           --> scalar_function_wrapper(args(:)) 
C               --> scalar_function(args(1), args(2), args(3), ...)
C
C     Given a 2-dimensional array `args` of arbitrary values and a
C     1-dimensional array `n_elems` containing the number of elements
C     used in `args` along `dim=1`, `var_args_r` recursively calls 
C     itself with arrays `args` of reduced size until `args` is
C     effectively a 1-dimensional array. In this case, `args` is passed
C     to the function `pntr_s`.
C   
C     Requirements:
C       - size(args,dim=2) = size(n_elems,dim=1)
C       - size(args,dim=1) >= max(n_elems)
C
C     By default, the function `pntr_s` is called with all 1-d
C     combinations of `args`, i.e. values are treated independent of
C     each other. If the number of neighboring elements of `args` along
C     `dim=1` is identical, they will be treated dependent on each
C     other.
C
C     Examples:
C     1. Independent mode:
C        Calling `var_args_r(pntr, reshape((/10._dp, 11._dp, 0._dp/), 
C           (/0._dp, 1._dp, 2._dp/)), (/3, 2/)), (/2, 3/))`
C        with different number of elements in `args` along `dim=1)        
C        will call function `pntr` in the following order:
C           `pntr((/ 10._dp, 0._dp /))
C           `pntr((/ 11._dp, 0._dp /))
C           `pntr((/ 10._dp, 1._dp /))
C           `pntr((/ 11._dp, 1._dp /))
C           `pntr((/ 10._dp, 2._dp /))
C           `pntr((/ 11._dp, 2._dp /))
C
C     2. Dependent mode:
C        Calling `var_args_r(pntr, reshape((/10._dp, 11._dp, 0._dp/), 
C           (/0._dp, 1._dp, 2._dp/)), (/3, 2/)), (/2, 2/))`
C        with identical number of elements in `args` along `dim=1)        
C        will call function `pntr` in the following order:
C           `pntr((/ 10._dp, 0._dp /))
C           `pntr((/ 11._dp, 1._dp /))
C
C     3. Mixed independent/dependent mode:
C        Calling `var_args_r(pntr, reshape((/100._dp, 101._dp/), 
C           (/10._dp, 0._dp/), (/0._dp, 1._dp/)), (/2, 3/)), 
C           (/2, 1, 2/))`
C        will call function `pntr` in the following order:
C           `pntr((/100._dp, 10._dp, 0._dp/))`
C           `pntr((/101._dp, 10._dp, 0._dp/))`
C           `pntr((/100._dp, 10._dp, 1._dp/))`
C           `pntr((/101._dp, 10._dp, 1._dp/))`
C
C     4. Mixed independent/dependent mode:
C        Switching 1st and 2nd element of `args` in previous example, 
C        calling `var_args_r(pntr, reshape((/10._dp, 0._dp/), 
C           (/100._dp, 101._dp/), (/0._dp, 1._dp/)), (/2, 3/)), 
C           (/1, 2, 2/))`
C        will call function `pntr` in the following order:
C           `pntr((/10._dp, 100._dp, 0._dp/))`
C           `pntr((/10._dp, 101._dp, 1._dp/))`
C
C----------------------------------------------------------------------|
        ! ----- Function arguments ------------------------------------|
      integer, dimension(:), intent(in) ::
     >  n_elems

      real(kind=dp), dimension(:,:), intent(in) ::
     >  args

      real(kind=dp), dimension(:), allocatable :: res

        ! ----- Function interface ------------------------------------|
      abstract interface
        function func(args)
            real(kind=kind(1.d0)) :: func
            real(kind=kind(1.d0)), dimension(:), intent(in) :: args
        end function func
      end interface
      procedure (func) :: pntr_s

        ! ----- Local variables ---------------------------------------|
      integer :: 
     >  i, ind_l, ind_u, n, n_args, n_sub

      integer, dimension(size(n_elems)) ::
     >  n_elems_sub

      real(kind=dp), dimension(size(args,dim=1),size(args,dim=2)) ::
     >  args_sub

C----------------------------------------------------------------------|
C     Check input dimensions
C----------------------------------------------------------------------|
      if (size(args, dim=2).ne.size(n_elems)) then
        write(*,'(a)') 'VAR_ARGS_R: ERROR: Dimension 1 of input ' //
     >      'arrays do not match.'
        stop
      else if (size(args, dim=1).lt.max_a(n_elems)) then
        write(*,'(a)') 'VAR_ARGS_R: ERROR: Dimension 2 of argument ' //
     >      'array smaller than specified in element array.'
        stop
      endif

C----------------------------------------------------------------------|
C     Scalar input received: Call target function        
C----------------------------------------------------------------------|
      if (all(n_elems(:).eq.1)) then
            ! Allocate array
        allocate(res(1))

            ! Call function with 1-dimensional argument
        res(1) = pntr_s((/ (args(1,i), i=1,size(n_elems)) /))

C----------------------------------------------------------------------|
C     Vector input received: Call target function        
C----------------------------------------------------------------------|
      else
            ! Determine indices of dependent elements and size of 
            ! corresponding array
        call get_dependent_indices_and_size(n_elems, ind_l, ind_u, n)

            ! Set up sub-array
        n_sub = n/n_elems(ind_l)
        args_sub = args
        n_elems_sub = n_elems
        n_elems_sub(ind_l:ind_u) = 1

            ! Allocate array
        allocate(res(1:n))

            ! Call function recursively on subset of input parameters
        do i=1,n_elems(ind_u)
            args_sub(1,ind_l:ind_u) = args(i,ind_l:ind_u)

            res(n_sub*(i-1)+1:n_sub*i) =
     >          var_args_r(pntr_s, args_sub, n_elems_sub)
        
        enddo
      endif

      return
      end function var_args_r
C======================================================================|


C======================================================================|
      subroutine get_dependent_indices_and_size(n_elems, ind_l, ind_u,
     >  n)
C----------------------------------------------------------------------|
C     Finds dependent indices and calculates size of corresponding
C     array.
C----------------------------------------------------------------------|
        ! ---- Subroutine arguments -----------------------------------|
      integer, dimension(:), intent(in) :: 
     >  n_elems

      integer, intent(out) :: 
     >  n, ind_l, ind_u

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i

C----------------------------------------------------------------------|
        ! Determine last argument with > 1 elements
      do i=size(n_elems),1,-1
        if (n_elems(i).gt.1) then
            ind_u=i
            exit
        endif
      enddo

        ! Search for argument directly below with identical size
      do i=ind_u,1,-1
        if (n_elems(i).eq.n_elems(ind_u)) then
            ind_l = i
        else
            exit
        endif
      enddo

        ! Size of corresponding array
      n = product(n_elems(1:ind_l))

      return
      end subroutine get_dependent_indices_and_size
C======================================================================|

      end module var_args
C######################################################################|
