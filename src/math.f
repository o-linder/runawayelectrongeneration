C######################################################################|
      module math

      use double

      implicit none

      interface max_a
          module procedure max_ai, max_ar
      end interface

      interface min_a
          module procedure min_ai, min_ar
      end interface

      public ::
     >  max_a,
     >  min_a

      private ::
     >  max_ai, max_ar,
     >  min_ai, min_ar,
     >
     >  dp

      contains

C======================================================================|
C     Function max_a
C----------------------------------------------------------------------|
C     Determines maximum value of an array
C----------------------------------------------------------------------|
      !----------------------------------------------------------------|
      integer function max_ai(args)
      !----------------------------------------------------------------|
      integer, dimension(:), intent(in) :: args
      integer :: i
      max_ai = args(1)
      do i=2,size(args)
        if (args(i).gt.max_ai) max_ai = args(i)
      enddo
      
      return
      end function max_ai
      !----------------------------------------------------------------|


      !----------------------------------------------------------------|
      real(kind=dp) function max_ar(args)
      !----------------------------------------------------------------|
      real(kind=dp), dimension(:), intent(in) :: args
      integer :: i
      max_ar = args(1)
      do i=2,size(args)
        if (args(i).gt.max_ar) max_ar = args(i)
      enddo
      
      return
      end function max_ar
      !----------------------------------------------------------------|
C======================================================================|


C======================================================================|
C     Function min_a
C----------------------------------------------------------------------|
C     Determines minimum value of an array
C----------------------------------------------------------------------|
      !----------------------------------------------------------------|
      integer function min_ai(args)
      !----------------------------------------------------------------|
      integer, dimension(:), intent(in) :: args
      min_ai = -1*max_a(-1*args)
      return
      end function min_ai
      !----------------------------------------------------------------|


      !----------------------------------------------------------------|
      real(kind=dp) function min_ar(args)
      !----------------------------------------------------------------|
      real(kind=dp), dimension(:), intent(in) :: args
      min_ar = -1*max_a(-1*args)
      return
      end function min_ar
      !----------------------------------------------------------------|
C======================================================================|


      end module math
C######################################################################|
