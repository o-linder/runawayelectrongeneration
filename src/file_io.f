C######################################################################|
      module file_io
C----------------------------------------------------------------------|
C     Provides subroutines for file input and output.
C----------------------------------------------------------------------|
      use double

      implicit none

      private

      public ::
     >  store_real_arr

      contains

C======================================================================|
      subroutine store_real_arr(un, args)
C----------------------------------------------------------------------|
C     Stores a variable number of elements located in real array `args`
C     in an opened file with unit number `un`.
C----------------------------------------------------------------------|
        ! ----- Subroutine arguments ----------------------------------|
      integer, intent(in) ::
     >  un

      real(kind=dp), dimension(:), intent(in) ::
     >  args

        ! ----- Local variables ---------------------------------------|
      character(len=12) ::
     >  frmt

C----------------------------------------------------------------------|
C     Create format string
C----------------------------------------------------------------------|
      write(unit=frmt, fmt='(a,i0.3,a)') 
     >  '(', size(args), 'es20.12)'

C----------------------------------------------------------------------|
C     Write to file
C----------------------------------------------------------------------|
      write(unit=un, fmt=frmt) args

      return
      end subroutine store_real_arr
C======================================================================|


      end module file_io
C######################################################################|
