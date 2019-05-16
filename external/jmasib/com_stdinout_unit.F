      module com_stdinout_unit 
!
      implicit none
!
      integer,save :: iunit_in  = 5 
      integer,save :: iunit_out = 6 
!
      CONTAINS 
!================================================================
      SUBROUTINE COM_STDINOUT_UNIT_INI ( iin , iout )
!
        integer,intent(in) :: iin 
        integer,intent(in) :: iout
#ifdef MS_WINDOWS
        logical:: exist
#endif
!
        iunit_in  = iin 
        iunit_out = iout 
!
        ! Windows の統合開発環境ではリダイレクトができないのでこうする 
#ifdef MS_WINDOWS
          inquire(file='STDIN.TXT', exist=exist)
          if (exist) then
            open(unit=iunit_in, file='STDIN.TXT', action='READ')
            inquire(file='STDOUT.TXT', exist=exist)
            if (exist) then
              open(unit=iunit_out, file='STDOUT.TXT', action='WRITE')
            endif
          endif
#endif
!
        RETURN 
      END SUBROUTINE COM_STDINOUT_UNIT_INI
      END MODULE com_stdinout_unit 
        

