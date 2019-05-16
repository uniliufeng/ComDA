      MODULE COM_DEBUG
!
      INTEGER,SAVE ::  istep  
      INTEGER,SAVE ::  istep_set  = 1000000  
!     INTEGER,SAVE ::  ipoint_ij  = 11
      INTEGER,SAVE ::  ipoint_ij  = 0
      INTEGER,SAVE ::  ipoint_idx
      INTEGER,SAVE ::  ipoint_new
      INTEGER,SAVE ::  ipmask_idx
      INTEGER,SAVE ::  ipmask_new
!
! ・ipoint_ij 
!    モデル空間での位置。海陸関係無し。ここで与える。
!    ij=1 は北西端、ij=imax*jmax は南東端。
!    ij が増える後とに東へ、東端に達すれば一つ南へ。
! ・ipoint_idx
!    L3SiB 等で用いられる。ipoint_ij に対応する I=1,NMAX の I 。
!    L3SiB では main_sib3.F で設定される。
!    ただし現状(00/08/17) JLMAX=1 以外には対応していない。 
! ・ipoint_new
!    NEWSiB で用いられる。ipoint_ij に対応する I=1,NMAX 等の I 。
!    sib0007_mj.F の INI_LOOP で設定される。
!    JLMAX=1 以外にも対応している。
!
! NEWSIB での使用例
!(1)  do i=nstr,nend ;  enddo の loop 中で使う場合
!      if ( i.eq.ipoint_new ) then 
!        write(6,*) x(i) 
!      endif
!     のようにする
!(2)  do loop の外で使う場合
!      if ( ipoint_new.ge.nstr .and. ipoint_new.le.nend ) then 
!        write(6,*) x(ipoint_new) 
!      endif
!     のようにする。
!
      contains 
      subroutine com_debug_ini
!
      namelist /nam_com_debug / ipoint_ij , istep_set 
!
!      read (5,nam_com_debug)
!      write(6,nam_com_debug)
!
      end subroutine com_debug_ini
      END MODULE COM_DEBUG
