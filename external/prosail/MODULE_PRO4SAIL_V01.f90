MODULE ANGLE
	REAL(KIND=8),SAVE :: pi,rd
END MODULE

!*********************

MODULE dataparm
	INTEGER(KIND=4),SAVE :: ns
	INTEGER,PARAMETER :: nw=2101
END MODULE

!*********************

MODULE lidf_parm
	INTEGER(KIND=4),SAVE :: na
	REAL(KIND=8),ALLOCATABLE,SAVE :: thetal(:)
END MODULE

!*********************

MODULE spectral
	REAL(KIND=8),ALLOCATABLE,SAVE :: rsoil(:,:),rsoil0(:)
	REAL(KIND=8),ALLOCATABLE,SAVE :: PARdiro(:),PARdifo(:)
END MODULE

!*********************

MODULE staticvar
	REAL(KIND=8),ALLOCATABLE,SAVE :: lidf(:)
	REAL(KIND=8),SAVE :: cts,cto,ctscto
	REAL(KIND=8),SAVE :: ddb,ddf,dob,dof,dso
	REAL(KIND=8),SAVE :: ko,ks,sdb,sdf
	REAL(KIND=8),SAVE :: sob,sof,sumint
	REAL(KIND=8),ALLOCATABLE,SAVE :: sb(:),sf(:),vb(:),vf(:),w(:)
	REAL(KIND=8),ALLOCATABLE,SAVE :: m(:),m2(:),att(:),sigb(:),rinf(:)
END MODULE

!*********************

MODULE output_PROSPECT
	REAL(KIND=8),ALLOCATABLE,SAVE :: LRT(:,:),rho(:),tau(:)
END	MODULE

!*********************

MODULE flag_util
	LOGICAL flag(9),init_completed,init_completed0
	LOGICAL delta_geom,delta_lai,delta_hot,delta_leaf,delta_skyl,delta_soil,delta_lidf
	REAL(KIND=8),SAVE :: N_old,Cab_old,Car_old,Cbrown_old,Cw_old,Cm_old
	REAL(KIND=8),SAVE :: lai_old,angl_old,q_old,psoil_old,skyl_old
	REAL(KIND=8),SAVE :: tts_old,tto_old,psi_old
END	MODULE

!*********************

MODULE rfltrn
	REAL(KIND=8),ALLOCATABLE,SAVE :: rsos(:),rsod(:),rddt(:),rsdt(:),rdot(:),rsodt(:),rsost(:),rsot(:)
	REAL(KIND=8),SAVE :: tsstoo
END MODULE

!*********************

MODULE outSAIL
	REAL(KIND=8),SAVE :: tss,too
	REAL(KIND=8),ALLOCATABLE,SAVE :: tso(:),tsd(:),tdd(:),tdo(:),rsd(:),rdd(:),rso(:),rdo(:)
END MODULE

!*********************
