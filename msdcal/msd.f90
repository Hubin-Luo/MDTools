        PROGRAM main
!=======================================================
! calculate root mean squared displacement of different 
! species within certain timesteps
! NOTE !! No PBC considered here. It must be modified
! if the PBC is not ignorable.
! msdrad is a radial msd for the particle.
!                                              H.B.LUO
!=======================================================
        USE dumpInfo
        IMPLICIT NONE
        INTEGER :: i,j,k,l,ifid,ifty,ifx,ify,ifz,itype
        INTEGER :: interval,n
        REAL(8) :: msd,msdrad,ri,rj
        REAL(8) :: vec(3)
        CHARACTER(50) :: fnm
!
        WRITE(*,*) "Input dumpfile, iType and frame interval(determine tau)"
        READ(*,*) fnm, itype, interval
        CALL readDump(fnm)
!
        DO i=1,nframe
           msd=0.d0; msdrad=0.d0; n=0
           IF(i+interval>nframe) EXIT
           DO j=1,SIZE(frame(i)%field(:))
              IF(TRIM(frame(i)%field(j))=="id") ifid=j
              IF(TRIM(frame(i)%field(j))=="type") ifty=j
              IF(TRIM(frame(i)%field(j))=="x") ifx=j
              IF(TRIM(frame(i)%field(j))=="y") ify=j
              IF(TRIM(frame(i)%field(j))=="z") ifz=j
           END DO
           DO k=1,frame(i)%natom
              IF(NINT(frame(i)%atomInfo(k,ifty))==itype) THEN
                 n=n+1
!                 DO l=1,frame(i+interval)%natom
!                    IF(NINT(frame(i+interval)%atomInfo(l,ifid))==&
!                      &NINT(frame(i)%atomInfo(k,ifid))) THEN
!                       vec(1)=frame(i+interval)%atomInfo(l,ifx)-&
!                             &frame(i)%atomInfo(k,ifx)
!                       vec(2)=frame(i+interval)%atomInfo(l,ify)-&
!                             &frame(i)%atomInfo(k,ify)
!                       vec(3)=frame(i+interval)%atomInfo(l,ifz)-&
!                             &frame(i)%atomInfo(k,ifz)
!                       msd=msd+DOT_PRODUCT(vec,vec)
                       vec(1)=frame(i)%atomInfo(k,ifx)-&
                             &(frame(i)%box(1,2)+frame(i)%box(1,1))/2
                       vec(2)=frame(i)%atomInfo(k,ify)-&
                             &(frame(i)%box(2,2)+frame(i)%box(2,1))/2
                       vec(3)=frame(i)%atomInfo(k,ifz)-&
                             &(frame(i)%box(3,2)+frame(i)%box(3,1))/2
                       ri=SQRT(DOT_PRODUCT(vec,vec))
                       vec(1)=frame(i+interval)%atomInfo(k,ifx)-&
                             &(frame(i+interval)%box(1,2)+frame(i+interval)%box(1,1))/2
                       vec(2)=frame(i+interval)%atomInfo(k,ify)-&
                             &(frame(i+interval)%box(2,2)+frame(i+interval)%box(2,1))/2
                       vec(3)=frame(i+interval)%atomInfo(k,ifz)-&
                             &(frame(i+interval)%box(3,2)+frame(i+interval)%box(3,1))/2
                       rj=SQRT(DOT_PRODUCT(vec,vec))
                       msdrad=msdrad+(rj-ri)*(rj-ri)
!                       EXIT
!                    END IF
!                 END DO
              END IF
           END DO
!           msd=msd/n
           msdrad=msdrad/n
!           WRITE(*,'(I6,2F20.4)') i+interval-1, msd, msdrad
           WRITE(*,'(I6,F20.4)') i+interval-1, msdrad
        END DO
!
        END PROGRAM
