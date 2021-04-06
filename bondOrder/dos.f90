        PROGRAM dos
!-----------------------------------------------------------
! calculate the dos using Gaussion function
! It reads the file with surrogate energies output from
! expan.exe calcualtion                           Hu-Bin Luo
!-----------------------------------------------------------
        IMPLICIT NONE
        INTEGER :: i,ie,nein,neout,itmp, ntstp, ierr, n
        INTEGER :: iq,nqout,nLiq, nCry                 ! numbers of atoms for Liquid and Crystal
        REAL(8) :: pi,tmp
        REAL(8),ALLOCATABLE :: eng(:),qLiq(:),qCry(:)
        REAL(8) :: sigma,emin,emax,de,eout,ediff,qmin,qmax
        REAL(8) :: tdos,qdos1,qdos2,qout,dq,qdiff,sigmaq,stroderthres
        CHARACTER(50) :: fnm, string
!
        pi=ACOS(-1.d0)
!        WRITE(*,*) "input the filename where energies stored..."
        READ(*,*) fnm
!        WRITE(*,*) "numbers of output energies..."
        READ(*,*) neout,nqout
!        WRITE(*,*) "input the smearing in Gaussion function..."
        READ(*,*) sigma,sigmaq,stroderthres
!        WRITE(*,*) "input emin, emax..."
        READ(*,*) emin,emax                    ! for S distribution
        READ(*,*) qmin,qmax                    ! for Q distribution
!
        OPEN(8,file=TRIM(fnm),status="old")
!
        ntstp=0
        DO WHILE(.TRUE.)
           READ(8,'(A50)',iostat=ierr) string
           IF(string(7:14)=="TIMESTEP") ntstp=ntstp+1
           IF(ierr/=0) EXIT
        END DO
        REWIND(8)
!
        n=0
        DO WHILE(n<ntstp)             !loop frames
        nLiq=0; nCry=0
        WRITE(string,*) n
        OPEN(9,file="dos-S-"//TRIM(ADJUSTL(string))//".txt",&
            &status="replace")
        OPEN(10,file="dos-Q-"//TRIM(ADJUSTL(string))//".txt",&
            &status="replace")
!
        DO i=1,9                      !skip the first 9 lines in dump file
           IF(i==4) THEN
              READ(8,*) nein
              ALLOCATE(eng(nein));ALLOCATE(qLiq(nein));ALLOCATE(qCry(nein))
           ELSE 
              READ(8,*)
           END IF
        END DO
        DO i=1,nein
           READ(8,*) tmp,tmp,tmp,tmp,tmp,eng(i),tmp !structural order para. in the 6th column
           IF(eng(i)<strOderThres) THEN
             nLiq=nLiq+1
             qLiq(nLiq)=tmp
           ELSE
             nCry=nCry+1
             qCry(nCry)=tmp
           END IF
        END DO
!
        de=(emax-emin)/(neout-1)
        DO ie=1,neout
           eout=emin+(ie-1)*de
           tdos=0.d0
           DO i=1,nein
              ediff=eng(i)-eout
              tdos=tdos+exp(-ediff*ediff/(2*sigma*sigma)) &
             &        /SQRT(2*pi*sigma*sigma)
           END DO
           tdos=tdos/nein
           WRITE(9,*) eout, tdos
        END DO
!
        dq=(qmax-qmin)/(nqout-1)
        DO iq=1,nqout
           qout=qmin+(iq-1)*dq
           qdos1=0.d0; qdos2=0.d0
           DO i=1,nLiq
              qdiff=qLiq(i)-qout
              qdos1=qdos1+exp(-qdiff*qdiff/(2*sigmaq*sigmaq)) &
             &        /SQRT(2*pi*sigmaq*sigmaq)
           END DO
           qdos1=qdos1/nein      ! =qdos1/nLiq*nLiq/nein (normalized then scaled by composition)
           DO i=1,nCry
              qdiff=qCry(i)-qout
              qdos2=qdos2+exp(-qdiff*qdiff/(2*sigmaq*sigmaq)) &
             &        /SQRT(2*pi*sigmaq*sigmaq)
           END DO
           qdos2=qdos2/nein
           WRITE(10,*) qout, qdos1, qdos2
        END DO
        CLOSE(9);CLOSE(10)
!
        n=n+1
        DEALLOCATE(eng,qLiq,qCry)
        END DO
!
        CLOSE(8)
!
        END PROGRAM
