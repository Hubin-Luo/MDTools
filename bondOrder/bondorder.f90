        PROGRAM bondorder
!==========================================================================
! Calculate the local structural order parameter and bond order parameter
! for each atom according to the paper:
! PCCP 17(2015) 27127. The bond order parameter can be used to distinguish
! solid-like and liquid-like atoms.
! The input file is the dump file.
!                                                                Hubin Luo
!==========================================================================
        USE linkedcell
        IMPLICIT NONE
!
        COMPLEX(8),EXTERNAL :: sphHarm
        INTEGER :: i,j,k,l,m,n,natom,iatom,nnn,ierr,ntstp
        REAL(8) :: xlo,xhi,ylo,yhi,zlo,zhi,pi,theta,phi,normSum
        REAL(8) :: pairdist,shift(3),vectmp(3),dist,lx,ly,lz
        REAL(8),ALLOCATABLE :: atominfo(:,:), shiftpos(:,:)
        COMPLEX(8),ALLOCATABLE :: qlm(:,:),qOrig(:,:),qtmp(:,:) !qOrig is unnormalized. qtmp is locally averaged over qOrig
        CHARACTER(50) :: ctemp,dumpFile, string
!
        pi=ACOS(-1.d0); ierr=0
!
! set the l for the spherical harmonics and cutoff for the neighbors and the dump file
!
!        WRITE(*,*) "---Input the l and cutoff---"
        READ(*,*) l, pairdist
!        WRITE(*,*) "---Input the name of dump file---"
        READ(*,*) dumpFile
        OPEN(8,file=TRIM(dumpFile),status="old")
!
! find out how many frames there are
!
        ntstp=0
        DO WHILE(.TRUE.)
           READ(8,'(A50)',iostat=ierr) string
           IF(string(7:14)=="TIMESTEP") ntstp=ntstp+1
           IF(ierr/=0) EXIT
        END DO
        REWIND(8)
!
! loop frames
!
        n=0
        DO WHILE(n<ntstp)
!
! read in the dump file
!
        READ(8,'(A50)') ctemp
        WRITE(*,'(A50)') ctemp
        READ(8,'(A50)') ctemp
        WRITE(*,'(A50)') ctemp
        READ(8,'(A50)') ctemp
        WRITE(*,'(A50)') ctemp
        READ(8,*) natom
        WRITE(*,*) natom
        READ(8,'(A50)') ctemp
        WRITE(*,'(A50)') ctemp
        READ(8,*) xlo, xhi
        READ(8,*) ylo, yhi
        READ(8,*) zlo, zhi
        WRITE(*,*) xlo, xhi
        WRITE(*,*) ylo, yhi
        WRITE(*,*) zlo, zhi
        READ(8,*) ctemp
        WRITE(*,'("ITEM: ATOMS id type x y z c_strOrder c_bondOrder")')
!
        ALLOCATE(atominfo(natom,7))     ! atomic information. The 6th and 7th columns are the structural order para.
        ALLOCATE(shiftpos(natom,3))     ! store positions shifted regarding that the box may be shifted
        ALLOCATE(qlm(natom,2*l+1))      ! store the complex vectors calculate using sphHarm
        ALLOCATE(qOrig(natom,2*l+1));ALLOCATE(qtmp(natom,2*l+1))
        atominfo=0.d0
        DO i=1, natom
           READ(8,*) (atominfo(i,j), j=1, 5)
        END DO
!
! get the information of linked cells, head and list
!
        shift(:)=(/xlo,ylo,zlo/)
        lx=xhi-xlo; ly=yhi-ylo; lz=zhi-zlo
        DO i=1, natom
           atominfo(i,3)=atominfo(i,3)-FLOOR((atominfo(i,3)-xlo)/lx)*lx  ! These 3 statements make sure that
           atominfo(i,4)=atominfo(i,4)-FLOOR((atominfo(i,4)-ylo)/ly)*ly  ! the atoms are in the box since
           atominfo(i,5)=atominfo(i,5)-FLOOR((atominfo(i,5)-zlo)/lz)*lz  ! lammps produces outbox atoms sometimes.
           shiftpos(i,:)=atominfo(i,3:5)-shift(:)
        END DO
        CALL cellmap(xhi-xlo, yhi-ylo, zhi-zlo, pairdist)
        CALL link(shiftpos, natom)
!
! loop the atoms to calculate the complex vectors
!
        qlm(:,:)=(0.d0,0.d0);qOrig(:,:)=(0.d0,0.d0)
        DO i=1,natom
           nnn=0
     L1:   DO j=1,27
              iatom=head(neigh(cellnum(i),j))
        L2:   DO WHILE (iatom/=0)
                 vectmp(1)=atominfo(iatom,3)+iSuper(cellnum(i),j,1)*lx  ! shift pos if out of box
                 vectmp(2)=atominfo(iatom,4)+iSuper(cellnum(i),j,2)*ly
                 vectmp(3)=atominfo(iatom,5)+iSuper(cellnum(i),j,3)*lz
                 vectmp(:)=vectmp(:)-atominfo(i,3:5)
                 dist=SQRT(DOT_PRODUCT(vectmp,vectmp))
!write(*,*) i, dist, iatom, (vectmp(k),k=1,3)
                 IF(dist<pairdist .AND. dist>0.00001) THEN
                    nnn=nnn+1
                    theta=ACOS(vectmp(3)/dist)
                    IF(ABS(SIN(theta))<0.00001) THEN          ! avoid sin(theta)=0 in denominator
                       phi=0.d0                          ! phi can be random according to sphHarm
                    ELSE
                       IF(ABS(vectmp(2))<0.00000001) THEN   ! when vectmp(2)~=0, ACOS may get error 
                          phi=0.d0                          ! if numerically it receives a number out of range
                       ELSE 
                          IF(vectmp(2)>0.d0)  phi=ACOS(vectmp(1)/(dist*SIN(theta)))      ! when phi is smaller than 180
                          IF(vectmp(2)<=0.d0) phi=2*pi-ACOS(vectmp(1)/(dist*SIN(theta))) ! when phi is larger than 180
                       END IF
                    END IF
!write(*,*) theta, phi
                    DO m=-l, l
                       qlm(i,m+l+1)=qlm(i,m+l+1)+sphHarm(l,m,theta,phi)
                    END DO
!write(*,*) (qlm(i,k),k=1,2*l+1)
                 END IF
                 IF(list(iatom)==0) THEN
                    EXIT L2
                 ELSE
                    iatom=list(iatom)
                 END IF
              END DO L2
           END DO L1
!write(*,*) i, nnn, (qlm(i,k),k=1,2*l+1)
!stop
           IF(nnn/=0) qlm(i,:)=qlm(i,:)/nnn               ! average over the number of neighbors
           qOrig(i,:)=qlm(i,:)
           normSum=0.d0                        ! normalize the vectors
           DO m=-l,l
              normSum=normSum+ABS(qlm(i,m+l+1))*ABS(qlm(i,m+l+1))
           END DO
           DO m=-l,l
              qlm(i,m+l+1)=qlm(i,m+l+1)/SQRT(normSum)
           END DO
!write(*,*) i,(qlm(i,k),k=1,2*l+1)
        END DO
!stop
!
! calculate the structural order parameters and bond order para.
!
        qtmp(:,:)=(0.d0,0.d0)
        DO i=1,natom
           nnn=0
      L3:  DO j=1,27
              iatom=head(neigh(cellnum(i),j))
         L4:  DO WHILE (iatom/=0)
                 vectmp(1)=atominfo(iatom,3)+iSuper(cellnum(i),j,1)*lx  ! shift pos if out of box
                 vectmp(2)=atominfo(iatom,4)+iSuper(cellnum(i),j,2)*ly
                 vectmp(3)=atominfo(iatom,5)+iSuper(cellnum(i),j,3)*lz
                 vectmp(:)=vectmp(:)-atominfo(i,3:5)
                 dist=SQRT(DOT_PRODUCT(vectmp,vectmp))
                 IF(dist<pairdist .AND. dist>0.000001) THEN
                    DO m=-l,l
                       atominfo(i,6)=atominfo(i,6)+REAL(qlm(i,m+l+1)*CONJG(qlm(iatom,m+l+1)))
                    END DO
                    nnn=nnn+1
                    qtmp(i,:)=qtmp(i,:)+qOrig(iatom,:)
                 END IF
                 IF(list(iatom)==0) THEN
                    EXIT L4
                 ELSE
                    iatom=list(iatom)
                 END IF
              END DO L4
           END DO L3
           qtmp(i,:)=(qtmp(i,:)+qOrig(i,:))/(nnn+1)
           DO m=-l,l
              atominfo(i,7)=atominfo(i,7)+REAL(qtmp(i,m+l+1)*CONJG(qtmp(i,m+l+1)))
           END DO
           atominfo(i,7)=SQRT(4*pi/(2*l+1)*atominfo(i,7))
        END DO

!
! print the atom information
!
        DO i=1,natom
           WRITE(*,'(I9,I4,3F15.8,F8.2,F8.2)') NINT(atominfo(i,1)),NINT(atominfo(i,2)),(atominfo(i,j),j=3,7)
        END DO
!
! Deallocate
!
        DEALLOCATE(atominfo); DEALLOCATE(qlm,qOrig,qtmp); DEALLOCATE(shiftpos)
        DEALLOCATE(neigh, iSuper, cellnum, head, list)
!
        n=n+1
        END DO
!
        END PROGRAM
