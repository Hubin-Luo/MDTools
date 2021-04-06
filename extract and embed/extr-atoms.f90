	PROGRAM main
!========================================================
! extract atoms from one frame and embed them into another
! frame 
! Presently, it works only for orthorhombic box
!                                                H.B. Luo
!=========================================================
	IMPLICIT NONE
	INTEGER :: i,j,k,natom1,natom2,ntype,n1,n2,nfin,nid,nclst,nmtx
	REAL(8) :: radius,vecmod,dist, theta
	REAL(8) :: celldim(3,2),abc(3),vectmp(3),vecout(3)
	REAL(8) :: exCenter(3),emCenter(3), axis(3),trsl(3)
	REAL(8),ALLOCATABLE :: coor1(:,:),coor2(:,:),clst(:,:),mtx(:,:)
	CHARACTER(20) :: fnm1, fnm2
!
	WRITE(*,*) "Dump Filename(extract from); Data Filename(embed into)"
	READ(*,*) fnm1, fnm2 
	WRITE(*,*) "The radius of the sphere; The distance of first coordinate"
	READ(*,*) radius, dist
        READ(*,*) (axis(i),i=1,3), theta
        trsl(:)=0.d0
!
!read in the coordinates of the cluster
!
	OPEN(8,file=TRIM(fnm1),status='old')
	READ(8,*); READ(8,*); READ(8,*)
	READ(8,*) natom1
	ALLOCATE(coor1(natom1,6))
	READ(8,*)
	DO i=1,3
	   READ(8,*) (celldim(i,j),j=1,2)
	   abc(i)=celldim(i,2)-celldim(i,1)
           IF(abc(i)<2*radius) THEN        ! Check if periodicity can be satisfied
              WRITE(*,*) "Too large cluster radius to extract!!"
              STOP
           END IF
           exCenter(i)=(celldim(i,2)+celldim(i,1))/2 
	END DO
        READ(8,*)
        n1=0
	DO i=1, natom1
	   READ(8,*) (coor1(i,j),j=1,5)
	   vectmp(:)=coor1(i,3:5)-exCenter(:)
	   DO k=1,3
	      IF(ABS(vectmp(k))>abc(k)/2) THEN
	         vectmp(k)=vectmp(k)-abc(k)*vectmp(k)/ABS(vectmp(k))
	      END IF
	   END DO
	   vecmod=SQRT(DOT_PRODUCT(vectmp,vectmp))
	   IF(vecmod>radius) THEN
	      coor1(i,:)=0.d0
	   ELSE
              CALL rotardaxis(axis, trsl, theta, vectmp, vecout, 0)
	      coor1(i,3:5)=vecout(:)
              coor1(i,6)=vecmod
              n1=n1+1
	   END IF
	END DO
        nclst=n1
        ALLOCATE(clst(nclst,6))
!
        OPEN(11,file="EXCOOR",status='replace')
        WRITE(11,*) "#Extracted from the dump file"
        WRITE(11,'(1X,I9,A6)') n1,"atoms"
        WRITE(11,'(1X,I2,A11)') NINT(MAXVAL(coor1(:,2))), "atom types"
        DO i=1,3
           IF(i==1) WRITE(11,'(1X,2F20.12,A8)') (celldim(i,j),j=1,2), "xlo xhi"
           IF(i==2) WRITE(11,'(1X,2F20.12,A8)') (celldim(i,j),j=1,2), "ylo yhi"
           IF(i==3) WRITE(11,'(1X,2F20.12,A8)') (celldim(i,j),j=1,2), "zlo zhi"
        END DO
        WRITE(11,*)
        WRITE(11,*) "Atoms"
        WRITE(11,*)
        n1=0
        DO i=1,natom1
           IF(MAXVAL(ABS(coor1(i,:)))<0.1) THEN
              CYCLE
           ELSE
              n1=n1+1
              coor1(i,1)=n1
              clst(n1,:)=coor1(i,:)
              WRITE(11,'(1X,2I9,3F15.6)') NINT(coor1(i,1)),NINT(coor1(i,2)),&
            &                             (coor1(i,j),j=3,5)
           END IF
        END DO
	CLOSE(8); CLOSE(11)
        CALL quicksort(clst,1,nclst,nclst,6,6)
!
!read in the coordinates in which the cluster is embeded
!
        n2=0
	OPEN(9,file=TRIM(fnm2),status='old')
	READ(9,*)
	READ(9,*) natom2
        ALLOCATE(coor2(natom2,6))
	READ(9,*) ntype
	DO i=1,3
           READ(9,*) (celldim(i,j),j=1,2)
           abc(i)=celldim(i,2)-celldim(i,1)
           IF(abc(i)<2*radius) THEN        ! Check if periodicity can be satisfied
              WRITE(*,*) "Matrix is not large enough to contain the cluster!!"
              STOP
           END IF
           emCenter(i)=(celldim(i,2)+celldim(i,1))/2
        END DO
	READ(9,*);READ(9,*);READ(9,*)
	DO i=1,natom2
	   READ(9,*) (coor2(i,j),j=1,5)
	   vectmp=coor2(i,3:5)-emCenter(:)
           DO k=1,3
              IF(ABS(vectmp(k))>abc(k)/2) THEN
                 vectmp(k)=vectmp(k)-abc(k)*vectmp(k)/ABS(vectmp(k))
              END IF
           END DO
           vecmod=SQRT(DOT_PRODUCT(vectmp,vectmp))
           IF(vecmod<radius) THEN
              coor2(i,:)=0.d0
           ELSE
              coor2(i,6)=vecmod
              n2=n2+1
           END IF   
        END DO
        nmtx=n2
        ALLOCATE(mtx(nmtx,6))
        CLOSE(9)
        n2=0
        DO i=1,natom2
           IF(MAXVAL(ABS(coor2(i,:)))>0.1) THEN
              n2=n2+1
              mtx(n2,:)=coor2(i,:)
           END IF
        END DO
!
!remove the atoms in mtx(:,:) that are too close to the cluster atoms
!
        DO i=1,nclst
           clst(i,3:5)=clst(i,3:5)+emCenter(:)
           DO j=3,5
              IF(clst(i,j)<celldim(j-2,1)) clst(i,j)=clst(i,j)+abc(j-2)
              IF(clst(i,j)>celldim(j-2,2)) clst(i,j)=clst(i,j)-abc(j-2)
           END DO
        END DO
        CALL quicksort(mtx,1,nmtx,nmtx,6,6)
        DO i=nclst,1,-1
           IF(clst(i,6)<radius-dist) EXIT
           DO j=1,nmtx
              IF(mtx(j,6)>radius+dist) EXIT
              vectmp=mtx(j,3:5)-clst(i,3:5)
              DO k=1,3
                 IF(ABS(vectmp(k))>abc(k)/2) THEN
                    vectmp(k)=vectmp(k)-abc(k)*vectmp(k)/ABS(vectmp(k))
                 END IF
              END DO
              vecmod=SQRT(DOT_PRODUCT(vectmp,vectmp))
              IF(vecmod<dist) mtx(j,:)=0.d0
           END DO
        END DO
        n2=0
        DO j=1,nmtx
           IF(MAXVAL(ABS(mtx(j,:)))<0.1) THEN
              CYCLE
           ELSE
              n2=n2+1
           END IF
        END DO
        nfin=nclst+n2
!
!write the head lines of the NEWXYZ
!
        OPEN(10,file="NEWXYZ",status='replace')
        WRITE(10,*)
        WRITE(10,*) nfin, "atoms"
        WRITE(10,*) ntype, "atom types"
        DO i=1,3
           IF(i==1) WRITE(10,'(1X,2F20.12,A8)') (celldim(i,j),j=1,2), "xlo xhi"
           IF(i==2) WRITE(10,'(1X,2F20.12,A8)') (celldim(i,j),j=1,2), "ylo yhi"
           IF(i==3) WRITE(10,'(1X,2F20.12,A8)') (celldim(i,j),j=1,2), "zlo zhi"
        END DO
        WRITE(10,*)
        WRITE(10,*) "Atoms"
        WRITE(10,*)
!
!write the coordinates and re-assign IDs to the atoms
!
        nid=0
        DO i=1,nclst
           nid=nid+1
           clst(i,1)=nid
           WRITE(10,'(1X,2I9,3F15.6)') NINT(clst(i,1)),NINT(clst(i,2)), &
         &                             (clst(i,j),j=3,5)
        END DO
        DO i=1,nmtx
           IF(MAXVAL(ABS(mtx(i,:)))<0.1) THEN
              CYCLE
           ELSE
              nid=nid+1
              mtx(i,1)=nid
              WRITE(10,'(1X,2I9,3F15.6)') NINT(mtx(i,1)),NINT(mtx(i,2)), &
            &                             (mtx(i,j),j=3,5)
           END IF
        END DO
        CLOSE(10); DEALLOCATE(coor1,coor2,clst,mtx)
!
        END PROGRAM

!=======================================================================================
        SUBROUTINE rotardaxis(axis, tranv, thetadeg, vecin, vecout, pov)
!---------------------------------------------------------------
! Rotate a vector around a given axis
! tranv is the vector to translate the axis to origin(not unique)
! The last parameter defines rotation of position(0) or vector(1)
!                                                      H. B. Luo
!---------------------------------------------------------------
        IMPLICIT NONE
        INTEGER :: i, j, k
        INTEGER,INTENT(IN) :: pov
        REAL(8),INTENT(IN) :: thetadeg
        REAL(8),INTENT(IN) :: axis(3), tranv(3), vecin(3)
        REAL(8),INTENT(OUT) :: vecout(3)
        REAL(8) :: vectmpin(4), vectmpout(4)
        REAL(8) :: rotxp(4,4), rotyp(4,4), rotzp(4,4), tranmp(4,4)
        REAL(8) :: rotxn(4,4), rotyn(4,4), rotfin(4,4),tranmn(4,4)
        REAL(8) :: sinalph, cosalph, sinbeta, cosbeta, theta, pi
!
        pi=ACOS(-1.d0)
!        axis=(/1.d0,1.d0,0.d0/); theta=pi
!        vecin=(/1.d0, 0.d0, 0.d0/)
        rotxp=0.d0; rotyp=0.d0; rotzp=0.d0; tranmp=0.d0; rotfin=0.d0
!
! set the rotation sin and cos of the axis around x and y. remember if the
! axis is x axis, we have to set sinalph and cosalph manually to avoid singularity
!
        IF(ABS(axis(2))<1.d-6 .AND. ABS(axis(3))<1.d-6) THEN
           cosalph=1.d0
           sinalph=0.d0
        ELSE
           sinalph=axis(2)/(SQRT(axis(2)*axis(2)+axis(3)*axis(3)))
           cosalph=axis(3)/(SQRT(axis(2)*axis(2)+axis(3)*axis(3)))
        END IF
        sinbeta=axis(1)/SQRT(axis(1)*axis(1)+axis(2)*axis(2)+axis(3)*axis(3))
        cosbeta=SQRT(axis(2)*axis(2)+axis(3)*axis(3))/ &
       &        SQRT(axis(1)*axis(1)+axis(2)*axis(2)+axis(3)*axis(3))
        theta=thetadeg/180*pi
!
        rotxp(1,1)=1.d0;
        rotxp(2,2)=cosalph
        rotxp(2,3)=-sinalph
        rotxp(3,2)=-rotxp(2,3)
        rotxp(3,3)=rotxp(2,2)
        rotxp(4,4)=1.d0
        rotxn=rotxp; rotxn(2,3)=-rotxn(2,3); rotxn(3,2)=-rotxn(3,2)
!
        rotyp(1,1)=cosbeta
        rotyp(1,3)=sinbeta
        rotyp(3,1)=-rotyp(1,3)
        rotyp(3,3)=rotyp(1,1)
        rotyp(2,2)=1.d0
        rotyp(4,4)=1.d0
        rotyn=rotyp; rotyn(1,3)=-rotyn(1,3); rotyn(3,1)=-rotyn(3,1)
!
        rotzp(1,1)=COS(theta)
        rotzp(1,2)=-SIN(theta)
        rotzp(2,1)=-rotzp(1,2)
        rotzp(2,2)=rotzp(1,1)
        rotzp(3,3)=1.d0
        rotzp(4,4)=1.d0
!
        DO i=1, 4
           tranmp(i, i)=1.d0
           tranmn(i, i)=1.d0
        END DO
        tranmn(1:3, 4)=-tranv
        tranmp(1:3, 4)=tranv
!
        vectmpin(1:3)=vecin
        IF (pov==0) THEN
            vectmpin(4)=1.d0     ! rotate position
        ELSE IF (pov==1) THEN
            vectmpin(4)=0.d0     ! rotate vector
        ELSE
            WRITE(*,*) "[ROTARDAXIS]: Please define the rotation, position or vector"
            STOP
        END IF
        DO i=1, 4
           rotfin(i, i)=1.d0
        END DO
        rotfin=MATMUL(rotfin, tranmp)
        rotfin=MATMUL(rotfin, rotxn)
        rotfin=MATMUL(rotfin, rotyp)
        rotfin=MATMUL(rotfin, rotzp)
        rotfin=MATMUL(rotfin, rotyn)
        rotfin=MATMUL(rotfin, rotxp)
        rotfin=MATMUL(rotfin, tranmn)
!
        vectmpout=MATMUL(rotfin, vectmpin)
        vecout=vectmpout(1:3)
!
        END SUBROUTINE
!
!===============================================================================
        RECURSIVE SUBROUTINE quicksort(a, first, last, nlow, ncolum, icolum)
        IMPLICIT NONE
        INTEGER :: first, last, nlow, ncolum, icolum
        INTEGER :: i, j
        REAL(8) :: x
        REAL(8) :: a(nlow,ncolum)
        REAL(8) :: t(ncolum)

        x = a( (first+last) / 2, icolum )
        i = first
        j = last
        DO
           DO WHILE (a(i,icolum) < x)
              i=i+1
           END DO
           DO WHILE (x < a(j,icolum))
              j=j-1
           END DO
           IF (i >= j) EXIT
           t(:) = a(i,:);  a(i,:) = a(j,:);  a(j,:) = t(:)
           i=i+1
           j=j-1
       END DO
       IF (first < i-1) CALL quicksort(a, first, i-1, nlow, ncolum, icolum)
       IF (j+1 < last)  CALL quicksort(a, j+1, last, nlow, ncolum, icolum)
       END SUBROUTINE quicksort
