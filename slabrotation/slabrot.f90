	PROGRAM main
!===========================================================================
! rotate the extracted slab from ovito (lammps data file) first around the 
! x-axis to make a smallest span of position.z, then rotate it around y-axis 
! to get a further smallest span of position.z, then rotate it around z-axis 
! to get a further smallest span of position.y. Further, we can get the slab 
! with orthogonal orientation of the box best fitting its volume.
!                                                                   H.B. Luo
!===========================================================================
	IMPLICIT NONE
	INTEGER :: i, j, k, ndeg,natom,ntype
	INTEGER :: iloc(1)                       ! to recieve a location in a vector
	INTEGER,ALLOCATABLE :: itypenum(:)
	REAL(8),ALLOCATABLE :: atominfo(:,:), postmp(:,:), yspan(:), zspan(:)
	REAL(8) :: degstp, rotdeg, xlo, xhi, ylo, yhi, zlo, zhi, shift(3)
	CHARACTER(50) :: fmt
!
        ndeg=360
	degstp=180.d0/ndeg
!
	READ(*,*)
	READ(*,*) natom
	READ(*,*) ntype
	DO i=1,6
	   READ(*,*)
	END DO
!
	ALLOCATE(itypenum(ntype)); ALLOCATE(atominfo(natom,5))
	ALLOCATE(postmp(natom,3))
	ALLOCATE(yspan(ndeg+1)); ALLOCATE(zspan(ndeg+1))
!
	itypenum(:)=0
	DO i=1,natom
	   READ(*,*) (atominfo(i,j),j=1,5)
	   DO k=1, ntype
	      IF(NINT(atominfo(i,2))==k) itypenum(k)=itypenum(k)+1
	   END DO
	END DO
!
! start the rotation
!
	yspan(:)=0.d0; zspan(:)=0.d0; postmp(:,:)=0.d0
   !
   !rotate around x axis
   !
	DO i=0, ndeg
	   DO j=1, natom
	      CALL rotardaxis((/1.d0,0.d0,0.d0/), (/0.d0,0.d0,0.d0/), degstp*i, atominfo(j,3:5), postmp(j,:), 0)
	   END DO
	   zspan(i+1)=MAXVAL(postmp(:,3))-MINVAL(postmp(:,3))
	END DO
     !
     !find out the degree resulting smallest zspan and rotate the positions accordingly
     !
	iloc=MINLOC(zspan)
	DO j=1, natom
	   CALL rotardaxis((/1.d0,0.d0,0.d0/), (/0.d0,0.d0,0.d0/), degstp*(iloc(1)-1), atominfo(j,3:5), postmp(j,:), 0)
	END DO
	atominfo(:,3:5)=postmp(:,:)
   !
   !rotate around y axis
   !
	DO i=0, ndeg
	   DO j=1, natom
	      CALL rotardaxis((/0.d0,1.d0,0.d0/), (/0.d0,0.d0,0.d0/), degstp*i, atominfo(j,3:5), postmp(j,:), 0)
           END DO
	   zspan(i+1)=MAXVAL(postmp(:,3))-MINVAL(postmp(:,3))
	END DO
     !
     !find out the degree resulting smallest zspan and rotate the positions accordingly
     !
	iloc=MINLOC(zspan)
	DO j=1, natom
	   CALL rotardaxis((/0.d0,1.d0,0.d0/), (/0.d0,0.d0,0.d0/), degstp*(iloc(1)-1), atominfo(j,3:5), postmp(j,:), 0)
	END DO
	atominfo(:,3:5)=postmp(:,:)
   !
   !rotate around z axis
   !
	DO i=0, ndeg
	   DO j=1, natom
	      CALL rotardaxis((/0.d0,0.d0,1.d0/), (/0.d0,0.d0,0.d0/), degstp*i, atominfo(j,3:5), postmp(j,:), 0)
	   END DO
	   yspan(i+1)=MAXVAL(postmp(:,2))-MINVAL(postmp(:,2))
	END DO
     !
     !find out the degree resulting smallest yspan and rotate the positions accordingly
     !
	iloc=MINLOC(yspan)
	DO j=1, natom
	   CALL rotardaxis((/0.d0,0.d0,1.d0/), (/0.d0,0.d0,0.d0/), degstp*(iloc(1)-1), atominfo(j,3:5), postmp(j,:), 0)
	END DO
	atominfo(:,3:5)=postmp(:,:)
!
! Now the rotations are done. Find out the vector to shift the positions to be all positive.
! Then, reshape the box to best fit the position.
!
	shift(1)=MINVAL(atominfo(:,3))-0.5       !introduce a little redundance
	shift(2)=MINVAL(atominfo(:,4))-0.5
	shift(3)=MINVAL(atominfo(:,5))-0.5
	DO i=1, natom
	   atominfo(i,3:5)=atominfo(i,3:5)-shift(:)
	END DO
	xlo=0.d0; xhi=MAXVAL(atominfo(:,3))+0.5
	ylo=0.d0; yhi=MAXVAL(atominfo(:,4))+0.5
	zlo=0.d0; zhi=MAXVAL(atominfo(:,5))+0.5
!
! print the new data file
!
	OPEN(9,file='newslab.xyz', status='replace')
	WRITE(fmt,*) ntype
	WRITE(9,'("#generated by SLABROT with number of atoms:",1X,'//TRIM(ADJUSTL(fmt))//'I10)') (itypenum(k),k=1,ntype)
	WRITE(9,*) natom, "atoms"
	WRITE(9,*) ntype, "atom types"
	WRITE(9,'(1X,2F16.8,1X,"xlo",1X,"xhi")') xlo, xhi
	WRITE(9,'(1X,2F16.8,1X,"ylo",1X,"yhi")') ylo, yhi
	WRITE(9,'(1X,2F16.8,1X,"zlo",1X,"zhi")') zlo, zhi
	WRITE(9,*)
	WRITE(9,'(A5)') "Atoms"
	WRITE(9,*)
	DO i=1, natom
	   WRITE(9,'(1X,I10,1X,I3,1X,3F16.7)') NINT(atominfo(i,1)), NINT(atominfo(i,2)), (atominfo(i,j),j=3,5)
	END DO
	CLOSE(9)
	DEALLOCATE(itypenum); DEALLOCATE(atominfo)
        DEALLOCATE(postmp)
        DEALLOCATE(yspan); DEALLOCATE(zspan)
!
	END PROGRAM
!
!===============================================================================
!

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
        rotxp=0.d0; rotyp=0.d0; rotzp=0.d0
        tranmp=0.d0; tranmn=0.d0; rotfin=0.d0
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
