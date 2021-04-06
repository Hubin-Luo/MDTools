    PROGRAM main
!==================================================================
! Averaging the stress tensor around a central atom within a certain
! distance and calculate the Von Mises stress. Also, directly averaged
! Von Mises stresses are given for comparison.
! (1) if the volume per atom is not output from Lammps, output it from
!     ovito by vonoroi analysis before doing the calculation
! (2) the order of properties of each atom in the dump file must be:
!     id type x y z von stress[1-6] vol
! (3) the second line of the dump file should be changed to
!     cutoff ntype
! 
! Note: the default lammps unit (metal) for stress is Bar*vol
!                                                        H.B. Luo
!==================================================================
    USE linkedcell
    IMPLICIT NONE
    INTEGER :: i, j, k, iatom, natom, ntype, itype, numLocal
    INTEGER, ALLOCATABLE :: ncount(:), nEach(:)       ! nEach stores numbers of each type
    REAL(8) :: xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz
    REAL(8) :: pairdist, vectmp(3),dist, shift(3)
    REAL(8),ALLOCATABLE :: atominfo(:,:), shiftpos(:,:)
    REAL(8),ALLOCATABLE :: avrg(:,:)
!
! read in the data
!
    READ(*,*)
    READ(*,*) pairdist, ntype   !set the threshold of the distance to include atoms in the second line of dump frame
    READ(*,*)
    READ(*,*) natom
    READ(*,*)
    READ(*,*) xlo, xhi
    READ(*,*) ylo, yhi
    READ(*,*) zlo, zhi
!
! allocate the arrays. The the columns are in the order: id,type,x,y,z,von,pa[1-6],vol,
! von', pa'[1-6], von''. The primed the quantities are averaged results as mentioned above
!
    ALLOCATE(nEach(ntype)); ALLOCATE(avrg(ntype, 7))   ! store the quantities: v_von and  v_pa[1-6]
    ALLOCATE(ncount(ntype))                            ! store the number of atom locally
    ALLOCATE(atominfo(natom,21))                       ! atomic information
    ALLOCATE(shiftpos(natom,3))                        ! store positions shifted regarding that the box may be shifted
!
    atominfo=0.d0; nEach=0
    READ(*,*)
    DO i=1, natom
       READ(*,*) (atominfo(i,j), j=1, 13)
       nEach(NINT(atominfo(i,2)))=nEach(NINT(atominfo(i,2)))+1
    END DO
!
!get the information of linked cells, head and list
!
    shift(:)=(/xlo,ylo,zlo/)
    lx=xhi-xlo; ly=yhi-ylo; lz=zhi-zlo
    DO i=1, natom
       shiftpos(i,:)=atominfo(i,3:5)-shift(:)
    END DO
    CALL cellmap(lx, ly, lz, pairdist)
    CALL link(shiftpos, natom)
!
!do the loops to calculate averaged stresses
!
    DO i=1, natom   ! set one by one as the center
       avrg=0.d0
       ncount=0
!
!loop all neighboring cells including itself
!
   L1: DO j=1, 27
          iatom=head(neigh(cellnum(i),j))
      L2: DO WHILE (iatom/=0)
             vectmp(1)=atominfo(iatom,3)+iSuper(cellnum(i),j,1)*lx  ! shift positions if out of box
             vectmp(2)=atominfo(iatom,4)+iSuper(cellnum(i),j,2)*ly
             vectmp(3)=atominfo(iatom,5)+iSuper(cellnum(i),j,3)*lz
             vectmp(:)=vectmp(:)-atominfo(i,3:5)
             dist=SQRT(DOT_PRODUCT(vectmp,vectmp))
             IF(dist<pairdist) THEN
                ncount(NINT(atominfo(iatom,2)))=ncount(NINT(atominfo(iatom,2)))+1
                avrg(NINT(atominfo(iatom,2)),:)=avrg(NINT(atominfo(iatom,2)),:)+atominfo(iatom,6:12)/atominfo(iatom,13)
             END IF
             IF(list(iatom)==0) THEN
                EXIT L2
             ELSE
                iatom=list(iatom)
             END IF
          END DO L2
       END DO L1
!
       numLocal=0
       DO itype=1,ntype
          numLocal=numLocal+ncount(itype)
       END DO
       DO itype=1, ntype
          IF(ncount(itype)==0) CYCLE    ! in case there are some species missing within the cutoff
          avrg(itype,:)=avrg(itype,:)/ncount(itype)   ! average over number of atoms of each type separately
          atominfo(i,14:20)=atominfo(i,14:20)+avrg(itype,:)*ncount(itype)/numLocal ! average according to LOCAL composition
       END DO
!
! calculate the von Mises stress according to the averaged stress tensor
!
       atominfo(i,21)=SQRT( 0.5*( (atominfo(i,15)-atominfo(i,16))**2+ &
                                & (atominfo(i,16)-atominfo(i,17))**2+ &
                                & (atominfo(i,17)-atominfo(i,15))**2 )+ &
                                & 3*(atominfo(i,18)*atominfo(i,18)+ &
                                &    atominfo(i,19)*atominfo(i,19)+ &
                                &    atominfo(i,20)*atominfo(i,20)) )
    END DO
!
! print the result
!
    WRITE(*,'("ITEM: TIMESTEP")')
    WRITE(*,'("0")')
    WRITE(*,'("ITEM: NUMBER OF ATOMS")')
    WRITE(*,'(I10)') natom
    WRITE(*,'("ITEM: BOX BOUNDS pp pp pp")')
    WRITE(*,'(2F25.15)') xlo, xhi
    WRITE(*,'(2F25.15)') ylo, yhi
    WRITE(*,'(2F25.15)') zlo, zhi
    WRITE(*,'("ITEM: ATOMS id type x y z v_von1 v_pa1 v_pa2 v_pa3 v_pa4 v_pa5 v_pa6 v_von2")')
    DO i=1, natom
       WRITE(*,'(I8,1X,I2,3F11.5,1X,F8.2,6F9.2,F9.2)') NINT(atominfo(i,1)), NINT(atominfo(i,2)), &
             & (atominfo(i,j), j=3,5), atominfo(i,14),(atominfo(i,j),j=15,20), atominfo(i,21)
    END DO
!
    DEALLOCATE(atominfo)
!
    END PROGRAM
