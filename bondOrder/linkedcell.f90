    MODULE linkedcell
!=======================================================================
! A code following linked cell list algorithm
! The cell number is counted from along x, then y, then z.
! e.g. the first two layers in xy plane has cell numbers as
!     1 2 3   10 11 12
!     4 5 6   13 14 15
!     7 8 9   16 17 18
! You should make sure the positions are in a simulation box with an 
! origin (0, 0, 0). So, if simulation box is shifted, shift the positions
! accordingly.
!                                                                  HBLUO
!=======================================================================
    IMPLICIT NONE
!
!variables used to contruct the cell map. You have to call cellmap first to 
!obtain the cell information stored in the variables below.
!
    INTEGER :: nx,ny,nz,ncell
    INTEGER,ALLOCATABLE :: neigh(:,:)       ! store 27 cell numbers with the 1st itself and 26 neigbors
    INTEGER,ALLOCATABLE :: iSuper(:,:,:)    ! which supercell the cell is in. Current: (0, 0, 0)
    REAL(8) :: dx,dy,dz        ! store the real cell dimensions
!
!variables used to get head and list arrays
!
    INTEGER,ALLOCATABLE :: cellnum(:) ! the cell number for each atom
    INTEGER,ALLOCATABLE :: head(:)    ! the largest position index of each cell
    INTEGER,ALLOCATABLE :: list(:)    ! one linked position in the same cell
!
!
    CONTAINS
!
    SUBROUTINE cellmap(lx,ly,lz,rcut)
    IMPLICIT NONE
    INTEGER :: i,k, ix,iy,iz, iNeighx, iNeighy, iNeighz
    INTEGER :: ishiftx,ishifty,ishiftz,iBox(3)
    REAL(8) :: lx,ly,lz,rcut
    
!
!number of cell along each direction
!
    nx=INT(lx/rcut)
    ny=INT(ly/rcut)
    nz=INT(lz/rcut)
    ncell=nx*ny*nz    ! the total number of cells
!
!cell dimenstions
!
    dx=lx/nx
    dy=ly/ny
    dz=lz/nz
!
!find out the neighoring cell of each cell. I am going to give all 26 explicitly
!
    IF(ALLOCATED(neigh)) DEALLOCATE(neigh)
    IF(ALLOCATED(iSuper)) DEALLOCATE(iSuper)
    ALLOCATE(neigh(ncell,27))
    ALLOCATE(iSuper(ncell,27,3))
    iSuper(:,:,:)=0
    DO i=1, ncell                          ! find neighbor cells for each cell
       neigh(i,1)=i                        ! the first is itself
       iz=(i-1)/nx/ny + 1                  ! ix, iy, iz are the cell indices of the central cell
       iy=MOD(i-1, nx*ny)/nx + 1
       ix=i - (iy-1)*nx - (iz-1)*nx*ny
!write(*,*) ix,iy,iz, i
       k=1                             ! initialize the count of neighbors so it actually start from 2
       DO ishiftx= -1, 1               ! shift the cell indices
          iNeighx=ix+ishiftx
          IF(iNeighx<=0) THEN
             iBox(1)=-1
             iNeighx=iNeighx+nx  ! if the index is outside the box, translate it
          ELSE IF(iNeighx>nx) THEN
             iBox(1)=1
             iNeighx=iNeighx-nx
          ELSE
             iBox(1)=0
          END IF
          DO ishifty= -1, 1
             iNeighy=iy+ishifty
             IF(iNeighy<=0) THEN
                iBox(2)=-1
                iNeighy=iNeighy+ny
             ELSE IF(iNeighy>ny) THEN
                iBox(2)=1
                iNeighy=iNeighy-ny
             ELSE
                iBox(2)=0
             END IF
             DO ishiftz= -1, 1
                iNeighz=iz+ishiftz
                IF(iNeighz<=0) THEN
                   iBox(3)=-1
                   iNeighz=iNeighz+nz
                ELSE IF(iNeighz>nz) THEN
                   iBox(3)=1
                   iNeighz=iNeighz-nz
                ELSE
                   iBox(3)=0
                END IF
                IF(ishiftx==0 .AND. ishifty==0 .AND. ishiftz==0) THEN      ! exclude the central cell
                   CYCLE
                ELSE
                   k=k+1
                   neigh(i,k)=iNeighx + (iNeighy-1)*nx + (iNeighz-1)*nx*ny ! calculate the cell number of the neighbor
                   iSuper(i,k,:)=iBox(:)
                END IF
             END DO
          END DO
       END DO
    END DO

    END SUBROUTINE
!
    SUBROUTINE link(pos, natom)
    IMPLICIT NONE
    INTEGER :: i, icell, natom
    REAL(8) :: pos(natom,3)
!
    IF(ALLOCATED(head)) DEALLOCATE(head)
    IF(ALLOCATED(list)) DEALLOCATE(list)
    IF(ALLOCATED(cellnum)) DEALLOCATE(cellnum)
    ALLOCATE(head(ncell))
    ALLOCATE(list(natom))
    ALLOCATE(cellnum(natom))
!
!initialize the head for each cell
!
    head(:)=0
!
!loop all the atoms to set up the head and list
!
    DO i=1, natom
       icell=1 + INT(pos(i,1)/dx) + INT(pos(i,2)/dy)*nx + INT(pos(i,3)/dz)*nx*ny
       cellnum(i)=icell
       list(i)=head(icell)
       head(icell)=i
    END DO
!
    END SUBROUTINE
!
!
    END MODULE
