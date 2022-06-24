module CDS_Mod

use Matrix_Base

implicit none

type, extends(t_matrix_base)  ::  t_cds  !A type designed to store a matrix in Compressed Diagonal storage fromat using the Intel MKL diagonal storage format
  integer, private                            ::  ndiag !An integer containing the number of diagonals with non-zero elements
  integer, dimension(:), allocatable, private ::  distance  !An array which contains the locations of the diagonals with non-zero elements
  real(kind=dp), dimension(:,:), allocatable  ::  values  !An array which contains the values found in each diagonal. The first dimension increases with increasing diagonal number, the second with increasing row number
contains
  procedure ::  construct => cds_construct
  procedure ::  destroy => cds_destroy
  procedure ::  get => cds_get
  procedure ::  find_diag_ref => cds_find_diag_ref
  procedure ::  set_values => cds_set_values
  procedure ::  return_size => cds_return_size
  procedure ::  set_ndiag => cds_set_ndiag
  procedure ::  return_ndiag => cds_return_ndiag
  procedure ::  set_distance => cds_set_distance
  procedure ::  return_distance => cds_return_distance
  procedure ::  set => cds_set
  procedure ::  remove_zero => cds_remove_zero
  procedure ::  check_symmetric => cds_check_symmetric
  procedure ::  max_entries_row => cds_max_entries_row
  procedure ::  operate => operate_cds

end type t_cds

contains

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

subroutine cds_construct(this, n_row, n_column) !A subroutine which takes nthe size of the matrix as an argument to create an instance of the type t_cds
  class(t_cds), intent(inout) :: this !The instance of t_cds to be created
  integer, intent(in)              :: n_row, n_column !the size of the matrix

  call this%destroy()

  if (n_column /= n_row)then
    Error Stop "In cds_construct, the number of rows and number of columns were not equal. Terminating."
  end if

  this%n_row=n_row
  this%ndiag=0

end subroutine

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

subroutine cds_destroy(this) !A subroutine which takes an instance of the t_cds type and deallocates all arrays and returns the equivalent of a matrix of size 0.
  class(t_cds), intent(inout) :: this !The instance of t_cds to be destroyed

  if (allocated(this%distance)) deallocate(this%distance)
  if (allocated(this%values)) deallocate(this%values)

end subroutine

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

integer function cds_return_size(this)
  class(t_cds), intent(in) ::  this !The matrix whose size is to be returned

  cds_return_size=this%n_row

end function

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

subroutine cds_set_ndiag(this, ndiag) !A subroutine which sets the number of non-zero diagonals. Sets all the distances and values to zero
  class(t_cds), intent(inout)  :: this !The instance of t_cds
  integer, intent(in)               :: ndiag !The value the size the matrix will be set to

  if (allocated(this%values)) deallocate(this%values)
  if (allocated(this%distance)) deallocate(this%distance)

  this%ndiag=ndiag

  allocate(this%values(this%ndiag,this%n_row), this%distance(this%ndiag))

  this%distance=0
  this%values=0.0_dp

end subroutine

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

integer function cds_return_ndiag(this)
  class(t_cds), intent(in) ::  this !The matrix for which the number of non-zero diagonals is to be returned

  cds_return_ndiag=this%ndiag

end function

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

subroutine cds_set_distance(this, distance) !A subroutine which sets the distance the non-zero diagonals are from the central diagonal.
  class(t_cds), intent(inout)  :: this !The instance of t_cds
  integer, dimension(:), intent(in) :: distance !The array which distances will be set to

  !Check that the number of entries in the supplied distance is equal to what is expected.
  if (size(distance) .ne. this%ndiag) then
    !If not, display a warning then reset the matrix with these new diagonals
    print*, "Warning: In cds class you have asked to set the distance parameter for the non-zero diagonals to a different "&
    //"number of values than there are diagonals. Deleting the current contents of matrix, setting the number of diagonals "&
    //"to match and setting all values to zero."
    call this%set_ndiag(size(distance))
  end if

  this%distance=distance

end subroutine

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

function cds_return_distance(this) !This function returns the distance variable
  class(t_cds), intent(in)       ::  this !The matrix for which the number of non-zero diagonals is to be returned
  integer, dimension(:), allocatable  ::  cds_return_distance

  allocate(cds_return_distance(this%ndiag))

  cds_return_distance=this%distance

end function

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

real(kind=dp) function cds_get(this, row, column) !A function which returns the value of a specific location in the matrix
  class(t_cds), intent(in) ::  this !The matric whose values are to be returned
  integer, intent(in)           ::  row  !The row of the value to be returned
  integer, intent(in)           ::  column  !The column of the value to be returned
  integer                       ::  diagonal !The diagonal of the value to be returned
  integer                       ::  location  !The location of the value to be returned within the diagonal
  integer                       ::  diagref  !The first array index of the value to be returned in the values array of the t_cds type
  integer                       ::  ii  !A generic counting variable
  logical                       ::  padded  !Will be set to false if the value is not padded

!print*, "get", this%ndiag

  !Check the called for value is within the matrix. If it's not then end the program
  if (row.gt.this%n_row .or. column.gt.this%n_row .or. row.lt.1 .or. column.lt.1) then
    print*, "The value asked to be retrieved by cds_get is outside the size of the matrix. Terminating program."
    STOP
  end if

  !Check the matrix is non-zero. If it is, return zero
  if (this%ndiag==0) then
    cds_get=0.0_dp
    return
  end if

  !Find the diagonal and the location in that diagonal of the value to be returned
  diagonal=column-row
  location=column+row-1

  call this%find_diag_ref(diagonal, diagref, padded)

  !If the value is not in a stored diagonal the return the value 0
  if (padded) then
    cds_get=0.0_dp
  !Otherwise return the value stored in the values array of the containing t_cds type
  else
    cds_get=this%values(diagref,row)
  end if

end function

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

subroutine cds_set_values(this, nvaluesin, rows, columns, valuesin) !A subroutine which sets specific entries in the matrix to a certain values
  class(t_cds), intent(inout)                      ::  this !The matrix to be modified
  integer, intent(in)                                   ::  nvaluesin  !The number of values to be inserted
  integer, dimension(:), intent(in)                     ::  rows !The row of the entry in the matrix to be modified
  integer, dimension(:), intent(in)                     ::  columns  !The column of the entry to be modified
  real(kind=dp), dimension(:), intent(in)               ::  valuesin !The value which is to be inserted
  integer, dimension(:), allocatable                    ::  diagonals !The diagonal of the new values
  integer, dimension(:), allocatable                    ::  diagref !The references of the diagonal in the initial matrix
  logical, dimension(:), allocatable                    ::  padded !The array telling us whether the current diagonals are padded or not
  integer, dimension(:), allocatable                    ::  newdiags  !The matrix which contains the values of the diagonals to be created
  integer, dimension(:), allocatable                    ::  newdiagstemp  !Temporary version of newdiags
  integer, dimension(:), allocatable                    ::  distancetemp  !Temporary version of distance array
  real(kind=dp), dimension(:,:), allocatable            ::  valuestemp  !Temporary version of the values array from the cds
  integer                                               ::  nnewdiags !The number of new diagonals added
  integer                                               ::  insertref !The entry in the distance matrix a new diagonal will take
  logical                                               ::  newdiag !A variables used to track if a diagonal has already been added
  integer                                               ::  ii,jj !Generic counting variables

  allocate(padded(nvaluesin), diagref(nvaluesin), diagonals(nvaluesin))

  !Find the diagonal of the matrix to be modified
  diagonals=columns-rows

  !Set the number of new diagonals to zero to begin with
  nnewdiags=0

  !Loop over the values to be added
  do ii=1, nvaluesin
    !Find out whether the diagnoals exist and, if so, what their references are
    call this%find_diag_ref(diagonals(ii), diagref(ii), padded(ii))
    if (padded(ii)) then
      !If it is padded then see if the newdiags variable is allocated
      if (allocated(newdiags) .eqv. .FALSE.) then
        !If it isn't then allocate it and store the value
        nnewdiags=1
        allocate(newdiags(nnewdiags))
        newdiags(1)=diagonals(ii)
      else
        !Check to see the diagonal hasn't already been marked to be added
        newdiag=.true.
        do jj=1, nnewdiags
          if (diagonals(ii)==newdiags(jj)) then
            newdiag=.false.
            exit
          end if
        end do
        !If it hasn't already been marked to be added then add it to the list of diagonals to be added, if the value to be added sin't zero
        if (newdiag .and. valuesin(ii).ne.0.0_dp) then
          !Otherwise, temporarily store the newdiags
          nnewdiags=nnewdiags+1
          allocate(newdiagstemp(nnewdiags))
          newdiagstemp(1:nnewdiags-1)=newdiags
          !Now transfer it back to the newdiags once it has a larger size
          deallocate(newdiags)
          allocate(newdiags(nnewdiags))
          newdiags=newdiagstemp
          !Now add the new diagonal and deallocate the temporary array
          newdiags(nnewdiags)=diagonals(ii)
          deallocate(newdiagstemp)
        end if
      end if
    end if
  end do

  !If the values array of the cds is allocated then store the value array of the cds into the temporary array
  if (allocated(this%values)) then
    allocate(valuestemp(this%ndiag, this%n_row))
    valuestemp=this%values
    !Reallocate the value array and put the values back in
    deallocate(this%values)
  end if
print*, "REALLOCATE", this%ndiag+nnewdiags, this%n_row
  !Change values so it is the correct value in previously define entries and zero in new diagonals
  allocate(this%values(this%ndiag+nnewdiags,this%n_row))
  if(allocated(valuestemp))this%values(1:this%ndiag,:)=valuestemp
  this%values(this%ndiag+1:this%ndiag+nnewdiags,:)=0.0_dp

  !Now do the same for distance
  if (allocated(this%distance)) then
    allocate(distancetemp(this%ndiag))
    distancetemp=this%distance
    !Reallocate the value array and put the values back in
    deallocate(this%distance)
  end if

  !Change distance so it is the correct value in previously defined entries and zero for new diagonals which are at the end
  allocate(this%distance(this%ndiag+nnewdiags))
  if(allocated(distancetemp))this%distance(1:this%ndiag)=distancetemp
  this%distance(this%ndiag+1:this%ndiag+nnewdiags)=0

  if (allocated(distancetemp)) deallocate(distancetemp)
  if (allocated(valuestemp)) deallocate(valuestemp)

  !Now insert the new diagonals
  do ii=1, nnewdiags
    !First, find where the diagonal needs to be
    do jj=1, this%ndiag+ii
      if (this%distance(jj).gt.newdiags(ii) .or. jj==this%ndiag+ii) then
        !Now we have found where it should be inserted, insert the diagonal in distance and move the values array appropriately. Set the new diagonal to a value of zero for now.
        if (jj.ne.this%ndiag+ii) then
          this%distance(jj+1:this%ndiag+ii)=this%distance(jj:this%ndiag+ii-1)
          !Transferring the values array to the right is necessary via the array valuestemp to avoid stack overflow in the case of a large number of matrix rows
          allocate(valuestemp(this%ndiag+ii-jj, this%n_row))
          valuestemp=this%values(jj:this%ndiag+ii-1,:)
          this%values(jj+1:this%ndiag+ii,:)=valuestemp
        end if
        this%distance(jj)=newdiags(ii)
        this%values(jj,:)=0.0_dp
        exit
      end if
    end do
  end do

  !Now set the ndiag value of the cds to be correct
  this%ndiag=this%ndiag+nnewdiags

  !Now, find the relevant diagonal for each new value and insert it unless the diagonal is not inserted
  do ii=1,nvaluesin
    call this%find_diag_ref(diagonals(ii), diagref(ii), padded(ii))
    if (padded(ii).eqv. .false.) this%values(diagref(ii),rows(ii))=valuesin(ii)
  end do

  !Remove any zero diagonals which have been created
  call this%remove_zero()

end subroutine cds_set_values

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

subroutine cds_set(this, row, column, value)
  class(t_cds), intent(inout)            ::  this !The matrix to be modified
  integer, intent(in)                         ::  row !The row of the entry to be modified
  integer, intent(in)                         ::  column !The column of the entry to be modified
  real(kind=dp), intent(in)                   ::  value !The value it is to be modifed to
  integer, dimension(1)                       ::  rowarray, columnarray
  real(kind=dp), dimension(1)                 ::  valuearray

  rowarray=row
  columnarray=column
  valuearray=value

  call this%set_values(1, rowarray, columnarray, valuearray)

end subroutine


!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

subroutine cds_find_diag_ref(this, diagonal, diagref, padded)!A subroutine which returns whether the diagonal is padded and, if not, what the reference of the diagonal is.

  class(t_cds), intent(in) ::  this !The matrix whose values are to be returned
  integer                       ::  diagonal !The location of the diagonal whose reference is to be found
  integer, intent(out)          ::  diagref  !The first array index of the value to be returned in the values array of the t_cds type
  integer                       ::  ii  !A generic counting variable
  logical, intent(out)          ::  padded  !Will be set to false if the value is not padded

  !Initially assume the variable is in a padded cell
  padded=.true.
  !Find out which diagonal stored in t_cds type the value to bre returned is and turn padded to false. If
  do ii=1, this%ndiag
    if (diagonal==this%distance(ii)) then
      padded=.false.
      diagref=ii
      exit
    else if (diagonal.lt.this%distance(ii)) then
      exit
    end if
  end do

  !If the diagnoal is padded then return the dummy value -666 as the diagonal value
  if(padded) diagref=-666

end subroutine

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

subroutine cds_remove_zero(this) !This subroutine removes any diagonals which contain only zeros
  class(t_cds), intent(inout)            ::  this !The matrix to be modified
  integer                                     ::  ndiagtemp !Temporary value of nubmer of non-zero diagonals
  logical                                     ::  allzero !False if a diagonal has a non-zero entry
  integer                                     ::  lower, upper  !Lower and upper locations of
  integer, dimension(:), allocatable          ::  distancetemp  !Temporary storage for distance variable
  real(kind=dp), dimension(:,:), allocatable  ::  valuestemp  !Temporary storage for values variable
  integer, dimension(:), allocatable          ::  diagdelref !The references of diagonals to be deleted
  integer                                     ::  ndel  !The number of diagonals to delete
  integer                                     ::  ii,jj !Generic counting variables

  allocate(diagdelref(this%ndiag))

  !Check if there are zero diagonals. If so, return.
  if (this%ndiag==0) return

  !Initially set ndiagtemp to ndiag and reduce it as zero diagonals are found
  ndiagtemp=this%ndiag
  diagdelref=0
  ndel=0

  !Set allzero to true and then loop through each diagonal, setting it to flase if a non-zero entry is found
  do ii=1, this%ndiag
    allzero=.true.
    if (this%distance(ii).ge.0) then
      lower=1
      upper=this%n_row-this%distance(ii)
    else
      lower=1-this%distance(ii)
      upper=this%n_row
    end if
      do jj=lower, upper
        if (this%values(ii,jj).ne.0.0_dp) then
          allzero=.false.
          exit
        end if
      end do
    !If the diagonal only contains zeros then make a note of it and reduce ndiagtemp
    if (allzero) then
      ndel=ndel+1
      diagdelref(ndel)=ii
      ndiagtemp=ndiagtemp-1
    end if
  end do

  !If a diagonal is to be deleted then shift the distance and values arrays into its position instead, if necessary
  do ii=1, ndel-1
    this%distance(diagdelref(ii)-ii+1:diagdelref(ii+1)-ii-1)=this%distance(diagdelref(ii)+1:diagdelref(ii+1)-1)
    this%values(diagdelref(ii)-ii+1:diagdelref(ii+1)-ii-1,:)=this%values(diagdelref(ii)+1:diagdelref(ii+1)-1,:)
  end do
  !Finally, do it for the final diagonal to be deleted, if it is not the final diagonal of the input matrix
  if (ndel.gt.0) then
    if(diagdelref(ndel).ne.this%ndiag) then
      this%distance(diagdelref(ndel)-ndel+1:this%ndiag-ndel)=this%distance(diagdelref(ndel)+1:this%ndiag)
      this%values(diagdelref(ndel)-ndel+1:this%ndiag-ndel,:)=this%values(diagdelref(ndel)+1:this%ndiag,:)
    end if
  end if

!  !Finally, if necessary, place the values and distance variables into temporary arrays before real(kind=dp)locating the main array and insering variables
  if (this%ndiag.ne.ndiagtemp)then
    this%ndiag=ndiagtemp
    allocate(distancetemp(ndiagtemp), valuestemp(ndiagtemp,this%n_row))
    distancetemp=this%distance(1:ndiagtemp)
    valuestemp=this%values(1:ndiagtemp,:)
    deallocate(this%distance, this%values)
    allocate(this%distance(this%ndiag), this%values(this%ndiag,this%n_row))
    this%distance=distancetemp
    this%values=valuestemp
  end if

end subroutine

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

subroutine cds_check_symmetric(this)

!This subroutine checks to see if the matrix is symmetric and tells the user the highest absolute and relative differences and their locations between the two halves

  class(t_cds), intent(inout)    ::  this !The matrix being checked
  real(kind=dp)                       ::  highest_diff !The highest absolute difference between M(ii,jj) and M(jj,ii)
  real(kind=dp)                       ::  highest_rel_diff !The highest relative difference between M(ii,jj) and M(jj,ii)
  integer                             ::  highest_diff_row, highest_diff_column, highest_rel_diff_row, highest_rel_diff_column !The locations of the highest relative and absolute differences
  logical, dimension(:), allocatable  ::  diagonal_checked !A variabel which checks if a particular diagonal has been checked
  integer                             ::  opposite_diag_ref !The location in the distance array of the diagonal symmetric to this one.
  logical                             ::  oppsoite_diag_padded !Is true if the opposite diagonal is padded
  integer, dimension(:), allocatable  ::  matrix_distances !An array storing the distance the different diagonals are away from the main diagonal
  integer                             ::  ii, jj !General counting variables

  if (this%return_ndiag()==0)then
    print*, "The matrix is entirely zero. This means it is trivaially symmetric."
  end if

  !Set up variables
  highest_diff=0.0_dp
  highest_rel_diff=0.0_dp
  highest_diff_row=0
  highest_diff_column=0
  highest_rel_diff_row=0
  highest_rel_diff_column=0

  allocate(diagonal_checked(this%return_ndiag()))

  diagonal_checked=.false.

  do ii=1, this%return_ndiag()
    if (diagonal_checked(ii) .or. this%distance(ii)==0) cycle
    call this%find_diag_ref(-this%distance(ii), opposite_diag_ref, oppsoite_diag_padded)
    do jj=max(1, -this%distance(ii)+1), min(this%return_size(), this%n_row-this%distance(ii))
      if (oppsoite_diag_padded)then
        highest_rel_diff=1.0_dp
        highest_rel_diff_row=jj
        highest_rel_diff_column=this%distance(ii)+jj
        if (abs(this%values(ii, jj))>highest_diff)then
          highest_diff=abs(this%values(ii, jj))
          highest_diff_row=jj
          highest_diff_column=this%distance(ii)+jj
        end if
      else
        if (abs(this%values(ii,jj)-this%values(opposite_diag_ref, this%distance(ii)+jj))/max(abs(this%values(ii,jj)), abs(this%values(opposite_diag_ref, this%distance(ii)+jj)))>highest_rel_diff)then
        highest_rel_diff=abs(this%values(ii,jj)-this%values(opposite_diag_ref, ii+jj))/max(abs(this%values(ii,jj)), abs(this%values(opposite_diag_ref, this%distance(ii)+jj)))
        highest_rel_diff_row=jj
        highest_rel_diff_column=this%distance(ii)+jj
        end if
        if (abs(this%values(ii,jj)-this%values(opposite_diag_ref, this%distance(ii)+jj))>highest_diff)then
          highest_diff=abs(this%values(ii,jj)-this%values(opposite_diag_ref, ii+jj))
          highest_diff_row=jj
          highest_diff_column=this%distance(ii)+jj
        end if
        diagonal_checked(opposite_diag_ref)=.true.
        diagonal_checked(ii)=.true.
      end if
    end do
  end do

  if (highest_diff==0.0_dp)then
    print*, "The matrix is exactly symmetric."
  else
    print*, "The highest absolute differences is between matrix elements (", highest_diff_row, ",", highest_diff_column, ") (", this%get(highest_diff_row, highest_diff_column) ,") and (", highest_diff_column, ",", highest_diff_row, ") (", this%get(highest_diff_column, highest_diff_row) ,")."
    !print*, "The highest relative differences is between matrix elements (", highest_rel_diff_row, ",", highest_rel_diff_column, ") (", this%get(highest_rel_diff_row, highest_rel_diff_column) ,") and (", highest_rel_diff_column, ",", highest_rel_diff_row, ") (", this%get(highest_rel_diff_column, highest_rel_diff_row) ,")."
  end if


end subroutine cds_check_symmetric

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

integer function cds_max_entries_row(this)result(n_rowentriesmax)
  !This function returns how many entries there are in the most populated row

  class(t_cds), intent(in)           ::  this !The matrix being examined
  integer                                 ::  ii, jj !Generic counting parameters

  n_rowentriesmax=0
  rowmaxouter: do ii=1, this%ndiag
    rowmaxinner: do jj=ii,this%ndiag
      if (this%distance(jj)-this%distance(ii)>=this%n_row) then
        n_rowentriesmax=max(n_rowentriesmax, jj-ii)
        cycle rowmaxouter
      end if
    end do rowmaxinner
    n_rowentriesmax=max(n_rowentriesmax, this%ndiag-ii+1)
  end do rowmaxouter

  !n_rowentriesmax=n_rowentriesmax+1

end function  cds_max_entries_row

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

subroutine operate_cds(this, vector_in, vector_out)
  class(t_cds), intent(in)               ::  this !The matrix to multiply the vector by
  real(kind=dp), dimension(:), intent(in)     ::  vector_in  !The vector to be multiplied
  real(kind=dp), dimension(:), intent(inout)  ::  vector_out  !The vector which results from the multiplication
  integer                                     ::  column !Temporarily stores the column of a variable
  integer                                     ::  ii, jj  !Generic counting variables

  !First, check if the vector to be multiplied has the same number of rows as the matrix
  if (size(vector_in).ne.this%n_row) then
    print*, "cds_multiply_vector has been given a vector of size ", size(vector_in), " to multiply which is of a different size to the matrix which is of size ", this%n_row, ". Terminating."
    stop
  end if

  !Set the output to zero start with, then add up all contributions
  vector_out=0.0_dp
  do ii=1, this%n_row
    do jj=1, this%ndiag
      !Caluculate column of element to see if the element is within the matrix
      column=this%distance(jj)+ii
      if (column.gt.0 .and. column .le. this%n_row)then
        !Perform the multiplication
        vector_out(ii)=vector_out(ii)+this%values(jj,ii)*vector_in(column)
      end if
    end do
  end do

end subroutine operate_cds

end module CDS_Mod