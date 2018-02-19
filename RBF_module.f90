module RBF_module
  implicit none
contains

  function distance(X_1, Y_1, X_2, Y_2) result(dist)
    real(kind=8), intent(in) :: X_1, Y_1, X_2, Y_2
    real(kind=8)             :: dist
    dist = (X_1 - X_2)**2 + (Y_1 - Y_2)**2
  end function distance

  function gaussian(dist, epsilon) result(output) !not currently in use
    real(kind=8), intent(in)    :: dist
    real(kind=8), intent(in)    :: epsilon
    real(kind=8)                :: output
    output = exp(-epsilon*dist)
  end function gaussian

  subroutine point_set_distance(points_1, points_2, dist_mat)
    !fill dist_mat with distances between point set 1, and 2      (To be optimized)
    real(kind=8), intent(in), dimension(:,:) :: points_1, points_2
    real(kind=8), intent(out), dimension(:,:):: dist_mat
    integer                                  :: i, j, i_max, j_max
    i_max = size(points_1, 1)
    j_max = size(points_2, 1)

    do i = 1,i_max
      do j = 1, j_max
        dist_mat(i,j) = distance(points_1(i, 1), points_1(i, 2), points_2(j, 1), points_2(j, 2))
      end do
    end do
  end subroutine

  subroutine matrix_eval(mat_in, mat_out, epsilon)
    !evaluate function for all elements of a matrix  (Not currently in use)
    real(kind=8), intent(in), dimension(:,:)  :: mat_in
    real(kind=8), intent(out), dimension(:,:) :: mat_out
    real(kind=8), intent(in)                  :: epsilon
    integer                                   :: i_max, j_max, i, j

    i_max = size(mat_in,1)
    j_max = size(mat_in,2)
    do i = 1,i_max
      do j = 1,j_max
        mat_out(i,j) = gaussian(mat_in(i,j), epsilon)
      end do
    end do
    print *, mat_out(1,1)
  end subroutine

  subroutine RBF_matrix(points, point_count, RF_mat, epsilon)
    !Determines proper coefficients for an RBF decomposition
    !over a set of 2d points
    real(kind=8), intent(in), dimension(:,:)                  :: points
    !real(kind=8), intent(in)                  :: RBF_func
    real(kind=8), intent(in)                                  :: epsilon
    real(kind=8), dimension(:,:), allocatable                 :: dist_mat !matrix for storing distances between nodes
    real(kind=8), intent(out)                                 :: RF_mat(point_count, point_count)   !matrix for storing values of RBF
!    real(kind=8), dimension(:)                   :: values   !function to be interpolated evaluated at various points (overwritten with coefficients)
    integer                                                   :: point_count
    integer                                                   :: SUCCESS_FLAG
!f2py intent(in), dimension(:,:) :: points
!f2py real(kind=8), intent(in)       :: epsilon
!f2py integer, intent(in)       :: point_count
!f2py real(kind=8), intent(out), dimension(:,:):: rf_mat
!f2py integer, intent(hide), depend(points)        :: point_count=shape(points, 0)
    allocate(dist_mat(point_count, point_count))
    !allocate(RF_mat(point_count, point_count))
    call point_set_distance(points, points, dist_mat)
    call matrix_eval(dist_mat, RF_mat, epsilon)
    call dpotrf('L', point_count, RF_mat, point_count, SUCCESS_FLAG) !lower triangular cholesky factorization
  end subroutine

  subroutine get_coefficients(rf_mat, values, point_count)
    real(kind=8), intent(in), dimension(:,:)                 :: rf_mat
    real(kind=8), intent(inout), dimension(:)                   :: values
    integer, intent(in)                                      :: point_count
    integer                                                  :: SUCCESS_FLAG
    integer                                                  :: i
!f2py intent(in), dimension(:,:)  :: rf_mat
!f2py real(kind=8), intent(in, out), dimension(:)    :: values
!f2py integer, intent(in), intent(hide)         :: point_count=shape(rf_mat,0)
  do i = 1, point_count
    print *, rf_mat(i,1)
  end do
  call dpotrs('L', point_count, 0,  rf_mat, point_count, values, point_count, SUCCESS_FLAG) !Solve a A*x = b where A is a symmetric matrix (uses the cholesky factorization)
  print *, SUCCESS_FLAG
  end subroutine

  subroutine eval_over_array(interp_points, point_count, coefficients, eval_points, values, epsilon)
    real(kind=8), intent(in), dimension(:,:)   :: interp_points
    real(kind=8), intent(in), dimension(:)     :: coefficients
    real(kind=8), intent(in), dimension(:,:)   :: eval_points
    integer,      intent(in)                   :: point_count
    integer                                    :: interp_point_count
    real(kind=8), intent(out), dimension(point_count)    :: values
    real(kind=8), dimension(:,:), allocatable            :: dist_mat
    !real(kind=8), dimension(:,:), allocatable            :: val_mat
    real(kind=8), intent(in)                             :: epsilon
    integer                                              :: i,j

!f2py intent(in), real(kind=8), dimension(:,:)  :: interp_points
!f2py intent(in), real(kind=8), dimension(:)  :: coefficients
!f2py intent(in), real(kind=8), dimension(:,:)  :: eval_points
!f2py intent(out), real(kind=8), dimension(:) :: values
!f2py integer, intent(hide), depend(eval_points)                             :: point_count=shape(eval_points,0)
!f2py real(kind=8), intent(in)                        :: epsilon

!The key function: evaluate the interpolant over a set of points. Currently two bottle-neck steps.

  interp_point_count = size(interp_points, 1)
  allocate(dist_mat(interp_point_count, point_count)) ! Could be cut down to a one or a few times allocation

  call point_set_distance(interp_points, eval_points, dist_mat) !slow
  
  dist_mat = -epsilon*dist_mat
  
  !calculate the exponential kernel
  dist_mat = exp(dist_mat) !slow    Could possibly be memoized if the same two grids are used several times

  do i = 1, interp_point_count
    dist_mat(i,:) = dist_mat(i,:)*coefficients(i)
  end do
  do j = 1, point_count
    values(j) = sum(dist_mat(:, j))
  end do
  end subroutine
end module RBF_module
