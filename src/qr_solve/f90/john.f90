program main

!*****************************************************************************80
!
!! MAIN is the main program for QR_SOLVE_PRB.
!
!  Discussion:
!
!    QR_SOLVE_PRB tests the QR_SOLVE library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(15)
  real ( kind = 8 ) b(5)
  integer ( kind = 4 ) jpvt(3)
  integer ( kind = 4 ) kr
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) qraux(3)
  real ( kind = 8 ) tol
  real ( kind = 8 ) work(3)
  real ( kind = 8 ) x(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QR_SOLVE_TEST'
  write ( *, '(a)' ) '  QR_SOLVE is a function with a simple interface which'
  write ( *, '(a)' ) '  solves a linear system A*x = b in the least squares sense.'
  write ( *, '(a)' ) '  Compare a tabulated solution X1 to the QR_SOLVE result X2.'

  a = (/ 1.0, 1.0, 1.0, 1.0, 1.0, &
         1.0, 2.0, 3.0, 4.0, 5.0, &
         1.0, 4.0, 9.0, 16.0, 25.0 /)
  lda = 5
  m = 5
  n = 3
  tol = 1.0D-10
  b = (/ 1.0, 2.3, 4.6, 3.1, 1.2 /)

  call dqrank ( a, lda, m, n, tol, kr, jpvt, qraux, work )

  write ( *, * ) kr
  write ( *, * ) jpvt
  write ( *, * ) qraux

  stop
end
subroutine dqrank ( a, lda, m, n, tol, kr, jpvt, qraux, work )

!*****************************************************************************80
!
!! DQRANK computes the QR factorization of a rectangular matrix.
!
!  Discussion:
!
!    This routine is used in conjunction with sqrlss to solve
!    overdetermined, underdetermined and singular linear systems
!    in a least squares sense.
!
!    DQRANK uses the LINPACK subroutine DQRDC to compute the QR
!    factorization, with column pivoting, of an M by N matrix A.
!    The numerical rank is determined using the tolerance TOL.
!
!    Note that on output, ABS ( A(1,1) ) / ABS ( A(KR,KR) ) is an estimate
!    of the condition number of the matrix of independent columns,
!    and of R.  This estimate will be <= 1/TOL.
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, the matrix whose
!    decomposition is to be computed.  On output, the information from DQRDC.
!    The triangular matrix R of the QR factorization is contained in the
!    upper triangle and information needed to recover the orthogonal
!    matrix Q is stored below the diagonal in A and in the vector QRAUX.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must
!    be at least M.
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input, real ( kind = 8 ) TOL, a relative tolerance used to determine the
!    numerical rank.  The problem should be scaled so that all the elements
!    of A have roughly the same absolute accuracy, EPS.  Then a reasonable
!    value for TOL is roughly EPS divided by the magnitude of the largest
!    element.
!
!    Output, integer ( kind = 4 ) KR, the numerical rank.
!
!    Output, integer ( kind = 4 ) JPVT(N), the pivot information from DQRDC.
!    Columns JPVT(1), ..., JPVT(KR) of the original matrix are linearly
!    independent to within the tolerance TOL and the remaining columns
!    are linearly dependent.
!
!    Output, real ( kind = 8 ) QRAUX(N), will contain extra information defining
!    the QR factorization.
!
!    Workspace, real ( kind = 8 ) WORK(N).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jpvt(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kr
  integer ( kind = 4 ) m
  real ( kind = 8 ) qraux(n)
  real ( kind = 8 ) tol
  real ( kind = 8 ) work(n)

  jpvt(1:n) = 0

  call dqrdc ( a, lda, m, n, qraux, jpvt, work, 1 )

  kr = 0
  k = min ( m, n )

  do j = 1, k
    if ( abs ( a(j,j) ) <= tol * abs ( a(1,1) ) ) then
      return
    end if
    kr = j
  end do

  return
end
subroutine dqrdc ( a, lda, n, p, qraux, jpvt, work, job )

!*****************************************************************************80
!
!! DQRDC computes the QR factorization of a real rectangular matrix.
!
!  Discussion:
!
!    DQRDC uses Householder transformations.
!
!    Column pivoting based on the 2-norms of the reduced columns may be
!    performed at the user's option.
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,P).  On input, the N by P matrix
!    whose decomposition is to be computed.  On output, A contains in
!    its upper triangle the upper triangular matrix R of the QR
!    factorization.  Below its diagonal A contains information from
!    which the orthogonal part of the decomposition can be recovered.
!    Note that if pivoting has been requested, the decomposition is not that
!    of the original matrix A but that of A with its columns permuted
!    as described by JPVT.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer ( kind = 4 ) N, the number of rows of the matrix A.
!
!    Input, integer ( kind = 4 ) P, the number of columns of the matrix A.
!
!    Output, real ( kind = 8 ) QRAUX(P), contains further information required
!    to recover the orthogonal part of the decomposition.
!
!    Input/output, integer ( kind = 4 ) JPVT(P).  On input, JPVT contains
!    integers that control the selection of the pivot columns.  The K-th
!    column A(*,K) of A is placed in one of three classes according to the
!    value of JPVT(K).
!      > 0, then A(K) is an initial column.
!      = 0, then A(K) is a free column.
!      < 0, then A(K) is a final column.
!    Before the decomposition is computed, initial columns are moved to
!    the beginning of the array A and final columns to the end.  Both
!    initial and final columns are frozen in place during the computation
!    and only free columns are moved.  At the K-th stage of the
!    reduction, if A(*,K) is occupied by a free column it is interchanged
!    with the free column of largest reduced norm.  JPVT is not referenced
!    if JOB == 0.  On output, JPVT(K) contains the index of the column of the
!    original matrix that has been interchanged into the K-th column, if
!    pivoting was requested.
!
!    Workspace, real ( kind = 8 ) WORK(P).  WORK is not referenced if JOB == 0.
!
!    Input, integer ( kind = 4 ) JOB, initiates column pivoting.
!    0, no pivoting is done.
!    nonzero, pivoting is done.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p

  real ( kind = 8 ) a(lda,p)
  integer ( kind = 4 ) jpvt(p)
  real ( kind = 8 ) qraux(p)
  real ( kind = 8 ) work(p)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lup
  integer ( kind = 4 ) maxj
  real ( kind = 8 ) maxnrm
  real ( kind = 8 ) nrmxl
  integer ( kind = 4 ) pl
  integer ( kind = 4 ) pu
  real ( kind = 8 ) ddot
  real ( kind = 8 ) dnrm2
  logical swapj
  real ( kind = 8 ) t
  real ( kind = 8 ) tt

  pl = 1
  pu = 0
!
!  If pivoting is requested, rearrange the columns.
!
  if ( job /= 0 ) then

    do j = 1, p

      swapj = 0 < jpvt(j)

      if ( jpvt(j) < 0 ) then
        jpvt(j) = - j
      else
        jpvt(j) = j
      end if

      if ( swapj ) then

        if ( j /= pl ) then
          call dswap ( n, a(1,pl), 1, a(1,j), 1 )
        end if

        jpvt(j) = jpvt(pl)
        jpvt(pl) = j
        pl = pl + 1

      end if

    end do

    pu = p

    do j = p, 1, -1

      if ( jpvt(j) < 0 ) then

        jpvt(j) = - jpvt(j)

        if ( j /= pu ) then
          call dswap ( n, a(1,pu), 1, a(1,j), 1 )
          jp = jpvt(pu)
          jpvt(pu) = jpvt(j)
          jpvt(j) = jp
        end if

        pu = pu - 1

      end if

    end do

  end if
!
!  Compute the norms of the free columns.
!
  do j = pl, pu
    qraux(j) = dnrm2 ( n, a(1,j), 1 )
  end do

  work(pl:pu) = qraux(pl:pu)
!
!  Perform the Householder reduction of A.
!
  lup = min ( n, p )

  do l = 1, lup
!
!  Bring the column of largest norm into the pivot position.
!
    if ( pl <= l .and. l < pu ) then

      maxnrm = 0.0D+00
      maxj = l
      do j = l, pu
        if ( maxnrm < qraux(j) ) then
          maxnrm = qraux(j)
          maxj = j
        end if
      end do

      if ( maxj /= l ) then
        call dswap ( n, a(1,l), 1, a(1,maxj), 1 )
        qraux(maxj) = qraux(l)
        work(maxj) = work(l)
        jp = jpvt(maxj)
        jpvt(maxj) = jpvt(l)
        jpvt(l) = jp
      end if

    end if
!
!  Compute the Householder transformation for column L.
!
    qraux(l) = 0.0D+00

    if ( l /= n ) then

      nrmxl = dnrm2 ( n-l+1, a(l,l), 1 )

      if ( nrmxl /= 0.0D+00 ) then

        if ( a(l,l) /= 0.0D+00 ) then
          nrmxl = sign ( nrmxl, a(l,l) )
        end if

        call dscal ( n-l+1, 1.0D+00 / nrmxl, a(l,l), 1 )
        a(l,l) = 1.0D+00 + a(l,l)
!
!  Apply the transformation to the remaining columns, updating the norms.
!
        do j = l + 1, p

          t = - ddot ( n-l+1, a(l,l), 1, a(l,j), 1 ) / a(l,l)
          call daxpy ( n-l+1, t, a(l,l), 1, a(l,j), 1 )

          if ( pl <= j .and. j <= pu ) then

            if ( qraux(j) /= 0.0D+00 ) then

              tt = 1.0D+00 - ( abs ( a(l,j) ) / qraux(j) )**2
              tt = max ( tt, 0.0D+00 )
              t = tt
              tt = 1.0D+00 + 0.05D+00 * tt * ( qraux(j) / work(j) )**2

              if ( tt /= 1.0D+00 ) then
                qraux(j) = qraux(j) * sqrt ( t )
              else
                qraux(j) = dnrm2 ( n-l, a(l+1,j), 1 )
                work(j) = qraux(j)
              end if

            end if

          end if

        end do
!
!  Save the transformation.
!
        qraux(l) = a(l,l)
        a(l,l) = - nrmxl

      end if

    end if

  end do

  return
end
subroutine dscal ( n, sa, x, incx )

!*****************************************************************************80
!
!! DSCAL scales a vector by a constant.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) SA, the multiplier.
!
!    Input/output, real ( kind = 8 ) X(*), the vector to be scaled.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) sa
  real ( kind = 8 ) x(*)

  if ( n <= 0 ) then

  else if ( incx == 1 ) then

    m = mod ( n, 5 )

    x(1:m) = sa * x(1:m)

    do i = m+1, n, 5
      x(i)   = sa * x(i)
      x(i+1) = sa * x(i+1)
      x(i+2) = sa * x(i+2)
      x(i+3) = sa * x(i+3)
      x(i+4) = sa * x(i+4)
    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    do i = 1, n
      x(ix) = sa * x(ix)
      ix = ix + incx
    end do

  end if

  return
end
subroutine dswap ( n, x, incx, y, incy )

!*****************************************************************************80
!
!! DSWAP interchanges two vectors.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input/output, real ( kind = 8 ) X(*), one of the vectors to swap.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
!    Input/output, real ( kind = 8 ) Y(*), one of the vectors to swap.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive 
!    elements of Y.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) temp
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    m = mod ( n, 3 )

    do i = 1, m
      temp = x(i)
      x(i) = y(i)
      y(i) = temp
    end do

    do i = m + 1, n, 3

      temp = x(i)
      x(i) = y(i)
      y(i) = temp

      temp = x(i+1)
      x(i+1) = y(i+1)
      y(i+1) = temp

      temp = x(i+2)
      x(i+2) = y(i+2)
      y(i+2) = temp

    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      temp = x(ix)
      x(ix) = y(iy)
      y(iy) = temp
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine daxpy ( n, da, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DAXPY computes constant times a vector plus a vector.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!    This routine uses unrolled loops for increments equal to one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in DX and DY.
!
!    Input, real ( kind = 8 ) DA, the multiplier of DX.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of DX.
!
!    Input/output, real ( kind = 8 ) DY(*), the second vector.
!    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive 
!    entries of DY.
!
  implicit none

  real ( kind = 8 ) da
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( n <= 0 ) then
    return
  end if

  if ( da == 0.0D+00 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dy(iy) = dy(iy) + da * dx(ix)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 4 )

    dy(1:m) = dy(1:m) + da * dx(1:m)

    do i = m+1, n, 4
      dy(i  ) = dy(i  ) + da * dx(i  )
      dy(i+1) = dy(i+1) + da * dx(i+1)
      dy(i+2) = dy(i+2) + da * dx(i+2)
      dy(i+3) = dy(i+3) + da * dx(i+3)
    end do

  end if

  return
end
function ddot ( n, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DDOT forms the dot product of two vectors.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!    This routine uses unrolled loops for increments equal to one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries in DX.
!
!    Input, real ( kind = 8 ) DY(*), the second vector.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive 
!    entries in DY.
!
!    Output, real ( kind = 8 ) DDOT, the sum of the product of the 
!    corresponding entries of DX and DY.
!
  implicit none

  real ( kind = 8 ) ddot
  real ( kind = 8 ) dtemp
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  ddot = 0.0D+00
  dtemp = 0.0D+00

  if ( n <= 0 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dtemp = dtemp + dx(ix) * dy(iy)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 5 )

    do i = 1, m
      dtemp = dtemp + dx(i) * dy(i)
    end do

    do i = m+1, n, 5

      dtemp = dtemp + dx(i  ) * dy(i  ) &
                    + dx(i+1) * dy(i+1) &
                    + dx(i+2) * dy(i+2) &
                    + dx(i+3) * dy(i+3) &
                    + dx(i+4) * dy(i+4)
    end do

  end if

  ddot = dtemp

  return
end
function dnrm2 ( n, x, incx )

!*****************************************************************************80
!
!! DNRM2 returns the euclidean norm of a vector.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!     DNRM2 ( X ) = sqrt ( X' * X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
!    David Kincaid, Fred Krogh.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector whose norm is to be computed.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
!    Output, real ( kind = 8 ) DNRM2, the Euclidean norm of X.
!
  implicit none

  real ( kind = 8 ) absxi
  real ( kind = 8 ) dnrm2
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm
  real ( kind = 8 ) scale
  real ( kind = 8 ) ssq
  real ( kind = 8 ) x(*)

  if ( n < 1 .or. incx < 1 ) then

    norm  = 0.0D+00

  else if ( n == 1 ) then

    norm  = abs ( x(1) )

  else

    scale = 0.0D+00
    ssq = 1.0D+00

    do ix = 1, 1 + ( n - 1 ) * incx, incx
      if ( x(ix) /= 0.0D+00 ) then
        absxi = abs ( x(ix) )
        if ( scale < absxi ) then
          ssq = 1.0D+00 + ssq * ( scale / absxi )**2
          scale = absxi
        else
          ssq = ssq + ( absxi / scale )**2
        end if
      end if
    end do
    norm  = scale * sqrt ( ssq )
  end if

  dnrm2 = norm

  return
end
