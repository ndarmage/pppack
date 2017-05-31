subroutine interv ( xt, lxt, x, left, mflag )

!*****************************************************************************80
!
!! INTERV brackets a real value in an ascending vector of values.
!
!  Discussion:
!
!    The XT array is a set of increasing values.  The goal of the routine
!    is to determine the largest index I so that 
!
!      XT(I) < XT(LXT)  and  XT(I) <= X.
!
!    The routine is designed to be efficient in the common situation
!    that it is called repeatedly, with X taken from an increasing
!    or decreasing sequence.
!
!    This will happen when a piecewise polynomial is to be graphed.
!    The first guess for LEFT is therefore taken to be the value
!    returned at the previous call and stored in the local variable ILO.
!
!    A first check ascertains that ILO < LXT.  This is necessary
!    since the present call may have nothing to do with the previous
!    call.  Then, if 
!      XT(ILO) <= X < XT(ILO+1), 
!    we set LEFT = ILO and are done after just three comparisons.
!
!    Otherwise, we repeatedly double the difference ISTEP = IHI - ILO
!    while also moving ILO and IHI in the direction of X, until
!      XT(ILO) <= X < XT(IHI)
!    after which we use bisection to get, in addition, ILO + 1 = IHI.
!    The value LEFT = ILO is then returned.
!
!    Thanks to Daniel Gloger for pointing out an important modification
!    to the routine, so that the piecewise polynomial in B-form is
!    left-continuous at the right endpoint of the basic interval,
!    17 April 2014.
!
!  Modified:
!
!    17 April 2014
!
!  Author:
!
!    Carl de Boor
!
!  Reference:
!
!    Carl de Boor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XT(LXT), a nondecreasing sequence of values.
!
!    Input, integer ( kind = 4 ) LXT, the dimension of XT.
!
!    Input, real ( kind = 8 ) X, the point whose location with 
!    respect to the sequence XT is to be determined.
!
!    Output, integer ( kind = 4 ) LEFT, the index of the bracketing value:
!      1     if             X  <  XT(1)
!      I     if   XT(I)  <= X  < XT(I+1), for all I, with
!    LXT-1   if   XT(I)  <= X == XT(I+1) == XT(LXT)
!
!    Output, integer ( kind = 4 ) MFLAG, indicates whether X lies within the
!    range of the data.
!    -1:            X  <  XT(1)
!     0: XT(I)   <= X  < XT(I+1) for all I, with XT(I==LXT) == X
!    +1: XT(LXT) <  X
!
  implicit none

  integer ( kind = 4 ), intent(in) :: lxt
  real ( kind = 8 ), intent(in) :: x, xt(lxt)
  integer ( kind = 4 ), intent(out) :: left, mflag
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ), save :: ilo = 1
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) middle

  mflag = 0
  ihi = ilo + 1

  if ( lxt <= ilo ) then

    if ( xt(lxt) <= x ) then
      go to 110
    end if

    if ( lxt <= 1 ) then
      mflag = -1
      ilo = 1
      left = ilo
      return
    end if

    ilo = lxt - 1
    ihi = lxt

    if ( xt(ilo) <= x ) then
      left = ilo
      return
    end if

  end if

  if ( xt(ihi) > x ) then

    if ( xt(ilo) <= x ) then
      left = ilo
      return
    end if
!
!  Now X < XT(ILO).  Decrease ILO to capture X.
!
    istep = 1
    ihi = ilo
    ilo = ihi - istep
    do while ( ilo > 1 )
      if ( xt(ilo) <= x ) then
        go to 50
      end if
      istep = istep * 2
      ihi = ilo
      ilo = ihi - istep
    end do
    ilo = 1
    if ( x < xt(1) ) then
      mflag = -1
      left = ilo
      return
    end if
    go to 50
!
!  Now XT(IHI) <= X.  Increase IHI to capture X.
!
  else

    istep = 1
    ilo = ihi
    ihi = ilo + istep
    do while ( ihi < lxt )
      if ( x < xt(ihi) ) then
        go to 50
      end if
      istep = istep * 2
      ilo = ihi
      ihi = ilo + istep
    end do
    if ( xt(lxt) <= x ) then
      go to 110
    end if
    ihi = lxt

  end if
!
!  Now XT(ILO) < = X < XT(IHI).  Narrow the interval.
!
50 continue

  middle = ( ilo + ihi ) / 2
  
  do while ( middle /= ilo )
!
!  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
!
    if ( xt(middle) <= x ) then
      ilo = middle
    else
      ihi = middle
    end if

    middle = ( ilo + ihi ) / 2

  end do

  left = ilo
  return

!
!  Set output and return (for the last interval).
!
110 continue

  mflag = 1
!
!  Enforce left continuity at the right endpoint of the basic interval
!
  if ( x == xt(lxt) ) then
    mflag = 0
  end if
  
  left = lxt - 1
  do while ( xt(left) >= xt(lxt) )
    left = left - 1
    if (left == 0) then
      exit
    end if
  end do
  ilo = left

  return
end
