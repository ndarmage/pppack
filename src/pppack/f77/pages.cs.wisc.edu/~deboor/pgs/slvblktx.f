      subroutine slvblk ( bloks, integs, nbloks, b, ipivot, x, iflag )
c    this program solves  the  linear system  a*x = b  where a is an
c  almost block diagonal matrix.  such almost block diagonal matrices
c  arise naturally in piecewise polynomial interpolation or approx-
c  imation and in finite element methods for two-point boundary value
c  problems.  the plu factorization method is implemented here to take
c  advantage of the special structure of such systems for savings in
c  computing time and storage requirements.
c
c                  parameters
c  bloks   a one-dimenional array, of length
c                   sum( integs(1,i)*integs(2,i) ; i = 1,nbloks )
c          on input, contains the blocks of the almost block diagonal
c          matrix  a  .  the array integs (see below and the example)
c          describes the block structure.
c          on output, contains correspondingly the plu factorization
c          of  a  (if iflag .ne. 0).  certain of the entries into bloks
c          are arbitrary (where the blocks overlap).
c  integs  integer array description of the block structure of  a .
c            integs(1,i) = no. of rows of block i        =  nrow
c            integs(2,i) = no. of colums of block i      =  ncol
c            integs(3,i) = no. of elim. steps in block i =  last
c                          i  = 1,2,...,nbloks
c          the linear system is of order
c                n  =  sum ( integs(3,i) , i=1,...,nbloks ),
c          but the total number of rows in the blocks is
c              nbrows = sum( integs(1,i) ; i = 1,...,nbloks)
c  nbloks  number of blocks
c  b       right side of the linear system, array of length nbrows.
c          certain of the entries are arbitrary, corresponding to
c          rows of the blocks which overlap (see block structure and
c          the example below).
c  ipivot  on output, integer array containing the pivoting sequence
c          used. length is nbrows
c  x       on output, contains the computed solution (if iflag .ne. 0)
c          length is n.
c  iflag   on output, integer
c            = (-1)**(no. of interchanges during factorization)
c                   if  a  is invertible
c            = 0    if  a  is singular
c
c                   auxiliary programs
c  fcblok (bloks,integs,nbloks,ipivot,scrtch,iflag)  factors the matrix
c           a , and is used for this purpose in slvblk. its arguments
c          are as in slvblk, except for
c              scrtch = a work array of length max(integs(1,i)).
c
c  sbblok (bloks,integs,nbloks,ipivot,b,x)  solves the system a*x = b
c          once  a  is factored. this is done automatically by slvblk
c          for one right side b, but subsequent solutions may be
c          obtained for additional b-vectors. the arguments are all
c          as in slvblk.
c
c  dtblok (bloks,integs,nbloks,ipivot,iflag,detsgn,detlog) computes the
c          determinant of  a  once slvblk or fcblok has done the fact-
c          orization.the first five arguments are as in slvblk.
c              detsgn  = sign of the determinant
c              detlog  = natural log of the determinant
c
c             ------ block structure of  a  ------
c  the nbloks blocks are stored consecutively in the array  bloks .
c  the first block has its (1,1)-entry at bloks(1), and, if the i-th
c  block has its (1,1)-entry at bloks(index(i)), then
c         index(i+1) = index(i)  +  nrow(i)*ncol(i) .
c    the blocks are pieced together to give the interesting part of  a
c  as follows.  for i = 1,2,...,nbloks-1, the (1,1)-entry of the next
c  block (the (i+1)st block ) corresponds to the (last+1,last+1)-entry
c  of the current i-th block.  recall last = integs(3,i) and note that
c  this means that
c      a. every block starts on the diagonal of  a .
c      b. the blocks overlap (usually). the rows of the (i+1)st block
c         which are overlapped by the i-th block may be arbitrarily de-
c         fined initially. they are overwritten during elimination.
c    the right side for the equations in the i-th block are stored cor-
c  respondingly as the last entries of a piece of  b  of length  nrow
c  (= integs(1,i)) and following immediately in  b  the corresponding
c  piece for the right side of the preceding block, with the right side
c  for the first block starting at  b(1) . in this, the right side for
c  an equation need only be specified once on input, in the first block
c  in which the equation appears.
c
c             ------ example and test driver ------
c    the test driver for this package contains an example, a linear
c  system of order 11, whose nonzero entries are indicated in the fol-
c  lowing schema by their row and column index modulo 10. next to it
c  are the contents of the  integs  arrray when the matrix is taken to
c  be almost block diagonal with  nbloks = 5, and below it are the five
c  blocks.
c
c                      nrow1 = 3, ncol1 = 4
c           11 12 13 14
c           21 22 23 24   nrow2 = 3, ncol2 = 3
c           31 32 33 34
c  last1 = 2      43 44 45
c                 53 54 55            nrow3 = 3, ncol3 = 4
c        last2 = 3         66 67 68 69   nrow4 = 3, ncol4 = 4
c                          76 77 78 79      nrow5 = 4, ncol5 = 4
c                          86 87 88 89
c                 last3 = 1   97 98 99 90
c                    last4 = 1   08 09 00 01
c                                18 19 10 11
c                       last5 = 4
c
c         actual input to bloks shown by rows of blocks of  a .
c      (the ** items are arbitrary, this storage is used by slvblk)
c
c  11 12 13 14  / ** ** **  / 66 67 68 69  / ** ** ** **  / ** ** ** **
c  21 22 23 24 /  43 44 45 /  76 77 78 79 /  ** ** ** ** /  ** ** ** **
c  31 32 33 34/   53 54 55/   86 87 88 89/   97 98 99 90/   08 09 00 01
c                                                           18 19 10 11
c
c  index = 1      index = 13  index = 22     index = 34     index = 46
c
c         actual right side values with ** for arbitrary values
c  b1 b2 b3 ** b4 b5 b6 b7 b8 ** ** b9 ** ** b10 b11
c
c  (it would have been more efficient to combine block 3 with block 4)
c
      integer integs(3,nbloks),ipivot(1),iflag
      real bloks(1),b(1),x(1)
c     in the call to fcblok,  x  is used for temporary storage.
      call fcblok(bloks,integs,nbloks,ipivot,x,iflag)
      if (iflag .eq. 0)                 return
      call sbblok(bloks,integs,nbloks,ipivot,b,x)
                                        return
      end
      subroutine fcblok ( bloks, integs, nbloks, ipivot, scrtch, iflag )
calls subroutines  f a c t r b  and  s h i f t b .
c
c   f c b l o k  supervises the plu factorization with pivoting of
c  scaled rows of the almost block diagonal matrix stored in the arrays
c   b l o k s  and  i n t e g s .
c
c   factrb = subprogram which carries out steps 1,...,last of gauss
c            elimination (with pivoting) for an individual block.
c   shiftb = subprogram which shifts the remaining rows to the top of
c            the next block
c
c parameters
c    bloks   an array that initially contains the almost block diagonal
c            matrix  a  to be factored, and on return contains the com-
c            puted factorization of  a .
c    integs  an integer array describing the block structure of  a .
c    nbloks  the number of blocks in  a .
c    ipivot  an integer array of dimension  sum (integs(1,n) ; n=1,
c            ...,nbloks) which, on return, contains the pivoting stra-
c            tegy used.
c    scrtch  work area required, of length  max (integs(1,n) ; n=1,
c            ...,nbloks).
c    iflag   output parameter;
c            = 0  in case matrix was found to be singular.
c            otherwise,
c            = (-1)**(number of row interchanges during factorization)
c
      integer integs(3,nbloks),ipivot(1),iflag, i,index,indexb,indexn,
     *        last,ncol,nrow
      real bloks(1),scrtch(1)
      iflag = 1
      indexb = 1
      indexn = 1
      i = 1
c                        loop over the blocks.  i  is loop index
   10    index = indexn
         nrow = integs(1,i)
         ncol = integs(2,i)
         last = integs(3,i)
c        carry out elimination on the i-th block until next block
c        enters, i.e., for columns 1,...,last  of i-th block.
         call factrb(bloks(index),ipivot(indexb),scrtch,nrow,ncol,last,
     *            iflag)
c         check for having reached a singular block or the last block
         if (iflag .eq. 0 .or. i .eq. nbloks)
     *                                  return
         i = i+1
         indexn = nrow*ncol + index
c              put the rest of the i-th block onto the next block
         call shiftb(bloks(index),ipivot(indexb),nrow,ncol,last,
     *            bloks(indexn),integs(1,i),integs(2,i))
         indexb = indexb + nrow
                                        go to 10
      end
      subroutine factrb ( w, ipivot, d, nrow, ncol, last, iflag )
c  adapted from p.132 of 'element.numer.analysis' by conte-de boor
c
c  constructs a partial plu factorization, corresponding to steps 1,...,
c   l a s t   in gauss elimination, for the matrix  w  of order
c   ( n r o w ,  n c o l ), using pivoting of scaled rows.
c
c  parameters
c    w       contains the (nrow,ncol) matrix to be partially factored
c            on input, and the partial factorization on output.
c    ipivot  an integer array of length nrow containing a record of the
c            pivoting strategy used; row ipivot(i) is used during the
c            i-th elimination step, i=1,...,last.
c    d       a work array of length nrow used to store row sizes
c            temporarily.
c    nrow    number of rows of w.
c    ncol    number of columns of w.
c    last    number of elimination steps to be carried out.
c    iflag   on output, equals iflag on input times (-1)**(number of
c            row interchanges during the factorization process), in
c            case no zero pivot was encountered.
c            otherwise, iflag = 0 on output.
c
      integer ipivot(nrow),ncol,last,iflag, i,ipivi,ipivk,j,k,kp1
      real w(nrow,ncol),d(nrow), awikdi,colmax,ratio,rowmax
c  initialize ipivot, d
      do 10 i=1,nrow
         ipivot(i) = i
         rowmax = 0.
         do 9 j=1,ncol
    9       rowmax = amax1(rowmax, abs(w(i,j)))
         if (rowmax .eq. 0.)            go to 999
   10    d(i) = rowmax
c gauss elimination with pivoting of scaled rows, loop over k=1,.,last
      k = 1
c        as pivot row for k-th step, pick among the rows not yet used,
c        i.e., from rows ipivot(k),...,ipivot(nrow), the one whose k-th
c        entry (compared to the row size) is largest. then, if this row
c        does not turn out to be row ipivot(k), redefine ipivot(k) ap-
c        propriately and record this interchange by changing the sign
c        of  i f l a g .
   11    ipivk = ipivot(k)
         if (k .eq. nrow)               go to 21
         j = k
         kp1 = k+1
         colmax = abs(w(ipivk,k))/d(ipivk)
c              find the (relatively) largest pivot
         do 15 i=kp1,nrow
            ipivi = ipivot(i)
            awikdi = abs(w(ipivi,k))/d(ipivi)
            if (awikdi .le. colmax)     go to 15
               colmax = awikdi
               j = i
   15       continue
         if (j .eq. k)                  go to 16
         ipivk = ipivot(j)
         ipivot(j) = ipivot(k)
         ipivot(k) = ipivk
         iflag = -iflag
   16    continue
c        if pivot element is too small in absolute value, declare
c        matrix to be noninvertible and quit.
         if (abs(w(ipivk,k))+d(ipivk) .le. d(ipivk))
     *                                  go to 999
c        otherwise, subtract the appropriate multiple of the pivot
c        row from remaining rows, i.e., the rows ipivot(k+1),...,
c        ipivot(nrow), to make k-th entry zero. save the multiplier in
c        its place.
         do 20 i=kp1,nrow
            ipivi = ipivot(i)
            w(ipivi,k) = w(ipivi,k)/w(ipivk,k)
            ratio = -w(ipivi,k)
            do 20 j=kp1,ncol
   20          w(ipivi,j) = ratio*w(ipivk,j) + w(ipivi,j)
         k = kp1
c        check for having reached the next block.
         if (k .le. last)               go to 11
                                        return
c     if  last  .eq. nrow , check now that pivot element in last row
c     is nonzero.
   21 if( abs(w(ipivk,nrow))+d(ipivk) .gt. d(ipivk) )
     *                                  return
c                   singularity flag set
  999 iflag = 0
                                        return
      end
      subroutine shiftb ( ai, ipivot, nrowi, ncoli, last,
     *                    ai1, nrowi1, ncoli1 )
c  shifts the rows in current block, ai, not used as pivot rows, if
c  any, i.e., rows ipivot(last+1),...,ipivot(nrowi), onto the first
c  mmax = nrow-last rows of the next block, ai1, with column last+j of
c  ai  going to column j , j=1,...,jmax=ncoli-last. the remaining col-
c  umns of these rows of ai1 are zeroed out.
c
c                             picture
c
c       original situation after         results in a new block i+1
c       last = 2 columns have been       created and ready to be
c       done in factrb (assuming no      factored by next factrb call.
c       interchanges of rows)
c                   1
c              x  x 1x  x  x           x  x  x  x  x
c                   1
c              0  x 1x  x  x           0  x  x  x  x
c  block i          1                       ---------------
c  nrowi = 4   0  0 1x  x  x           0  0 1x  x  x  0  01
c  ncoli = 5        1                       1             1
c  last = 2    0  0 1x  x  x           0  0 1x  x  x  0  01
c  -------------------------------          1             1   new
c                   1x  x  x  x  x          1x  x  x  x  x1  block
c                   1                       1             1   i+1
c  block i+1        1x  x  x  x  x          1x  x  x  x  x1
c  nrowi1= 5        1                       1             1
c  ncoli1= 5        1x  x  x  x  x          1x  x  x  x  x1
c  -------------------------------          1-------------1
c                   1
c
      integer ipivot(nrowi),last, ip,j,jmax,jmaxp1,m,mmax
      real ai(nrowi,ncoli),ai1(nrowi1,ncoli1)
      mmax = nrowi - last
      jmax = ncoli - last
      if (mmax .lt. 1 .or. jmax .lt. 1) return
c              put the remainder of block i into ai1
      do 10 m=1,mmax
         ip = ipivot(last+m)
         do 10 j=1,jmax
   10       ai1(m,j) = ai(ip,last+j)
      if (jmax .eq. ncoli1)             return
c              zero out the upper right corner of ai1
      jmaxp1 = jmax + 1
      do 20 j=jmaxp1,ncoli1
         do 20 m=1,mmax
   20       ai1(m,j) = 0.
                                        return
      end
      subroutine sbblok ( bloks, integs, nbloks, ipivot, b, x )
calls subroutines  s u b f o r  and  s u b b a k .
c
c  supervises the solution (by forward and backward substitution) of
c  the linear system  a*x = b  for x, with the plu factorization of  a
c  already generated in  f c b l o k .  individual blocks of equations
c  are solved via  s u b f o r  and  s u b b a k .
c
c parameters
c    bloks, integs, nbloks, ipivot    are as on return from fcblok.
c    b       the right side, stored corresponding to the storage of
c            the equations. see comments in  s l v b l k  for details.
c    x       solution vector
c
      integer integs(3,nbloks),ipivot(1), i,index,indexb,indexx,j,last,
     *        nbp1,ncol,nrow
      real bloks(1),b(1),x(1)
c
c      forward substitution pass
c
      index = 1
      indexb = 1
      indexx = 1
      do 20 i=1,nbloks
         nrow = integs(1,i)
         last = integs(3,i)
         call subfor(bloks(index),ipivot(indexb),nrow,last,b(indexb),
     *               x(indexx))
         index = nrow*integs(2,i) + index
         indexb = indexb + nrow
   20    indexx = indexx + last
c
c     back substitution pass
c
      nbp1 = nbloks + 1
      do 30 j=1,nbloks
         i = nbp1 - j
         nrow = integs(1,i)
         ncol = integs(2,i)
         last = integs(3,i)
         index = index - nrow*ncol
         indexb = indexb - nrow
         indexx = indexx - last
   30    call subbak(bloks(index),ipivot(indexb),nrow,ncol,last,
     *               x(indexx))
                                        return
      end
      subroutine subfor ( w, ipivot, nrow, last, b, x )
c  carries out the forward pass of substitution for the current block,
c  i.e., the action on the right side corresponding to the elimination
c  carried out in  f a c t r b  for this block.
c     at the end, x(j) contains the right side of the transformed
c  ipivot(j)-th equation in this block, j=1,...,nrow. then, since
c  for i=1,...,nrow-last, b(nrow+i) is going to be used as the right
c  side of equation  i  in the next block (shifted over there from
c  this block during factorization), it is set equal to x(last+i) here.
c
c parameters
c    w, ipivot, nrow, last  are as on return from factrb.
c    b(j)   is expected to contain, on input, the right side of j-th
c           equation for this block, j=1,...,nrow.
c    b(nrow+j)   contains, on output, the appropriately modified right
c           side for equation j in next block, j=1,...,nrow-last.
c    x(j)   contains, on output, the appropriately modified right
c           side of equation ipivot(j) in this block, j=1,...,last (and
c           even for j=last+1,...,nrow).
c
      integer ipivot(nrow), ip,jmax,k
c     dimension b(nrow + nrow-last)
      real w(nrow,last),b(1),x(nrow)
      ip = ipivot(1)
      x(1) = b(ip)
      if (nrow .eq. 1)                  go to 99
      do 15 k=2,nrow
         ip = ipivot(k)
         jmax = amin0(k-1,last)
         sum = 0.
         do 14 j=1,jmax
   14       sum = w(ip,j)*x(j) + sum
   15    x(k) = b(ip) - sum
c
c     transfer modified right sides of equations ipivot(last+1),...,
c     ipivot(nrow) to next block.
      nrowml = nrow - last
      if (nrowml .eq. 0)                go to 99
      lastp1 = last+1
      do 25 k=lastp1,nrow
   25    b(nrowml+k) = x(k)
   99                                   return
      end
      subroutine subbak ( w, ipivot, nrow, ncol, last, x )
c  carries out backsubstitution for current block.
c
c parameters
c    w, ipivot, nrow, ncol, last  are as on return from factrb.
c    x(1),...,x(ncol)  contains, on input, the right side for the
c            equations in this block after backsubstitution has been
c            carried up to but not including equation ipivot(last).
c            means that x(j) contains the right side of equation ipi-
c            vot(j) as modified during elimination, j=1,...,last, while
c            for j .gt. last, x(j) is already a component of the solut-
c            ion vector.
c    x(1),...,x(ncol) contains, on output, the components of the solut-
c            ion corresponding to the present block.
c
      integer ipivot(nrow),last,  ip,j,k,kp1
      real w(nrow,ncol),x(ncol), sum
      k = last
      ip = ipivot(k)
      sum = 0.
      if (k .eq. ncol)                  go to 4
      kp1 = k+1
    2    do 3 j=kp1,ncol
    3       sum = w(ip,j)*x(j) + sum
    4    x(k) = (x(k) - sum)/w(ip,k)
         if (k .eq. 1)                  return
         kp1 = k
         k = k-1
         ip = ipivot(k)
         sum = 0.
                                        go to 2
      end
      subroutine dtblok ( bloks, integs, nbloks, ipivot, iflag,
     *                    detsgn, detlog )
c  computes the determinant of an almost block diagonal matrix whose
c  plu factorization has been obtained previously in fcblok.
c  *** the logarithm of the determinant is computed instead of the
c  determinant itself to avoid the danger of overflow or underflow
c  inherent in this calculation.
c
c parameters
c    bloks, integs, nbloks, ipivot, iflag  are as on return from fcblok.
c            in particular, iflag = (-1)**(number of interchanges dur-
c            ing factorization) if successful, otherwise iflag = 0.
c    detsgn  on output, contains the sign of the determinant.
c    detlog  on output, contains the natural logarithm of the determi-
c            nant if determinant is not zero. otherwise contains 0.
c
      integer integs(3,nbloks),ipivot(1),iflag, i,indexp,ip,k,last
      real bloks(1),detsgn,detlog
c
      detsgn = iflag
      detlog = 0.
      if (iflag .eq. 0)                 return
      index = 0
      indexp = 0
      do 2 i=1,nbloks
         nrow = integs(1,i)
         last = integs(3,i)
         do 1 k=1,last
            ip = index + nrow*(k-1) + ipivot(indexp+k)
            detlog = detlog + alog(abs(bloks(ip)))
    1       detsgn = detsgn*sign(1.,bloks(ip))
         index = nrow*integs(2,i) + index
    2    indexp = indexp + nrow
                                        return
      end
