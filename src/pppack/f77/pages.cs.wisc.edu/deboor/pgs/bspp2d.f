      subroutine bspp2d ( t, bcoef, n, k, m, scrtch, break, coef)       ********    
c  from  * a practical guide to splines *  by c. de Boor    
calls  bsplvb     
c  this is an extended version of  bsplpp  for use with tensor products --------    
c     
converts the b-representation  t, bcoef(.,j), n, k  of some spline into --------    
c  its pp-representation  break, coef(j,.,.), l, k , j=1, ..., m  .     --------    
c     
c******  i n p u t  ******    
c  t.....knot sequence, of length  n+k    
c  bcoef(.,j) b-spline coefficient sequence, of length  n ,j=1,...,m    --------    
c  n.....length of  bcoef  and  dimension of spline space  spline(k,t)  
c  k.....order of the spline  
c     
c  w a r n i n g  . . .  the restriction   k .le. kmax (= 20)   is impo-
c        sed by the arbitrary dimension statement for  biatx  below, but
c        is  n o w h e r e   c h e c k e d   for.     
c     
c  m     number of data sets                                            ********    
c     
c******  w o r k   a r e a  ******  
c  scrtch   of size  (k,k,m), needed to contain bcoeffs of a piece of   ********    
c        the spline and its  k-1  derivatives   for each of the m sets  --------    
c     
c******  o u t p u t  ******  
c  break.....breakpoint sequence, of length  l+1, contains (in increas- 
c        ing order) the distinct points in the sequence  t(k),...,t(n+1)
c  coef(mm,.,.)  array of size (k,n), with  coef(mm,i,j) = (i-1)st der- ********    
c        ivative of  mm-th  spline at break(j) from the right, mm=1,.,m ********    
c     
c******  m e t h o d  ******  
c     for each breakpoint interval, the  k  relevant b-coeffs of the    
c  spline are found and then differenced repeatedly to get the b-coeffs 
c  of all the derivatives of the spline on that interval. the spline and
c  its first  k-1  derivatives are then evaluated at the left end point 
c  of that interval, using  bsplvb  repeatedly to obtain the values of  
c  all b-splines of the appropriate order at that point.    
c     
      integer k,m,n,   i,j,jp1,kmax,kmj,left,lsofar,mm                  ********
      parameter (kmax = 20)     
      real bcoef(n,m),break(n+2-k),coef(m,k,n+1-k),t(n+k)               ********
     *                         ,scrtch(k,k,m),biatx(kmax),diff,fkmj,sum
c
      lsofar = 0  
      break(1) = t(k)   
      do 50 left=k,n    
c                                find the next nontrivial knot interval.
         if (t(left+1) .eq. t(left))    go to 50
         lsofar = lsofar + 1  
         break(lsofar+1) = t(left+1)
         if (k .gt. 1)                  go to 9 
         do 5 mm=1,m                                                    ********    
    5       coef(mm,1,lsofar) = bcoef(left,mm)                          ********    
                                        go to 50
c        store the k b-spline coeff.s relevant to current knot interval 
c                             in  scrtch(.,1) . 
    9    do 10 i=1,k    
            do 10 mm=1,m                                                ********    
   10          scrtch(i,1,mm) = bcoef(left-k+i,mm)                      ********    
c     
c        for j=1,...,k-1, compute the  k-j  b-spline coeff.s relevant to
c        current knot interval for the j-th derivative by differencing  
c        those for the (j-1)st derivative, and store in scrtch(.,j+1) . 
         do 20 jp1=2,k  
            j = jp1 - 1 
            kmj = k - j 
            fkmj = float(kmj) 
            do 20 i=1,kmj     
               diff = (t(left+i) - t(left+i - kmj))/fkmj                --------    
               if (diff .le. 0.)         go to 20                       --------    
               do 15 mm=1,m                                             ********    
   15             scrtch(i,jp1,mm) =                                    ********    
     *            (scrtch(i+1,j,mm) - scrtch(i,j,mm))/diff              ********    
   20          continue 
c     
c        for  j = 0, ..., k-1, find the values at  t(left)  of the  j+1 
c        b-splines of order  j+1  whose support contains the current    
c        knot interval from those of order  j  (in  biatx ), then comb- 
c        ine with the b-spline coeff.s (in scrtch(.,k-j) ) found earlier
c        to compute the (k-j-1)st derivative at  t(left)  of the given  
c        spline.  
c           note. if the repeated calls to  bsplvb  are thought to gene-
c        rate too much overhead, then replace the first call by   
c           biatx(1) = 1.     
c        and the subsequent call by the statement     
c           j = jp1 - 1 
c        followed by a direct copy of the lines 
c           deltar(j) = t(left+j) - x     
c                  ......     
c           biatx(j+1) = saved
c        from  bsplvb . deltal(kmax)  and  deltar(kmax)  would have to  
c        appear in a dimension statement, of course.  
c     
         call bsplvb ( t, 1, 1, t(left), left, biatx )
         do 25 mm=1,m                                                   ********    
   25       coef(mm,k,lsofar) = scrtch(1,k,mm)                          ********    
         do 30 jp1=2,k  
            call bsplvb ( t, jp1, 2, t(left), left, biatx ) 
            kmj = k+1 - jp1   
            do 30 mm=1,m                                                ********    
               sum = 0.                                                 --------    
               do 28 i=1,jp1                                            --------    
   28             sum = biatx(i)*scrtch(i,kmj,mm) + sum                 ********    
   30          coef(mm,kmj,lsofar) = sum                                ********    
   50    continue 
                                        return  
      end   
