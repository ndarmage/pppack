chapter xiii, example 3 . test of optimal spline interpolation routine
c                         on titanium heat data .
c  from  * a practical guide to splines *  by c. de boor    
calls titand,splopt(bsplvb,banfac/slv),splint(*),bvalue(interv)
c
C     data n,k /12,5/   
c     lenscr = (n-k)(2k+3)+5k+3  is the length of scrtch required in
c                                splopt .
      integer k,lenscr,n,ntitan
      parameter (n=12,ntitan=49,k=5,lenscr=(n-k)*(2*k+3)+5*k+3)
      integer i,iflag,ipick(n),ipicki,lx,nmk
      real a(n),gtitan(ntitan),gtau(ntitan),scrtch(lenscr),t(n+k),tau(n)
     *    ,x(ntitan)
C     integer i,iflag,ipick(12),ipicki,lx,nmk
C     real a(12),gtitan(49),gtau(49),scrtch(119),t(17),tau(12)
C    *    ,x(49)
      data ipick /1,5,11,21,27,29,31,33,35,40,45,49/
      call titand ( x, gtitan, lx )
      do 10 i=1,n
         ipicki = ipick(i)
         tau(i) = x(ipicki)
   10    gtau(i) = gtitan(ipicki)
      call splopt ( tau, n, k, scrtch, t, iflag )
      if (iflag .gt. 1)                 stop
      call splint ( tau, gtau, t, n, k, scrtch, a, iflag )
      if (iflag .gt. 1)                 stop
      do 20 i=1,lx
         gtau(i) = bvalue ( t, a, n, k, x(i), 0 )
   20    scrtch(i) = gtitan(i) - gtau(i)
      print 620,(i,x(i),gtitan(i),gtau(i),scrtch(i),i=1,lx)
  620 format(41h  i, data point, data, interpolant, error//
     2         (i3,f8.0,f10.4,f9.4,e11.3))
      nmk = n-k
      print 621,(i,t(k+i),i=1,nmk)
  621 format(///16h optimal knots =/(i5,f15.9))
                                        stop
      end
