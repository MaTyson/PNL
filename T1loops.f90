module globais
    integer np, dm
    double precision M
end module globais
!Trabalho 1
program t1
    use globais
    implicit none
    integer nn
    parameter(nn=1.d7)
    double precision l(nn), u(nn)
    double precision x(nn), g(nn), xn(nn), gn(nn), fm(nn)
    double precision alfa, eps, f, seed, z, melhorf
    integer n, itmax, nafmax, j, i, k, Q

    write(*,*)'Entre com dimensao e numero de pontos'
    read(*,*) dm, np

    n=dm*np
    seed=18169803
    do j=1,n
        call rando(seed,z)
        x(j)=2*z-1
    end do

    eps=1.d-10
    itmax=1.d7
    nafmax=1.d7
    alfa=1.d-2
    M=1.d99
    Q=10

    do j=1,n
        l(j)=-1
        u(j)=1
    end do
    call spg(n,x,xn,g,gn,l,u,alfa,eps,itmax,nafmax,f,fm,Q)
    melhorf=f
    do j=1,10000
    do i=1,n
        call rando(seed,z)
        x(i)=2*z-1
    end do
    call spg(n,x,xn,g,gn,l,u,alfa,eps,itmax,nafmax,f,fm,Q)
    melhorf= dmin1(melhorf,f)
    end do
    write(*,*)'o melhor f=',melhorf
    stop
end program t1

!subrotina função inverso da distancia ao quadrado
subroutine fun(n,x,f)
    use globais
    implicit none
    integer n, i, j, k
    double precision x(n), f, dist

    f=0.d0

    do i=1,np-1
        do j=i+1,np
            dist=0.d0
            do k=1,dm
                dist=dist+(x((i-1)*dm+k)-x((j-1)*dm+k))**2
            end do
            if(dist.le.(1.d-12))then
                f=f+M
            else
                f=f+sqrt(1/dist)
            end if
        end do
    end do
    return
    end

!subrotina para o gradiente
subroutine grad(n,x,g)
    use globais
    implicit none
    integer i, j, k, n
    double precision x(n), g(n), dist

    do i=1,n
        g(i)=0.d0
    end do

    do i=1,np-1
        do j=i+1,np
            dist=0.d0
            do k=1,dm
                dist=dist+(x((i-1)*dm+k)-x((j-1)*dm+k))**2
            end do
            do k=1,dm
                if(dist.le.(1.d-12))then
                    g((i-1)*dm+k)=g((i-1)*dm+k)-M
                    g((j-1)*dm+k)=g((j-1)*dm+k)+M
                else
                    g((i-1)*dm+k)=g((i-1)*dm+k)-(x((i-1)*dm+k)-x((j-1)*dm+k))/(dist**(3/2))
                    g((j-1)*dm+k)=g((j-1)*dm+k)+(x((i-1)*dm+k)-x((j-1)*dm+k))/(dist**(3/2))
                end if
            end do
        end do
    end do
    return
    end

!subrotina SPG com decréscimo nao monótono
subroutine spg(n,x,xn,g,gn,l,u,alfa,eps,itmax,nafmax,f,fm,Q)
    implicit none
    integer k, n, j, itmax, nafmax, naf, Q
    double precision eps, alfa
    double precision gnor, num, den, p, a, f, fn, t, fmor
    double precision x(n), xn(n), g(n), gn(n), u(n), l(n), fm(Q)
    
    !inicializando contagem
    k=0
    naf=1

    do j=1,Q
        fm(j)=0
    end do

    call fun(n,x,f)
    call grad(n,x,g)
    !iniciamos t0
    t=1.d0
    fmor=f
    

    !Calculamos norma 1 do gradiente projetado
10  do j=1,n 
        xn(j) = x(j) - t*g(j)
        !projetamos cada coordenada na caixa
        xn(j) = dmax1(l(j), dmin1(xn(j),u(j)))
    end do 

    gnor = 0.0d0
    do j=1,n
        gnor = dmax1(gnor, dabs(xn(j)-x(j)))
    end do
 
    !testamos se gnor menor que eps
    if(gnor.le.eps) then
        !write(*,*)'paro feliz, norma do grad projetado pequena'
        !write(*,*)'foram feitas',k,'iteracoes'
        !write(*,*)'x=',x(1),'...',x(n)
        !write(*,*)'f=',f
        return
    end if

20  do j=1,n
        xn(j)=x(j)-t*g(j)
        xn(j)=dmax1(l(j),dmin1(xn(j),u(j)))
    end do
    
    !calculamos o produto <g,xn-x>
    p=0
    do j=1,n
        p=p+g(j)*(xn(j)-x(j))
    end do
    
    !testamos f(xn)
    call fun(n,xn,fn)
    naf=naf+1

    if(k.gt.1)then
      if(k.le.Q) then
        do j=1,k
            fmor=dmax1(fmor,fm(j))
        end do
      else
        do j=1,Q
            fmor=dmax1(fmor,fm(j))
        end do
      end if
    end if
    
    !Checamos avaliacoes da f
    if(naf.ge.nafmax) then
        !write(*,*)'paro triste, avaliamos f demais'
        !write(*,*)'foram feitas',k,'iteracoes e',naf,'avaliacoes'
        !write(*,*)'x=',x(1),'...',x(n)
        !write(*,*)'f=',f
        return
    end if

    if(fn.le.(fmor + t*alfa*p))then
        if(k.eq.0)then
            t=1.d0  
        else
            num=0
            den=0  
            do j=1,n
                num=num+(xn(j)-x(j))**2
                den=den+(xn(j)-x(j))*(g(j)-gn(j))
            end do
            if(den.eq.(0.d0)) then 
                t=1.d0
            else 
                a=num/den            
                t=dabs(a)
            end if
            if(t.gt.(1.d20)) t=1.d0 
            if(t.lt.(1.d-2)) t=1.d0
        end if       
        do j=1,n
            x(j)=xn(j)
        end do
        do j=1, n
            gn(j)=g(j)
        end do            
        f=fn
        call grad(n,x,g)
        k=k+1

        if(k.le.Q) then 
            fm(k)=fn 
        else
            do j=1,Q-1
                fm(j)=fm(j+1)
            end do
            fm(Q)=fn 
        end if 
        !Checamos iteracoes    
        if(k.ge.itmax) then
            !write(*,*)'paro triste, iteracoes demais'
            !write(*,*)'foram feitas',k,'iteracoes'
            !write(*,*)'x=',x(1),'...',x(n)
            !write(*,*)'f=',f
            return
        end if

        go to 10

    end if

    t=t/2

    go to 20

    end

subroutine rando(seed, x)

!     This is the random number generator of Schrage:
!
!     L. Schrage, A more portable Fortran random number generator, ACM
!     Transactions on Mathematical Software 5 (1979), 132-138.

      double precision seed, x

      double precision a,p,b15,b16,xhi,xalo,leftlo,fhi,k
      data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/

      xhi= seed/b16
      xhi= xhi - dmod(xhi,1.d0)
      xalo= (seed-xhi*b16)*a
      leftlo= xalo/b16
      leftlo= leftlo - dmod(leftlo,1.d0)
      fhi= xhi*a + leftlo
      k= fhi/b15
      k= k - dmod(k,1.d0)
      seed= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if (seed.lt.0) seed = seed + p
      x = seed*4.656612875d-10

      return

      end




