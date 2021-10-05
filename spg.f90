!subrotina SPG com decréscimo monótono
subroutine spg(n,x,xn,g,gn,l,u,alfa,eps,itmax,nafmax)
    implicit none
    integer k, n, j, itmax, nafmax, naf
    double precision eps, alfa
    double precision gnor, num, den, p, a, f, fn, t
    double precision x(n), xn(n), g(n), gn(n), u(n), l(n)
    
    !inicializando contagem
    k=0
    naf=1

    call fun(n,x,f)
    call grad(n,x,g)
    !iniciamos t0
    t=1.d0

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
        write(*,*)'f=',f
        return
    end if

    if(f.le.(1.d-8))then
        write(*,*)'f=',f        
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
    
    !Checamos avaliacoes da f
    if(naf.ge.nafmax) then
        !write(*,*)'paro triste, avaliamos f demais'
        !write(*,*)'foram feitas',k,'iteracoes e',naf,'avaliacoes'
        !write(*,*)'x=',x(1),'...',x(n)
        write(*,*)'f=',f
        return
    end if

    if(fn.le.(f + t*alfa*p))then
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
        !Checamos iteracoes    
        if(k.ge.itmax) then
            !write(*,*)'paro triste, iteracoes demais'
            !write(*,*)'foram feitas',k,'iteracoes'
            !write(*,*)'x=',x(1),'...',x(n)
            write(*,*)'f=',f
            return
        end if

        go to 10

    end if

    t=t/2

    go to 20

    end
