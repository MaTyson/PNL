!definimos as variaveis globais np numero de pontos, dm dimensao da caixa
!M grande para adicionar à função em caso de zero no denominador
module globais  
    integer np, dm
    double precision M
end module globais
!Trabalho 3
program t3
    use globais
    implicit none
    integer nn
    parameter(nn=1.d7)
!iniciamos variaveis
    double precision l(nn), u(nn)
    double precision x(nn), g(nn), xn(nn), gn(nn), fm(nn), lambda(nn), h(nn)
    double precision alfa, eps, f, seed, z, rho
    integer n, itmax, nafmax,i, j, Q
!lemos a dimensao e o numero de pontos
    write(*,*)'Entre com dimensao e o numero de pontos '
    read(*,*) dm, np
!estabelecemos o tamanho de x, que e o produto da dimensao
!com o numero de pontos mais o numero de restriçoes que e np   
    n=dm*np+np
!inicializamos uma semente aleatória qualquer, mas na prática 
!o programa foi rodado em loop para testar varios pontos iniciais
    seed=983456
!inicializamos o ponto inicial x 
    do j=1,dm*np
        call rando(seed,z)
        x(j)=z*2 -1
    end do
    do j=1,np
        x(dm*np+j)=0.d0
    end do
!definimos eps tolerancia de parada, itmax maximo de iteracoes
!nafmax maximo de avaliacoes da funcao objetivo
!alfa parametro de decrescimo para armijo
!M grande em caso de zeros no denominador de f
!Q parametro de decrescimo nao monotono
!rho parametro de penalidade
    eps=1.d-8
    itmax=1.d7
    nafmax=1.d7
    alfa=1.d-2
    M=1.d99
    Q=10
    rho=10
!lambda multiplicador de lagrange
    do i=1,np
      lambda(i)=0.d0
    end do
!definimos a caixa como limitante para os pontos
![-1,1] para os pontos na bola
![0,1] para as folgas
    do j=1,dm*np
        l(j)=-1
        u(j)=1
    end do
    do j=dm*np+1,(dm+1)*np
        l(j)=0
        u(j)=1
    end do
!chamamos lagaum, escrevemos na tela o valor da funcao objetivo
!escrevemos na tela x com os pontos encontrados    
    call lagaum(n,x,xn,g,gn,l,u,alfa,eps,itmax,nafmax,f,fm,Q,lambda,rho,h)  
    stop
end program t3

subroutine fun(n,x,f,h,lambda,rho)
use globais
implicit none
integer n
double precision  x(n),a, hnorm, f, h(np),lambda(np),rho
f=0.d0
call funf(n,x,a)
f=a
call funh(n,x,h,lambda,rho,a)
f=f+a
return
end subroutine

subroutine funh(n,x,h,lambda,rho,v)
use globais
implicit none
integer n, i, j
double precision x(n),lambda(np),rho,v, h(np), hnorm
v=0.d0
do i=1, np
h(i)=0.d0
end do
do i=1, np
  do j=1, dm
    h(i)=h(i)+x((i-1)*dm+j)**2
  end do
  h(i)=h(i)+x(np*dm+i)-1
end do
call nfunh(n,x,h,lambda,rho,hnorm)
v=v+(rho/2)*hnorm**2
do i=1, np
v=v+lambda(i)*h(i)
end do
return
end 

subroutine funf(n,x,f)
    use globais    
    implicit none
    integer n, i, j, k
!declaramos x, f e dist, variavel que auxilia no calculo da distancia
    double precision x(n), f, dist
!iniciamos f=0
    f=0.d0
!calculamos dist    
    do i=1,np-1
        do j=i+1,np
!sempre reinicializamos como zero para cada novos 2 pontos
            dist=0.d0
            do k=1,dm
                dist=dist+(x((i-1)*dm+k)-x((j-1)*dm+k))**2
            end do
!se a distancia ao for muito proxima de zero
!como teriamos um zero no denominador de f
!adicionamos a f o valor M
            if(dist.le.(1.d-12))then
                f=f+M
            else
!senao, adicionamos a f o inverso da raiz quadrada de dist
                f=f+sqrt(1/dist)
            end if
        end do
    end do
    return
end subroutine

subroutine grad(n,x,g,h,lambda,rho)   
    use globais    
    implicit none
    integer i, j, k, n
!novamente usamos dist para auxiliar na distancia
    double precision x(n), g(n), dist, h(np),rho, lambda(np)
!inicializamos g=0    
    do i=1,n
        g(i)=0.d0
    end do
!armazenamos a distancia entre dois pontos pois esta aparece
!no denominador do gradiente
    do i=1,np-1
        do j=i+1,np
!novamente inicializando dist a cada loop como 0
            dist=0.d0
            do k=1,dm
                dist=dist+(x((i-1)*dm+k)-x((j-1)*dm+k))**2
            end do
!se dist for muito proximo de zero, em uma das componentes
!do gradiente subtraimos M, e na outra adicionamos M
            do k=1,dm
                if(dist.le.(1.d-12))then
                    g((i-1)*dm+k)=g((i-1)*dm+k)-M
                    g((j-1)*dm+k)=g((j-1)*dm+k)+M
!senao, adicionamos ao gradiente esta expressao, conforme calculado
                else
                    g((i-1)*dm+k)=g((i-1)*dm+k)-(x((i-1)*dm+k)-x((j-1)*dm+k))/(dist**(3/2))
                    g((j-1)*dm+k)=g((j-1)*dm+k)+(x((i-1)*dm+k)-x((j-1)*dm+k))/(dist**(3/2))
                end if
            end do
        end do
    end do
    do i=1, np
      do j=1, dm
        g((i-1)*dm+j)=g((i-1)*dm+j)+2*lambda(i)*x((i-1)*dm+j)
		  end do
      g(np*dm+i)=g(np*dm+i)+lambda(i)
    end do
    do i=1, np
      dist=0.d0
      do j=1,dm
        dist=dist+x((i-1)*dm+j)**2
      end do
      dist=dist+x(np*dm+i)-1
      do j=1, dm
        g((i-1)*dm+j)=g((i-1)*dm+j)+2*rho*x((i-1)*dm+j)*dist
      end do
    g(np*dm+i)=g(np*dm+i)+rho*dist
    end do
    return
end subroutine

subroutine lagaum(n,x,xn,g,gn,l,u,alfa,eps,itmax,nafmax,f,fm,Q,lambda,rho,h)
use globais
implicit none
integer i,j, n, k, itmax, nafmax, Q
double precision eps, alfa, hnorm, hnormn,rho,norma
double precision x(n),h(np),lambda(np), xn(n),g(n),gn(n),l(n),u(n)
double precision f, fm(Q) 
k=0.d0
10  call funh(n,x,h,lambda,rho,norma)
    call nfunh(n,x,h,lambda,rho,hnorm)
    call spg(n,x,xn,g,gn,l,u,alfa,eps,itmax,nafmax,f,fm,Q,h,lambda,rho)
    call funh(n,x,h,lambda,rho,norma)
    call nfunh(n,x,h,lambda,rho,hnormn)
    if(hnormn.le.((0.1)*hnorm)) then
    do i=1,np
      lambda(i)=lambda(i)+rho*h(i)
    end do
    else
      rho=10*rho
    end if

    if(hnormn.le.eps) then
      write(*,*)'lambda=',lambda
      write(*,*)'norma de h=',hnormn
      write(*,*)'foram feitas',k,'iteracoes'
      write(*,*)'x sem projetar=',(x(j), j=1,dm*np)
      call funf(n,x,f)
      write(*,*)'f sem projetar=',f
      !do i=1,np
        !norma=0.d0
        !do j=1,dm
        !  norma=norma+x((i-1)*dm+j)**2
        !end do
        !norma=sqrt(norma)
        !if(norma.gt.(1.d0))then
          !do j=1,dm
         !   x((i-1)*dm+j)=x((i-1)*dm+j)/norma
        !  end do
       ! end if
      !end do
      !write(*,*)''
      !write(*,*)'x projetado=',(x(j), j=1,dm*np)
     ! call funf(n,x,f)
    !  write(*,*)'f projetado=',f
      return
    else
      k=k+1
	  if(k.ge.itmax)return
      go to 10
    end if
end subroutine

subroutine nfunh(n,x,h,lambda,rho,hnorm)
use globais
integer i, n
double precision x(n), hnorm, h(np), rho, lambda(np)
hnorm=0.d0
do i=1,n
  hnorm=hnorm+h(i)**2
end do
hnorm=sqrt(hnorm)
end subroutine

!subrotina SPG com decréscimo nao monótono
subroutine spg(n,x,xn,g,gn,l,u,alfa,eps,itmax,nafmax,f,fm,Q,h,lambda,rho)
    implicit none
    integer k, n, j, itmax, nafmax, naf, Q
    double precision eps, alfa
    double precision gnor, num, den, p, a, f, fn, t, fmor
    double precision x(n), xn(n), g(n), gn(n), u(n), l(n), fm(Q),h(n),rho, lambda(n)
    
    !inicializando contagem
    k=0
    naf=1

    do j=1,Q
        fm(j)=0
    end do

    call fun(n,x,f,h,lambda,rho)
    call grad(n,x,g,h,lambda,rho)
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
    call fun(n,xn,fn,h,lambda,rho)
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
        !write(*,*)'foram feitas',k,'iteracoes '
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
        call grad(n,x,g,h,lambda,rho)
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
