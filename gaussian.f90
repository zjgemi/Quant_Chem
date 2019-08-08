module gauss_integration

    use constant
    use basic
    use gamma_func
    contains


    function overlap(a,RA,Ia,b,RB,Ib)
        implicit none
        ! calculate the overlap function for gauss function, center at Ra and rb
        ! a for exp{-a*(r-Ra)^2}
        ! Ia means (0,0,0) for s type function,(1,0,0) for px-type
        ! the prefactor always 1
        real*8,intent(in)::a,b,Ra(3),Rb(3)
        integer,intent(in)::Ia(3),Ib(3)
        integer::i
        real*8:: Rp(3),I1(3),overlap,K

        K=exp(-(a*b)/(a+b)*length(Ra-Rb))

        Rp=(b/(a+b))*Rb+(a/(a+b))*Ra

        do i=1,3


            if ((Ia(i).eq.0).and.(Ib(i).eq.0)) then

             I1(i)= sqrt(Pi/(a+b))


            end if

            if ((Ia(i).eq.1).and.(Ib(i).eq.1)) then

             I1(i)=sqrt(Pi)/(2*(a+b)**(1.5_8))+(Rp(i)-Ra(i))*(Rp(i)-Rb(i))*sqrt(Pi/(a+b))

            end if

            if ((Ia(i).eq.1).and.(Ib(i).eq.0)) then

            I1(i)=(Rp(i)-Ra(i))*sqrt(Pi/(a+b))

            end if


            if ((Ia(i).eq.0).and.(Ib(i).eq.1)) then

            I1(i)=(Rp(i)-Rb(i))*sqrt(Pi/(a+b))

            end if

        end do

        overlap=K*I1(1)*I1(2)*I1(3)

    end function

    function overlap_v2(a,RA,Ia,b,RB,Ib)
        implicit none
        ! calculate the overlap function for gauss function, center at Ra and rb,at one direction
        ! a for exp{-a*(r-Ra)^2}
        ! Ia means (0,0,0) for s type function,(1,0,0) for px-type
        ! the prefactor always 1
        real*8,intent(in)::a,b,Ra,Rb
        integer,intent(in)::Ia,Ib
        real*8:: Rp,I,overlap_v2,K

        K=exp(-(a*b)/(a+b)*(Ra-Rb)**2)

        Rp=(b/(a+b))*Rb+(a/(a+b))*Ra




            if ((Ia.eq.0).and.(Ib.eq.0)) then

             I= sqrt(Pi/(a+b))


            end if

            if ((Ia.eq.1).and.(Ib.eq.1)) then

             I=sqrt(Pi)/(2*(a+b)**(1.5_8))+(Rp-Ra)*(Rp-Rb)*sqrt(Pi/(a+b))

            end if

            if ((Ia.eq.1).and.(Ib.eq.0)) then

             I=(Rp-Ra)*sqrt(Pi/(a+b))

            end if


            if ((Ia.eq.0).and.(Ib.eq.1)) then

             I=(Rp-Rb)*sqrt(Pi/(a+b))

            end if



        overlap_v2=K*I

    end function

 function kinetic(a,RA,Ia,b,RB,Ib)
    implicit none
        ! calculate the kinetic function for gauss function, center at Ra and rb
        ! a for exp{-a*(r-Ra)^2}
        ! Ia means (0,0,0) for s type function,(1,0,0) for px-type
        ! the prefactor always 1
        real*8,intent(in)::a,b,Ra(3),Rb(3)
        integer,intent(in)::Ia(3),Ib(3)
        integer::i
        real*8:: Rp(3),I1(3),kinetic,K,S(3),aa

        K=exp(-(a*b)/(a+b)*length(Ra-Rb))

        Rp=(b/(a+b))*Rb+(a/(a+b))*Ra

        S(1)=overlap_v2(a,Ra(2),Ia(2),b,Rb(2),Ib(2))*overlap_v2(a,Ra(3),Ia(3),b,Rb(3),Ib(3))

        S(2)=overlap_v2(a,Ra(3),Ia(3),b,Rb(3),Ib(3))*overlap_v2(a,Ra(1),Ia(1),b,Rb(1),Ib(1))

        S(3)=overlap_v2(a,Ra(2),Ia(2),b,Rb(2),Ib(2))*overlap_v2(a,Ra(1),Ia(1),b,Rb(1),Ib(1))


        do i=1,3

            if ((Ia(i).eq.0).and.(Ib(i).eq.0)) then


            I1(i)=-2*b*overlap_v2(a,ra(i),0,b,rb(i),0)+4*b**2*K*(sqrt(Pi)/(2*(a+b)**(1.5_8))+(Rp(i)-Rb(i))**2*sqrt(Pi/(a+b)))

            end if

            if ((Ia(i).eq.1).and.(Ib(i).eq.1)) then

            aa=3*sqrt(Pi)/(4*(a+b)**(2.5_8))+3*(rp(i)-rb(i))*(2*rp(i)-rb(i)-ra(i))*sqrt(Pi)/(2*(a+b)**(1.5_8))
            aa=aa+(rp(i)-ra(i))*(rp(i)-rb(i))**3*sqrt(Pi/(a+b))

            I1(i)=-6*b*overlap_v2(a,ra(i),1,b,rb(i),1)+4*b**2*K*aa

            end if

            if ((Ia(i).eq.1).and.(Ib(i).eq.0)) then


            aa=(3*rp(i)-2*rb(i)-ra(i))*sqrt(Pi)/(2*(a+b)**(1.5_8))+(rp(i)-ra(i))*(rp(i)-rb(i))**2*sqrt(Pi/(a+b))

            I1(i)=-2*b*overlap_v2(a,ra(i),1,b,rb(i),0)+4*b**2*K*aa

            end if


            if ((Ia(i).eq.0).and.(Ib(i).eq.1)) then


            aa=3*(rp(i)-rb(i))*sqrt(Pi)/(2*(a+b)**(1.5_8))+(rp(i)-rb(i))**3*sqrt(Pi/(a+b))


            I1(i)=-6*b*overlap_v2(a,ra(i),0,b,rb(i),1)+4*b**2*K*aa

            end if



        end do


        kinetic=-0.5*(I1(1)*S(1)+I1(2)*S(2)+I1(3)*S(3))

 end function

  function e_attract(a,RA,Ia,b,RB,Ib,Rc)
        implicit none
        ! calculate the eletron attract function for gauss function, center at Ra and rb,towards rc
        ! a for exp{-a*(r-Ra)^2}
        ! Ia means (0,0,0) for s type function,(1,0,0) for px-type
        ! the prefactor always 1

        real*8,intent(in)::a,b,Ra(3),Rb(3),Rc(3)
        integer,intent(in)::Ia(3),Ib(3)
        integer::i,sum1,sum2
        real*8:: Rp(3),e_attract,K,a1,a2


        K=exp(-(a*b)/(a+b)*length(Ra-Rb))

        Rp=(b/(a+b))*Rb+(a/(a+b))*Ra

        sum1=ia(1)+ia(2)+ia(3)
        sum2=ib(1)+ib(2)+ib(3)

        if ((sum1.eq.0).and.(sum2.eq.0)) then

        e_attract=K*2*Pi/(a+b)*F_func(0,(a+b)*length(Rp-Rc))

        end if



        if ((sum1.eq.1).and.(sum2.eq.0)) then

        do i=1,3

            if(ia(i).eq.1) then
                sum1=ia(i)
            end if

        end do

        a1=2*a*(rc(sum1)-rp(sum1))*F_func(1,(a+b)*length(rp-rc))
        a2=2*a*b/(a+b)*(ra(sum1)-rb(sum1))*F_func(0,(a+b)*length(rp-rc))
        e_attract=Pi/(a*(a+b))*(a1-a2)*K



        end if



        if ((sum1.eq.0).and.(sum2.eq.1)) then
        do i=1,3

            if(ib(i).eq.1) then
                sum1=ib(i)
            end if

        end do

        a1=2*b*(rc(sum1)-rp(sum1))*F_func(1,(a+b)*length(rp-rc))
        a2=2*a*b/(a+b)*(-ra(sum1)+rb(sum1))*F_func(0,(a+b)*length(rp-rc))
        e_attract=Pi/(a*(a+b))*(a1-a2)*K
        end if




        if ((sum1.eq.1).and.(sum2.eq.1)) then
        e_attract=K*2*Pi/(a+b)*F_func(0,(a+b)*length(Rp-Rc))
        end if
  end function

end module
