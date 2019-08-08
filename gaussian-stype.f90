module gauss_integration_s

    use gauss_integration
    use constant
    use basic
    use gamma_func
    contains


    function overlap_s(a,RA,b,RB)
        implicit none
        ! calculate the overlap function for gauss function, center at Ra and rb
        ! a for exp{-a*(r-Ra)^2}
        ! Ia means (0,0,0) for s type function,(1,0,0) for px-type
        ! the prefactor always 1
        real*8,intent(in)::a,b,Ra(3),Rb(3)
        integer::i
        real*8:: Rp(3),I1(3),overlap_s,K

        K=exp(-(a*b)/(a+b)*length(Ra-Rb))

        !Rp=(b/(a+b))*Rb+(a/(a+b))*Ra

        do i=1,3
             I1(i)= sqrt(Pi/(a+b))
        end do

        overlap_s=K*I1(1)*I1(2)*I1(3)

    end function

    function overlap_v2_s(a,RA,b,RB)
        implicit none
        ! calculate the overlap function for gauss function, center at Ra and rb,at one direction
        ! a for exp{-a*(r-Ra)^2}
        ! Ia means (0,0,0) for s type function,(1,0,0) for px-type
        ! the prefactor always 1
        real*8,intent(in)::a,b,Ra,Rb
        real*8:: Rp,I,overlap_v2_s,K

        K=exp(-(a*b)/(a+b)*(Ra-Rb)**2)


        !Rp=(b/(a+b))*Rb+(a/(a+b))*Ra

        I= sqrt(Pi/(a+b))

        overlap_v2_s=K*I

    end function

 function kinetic_s(a,RA,b,RB)
    implicit none
        ! calculate the kinetic function for gauss function, center at Ra and rb
        ! a for exp{-a*(r-Ra)^2}
        ! Ia means (0,0,0) for s type function,(1,0,0) for px-type
        ! the prefactor always 1
        real*8,intent(in)::a,b,Ra(3),Rb(3)
        integer::i
        real*8:: Rp(3),I1(3),kinetic_s,K

        K=exp(-(a*b)/(a+b)*length(Ra-Rb))

        Rp=(b/(a+b))*Rb+(a/(a+b))*Ra

        kinetic_s=a*b/(a+b)*(3-2*a*b/(a+b)*length(ra-rb))*K*(Pi/(a+b))**(1.5_8)

 end function

  function e_attract_s(a,RA,b,RB,Rc)
        implicit none
        ! calculate the eletron attract function for gauss function, center at Ra and rb,towards rc
        ! a for exp{-a*(r-Ra)^2}
        ! Ia means (0,0,0) for s type function,(1,0,0) for px-type
        ! the prefactor always 1

        real*8,intent(in)::a,b,Ra(3),Rb(3),Rc(3)
        integer::i
        real*8:: Rp(3),e_attract_s,K


        K=exp(-(a*b)/(a+b)*length(Ra-Rb))

        Rp=(b/(a+b))*Rb+(a/(a+b))*Ra

        e_attract_s=K*2*Pi/(a+b)*F_func(0,(a+b)*length(Rp-Rc))

  end function

  function e_repulsion_s(a,RA,b,RB,c,Rc,d,rd)
        implicit none
        ! calculate the eletron repulsion function for gauss function, center at Ra and rb,towards rc
        ! a for exp{-a*(r-Ra)^2}
        ! Ia means (0,0,0) for s type function,(1,0,0) for px-type
        ! the prefactor always 1

        real*8,intent(in)::a,b,c,d,Ra(3),Rb(3),Rc(3),Rd(3)
        integer::i
        real*8:: Rp1(3),e_repulsion_s,K1,Rp2(3),K2,K,aa


        K1=exp(-(a*b)/(a+b)*length(Ra-Rb))
        K2=exp(-(c*d)/(c+d)*length(Rc-Rd))
        Rp1=(b/(a+b))*Rb+(a/(a+b))*Ra
        Rp2=(c/(c+d))*Rc+(d/(c+d))*Rd
        K=2*Pi**(2.5_8)/((a+b)*(c+d)*sqrt(a+b+c+d))
        aa=length(rp1-rp2)*(a+b)*(c+d)/(a+b+c+d)
        e_repulsion_s=K*K1*K2*F_func(0,aa)

  end function

end module
