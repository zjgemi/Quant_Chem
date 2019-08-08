


module basic



    contains

    function double_factorial(n)
        ! given a number n,return (2n-1)!1

        integer,intent(in)::n
        real*8::double_factorial
        integer::i

        double_factorial=1.0_8

        do i=1,n
            double_factorial=double_factorial*(2*i-1)
        end do

    end function

    function length(b)
        real*8,intent(in)::b(3)
        real*8::length
        length=b(1)**2+b(2)**2+b(3)**2
    end function


end module

module constant
    use basic
    implicit none
    real*8,parameter::Pi=3.141592653589793_8
    real*8::w(0:5)
    real*8::F(0:5,2,0:6)
end module

module gamma_func
    use basic
    use constant

    contains

    function F_func(m,w1)
    integer,intent(in)::m
    real*8,intent(in)::w1
    real*8::F_func,a,b,ww
    integer::i

    ww=1.0_8
    a=0.0_8
    b=0.0_8

    if (w1<w(m)) then
    do i=0,5
        a=a+F(m,1,i)*ww
        b=b+F(m,2,i)*ww
        ww=ww*w1
    end do

    b=b+F(m,2,6)*ww

    F_func=(a/b)**(m+0.5_8)
    else

    F_func=(double_factorial(m)/(2*w1)**(m+0.5_8))*(Pi/2)**(0.5_8)

    end if


    end function


    subroutine init_param
    implicit none
    w(0)=16.3578_8
    w(1)=17.4646_8
    w(2)=15.2368_8
    w(3)=16.0419_8
    w(4)=16.8955_8
    w(5)=17.7822_8

    F(0,1,0)=1.0_8
    F(0,1,1)=0.213271302431420_8
    F(0,1,2)=0.0629344460255614_8
    F(0,1,3)=0.00769838037756759_8
    F(0,1,4)=0.000758433197127160_8
    F(0,1,5)=0.0000564691197633667_8
    F(0,2,0)=1.0_8
    F(0,2,1)=0.879937801660182_8
    F(0,2,2)=0.338450368470103_8
    F(0,2,3)=0.0738522953299624_8
    F(0,2,4)=0.0101431553402629_8
    F(0,2,5)=0.000955528842975585_8
    F(0,2,6)=0.0000720266520392572_8

    F(1,1,0)=0.4807498567691361_8
    F(1,1,1)=0.0295195994716045_8
    F(1,1,2)=0.0128790985465415_8
    F(1,1,3)=0.000998165499553218_8
    F(1,1,4)=0.0000970927983276419_8
    F(1,1,5)=0.00000493839847029699_8
    F(1,2,0)=1.0_8
    F(1,2,1)=0.461403194579124_8
    F(1,2,2)=0.108494164372449_8
    F(1,2,3)=0.0171462934845042_8
    F(1,2,4)=0.00196918657845508_8
    F(1,2,5)=0.000160138863265254_8
    F(1,2,6)=0.00000857708713007233_8

    F(2,1,0)=0.5253055608807534_8
    F(2,1,1)=-0.00575763488635418_8
    F(2,1,2)=0.00731474973333076_8
    F(2,1,3)=0.000251276149443393_8
    F(2,1,4)=0.0000264336244559094_8
    F(2,1,5)=0.0_8
    F(2,2,0)=1.0_8
    F(2,2,1)=0.274754154712841_8
    F(2,2,2)=0.0425364830353043_8
    F(2,2,3)=0.00493902790955943_8
    F(2,2,4)=0.000437251500927601_8
    F(2,2,5)=0.0000288914662393981_8
    F(2,2,6)=0.0_8
    end subroutine init_param


end module
