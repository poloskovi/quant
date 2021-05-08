//! Библиотека работы с комплексными числами

use std::ops::{Add, Mul, Sub, Div};
use std::fmt;
use std::fmt::Display;
//use num_traits::identities::One;

pub trait One<T>{
    fn one()->T;
}

impl One<f64> for f64{
    fn one() -> f64{
        1.0
    }
}
impl One<f32> for f32{
    fn one() -> f32{
        1.0
    }
}
impl One<i32> for i32{
    fn one() -> i32{
        1
    }
}

#[derive(Debug,Copy,Clone)]
pub struct Complex<T>{
    pub re: T,
    pub im: T
}

impl<T> Complex<T>{
    pub fn new(re: T, im: T) -> Complex<T>{
        Complex{
            re,
            im,
        }
    }
    /// квадрат модуля комплексного числа
    pub fn sqrm(self) -> T
    where T: Add<Output=T> + Mul<Output=T> + Copy{
        self.re*self.re + self.im*self.im
    }
    pub fn one() -> Self
    where T: One<T> + Default{
        Self{
            re: T::one(),
            im: T::default(),
        }
    }
    pub fn zero() -> Self
    where T: Default{
        Self{
            re: T::default(),
            im: T::default(),
        }
    }
}

impl<T> Add for Complex<T>
    where T: Add<Output=T> + Copy{
    type Output = Self;
    fn add(self, other: Self) -> Self{
        Self{
            re: self.re + other.re,
            im: self.im + other.im
        }
    }
}

impl<T:Add> Sub for Complex<T>
    where T: Sub<Output=T> + Copy{
    type Output = Self;
    fn sub(self, other: Self) -> Self{
        Self{
            re: self.re - other.re,
            im: self.im - other.im
        }
    }
}

impl<T> Mul for Complex<T>
    where T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Copy{
    type Output = Self;
    fn mul(self, other: Self) -> Self{
        Self{
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re
        }
    }
}

impl<T> Div for Complex<T>
    where T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Default + Copy{
    type Output = Self;
    fn div(self, other: Self) -> Self{
        let x = Self{
            re: other.re,
            im: T::default() - other.im,
        };
        let y = other.re * other.re + other.im * other.im;
        let z = self * x;
        Self{
            re: z.re / y,
            im: z.im / y,
        }
    }
}

impl<T> Default for Complex<T>
    where T: Default{
    fn default() -> Complex<T>{
        Complex{
            re: T::default(),
            im: T::default(),
        }
    }
}

impl<T> fmt::Display for Complex<T>
    where T: Default + PartialEq + PartialOrd + Display + Sub<Output=T> + Copy + One<T>{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let zero = T::default();
        if self.re == zero && self.im == zero{
            write!(f, "{:12}", zero)
        }else if self.im == zero{
           write!(f, "{:12.2}", self.re)
        }else if self.re == zero{
           write!(f, "{:11.2}i", self.im)
        }else if self.im < zero{
            write!(f, "{:5.2}-{:5.2}i", self.re, zero-self.im)
        }else{
            write!(f, "{:5.2}+{:5.2}i", self.re, self.im)
        }
    }
}
