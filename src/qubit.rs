//! Библиотека работы с кубитами

use std::ops::{Add, Mul, Sub, Div};
use crate::complex::{Complex, One};
use crate::matrix::Matrix;
use std::fmt;
use std::fmt::Display;

// Типы с плавающей точкой, на которых построены вычисления кубита.
// Определяем для этих типов константы и функции, испольуземые в унитарных операторах эволюции:
pub trait Floating<T>{
    // квадратный корень из 2
    fn sqrt_2()->T;
    //некое маленькое число, используется как граница погрешности в определении корректности кубита
    fn epsilon()->T;
}
impl Floating<f64> for f64{
    fn sqrt_2() -> f64{
        std::f64::consts::SQRT_2
    }
    fn epsilon() -> f64{
        0.001
    }
}
impl Floating<f32> for f32{
    fn sqrt_2() -> f32{
        std::f32::consts::SQRT_2
    }
    fn epsilon() -> f32{
        0.001
    }
}

// Кубит в состоянии [alfa|0> + beta|1>]
#[derive(Debug,Copy,Clone)]
pub struct Qubit<T>{
    pub alfa: Complex<T>,
    pub beta: Complex<T>,
}

impl<T> Qubit<T>{
    pub fn new(alfa: Complex<T>, beta: Complex<T>) -> Qubit<T>{
        Qubit{
            alfa,
            beta,
        }
    }
    // вероятность нахождения в состоянии в заданном состоянии (0 или 1)
    pub fn probability(self, state:u8) -> T
    where T: Add<Output=T> + Mul<Output=T> + Copy{
        if state == 0{
            self.alfa.sqrm()
        }else{
            self.beta.sqrm()
        }
    }
    // проверка, что это корректный кубит: сумма вероятностей в состоянии 0 и 1 равна единице
    pub fn is_correct(self) -> bool
    where T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Copy + PartialOrd + Default + One<T> + Floating<T>{
        // полная вероятность
        let full_probability = self.probability(0) + self.probability(1);
        // разность между полной вероятностью и единицей
        let delta = full_probability - T::one();
        let zero = T::default();//нуль
        if delta < zero{
            (zero - delta) < T::epsilon()
        }else{
            (delta - zero) < T::epsilon()
        }
    }
    pub fn as_matrix(self) -> Matrix<Complex<T>>
    where T: Default + Copy + Clone{
        let nrow = 2;
        let mut res = Matrix::new(nrow,1);
        res.set(0,0,self.alfa);
        res.set(1,0,self.beta);
        res
    }
    pub fn from_matrix(m: Matrix<Complex<T>>) -> Self
    where T: Copy{
        Self::new(m.get(0,0), m.get(0,1))
    }
}

impl<T> fmt::Display for Qubit<T>
where T: Default + PartialEq + PartialOrd + Display + Sub<Output=T> + Copy + One<T>{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}|0> + {}|1>", self.alfa, self.beta)
    }
}

// ОПЕРАТОРЫ ЭВОЛЮЦИИ

// Оператор Адамара H
pub fn hadamar<T>(q: Qubit<T>) -> Qubit<T>
where T: Default + Sub<Output=T> + Add<Output=T> + Mul<Output=T> + Div<Output=T> + Copy + One<T>+Floating<T>{

    let mut hadamar = Matrix::new(2,2);

    let one_by_sqrt_2 = Complex{
        re: T::one() / T::sqrt_2(),
        im: T::default(),
    };

    hadamar.set(0,0,one_by_sqrt_2);
    hadamar.set(0,1,one_by_sqrt_2);
    hadamar.set(1,0,one_by_sqrt_2);
    hadamar.set(1,1,Complex::zero() - one_by_sqrt_2);

    let m = Matrix::mul(&hadamar, &q.as_matrix());
    Qubit::from_matrix(m)

}

// Гейт X (оператор квантового NOT)
pub fn q_not<T>(q: Qubit<T>) -> Qubit<T>
where T: Default + Sub<Output=T> + Add<Output=T> + Mul<Output=T> + Div<Output=T> + Copy + One<T>{

    let mut q_not = Matrix::new(2,2);

    let one = Complex::one();
    let zero = Complex::zero();

    q_not.set(0,0,zero);
    q_not.set(0,1,one);
    q_not.set(1,0,one);
    q_not.set(1,1,zero);

    let m = Matrix::mul(&q_not, &q.as_matrix());
    Qubit::from_matrix(m)

}
