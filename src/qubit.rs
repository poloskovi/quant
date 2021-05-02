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
    fn sqrt(self)->T;
}
impl Floating<f64> for f64{
    fn sqrt_2() -> f64{
        std::f64::consts::SQRT_2
    }
    fn epsilon() -> f64{
        0.001
    }
    fn sqrt(self) -> f64{
        self.sqrt()
    }
}
impl Floating<f32> for f32{
    fn sqrt_2() -> f32{
        std::f32::consts::SQRT_2
    }
    fn epsilon() -> f32{
        0.001
    }
    fn sqrt(self) -> f32{
        self.sqrt()
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
    // нормализация, то есть приведение к вектору единичной длины
    pub fn normalize(&mut self)
    where T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Default + Copy + Floating<T>{
        let full_probability = (self.probability(0) + self.probability(1)).sqrt();
        let c = Complex{re:full_probability, im:T::default()};
        self.alfa = self.alfa / c;
        self.beta = self.beta / c;
    }
    // вероятность нахождения в заданном состоянии (0 или 1)
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
    // здесь сделано неправильно. Нужно делать тензорное произведение. Пока такая функция не требуется
//    // преобразовать вектор кубитов в матрицу
//    pub fn vector_to_matrix(vector: &Vec<Qubit<T>>) -> Matrix<Complex<T>>
//    where T: Default + Copy + Clone{
//        let nrow = 2*vector.len();
//        let mut matrix = Matrix::new(nrow,1);
//        for (i, value) in vector.iter().enumerate() {
//            matrix.set(i*2, 0, value.alfa);
//            matrix.set(i*2+1, 0, value.beta);
//        }
//        matrix
//    }
}

impl<T> fmt::Display for Qubit<T>
where T: Default + PartialEq + PartialOrd + Display + Sub<Output=T> + Copy + One<T>{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}) |0> + ({}) |1>", self.alfa, self.beta)
    }
}

// ОПЕРАТОРЫ ЭВОЛЮЦИИ

// Оператор Адамара H
pub fn hadamar<T>(q: Qubit<T>) -> Qubit<T>
where T: Default + Sub<Output=T> + Add<Output=T> + Mul<Output=T> + Div<Output=T> + Copy + One<T>+Floating<T>{

    let matrix_hadamar = matrix_hadamar();
    let m = Matrix::mul(&matrix_hadamar, &q.as_matrix());
    Qubit::from_matrix(m)

}

// Гейт X (оператор квантового NOT)
pub fn X<T>(q: Qubit<T>) -> Qubit<T>
where T: Default + Sub<Output=T> + Add<Output=T> + Mul<Output=T> + Div<Output=T> + Copy + One<T>{

    let matrix_X = matrix_X();
    let m = Matrix::mul(&matrix_X, &q.as_matrix());
    Qubit::from_matrix(m)

}

// Матрица оператора Адамара
pub fn matrix_hadamar<T>() -> Matrix<Complex<T>>
where T: Default + Sub<Output=T> + Add<Output=T> + Mul<Output=T> + Div<Output=T> + Copy + One<T>+Floating<T>{

    let mut matrix = Matrix::new(2,2);

    let one_by_sqrt_2 = Complex{
        re: T::one() / T::sqrt_2(),
        im: T::default(),
    };

    matrix.set(0,0,one_by_sqrt_2);
    matrix.set(0,1,one_by_sqrt_2);
    matrix.set(1,0,one_by_sqrt_2);
    matrix.set(1,1,Complex::zero() - one_by_sqrt_2);

    matrix

}

// Матрица гейта X (квантовый NOT)
pub fn matrix_X<T>() -> Matrix<Complex<T>>
where T: Default + Copy + One<T>{

    let mut matrix = Matrix::new(2,2);

    let one = Complex::one();
    let zero = Complex::zero();

    matrix.set(0,0,zero);
    matrix.set(0,1,one);
    matrix.set(1,0,one);
    matrix.set(1,1,zero);

    matrix

}

// Матрица гейта CNOT
pub fn matrix_CNOT<T>() -> Matrix<Complex<T>>
where T: Default + Copy + One<T>+Floating<T>{

    let mut matrix = Matrix::new(4,4);

    let one = Complex::one();
    let zero = Complex::zero();

    matrix.set(0,0,one);
    matrix.set(0,1,zero);
    matrix.set(0,2,zero);
    matrix.set(0,3,zero);

    matrix.set(1,0,zero);
    matrix.set(1,1,one);
    matrix.set(1,2,zero);
    matrix.set(1,3,zero);

    matrix.set(2,0,zero);
    matrix.set(2,1,zero);
    matrix.set(2,2,zero);
    matrix.set(2,3,one);

    matrix.set(3,0,zero);
    matrix.set(3,1,zero);
    matrix.set(3,2,one);
    matrix.set(3,3,zero);

    matrix

}

// Единичная матрица 2х2
pub fn ed_matrix<T>() -> Matrix<Complex<T>>
where T: Default + Copy + One<T>{

    let mut ed_matrix = Matrix::new(2,2);

    ed_matrix.set(0,0,Complex::<T>::one());
    ed_matrix.set(0,1,Complex::<T>::zero());
    ed_matrix.set(1,0,Complex::<T>::zero());
    ed_matrix.set(1,1,Complex::<T>::one());

    ed_matrix

}

