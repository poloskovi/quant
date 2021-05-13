//! Библиотека работы с кубитами

use std::ops::{Add, Mul, Sub, Div, Neg};
use crate::complex::{Complex, One};
use crate::matrix::Matrix;
use std::fmt;
use std::fmt::Display;

// Типы с плавающей точкой, на которых построены вычисления кубита.
// Определяем для этих типов константы и функции, испольуземые в унитарных операторах эволюции:
pub trait Floating<T>{
    //некое маленькое число, используется как граница погрешности в определении корректности кубита
    fn epsilon()->T;
    // квадратный корень из 2
    fn sqrt(self)->T;
}
impl Floating<f64> for f64{
    fn epsilon() -> f64{
        0.001
    }
    fn sqrt(self) -> f64{
        self.sqrt()
    }
}
impl Floating<f32> for f32{
    fn epsilon() -> f32{
        0.001
    }
    fn sqrt(self) -> f32{
        self.sqrt()
    }
}

#[allow(non_camel_case_types)]
pub trait Sqrt_2<T>{
    fn sqrt_2()->T;
}

impl Sqrt_2<f64> for f64{
    fn sqrt_2() -> f64{
        std::f64::consts::SQRT_2
    }
}
impl Sqrt_2<f32> for f32{
    fn sqrt_2() -> f32{
        std::f32::consts::SQRT_2
    }
}
impl<T> Sqrt_2<Complex<T>> for Complex<T>
where T: Sqrt_2<T> + Default{
    fn sqrt_2()->Complex<T>{
        Complex{
            re: T::sqrt_2(),
            im: T::default(),
        }
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
    // Состояние |0>
    pub fn zero() -> Qubit<T>
    where T: One<T> + Default{
        Qubit{
            alfa: Complex::one(),
            beta: Complex::zero()
        }
    }
    // Состояние |1>
    pub fn one() -> Qubit<T>
    where T: One<T> + Default{
        Qubit{
            alfa: Complex::zero(),
            beta: Complex::one()
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
}

impl<T> fmt::Display for Qubit<T>
where T: Default + PartialEq + PartialOrd + Display + Sub<Output=T> + Copy + One<T>{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}) |0> + ({}) |1>", self.alfa, self.beta)
    }
}

// ОПЕРАТОРЫ ЭВОЛЮЦИИ

// Оператор Адамара H
pub fn hadamard<T>(q: Qubit<T>) -> Qubit<T>
where T: Default + Sub<Output=T> + Add<Output=T> + Mul<Output=T> + Div<Output=T> + Copy + One<T> + Sqrt_2<T>{

    let matrix_hadamard = matrix_hadamard();
    let m = MatrixZeroOne::mul(&matrix_hadamard, &q.as_matrix());
    Qubit::from_matrix(m)

}

// Гейт X (оператор квантового "NOT")
#[allow(non_snake_case)]
#[allow(dead_code)]
pub fn X<T>(q: Qubit<T>) -> Qubit<T>
where T: Default + Sub<Output=T> + Add<Output=T> + Mul<Output=T> + Div<Output=T> + Copy + One<T>{

    let matrix_X = matrix_X();
    let m = MatrixZeroOne::mul(&matrix_X, &q.as_matrix());
    Qubit::from_matrix(m)

}

// Матрица оператора Адамара H
pub fn matrix_hadamard<C>() -> MatrixZeroOne<C>
where C: Div<Output=C> + One<C> + Sqrt_2<C>{

    let mut matrix = Matrix::new(2,2);

    let one = ZeroOne::One;
    let neg_one = ZeroOne::NegOne;

    matrix.set(0,0,one);
    matrix.set(0,1,one);
    matrix.set(1,0,one);
    matrix.set(1,1,neg_one);

    MatrixZeroOne{
        matrix: matrix,
        multiplier: C::one() / C::sqrt_2(),
    }

}

// Матрица n-мерного оператора Адамара Hn
pub fn matrix_hadamard_n<C>(n: u8) -> MatrixZeroOne<C>
where C: Mul<Output=C> + Div<Output=C> + Neg<Output=C> + One<C> + Sqrt_2<C> + Copy{

    let mut zo_matrix = matrix_hadamard();
    for _i in 1..n {
        zo_matrix = MatrixZeroOne::kroneker_product_zo(
            &zo_matrix,
            &matrix_hadamard()
        );
    };

    zo_matrix

}

// Матрица гейта X (квантовый NOT)
#[allow(non_snake_case)]
pub fn matrix_X<C>() -> MatrixZeroOne<C>
where C: One<C>{

    let mut matrix = Matrix::new(2,2);

    let zero = ZeroOne::Zero;
    let one = ZeroOne::One;

    matrix.set(0,0,zero);
    matrix.set(0,1,one);
    matrix.set(1,0,one);
    matrix.set(1,1,zero);

    MatrixZeroOne{
        matrix: matrix,
        multiplier: C::one(),
    }

}

// Матрица гейта CNOT
#[allow(non_snake_case)]
pub fn matrix_CNOT<C>() -> MatrixZeroOne<C>
where C: One<C>{
    matrix_CNOT_n(2, 1, 0)
}

// Матрица гейта CNOT для системы n кубитов
// control - номер контролирующего кубита
// target - номер контролируемого кубита
#[allow(non_snake_case)]
pub fn matrix_CNOT_n<C>(n:u32, control:u32, target: u32) -> MatrixZeroOne<C>
where C: One<C>{

    assert!(n >= 2);
    assert!(n > control);
    assert!(n > target);
    assert!(control != target);

    let size = 2_u32.pow(n) as usize;
    let mut matrix = Matrix::new(size, size);

    let p_control = 2_u32.pow(control) as usize;
    let p_target = 2_u32.pow(target) as usize;

    for column in 0..size{
        let mut row = column;
        let is_control_1 = (column & p_control) != 0;
        if is_control_1{
            // контролирующий кубит = 1
            let is_target_1 = (column & p_target) != 0;
            if is_target_1{
                // контролируемый кубит = 1
                row = column - p_target; // строго не доказал
            }else{
                // контролируемый кубит = 0
                row = column + p_target;
            }
        }
        matrix.set(row, column, ZeroOne::One);
    }

    MatrixZeroOne{
        matrix: matrix,
        multiplier: C::one(),
    }
}

// Единичная матрица 2х2
pub fn ed_matrix<C>() -> MatrixZeroOne<C>
where C: One<C>{

    let mut matrix = Matrix::new(2,2);

    matrix.set(0, 0, ZeroOne::One);
    matrix.set(0, 1, ZeroOne::Zero);
    matrix.set(1, 0, ZeroOne::Zero);
    matrix.set(1, 1, ZeroOne::One);

    MatrixZeroOne{
        matrix: matrix,
        multiplier: C::one(),
    }

}

#[derive(Debug,Copy,Clone,PartialEq)]
pub enum ZeroOne{
    Zero,   // 0
    One,    // 1
    NegOne, // -1
}
impl Default for ZeroOne{
    fn default() -> ZeroOne{
        ZeroOne::Zero
    }
}
impl Neg for ZeroOne{
    type Output = Self;
    fn neg(self) -> Self{
        match self{
            ZeroOne::Zero => ZeroOne::Zero,
            ZeroOne::One => ZeroOne::NegOne,
            ZeroOne::NegOne => ZeroOne::One,
        }
    }
}
impl<T> Mul<T> for ZeroOne
where T: Neg<Output=T> + Default{
    type Output = T;
    fn mul(self, other:T) -> T{
        match self{
            ZeroOne::Zero => T::default(),
            ZeroOne::One => other,
            ZeroOne::NegOne => -other,
        }
    }
}
impl fmt::Display for ZeroOne{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self{
            ZeroOne::Zero => write!(f, " 0"),
            ZeroOne::One => write!(f, " 1"),
            ZeroOne::NegOne => write!(f, "-1")
        }
    }
}

// матрица, состоящая из нулей и единиц, с сомножителем
//#[derive(Debug)]
pub struct MatrixZeroOne<C>{
    // множитель перед матрицей
    pub multiplier: C,
    // сама матрица
    pub matrix: Matrix<ZeroOne>,
}

impl<C> MatrixZeroOne<C>{

    pub fn new(nrow: usize, ncol: usize, multiplier: C) -> MatrixZeroOne<C>{
        MatrixZeroOne{
            matrix: Matrix::new(nrow, ncol),
            multiplier: multiplier
        }
    }

    /// Произведение Кронекера (тензорное умножение)
    /// Возвращает обычную матрицу
    #[allow(dead_code)]
    pub fn kroneker_product(m1: &Matrix<C>, m_zo: &MatrixZeroOne<C>) -> Matrix<C>
    where C: Mul<Output=C> + Neg<Output=C> + Copy + Default{

        let m2 = &m_zo.matrix;
        let mut result = Matrix::new(m1.nrow*m2.nrow, m1.ncol*m2.ncol);

        for i1 in 0..m1.nrow{
            for j1 in 0..m1.ncol{
                for i2 in 0..m2.nrow{
                    for j2 in 0..m2.ncol{
                        let m1_ij = m1.get(i1, j1);// значение ячейки m1
                        let m2_ij = m2.get(i2, j2);// значение ячейки m2
                        let row = i1 * m2.nrow + i2;// строка матрици результата
                        let col = j1 * m2.ncol + j2;// колонка матрицы результата
                        result.set(row,col,
                            m2_ij * m1_ij * m_zo.multiplier);
                    }
                }
            }
        }
        result
    }

    /// Произведение Кронекера (тензорное умножение)
    /// Возвращает матрицу ZeroOne
    pub fn kroneker_product_zo(m_zo1: &MatrixZeroOne<C>, m_zo2: &MatrixZeroOne<C>) -> MatrixZeroOne<C>
    where C: Mul<Output=C> + Neg<Output=C> + Copy{

        let m1 = &m_zo1.matrix;
        let m2 = &m_zo2.matrix;

        let mut result = MatrixZeroOne::new(
            m1.nrow*m2.nrow,
            m1.ncol*m2.ncol,
            m_zo1.multiplier * m_zo2.multiplier
        );

        for i1 in 0..m1.nrow{
            for j1 in 0..m1.ncol{
                for i2 in 0..m2.nrow{
                    for j2 in 0..m2.ncol{
                        let m1_ij = m1.get(i1, j1);// значение ячейки m1
                        let m2_ij = m2.get(i2, j2);// значение ячейки m2
                        let row = i1 * m2.nrow + i2;// строка матрици результата
                        let col = j1 * m2.ncol + j2;// колонка матрицы результата
                        result.matrix.set(row,col,
                            m2_ij * m1_ij);
                    }
                }
            }
        }
        result
    }

    /// Умножение
    pub fn mul(m_zo: &MatrixZeroOne<C>, m2: &Matrix<C>) -> Matrix<C>
    where C: Add<Output=C> + Sub<Output=C> + Mul<Output=C> + Default + Copy{

        let m1 = &m_zo.matrix;
        assert_eq!(m1.ncol, m2.nrow);
        let mut result = Matrix::new(m1.nrow, m2.ncol);

        for i in 0..m1.nrow {
            for j in 0..m2.ncol {
                let mut cij = C::default();
                for r in 0..m1.ncol {
                    let m1_ir = m1.get(i,r);
                    let m2_rj = m2.get(r,j);
                    if m1_ir == ZeroOne::One{
                        cij = cij + m2_rj;
                    }else if m1_ir == ZeroOne::NegOne{
                        cij = cij - m2_rj;
                    }
                };
                result.set(i,j,
                    cij * m_zo.multiplier);
            }
        }
        result
    }

    pub fn to_matrix(&self) -> Matrix<C>
    where C: Neg<Output=C> + Default + Copy{
        let m1 = &self.matrix;
        let mut result = Matrix::new(m1.nrow, m1.ncol);
        for i in 0..m1.nrow {
            for j in 0..m1.ncol {
                let m1_ij = m1.get(i,j);
                result.set(i, j, m1_ij * self.multiplier);
            }
        }
        result
    }
}

impl<C> fmt::Display for MatrixZeroOne<C>
    where C :Display + Copy{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{:.5} * ", self.multiplier)?;
        write!(f, "{}", self.matrix)
    }
}



