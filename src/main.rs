// Решение упражнений из книги С.С.Сысоев "Введение в квантовые вычисления. Квантовые алгоритмы"

mod complex;
use complex::Complex;
mod qubit;
use qubit::Qubit;

extern crate matrix;
use matrix::Matrix;

fn print_correct(qb: Qubit<f64>){
    if qb.is_correct(){
        println!("Это корректный кубит {}", qb);
    }else{
        println!("Кубит некорректен! {}", qb);
    }
}

// Состояние |φ> = (1-√2i)/4 |0> - (3-2i)/4 |1>
// измеряют в стандартном базисе |0>, |1>.
// Какова вероятность получить результат |1> ?
fn exercise_2_6(){

    println!("");
    println!("Упражнение 2.6");

    let sqrt2 = std::f64::consts::SQRT_2;

    let c1:Complex<f64> = Complex::new(1.0/4.0, -sqrt2/4.0);
    let c2:Complex<f64> = Complex::new(-3.0/4.0, 2.0/4.0);

    let q1 = Qubit::new(c1, c2);
    print_correct(q1);
    println!("вероятность нахождения в состоянии '1': {:.4}", q1.probability(1));
}

// Состояние |φ> = (1-√2i)/4 |0> - (3-2i)/4 |1>
// измеряют в базисе Адамара |+>, |->.
// Какова вероятность получить результат |+> ?
fn exercise_2_7(){

    println!("");
    println!("Упражнение 2.7");

    let sqrt2 = std::f64::consts::SQRT_2;

    let c1:Complex<f64> = Complex::new(1.0/4.0, -sqrt2/4.0);
    let c2:Complex<f64> = Complex::new(-3.0/4.0, 2.0/4.0);

    let q1 = Qubit::new(c1, c2);
    print_correct(q1);
    let qh1 = qubit::hadamar(q1);
    println!("после применения оператора Адамара: {}", qh1);
    println!("вероятность нахождения в состоянии '+': {:.4}", qh1.probability(0));

}

// Первый кубит состояния |φ> = 1/2 |00> - 1/2 |01> + 1/2 |10> + 1/2 |11>
// измерили в стандартном базисе и получили вектор |0>
// Какова вероятность при измерении второго кубита в базисе Адамара получить вектор |+> ?
fn exercise_2_10(){

    println!("");
    println!("Упражнение 2.10");

    // Решение: При измерении первого кубита исходное состояние переходит в [1/2 |00> - 1/2 |01>]
    // Домножаем каждую вероятность на 1/√2, чтобы получить корректное состояние (с полной вероятностью = 1)

    let sqrt2 = std::f64::consts::SQRT_2;

    let c1:Complex<f64> = Complex::new(1.0/sqrt2, 0.0);
    let c2:Complex<f64> = Complex::new(-1.0/sqrt2, 0.0);
    let q1 = Qubit::new(c1, c2);
    print_correct(q1);
    let qh1 = qubit::hadamar(q1);
    println!("после применения оператора Адамара: {}", qh1);
    println!("вероятность нахождения в состоянии '+': {:.4}", qh1.probability(0));

}

// Первый кубит состояния |φ> = (1+√3i)/4 |00> + (1-√3i)/4 |01> + (1-√3i)/4 |10> + (1+√3i)/4 |11>
// измерили в базисе Адамара и получили вектор |->
// Какова вероятность при измерении второго кубита в стандартном базисе, получить вектор |0> ?
fn exercise_2_11(){

    println!("");
    println!("Упражнение 2.11");

    let sqrt3 = 3_f64.sqrt();
    let c1:Complex<f64> = Complex::new(1.0/4.0, sqrt3/4.0);
    let c2:Complex<f64> = Complex::new(1.0/4.0, -sqrt3/4.0);

    let mut state = Matrix::new(4,1);
    state.set(0,0,c1);
    state.set(0,1,c2);
    state.set(0,2,c2);
    state.set(0,3,c1);

    println!("Исходное состояние:");
    println!("{}", state);

    //к первому кубиту применяем оператор Адамара, ко второму - единичный оператор
    let operator = Matrix::kroneker_product(
        &qubit::matrix_hadamar(),
        &qubit::ed_matrix()
    );
    println!("Оператор H⊗I:");
    println!("{}",operator);

    let state_2 = Matrix::mul(&operator, &state);
    println!("Cостояние системы после воздействия оператором H⊗I:");
    println!("{}", state_2);

    // Измерение показало, что первый вектор находится в состоянии |->. Следовательно, в состоянии
    // а |+0> + б |+1> + г |-0> + д |-1> реализованы два последних слагаемых г |-0> + д |-1>:
    let mut q2 = Qubit::new(
        state_2.get(2,0),
        state_2.get(3,0)
    );
    println!("вектор 2: {}", q2);
    q2.normalize();
    println!("после нормализации : {}", q2);
    println!("вероятность нахождения в состоянии '0': {:.4}", q2.probability(0));

}

// Задача Дойча
// Задача Дойча — Йожи заключается в определении, является ли функция нескольких двоичных переменных
// f ( x 1 , x 2 , … , x n ) постоянной (принимает либо значение 0, либо 1 при любых аргументах)
// или сбалансированной (для половины области определения принимает значение 0, для другой половины 1)
fn algorithm_deutsch(){

    println!("");
    println!("Задача Дойча");

    let qx = Qubit::<f64>::new(Complex::one(), Complex::zero());
    let qy = Qubit::<f64>::new(Complex::zero(), Complex::one());

    let state = Matrix::kroneker_product(
        &qx.as_matrix(),
        &qy.as_matrix()
    );
    println!("Начальное состояние:");
    println!("{}", state);

    for variant in 0..4{

        println!("вариант Оракула: {}", variant);
        //к обоим кубитам применяем оператор Адамара
        let operator = Matrix::kroneker_product(
            &qubit::matrix_hadamar(),
            &qubit::matrix_hadamar()
        );

        let state_2 = Matrix::mul(&operator, &state);

        // состояние после воздействия неизвестного оператора
        let unknown_operator = unknown_operator(variant);
        let state_3 = Matrix::mul(&unknown_operator, &state_2);

        // к кубиту x применяем оператор Адамара
        let operator_hadamar_ed = Matrix::kroneker_product(
            &qubit::matrix_hadamar(),
            &qubit::ed_matrix()
        );
        let state_4 = Matrix::mul(&operator_hadamar_ed, &state_3);
        println!("Cостояние системы в конце:");
        println!("{}", state_4);

    }

}

fn unknown_operator(variant: u8) -> Matrix<Complex<f64>>{
    match variant{
        // константа f(x) = 0
        0 => Matrix::kroneker_product(
                &qubit::ed_matrix(),
                &qubit::ed_matrix()
            ),
        // константа f(x) = 1
        1 => Matrix::kroneker_product(
                &qubit::ed_matrix(),
                &qubit::matrix_X()
            ),
        // сбалансирована f(x) = x
        2 => qubit::matrix_CNOT(),
        // сбалансирована f(x) = -x
        3 => Matrix::mul(
                &Matrix::mul(
                    &Matrix::kroneker_product(
                        &qubit::matrix_X(),
                        &qubit::ed_matrix()
                    ),
                    &qubit::matrix_CNOT()
                ),
                &Matrix::kroneker_product(
                    &qubit::matrix_X(),
                    &qubit::ed_matrix()
                )
            ),
        _=> panic!("")
    }
}

fn main() {

    exercise_2_6();
    exercise_2_7();
    exercise_2_10();
    exercise_2_11();

    algorithm_deutsch();

}
