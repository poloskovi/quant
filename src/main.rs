// Решение упражнений из книги С.С.Сысоев "Введение в квантовые вычисления. Квантовые алгоритмы"

mod complex;
use complex::Complex;
mod qubit;
use qubit::{Qubit, MatrixZeroOne};

type Tfloat = f64;
type Tcomplex = Complex<Tfloat>;

extern crate matrix;
use matrix::Matrix;

fn print_correct(qb: Qubit<Tfloat>){
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
    let qh1 = qubit::hadamard(q1);
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
    let qh1 = qubit::hadamard(q1);
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
    let c1:Tcomplex = Complex::new(1.0/4.0, sqrt3/4.0);
    let c2:Tcomplex = Complex::new(1.0/4.0, -sqrt3/4.0);

    let mut state = Matrix::new(4,1);
    state.set(0,0,c1);
    state.set(0,1,c2);
    state.set(0,2,c2);
    state.set(0,3,c1);

    println!("Исходное состояние:");
    println!("{}", state);

    //к первому кубиту применяем оператор Адамара, ко второму - единичный оператор
    let operator = MatrixZeroOne::kroneker_product_zo(
        &qubit::matrix_hadamard::<Tcomplex>(),
        &qubit::ed_matrix::<Tcomplex>()
    );
    println!("Оператор H⊗I:");
    println!("{}",operator);

    let state_2 = MatrixZeroOne::mul(&operator, &state);
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
// Задача Дойча заключается в определении, является ли функция f (x) постоянной
// (принимает либо значение 0, либо 1 при любых аргументах)
// или сбалансированной (для половины области определения принимает значение 0, для другой половины 1)
fn algorithm_deutsch(){

    println!("");
    println!("Задача Дойча");

    let qx = Qubit::<Tfloat>::zero();
    let qy = Qubit::<Tfloat>::one();

    let state = Matrix::kroneker_product(
        &qx.as_matrix(),
        &qy.as_matrix()
    );

    for variant in 0..4{

        println!("");
        println!("вариант Оракула: {}", variant);

        //к обоим кубитам применяем оператор Адамара
        let operator = qubit::matrix_hadamard_n(2);
        let state_2 = MatrixZeroOne::mul(&operator, &state); // Matrix

        // к обоим кубитам применяем неизвестный оператор
        let unknown_operator = unknown_operator(variant); // MatrixZeroOne
        let state_3 = Matrix::mul(&unknown_operator, &state_2); // Matrix

        // к кубиту x применяем оператор Адамара
        let operator_hadamard_ed = MatrixZeroOne::kroneker_product_zo(
            &qubit::matrix_hadamard::<Tcomplex>(),
            &qubit::ed_matrix()
        );
        let state_4 = MatrixZeroOne::mul(&operator_hadamard_ed, &state_3);

        // измерение кубита x:
        let probability_const = state_4.get(0,0).sqrm() + state_4.get(1,0).sqrm();
        let probability_balans = state_4.get(2,0).sqrm() + state_4.get(3,0).sqrm();

        if probability_const>probability_balans{
            println!("Оракул - константа с вероятностью {:.2}%", probability_const*100.0);
        }else{
            println!("Оракул - сбалансирован с вероятностью {:.2}%", probability_balans*100.0);
        }
    }

}

// поменять результат на Matrix, так как в последнем варианте при перемножении матриц
// законно могут появиться "2" в результате.
// Сделать функцию пребразования MatrixZeroOne в Matrix.
fn unknown_operator(variant: u8) -> Matrix<Tcomplex>{
    match variant{
        // константа f(x) = 0
        0 => MatrixZeroOne::kroneker_product_zo(
                &qubit::ed_matrix(),
                &qubit::ed_matrix()
            ).to_matrix(),
        // константа f(x) = 1
        1 => MatrixZeroOne::kroneker_product_zo(
                &qubit::ed_matrix(),
                &qubit::matrix_X()
            ).to_matrix(),
        // сбалансирована f(x) = x
        2 => qubit::matrix_CNOT().to_matrix(),
        // сбалансирована f(x) = -x
        3 => Matrix::mul(
                &MatrixZeroOne::mul(
                    &MatrixZeroOne::kroneker_product_zo(
                        &qubit::matrix_X(),
                        &qubit::ed_matrix()
                    ),
                    &qubit::matrix_CNOT().to_matrix()
                ),
                &MatrixZeroOne::kroneker_product_zo(
                    &qubit::matrix_X(),
                    &qubit::ed_matrix()
                ).to_matrix()
            ),
        _=> unreachable!()
    }
}

// Алгоритм Бернштейна - Вазирани.
// квантовый алгоритм, решающий задачу нахождения n-битного числа, скрытого в черном ящике
fn algoritm_bernstein(){

    println!("");
    println!("Алгоритм Бернштейна - Вазирани");

    // количество бит в числе
    let digits = 3;

    // генерируем матрицу начального состояния: [0 0 0 0 .... 0 1],
    // где 0 - на месте разрядов определяемого числа, а 1 - в последнем столбце

    let qubit_zero = Qubit::<Tfloat>::zero();
    let qubit_one = Qubit::<Tfloat>::one();

    let mut startstate = qubit_zero.as_matrix();
    for _i in 1..digits{
        startstate = Matrix::kroneker_product(
            &startstate,
            &qubit_zero.as_matrix()
        );
    }
    startstate = Matrix::kroneker_product(
        &startstate,
        &qubit_one.as_matrix()
    );

    println!("Начальное состояние:");
    println!("{}", startstate);

    // в результате получается столбец из 2**(digits+1) строк, в котором 2 строка = 1, а остальные = 0
    // поэтому при желании можно эту часть легко оптимизивать, не вычисляя матрицу начального состояния, а генерируя ее.

    let hadamard_n = qubit::matrix_hadamard_n(digits + 1);

    println!("Матрица оператора Адамара Hn:");
    println!("{}", hadamard_n);

    let state_2 = MatrixZeroOne::mul(&hadamard_n, &startstate);

    println!("Состояние после воздействия оператора Адамара:");
    println!("{}", state_2);

    // не доделано

}

fn test_cnot_n(){

    // проверка работы оператора CNOT для n кубитов

    let operator: MatrixZeroOne::<Tcomplex> = qubit::matrix_CNOT_n(2,1,0);
    println!("{}", operator);

    let operator: MatrixZeroOne::<Tcomplex> = qubit::matrix_CNOT();
    println!("{}", operator);
    println!("I x CNOT = {}", MatrixZeroOne::kroneker_product_zo(&qubit::ed_matrix(), &operator));

    let operator: MatrixZeroOne::<Tcomplex> = qubit::matrix_CNOT_n(3,1,0);
    println!("1 - 0 (I x CNOT): {}", operator);

    let operator: MatrixZeroOne::<Tcomplex> = qubit::matrix_CNOT_n(3,2,0);
    println!("2 - 0: {}", operator);

    let operator: MatrixZeroOne::<Tcomplex> = qubit::matrix_CNOT_n(3,2,1);
    println!("2 - 1 (CNOT x I): {}", operator);

}

fn main() {

    exercise_2_6();
    exercise_2_7();
    exercise_2_10();
    exercise_2_11();

    algorithm_deutsch();

    algoritm_bernstein();

    test_cnot_n();

}
