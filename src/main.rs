// Решение упражнений из книги С.С.Сысоев "Введение в квантовые вычисления. Квантовые алгоритмы"

mod complex;
use complex::{Complex,One};
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
        let unknown_operator = unknown_operator_deutsch(variant); // MatrixZeroOne
        let state_3 = MatrixZeroOne::mul(&unknown_operator, &state_2); // Matrix

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

fn unknown_operator_deutsch(variant: u8) -> MatrixZeroOne<Tcomplex>{
    match variant{
        // константа f(x) = 0
        0 => MatrixZeroOne::kroneker_product_zo(
                &qubit::ed_matrix(),
                &qubit::ed_matrix()
            ),
        // константа f(x) = 1
        1 => MatrixZeroOne::kroneker_product_zo(
                &qubit::ed_matrix(),
                &qubit::matrix_X()
            ),
        // сбалансирована f(x) = x
        2 => qubit::matrix_CNOT(),
        // сбалансирована f(x) = -x
        3 => MatrixZeroOne::mul_zo(
                &MatrixZeroOne::mul_zo(
                    &MatrixZeroOne::kroneker_product_zo(
                        &qubit::matrix_X(),
                        &qubit::ed_matrix()
                    ),
                    &qubit::matrix_CNOT()
                ),
                &MatrixZeroOne::kroneker_product_zo(
                    &qubit::matrix_X(),
                    &qubit::ed_matrix()
                )
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
    let digits :u8 = 8;
    let unknown_number = unknown_number();

    // генерируем матрицу начального состояния: [0 0 0 0 .... 0 1],
    // где 0 - на месте разрядов определяемого числа, а 1 - в последнем кубите

    // в результате получается столбец из 2**(digits+1) строк, в котором 2 строка = 1, а остальные = 0
    // поэтому при желании можно эту часть легко оптимизивать, не вычисляя матрицу начального состояния, а генерируя ее.
//    let qubit_zero = Qubit::<Tfloat>::zero();
//    let qubit_one = Qubit::<Tfloat>::one();

//    let mut startstate = qubit_zero.as_matrix();
//    for _i in 1..digits{
//        startstate = Matrix::kroneker_product(
//            &startstate,
//            &qubit_zero.as_matrix()
//        );
//    }
//    startstate = Matrix::kroneker_product(
//        &startstate,
//        &qubit_one.as_matrix()
//    );

    //оптимизировано: матрица начального состояния
    let mut startstate = Matrix::<Tcomplex>::new(
        2_u32.pow(digits as u32 + 1) as usize,
        1
    );
    startstate.set(1,0,Complex::one());

    println!("Это очень мощное колдунство. Поэтому могу считать долго...");


    let hadamard_n = qubit::matrix_hadamard_n(digits + 1);

    let state_2 = MatrixZeroOne::mul(&hadamard_n, &startstate);

    let unknown_operator = operator_bernstein(unknown_number, digits);
    let state_3;
    match unknown_operator{
        None => state_3 = state_2,
        //Some(operator) => state_3 = Matrix::mul(&operator, &state_2),
        Some(operator) => state_3 = MatrixZeroOne::mul(&operator, &state_2),
    }

    let state_final = MatrixZeroOne::mul(&hadamard_n, &state_3);

    let result = find_one_in_state(&state_final);
    match result{
        Some(x) => {

            let mut y = x >> 1;// младший разряд - от контролирующего кубита, отбрасываем

            // биты расположены в обратном порядке
            // (так как клали мы их в орукул тоже в обратном порядке, от младших к старшим)
            let mut z = 0;
            for _ in 0..digits{
                let ostatok = y % 2;
                y = y >> 1;
                z = z << 1;
                z = z + ostatok;
            }
            println!("Загаданное число: {}", z);
        },
        None => println!("Какая-то ошибка"),
    }

}

fn find_one_in_state(m: &Matrix<Tcomplex>) -> Option<usize>{
    for i in 0..m.nrow{
        if m.get(i,0).sqrm()>0.9{
            return Some(i)
        }
    }
    None
}

fn operator_bernstein(unknown_number: u8, digits: u8) -> Option<MatrixZeroOne::<Tcomplex>>{

    let mut operator = None;

    let mut number = unknown_number;

    let mut position = digits + 1;
    while number != 0{

        let ostatok = number % 2;
        number = number >> 1;
        position -= 1;

        if ostatok == 1 {
            let tek_operator = qubit::matrix_CNOT_n(
                (digits+1).into(), position.into() ,0);
            match operator{
                None =>
                    operator = Some(tek_operator),
                Some(old_operator) =>
                    operator = Some(MatrixZeroOne::mul_zo(&old_operator, &tek_operator)),
            }
        }
    };

    operator
}

fn unknown_number() ->u8{
    163
}

fn main() {

    exercise_2_6();
    exercise_2_7();
    exercise_2_10();
    exercise_2_11();

    algorithm_deutsch();

    algoritm_bernstein();

    //test_cnot_n();

}
