mod complex;
use complex::Complex;
mod qubit;
use qubit::Qubit;

extern crate matrix;

fn main() {

    let sqrt2 = std::f64::consts::SQRT_2;

//    let c1:Complex<f64> = Complex::new(1.0/sqrt2, 0.0);
//    let c2:Complex<f64> = Complex::new(-1.0/sqrt2, 0.0);

    let c1:Complex<f64> = Complex::new(1.0/4.0, -sqrt2/4.0);
    let c2:Complex<f64> = Complex::new(-3.0/4.0, 2.0/4.0);

    let q1 = Qubit::new(c1, c2);
    println!("q1: {}", q1);
    if q1.is_correct(){
        println!("Это корректный кубит");
    }else{
        println!("Кубит некорректен!");
    }
    println!("вероятность нахождения в состоянии '0': {:.4}", q1.probability(0));
    println!("вероятность нахождения в состоянии '1': {:.4}", q1.probability(1));
    let qh1 = qubit::hadamar(q1);
    println!("после применения оператора Адамара: {}", qh1);
    println!("вероятность нахождения в состоянии '0': {:.4}", qh1.probability(0));
    println!("вероятность нахождения в состоянии '1': {:.4}", qh1.probability(1));
}
