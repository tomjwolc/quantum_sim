extern crate quantum_sim;
use quantum_sim::*;
use ga_macros::*;

#[test]
fn factor() {
    for i in 0..2 {
        for j in 0..2 {
            let mut v1 = (if i == 0 { ZERO.clone() } else { ONE.clone() }).to_vec();
            let v2 = (if j == 0 { ZERO.clone() } else { ONE.clone() }).to_vec();

            pass_gate(&mut v1, &HADAMARD);

            let p = tensor_product_vector(vec![&v1, &v2]);
            
            println!(
                "{} ⊗ {} = {} -> factored -> {}", 
                print_tensor(&v1),
                print_tensor(&v2),
                print_tensor(&p),
                print_matrix(&tensor_factor(p.clone()))
            );
        }
    }
}

#[test]
fn factor2() {
    let qubit_tensor = vec![
        eq!(1/2),
        eq!(0),
        eq!(0),
        eq!(1/2),
        eq!(1/2),
        eq!(0),
        eq!(0),
        eq!(1/2),
    ];

    println!("{}", print_matrix(&tensor_factor(qubit_tensor)));
}

#[test]
fn tests() {
    println!("{}", print_matrix(&tensor_product_matrix(vec![&CNOT, &IDENTITY])))
}

#[test]
fn gate_test() {
    let mut x = ZERO.clone();
    print!("{} -> ", print_tensor(&x));
    pass_gate(&mut x, &HADAMARD);
    println!("{}", print_tensor(&x));

    let mut x = ONE.clone();
    print!("{} -> ", print_tensor(&x));
    pass_gate(&mut x, &HADAMARD);
    println!("{}", print_tensor(&x));
}