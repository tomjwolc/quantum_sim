use ga_macros::*;
use lazy_static::lazy_static;
use colored::*;

lazy_static! {
    pub static ref ZERO: Vec<t!()> = vec![eq!(1), eq!(0)];
    pub static ref ONE: Vec<t!()> = vec![eq!(0), eq!(1)];

    pub static ref IDENTITY: Vec<Vec<t!()>> = vec![
        vec![eq!(1), eq!(0)],
        vec![eq!(0), eq!(1)]
    ];

    pub static ref PAULIX: Vec<Vec<t!()>> = vec![
        vec![eq!(0), eq!(1)],
        vec![eq!(1), eq!(0)]
    ];

    pub static ref PAULIY: Vec<Vec<t!()>> = vec![
        vec![eq!(0), eq!(-i)],
        vec![eq!(i), eq!(0)]
    ];

    pub static ref PAULIZ: Vec<Vec<t!()>> = vec![
        vec![eq!(1), eq!(0)],
        vec![eq!(0), eq!(-1)]
    ];

    pub static ref HADAMARD: Vec<Vec<t!()>> = vec![
        vec![eq!(1/2^0.5), eq!(1/2^0.5)],
        vec![eq!(1/2^0.5), eq!(-1/2^0.5)]
    ];

    pub static ref CNOT: Vec<Vec<t!()>> = vec![
        vec![eq!(1), eq!(0), eq!(0), eq!(0)],
        vec![eq!(0), eq!(1), eq!(0), eq!(0)],
        vec![eq!(0), eq!(0), eq!(0), eq!(1)],
        vec![eq!(0), eq!(0), eq!(1), eq!(0)]
    ];
}

pub fn tensor_product_vector(vectors: Vec<&Vec<t!()>>) -> Vec<t!()> {
    let mut accum = vectors[0].clone();

    for vector in vectors[1..].iter() {
        accum = tensor_product_vector_single(&accum, vector);
    }

    accum
}

fn tensor_product_vector_single(v1: &Vec<t!()>, v2: &Vec<t!()>) -> Vec<t!()> {
    let mut result = Vec::new();

    for n1 in v1.iter() {
        for n2 in v2.iter() {
            result.push(eq!(n1 * n2));
        }
    }

    result
}

pub fn tensor_product_matrix(matrices: Vec<&Vec<Vec<t!()>>>) -> Vec<Vec<t!()>> {
    let mut accum = matrices[0].clone();

    for matrix in matrices[1..].iter() {
        accum = tensor_product_matrix_single(&accum, matrix);
    }

    accum
}

fn tensor_product_matrix_single(m1: &Vec<Vec<t!()>>, m2: &Vec<Vec<t!()>>) -> Vec<Vec<t!()>> {
    let mut result = Vec::new();

    for v1 in m1.iter() {
        for v2 in m2.iter() {
            result.push(tensor_product_vector_single(v1, v2));
        }
    }

    result
}

// Only the magnitudes will be correct! The phases might be wrong
pub fn tensor_factor(mut tensor: Vec<t!()>) -> Vec<Vec<t!()>> {
    let mut base_components = Vec::new();
    let mut component_size = 2;
    let mut component2_size = tensor.len() / 2;

    let t = tensor.clone();

    while component_size < tensor.len() {
        let mut can_be_factored = true;

        let mut i = 0;

        // Checks that all the pairs of pairs are equal
        while i < tensor.len() {
            for j in 0..tensor.len() {
                let ia = i / component2_size; // the part of tensor[i] that was made up from the first component
                let ib = i % component2_size; // the part of tensor[i] that was made up from the second component
                let ja = j / component2_size;
                let jb = j % component2_size;

                if ia == ja || ib == jb { continue };

                can_be_factored = can_be_factored && (
                    eq!(tensor[ia * component2_size + ib] * tensor[ja * component2_size + jb]) == 
                    eq!(tensor[ia * component2_size + jb] * tensor[ja * component2_size + ib])
                );

                if !can_be_factored { break };

            }

            if !can_be_factored { break };

            i += 1;
        }

        if can_be_factored {
            let divisor_index = tensor.iter().position(|z| z != &[0.0, 0.0]).expect(format!("Could not find a non-zero part of the tensor: {}", print_tensor(&tensor)).as_str());
            let divisor = (tensor[divisor_index][0].powf(2.0) + tensor[divisor_index][1].powf(2.0)).powf(0.5);

            tensor[divisor_index][1] *= -1.0; 

            let mut component: Vec<t!()> = (0..component_size).map(|i| {
                let index: usize = divisor_index % component2_size + i * component2_size;

                eq!(tensor[index] / #divisor)
            }).collect();

            tensor = (0..component2_size).map(|i| {
                let index: usize = divisor_index / component2_size * component2_size + i;

                eq!(tensor[index] * tensor[divisor_index] / (#divisor ^ 2))
            }).collect();

            norm(&mut component);
            norm(&mut tensor);

            base_components.push(component);
            component_size = 1;
        }

        component_size *= 2;
        component2_size = tensor.len() / component_size;
    }

    base_components.push(tensor);

    if print_tensor(&tensor_product_vector(base_components.iter().collect())) != print_tensor(&t) {
        panic!(
            "tensor factor failed! result != original: \n    {}\n    or\n    {} \n    != \n    {}", 
            base_components.iter().map(|x| print_tensor(x)).collect::<Vec<String>>().join(" ⊗ "),
            print_tensor(&tensor_product_vector(base_components.iter().collect())), 
            print_tensor(&t)
        )
    }

    base_components
}

pub fn print_tensor(tensor: &Vec<[f64; 2]>) -> String {
    let t = tensor.iter().map(|n| [format!("{:.3}", n[0]), format!("{:.3}", n[1])]).collect::<Vec<[String; 2]>>();

    match *t.iter().map(|n| [correct(n[0].as_str()), correct(n[1].as_str())]).collect::<Vec<[&str; 2]>>() {
        [["1.000", "0.000"], ["0.000", "0.000"]] => {return String::from("      |0⟩     ");},
        [["0.000", "0.000"], ["1.000", "0.000"]] => {return String::from("      |1⟩     ");},
        _ => {}
    }

    format!("[{}]",
    tensor.iter().map(|n| {
        match [correct(format!("{:.3}", n[0]).as_str()), correct(format!("{:.3}", n[1]).as_str())] {
            ["0.000", "0.000"] => String::from("  0  ").bold(),
            ["1.000", "0.000"] => String::from("  1  ").bold(),
            ["0.000", "1.000"] => String::from("  i  ").bold(),
            [a, "0.000"] => format!("{}", a).normal(),
            ["0.000", b] => format!("{}i", b).normal(),
            [a, b] => format!("{} + {}i", a, b).normal()
        }.to_string()
    }).collect::<Vec<String>>().join(", "))
}

pub fn correct(str: &str) -> &str {
    if str == "-0.000" { "0.000" } else { str }
}

pub fn print_matrix(matrix: &Vec<Vec<[f64; 2]>>) -> String {
    format!("[\n    {}\n]",
    matrix.iter().map(|tensor| {
        print_tensor(tensor)
    }).collect::<Vec<String>>().join(", \n    "))
}

pub fn print_matrix_without_breaklines(matrix: &Vec<Vec<[f64; 2]>>) -> String {
    format!("[{}]",
    matrix.iter().map(|tensor| {
        print_tensor(tensor)
    }).collect::<Vec<String>>().join(", "))
}

pub fn pass_gate(tensor: &mut Vec<t!()>, gate: &Vec<Vec<t!()>>) {
    if tensor.len() != gate.len() { panic!("tensor and gate should be the same length") };

    let mut result: Vec<t!()> = vec![eq!(0); tensor.len()];

    for i in 0..tensor.len() {
        for j in 0..tensor.len() {
            result[i] = eq!(result[i] + tensor[j] * gate[i][j]);
        }
    }

    *tensor = result;
}

pub fn norm(tensor: &mut Vec<t!()>) {
    let len = length(tensor);

    for z in tensor.iter_mut() {
        *z = eq!(z / #len);
    }
}

fn length(tensor: &Vec<t!()>) -> f64 {
    let mut len = 0.0;

    for z in tensor.iter() {
        len += z[0].powf(2.0) + z[1].powf(2.0);
    };

    len.powf(0.5)
}

pub fn complement(v: t!()) -> t!() {
    [v[0], -1.0 * v[1]]
}

/* 
        I don't know how to name this function, it takes the index at which your tensor is 
    located in the factored tensor product: t1 ⊗ t2 ⊗ t3, and the index within the product
    state in order to return wether the product state index was calculate with the bottom index
    of the factored tensor you specified
*/  
pub fn is_from_one(tensor_order_index: usize, loc_in_tensor: usize, tensor_len: usize) -> bool {
    loc_in_tensor % (tensor_len / ((2 as f64).powf(tensor_order_index as f64) as usize)) >= tensor_len / ((2 as f64).powf(tensor_order_index as f64 + 1.0) as usize)
}