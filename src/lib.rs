use std::collections::HashMap;

use ga_macros::*;
use rand::Rng;
// use colored::*;

pub mod tensor;

pub use tensor::*;

#[derive(Clone)]
pub enum Instruction {
    Measure(usize, bool),
    MeasureAtAngle(f64, usize, bool),
    MeasureAtSpinVector(Vec<t!()>, usize, bool),
    Circuit(Vec<Instruction>, usize),
    Gate(Vec<Vec<t!()>>, usize),
    Dependent(Box<Instruction>, usize)
}

impl std::fmt::Debug for Instruction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Instruction::Measure(index, _) => write!(f, "{{Measure {} @ |1⟩}}", index),
            Instruction::MeasureAtAngle(angle, index, _) => write!(f, "{{Measure {} @ {:.0}}}", index, angle * 180.0 / std::f64::consts::PI),
            Instruction::MeasureAtSpinVector(spin_vector, index, _) => write!(f, "{{Measure {} @ {}}}", index, print_tensor(spin_vector)),
            Instruction::Gate(gate, index) => write!(f, "{{Gate: {} @ {}}}", print_matrix_without_breaklines(gate), index),
            Instruction::Dependent(circuit, measurement_index) => write!(f, "{{{:?} depending on measurement #{}}}", *circuit, measurement_index),
            Instruction::Circuit(circuit, index) => write!(f, "{{circuit: {:?} @ {}}}", *circuit, index)
        }
    }
}

fn simplify_instructions(instructions: &mut Vec<Instruction>, measurements: usize, vertical_shift: usize) {
    let mut i = 0;
    let mut new_measurements = 0;
    let mut skip_increment = false;

    while i < instructions.len() {
        match instructions.remove(i) {
            Instruction::Circuit(mut circuit, index) => {
                simplify_instructions(&mut circuit, new_measurements + measurements, vertical_shift + index);

                let len = circuit.len();

                while circuit.len() > 0 {
                    instructions.insert(i, circuit.pop().unwrap());
                }

                i += len;
                skip_increment = true;
            },
            Instruction::Measure(index, b) => {
                new_measurements += 1;

                instructions.insert(i, Instruction::Measure(index + vertical_shift, b));
            },
            Instruction::MeasureAtAngle(angle, index, b) => {
                skip_increment = true;

                instructions.insert(
                    i, 
                    Instruction::MeasureAtSpinVector(
                        vec![[(angle / 2.0).sin(), 0.0], [(angle / 2.0).cos(), 0.0]], 
                        index, 
                        b
                    )
                )
            },
            Instruction::MeasureAtSpinVector(spin_vector, index, display) => {
                skip_increment = true;

                let (a, b) = (spin_vector[0], spin_vector[1]);
                let a_ = [a[0], -1.0 * a[1]];
                let b_ = [b[0], -1.0 * b[1]];

                // Rotates the measurement to be vertical
                instructions.insert(
                    i, 
                    Instruction::Gate(
                        vec![
                            vec![  b   ,  eq!(-1 * a)  ],
                            vec![  a_  ,       b_      ]
                        ], 
                        index
                    )
                );

                // Measure
                instructions.insert(
                    i + 1, 
                    Instruction::Measure(index, display)
                ); 

                // Rotates back
                instructions.insert(
                    i + 2, 
                    Instruction::Gate(
                        vec![
                            vec![       b_       ,  a  ],
                            vec![  eq!(-1 * a_)  ,  b  ]
                        ], 
                        index
                    )
                );
            },
            Instruction::Gate(gate, index) => {
                instructions.insert(i, Instruction::Gate(gate, index + vertical_shift));
            },
            Instruction::Dependent(boxed_instruction, measurement_index) => {
                let mut circuit = vec![*boxed_instruction];

                simplify_instructions(&mut circuit, new_measurements + measurements, vertical_shift);

                instructions.insert(
                    i,  
                    Instruction::Dependent(
                        Box::new(Instruction::Circuit(circuit, vertical_shift)), 
                        measurements + measurement_index
                    )
                );
            }
        }

        if !skip_increment { i += 1 };
        skip_increment = false;
    }
}

pub struct Measurement(pub usize, pub bool);

impl std::fmt::Debug for Measurement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}{}", self.0, if self.1 { " " } else { "" }, self.1)
    }
}

pub struct Experiment {
    qbits: Vec<Vec<t!()>>,
    instructions: Vec<Instruction>
}

impl Experiment {
    pub fn new(qbits: Vec<Vec<t!()>>, mut instructions: Vec<Instruction>) -> Self {
        simplify_instructions(&mut instructions, 0, 0);

        // println! ("[\n    {}\n]", instructions.iter().map(|ins| format!("{:?}", ins)).collect::<Vec<String>>().join("\n    "));

        Experiment { qbits, instructions }
    }

    pub fn average_out(&self, sims: usize) -> Vec<(String, Vec<t!()>, Vec<Measurement>, f64)> {
        let mut result_map: HashMap<String, (Vec<t!()>, Vec<Measurement>, f64)> = HashMap::new();

        for _ in 0..sims {
            let run = self.run();

            match result_map.get_mut(&format!("{:?}", run.1)) {
                Some((_, _, num)) => *num += 1.0,
                None => { result_map.insert(format!("{:?}", run.1), (run.0, run.1, 1.0)); }
            };
        }

        result_map.into_iter().map(|(x, (y, z, w))| (x, y, z, w / (sims as f64))).collect()
    }

    pub fn average_out_pretty(&self, sims: usize) -> String {
        self
            .average_out(sims)
            .iter()
            .map(|(_, _, measurements, num)| {
                format!(
                    "Measurements: [ {} ] => {:.2}%", 
                    measurements.iter().map(|m| format!("{:?}", m)).collect::<Vec<String>>().join(", "),
                    100.0 * num
                )
            })
            .collect::<Vec<String>>()
            .join("\n")
    }
    
    pub fn run(&self) -> (Vec<t!()>, Vec<Measurement>) {
        let mut measurements: Vec<(bool, Measurement)> = Vec::new();
        let mut qbits_tensor = tensor_product_vector(self.qbits.iter().collect());
        let mut instructions: Vec<&Instruction> = self.instructions.iter().collect();
        let mut i = 0;
    
        while i < instructions.len() {
            match instructions[i] {
                Instruction::Measure(index, display) => {
                    let measurement = choose_probability(
                        &qbits_tensor.iter().map(|z| z[0].powf(2.0) + z[1].powf(2.0)).collect()
                    ).expect("choose_probability couldn't choose a probability");

                    let result = is_from_one(*index, measurement, qbits_tensor.len());

                    for i in 0..qbits_tensor.len() {
                        if is_from_one(*index, i, qbits_tensor.len()) != result {
                            qbits_tensor[i] = eq!(0);
                        }
                    }

                    measurements.push((*display, Measurement(*index, result)));

                    norm(&mut qbits_tensor);
                },
                Instruction::MeasureAtAngle(_, _, _) => {}, // Should have been removed
                Instruction::MeasureAtSpinVector(_, _, _) => {}, // Should have been removed
                Instruction::Circuit(circuit, _) => { // Should only activate after dependent
                    for instruction in circuit.iter() {
                        instructions.insert(i + 1, instruction);
                    }
                },
                Instruction::Gate(gate, index) => {
                    let mut matrices: Vec<&Vec<Vec<t!()>>> = vec![&*IDENTITY; self.qbits.len()];

                    matrices[*index] = &gate; 

                    let excess_amount = (gate.len() as f64).log2() as usize - 1;

                    for _ in 0..excess_amount { matrices.remove(*index + 1); }

                    pass_gate(&mut qbits_tensor, &tensor_product_matrix(matrices));
                },
                Instruction::Dependent(instruction, measurement_index) => {
                    if measurements
                        .get(*measurement_index)
                        .expect(format!(
                            "measurement_index: {} should not be greater than the size of measurements: {:?}", 
                            *measurement_index, 
                            measurements
                        ).as_str()).1.1
                    {
                        instructions.insert(i + 1, &*instruction)
                    }
                }
            }

            i += 1;
        };
    
        (qbits_tensor, measurements.into_iter().filter(|m| m.0 == true).map(|(_, m)| m).collect())
    }
}

fn choose_probability(probabilities: &Vec<f64>) -> Option<usize> {
    let mut rng = rand::thread_rng();
    let mut r = rng.gen::<f64>();

    for i in 0..probabilities.len() {
        if r < probabilities[i] { return Some(i) };

        r -= probabilities[i];
    }

    None
}