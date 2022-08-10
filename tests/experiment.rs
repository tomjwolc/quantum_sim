extern crate quantum_sim;

use quantum_sim::*;
use lazy_static::lazy_static;
use rand::Rng;

lazy_static! {
    static ref ENTANGLE: Vec<Instruction> = vec![
        Instruction::Gate(HADAMARD.clone(), 0),
        Instruction::Gate(CNOT.clone(), 0),
    ];

    static ref TELEPORT: Vec<Instruction> = vec![
        Instruction::Circuit(ENTANGLE.clone(), 1),
        Instruction::Gate(CNOT.clone(), 0),
        Instruction::Gate(HADAMARD.clone(), 0),
        Instruction::Measure(1, false),
        Instruction::Measure(0, false),
        Instruction::Dependent(Box::new(Instruction::Gate(PAULIX.clone(), 2)), 0),
        Instruction::Dependent(Box::new(Instruction::Gate(PAULIZ.clone(), 2)), 1),
    ];
}


#[cfg(test)]
mod experiment {
    use super::*;

    #[test]
    fn entanglement() {
        let experiment = Experiment::new(
            vec![
                ZERO.clone(),
                ZERO.clone()
            ],
            vec![
                Instruction::Gate(HADAMARD.clone(), 0),
                Instruction::Gate(CNOT.clone(), 0),
                Instruction::Measure(0, true),
                Instruction::Measure(1, true)
            ]
        );

        println!("{}", &experiment.average_out_pretty(10000));
    }

    #[test]
    fn independent_entanglement() {
        let experiment = Experiment::new(
            vec![
                ZERO.clone(),
                ZERO.clone(),
                ZERO.clone(),
                ZERO.clone(),
            ],
            vec![
                Instruction::Gate(HADAMARD.clone(), 0),
                Instruction::Gate(CNOT.clone(), 0),
                Instruction::Gate(HADAMARD.clone(), 2),
                Instruction::Gate(CNOT.clone(), 2),
                Instruction::Measure(0, true),
                Instruction::Measure(1, true),
                Instruction::Measure(2, true),
                Instruction::Measure(3, true)
            ]
        );

        println!("{}", &experiment.average_out_pretty(1000));
    }

    #[test]
    fn deutsch_oracle() {
        let circuits = vec![
            ("Constant-0", Instruction::Circuit(vec![], 0)),
            ("Constant-1", Instruction::Circuit(vec![
                Instruction::Gate(PAULIX.clone(), 1)
            ], 0)),
            ("Identity", Instruction::Circuit(vec![
                Instruction::Gate(CNOT.clone(), 0)
            ], 0)),
            ("Negation", Instruction::Circuit(vec![
                Instruction::Gate(CNOT.clone(), 0),
                Instruction::Gate(PAULIX.clone(), 1)
            ], 0)),
        ];

        for circuit in circuits.into_iter() {
            let title = circuit.0;

            let experiment = Experiment::new(vec![
                ZERO.clone(),
                ZERO.clone()
            ], vec![
                Instruction::Gate(PAULIX.clone(), 0),
                Instruction::Gate(PAULIX.clone(), 1),
                Instruction::Gate(HADAMARD.clone(), 0),
                Instruction::Gate(HADAMARD.clone(), 1),
                circuit.1,
                Instruction::Gate(HADAMARD.clone(), 0),
                Instruction::Gate(HADAMARD.clone(), 1),
                Instruction::Measure(0, true),
                Instruction::Measure(1, true),
            ]);

            let measurements = experiment.run().1;

            println!("{} is {}", title, if measurements[0].1 { "constant" } else { "variable" });
        }
    }

    #[test]
    fn teleportation() {
        let qubit = ONE.clone(); //vec![eq!(1/2), eq!((3/4)^0.5)];

        // norm(&mut qubit);

        println!("Initial qubit: {}", print_tensor(&qubit));

        let experiment = Experiment::new(
            vec![
                qubit,
                ZERO.clone(),
                ZERO.clone()
            ],
            vec![
                Instruction::Gate(HADAMARD.clone(), 1),
                Instruction::Gate(CNOT.clone(), 1),
                Instruction::Gate(CNOT.clone(), 0),
                Instruction::Gate(HADAMARD.clone(), 0),
                Instruction::Measure(1, false),
                Instruction::Measure(0, false),
                Instruction::Dependent(Box::new(Instruction::Gate(PAULIX.clone(), 2)), 0),
                Instruction::Dependent(Box::new(Instruction::Gate(PAULIZ.clone(), 2)), 1),
                Instruction::Measure(2, true)
            ]
        );

        println!("{}", &experiment.average_out_pretty(100));
    }

    #[test]
    fn teleport_with_circuit() {
        let mut qubit = ONE.clone();

        pass_gate(&mut qubit, &HADAMARD);

        println!("Input: {}", print_tensor(&qubit));

        let experiment = Experiment::new(vec![
            qubit.clone(),
            ZERO.clone(),
            ZERO.clone()
        ], vec![
            Instruction::Circuit(TELEPORT.clone(), 0),
            Instruction::Measure(2, true)
        ]);

        println!("{}", experiment.average_out_pretty(1000));
    }

    #[test]
    fn teleport_entangled() {
        let experiment = Experiment::new(vec![
            ZERO.clone(),
            ZERO.clone(),
            ZERO.clone(),
            ZERO.clone(),
        ], vec![
            Instruction::Circuit(ENTANGLE.clone(), 0),
            Instruction::Circuit(TELEPORT.clone(), 1),
            Instruction::Measure(0, true),
            Instruction::Measure(3, true)
        ]);

        println!("{}", experiment.average_out_pretty(1000));
    }

    #[test]
    fn teleport_angled() {
        let experiment = Experiment::new(vec![
            ZERO.clone(),
            ZERO.clone(),
            ZERO.clone(),
        ], vec![
            Instruction::MeasureAtAngle((90.0 as f64).to_radians(), 0, false),
            Instruction::Circuit(TELEPORT.clone(), 0),
            Instruction::Measure(2, true)
        ]);

        println!("{}", experiment.average_out_pretty(10000));
    }

    #[test]
    fn measure_at_angle() {
        let angle: f64 = 0.0;

        let experiment = Experiment::new(vec![
            ZERO.clone()
        ], vec![
            Instruction::MeasureAtAngle(angle.to_radians(), 0, true)
        ]);

        println!("{}", experiment.average_out_pretty(10000));
    }

    #[test]
    fn entangle_angle() {
        let angle1: f64 = 45.0;
        let angle2: f64 = 90.0;

        let experiment = Experiment::new(vec![
            ZERO.clone(),
            ZERO.clone(),
        ], vec![
            Instruction::Circuit(ENTANGLE.clone(), 0),
            Instruction::MeasureAtAngle(angle1.to_radians(), 0, true),
            Instruction::MeasureAtAngle(angle2.to_radians(), 0, true),
        ]);

        println!("{}", experiment.average_out_pretty(10000));
    }

    #[test]
    fn gamification_of_bells_theorem() { // youtube.com/watch?v=v7jctqKsUMA
        let mut rng = rand::thread_rng();
        let mut correct = 0.0;
        let total = 10000;

        for _ in 0..total {
            let a = rng.gen::<bool>();
            let b = rng.gen::<bool>();

            let angle1: f64 = if a { 90.0  } else { 0.0   };
            let angle2: f64 = if b { 135.0 } else { 215.0 };

            let measurements = Experiment::new(vec![
                ONE.clone(),
                ONE.clone(),
            ], vec![
                Instruction::Circuit(ENTANGLE.clone(), 0),
                Instruction::MeasureAtAngle(angle1.to_radians(), 0, true),
                Instruction::MeasureAtAngle(angle2.to_radians(), 1, true),
            ]).run().1;

            match (a, b, measurements[0].1, measurements[1].1) {
                (true, true, x, y) => if x != y { correct += 1.0 },
                (_   , _   , x, y) => if x == y { correct += 1.0 }
            }
        }

        println!("win %: {:.2}%", 100.0 * correct / (total as f64));
    }
}