use std::fmt::Display;
use genalg::{Chromosome, Run};
use rand::{thread_rng, Rng};

const N_KEEP: usize = 6;
const N_POP: usize = 12;

// A two dimensional point
#[derive(Clone)]
struct Point {
    x: f64,
    y: f64,
}

impl Chromosome for Point {
    fn crossover(x: Self, y: Self) -> (Self, Self) {
        let mut rng = thread_rng();
        // Do we change x or y?
        let param = rng.gen_range(0..2);

        match param {
            // Change x coordinate
            0 => {
                let child1 = Point {
                    x: y.x,
                    y: x.y,
                };
                let child2 = Point {
                    x: x.x,
                    y: y.y,
                };

                (child1, child2)
            }
            // Change y coordinate
            1 => {
                let child1 = Point {
                    x: x.x,
                    y: y.y,
                };
                let child2 = Point {
                    x: y.x,
                    y: x.y,
                };

                (child1, child2)
            }
            _ => {
                println!("An error occured!!!");

                (x, y)
            }
        }
    }

    fn fitness(&self) -> f64 {
        // The function we want to maximize
        0.8 * self.x.sin() * (2f64 * self.y).cos() -3f64 * (self.x - self.y).cos()
    }
    // Mutate each of the parameters of the chromosome
    fn mutate(&self, mutations: &[bool]) -> Self {
        let mut rng = thread_rng();

        // If a mutation happened, we replace the mutated coordinate with a random number from 0 to 10
        Point {
            x: if mutations[0] == true {
                rng.gen_range(0f64..10f64)
            } else {
                self.x
            },
            y: if mutations[1] == true {
                rng.gen_range(0f64..10f64)
            } else {
                self.y
            },
        }
    }

    fn random() -> Self {
        let mut rng = thread_rng();
        let x = rng.gen_range(0f64..10f64);
        let y = rng.gen_range(0f64..10f64);

        Point { x, y }
    }
}

impl Display for Point {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
    }
}

fn main() {
    let mut run = Run::<Point>::new(200, 0.2, N_KEEP, N_POP, 2, genalg::Selection::UpperItems);
    let history = run.begin();

    for i in 0..200 {
        println!("{}\t{}", history[i].0, history[i].1);
    }
}
