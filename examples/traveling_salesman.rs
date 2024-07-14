use std::{fmt::Display, fs::File, io::Read};

use genalg::{Chromosome, Run};
use rand::{seq::SliceRandom, thread_rng, Rng};

//const NUM_CITIES: usize = 20;
const NUM_CITIES: usize = 6;
const N_POP: usize = 20;
const N_KEEP: usize = 10;

#[derive(Debug, Clone)]
struct Permutation {
    // Permutation of cities 1..n - the order of cities to visit
    perm: [usize; NUM_CITIES],
}

#[derive(Debug, Clone)]
struct DistMatrix {
    dist_mat: [f64; NUM_CITIES * NUM_CITIES],
}

impl DistMatrix {
    fn new() -> DistMatrix {
        // In case we want to read the distance matrix from a file
        //let mut dist_mat = vec![];
        //let mut dist_mat_s = String::new();
        //let mut dist_mat_f = File::open("./dist_mat").unwrap();

        //dist_mat_f.read_to_string(&mut dist_mat_s).unwrap();

        //for line in dist_mat_s.lines() {
        //    // Each line is composed of NUM_CITIES=20 numbers
        //    for number in line.split(" ") {
        //        if !number.is_empty() {
        //            dist_mat.push(number.parse::<f64>().unwrap());
        //        }
        //    }
        //}
        let dist_mat = vec![1f64,2f64,3f64,9f64,27f64,8f64,23f64,15f64,11f64,18f64,13f64,24f64,34f64,2f64,29f64,33f64,14f64,32f64,25f64,35f64,20f64,31f64,6f64,10f64,12f64,17f64,5f64,16f64,3f64,21f64,4f64,30f64,36f64,22f64,26f64,1f64];

        println!("{}", dist_mat.len());

        DistMatrix {
            dist_mat: dist_mat.as_slice().try_into().unwrap(),
        }
    }

    // Calculate the distance from city a to city b (and vice versa)
    // given the NUM_CITIESxNUM_CITIES square distance matrix
    fn dist(&self, city_a: usize, city_b: usize) -> f64 {
        self.dist_mat[city_a * NUM_CITIES + city_b]
    }
}

impl<'a> Chromosome for Permutation {
    fn fitness(&self) -> f64 {
        let mut cost = 0f64;
        let dist_mat = DistMatrix::new();

        for i in 0..NUM_CITIES - 1 {
            let city_a = self.perm[i];
            let city_b = self.perm[i + 1];

            cost += dist_mat.dist(city_a, city_b);
        }

        // Getting from the last city to the first city
        cost += dist_mat.dist(self.perm[NUM_CITIES - 1], self.perm[0]);

        -1f64 * cost as f64
    }

    fn crossover(x: Self, y: Self) -> (Self, Self) {
        let mut first_offspring = x.perm;
        let mut second_offspring = y.perm;
        // Cycle crossover
        let mut rng = thread_rng();
        let mut visited = vec![];
        // Exchange a random index
        let random_idx = rng.gen_range(0..NUM_CITIES);
        let tmp = first_offspring[random_idx];
        // The current element which is a duplicate in the first offspring
        let mut dup_elem;
        first_offspring[random_idx] = second_offspring[random_idx];
        dup_elem = second_offspring[random_idx];
        second_offspring[random_idx] = tmp;
        visited.push(random_idx);

        loop {
            for i in 0..NUM_CITIES {
                // If this is the duplicate and we haven't touched it yet
                if first_offspring[i] == dup_elem && !visited.contains(&i) {
                    let tmp = first_offspring[i];
                    first_offspring[i] = second_offspring[i];
                    dup_elem = second_offspring[i];
                    second_offspring[i] = tmp;
                    visited.push(i);
                }
            }

            // If there are no duplicates
            if first_offspring.iter().sum::<usize>() == (0..NUM_CITIES).sum() {
                break;
            }
        }

        let first_offspring = Permutation {
            perm: first_offspring,
        };
        let second_offspring = Permutation {
            perm: second_offspring,
        };

        (first_offspring, second_offspring)
    }

    fn mutate(&self, mutations: &[bool]) -> Self {
        let is_mut = mutations[0];
        let mut rng = thread_rng();
        let perm = self.perm;

        if is_mut {
            let mut mutated_perm = perm;
            // Generate two random indexes
            let rand_idx_1 = rng.gen_range(0..NUM_CITIES);
            let rand_idx_2 = rng.gen_range(0..NUM_CITIES);
            // Swap between them
            let tmp = mutated_perm[rand_idx_1];
            mutated_perm[rand_idx_1] = mutated_perm[rand_idx_2];
            mutated_perm[rand_idx_2] = tmp;

            Permutation { perm: mutated_perm }
        } else {
            Permutation { perm }
        }
    }

    fn random() -> Self {
        let mut rng = thread_rng();
        // Create a random permutation
        let mut perm: Vec<usize> = (0..NUM_CITIES).collect();
        perm.shuffle(&mut rng);
        // Return the Permutation
        Permutation {
            perm: perm.try_into().unwrap(),
        }
    }
}

impl Display for Permutation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for x in self.perm {
            write!(f, "{} ", x).unwrap()
        }

        return Ok(());
    }
}

fn main() {
    let mut run =
        Run::<Permutation>::new(200, 0.2, N_KEEP, N_POP, 1, genalg::Selection::UpperItems);

    let history = run.begin();

    for (perm, fitness) in history {
        println!("{}\t{}", perm, fitness);
    }
}
