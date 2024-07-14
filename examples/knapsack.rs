use std::fmt::Display;

use genalg::{Chromosome, Run};
use rand::{
    distributions::{Bernoulli, Distribution},
    thread_rng
};

// Total number of items
const NUM_ITEMS: usize = 12;
// The maximum weight we can carry (kg)
const MAX_WEIGHT: f64 = 20f64;
// Size of population
const N_POP: usize = 20;

#[derive(Clone)]
struct Knapsack {
    // True if we include the item in the knapsack and
    // false otherwise
    items: [bool; NUM_ITEMS],
}

struct KnapsackData {
    // The value of each item
    values: [f64; NUM_ITEMS],
    // The weight of each item
    weights: [f64; NUM_ITEMS],
}

impl KnapsackData {
    fn new() -> KnapsackData {
        let values = [
            3f64, 2f64, 4f64, 3.5f64, 7.2f64, 5.9f64, 0.7f64, 8f64, 6f64, 10f64, 1f64, 3.1f64,
        ];
        let weights = [
            8f64, 2f64, 1f64, 3.5f64, 5.7f64, 3f64, 8f64, 9f64, 2.3f64, 5.5f64, 6f64, 0.2f64,
        ];

        KnapsackData { values, weights }
    }
}

impl Chromosome for Knapsack {
    fn crossover(x: Self, y: Self) -> (Self, Self) {
        let mut new_x = Knapsack { items: x.items };
        let mut new_y = Knapsack { items: y.items };

        // Exchange 2/3 of the bits
        for i in 0..2 * (NUM_ITEMS / 3) {
            let tmp = new_x.items[i];
            new_x.items[i] = new_y.items[i];
            new_y.items[i] = tmp;
        }

        (new_x, new_y)
    }

    fn fitness(&self) -> f64 {
        // Sum all of the weight values
        // Stop summing and return a negative value if we exceed the maximum weight allowed
        let mut value_sum = 0f64;
        let mut weight_sum = 0f64;
        // Knapsack data (weights and values)
        let data = KnapsackData::new();
        let (values, weights) = (data.values, data.weights);

        for i in 0..NUM_ITEMS {
            // If the knapsack contains the item
            if self.items[i] {
                value_sum += values[i];
                weight_sum += weights[i];
            }

            // We exceeded the maximum allowed weight
            if weight_sum > MAX_WEIGHT {
                return -1f64;
            }
        }

        value_sum
    }

    fn mutate(&self, mutations: &[bool]) -> Self {
        let mut new_knapsack = Knapsack { items: self.items };

        for i in 0..NUM_ITEMS {
            // If this bit has mutated, flip it
            if mutations[i] {
                new_knapsack.items[i] = !new_knapsack.items[i];
            }
        }

        new_knapsack
    }

    // Generate a random (possibly weight-exceeding) knapsack
    fn random() -> Self {
        let mut rng = thread_rng();
        let dist = Bernoulli::new(0.5).unwrap();
        // Generate a vector of NUM_ITEMS random booleans
        let items: Vec<bool> = dist.sample_iter(&mut rng).take(NUM_ITEMS).collect();

        Knapsack {
            items: items.try_into().unwrap(),
        }
    }
}

impl Display for Knapsack {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..NUM_ITEMS {
            if self.items[i] {
                write!(f, "{} ", i).unwrap();
            }
        }

        Ok(())
    }
}

fn main() {
    let mut run = Run::<Knapsack>::new(200, 0.05, 0, N_POP, NUM_ITEMS, genalg::Selection::Softmax);
    let history = run.begin();

    for (knapsack, value) in history {
        println!("{}\t{}", knapsack, value);
    }
}
