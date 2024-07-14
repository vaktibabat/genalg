use rand::{
    distributions::{Bernoulli, Distribution, Uniform},
    rngs::ThreadRng,
    thread_rng,
};

// Types of ways to generate the probability distribution that controls
// which chromosomes are more likely to reproduce
pub enum Selection {
    // Fitness of the current chromosome / Sum of fitnesses
    Fraction,
    // Only keep the upper N_KEEP items of the population and let them reproduce
    // with probabilities proportional to their place (i.e. first place is most likely to reproduce)
    UpperItems,
    // Softmax function
    Softmax,
}

// How to treat the structure as a chromosome? i.e. how should it reproduce
// and how likely is it to be selected for reproduction
pub trait Chromosome
where
    Self: Sized,
{
    // How should two chromosomes of this type be crossed together to produce two new offspring?
    fn crossover(x: Self, y: Self) -> (Self, Self);
    // Create a new random chromosome
    fn random() -> Self;
    // How fit is this chromosome?
    fn fitness(&self) -> f64;
    // How should this chromosome mutate?
    fn mutate(&self, mutations: &[bool]) -> Self;
}

// A run of a genetic algorithm
pub struct Run<T: Chromosome> {
    // How many generations do we want this algorithm to go on for?
    num_iter: usize,
    // What is the probability of mutation (i.e. randomly bit flipping after a crossover)?
    p_mutation: f64,
    // Do we want to always keep the fitter part of the population?
    n_keep: usize,
    // What is the population size?
    n_pop: usize,
    // How many parameters does each chromosome have? For example, 2D points have two parameters: x and y
    n_params: usize,
    // How do we select which chromosomes are likely to reproduce?
    selection: Selection,
    // The current generation
    gen: Vec<T>,
    // Best-so-far vector of the optimal values of the cost function along with the chromosomes
    // that yielded them
    best_so_far: Vec<(T, f64)>,
}

impl<T: Chromosome + Clone> Run<T> {
    // Initialize the run
    pub fn new(
        num_iter: usize,
        p_mutation: f64,
        n_keep: usize,
        n_pop: usize,
        n_params: usize,
        selection: Selection,
    ) -> Run<T> {
        let gen = vec![];
        let best_so_far = vec![];

        Run {
            num_iter,
            p_mutation,
            n_keep,
            n_pop,
            n_params,
            selection,
            gen,
            best_so_far,
        }
    }
    // Run the algorithm for the specified number of iterations
    pub fn begin(&mut self) -> Vec<(T, f64)> {
        // Initialize the first generation by generating N_pop random chromosomes
        for _ in 0..self.n_pop {
            self.gen.push(T::random());
        }
        // Used later for sampling random chromosomes
        let mut rng = thread_rng();
        let uni = Uniform::new(0f64, 1f64);
        let mut_dist = Bernoulli::new(self.p_mutation).unwrap();

        for _ in 0..self.num_iter {
            // Sort the chromosomes in the current generation by fitness
            let mut chromo_by_fit: Vec<(&T, f64)> =
                self.gen.iter().map(|c| (c, c.fitness())).collect();
            chromo_by_fit.sort_by(|a, b| a.1.total_cmp(&b.1));
            // Put them in descending order
            chromo_by_fit.reverse();
            // Push the fittest chromosome in the current generation to the best_so_far vector
            let best_clone = chromo_by_fit.first().unwrap();
            self.best_so_far.push((best_clone.0.clone(), best_clone.1));
            // Generate the new generation:
            // 1. Take the fittest N_keep chromosomes in the current generation and add them to the next one
            // 2. Generate the rest through natural selection
            let mut next_gen = vec![];
            // Step 1 - The fittest N_keep chromosomes are the last N_keep in chromo_by_fit
            for (c, _) in &chromo_by_fit[0..self.n_keep] {
                next_gen.push((*c).clone());
            }
            // Step 2 - The rest of the chromosomes are generated through natural selection
            for _ in 0..(self.n_pop - self.n_keep) / 2 {
                let parent1 = self
                    .random_parent(uni, &mut rng, &chromo_by_fit)
                    .expect("Cannot sample a random parent from an empty population!");
                let parent2 = self
                    .random_parent(uni, &mut rng, &chromo_by_fit)
                    .expect("Cannot sample a random parent from an empty population!");
                // Perform a crossover
                let (child1, child2) = T::crossover(parent1, parent2);
                // Introduce mutations randomly
                let mutations: Vec<bool> = mut_dist.sample_iter(&mut rng).take(self.n_params).collect();
                let child1 = child1.mutate(mutations.as_slice());
                let child2 = child2.mutate(mutations.as_slice());

                // Push to the next generation
                next_gen.push(child1);
                next_gen.push(child2);
            }

            // Assign it to the next gen
            self.gen = next_gen;
        }

        self.best_so_far.clone()
    }

    // Sample a random chromosome from this generation WRT the selection distribution
    // This is done by sampling a random number from 0 to 1 from a uniform distribution
    // and then picking the chromosome corresponding to that interval
    // for example, if we have P(C1) = 0.5 P(C2) = 0.25 and P(C3) = 0.25
    // If the random number y is smaller than 0.5, we pick C1
    // If it's between 0.5 and 0.5 + 0.25 = 0.75, we pick C2
    // and if it's larger than 0.75, we pick C3
    fn random_parent(
        &self,
        dist: Uniform<f64>,
        rng: &mut ThreadRng,
        chromo_by_fit: &Vec<(&T, f64)>,
    ) -> Option<T> {
        // The chromosomes by probability to be reproduced
        let mut chromo_by_prob: Vec<(&T, f64)> = chromo_by_fit
            .iter()
            .enumerate()
            .map(|(i, (a, _))| (*a, match self.selection {
                // Different ways to compute the selection probability
                Selection::Fraction => {
                    let fit_sum = chromo_by_fit.iter().map(|(_, fit)| fit).sum::<f64>();
                    
                    a.fitness() / fit_sum
                },
                Selection::UpperItems => {
                    let place = i + 1;
                    
                    if 1 <= place && place <= self.n_keep {
                        (self.n_keep as f64 - place as f64 + 1f64)
                            / ((1..self.n_keep + 1).sum::<usize>() as f64)
                    } else {
                        0f64
                    }
                }
                Selection::Softmax => {
                    let fit_exp_sum = chromo_by_fit.iter().map(|(_, fit)| fit.exp()).sum::<f64>();

                    a.fitness().exp() / fit_exp_sum
                }
            }))
            .collect();

        chromo_by_prob.sort_by(|a, b| a.1.total_cmp(&b.1));
        chromo_by_prob.reverse();
        // Generate a random number y
        let y = dist.sample(rng);
        let mut prev = 0f64;
        let mut curr = 0f64;

        for chromo in chromo_by_prob {
            curr += chromo.1;

            if prev <= y && y < curr {
                return Some(chromo.0.clone());
            }

            prev = curr;
        }

        None
    }
}
