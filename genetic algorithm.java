import java.util.*;
import java.io.*;

//genetic Algorithm
class Gene {
	private int number;

	public Gene(int number) {
		this.number = number;
	}

	public void setNumber(int number) {
		this.number = number;
	}

	public int getNumber() {
		return number;
	}

	@Override
	public String toString() {
		return String.valueOf(this.getNumber());
	}
}

class Individual {
	private Gene[] chromosome;
	private double fitness = -1;

	public Individual(Gene[] chromosome) {
		this.chromosome = chromosome;
	}

	public Individual(int chromosomeSize) {
		this.chromosome = new Gene[chromosomeSize];
		this.initRandomIndividual();
	}

	public Gene[] getChromosome() {
		return this.chromosome;
	}

	public void setFitness(double fitness) {
		this.fitness = fitness;
	}

	public double getFitness() {
		return fitness;
	}

	public void setGene(int index, Gene gene) {
		this.chromosome[index] = gene;
	}

	public Gene getGene(int index) {
		return this.chromosome[index];
	}

	private void initRandomIndividual() {
		if (this.chromosome == null)
			throw new NullPointerException("Chromosome is null !");
		for (int i = 0; i < this.chromosome.length; i++) {
			double rand = Math.random();
			this.setGene(i, new Gene(rand > 0.5 ? 1 : 0));
		}
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (Gene gene : this.chromosome) {
			sb.append(gene.toString());
		}
		return sb.toString();
	}
}

class Population {
	private Individual[] individuals;
	private double fitness = -1;

	public Population(int populationSize) {
		this.individuals = new Individual[populationSize];
	}

	public Population(int populationSize, int chromosomeSize) {
		this.individuals = new Individual[populationSize];

		for (int i = 0; i < populationSize; i++) {
			this.individuals[i] = new Individual(chromosomeSize);
		}
	}

	public Individual[] getIndividuals() {
		return individuals;
	}

	public double getFitness() {
		return fitness;
	}

	public void setFitness(double fitness) {
		this.fitness = fitness;
	}

	public Individual getIndividual(int index) {
		return this.individuals[index];
	}

	public void setIndividual(int index, Individual individual) {
		this.individuals[index] = individual;
	}

	public Individual getFittestIndividual(int offset) {
		// sort the population so that first individual has highest fitness score
		Arrays.sort(this.individuals, new Comparator<Individual>() {
			@Override
			public int compare(Individual a, Individual b) {
				if (a.getFitness() < b.getFitness())
					return 1;
				if (a.getFitness() > b.getFitness())
					return -1;
				return 0;
			}
		});
		return this.individuals[offset];
	}

	public void sufflePopulation() {
		List<Individual> list = Arrays.<Individual>asList(this.individuals);
		Collections.shuffle(list);
		list.toArray(this.individuals);
	}

}

class Algorithm {
	private int popolationSize;
	private double mutationRate, crossoverRate;
	private int elitismCount;

	public Algorithm(int populationSize, double mutationRate, double crossoverRate, int elitismCount) {
		this.popolationSize = populationSize;
		this.mutationRate = mutationRate;
		this.crossoverRate = crossoverRate;
		this.elitismCount = elitismCount;
	}

	public Population initPopulation(int chromosomeLength) {
		return new Population(this.popolationSize, chromosomeLength);
	}

	public void evaluateFitness(Individual individual) {
		int correct = 0;
		for (int i = 0; i < individual.getChromosome().length; i++) {
			if (individual.getGene(i).getNumber() == 1) {
				correct++;
			}
		}
		double fitness = (double) correct / individual.getChromosome().length;
		individual.setFitness(fitness);
	}

	public void evaluateFitness(Population population) {
		double fitness = 0.0;
		for (Individual individual : population.getIndividuals()) {
			evaluateFitness(individual);
			fitness += individual.getFitness();
		}
		population.setFitness(fitness);
	}

	public boolean solutionFound(Population population) {
		return population.getFittestIndividual(0).getFitness() == 1;
	}

	public Individual selectParentUsingRouletteWheel(Population population) {
		Individual[] individuals = population.getIndividuals();
		double wheelPosition = population.getFitness() * Math.random();
		double wheel = 0;
		for (Individual individual : individuals) {
			wheel += individual.getFitness();
			if (wheel >= wheelPosition) {
				return individual;
			}
		}
		return individuals[individuals.length - 1];
	}

	public Population crossover(Population population) {
		Population newPopulation = new Population(population.getIndividuals().length);
		for (int i = 0; i < population.getIndividuals().length; i++) {
			Individual firstParent = population.getFittestIndividual(i);
			// do crossover if not elite
			if (crossoverRate > Math.random() && i > elitismCount) {
				// perform crossover
				Individual secondParent = selectParentUsingRouletteWheel(population);
				Individual offspring = new Individual(firstParent.getChromosome().length);

				for (int j = 0; j < firstParent.getChromosome().length; j++) {
					offspring.setGene(j, Math.random() < 0.5 ? firstParent.getGene(j) : secondParent.getGene(j));
				}
				newPopulation.setIndividual(i, offspring);
			} else {
				// no crossover
				newPopulation.setIndividual(i, firstParent);
			}
		}
		return newPopulation;
	}

	public Population mutate(Population population) {
		Population newPopulation = new Population(population.getIndividuals().length);
		for (int i = 0; i < population.getIndividuals().length; i++) {
			Individual individual = population.getFittestIndividual(i);
			for (int j = 0; j < individual.getChromosome().length; j++) {
				// skip for elite individual
				if (i > elitismCount && Math.random() > mutationRate) {
					Gene mutatedGene = new Gene(1);
					if (individual.getGene(j).getNumber() == 1) {
						mutatedGene.setNumber(0);
					}
					individual.setGene(j, mutatedGene);
				}
			}
			newPopulation.setIndividual(i, individual);
		}
		return newPopulation;
	}
}

public class Main {

	public static void main(String[] args) {
		Algorithm alg = new Algorithm(100, 0.001, 0.95, 2);
		Population population = alg.initPopulation(25);
		alg.evaluateFitness(population);
		int generation = 1;

		while (!alg.solutionFound(population)) {
			// evaluate fitness
			alg.evaluateFitness(population);

			System.out.println(generation + ": " + population.getFittestIndividual(0) + ": "
					+ population.getFittestIndividual(0).getFitness());

			// crossover
			population = alg.crossover(population);
			// mutation
			population = alg.mutate(population);

			generation++;
		}

		System.out.println("-------------------------------");
		System.out.printf("FOUND SOLUTION IN %d GENERATIONS !\n", generation);
		System.out.println("SOLUTION: " + population.getFittestIndividual(0).toString() + ": "
				+ population.getFittestIndividual(0).getFitness());
	}
}
