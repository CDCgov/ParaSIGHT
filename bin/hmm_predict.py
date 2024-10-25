#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import json

class HMM:
    def __init__(self, X, A_prime, B_prime, start_probs, obs_states):
        self.X = X  # Hidden states (species, subtype, host, etc.)
        self.A_prime = A_prime  # Transition matrix
        self.B_prime = B_prime  # Emission matrix
        self.start_probs = start_probs  # Initial probabilities
        self.obs_states = obs_states  # Set of observed states (alleles)

    # Function to construct transition matrix with pseudocounts
    @staticmethod
    def construct_A_prime(allele_transitions, pseudocount=1):
        num_hidden_states = len(allele_transitions)
        A_prime = np.zeros((num_hidden_states, num_hidden_states))

        for i in range(num_hidden_states):
            row_total = sum(allele_transitions[i]) + pseudocount * num_hidden_states
            A_prime[i] = [(count + pseudocount) / row_total for count in allele_transitions[i]]

        return A_prime

    # Function to construct emission matrix from observed allele frequencies
    @staticmethod
    def construct_B_prime(allele_freqs):
        num_hidden_states = len(allele_freqs)
        num_obs_states = len(allele_freqs[0])
        B_prime = np.zeros((num_hidden_states, num_obs_states))

        for i in range(num_hidden_states):
            total = sum(allele_freqs[i])
            B_prime[i] = [freq / total for freq in allele_freqs[i]]

        return B_prime

    # Viterbi algorithm for computing the most likely sequence of hidden states
    def viterbi(self, observed_seq):
        T = len(observed_seq)
        N = len(self.X)

        viterbi = np.zeros((N, T))
        backpointer = np.zeros((N, T), dtype=int)

        # Initialize base cases (t == 0)
        for s in range(N):
            viterbi[s, 0] = self.start_probs[s] * self.B_prime[s, self.obs_states.index(observed_seq[0])]
            backpointer[s, 0] = 0

        # Run Viterbi for t > 0
        for t in range(1, T):
            for s in range(N):
                prob_state = [viterbi[s_prev, t-1] * self.A_prime[s_prev, s] * self.B_prime[s, self.obs_states.index(observed_seq[t])] for s_prev in range(N)]
                viterbi[s, t] = max(prob_state)
                backpointer[s, t] = np.argmax(prob_state)

        # Backtrack to find the best path
        best_path = np.zeros(T, dtype=int)
        best_path[-1] = np.argmax(viterbi[:, T-1])
        for t in range(T-2, -1, -1):
            best_path[t] = backpointer[best_path[t+1], t+1]

        return [self.X[state] for state in best_path], viterbi.max(axis=0)

# Example input data for constructing A_prime and B_prime
allele_transitions = [
    [5, 1, 3],  # Transitions for species
    [2, 3, 1],  # Transitions for subtype
    [1, 3, 5],  # Transitions for host
]

allele_freqs = [
    [0.6, 0.4],  # Emissions for species
    [0.7, 0.3],  # Emissions for subtype
    [0.2, 0.8],  # Emissions for host
]

# Sets of hidden states (e.g., species, subtype, host, etc.)
X_species   = ['species', 'subtype', 'host']
X_subtype   = []
X_host      = [ 'Human' ]
X_geostate  = [ 'AL','AK', 'AZ', 'AR', 'CA', 'CO', 'CT', 'DE',
                'FL', 'GA', 'HI', 'ID', 'IL', 'IN', 'IA', 'KS',
                'KY', 'LA', 'ME', 'MD', 'MA', 'MI', 'MN', 'MS',
                'MO', 'MT', 'NE', 'NV', 'NH', 'NJ', 'NM', 'NY',
                'NC', 'ND', 'OH', 'OK', 'OR', 'PA', 'RI', 'SC',
                'SD', 'TN', 'TX', 'UT', 'VT', 'VA', 'WA', 'WV',
                'WI', 'WY', 'DC'
            ]

# Starting probabilities (allele frequencies for the first gene)
start_probs = [0.5, 0.3, 0.2]

# Observed states (possible alleles)
obs_states = ['allele1', 'allele2']

# Construct the transition (A_prime) and emission (B_prime) matrices
A_prime = HMM.construct_A_prime(allele_transitions, pseudocount=1)
B_prime = HMM.construct_B_prime(allele_freqs)

# Save the matrices to disk
with open('A_prime.json', 'w') as f:
    json.dump(A_prime.tolist(), f)
with open('B_prime.json', 'w') as f:
    json.dump(B_prime.tolist(), f)

# Define a set of observed alleles from an unknown genome
observed_seq = ['allele1', 'allele2', 'allele1']

# Create HMM model
hmm = HMM(X, A_prime, B_prime, start_probs, obs_states)

# Run the Viterbi algorithm to predict hidden states
best_path, viterbi_probs = hmm.viterbi(observed_seq)

print("Most likely hidden states:", best_path)
print("Viterbi probabilities:", viterbi_probs)


if __name__ == "__main__":
    # Argument parser for input and output files
    parser = argparse.ArgumentParser(description="Run Viterbi algorithm over a state transition and emission matrix to determine most probable labels.")
    parser.add_argument("-a", "--transitions", required=True, help="Input hidden state transitions CSV file")
    parser.add_argument("-b", "--emissions", required=True, help="Input observed/emission probabilities CSV file")
    args = parser.parse_args()

    # Run the main function
    main(args.input, args.output)